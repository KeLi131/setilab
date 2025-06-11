#define _GNU_SOURCE
#include <assert.h>
#include <ctype.h>
#include <pthread.h> // pthread api
#include <sched.h>   // for processor affinity
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // unix standard apis

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40 // Max number of printed *
#define THRESHOLD 2.0
#define ALIENS_LOW 50000.0
#define ALIENS_HIGH 150000.0

typedef struct {
  double Fs;
  double Fcl;
  double Fch;
  int order;
  double *coeffs;
} generate_band_pass_args_t;

typedef struct {
  int order;
  double *coeffs;
} hamming_window_args_t;

typedef struct {
  int length;
  double *input_signal;
  int order;
  double *coeffs;
  double *power;
} convolve_and_compute_power_args_t;

typedef struct {
  // generate_band_pass_args_t *generate_band_pass_args;
  // hamming_window_args_t *hamming_window_args;
  // convolve_and_compute_power_args_t *convolve_and_compute_power_args;
  long myid;
  signal *sig;
  int filter_order;
  double bandwidth;
  double *band_power;
  int num_band_threads;
  int num_conv_threads;
} band_thread_args_t;

typedef struct {
  double *input_signal;
  int start_sample;
  int end_sample;
  int order;
  double *coeffs;
  double *partial_power_sum;
  int total_length;
} conv_thread_args_t;

int num_threads;    // number of threads we will use
int num_processors; // number of processors we will use
int num_bands;
double *band_power;
pthread_t *tid; // array of thread ids

void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands "
         "num_threads num_processors\n");
}

double avg_power(double *data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double *data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double *data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double *data, int num) {

  double dc = avg_of(data, num);

  printf("Removing DC component of %lf\n", dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}

void *convolution_worker(void *arg) {
  conv_thread_args_t *conv_thread_args = (conv_thread_args_t *)arg;

  double local_pow_sum = 0;

  for (int i = conv_thread_args->start_sample; i < conv_thread_args->end_sample;
       i++) {
    double cur_sum = 0;
    for (int j = conv_thread_args->order; j >= 0; j--) {
      if ((i - j) >= 0 && (i - j) < conv_thread_args->total_length) {
        cur_sum +=
            conv_thread_args->input_signal[i - j] * conv_thread_args->coeffs[j];
      }
    }
    local_pow_sum += cur_sum * cur_sum;
  }

  *(conv_thread_args->partial_power_sum) = local_pow_sum;
  free(conv_thread_args);
  pthread_exit(NULL);
}

int p_convolve_and_compute_power(int length, double input_signal[], int order,
                                 double coeffs[], double *power,
                                 int num_conv_threads) {

  // if 1 or 0 threads available, reuse old solution
  if (num_conv_threads <= 1) {
    printf("parallel bands; sequential convolutions\n");
    convolve_and_compute_power(length, input_signal, order, coeffs, power);
    return 0;
  }

  // to parallelize

  pthread_t *tid = malloc(num_conv_threads * sizeof(pthread_t));
  double *partial_sums = malloc(num_conv_threads * sizeof(double));

  int samples_per_thread = length / num_conv_threads;

  for (int t = 0; t < num_conv_threads; t++) {
    conv_thread_args_t *args = malloc(sizeof(conv_thread_args_t));

    args->input_signal = input_signal;
    args->start_sample = t * samples_per_thread;
    args->end_sample =
        (t == num_conv_threads - 1) ? length : (t + 1) * samples_per_thread;
    args->order = order;
    args->coeffs = coeffs;
    args->partial_power_sum = &partial_sums[t];
    args->total_length = length;

    int ret = pthread_create(&tid[t], NULL, convolution_worker, args);

    if (ret != 0) {
      perror("failed to start thread (conv)");
      exit(-1);
    }
  }

  double total_power = 0;

  for (int t = 0; t < num_conv_threads; t++) {
    int ret = pthread_join(tid[t], NULL);
    if (ret != 0) {
      perror("join failed (conv threads)");
      exit(-1);
    }
    total_power += partial_sums[t];
  }

  *power = total_power / length;

  free(tid);
  free(partial_sums);
  return 0;
}

// Function run by each thread calculating each band
void *band_worker(void *arg) {
  band_thread_args_t *band_thread_args = (band_thread_args_t *)arg;

  long myid = band_thread_args->myid;
  signal *sig = band_thread_args->sig;
  double bandwidth = band_thread_args->bandwidth;
  double *band_power = band_thread_args->band_power;
  int filter_order = band_thread_args->filter_order;
  int num_conv_threads = band_thread_args->num_conv_threads;
  int num_band_threads = band_thread_args->num_band_threads;
  double filter_coeffs[filter_order + 1];

  // put ourselves on the desired processor
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(myid % num_processors, &set);
  if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
    perror("Can't setaffinity");                     // hopefully doesn't fail
    exit(-1);
  }

  // This figures out the chunk of the bin_power I should
  // work on based on my id
  int blocksize = num_bands / num_band_threads;

  // If there are more threads than bands
  if (blocksize < 1)
    blocksize = 1;
  int mystart = myid * blocksize;
  int myend = 0;
  if (myid >= (num_bands - 1)) { // last thread
    // the last thread will take care of the leftover
    // elements of the vector, in case num_threads doesn't
    // divide vector_len
    // WARNING: this is a suboptimal solution. It means that the last thread
    // might do much more work than the other threads (up to almost double)
    // which will slow down the entire job. A better solution would split up
    // remainder work equally between threads...
    myend = num_bands;
  } else {
    myend = (myid + 1) * blocksize;
  }

  if (myend > mystart) {
    for (int i = mystart; i < myend; i++) {

      // Make the filter
      generate_band_pass(band_thread_args->sig->Fs,
                         i * bandwidth + 0.0001, // keep within limits
                         (i + 1) * bandwidth - 0.0001, filter_order,
                         filter_coeffs);

      hamming_window(filter_order, filter_coeffs);

      // Convolve
      p_convolve_and_compute_power(sig->num_samples, sig->data, filter_order,
                                   filter_coeffs, &(band_power[i]),
                                   num_conv_threads);
    }
  }
  // Done.  The master thread will sum up the partial sums
  free(band_thread_args);
  pthread_exit(NULL); // finish - no return value
}

int analyze_signal(signal *sig, int filter_order, int num_bands, double *lb,
                   double *ub) {

  double Fc = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data, sig->num_samples);

  double signal_power = avg_power(sig->data, sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart, THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  band_power = (double *)malloc(num_bands * sizeof(double));

  // distribute threads btwn bands and convolutions?

  int num_band_threads, num_conv_threads;

  num_band_threads = num_threads / 2;
  if (num_band_threads < 1)
    num_band_threads = 1;
  num_conv_threads = num_threads - num_band_threads;

  for (int i = 0; i < num_threads; i++) {

    // generate_band_pass_args_t *generate_band_pass_args =
    // (generate_band_pass_args_t *)malloc(sizeof(generate_band_pass_args_t));
    // generate_band_pass_args->Fs = sig->Fs;
    // generate_band_pass_args->Fcl = band * bandwidth + 0.0001; // keep within
    // limits generate_band_pass_args->Fch = (band + 1) * bandwidth - 0.0001;
    // generate_band_pass_args->order = filter_order;
    // generate_band_pass_args->coeffs = filter_coeffs;

    // band_thread_args->generate_band_pass_args = generate_band_pass_args;

    // hamming_window_args_t *hamming_window_args = (hamming_window_args_t
    // *)malloc(sizeof(hamming_window_args_t)); hamming_window_args->order =
    // filter_order; hamming_window_args->coeffs = filter_coeffs;

    // band_thread_args->hamming_window_args = hamming_window_args;

    // convolve_and_compute_power_args_t *convolve_and_compute_power_args =
    // (convolve_and_compute_power_args_t
    // *)malloc(sizeof(convolve_and_compute_power_args_t));
    // convolve_and_compute_power_args->length = sig->num_samples;
    // convolve_and_compute_power_args->input_signal = sig->data;
    // convolve_and_compute_power_args->order = filter_order;
    // convolve_and_compute_power_args->coeffs = filter_coeffs;
    // convolve_and_compute_power_args->power = &(band_power[band]);

    // band_thread_args->convolve_and_compute_power_args =
    // convolve_and_compute_power_args;

    // Create the input arguments to band_worker function
    band_thread_args_t *band_thread_args =
        (band_thread_args_t *)malloc(sizeof(band_thread_args_t));

    band_thread_args->band_power = band_power;
    band_thread_args->bandwidth = bandwidth;
    band_thread_args->filter_order = filter_order;
    band_thread_args->myid = i;
    band_thread_args->sig = sig;
    band_thread_args->num_conv_threads = num_conv_threads;
    band_thread_args->num_band_threads = num_band_threads;

    int returncode =
        pthread_create(&(tid[i]),   // thread id gets put here
                       NULL,        // use default attributes
                       band_worker, // thread will begin in this function
                       (void *)band_thread_args);
    if (returncode != 0) {
      perror("Failed to start thread");
      exit(-1);
    }
  }

  int i_max = (num_threads < num_bands) ? num_threads : num_bands;
  // now we will join all the threads
  for (int i = 0; i < i_max; i++) {
    int returncode = pthread_join(tid[i], NULL);
    if (returncode != 0) {
      perror("join failed");
      exit(-1);
    }
  }

  unsigned long long tend = get_cycle_count();
  double end = get_seconds();

  resources rend;
  get_resources(&rend, THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power, num_bands);
  double avg_band_power = avg_of(band_power, num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++) {
    double band_low = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ", band, band_low, band_high,
           band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n\
User time        %lf seconds\n\
System time      %lf seconds\n\
Page faults      %ld\n\
Page swaps       %ld\n\
Blocks of I/O    %ld\n\
Signals caught   %ld\n\
Context switches %ld\n",
         rdiff.usertime, rdiff.systime, rdiff.pagefaults, rdiff.pageswaps,
         rdiff.ioblocks, rdiff.sigs, rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing "
         "overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one "
         "core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

int main(int argc, char *argv[]) {

  if (argc != 8) {
    usage();
    return -1;
  }

  char sig_type = toupper(argv[1][0]);
  char *sig_file = argv[2];
  double Fs = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  num_bands = atoi(argv[5]);
  num_threads = atoi(argv[6]);
  num_processors = atoi(argv[7]);
  tid = (pthread_t *)malloc(sizeof(pthread_t) * num_threads);

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n\
file:     %s\n\
Fs:       %lf Hz\n\
order:    %d\n\
bands:    %d\n",
         sig_type == 'T'
             ? "Text"
             : (sig_type == 'B'
                    ? "Binary"
                    : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file, Fs, filter_order, num_bands);

  printf("Load or map file\n");

  signal *sig;
  switch (sig_type) {
  case 'T':
    sig = load_text_format_signal(sig_file);
    break;

  case 'B':
    sig = load_binary_format_signal(sig_file);
    break;

  case 'M':
    sig = map_binary_format_signal(sig_file);
    break;

  default:
    printf("Unknown signal type\n");
    return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end,
           (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}
