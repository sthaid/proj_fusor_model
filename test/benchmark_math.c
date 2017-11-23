#if 0  
// XXX description

// RESULTS 
$ ./benchmark  1
      int32_test  threads=1  million-ops/s=346.253538
      int64_test  threads=1  million-ops/s=116.898357
     double_test  threads=1  million-ops/s=249.112458
$ ./benchmark  2
      int32_test  threads=2  million-ops/s=687.078873
      int64_test  threads=2  million-ops/s=238.941244
     double_test  threads=2  million-ops/s=498.122757
$ ./benchmark  3
      int32_test  threads=3  million-ops/s=766.516162
      int64_test  threads=3  million-ops/s=274.279253
     double_test  threads=3  million-ops/s=497.769643
$ ./benchmark  4
      int32_test  threads=4  million-ops/s=840.779316
      int64_test  threads=4  million-ops/s=309.315182
     double_test  threads=4  million-ops/s=497.681311
#endif

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#include "util_misc.h"

// variables

volatile bool terminate;
int32_t terminated_count;
pthread_barrier_t barrier;

// prototypes

void run_benchmark(char * name, void * (*proc)(void *), int32_t num_threads);
void use_the_result(void * volatile x);
void * int32_test(void * cx);
void * int64_test(void * cx);
void * double_test(void * cx);

// -----------------  MAIN  --------------------------------------------

int main(int argc, char ** argv)
{
    int32_t num_threads;

    if (argc != 2 || sscanf(argv[1], "%d", &num_threads) != 1) {
        printf("usage: benchmark <num_threads>\n");
        return -1;
    }

    run_benchmark("int32_test", int32_test, num_threads);
    run_benchmark("int64_test", int64_test, num_threads);
    run_benchmark("double_test", double_test, num_threads);
}

void run_benchmark(char * name, void * (*proc)(void *), int32_t num_threads)
{
    int64_t ops_count[1000], total_ops_count, start_us, end_us;
    int32_t i;
    pthread_t thread_id;

    // init, and create threads
    terminate = false;
    terminated_count = 0;
    memset(ops_count, 0, sizeof(ops_count));
    total_ops_count = 0;
    pthread_barrier_init(&barrier, NULL, num_threads+1);
    for (i = 0; i < num_threads; i++) {
        pthread_create(&thread_id, NULL, proc, ops_count+i);
    }

    // let the threads run for 1 second
    pthread_barrier_wait(&barrier);
    start_us = microsec_timer();
    sleep(1);
    end_us = microsec_timer();
    terminate = true;

    // wait for all threads to indicate they have terminated
    while (terminated_count != num_threads) {
        usleep(1000);
    }

    // print results
    for (i = 0; i < num_threads; i++) {
        total_ops_count += ops_count[i];
    }
    printf("%16s  threads=%d  million-ops/s=%f\n",
        name, num_threads,
        (double)total_ops_count * (1000000.0/(end_us-start_us)) / 1000000.0);

    // cleanup
    pthread_barrier_destroy(&barrier);
}

void use_the_result(void * volatile x) 
{
    // call to this prevents the compiler from optimizing 
    // out the benchmark, which it will do because it determines
    // that the result of the computation is not used
}

// -----------------  BENCHMARKS  --------------------------------------

void * int32_test(void * cx)
{
    int64_t * ops_count = cx;
    int32_t sum, a, b, c;

    pthread_barrier_wait(&barrier);

    sum = a = b = 0; c = 1;
    while (true) {
        sum = sum + a * b / c;
        a += 1;
        b += 2;
        c += 1;

        (*ops_count)++;
        if (terminate) {
            __sync_fetch_and_add(&terminated_count,1);
            use_the_result(&sum);
            break;
        }
    }

    return NULL;
}

void * int64_test(void * cx)
{
    int64_t * ops_count = cx;
    int64_t sum, a, b, c;

    pthread_barrier_wait(&barrier);

    sum = a = b = 0; c = 1;
    while (true) {
        sum = sum + a * b / c;
        a += 1;
        b += 2;
        c += 1;

        (*ops_count)++;
        if (terminate) {
            __sync_fetch_and_add(&terminated_count,1);
            use_the_result(&sum);
            break;
        }
    }

    return NULL;
}

void * double_test(void * cx)
{
    int64_t * ops_count = cx;
    double sum, a, b, c;

    pthread_barrier_wait(&barrier);

    sum = a = b = 0; c = 1;
    while (true) {
        sum = sum + a * b / c;
        a += 1;
        b += 2;
        c += 1;

        (*ops_count)++;
        if (terminate) {
            __sync_fetch_and_add(&terminated_count,1);
            use_the_result(&sum);
            break;
        }
    }

    return NULL;
}
