// XXX comments

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "util_misc.h"

#define MAX_GEN   100000000   // 100 million
#define MAX_CACHE 1000000     // 1 million

#define MIN_VALUE  -500000
#define MAX_VALUE   500000

void * random_range_cache_init(int32_t min, int32_t max, int32_t max_cache);
void random_range_cache_destroy(void * handle);
int32_t random_range_cache(void * handle);

// -----------------  MAIN  ---------------------------------------


int32_t main(int argc, char **argv)
{
    void * rrc_handle;
    int64_t start_us, end_us, i;

    // init random_range_cache with MAX_CACHE entries
    rrc_handle = random_range_cache_init(MIN_VALUE, MAX_VALUE, MAX_CACHE);

    // time generating MAX_GEN random numbers by calling random_range
    start_us = microsec_timer();
    for (i = 0; i < MAX_GEN; i++) {
        random_range(MIN_VALUE, MAX_VALUE);
    }
    end_us = microsec_timer();
    INFO("random_range       : %lf ns\n", (end_us - start_us) * 1000.0 / MAX_GEN);

    // time generating MAX_GEN random numbers by calling random_range_cache
    start_us = microsec_timer();
    for (i = 0; i < MAX_GEN; i++) {
        random_range_cache(rrc_handle);
    }
    end_us = microsec_timer();
    INFO("random_range_cache : %lf ns\n", (end_us - start_us) * 1000.0 / MAX_GEN);

    // cleanup and exit
    random_range_cache_destroy(rrc_handle);
    return 0;
}

// -----------------  RANDOM RANGE CACHE  ------------------------------

typedef struct {
    int32_t max_cache;
    int32_t index;
    int32_t value[0];
} random_range_cache_t;

void * random_range_cache_init(int32_t min, int32_t max, int32_t max_cache)
{
    random_range_cache_t * rrc;
    int32_t i;

    // allocate
    rrc = malloc(sizeof(random_range_cache_t) + max_cache * sizeof(int32_t));
    assert(rrc);

    // init
    rrc->max_cache = max_cache;
    rrc->index = 0;
    for (i = 0; i < max_cache; i++) {
        rrc->value[i] = random_range(min,max);
    }

    // return handle
    return rrc;
}

void random_range_cache_destroy(void * handle)
{
    free(handle);
}

int32_t random_range_cache(void * handle)
{
    random_range_cache_t * rrc = handle;
    int32_t index;

    index = rrc->index;
    rrc->index = (index + 1) < rrc->max_cache ? (index + 1) : 0;
    return rrc->value[index];
}
