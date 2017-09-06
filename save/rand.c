#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#define MAX 1000000000L

int main()
{
    int64_t i, count=0;
    double p;

    printf("RAND_MAX %d\n", RAND_MAX);

    p = 1e-8;

    for (i = 0; i < MAX; i++) {
        if ((double)random()/RAND_MAX <= p) {
            count++;
        }
    }
    printf("p = %lg\n", p);
    printf("MAX = %ld\n", MAX);
    printf("count = %ld\n", count);
    printf("result = %lg\n", (double)count/MAX);

    printf("compare = %lf\n", ((double)count/MAX) / p);
    return 0;
}
