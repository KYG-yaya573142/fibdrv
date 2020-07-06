#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bn.h"

#define ITH 1000

int main(int argc, char const *argv[])
{
    bn *test = bn_alloc(1);
    for (int i = 0; i < ITH + 1; i++) {
        bn_fib_fdoubling(test, i);
        char *p = bn_to_string(test);
        printf("%d,%s\n", i, p);
        free(p);
    }
    bn_free(test);
    return 0;
}