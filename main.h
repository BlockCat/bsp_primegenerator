#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multicore-bsp/include/mcbsp.h"
#include "bitarray.c"

void spmd();
Bitarray preProcessingPrimes(int upper);
Bitarray crossOutPrimes(Bitarray bitarray, int rangeStart, int rangeEnd);
int countPrimes(Bitarray bitarray, int bound);


