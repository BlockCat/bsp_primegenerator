#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "bitarray.c"

#define MAX_PRIMES 1000
#define MAX_PRIMES_ROOT sqrt(MAX_PRIMES)

#define DEBUG

int main( int argc, char ** argv ) {	
	// Initiate bitarray
	int int_blocks;
	Bitarray bitarray = bitarray_create(&int_blocks, MAX_PRIMES);
	for (int i = 0; i < int_blocks; i++) {
		bitarray[i] = 0;
	}

	clock_t start = clock();
	int cost = 0;
	for (int j = 2; 2 * j <= MAX_PRIMES; j++) {
		bitarray_set(&bitarray, 2*j);
		cost++;
	}

	for (int i = 3; i < MAX_PRIMES_ROOT; i+=2) {
		if (bitarray_get(bitarray, i)) continue;

		int mmm = MAX_PRIMES / i;
		for (int j = i; j <= mmm; j+=2) {
			bitarray_set(&bitarray, i*j);
			cost++;
		}
	}

	clock_t time = clock() - start;

	int sum = 0;
	for (int i = 2; i <= MAX_PRIMES; i++) {
		if (bitarray_get(bitarray, i)) continue;
		sum++;
	}

	printf("Total primes: %d in seconds %ld,%ld\n", sum, (time / CLOCKS_PER_SEC),(time%CLOCKS_PER_SEC));
	printf("Total costs: %d\n", cost);	
	printf("%lu", sizeof(Bitarray));
    printf("\nPress enter to continue...");
    getchar();

    return EXIT_SUCCESS;
}


