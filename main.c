#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multicore-bsp/include/mcbsp.h"
#include "main.h"

#ifndef MAX_PRIMES
#define MAX_PRIMES 1000000000
#endif // !MAX_PRIMES


#define MAX_PRIMES_ROOT sqrt(MAX_PRIMES)


int main( int argc, char ** argv ) {
    printf("cores: %d\n", bsp_nprocs());

    bsp_init( &spmd, argc, argv );	
    spmd();

    printf("\nPress enter to continue...");
    getchar();

    return EXIT_SUCCESS;
}

void calculateRange(int pid, int cores, int * rangeStart, int * rangeEnd) {
	(*rangeStart) = (MAX_PRIMES_ROOT + (pid * (MAX_PRIMES - MAX_PRIMES_ROOT) / cores));
	(*rangeEnd) = (MAX_PRIMES_ROOT + ceil((pid + 1) * (MAX_PRIMES - MAX_PRIMES_ROOT) / cores));
}

void spmd() {
	int cores = bsp_nprocs();	
	cores = 4;
	bsp_begin(cores);		
	// Get the first primes up to MAX_PRIMES_ROOT
	Bitarray preprocess = preProcessingPrimes(MAX_PRIMES_ROOT);	
	int pid = bsp_pid();

	int rangeStart;// = MAX_PRIMES_ROOT + (pid * (MAX_PRIMES - MAX_PRIMES_ROOT) / cores);
	int rangeEnd;// = MAX_PRIMES_ROOT + ceil((pid + 1) * (MAX_PRIMES - MAX_PRIMES_ROOT) / cores);

	calculateRange(pid, cores, &rangeStart, &rangeEnd);

	printf("(%d) Range: [%d,%d]\n", pid, rangeStart, rangeEnd);
	double start = bsp_time();
	Bitarray primesArray = crossOutPrimes(preprocess, rangeStart, rangeEnd);
	printf("Time crossed out: %f\n", bsp_time() - start);
	// TODO: SYNC
	bsp_sync();	

	// Start merging
	int blockSize = bitarray_blocks(MAX_PRIMES);	
	Bitarray vector = (Bitarray)malloc(blockSize * cores * (int)sizeof(BitBlock));
	
	bsp_push_reg(vector, blockSize * cores * (int)sizeof(BitBlock));
	bsp_sync();	
	for (int i = 0; i < cores; i++) {
		bsp_put(i, primesArray, vector, pid * blockSize * sizeof(BitBlock), blockSize * sizeof(BitBlock));
	}
	bsp_sync();
	
	//Merge preprocessed primes
	for (int i = 0; i < bitarray_blocks(MAX_PRIMES_ROOT); i++) {		
		primesArray[i] |= preprocess[i];
	}	

	for (int core = 0; core < cores; core++) {
		int rl, rr;
		calculateRange(core, cores, &rl, &rr);
		int leftBlock = rl / BLOCK_SIZE;
		int rightBlock = rr / BLOCK_SIZE;

		for (int i = leftBlock; i <= rightBlock && i < blockSize; i++) {
			primesArray[i] |= vector[i + blockSize * core];
		}

	}

	int amount = countPrimes(primesArray, MAX_PRIMES);
	printf("Amount of primes: %d\n", amount);

	/*if (pid == 0) {
		for (int i = 2; i < MAX_PRIMES; i++) {
			if (bitarray_get(primesArray, i) == 0) {
				printf("%d, ", i);
			}
		}
	}*/

	// Clean up memory	
	free(primesArray);
	free(preprocess);
	bsp_end();
}

///
/// Calculate the first primes till MAX_PRIMES_ROOT
///
Bitarray preProcessingPrimes(int upper) {
	int blocks;
	Bitarray bitarray = bitarray_create(&blocks, upper);

	// Cross out	
	for (int i = 2*2; i <= upper; i+=2) {
		bitarray_set(&bitarray, i);
	}

	// Cross out other primes, skip even numbers.
	int sqrtUpper = ceil(sqrt(upper));
	for (int prime = 3; prime < sqrtUpper; prime += 2) {			
		for (int j = 2*prime; j <= upper; j += prime) {						
			bitarray_set(&bitarray, j);
		}
	}	
	return bitarray;
}

Bitarray crossOutPrimes(Bitarray bitarray, int rangeStart, int rangeEnd) {
	int range = rangeEnd - rangeStart;
	int blocks;
	Bitarray blockPrimes = bitarray_create(&blocks, MAX_PRIMES);
	
	int start_2 = ((int)(rangeStart / 2) + 1) * 2; // Find the first multiplication of 2 after rangeStart, example: ceil(1001 / 2) * 2 = 1002	
	for (int crossOut = start_2; crossOut <= rangeEnd; crossOut+=2) {
		bitarray_set(&blockPrimes, crossOut);
	}

	for (int prime = 3; prime <= MAX_PRIMES_ROOT; prime += 2) {		
		if (bitarray_get(bitarray, prime) == 0) {			
			int start_prime = ((int)(rangeStart / prime) + 1) * prime;			
			for (int crossOut = start_prime; crossOut <= rangeEnd; crossOut += prime) {				
				bitarray_set(&blockPrimes, crossOut);
			}
		}
	}

	return blockPrimes;
}

int countPrimes(Bitarray bitarray, int bound) {
	int sum = 0;
	for (int i = 2; i <= bound; i++) {
		if (bitarray_get(bitarray, i) == 0) {
			sum += 1;
		}
	}
	return sum;
}