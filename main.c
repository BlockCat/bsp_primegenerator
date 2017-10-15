#include "main.h"

// If build has not defined MAX_PRIMES
#ifndef MAX_PRIMES
#define MAX_PRIMES 1000000000
#endif // !MAX_PRIMES
#define MAX_PRIMES_ROOT sqrt(MAX_PRIMES)

// If build has not defined CORES
#ifndef CORES
#define CORES 4
#endif

int main( int argc, char ** argv ) {
    printf("cores: %d\n", bsp_nprocs());

    bsp_init( &spmd, argc, argv );	
    spmd();

    return EXIT_SUCCESS;
}

// Calculate the range of {processor} given amount of {cores}
void calculateRange(int pid, int cores, int * rangeStart, int * rangeEnd) {
	(*rangeStart) = (MAX_PRIMES_ROOT + (pid * (MAX_PRIMES - MAX_PRIMES_ROOT) / cores));
	(*rangeEnd) = (MAX_PRIMES_ROOT + ceil((pid + 1) * (MAX_PRIMES - MAX_PRIMES_ROOT) / cores));
}

void spmd() {	
	int cores = CORES;
	bsp_begin(cores);		
	double start = bsp_time();
	// Get the first primes up to MAX_PRIMES_ROOT
	Bitarray preprocess = preProcessingPrimes(MAX_PRIMES_ROOT);	

	int pid = bsp_pid();

	int rangeStart;
	int rangeEnd;

	calculateRange(pid, cores, &rangeStart, &rangeEnd);

	printf("(%d) Range: [%d,%d]\n", pid, rangeStart, rangeEnd);
	
	Bitarray primesArray = crossOutPrimes(preprocess, rangeStart, rangeEnd);
	printf("Time crossed out: %f\n", bsp_time() - start);		

	// Start merging
	int blockSize = bitarray_blocks(MAX_PRIMES);	
	//Bitarray vector = (Bitarray)malloc(blockSize * cores * (int)sizeof(BitBlock));
	//bsp_push_reg(vector, blockSize * cores * (int)sizeof(BitBlock));

	Bitarray vector = (Bitarray)malloc(blockSize * (int)sizeof(BitBlock));
	bsp_push_reg(vector, blockSize * (int)sizeof(BitBlock));

	bsp_sync();	

	// Send data to other processors		
	// Calculate range of blocks
	int leftBlock = rangeStart / BLOCK_SIZE;
	int rightBlock = ceil(rangeEnd / (double)BLOCK_SIZE);
	
	//Get the relevant array
	int* partial = &primesArray[leftBlock];

	for (int i = 0; i < cores; i++) {		
		bsp_put(i, partial, vector, leftBlock * sizeof(BitBlock), (rightBlock - leftBlock) * sizeof(BitBlock));				
	}
	bsp_sync();
	
	//Merge preprocessed primes into the result
	for (int i = 0; i < bitarray_blocks(MAX_PRIMES_ROOT); i++) {						
		vector[i] |= preprocess[i];
	}	

	printf("Total time: %f\n", bsp_time() - start);
	int amount = countPrimes(vector, MAX_PRIMES);
	printf("Amount of primes: %d\n", amount);
	
	// Print out twin primes
	if (pid == 0) {
		for (int i = 2; i < MAX_PRIMES-2; i++) {
			if (bitarray_get(vector, i) == 0 && bitarray_get(vector, i + 2) == 0 ) {
				printf("%d:%d, ", i, i+2);
			}
		}
	}

	// Print out goldbach primes
	if (pid == 0) {
		struct GoldBach* bacharray = createGoldBachPairs(vector, MAX_PRIMES);
		printGoldBachArray(bacharray, MAX_PRIMES);
	}

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
		for (int j = prime*prime; j <= upper; j += prime+prime) {						
			bitarray_set(&bitarray, j);
		}
	}	
	return bitarray;
}

Bitarray crossOutPrimes(Bitarray bitarray, int rangeStart, int rangeEnd) {
	//Make the ranges dividable by block_size for faster merging.
	rangeStart = (rangeStart / BLOCK_SIZE)*BLOCK_SIZE - 1;
	rangeEnd = rangeEnd /BLOCK_SIZE*BLOCK_SIZE + 1;
	rangeStart = (rangeStart > 0) ? rangeStart : 0;
	rangeStart = (rangeStart < MAX_PRIMES_ROOT) ? MAX_PRIMES_ROOT : rangeStart;
	rangeEnd = (rangeEnd > MAX_PRIMES) ? rangeEnd : MAX_PRIMES;

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

struct GoldBach* createGoldBachPairs(Bitarray primes, int upperBound) {
	
	struct GoldBach* bacharray = (struct GoldBach*)malloc(sizeof(struct GoldBach) * upperBound / 2);
	bacharray[4 / 2] = (struct GoldBach){ 2, 2 };

	for (int i = 2; i < upperBound / 2; i++) {
		if (bitarray_get(primes, i)) continue;
		for (int j = i; i + j < upperBound; j++) {
			if (bitarray_get(primes, j)) continue;

			//i and j are both primes,
			bacharray[(i + j) / 2] = (struct GoldBach) { i, j };
		}
	}
	return bacharray;
}

void printGoldBachArray(struct GoldBach* bacharray, int upperbound) {
	for (int i = 2; i < upperbound / 2; i++) {
		struct GoldBach gb = bacharray[i];
		printf("%d: [%d,%d]\n", i * 2, gb.prime1, gb.prime2);
	}
}