gcc main.c -o main.out multicore-bsp/lib/libmcbsp1.2.0.a -pthread -lrt -lm -DMAX_PRIMES=$1 -DCORES=$2
./main.out
