#ifndef GENRAND_H
#define GENRAND_H
#include <curand_kernel.h>
__device__ void genrand(float * numbers, const size_t numbersLen, curandState * state);
#endif // GENRAND_H
