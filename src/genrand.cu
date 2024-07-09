#include "genrand.h"
#include <cuda.h>
#include <curand_kernel.h>
__device__ void genrand(float * numbers, const size_t numbersLen, curandState * state)
{
    for(int i = 0; i < numbersLen; ++i)
    {
        numbers[i] = curand_uniform(state);
    }
}
