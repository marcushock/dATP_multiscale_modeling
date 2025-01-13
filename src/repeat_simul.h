#ifndef REPEAT_SIMUL_H
#define REPEAT_SIMUL_H
#include "problemDefines.h"
__global__ void repeat_simul(float lambda,
                             const unsigned long randSeed,
                             float *k4_plus_drug,
                             float *k4_plus_baseline,
                             float *k4_minus,
                             float k3_plus_drug,
                             float k3_plus_baselineeline,
                             float k3_minus,
                             float *k2_plus_drug,
                             float *k2_plus_baseline,
                             float *k2_minus,
                             float *kB_plus,
                             float *kB_minus,
                             float kCa_plus_ref,
                             float kCa_minus_ref,
                             float percent_drug,
                             float k_force_drug,
                             float k_force_baseline,
                            float k_plus_SR_drug,
                            float k_plus_SR_baseline,
                            float k_minus_SR,
                             float * M3,
                             float * M1,
                             float * C,
                             float * B,
                             float * SR,
                             int cc, 
                             float protocol, 
                             float Calc_conc_exp 
                             );
#endif // REPEAT_SIMUL_H
