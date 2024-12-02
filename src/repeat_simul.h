#ifndef REPEAT_SIMUL_H
#define REPEAT_SIMUL_H
#include "problemDefines.h"
__global__ void repeat_simul(float lambda,
                             const unsigned long randSeed,
                             float *k4_plus_dATP,
                             float *k4_plus_ATP,
                             float *k4_minus,
                             float k3_plus_dATP,
                             float k3_plus_ATP,
                             float k3_minus,
                             float *k2_plus_dATP,
                             float *k2_plus_ATP,
                             float *k2_minus,
                             float *kB_plus,
                             float *kB_minus,
                             float kCa_plus_ref,
                             float kCa_minus_ref,
                             float percent_dATP,
                             float k_force_dATP,
                             float k_force_ATP,
			     float k_plus_SR_ref_dATP,
			     float k_plus_SR_ref_ATP,
			     float k_minus_SR_ref,
                             float * Force,
                             float * Mc,
                             float * C,
                             float * B,
                             float * SR,
                             int cc, 
                             int protocol, 
                             float Calc_conc_exp 
                             );
#endif // REPEAT_SIMUL_H
