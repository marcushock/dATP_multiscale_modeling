#ifndef UPDATE_RUS_H
#define UPDATE_RUS_H
#include "problemDefines.h"
__device__ void update_RUs(float lambda,
                float dt,
                float kCa_plus,
                float kCa_minus,
                float randNum[N_RU],
                float rand_dATP[N_RU],
                int   RU[N_RU],
                bool caRU[N_RU],
                float * kB_plus,
                float * kB_minus,
                float * k2_plus_dATP,
                float * k2_plus_ATP,
                float * k2_minus,
                float k3_plus_dATP,
                float k3_plus_ATP,
                float k3_minus,
                float *k4_plus_dATP,
                float *k4_plus_ATP,
                float *k4_minus,
                float percent_dATP,
                float k_force_dATP,
                float k_force_ATP, 
                float k_plus_SR_ref_dATP,
                float k_plus_SR_ref_ATP, 
                float k_minus_SR_ref,
                float f
               );
#endif // UPDATE_RUS_H
