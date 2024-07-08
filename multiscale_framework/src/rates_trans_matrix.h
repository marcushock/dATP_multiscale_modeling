#ifndef RATES_TRANS_MATRIX_H
#define RATES_TRANS_MATRIX_H
#include "problemDefines.h"
void rates_trans_matrix(const int n_s,
                        const float kB_plus_ref,
                        const float kB_minus_ref,
                        const float k2_plus_ref_dATP,
                        const float k2_plus_ref_ATP,
                        const float k2_minus_ref,
                        const float k4_plus_ref_dATP,
                        const float k4_plus_ref_ATP,
                        const float k4_minus_ref,
                        const float gamma_B,
                        const float gamma_M,
                        const float mu_B,
                        const float mu_M,
                        const float r,
                        const float q,
                        float *kB_plus,
                        float *kB_minus,
                        float *k2_plus_dATP,
                        float *k2_plus_ATP,
                        float *k2_minus,
                        float *k4_plus_dATP,
                        float *k4_plus_ATP,
                        float *k4_minus
                       );
#endif // RATES_TRANS_MATRIX_H
