// Header file for compute_coop_factor
#ifndef COMPUTE_COOP_FACTOR_H
#define COMPUTE_COOP_FACTOR_H
#include "problemDefines.h"
void compute_coop_factor(const int n_s,
                        const float parameter_reference_value,
                        const float B_coef, // like gamma_B or mu_B 
                        const float M_coef, // like gamma_B or mu_M
                        float *parameter_array_out, // like kB_plus and is length N_S*N_S this is the value that will be updated
                        const float final_exp // Like q = 1 or q = 1 -1 (etc) or r = 1
);
#endif // COMPUTE_COOP_FACTOR_H