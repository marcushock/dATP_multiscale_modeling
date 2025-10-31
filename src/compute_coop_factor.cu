#include "rates_trans_matrix.h"
#include "problemDefines.h"

#include <math.h>
// arr[x*N_S+y] == arr[x * row_len + y]
void compute_coop_factor(const int n_s,
    const float parameter_reference_value,
    const float B_coef, // like gamma_B or mu_B 
    const float M_coef, // like gamma_B or mu_M
    float *parameter_array_out, // like kB_plus and is length N_S*N_S this is the value that will be updated
    const float final_exp // Like q = 1 or q = 1 -1 (etc) or r = 1 
)
    {
        int B_states[] = {0, 2, 7};
        int M_states[] = {5, 6}; // Normally 4,5,6 but 4 is not used in this iteration. 
        // int C_states[3] = {1, 3, 8}; // not used in this function

        // Number of elements in the array
        int N_B_states = sizeof(B_states) / sizeof(B_states[0]);
        int N_M_states = sizeof(M_states) / sizeof(M_states[0]);

        for (int row = 0; row < n_s; row++) {
            for (int col = 0; col < n_s; col++){
                int B_count = 0;
                int M_count = 0;

                // Check if row is in B_states or M_states
                for (int i = 0; i < N_B_states; i++) {
                    if (row == B_states[i]) {
                        B_count++;
                    }
                    if (col == B_states[i]) {
                        B_count++;
                    }
                }
                // Check if col is in B_states
                for (int i = 0; i < N_M_states; i++) {
                    if (row == M_states[i]) {
                        M_count++;
                    }
                    else if (col == M_states[i]) {
                        M_count++;
                    }
                }
                float cooperative_value = pow(B_coef, -B_count) * pow(M_coef, M_count);
                parameter_array_out[row * n_s + col] = parameter_reference_value * pow(cooperative_value, final_exp);
            }
        } 
        return; 
    }

// B_states = [0, 2, 7]
// M_states = [4, 5, 6]
// C_states = [1, 3, 8]

// void compute_coop_factor(parameter, N_S, gamma_B, gamma_M):
// for row in N_S:
//     for col in N_S:
//         B_count = 0
//         M_count = 0
//         if row in B_states: 
//             B_count +=1 
//         elif row in M_states: 
//             M_count +=1
//         if col in B_states: 
//             B_count +=1
//         elif col in M_states: 
//             M_count +=1
//         parameter[row, col] = gamma_B ** B_count * gamma_M ** M_count
