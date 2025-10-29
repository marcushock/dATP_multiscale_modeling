//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |          Function Name                |                         //
//                    |        rates_trans_matrix ()          |                         //
//                    |                                       |                         //
//--------------------------------------------------------------------------------------//
// Inputs   |
//---------
// parameters / reference values for transition rates
//--------
// Outputs|
//---------
// kB_plus,kB_minus,f,g:  Matrices-->" Tables " of coefficients that depend on Neighboring states (X,Y)
//--------------------------------------------------------------------------
// Notation
//----------
//
//            |-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//            |  B* |  C* |   B |   C |  M1 |  M2 |  M3 | B** | C** |  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     |  B*  |(0,0)|(0,1)|(0,2)|(0,3)|(0,4)|(0,5)|(0,6)|(0,7)|(0,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     |  C*  |(1,0)|(1,1)|(1,2)|(1,3)|(1,4)|(1,5)|(1,6)|(1,7)|(1,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     |   B  |(2,0)|(2,1)|(2,2)|(2,3)|(2,4)|(2,5)|(2,6)|(2,7)|(2,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     |   C  |(3,0)|(3,1)|(3,2)|(3,3)|(3,4)|(3,5)|(3,6)|(3,7)|(3,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     |  M1  |(4,0)|(4,1)|(4,2)|(4,3)|(4,4)|(4,5)|(4,6)|(4,7)|(4,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     |  M2  |(5,0)|(5,1)|(5,2)|(5,3)|(5,4)|(5,5)|(5,6)|(5,7)|(5,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     |  M3  |(6,0)|(6,1)|(6,2)|(6,3)|(6,4)|(6,5)|(6,6)|(6,7)|(6,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     | B**  |(7,0)|(7,1)|(7,2)|(7,3)|(7,4)|(7,5)|(7,6)|(7,7)|(7,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  
//     | C**  |(8,0)|(8,1)|(8,2)|(8,3)|(8,4)|(8,5)|(8,6)|(8,7)|(8,8)|  
//     |------|-----|-----|-----|-----|-----|-----|-----|-----|-----|  

//--------------------------------------------------------------------------------
    #include "rates_trans_matrix.h"
    #include "problemDefines.h"
    # include "compute_coop_factor.h"
    #include <math.h>
    // arr[x*N_S+y] == arr[x * row_len + y]
    void rates_trans_matrix(const int n_s,
                            const float kB_plus_ref,
                            const float kB_minus_ref,
                            const float k1_plus_ref_drug,
                            const float k1_plus_ref_baseline,
                            const float k1_minus_ref,
                            const float k4_plus_ref_drug,
                            const float k4_plus_ref_baseline,
                            const float k4_minus_ref,
                            const float gamma_B,
                            const float gamma_M,
                            const float mu_B,
                            const float mu_M,
                            const float r,
                            const float q,
                            float *kB_plus,
                            float *kB_minus,
                            float *k1_plus_drug,
                            float *k1_plus_baseline,
                            float *k1_minus,
                            float *k4_plus_drug,
                            float *k4_plus_baseline,
                            float *k4_minus
                        )

    {
    //--------------------------------------------------
    // Step 1: Build the kB_plus [ns*N_S+ns] matrix
    //--------------------------------------------------
    compute_coop_factor(n_s, kB_plus_ref, gamma_B, gamma_M, kB_plus, q);


    //---------------------------------------------------
    // Step 2: Build the kB_minus [ns*N_S+ns] matrix
    //---------------------------------------------------
    compute_coop_factor(n_s, kB_minus_ref, gamma_B, gamma_M, kB_minus, q-1);
    

    //--------------------------------------------------
    // Step 3: Build the k1_plus_reference [ns*N_S+ns] matrix
    //--------------------------------------------------
    compute_coop_factor(n_s, k1_plus_ref_baseline, mu_B, mu_M, k1_plus_baseline, r);
    

    //--------------------------------------------------
    // Step 3: Build the k1_plus_drug [ns*N_S+ns] matrix
    //--------------------------------------------------
    compute_coop_factor(n_s, k1_plus_ref_drug, mu_B, mu_M, k1_plus_drug, r);

    //--------------------------------------------------
    // Step 4: Build the k1_minus [ns*N_S+ns] matrix
    //--------------------------------------------------
    compute_coop_factor(n_s, k1_minus_ref, mu_B, mu_M, k1_minus, r-1);


    //---------------------------------------------------
    // Step 5: Build the k4_plus_reference [ns*N_S+ns] matrix
    //---------------------------------------------------
    compute_coop_factor(n_s, k4_plus_ref_baseline, mu_B, mu_M, k4_plus_baseline, r-1);
    

    //---------------------------------------------------
    // Step 5: Build the k4_plus_drug [ns*N_S+ns] matrix
    //---------------------------------------------------
    compute_coop_factor(n_s, k4_plus_ref_drug, mu_B, mu_M, k4_plus_drug, r-1);

    
    //--------------------------------------------------
    // Step 6: Build the k4_minus [ns*N_S+ns] matrix
    //--------------------------------------------------
    compute_coop_factor(n_s, k4_minus_ref, mu_B, mu_M, k4_minus, r);
    }
