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
//          |-----|-----|-----|-----|-----|-----|-----|  
//          |  B* |  C* |   B |   C |  M1 |  M2 |  M3 |  
//     |----|-----|-----|-----|-----|-----|-----|-----|  
//     | B* |(0,0)|(0,1)|(0,2)|(0,3)|(0,4)|(0,5)|(0,6)|  
//     |----|-----|-----|-----|-----|-----|-----|-----|  
//     | C* |(1,0)|(1,1)|(1,2)|(1,3)|(1,4)|(1,5)|(1,6)|  
//     |----|-----|-----|-----|-----|-----|-----|-----|  
//     |  B |(2,0)|(2,1)|(2,2)|(2,3)|(2,4)|(2,5)|(2,6)|  
//     |----|-----|-----|-----|-----|-----|-----|-----|  
//     |  C |(3,0)|(3,1)|(3,2)|(3,3)|(3,4)|(3,5)|(3,6)|  
//     |----|-----|-----|-----|-----|-----|-----|-----|  
//     | M1 |(4,0)|(4,1)|(4,2)|(4,3)|(4,4)|(4,5)|(4,6)|  
//     |----|-----|-----|-----|-----|-----|-----|-----|  
//     | M2 |(5,0)|(5,1)|(5,2)|(5,3)|(5,4)|(5,5)|(5,6)|  
//     |----|-----|-----|-----|-----|-----|-----|-----|  
//     | M3 |(6,0)|(6,1)|(6,2)|(6,3)|(6,4)|(6,5)|(6,6)|  
//     |----|-----|-----|-----|-----|-----|-----|-----|
//--------------------------------------------------------------------------------
    #include "rates_trans_matrix.h"
    #include "problemDefines.h"

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
    kB_plus[0*N_S+0] = kB_plus_ref*pow(gamma_B,-2*q);
    kB_plus[0*N_S+1] = kB_plus_ref*pow(gamma_B,-q);
    kB_plus[0*N_S+2] = kB_plus_ref*pow(gamma_B,-2*q);
    kB_plus[0*N_S+3] = kB_plus_ref*pow(gamma_B,-q);
    kB_plus[0*N_S+4] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[0*N_S+5] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[0*N_S+6] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    //-------
    kB_plus[1*N_S+0] = kB_plus_ref*pow(gamma_B,-q);
    kB_plus[1*N_S+1] = kB_plus_ref;
    kB_plus[1*N_S+2] = kB_plus_ref*pow(gamma_B,-q);
    kB_plus[1*N_S+3] = kB_plus_ref;
    kB_plus[1*N_S+4] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[1*N_S+5] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[1*N_S+6] = kB_plus_ref*pow(gamma_M,q);
    //-------
    kB_plus[2*N_S+0] = kB_plus_ref*pow(gamma_B,-2*q);
    kB_plus[2*N_S+1] = kB_plus_ref*pow(gamma_B,-q);
    kB_plus[2*N_S+2] = kB_plus_ref*pow(gamma_B,-2*q);
    kB_plus[2*N_S+3] = kB_plus_ref*pow(gamma_B,-q);
    kB_plus[2*N_S+4] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[2*N_S+5] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[2*N_S+6] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    //-------
    kB_plus[3*N_S+0] = kB_plus_ref*pow(gamma_B,-q);
    kB_plus[3*N_S+1] = kB_plus_ref;
    kB_plus[3*N_S+2] = kB_plus_ref*pow(gamma_B,-q);
    kB_plus[3*N_S+3] = kB_plus_ref;
    kB_plus[3*N_S+4] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[3*N_S+5] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[3*N_S+6] = kB_plus_ref*pow(gamma_M,q);
    //-------
    kB_plus[4*N_S+0] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[4*N_S+1] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[4*N_S+2] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[4*N_S+3] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[4*N_S+4] = kB_plus_ref*pow(gamma_M,2*q);
    kB_plus[4*N_S+5] = kB_plus_ref*pow(gamma_M,2*q);
    kB_plus[4*N_S+6] = kB_plus_ref*pow(gamma_M,2*q);
    //-------
    kB_plus[5*N_S+0] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[5*N_S+1] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[5*N_S+2] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[5*N_S+3] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[5*N_S+4] = kB_plus_ref*pow(gamma_M,2*q);
    kB_plus[5*N_S+5] = kB_plus_ref*pow(gamma_M,2*q);
    kB_plus[5*N_S+6] = kB_plus_ref*pow(gamma_M,2*q);
    //-------
    kB_plus[6*N_S+0] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[6*N_S+1] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[6*N_S+2] = kB_plus_ref*pow((gamma_M/gamma_B),q);
    kB_plus[6*N_S+3] = kB_plus_ref*pow(gamma_M,q);
    kB_plus[6*N_S+4] = kB_plus_ref*pow(gamma_M,2*q);
    kB_plus[6*N_S+5] = kB_plus_ref*pow(gamma_M,2*q);
    kB_plus[6*N_S+6] = kB_plus_ref*pow(gamma_M,2*q);
    //---------------------------------------------------

    //---------------------------------------------------
    // Step 2: Build the kB_minus [ns*N_S+ns] matrix
    //---------------------------------------------------
    kB_minus[0*N_S+0] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
    kB_minus[0*N_S+1] = kB_minus_ref*pow(1/gamma_B,q-1);
    kB_minus[0*N_S+2] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
    kB_minus[0*N_S+3] = kB_minus_ref*pow(1/gamma_B,q-1);
    kB_minus[0*N_S+4] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[0*N_S+5] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[0*N_S+6] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    //-------
    kB_minus[1*N_S+0] = kB_minus_ref*pow(1/gamma_B,q-1);
    kB_minus[1*N_S+1] = kB_minus_ref;
    kB_minus[1*N_S+2] = kB_minus_ref*pow(1/gamma_B,q-1);
    kB_minus[1*N_S+3] = kB_minus_ref;
    kB_minus[1*N_S+4] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[1*N_S+5] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[1*N_S+6] = kB_minus_ref*pow(gamma_M,q-1);
    //-------
    kB_minus[2*N_S+0] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
    kB_minus[2*N_S+1] = kB_minus_ref*pow(1/gamma_B,q-1);
    kB_minus[2*N_S+2] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
    kB_minus[2*N_S+3] = kB_minus_ref*pow(1/gamma_B,q-1);
    kB_minus[2*N_S+4] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[2*N_S+5] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[2*N_S+6] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    //-------
    kB_minus[3*N_S+0] = kB_minus_ref*pow(1/gamma_B,q-1);
    kB_minus[3*N_S+1] = kB_minus_ref;
    kB_minus[3*N_S+2] = kB_minus_ref*pow(1/gamma_B,q-1);
    kB_minus[3*N_S+3] = kB_minus_ref;
    kB_minus[3*N_S+4] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[3*N_S+5] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[3*N_S+6] = kB_minus_ref*pow(gamma_M,q-1);
    //-------
    kB_minus[4*N_S+0] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[4*N_S+1] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[4*N_S+2] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[4*N_S+3] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[4*N_S+4] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
    kB_minus[4*N_S+5] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
    kB_minus[4*N_S+6] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
    //-------
    kB_minus[5*N_S+0] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[5*N_S+1] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[5*N_S+2] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[5*N_S+3] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[5*N_S+4] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
    kB_minus[5*N_S+6] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
    //-------
    kB_minus[6*N_S+0] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[6*N_S+1] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[6*N_S+2] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
    kB_minus[6*N_S+3] = kB_minus_ref*pow(gamma_M,q-1);
    kB_minus[6*N_S+4] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
    kB_minus[6*N_S+6] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
    //-----------------------------------

    //--------------------------------------------------
    // Step 3: Build the k1_plus [ns*N_S+ns] matrix
    //--------------------------------------------------
    k1_plus_baseline[0*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[0*N_S+1] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[0*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[0*N_S+3] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[0*N_S+4] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[0*N_S+5] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[0*N_S+6] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    //-------
    k1_plus_baseline[1*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[1*N_S+1] = k1_plus_ref_baseline;
    k1_plus_baseline[1*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[1*N_S+3] = k1_plus_ref_baseline;
    k1_plus_baseline[1*N_S+4] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[1*N_S+5] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[1*N_S+6] = k1_plus_ref_baseline*pow(mu_M,r);
    //-------
    k1_plus_baseline[2*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[2*N_S+1] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[2*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[2*N_S+3] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[2*N_S+4] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[2*N_S+5] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[2*N_S+6] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    //-------
    k1_plus_baseline[3*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[3*N_S+1] = k1_plus_ref_baseline;
    k1_plus_baseline[3*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[3*N_S+3] = k1_plus_ref_baseline;
    k1_plus_baseline[3*N_S+4] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[3*N_S+5] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[3*N_S+6] = k1_plus_ref_baseline*pow(mu_M,r);
    //-------
    k1_plus_baseline[4*N_S+0] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[4*N_S+1] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[4*N_S+2] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[4*N_S+3] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[4*N_S+4] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[4*N_S+5] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[4*N_S+6] = k1_plus_ref_baseline*pow(mu_M,2*r);
    //-------
    k1_plus_baseline[5*N_S+0] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[5*N_S+1] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[5*N_S+2] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[5*N_S+3] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[5*N_S+4] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[5*N_S+5] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[5*N_S+6] = k1_plus_ref_baseline*pow(mu_M,2*r);
    //-------
    k1_plus_baseline[6*N_S+0] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[6*N_S+1] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[6*N_S+2] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[6*N_S+3] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[6*N_S+4] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[6*N_S+5] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[6*N_S+6] = k1_plus_ref_baseline*pow(mu_M,2*r);
    //-------



    k1_plus_drug[0*N_S+0] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[0*N_S+1] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[0*N_S+2] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[0*N_S+3] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[0*N_S+4] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[0*N_S+5] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[0*N_S+6] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    //-------
    k1_plus_drug[1*N_S+0] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[1*N_S+1] = k1_plus_ref_drug;
    k1_plus_drug[1*N_S+2] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[1*N_S+3] = k1_plus_ref_drug;
    k1_plus_drug[1*N_S+4] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[1*N_S+5] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[1*N_S+6] = k1_plus_ref_drug*pow(mu_M,r);
    //-------
    k1_plus_drug[2*N_S+0] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[2*N_S+1] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[2*N_S+2] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[2*N_S+3] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[2*N_S+4] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[2*N_S+5] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[2*N_S+6] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    //-------
    k1_plus_drug[3*N_S+0] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[3*N_S+1] = k1_plus_ref_drug;
    k1_plus_drug[3*N_S+2] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[3*N_S+3] = k1_plus_ref_drug;
    k1_plus_drug[3*N_S+4] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[3*N_S+5] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[3*N_S+6] = k1_plus_ref_drug*pow(mu_M,r);
    //-------
    k1_plus_drug[4*N_S+0] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[4*N_S+1] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[4*N_S+2] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[4*N_S+3] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[4*N_S+4] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[4*N_S+5] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[4*N_S+6] = k1_plus_ref_drug*pow(mu_M,2*r);
    //-------
    k1_plus_drug[5*N_S+0] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[5*N_S+1] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[5*N_S+2] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[5*N_S+3] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[5*N_S+4] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[5*N_S+5] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[5*N_S+6] = k1_plus_ref_drug*pow(mu_M,2*r);
    //-------
    k1_plus_drug[6*N_S+0] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[6*N_S+1] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[6*N_S+2] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[6*N_S+3] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[6*N_S+4] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[6*N_S+5] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[6*N_S+6] = k1_plus_ref_drug*pow(mu_M,2*r);
    //---------------------------------------------------

    //---------------------------------------------------
    // Step 4: Build the k1_minus [ns*N_S+ns] matrix
    //---------------------------------------------------
    k1_minus[0*N_S+0] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[0*N_S+1] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[0*N_S+2] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[0*N_S+3] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[0*N_S+4] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[0*N_S+5] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[0*N_S+6] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    //-------
    k1_minus[1*N_S+0] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[1*N_S+1] = k1_minus_ref;
    k1_minus[1*N_S+2] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[1*N_S+3] = k1_minus_ref;
    k1_minus[1*N_S+4] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[1*N_S+5] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[1*N_S+6] = k1_minus_ref*pow(mu_M,r-1);
    //-------
    k1_minus[2*N_S+0] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[2*N_S+1] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[2*N_S+2] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[2*N_S+3] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[2*N_S+4] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[2*N_S+5] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[2*N_S+6] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    //-------
    k1_minus[3*N_S+0] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[3*N_S+1] = k1_minus_ref;
    k1_minus[3*N_S+2] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[3*N_S+3] = k1_minus_ref;
    k1_minus[3*N_S+4] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[3*N_S+5] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[3*N_S+6] = k1_minus_ref*pow(mu_M,r-1);
    //-------
    k1_minus[4*N_S+0] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[4*N_S+1] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[4*N_S+2] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[4*N_S+3] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[4*N_S+4] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[4*N_S+5] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[4*N_S+6] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    //-------
    k1_minus[5*N_S+0] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[5*N_S+1] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[5*N_S+2] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[5*N_S+3] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[5*N_S+4] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[5*N_S+5] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[5*N_S+6] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    //-------
    k1_minus[6*N_S+0] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[6*N_S+1] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[6*N_S+2] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[6*N_S+3] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[6*N_S+4] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[6*N_S+5] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[6*N_S+6] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    //-----------------------------------

    //---------------------------------------------------
    // Step 5: Build the k4_plus [ns*N_S+ns] matrix
    //---------------------------------------------------
    k4_plus_baseline[0*N_S+0] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[0*N_S+1] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[0*N_S+2] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[0*N_S+3] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[0*N_S+4] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[0*N_S+5] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[0*N_S+6] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    //-------
    k4_plus_baseline[1*N_S+0] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[1*N_S+1] = k4_plus_ref_baseline;
    k4_plus_baseline[1*N_S+2] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[1*N_S+3] = k4_plus_ref_baseline;
    k4_plus_baseline[1*N_S+4] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[1*N_S+5] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[1*N_S+6] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    //-------
    k4_plus_baseline[2*N_S+0] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[2*N_S+1] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[2*N_S+2] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[2*N_S+3] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[2*N_S+4] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[2*N_S+5] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[2*N_S+6] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    //-------
    k4_plus_baseline[3*N_S+0] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[3*N_S+1] = k4_plus_ref_baseline;
    k4_plus_baseline[3*N_S+2] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[3*N_S+3] = k4_plus_ref_baseline;
    k4_plus_baseline[3*N_S+4] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[3*N_S+5] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[3*N_S+6] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    //-------
    k4_plus_baseline[4*N_S+0] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[4*N_S+1] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[4*N_S+2] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[4*N_S+3] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[4*N_S+4] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[4*N_S+5] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[4*N_S+6] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    //-------
    k4_plus_baseline[5*N_S+0] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[5*N_S+1] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[5*N_S+2] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[5*N_S+3] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[5*N_S+4] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[5*N_S+5] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[5*N_S+6] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    //-------
    k4_plus_baseline[6*N_S+0] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[6*N_S+1] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[6*N_S+2] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[6*N_S+3] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[6*N_S+4] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[6*N_S+5] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[6*N_S+6] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    //-------




    k4_plus_drug[0*N_S+0] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[0*N_S+1] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[0*N_S+2] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[0*N_S+3] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[0*N_S+4] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[0*N_S+5] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[0*N_S+6] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    //-------
    k4_plus_drug[1*N_S+0] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[1*N_S+1] = k4_plus_ref_drug;
    k4_plus_drug[1*N_S+2] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[1*N_S+3] = k4_plus_ref_drug;
    k4_plus_drug[1*N_S+4] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[1*N_S+5] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[1*N_S+6] = k4_plus_ref_drug*pow(mu_M,(r-1));
    //-------
    k4_plus_drug[2*N_S+0] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[2*N_S+1] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[2*N_S+2] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[2*N_S+3] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[2*N_S+4] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[2*N_S+5] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[2*N_S+6] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    //-------
    k4_plus_drug[3*N_S+0] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[3*N_S+1] = k4_plus_ref_drug;
    k4_plus_drug[3*N_S+2] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[3*N_S+3] = k4_plus_ref_drug;
    k4_plus_drug[3*N_S+4] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[3*N_S+5] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[3*N_S+6] = k4_plus_ref_drug*pow(mu_M,(r-1));
    //-------
    k4_plus_drug[4*N_S+0] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[4*N_S+1] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[4*N_S+2] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[4*N_S+3] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[4*N_S+4] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[4*N_S+5] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[4*N_S+6] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    //-------
    k4_plus_drug[5*N_S+0] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[5*N_S+1] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[5*N_S+2] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[5*N_S+3] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[5*N_S+4] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[5*N_S+5] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[5*N_S+6] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    //-------
    k4_plus_drug[6*N_S+0] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[6*N_S+1] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[6*N_S+2] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[6*N_S+3] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[6*N_S+4] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[6*N_S+5] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[6*N_S+6] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    //-----------------------------------

    //--------------------------------------------------
    // Step 6: Build the k4_minus [ns*N_S+ns] matrix
    //--------------------------------------------------
    k4_minus[0*N_S+0] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[0*N_S+1] = k4_minus_ref*pow(mu_B,-1*r);
    k4_minus[0*N_S+2] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[0*N_S+3] = k4_minus_ref*pow(mu_B,-1*r);
    k4_minus[0*N_S+4] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[0*N_S+5] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[0*N_S+6] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    //-------
    k4_minus[1*N_S+0] = k4_minus_ref*pow(mu_B,-1*r);
    k4_minus[1*N_S+1] = k4_minus_ref;
    k4_minus[1*N_S+2] = k4_minus_ref*pow(mu_B,-1*r);
    k4_minus[1*N_S+3] = k4_minus_ref;
    k4_minus[1*N_S+4] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[1*N_S+5] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[1*N_S+6] = k4_minus_ref*pow(mu_M,1*r);
    //-------
    k4_minus[2*N_S+0] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[2*N_S+1] = k4_minus_ref*pow(mu_B,-1*r);
    k4_minus[2*N_S+2] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[2*N_S+3] = k4_minus_ref*pow(mu_B,-1*r);
    k4_minus[2*N_S+4] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[2*N_S+5] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[2*N_S+6] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    //-------
    k4_minus[3*N_S+0] = k4_minus_ref*pow(mu_B,-1*r);
    k4_minus[3*N_S+1] = k4_minus_ref;
    k4_minus[3*N_S+2] = k4_minus_ref*pow(mu_B,-1*r);
    k4_minus[3*N_S+3] = k4_minus_ref;
    k4_minus[3*N_S+4] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[3*N_S+5] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[3*N_S+6] = k4_minus_ref*pow(mu_M,1*r);
    //-------
    k4_minus[4*N_S+0] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[4*N_S+1] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[4*N_S+2] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[4*N_S+3] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[4*N_S+4] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[4*N_S+5] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[4*N_S+6] = k4_minus_ref*pow(mu_M,2*r);
    //-------
    k4_minus[5*N_S+0] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[5*N_S+1] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[5*N_S+2] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[5*N_S+3] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[5*N_S+4] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[5*N_S+5] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[5*N_S+6] = k4_minus_ref*pow(mu_M,2*r);
    //-------
    k4_minus[6*N_S+0] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[6*N_S+1] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[6*N_S+2] = k4_minus_ref*pow((mu_M/mu_B),1*r);
    k4_minus[6*N_S+3] = k4_minus_ref*pow(mu_M,1*r);
    k4_minus[6*N_S+4] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[6*N_S+5] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[6*N_S+6] = k4_minus_ref*pow(mu_M,2*r);

    }
