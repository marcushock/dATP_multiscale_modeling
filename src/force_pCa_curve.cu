//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |          Function Name                |                         //
//                    |          force_pCa_curve ()           |                         //
//                    |                                       |                         //
//--------------------------------------------------------------------------------------//
#include "force_pCa_curve.h"
#include "problemDefines.h"
//#include "experimentalDataHost.h"
#include "gpuErrchk.h"
#include <math.h>
#include <stdio.h>
#include <boost/thread.hpp>
#include <boost/atomic.hpp>
//------------------
// Functions used
//-----------------
#include "rates_trans_matrix.h"
#include "update_RUs.h"
#include "repeat_simul.h"
#include "setGPU.h"
//-------------------

static boost::atomic<int> totalThreadsFinishedMallocing(0);
static boost::mutex lock;

//------------------------------------------
// The force_pCa_Curve function definition:
//------------------------------------------

void force_pCa_curve(initParticleArgs & args,
                     unsigned long randSeed,
                     float * Fss,
                     float * M1Arrays,
                     float * M2Arrays,
                     float * M3Arrays,
                     float * CArrays,
                     float * BArrays,
                     float * SRArrays,
                     float * ATPaseArrays,
                     int cc
                    )
{

// grab a new GPU to balance load
int GPUid = getGPU();
setGPU(GPUid);
// select beginning of this loop's Force array
float * M1 = &(M1Arrays[cc * MAX_TSTEPS]);
float * M2 = &(M2Arrays[cc * MAX_TSTEPS]);
float * M3 = &(M3Arrays[cc * MAX_TSTEPS]);
float * C = &(CArrays[cc * MAX_TSTEPS]);
float * B  = &(BArrays[cc * MAX_TSTEPS]);
float * SR  = &(SRArrays[cc * MAX_TSTEPS]);
float * ATPase = &(ATPaseArrays[cc * MAX_TSTEPS]);
float * kB_plus;
gpuErrchk(cudaMallocManaged(&kB_plus, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(kB_plus, 0, sizeof(float)*N_S*N_S));
float * kB_minus;
gpuErrchk(cudaMallocManaged(&kB_minus, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(kB_minus, 0, sizeof(float)*N_S*N_S));
float * k1_plus_drug;
gpuErrchk(cudaMallocManaged(&k1_plus_drug, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k1_plus_drug, 0, sizeof(float)*N_S*N_S));
float * k1_plus_baseline;
gpuErrchk(cudaMallocManaged(&k1_plus_baseline, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k1_plus_baseline, 0, sizeof(float)*N_S*N_S));
float * k1_minus;
gpuErrchk(cudaMallocManaged(&k1_minus, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k1_minus, 0, sizeof(float)*N_S*N_S));
float * k4_plus_drug;
gpuErrchk(cudaMallocManaged(&k4_plus_drug, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k4_plus_drug, 0, sizeof(float)*N_S*N_S));
float * k4_plus_baseline;
gpuErrchk(cudaMallocManaged(&k4_plus_baseline, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k4_plus_baseline, 0, sizeof(float)*N_S*N_S));
float * k4_minus;
gpuErrchk(cudaMallocManaged(&k4_minus, sizeof(float)*N_S*N_S));
gpuErrchk(cudaMemset(k4_minus, 0, sizeof(float)*N_S*N_S));

float gamma_B = args.gamma_B;
float gamma_M = args.gamma_M;
float mu_B = args.gamma_M;
float mu_M = args.mu_M;
float kB_plus_ref = args.kB_plus_ref;
float kB_minus_ref = args.kB_minus_ref;

float k1_plus_ref_baseline = args.k1_plus_ref_baseline;
float k1_plus_ref_drug = args.k1_plus_ref_drug;

float k2_plus_baseline = args.k2_plus_baseline; 
float k2_plus_drug   = args.k2_plus_drug;

float k3_plus_baseline = args.k3_plus_baseline; 
float k3_plus_drug   = args.k3_plus_drug;

float k4_plus_ref_baseline = args.k4_plus_ref_baseline; 
float k4_plus_ref_drug = args.k4_plus_ref_drug;

float percent_drug = args.percent_drug;
float kCa_plus_ref = args.kCa_plus_ref;
float kCa_minus_ref = args.kCa_minus_ref;

float k_force_baseline = args.k_force_baseline;
float k_force_drug = args.k_force_drug;

float k_plus_SR_baseline = args.k_plus_SR_baseline;
float k_plus_SR_drug = args.k_plus_SR_drug;

float k_minus_SR = args.k_minus_SR;
float protocol = args.protocol;



//-------------------------------
//   Set rates using the input arguments
//-------------------------------
float r = args.r; // parameter defined here
float q = args.q; // parameter defined here
// float lambda = 0;
float lambda = args.lambda; 
// calculating rates for XB cycling - use Tanner 2007/ Daniel 1998/ Pate & Cooke 1989
float k1_minus_ref, k2_minus, k3_minus, k4_minus_ref;
float conc_ADP,conc_Pi, conc_ATP, x_preR, g_Ca, g_Cb, g_M1, g_M2, g_M3, delta_G_ATP, delta_G, k_xb, x_xb;
//float  A, B, C, D, M, N, P, x_b0;
//metabolite concentrations in cytosol
conc_ADP    = args.conc_ADP;        //uM, Dawson et al 1978/ Kushmerick et al 1969 (frog)
conc_ATP    = args.conc_ATP;        //uM
conc_Pi     = args.conc_Pi;         //uM
// parameter defined here
//thermodynamic parameters
//r_gas         = 8.314;      // Gas constant, J/mol*K
//tc            = 15;
//temp          = tc + 273;   //temperature in Kelvin


// other constants
float alpha = args.alpha; // parameter defined here
float beta = args.beta; // parameter defined here
float eta = args.eta; // parameter defined here
//A = 2000; 
//B = 100; // all from Tanner et al, 2007.
//C = 1;
//D = 1;
//M = 3600;
//N = 40;
//P = 20;
k_xb = args.k_xb; // parameter defined here

delta_G_ATP = args.delta_G_ATP; // units = RT
delta_G = delta_G_ATP + log(conc_ATP/(conc_ADP*conc_Pi)); // units = RT (Changing to be +, based on delta_G as an input being negative)
x_preR      = args.x_preR; // 0; XB distortion when pre-rotated.
x_xb        = args.x_xb;        // 0.075; nm, XB distortion
//x_b0        = eta * delta_G / k_xb; // xb distortion due to ATP hydrolysis



g_Cb    =  args.g_Cb                                    ;//free energy of XB state Cb
g_M1    = alpha * delta_G + k_xb * pow(x_preR,2 )     ;//free energy of XB state M1

// Note that this term could be manipulated further and a new term for k could be defined as well 
g_M2    = beta * delta_G + k_xb*pow(x_xb,2)       ;//free energy of XB state M2

g_M3    = eta* delta_G + k_xb*pow(x_xb,2)       ;//free energy of XB state M3
g_Ca    =   args.g_Ca;                                ;//free energy of XB state Ca


// to get reverse values, keep in mind that rij/rji = e^(gi - gj)
// so, to find r21 = r12/e^(g1 - g2)

//kCa_plus_ref    = 0.09;
//kCa_minus_ref   = 0.113;                    //X_kCa_minus_ref_PSO[i];
//kB_minus_ref    = 0.327;                    //X_kB_minus_ref_PSO[i];
//k1_plus_ref     = A * pow(k_xb/2*M_PI,0.5)*exp(-k_xb*pow(x_preR-x_b0,2)/2); // from tanner 2007

// Mc = M2, 

k1_minus_ref    =  k1_plus_ref_baseline/ exp(g_Cb - g_M1);//0.5 / exp(g_Cb - g_Mc);    //using vals from optimization_0227 (k1_plus = 0.615440)
// 1;//

// NEED TO DEFINE K2_MINUS!!! 
k2_minus        =  k2_plus_baseline / exp(g_M1 - g_M2); //0.5 / exp(g_Mc - g_Md);    //using vals from optimization_0227 (k1_plus = 0.615440)
// 0.000477;//
//k3_plus         = (B/pow(k_xb,.5))*(1-tanh(C*pow(k_xb,.5)*(x_xb-x_b0)))+D;        //X_k3_plus_PSO[i];
k3_minus        =  k3_plus_baseline / exp(g_M2 - g_M3) ;//0.3 / exp(g_Mc - g_Md);  //
//0.834; //

//k4_plus_ref     = pow(k_xb,0.5)*(pow(M*pow(x_xb,2),0.5)-N*x_xb)+ P;                 //X_k4_plus_PSO[i];
k4_minus_ref    =  k4_plus_ref_baseline / exp(g_M3 - delta_G); // Changing terms based on the fact that delta_G is negative
// 3.349;//

//-------------------------------------
// Call the transition rates function:
//-------------------------------------

rates_trans_matrix(N_S,
kB_plus_ref,
kB_minus_ref,
k1_plus_ref_drug,
k1_plus_ref_baseline,
k1_minus_ref,
k4_plus_ref_drug,
k4_plus_ref_baseline,
k4_minus_ref,
gamma_B,
gamma_M,
mu_B,
mu_M,
r,
q,
kB_plus,
kB_minus,
k1_plus_drug,
k1_plus_baseline,
k1_minus,
k4_plus_drug,
k4_plus_baseline,
k4_minus
);

    //-----------------------------------------------------
    // start Ca- loop i.e., to get the entire F-Ca curve:
    //-----------------------------------------------------

    float Ftemp        = 0.0;                      // is used to calculate the steady-state force at the end
    float Calc_conc_exp           = pow(10.0f,-(args.experimentalData[cc].first-6));     // Ca2+ concentration in uM
    // float kCa_plus     = Cal_conc*kCa_plus_ref;
    // float kCa_plus     = kCa_plus_ref; Removing becuase not refenced
    // float kCa_minus    = kCa_minus_ref; Removing because not referenced in this part yet 
    const int n_pCa = args.experimentalData.size();
    
    
    //---------------------------------
    // Call the repeat_simul function:
    //----------------------------------
    cudaStream_t s;
    cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking);
    totalThreadsFinishedMallocing++;
    while(totalThreadsFinishedMallocing < n_pCa){
        boost::thread::yield();
    }
repeat_simul<<<MAX_REPS/32, 32, 0, s>>>(lambda,
                                        randSeed,
                                        k1_plus_drug,
                                        k1_plus_baseline,
                                        k1_minus,
                                        k2_plus_drug,
                                        k2_plus_baseline,
                                        k2_minus,
                                        k3_plus_drug,
                                        k3_plus_baseline,
                                        k3_minus,
                                        k4_plus_drug,
                                        k4_plus_baseline,
                                        k4_minus,
                                        kB_plus,
                                        kB_minus,
                                        kCa_plus_ref,
                                        kCa_minus_ref,
                                        percent_drug,
                                        k_force_drug,
                                        k_force_baseline,
                                        k_plus_SR_drug,
                                        k_plus_SR_baseline,
                                        k_minus_SR,
                                        M1,
                                        M2,
                                        M3,
                                        C,
                                        B,
                                        SR,
                                        ATPase,
                                        cc, 
                                        protocol, 
                                        Calc_conc_exp
                                        );

    gpuErrchk(cudaStreamSynchronize(s)); // wait for device to finish repeat_simul
    gpuErrchk(cudaStreamDestroy(s));
    //--------------------------------------------------------------------------------------
    // Calculate The Steady-State Force using Impluse using data from the last 5 sec (was previously 0.5)
    // (i.e., just 100000 time steps) only using numerical trapaziodal integration
    //--------------------------------------------------------------------------------------

    for (int n = MAX_TSTEPS-1000000; n < MAX_TSTEPS-1; n++)  // time marching Originally was set to 100000
    {
        Ftemp = Ftemp+M3[n];
    }

    Fss[cc] = (Ftemp + (0.5f * M3[MAX_TSTEPS-1000001]) + (0.5f * M3[MAX_TSTEPS-1])) / 1000000.0f / MAX_REPS;    //Fss[cc] = 1;

    //--------------------------------

    // free allocated memory
gpuErrchk(cudaFree(kB_plus));
gpuErrchk(cudaFree(kB_minus));
gpuErrchk(cudaFree(k1_plus_drug));
gpuErrchk(cudaFree(k1_plus_baseline));
gpuErrchk(cudaFree(k1_minus));
gpuErrchk(cudaFree(k4_plus_drug));
gpuErrchk(cudaFree(k4_plus_baseline));
gpuErrchk(cudaFree(k4_minus));

    lock.lock();
    if(totalThreadsFinishedMallocing == n_pCa){
        totalThreadsFinishedMallocing = 0;
    }
    lock.unlock();

} // end main function
