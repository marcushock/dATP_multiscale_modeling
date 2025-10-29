#include "../src/compute_coop_factor.h" 
#include <math.h>
#include "../src/problemDefines.h"

#include "../src/gpuErrchk.h"
#include "../src/setGPU.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//--------------------------
// Function to be called
//--------------------------
#include "../src/particles.h"
#include "../src/csvReader.h"
#include "../src/argumentReader.h"

int main(int argc, const char *argv[])
{
if(argc < 2){
std::cerr << "Experimental data CSV needed as first argument." << std::endl;
return 1;
}
if(argc < 3){
std::cerr << "Argument list CSV needed as second argument." << std::endl;
return 1;
}
// std::vector< std::pair<float, float> > experimentalData = csvReader(argv[1]); // Defining a vector using std::vector < std::pair<flaot,float> > experimentalData. This is based on namespaces, and the standard library includes basic structures for vectors. 
// std::vector< std::vector<float> > argsArray = argumentReader(argv[2]); // Defining a vector using std::vector < std::vector<float> > argsArray. This is based on namespaces, and the standard library includes basic structures for vectors.
srand(SEED); //Random-Seed initialization (must be outside any loop)
std::cout << "SEED: " << SEED << std::endl;
long long startTime;
// make sure new GPUs are used
// initGPUSelection();


const int n_s = 9; 
const float kB_plus_ref = 8.9;
float kB_plus[n_s*n_s]; 
float kB_plus_new[n_s*n_s];
const float gamma_B = 45;
const float gamma_M = 21;
const float mu_B = 21;
const float mu_M = 3; 
float q = 0.5;

compute_coop_factor(n_s, kB_plus_ref, gamma_B, gamma_M, kB_plus_new, q);

kB_plus[0*n_s+0] = kB_plus_ref*pow(gamma_B,-2*q); // 0 and 0 which is B* and B* 
kB_plus[0*n_s+1] = kB_plus_ref*pow(gamma_B,-q); // 0 and 1 which is B* and C* 
kB_plus[0*n_s+2] = kB_plus_ref*pow(gamma_B,-2*q); // 0 and 2 which is B* and B 
kB_plus[0*n_s+3] = kB_plus_ref*pow(gamma_B,-q); // 0 and 3 which is B* and C 
kB_plus[0*n_s+4] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[0*n_s+5] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[0*n_s+6] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[0*n_s+7] = kB_plus_ref*pow(gamma_B,-2*q); // 0 and 7 which is B* and B** 
kB_plus[0*n_s+8] = kB_plus_ref*pow(gamma_B,-q); // 0 and 8 which is B* and C** 
//-------
kB_plus[1*n_s+0] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[1*n_s+1] = kB_plus_ref;
kB_plus[1*n_s+2] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[1*n_s+3] = kB_plus_ref;
kB_plus[1*n_s+4] = kB_plus_ref*pow(gamma_M,q);
kB_plus[1*n_s+5] = kB_plus_ref*pow(gamma_M,q);
kB_plus[1*n_s+6] = kB_plus_ref*pow(gamma_M,q);
kB_plus[1*n_s+7] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[1*n_s+8] = kB_plus_ref;
//-------
kB_plus[2*n_s+0] = kB_plus_ref*pow(gamma_B,-2*q);
kB_plus[2*n_s+1] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[2*n_s+2] = kB_plus_ref*pow(gamma_B,-2*q);
kB_plus[2*n_s+3] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[2*n_s+4] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[2*n_s+5] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[2*n_s+6] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[2*n_s+7] = kB_plus_ref*pow(gamma_B,-2*q);
kB_plus[2*n_s+8] = kB_plus_ref*pow(gamma_B,-q);
//-------
kB_plus[3*n_s+0] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[3*n_s+1] = kB_plus_ref;
kB_plus[3*n_s+2] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[3*n_s+3] = kB_plus_ref;
kB_plus[3*n_s+4] = kB_plus_ref*pow(gamma_M,q);
kB_plus[3*n_s+5] = kB_plus_ref*pow(gamma_M,q);
kB_plus[3*n_s+6] = kB_plus_ref*pow(gamma_M,q);
kB_plus[3*n_s+7] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[3*n_s+8] = kB_plus_ref;
//-------
kB_plus[4*n_s+0] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[4*n_s+1] = kB_plus_ref*pow(gamma_M,q);
kB_plus[4*n_s+2] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[4*n_s+3] = kB_plus_ref*pow(gamma_M,q);
kB_plus[4*n_s+4] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[4*n_s+5] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[4*n_s+6] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[4*n_s+7] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[4*n_s+8] = kB_plus_ref*pow(gamma_M,q);
//-------
kB_plus[5*n_s+0] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[5*n_s+1] = kB_plus_ref*pow(gamma_M,q);
kB_plus[5*n_s+2] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[5*n_s+3] = kB_plus_ref*pow(gamma_M,q);
kB_plus[5*n_s+4] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[5*n_s+5] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[5*n_s+6] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[5*n_s+7] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[5*n_s+8] = kB_plus_ref*pow(gamma_M,q);
//-------
kB_plus[6*n_s+0] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[6*n_s+1] = kB_plus_ref*pow(gamma_M,q);
kB_plus[6*n_s+2] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[6*n_s+3] = kB_plus_ref*pow(gamma_M,q);
kB_plus[6*n_s+4] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[6*n_s+5] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[6*n_s+6] = kB_plus_ref*pow(gamma_M,2*q);
kB_plus[6*n_s+7] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[6*n_s+8] = kB_plus_ref*pow(gamma_M,q);
//-------
kB_plus[7*n_s+0] = kB_plus_ref*pow(gamma_B,-2*q);
kB_plus[7*n_s+1] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[7*n_s+2] = kB_plus_ref*pow(gamma_B,-2*q);
kB_plus[7*n_s+3] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[7*n_s+4] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[7*n_s+5] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[7*n_s+6] = kB_plus_ref*pow((gamma_M/gamma_B),q);
kB_plus[7*n_s+7] = kB_plus_ref*pow(gamma_B,-2*q);
kB_plus[7*n_s+8] = kB_plus_ref*pow(gamma_B,-q);
//-------
kB_plus[8*n_s+0] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[8*n_s+1] = kB_plus_ref;
kB_plus[8*n_s+2] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[8*n_s+3] = kB_plus_ref;
kB_plus[8*n_s+4] = kB_plus_ref*pow(gamma_M,q);
kB_plus[8*n_s+5] = kB_plus_ref*pow(gamma_M,q);
kB_plus[8*n_s+6] = kB_plus_ref*pow(gamma_M,q);
kB_plus[8*n_s+7] = kB_plus_ref*pow(gamma_B,-q);
kB_plus[8*n_s+8] = kB_plus_ref;

std::cout << "Checking the kB_plus values" << std::endl;
for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_s; j++){
        float diff = fabs(kB_plus[i*n_s + j] - kB_plus_new[i*n_s + j]);
        if (diff > 1e-5){
            std::cout << "kB_plus[" << i << "][" << j << "] = " << kB_plus[i*n_s + j] << ", kB_plus_new[" << i << "][" << j << "] = " << kB_plus_new[i*n_s + j] << ", Diff = " << fabs(kB_plus[i*n_s + j] - kB_plus_new[i*n_s + j]) << std::endl;
        }
    }
}

// Testing kB_minus 
const float kB_minus_ref = 0.1;
float kB_minus[n_s*n_s];
float kB_minus_new[n_s*n_s];
q = 1;

compute_coop_factor(n_s, kB_minus_ref, gamma_B, gamma_M, kB_minus_new, q-1);

kB_minus[0*N_S+0] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[0*N_S+1] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[0*N_S+2] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[0*N_S+3] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[0*N_S+4] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[0*N_S+5] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[0*N_S+6] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[0*N_S+7] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[0*N_S+8] = kB_minus_ref*pow(1/gamma_B,q-1);
//-------
kB_minus[1*N_S+0] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[1*N_S+1] = kB_minus_ref;
kB_minus[1*N_S+2] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[1*N_S+3] = kB_minus_ref;
kB_minus[1*N_S+4] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[1*N_S+5] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[1*N_S+6] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[1*N_S+7] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[1*N_S+8] = kB_minus_ref;
//-------
kB_minus[2*N_S+0] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[2*N_S+1] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[2*N_S+2] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[2*N_S+3] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[2*N_S+4] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[2*N_S+5] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[2*N_S+6] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[2*N_S+7] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[2*N_S+8] = kB_minus_ref*pow(1/gamma_B,q-1);
//-------
kB_minus[3*N_S+0] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[3*N_S+1] = kB_minus_ref;
kB_minus[3*N_S+2] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[3*N_S+3] = kB_minus_ref;
kB_minus[3*N_S+4] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[3*N_S+5] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[3*N_S+6] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[3*N_S+7] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[3*N_S+8] = kB_minus_ref;
//-------
kB_minus[4*N_S+0] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[4*N_S+1] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[4*N_S+2] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[4*N_S+3] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[4*N_S+4] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[4*N_S+5] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[4*N_S+6] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[4*N_S+7] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[4*N_S+8] = kB_minus_ref*pow(gamma_M,q-1);
//-------
kB_minus[5*N_S+0] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[5*N_S+1] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[5*N_S+2] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[5*N_S+3] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[5*N_S+4] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[5*N_S+5] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[5*N_S+6] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[5*N_S+7] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[5*N_S+8] = kB_minus_ref*pow(gamma_M,q-1);
//-------
kB_minus[6*N_S+0] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[6*N_S+1] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[6*N_S+2] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[6*N_S+3] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[6*N_S+4] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[6*N_S+5] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[6*N_S+6] = kB_minus_ref*pow(pow(gamma_M,2),q-1);
kB_minus[6*N_S+7] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[6*N_S+8] = kB_minus_ref*pow(gamma_M,q-1);
//-------
kB_minus[7*N_S+0] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[7*N_S+1] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[7*N_S+2] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[7*N_S+3] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[7*N_S+4] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[7*N_S+5] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[7*N_S+6] = kB_minus_ref*pow(gamma_M/gamma_B,q-1);
kB_minus[7*N_S+7] = kB_minus_ref*pow(pow(gamma_B,-2),q-1);
kB_minus[7*N_S+8] = kB_minus_ref*pow(1/gamma_B,q-1);
//-------
kB_minus[8*N_S+0] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[8*N_S+1] = kB_minus_ref;
kB_minus[8*N_S+2] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[8*N_S+3] = kB_minus_ref;
kB_minus[8*N_S+4] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[8*N_S+5] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[8*N_S+6] = kB_minus_ref*pow(gamma_M,q-1);
kB_minus[8*N_S+7] = kB_minus_ref*pow(1/gamma_B,q-1);
kB_minus[8*N_S+8] = kB_minus_ref;


std::cout<<"Checking the kB_minus values" << std::endl;
// Check the kB_minus values
for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_s; j++){
        float diff = fabs(kB_minus[i*n_s + j] - kB_minus_new[i*n_s + j]);
        if (diff > 1e-5){
            std::cout << "kB_minus[" << i << "][" << j << "] = " << kB_minus[i*n_s + j] << ", kB_minus_new[" << i << "][" << j << "] = " << kB_minus_new[i*n_s + j] << ", Diff = " << fabs(kB_minus[i*n_s + j] - kB_minus_new[i*n_s + j]) << std::endl;
        }
    }
}

float k1_plus_ref_baseline = 0.0025;
float k1_plus_baseline[n_s*n_s];
float k1_plus_baseline_new[n_s*n_s];
float r = 0.5;
compute_coop_factor(n_s, k1_plus_ref_baseline, mu_B, mu_M, k1_plus_baseline_new, r);

float k1_plus_ref_drug = 0.0056;
float k1_plus_drug[n_s*n_s];
float k1_plus_drug_new[n_s*n_s];
compute_coop_factor(n_s, k1_plus_ref_drug, mu_B, mu_M, k1_plus_drug_new, r);

float k1_minus_ref = 0.0005;
float k1_minus[n_s*n_s];
float k1_minus_new[n_s*n_s];
compute_coop_factor(n_s, k1_minus_ref, mu_B, mu_M, k1_minus_new, r-1);

float k4_minus_ref = 0.05;
float k4_minus[n_s*n_s];
float k4_minus_new[n_s*n_s];
compute_coop_factor(n_s, k4_minus_ref, mu_B, mu_M, k4_minus_new, r);

float k4_plus_ref_baseline = 0.01;
float k4_plus_baseline[n_s*n_s];
float k4_plus_baseline_new[n_s*n_s];
compute_coop_factor(n_s, k4_plus_ref_baseline, mu_B, mu_M, k4_plus_baseline_new, r -1);

float k4_plus_ref_drug = 0.02;
float k4_plus_drug[n_s*n_s];
float k4_plus_drug_new[n_s*n_s];
compute_coop_factor(n_s, k4_plus_ref_drug, mu_B, mu_M, k4_plus_drug_new, r -1);

//--------------------------------------------------
    // Step 3: Build the k1_plus_reference [ns*N_S+ns] matrix
    //--------------------------------------------------
    k1_plus_baseline[0*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[0*N_S+1] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[0*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[0*N_S+3] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[0*N_S+4] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[0*N_S+5] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[0*N_S+6] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[0*N_S+7] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[0*N_S+8] = k1_plus_ref_baseline*pow(mu_B,-r);
    //-------
    k1_plus_baseline[1*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[1*N_S+1] = k1_plus_ref_baseline;
    k1_plus_baseline[1*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[1*N_S+3] = k1_plus_ref_baseline;
    k1_plus_baseline[1*N_S+4] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[1*N_S+5] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[1*N_S+6] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[1*N_S+7] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[1*N_S+8] = k1_plus_ref_baseline;
    //-------
    k1_plus_baseline[2*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[2*N_S+1] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[2*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[2*N_S+3] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[2*N_S+4] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[2*N_S+5] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[2*N_S+6] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[2*N_S+7] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[2*N_S+8] = k1_plus_ref_baseline*pow(mu_B,-r);
    //-------
    k1_plus_baseline[3*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[3*N_S+1] = k1_plus_ref_baseline;
    k1_plus_baseline[3*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[3*N_S+3] = k1_plus_ref_baseline;
    k1_plus_baseline[3*N_S+4] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[3*N_S+5] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[3*N_S+6] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[3*N_S+7] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[3*N_S+8] = k1_plus_ref_baseline;
    //-------
    k1_plus_baseline[4*N_S+0] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[4*N_S+1] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[4*N_S+2] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[4*N_S+3] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[4*N_S+4] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[4*N_S+5] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[4*N_S+6] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[4*N_S+7] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[4*N_S+8] = k1_plus_ref_baseline*pow(mu_M,r);
    //-------
    k1_plus_baseline[5*N_S+0] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[5*N_S+1] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[5*N_S+2] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[5*N_S+3] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[5*N_S+4] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[5*N_S+5] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[5*N_S+6] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[5*N_S+7] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[5*N_S+8] = k1_plus_ref_baseline*pow(mu_M,r);
    //-------
    k1_plus_baseline[6*N_S+0] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[6*N_S+1] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[6*N_S+2] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[6*N_S+3] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[6*N_S+4] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[6*N_S+5] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[6*N_S+6] = k1_plus_ref_baseline*pow(mu_M,2*r);
    k1_plus_baseline[6*N_S+7] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[6*N_S+8] = k1_plus_ref_baseline*pow(mu_M,r);
    //-------
    k1_plus_baseline[7*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[7*N_S+1] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[7*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[7*N_S+3] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[7*N_S+4] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[7*N_S+5] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[7*N_S+6] = k1_plus_ref_baseline*pow((mu_M/mu_B),r);
    k1_plus_baseline[7*N_S+7] = k1_plus_ref_baseline*pow(mu_B,-2*r);
    k1_plus_baseline[7*N_S+8] = k1_plus_ref_baseline*pow(mu_B,-r);
    //-------
    k1_plus_baseline[8*N_S+0] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[8*N_S+1] = k1_plus_ref_baseline;
    k1_plus_baseline[8*N_S+2] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[8*N_S+3] = k1_plus_ref_baseline;
    k1_plus_baseline[8*N_S+4] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[8*N_S+5] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[8*N_S+6] = k1_plus_ref_baseline*pow(mu_M,r);
    k1_plus_baseline[8*N_S+7] = k1_plus_ref_baseline*pow(mu_B,-r);
    k1_plus_baseline[8*N_S+8] = k1_plus_ref_baseline;
    //--------------------------------------------------
    // Step 3: Build the k1_plus_drug [ns*N_S+ns] matrix
    //--------------------------------------------------
    k1_plus_drug[0*N_S+0] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[0*N_S+1] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[0*N_S+2] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[0*N_S+3] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[0*N_S+4] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[0*N_S+5] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[0*N_S+6] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[0*N_S+7] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[0*N_S+8] = k1_plus_ref_drug*pow(mu_B,-r);
    //-------
    k1_plus_drug[1*N_S+0] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[1*N_S+1] = k1_plus_ref_drug;
    k1_plus_drug[1*N_S+2] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[1*N_S+3] = k1_plus_ref_drug;
    k1_plus_drug[1*N_S+4] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[1*N_S+5] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[1*N_S+6] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[1*N_S+7] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[1*N_S+8] = k1_plus_ref_drug;
    //-------
    k1_plus_drug[2*N_S+0] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[2*N_S+1] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[2*N_S+2] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[2*N_S+3] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[2*N_S+4] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[2*N_S+5] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[2*N_S+6] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[2*N_S+7] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[2*N_S+8] = k1_plus_ref_drug*pow(mu_B,-r);
    //-------
    k1_plus_drug[3*N_S+0] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[3*N_S+1] = k1_plus_ref_drug;
    k1_plus_drug[3*N_S+2] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[3*N_S+3] = k1_plus_ref_drug;
    k1_plus_drug[3*N_S+4] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[3*N_S+5] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[3*N_S+6] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[3*N_S+7] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[3*N_S+8] = k1_plus_ref_drug;
    //-------
    k1_plus_drug[4*N_S+0] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[4*N_S+1] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[4*N_S+2] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[4*N_S+3] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[4*N_S+4] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[4*N_S+5] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[4*N_S+6] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[4*N_S+7] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[4*N_S+8] = k1_plus_ref_drug*pow(mu_M,r);
    //-------
    k1_plus_drug[5*N_S+0] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[5*N_S+1] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[5*N_S+2] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[5*N_S+3] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[5*N_S+4] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[5*N_S+5] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[5*N_S+6] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[5*N_S+7] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[5*N_S+8] = k1_plus_ref_drug*pow(mu_M,r);
    //-------
    k1_plus_drug[6*N_S+0] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[6*N_S+1] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[6*N_S+2] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[6*N_S+3] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[6*N_S+4] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[6*N_S+5] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[6*N_S+6] = k1_plus_ref_drug*pow(mu_M,2*r);
    k1_plus_drug[6*N_S+7] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[6*N_S+8] = k1_plus_ref_drug*pow(mu_M,r);
    //-------
    k1_plus_drug[7*N_S+0] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[7*N_S+1] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[7*N_S+2] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[7*N_S+3] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[7*N_S+4] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[7*N_S+5] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[7*N_S+6] = k1_plus_ref_drug*pow((mu_M/mu_B),r);
    k1_plus_drug[7*N_S+7] = k1_plus_ref_drug*pow(mu_B,-2*r);
    k1_plus_drug[7*N_S+8] = k1_plus_ref_drug*pow(mu_B,-r);
    //-------
    k1_plus_drug[8*N_S+0] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[8*N_S+1] = k1_plus_ref_drug;
    k1_plus_drug[8*N_S+2] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[8*N_S+3] = k1_plus_ref_drug;
    k1_plus_drug[8*N_S+4] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[8*N_S+5] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[8*N_S+6] = k1_plus_ref_drug*pow(mu_M,r);
    k1_plus_drug[8*N_S+7] = k1_plus_ref_drug*pow(mu_B,-r);
    k1_plus_drug[8*N_S+8] = k1_plus_ref_drug;
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
    k1_minus[0*N_S+7] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[0*N_S+8] = k1_minus_ref*pow(1/mu_B,r-1);
    //-------
    k1_minus[1*N_S+0] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[1*N_S+1] = k1_minus_ref;
    k1_minus[1*N_S+2] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[1*N_S+3] = k1_minus_ref;
    k1_minus[1*N_S+4] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[1*N_S+5] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[1*N_S+6] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[1*N_S+7] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[1*N_S+8] = k1_minus_ref;
    //-------
    k1_minus[2*N_S+0] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[2*N_S+1] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[2*N_S+2] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[2*N_S+3] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[2*N_S+4] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[2*N_S+5] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[2*N_S+6] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[2*N_S+7] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[2*N_S+8] = k1_minus_ref*pow(1/mu_B,r-1);
    //-------
    k1_minus[3*N_S+0] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[3*N_S+1] = k1_minus_ref;
    k1_minus[3*N_S+2] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[3*N_S+3] = k1_minus_ref;
    k1_minus[3*N_S+4] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[3*N_S+5] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[3*N_S+6] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[3*N_S+7] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[3*N_S+8] = k1_minus_ref;
    //-------
    k1_minus[4*N_S+0] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[4*N_S+1] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[4*N_S+2] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[4*N_S+3] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[4*N_S+4] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[4*N_S+5] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[4*N_S+6] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[4*N_S+7] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[4*N_S+8] = k1_minus_ref*pow(mu_M,r-1);
    //-------
    k1_minus[5*N_S+0] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[5*N_S+1] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[5*N_S+2] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[5*N_S+3] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[5*N_S+4] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[5*N_S+5] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[5*N_S+6] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[5*N_S+7] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[5*N_S+8] = k1_minus_ref*pow(mu_M,r-1);
    //-------
    k1_minus[6*N_S+0] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[6*N_S+1] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[6*N_S+2] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[6*N_S+3] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[6*N_S+4] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[6*N_S+5] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[6*N_S+6] = k1_minus_ref*pow(pow(mu_M,2),r-1);
    k1_minus[6*N_S+7] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[6*N_S+8] = k1_minus_ref*pow(mu_M,r-1);
    //-------
    k1_minus[7*N_S+0] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[7*N_S+1] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[7*N_S+2] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[7*N_S+3] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[7*N_S+4] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[7*N_S+5] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[7*N_S+6] = k1_minus_ref*pow(mu_M/mu_B,r-1);
    k1_minus[7*N_S+7] = k1_minus_ref*pow(pow(mu_B,-2),r-1);
    k1_minus[7*N_S+8] = k1_minus_ref*pow(1/mu_B,r-1);
    //-------
    k1_minus[8*N_S+0] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[8*N_S+1] = k1_minus_ref;
    k1_minus[8*N_S+2] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[8*N_S+3] = k1_minus_ref;
    k1_minus[8*N_S+4] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[8*N_S+5] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[8*N_S+6] = k1_minus_ref*pow(mu_M,r-1);
    k1_minus[8*N_S+7] = k1_minus_ref*pow(1/mu_B,r-1);
    k1_minus[8*N_S+8] = k1_minus_ref;
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
    k4_plus_baseline[0*N_S+7] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[0*N_S+8] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    //-------
    k4_plus_baseline[1*N_S+0] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[1*N_S+1] = k4_plus_ref_baseline;
    k4_plus_baseline[1*N_S+2] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[1*N_S+3] = k4_plus_ref_baseline;
    k4_plus_baseline[1*N_S+4] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[1*N_S+5] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[1*N_S+6] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[1*N_S+7] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[1*N_S+8] = k4_plus_ref_baseline;
    //-------
    k4_plus_baseline[2*N_S+0] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[2*N_S+1] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[2*N_S+2] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[2*N_S+3] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[2*N_S+4] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[2*N_S+5] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[2*N_S+6] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[2*N_S+7] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[2*N_S+8] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    //-------
    k4_plus_baseline[3*N_S+0] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[3*N_S+1] = k4_plus_ref_baseline;
    k4_plus_baseline[3*N_S+2] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[3*N_S+3] = k4_plus_ref_baseline;
    k4_plus_baseline[3*N_S+4] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[3*N_S+5] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[3*N_S+6] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[3*N_S+7] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[3*N_S+8] = k4_plus_ref_baseline;
    //-------
    k4_plus_baseline[4*N_S+0] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[4*N_S+1] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[4*N_S+2] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[4*N_S+3] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[4*N_S+4] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[4*N_S+5] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[4*N_S+6] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[4*N_S+7] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[4*N_S+8] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    //-------
    k4_plus_baseline[5*N_S+0] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[5*N_S+1] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[5*N_S+2] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[5*N_S+3] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[5*N_S+4] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[5*N_S+5] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[5*N_S+6] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[5*N_S+7] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[5*N_S+8] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    //-------
    k4_plus_baseline[6*N_S+0] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[6*N_S+1] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[6*N_S+2] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[6*N_S+3] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[6*N_S+4] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[6*N_S+5] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[6*N_S+6] = k4_plus_ref_baseline*pow(pow(mu_M,2),(r-1));
    k4_plus_baseline[6*N_S+7] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[6*N_S+8] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    //-------
    k4_plus_baseline[7*N_S+0] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[7*N_S+1] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[7*N_S+2] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[7*N_S+3] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[7*N_S+4] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[7*N_S+5] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[7*N_S+6] = k4_plus_ref_baseline*pow(mu_M/mu_B,(r-1));
    k4_plus_baseline[7*N_S+7] = k4_plus_ref_baseline*pow(pow(mu_B,-2),(r-1));
    k4_plus_baseline[7*N_S+8] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    //-------
    k4_plus_baseline[8*N_S+0] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[8*N_S+1] = k4_plus_ref_baseline;
    k4_plus_baseline[8*N_S+2] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[8*N_S+3] = k4_plus_ref_baseline;
    k4_plus_baseline[8*N_S+4] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[8*N_S+5] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[8*N_S+6] = k4_plus_ref_baseline*pow(mu_M,(r-1));
    k4_plus_baseline[8*N_S+7] = k4_plus_ref_baseline*pow(1/mu_B,(r-1));
    k4_plus_baseline[8*N_S+8] = k4_plus_ref_baseline;



    k4_plus_drug[0*N_S+0] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[0*N_S+1] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[0*N_S+2] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[0*N_S+3] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[0*N_S+4] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[0*N_S+5] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[0*N_S+6] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[0*N_S+7] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[0*N_S+8] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    //-------
    k4_plus_drug[1*N_S+0] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[1*N_S+1] = k4_plus_ref_drug;
    k4_plus_drug[1*N_S+2] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[1*N_S+3] = k4_plus_ref_drug;
    k4_plus_drug[1*N_S+4] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[1*N_S+5] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[1*N_S+6] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[1*N_S+7] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[1*N_S+8] = k4_plus_ref_drug;
    //-------
    k4_plus_drug[2*N_S+0] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[2*N_S+1] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[2*N_S+2] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[2*N_S+3] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[2*N_S+4] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[2*N_S+5] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[2*N_S+6] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[2*N_S+7] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[2*N_S+8] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    //-------
    k4_plus_drug[3*N_S+0] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[3*N_S+1] = k4_plus_ref_drug;
    k4_plus_drug[3*N_S+2] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[3*N_S+3] = k4_plus_ref_drug;
    k4_plus_drug[3*N_S+4] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[3*N_S+5] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[3*N_S+6] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[3*N_S+7] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[3*N_S+8] = k4_plus_ref_drug;
    //-------
    k4_plus_drug[4*N_S+0] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[4*N_S+1] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[4*N_S+2] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[4*N_S+3] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[4*N_S+4] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[4*N_S+5] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[4*N_S+6] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[4*N_S+7] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[4*N_S+8] = k4_plus_ref_drug*pow(mu_M,(r-1));
    //-------
    k4_plus_drug[5*N_S+0] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[5*N_S+1] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[5*N_S+2] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[5*N_S+3] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[5*N_S+4] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[5*N_S+5] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[5*N_S+6] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[5*N_S+7] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[5*N_S+8] = k4_plus_ref_drug*pow(mu_M,(r-1));
    //-------
    k4_plus_drug[6*N_S+0] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[6*N_S+1] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[6*N_S+2] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[6*N_S+3] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[6*N_S+4] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[6*N_S+5] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[6*N_S+6] = k4_plus_ref_drug*pow(pow(mu_M,2),(r-1));
    k4_plus_drug[6*N_S+7] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[6*N_S+8] = k4_plus_ref_drug*pow(mu_M,(r-1));
    //-------
    k4_plus_drug[7*N_S+0] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[7*N_S+1] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[7*N_S+2] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[7*N_S+3] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[7*N_S+4] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[7*N_S+5] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[7*N_S+6] = k4_plus_ref_drug*pow(mu_M/mu_B,(r-1));
    k4_plus_drug[7*N_S+7] = k4_plus_ref_drug*pow(pow(mu_B,-2),(r-1));
    k4_plus_drug[7*N_S+8] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    //-------
    k4_plus_drug[8*N_S+0] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[8*N_S+1] = k4_plus_ref_drug;
    k4_plus_drug[8*N_S+2] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[8*N_S+3] = k4_plus_ref_drug;
    k4_plus_drug[8*N_S+4] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[8*N_S+5] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[8*N_S+6] = k4_plus_ref_drug*pow(mu_M,(r-1));
    k4_plus_drug[8*N_S+7] = k4_plus_ref_drug*pow(1/mu_B,(r-1));
    k4_plus_drug[8*N_S+8] = k4_plus_ref_drug;

    //-----------------------------------

    //--------------------------------------------------
    // Step 6: Build the k4_minus [ns*N_S+ns] matrix
    //--------------------------------------------------
    k4_minus[0*N_S+0] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[0*N_S+1] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[0*N_S+2] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[0*N_S+3] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[0*N_S+4] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[0*N_S+5] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[0*N_S+6] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[0*N_S+7] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[0*N_S+8] = k4_minus_ref*pow(mu_B,-r);
    //-------
    k4_minus[1*N_S+0] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[1*N_S+1] = k4_minus_ref;
    k4_minus[1*N_S+2] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[1*N_S+3] = k4_minus_ref;
    k4_minus[1*N_S+4] = k4_minus_ref*pow(mu_M,r);
    k4_minus[1*N_S+5] = k4_minus_ref*pow(mu_M,r);
    k4_minus[1*N_S+6] = k4_minus_ref*pow(mu_M,r);
    k4_minus[1*N_S+7] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[1*N_S+8] = k4_minus_ref;
    //-------
    k4_minus[2*N_S+0] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[2*N_S+1] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[2*N_S+2] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[2*N_S+3] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[2*N_S+4] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[2*N_S+5] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[2*N_S+6] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[2*N_S+7] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[2*N_S+8] = k4_minus_ref*pow(mu_B,-r);
    //-------
    k4_minus[3*N_S+0] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[3*N_S+1] = k4_minus_ref;
    k4_minus[3*N_S+2] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[3*N_S+3] = k4_minus_ref;
    k4_minus[3*N_S+4] = k4_minus_ref*pow(mu_M,r);
    k4_minus[3*N_S+5] = k4_minus_ref*pow(mu_M,r);
    k4_minus[3*N_S+6] = k4_minus_ref*pow(mu_M,r);
    k4_minus[3*N_S+7] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[3*N_S+8] = k4_minus_ref;
    //-------
    k4_minus[4*N_S+0] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[4*N_S+1] = k4_minus_ref*pow(mu_M,r);
    k4_minus[4*N_S+2] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[4*N_S+3] = k4_minus_ref*pow(mu_M,r);
    k4_minus[4*N_S+4] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[4*N_S+5] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[4*N_S+6] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[4*N_S+7] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[4*N_S+8] = k4_minus_ref*pow(mu_M,r);
    //-------
    k4_minus[5*N_S+0] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[5*N_S+1] = k4_minus_ref*pow(mu_M,r);
    k4_minus[5*N_S+2] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[5*N_S+3] = k4_minus_ref*pow(mu_M,r);
    k4_minus[5*N_S+4] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[5*N_S+5] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[5*N_S+6] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[5*N_S+7] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[5*N_S+8] = k4_minus_ref*pow(mu_M,r);
    //-------
    k4_minus[6*N_S+0] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[6*N_S+1] = k4_minus_ref*pow(mu_M,r);
    k4_minus[6*N_S+2] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[6*N_S+3] = k4_minus_ref*pow(mu_M,r);
    k4_minus[6*N_S+4] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[6*N_S+5] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[6*N_S+6] = k4_minus_ref*pow(mu_M,2*r);
    k4_minus[6*N_S+7] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[6*N_S+8] = k4_minus_ref*pow(mu_M,r);
    //-------
    k4_minus[7*N_S+0] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[7*N_S+1] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[7*N_S+2] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[7*N_S+3] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[7*N_S+4] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[7*N_S+5] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[7*N_S+6] = k4_minus_ref*pow((mu_M/mu_B),r);
    k4_minus[7*N_S+7] = k4_minus_ref*pow(mu_B,-2*r);
    k4_minus[7*N_S+8] = k4_minus_ref*pow(mu_B,-r);
    //------
    k4_minus[8*N_S+0] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[8*N_S+1] = k4_minus_ref;
    k4_minus[8*N_S+2] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[8*N_S+3] = k4_minus_ref;
    k4_minus[8*N_S+4] = k4_minus_ref*pow(mu_M,r);
    k4_minus[8*N_S+5] = k4_minus_ref*pow(mu_M,r);
    k4_minus[8*N_S+6] = k4_minus_ref*pow(mu_M,r);
    k4_minus[8*N_S+7] = k4_minus_ref*pow(mu_B,-r);
    k4_minus[8*N_S+8] = k4_minus_ref;



std::cout<<"Checking the k1_plus_drug values" << std::endl;
// Check the kB_minus values
for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_s; j++){
        float diff = fabs(k1_plus_baseline[i*n_s + j] - k1_plus_baseline_new[i*n_s + j]);
        if (diff > 1e-5){
            std::cout << "k1_plus_baseline[" << i << "][" << j << "] = " << k1_plus_baseline[i*n_s + j] << ", k1_plus_baseline_new[" << i << "][" << j << "] = " << k1_plus_baseline_new[i*n_s + j] << ", Diff = " << fabs(k1_plus_baseline[i*n_s + j] - k1_plus_baseline_new[i*n_s + j]) << std::endl;
        }
    }
}

std::cout<<"Checking the k1_plus_drug values" << std::endl;
// Check the k1_plus_drug values
for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_s; j++){
        float diff = fabs(k1_plus_drug[i*n_s + j] - k1_plus_drug_new[i*n_s + j]);
        if (diff > 1e-5){
            std::cout << "k1_plus_drug[" << i << "][" << j << "] = " << k1_plus_drug[i*n_s + j] << ", k1_plus_drug_new[" << i << "][" << j << "] = " << k1_plus_drug_new[i*n_s + j] << ", Diff = " << fabs(k1_plus_drug[i*n_s + j] - k1_plus_drug_new[i*n_s + j]) << std::endl;
        }
    }
}

std::cout<<"Checking the k1_minus values" << std::endl;
// Check the k1_minus values
for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_s; j++){
        float diff = fabs(k1_minus[i*n_s + j] - k1_minus_new[i*n_s + j]);
        if (diff > 1e-5){
            std::cout << "k1_minus[" << i << "][" << j << "] = " << k1_minus[i*n_s + j] << ", k1_minus_new[" << i << "][" << j << "] = " << k1_minus_new[i*n_s + j] << ", Diff = " << fabs(k1_minus[i*n_s + j] - k1_minus_new[i*n_s + j]) << std::endl;
        }
    }
}

std::cout<<"Checking the k4_plus_baseline values" << std::endl;
// Check the k4_plus_baseline values
for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_s; j++){
        float diff = fabs(k4_plus_baseline[i*n_s + j] - k4_plus_baseline_new[i*n_s + j]);
        if (diff > 1e-5){
            std::cout << "k4_plus_baseline[" << i << "][" << j << "] = " << k4_plus_baseline[i*n_s + j] << ", k4_plus_baseline_new[" << i << "][" << j << "] = " << k4_plus_baseline_new[i*n_s + j] << ", Diff = " << fabs(k4_plus_baseline[i*n_s + j] - k4_plus_baseline_new[i*n_s + j]) << std::endl;
        }
    }
}

std::cout<<"Checking the k4_plus_drug values" << std::endl;
// Check the k4_plus_drug values
for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_s; j++){
        float diff = fabs(k4_plus_drug[i*n_s + j] - k4_plus_drug_new[i*n_s + j]);
        if (diff > 1e-5){
            std::cout << "k4_plus_drug[" << i << "][" << j << "] = " << k4_plus_drug[i*n_s + j] << ", k4_plus_drug_new[" << i << "][" << j << "] = " << k4_plus_drug_new[i*n_s + j] << ", Diff = " << fabs(k4_plus_drug[i*n_s + j] - k4_plus_drug_new[i*n_s + j]) << std::endl;
        }
    }
}

std::cout<<"Checking the k4_minus values" << std::endl;
// Check the k4_minus values
for (int i = 0; i < n_s; i++){
    for (int j = 0; j < n_s; j++){
        float diff = fabs(k4_minus[i*n_s + j] - k4_minus_new[i*n_s + j]);
        if (diff > 1e-5){
            std::cout << "k4_minus[" << i << "][" << j << "] = " << k4_minus[i*n_s + j] << ", k4_minus_new[" << i << "][" << j << "] = " << k4_minus_new[i*n_s + j] << ", Diff = " << fabs(k4_minus[i*n_s + j] - k4_minus_new[i*n_s + j]) << std::endl;
        }
    }
}

return 0;

}