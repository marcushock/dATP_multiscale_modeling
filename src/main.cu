//--------------------------------------------------------------------------------------//
//                     University of California San Diego                               //
//                         Dept. of Bioengineering                                      //
//                      Cardiac Mechanics Research Group                                //
//                         PI: Prof. Andrew McCulloch                                   //
//--------------------------------------------------------------------------------------//
// Authors: Kim McCabe, Yasser Aboelkassem,  & Andrew McCulloch                         //
// Year  :  08/2017    									//
// Updated: Abby Teitgen 05/2023                                                        //
//-----------------------------                                                         //
//          This code uses the "Markov Chain Monte Carlo" algorithm to solve            //
//           cardiac thin filament activation / crossbridge cycling problem.            //
//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |      The TF/XB 10 State Model         |                         //
//                    |                                       |                         //
//--------------------------------------------------------------------------------------//
#include "problemDefines.h"

#include "gpuErrchk.h"
#include "setGPU.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//--------------------------
// Function to be called
//--------------------------
#include "particles.h"
#include "csvReader.h"
#include "argumentReader.h"

//-------------------------
// Main body code
//------------------------
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
std::vector< std::pair<float, float> > experimentalData = csvReader(argv[1]); // Defining a vector using std::vector < std::pair<flaot,float> > experimentalData. This is based on namespaces, and the standard library includes basic structures for vectors. 
std::vector< std::vector<float> > argsArray = argumentReader(argv[2]); // Defining a vector using std::vector < std::vector<float> > argsArray. This is based on namespaces, and the standard library includes basic structures for vectors.
srand(SEED); //Random-Seed initialization (must be outside any loop)
std::cout << "SEED: " << SEED << std::endl;
long long startTime;
// make sure new GPUs are used
initGPUSelection();

//---------------------------------------------
// Model reference parameters that we need to optimize
// Note: atof converts strings into doubles 
// The input arguments will define whether a twitch or 
// force pCa is ran. 
// Twitch = 0 
// Force pca = 1
//--------------------------------------------

// Loop through the vectors in the parametersArray and carry out simulations with each parameter set
std::vector<float> argsVector;
for(int i = 0; i < argsArray.size(); ++i){
    startTime = time(NULL);
    argsVector = argsArray[i];
    initParticleArgs args = initParticleArgs(experimentalData, argsVector);
    init_particle(args, i);
    std::cout << "One iteration runtime: " << (time(NULL)-startTime) << " second(s)" << std::endl;
}

// float gamma_B = atof(argv[2]); // [unitless] - RU-RU cooperative coefficient
// float gamma_M = atof(argv[3]); // [unitless] - XB-RU/RU-XB coopcoefficient (Note: gamma_M = mu_B)
// float mu_M = atof(argv[4]); // [unitless] - Inter-RU XB-XB cooperative coefficient
// float k1_plus_ref = atof(argv[5]);
// float k3_plus = atof(argv[6]);
// float k4_plus_ref = atof(argv[7]);
// float kB_plus_ref = atof(argv[8]); //               - (p = plus)
// float kB_minus_ref = atof(argv[9]); //               - (m = minus)
// float lambda = atof(argv[10]); // [unitless] must be between (0,1)
// float kCa_plus_ref = atof(argv[11]);
// float kCa_minus_ref = atof(argv[12]);
// float percent_dATP = atof(argv[13]);
// float k_force = atof(argv[14]);
// float k_plus_SR_ref = atof(argv[15]);
// float k_minus_SR_ref = atof(argv[16]);
// float protocol = atoi(argv[17]);


// initParticleArgs args = initParticleArgs(experimentalData,
// gamma_B,
// gamma_M,
// mu_M,
// k1_plus_ref,
// k3_plus,
// k4_plus_ref,
// kB_plus_ref,
// kB_minus_ref,
// lambda,
// kCa_plus_ref,
// kCa_minus_ref,
// percent_dATP,
// k_force,
// k_plus_SR_ref,
// k_minus_SR_ref, 
// protocol 
// );
// init_particle(args);
// std::cout << "One iteration runtime: " << (time(NULL)-startTime) << " second(s)" << std::endl;

return 0;

} // end main function
