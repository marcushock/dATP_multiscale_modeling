//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |          Function Name                |                         //
//                    |          force_pCa_curve ()           |                         //
//                    |                                       |                         //
//--------------------------------------------------------------------------------------//
#include "particles.h"
#include "problemDefines.h"
#include "setGPU.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <boost/thread.hpp>

#include "force_pCa_curve.h"
//#include "experimentalDataHost.h"
#include "gpuErrchk.h"

inline float getRand()
{
    return (float) (RANDVAL) / (float) (RAND_MAX);
}

void init_particle(initParticleArgs & args, int replicate_number)
{
    const int n_pCa = args.experimentalData.size();


    //-----------------------------------------------------------------------------------------
    // Call the force_pCa_curve function to get the force as a function of Ca++ concentrations.
    // NB: this function implicitly calls the other functions.
    //-----------------------------------------------------------------------------------------

    // allocate arrays for Force outside of kernel to use global GPU memory
    float * M3Arrays_return;
    gpuErrchk(cudaMallocManaged(&M3Arrays_return, sizeof(float)*n_pCa*MAX_TSTEPS));
    gpuErrchk(cudaMemset(M3Arrays_return, 0, sizeof(float)*n_pCa*MAX_TSTEPS));
    float * Fss_return;
    gpuErrchk(cudaMallocManaged(&Fss_return, sizeof(float)*n_pCa));
    gpuErrchk(cudaMemset(Fss_return, 0, sizeof(float)*n_pCa));
    float * M1Arrays_return;
    gpuErrchk(cudaMallocManaged(&M1Arrays_return, sizeof(float)*n_pCa*MAX_TSTEPS));
    gpuErrchk(cudaMemset(M1Arrays_return, 0, sizeof(float)*n_pCa*MAX_TSTEPS));
    float * M2Arrays_return;
    gpuErrchk(cudaMallocManaged(&M2Arrays_return, sizeof(float)*n_pCa*MAX_TSTEPS));
    gpuErrchk(cudaMemset(M2Arrays_return, 0, sizeof(float)*n_pCa*MAX_TSTEPS));
    float * CArrays_return;
    gpuErrchk(cudaMallocManaged(&CArrays_return, sizeof(float)*n_pCa*MAX_TSTEPS));
    gpuErrchk(cudaMemset(CArrays_return, 0, sizeof(float)*n_pCa*MAX_TSTEPS));
    float * BArrays_return;
    gpuErrchk(cudaMallocManaged(&BArrays_return, sizeof(float)*n_pCa*MAX_TSTEPS));
    gpuErrchk(cudaMemset(BArrays_return, 0, sizeof(float)*n_pCa*MAX_TSTEPS));
    float * SRArrays_return;
    gpuErrchk(cudaMallocManaged(&SRArrays_return, sizeof(float)*n_pCa*MAX_TSTEPS));
    gpuErrchk(cudaMemset(SRArrays_return, 0, sizeof(float)*n_pCa*MAX_TSTEPS));
    float * SSArrays_return; // State array to track the super slow state
    gpuErrchk(cudaMallocManaged(&SSArrays_return, sizeof(float)*n_pCa*MAX_TSTEPS));
    gpuErrchk(cudaMemset(SSArrays_return, 0, sizeof(float)*n_pCa*MAX_TSTEPS));
    float * ATPaseArrays_return;
    gpuErrchk(cudaMallocManaged(&ATPaseArrays_return, sizeof(float)*n_pCa*MAX_TSTEPS));
    gpuErrchk(cudaMemset(ATPaseArrays_return, 0, sizeof(float)*n_pCa*MAX_TSTEPS));

    std::cout << "Running force_pCa_curve" << std::endl;
    boost::thread_group pcaThreadGroup;
    for(int cc = 0; cc < n_pCa; ++cc)
    {
        pcaThreadGroup.add_thread(new boost::thread([&, cc]() {
            force_pCa_curve(
            args,
            RANDVAL,
            Fss_return,
            M1Arrays_return,
            M2Arrays_return,
            M3Arrays_return,
            CArrays_return,
            BArrays_return,
            SRArrays_return,
            SSArrays_return, 
            ATPaseArrays_return,
            cc);
        }));
    }
    pcaThreadGroup.join_all();

    float Fss_max_temp = 0.0;                         // is used to calculate the max s-s force
    //------------------------------------------------------------
    // Obtain the maximum vale of all the steady-state force values
    //-------------------------------------------------------------
    for(int j = 0; j < n_pCa; ++j)
    {
        if (Fss_return[j] > Fss_max_temp)
        {
            Fss_max_temp = Fss_return[j];
        }
    }

    //-------------------------------------
    // Formulate Residual/Cost Function :
    //--------------------------------------
    float residual_temp = 0.0;
    for (int j = 0; j < n_pCa; ++j)  // Ca-loop
    {
        // experimentalData[j].second == F_exp[j]
        residual_temp += pow((args.experimentalData[j].second - Fss_return[j] / Fss_max_temp),2); // normalized force
    }
    // set return residual
    float force_pCa_residual = pow(residual_temp,0.5f);
    // Note: All of these are looking at the drug parameter changes in the printed file name. 
    std::cout << "Residual: " << force_pCa_residual << std::endl;

        // std::string dataAppend =
        // //" percent " + std::to_string(args.percent_drug) +
        // //" gammaB " + std::to_string(args.gamma_B) +
        // //" gammaM " + std::to_string(args.gamma_M) +
        // //" mu_M " + std::to_string(args.mu_M) +
        // " k1_plus_ref " + std::to_string(args.k1_plus_ref_drug) +
        // " k3_plus " + std::to_string(args.k3_plus_drug) +
        // " k4_plus_ref " + std::to_string(args.k4_plus_ref_drug) +
        // " kB_plus_ref " + std::to_string(args.kB_plus_ref) +
        // " kB_minus_ref " + std::to_string(args.kB_minus_ref) +
        // //" lambda " + std::to_string(args.lambda) +
        // " kCa_plus_ref " + std::to_string(args.kCa_plus_ref) +
        // " drug " + std::to_string(args.percent_drug) +
        // " k_force " + std::to_string(args.k_force_drug) +
        // " k_plus_SR_ref " + std::to_string(args.k_plus_SR_drug) +
        // " k_minus_SR_ref " + std::to_string(args.k_minus_SR);
        std::string Force_out_Filename = ("MCMC_simulation_results/General_results/Rep_" + std::to_string(replicate_number) + "Force_out"+".csv");
        std::string ATP_out_Filename = ("MCMC_simulation_results/General_results/Rep_" + std::to_string(replicate_number) + "ATP_out"+".csv");
        std::string States_out_Filename = ("MCMC_simulation_results/General_results/Rep_" + std::to_string(replicate_number) + "States_out"+".csv");
        std::string Force_pCa_out_Filename = ("MCMC_simulation_results/General_results/Rep_" + std::to_string(replicate_number) + "Force_pCa_Optmz"+".csv");
        std::string Force_pCa_normalized_out_Filename = ("MCMC_simulation_results/General_results/Rep_" + std::to_string(replicate_number) + "Force_pCa_Optmz_Normalized"+".csv");

        /* raw force out */
        const int skipFactor = 1000;
        std::ofstream Force_out(Force_out_Filename); //opening an output stream for file *.csv
        for(int j = 0; j < MAX_TSTEPS; j+=skipFactor)
        {
            Force_out << DT*j;
            for (int cc = 0; cc < n_pCa; cc++)  // Ca-loop
            {
                Force_out << "," << M3Arrays_return[cc * MAX_TSTEPS + j]/MAX_REPS;

            }
            Force_out << std::endl;
        }
        Force_out.close();
        std::cout << "data successfully saved into the file name: " << Force_out_Filename << std::endl;

        /* raw ATP out */
        float ATPcumulative = 0.0;
        std::ofstream ATP_out(ATP_out_Filename); //opening an output stream for file *.csv
        // First timestep out will have 0 ATPase 
        ATP_out << 0;
        for (int cc = 0; cc < n_pCa; cc++){  // Ca-loop{ of all 0s to start
            ATP_out << "," << 0;
        }
        ATP_out << std::endl;
        // Then loop through the rest of the array to get the cumulative ATPases

        // Iteative through the ATPaseArray_return and print out the ATPaseArrays if it is not 0
        for (int x = 0; x < n_pCa*MAX_REPS; x++)
        {
            if (ATPaseArrays_return[x] != 0.0)
            {
                std::cout << "ATPaseArrays_return[" << x << "] = " << ATPaseArrays_return[x] << std::endl;
            }
        }


        for(int j = skipFactor; j < MAX_TSTEPS; j+=skipFactor)
        {
            ATP_out << DT*j;
            for (int cc = 0; cc < n_pCa; cc++)  // Ca-loop
            {
                // ATP_out << "," << ATPaseArrays_return[cc * MAX_TSTEPS + j]/MAX_REPS;
                ATPcumulative = 0.0;
                for (int k = 0; k < skipFactor; k+=1)
                {
                    ATPcumulative += (float) ATPaseArrays_return[cc * MAX_TSTEPS + j - k];
                    // printf("Test");
                }
                ATP_out << "," << ATPcumulative/MAX_REPS;
            }
            ATP_out << std::endl;
        }
        ATP_out.close();
        std::cout << "ATP data successfully saved into the file name: " << ATP_out_Filename << std::endl;

        std::ofstream States_out(States_out_Filename); //opening an output stream for file *.csv
        for(int j = 0; j < MAX_TSTEPS; j+=skipFactor)
        {
            States_out << DT*j;
            for (int cc = 0; cc < n_pCa; cc++)  // Ca-loop
            {
                States_out << "," << M3Arrays_return[cc * MAX_TSTEPS + j]/MAX_REPS << "," << M2Arrays_return[cc * MAX_TSTEPS + j]/MAX_REPS << "," << M1Arrays_return[cc * MAX_TSTEPS + j]/MAX_REPS << "," << CArrays_return[cc * MAX_TSTEPS + j]/MAX_REPS << "," << BArrays_return[cc * MAX_TSTEPS + j]/MAX_REPS << "," << SRArrays_return[cc * MAX_TSTEPS + j]/MAX_REPS<< "," << SSArrays_return[cc * MAX_TSTEPS + j]/MAX_REPS;

            }
            States_out << std::endl;
        }
        States_out.close();
        std::cout << "data successfully saved into the file name: " << States_out_Filename << std::endl;





        /* end raw force out */

        /* pca + pca normalized out */
        std::ofstream Force_pCa_out(Force_pCa_out_Filename); //opening an output stream for file *.csv
        std::ofstream Force_pCa_normalized_out(Force_pCa_normalized_out_Filename); //opening an output stream for file *.csv
        for (int cc = 0; cc < n_pCa; cc++)  // Ca-loop
        {
            // experimentalData[cc].first == pCa[cc]
            Force_pCa_out << args.experimentalData[cc].first << "," << Fss_return[cc] << std::endl; // write the average force at each time in a file
            Force_pCa_normalized_out << args.experimentalData[cc].first << "," << Fss_return[cc] / Fss_max_temp << std::endl; // normalized force
        }
        Force_pCa_out.close();
        std::cout << " data successfully saved into the file name: " << Force_pCa_out_Filename << std::endl;
        Force_pCa_normalized_out.close();
        std::cout << " data successfully saved into the file name: " << Force_pCa_normalized_out_Filename << std::endl;
        /* end pca + pca normalized out */

    gpuErrchk(cudaFree(Fss_return));
    gpuErrchk(cudaFree(M1Arrays_return));
    gpuErrchk(cudaFree(M2Arrays_return));
    gpuErrchk(cudaFree(M3Arrays_return));
    gpuErrchk(cudaFree(CArrays_return));
    gpuErrchk(cudaFree(BArrays_return));
    gpuErrchk(cudaFree(SRArrays_return));
    gpuErrchk(cudaFree(SSArrays_return));
    gpuErrchk(cudaFree(ATPaseArrays_return));

}
