//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |          Function Name                |                         //
//                    |         repeat_simul()                |                         //
//                    |                                       |                         //
//--------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------
// This function is used to update the states of each RUs based on the Markov step
//---------------------------------------------------------------------------------
// Input|
//-------
// as shown below
//--------
// Output|
//---------
// Force-Ca curve at at a given Ca value
//--------------------------------------------------------------------------%
#include "repeat_simul.h"
#include "problemDefines.h"
#include "lin_interp_ca.h"
//--------------------------------------
// functions used within this function
//--------------------------------------
#include "rates_trans_matrix.h"
#include "update_RUs.h"
#include "genrand.h"
#include <stdio.h>
//-----------------------------------------------
// This function definition
//-----------------------------------------------

__global__ void repeat_simul(float lambda,
const unsigned long randSeed,
float * k4_plus_dATP,
float * k4_plus_ATP,
float * k4_minus,
float k3_plus_dATP,
float k3_plus_ATP,
float k3_minus,
float * k2_plus_dATP,
float * k2_plus_ATP,
float * k2_minus,
float * kB_plus,
float * kB_minus,
float kCa_plus_ref,
float kCa_minus,
float percent_dATP,
float k_force_dATP,
float k_force_ATP,
float k_plus_SR_ref_dATP,
float k_plus_SR_ref_ATP,
float k_minus_SR_ref,
float * Force,
float * Mc,
float * C,
float * B,
float * SR,
int cc
)

{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    /* initialize random number generation per thread */
    float randNum[N_RU];
    float rand_dATP[N_RU];
    int RU[N_RU];
    bool caRU[N_RU];
    curandState_t state;
    float f;
    int flag = 0;
    float f_prev;
    float SR_prev;
    float kCa_plus;
    float calcium;
    float current_time;


    curand_init(randSeed, index, 0, &state);    
//--------------------------------------
    // start repeat loop i.e., using r-index
    //--------------------------------------

    //reset RUs again to B0
    memset(RU, 0, sizeof(int)*N_RU);
    memset(caRU, 0, sizeof(bool)*N_RU);
    RU[0]=2;
    RU[N_RU-1]=2;
    
    for(int i = 1; i < N_RU-1; ++i)
    {
    	RU[i]=0;
    }
    //------------------------------------
    // start time loop i.e., using n-index
    //------------------------------------
    for (int n = 0; n < MAX_TSTEPS; ++n)  // time marching
    {
        // begin n-loop for time marching
        int count_Md_state  = 0;                // used to find how many Md-state in each iteration
        int count_Mc_state  = 0;
        int count_C_state  = 0;
        int count_B_state   = 0;
        int count_SR_state = 0;
        genrand(randNum, N_RU, &state); // fills array with random numbers
        genrand(rand_dATP, N_RU, &state); // fills array with random numbers
        current_time = n*DT;
        calcium = lin_interp_ca(current_time);
        kCa_plus = kCa_plus_ref * calcium;
        //-----------------------------------
        // call the updated RUs
        //-----------------------------------
        
        // Ktr protocol
	//if (n == 3000001)
	//	{
	//	for (int y = 0; y<N_RU; ++y)
	//		{
	//		RU[y] = 0;
	//		}
	//	}

        
        
        if(n==0)
        {
        	f = 0;
        }
        // NOTE: This has been commented out, because I believe that this was the cause of the max_repeats issue. 
        // Force is eventually normalized when it saved, however, at this point, with all of the repeats 
        // running simultaneously, the Force array is inflated when there are more repeats running. 
        // Instead, we are getting the previous fraction of force states from the filament via the 
        // code below after counting the states (f = forceValue;)


        // else
        // {
        // 	f = (float)Force[n-1];
        // }
        // float current_max=0.0;
        // if (current_max < f){
        //     current_max = f;
        //     printf("New_max = %f, %i\n",f, cc);
        // }
        
        float k_plus_SR_ATP = k_plus_SR_ref_ATP; //*(1+k_force_ATP*f);
        float k_plus_SR_dATP = k_plus_SR_ref_dATP; //*(1+k_force_dATP*f);
        float k_minus_SR = k_minus_SR_ref;
        //printf("%f\n",k_plus_SR);
        //printf("%f\n",k_minus_SR);

       update_RUs(lambda, DT, kCa_plus, kCa_minus, randNum, rand_dATP, RU, caRU, kB_plus, kB_minus, k2_plus_dATP, k2_plus_ATP, k2_minus, k3_plus_dATP, k3_plus_ATP, k3_minus, k4_plus_dATP, k4_plus_ATP, k4_minus, percent_dATP, k_force_dATP, k_force_ATP, k_plus_SR_dATP, k_plus_SR_ATP, k_minus_SR,f);

        //--------------------------------------------
        // Obtain Force estimate based on the M-state
        //--------------------------------------------
        for(int i = 0; i < N_RU; ++i)
        {
            if(RU[i]==5) // this represents M2
            {
                ++count_Md_state;
            }
            else if(RU[i]==4) // this represents M1
            {
                ++count_Mc_state;
            }
            else if(RU[i]==3) // this represents C
            {
                ++count_C_state;
            }
            else if (RU[i]==2) // this represents B
            {
                ++count_B_state;
            }
            else if (RU[i]==1) // this represents C* (SRX)
            {
                ++count_SR_state;
            }
            else if (RU[i]==0) // this represents B* (SRX)
            {
                ++count_SR_state;
            }
        }
        float forceValue = (float)count_Md_state / (N_RU); // Type casting because count_Md_state is defined as an int 
        float McValue = (float)count_Mc_state / (N_RU);
        float CValue = (float)count_C_state / (N_RU);
        float BValue = (float)count_B_state / (N_RU);
        float SRValue = (float)count_SR_state / (N_RU);
        
        f = forceValue;
        // float current_max = 0;
        // // This is to look at what happens after the there is an instance where there is at least one state in the force producing state 
        // // It also looks at the following state to see if everything transitions out. 

        // if (flag == 1){
        //     // % ['Count', 'f_prev','SR_prev','f','SRValue']
        //     printf("%i, %f, %f, %f, %f\n",n, f_prev, SR_prev, f, SRValue);
        //     // printf("current SRValue = %f, %i\n", f, n);
        // }

        // if (f > 0) {
        //     flag = 1;
        //     // printf("current f = %f, %i\n",f, n);
        //     // printf("current SRValue = %f, %i\n", f, n);
        //     f_prev = f;
        //     SR_prev = SRValue;

        // }
        // else {
        //     flag = 0;
        // }

        if ( n % 10000 == 0){
            printf("current ca = %f, %f\n",calcium, current_time);
        }

        


        atomicAdd(&(Force[n]), forceValue); // add results every repeat
        atomicAdd(&(Mc[n]), McValue); // add results every repeat
        atomicAdd(&(C[n]), CValue); // add results every repeat
        atomicAdd(&(B[n]), BValue); // add results every repeat
        atomicAdd(&(SR[n]), SRValue); // add results every repeat
    } // end the (n-loop) of the time marching
}
