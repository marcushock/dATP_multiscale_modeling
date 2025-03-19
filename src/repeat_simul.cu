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
float * k1_plus_drug,
float * k1_plus_baseline,
float * k1_minus,
float k2_plus_drug,
float k2_plus_baseline,
float k2_minus,
float k3_plus_drug,
float k3_plus_baseline,
float k3_minus,
float * k4_plus_drug,
float * k4_plus_baseline,
float * k4_minus,
float * kB_plus,
float * kB_minus,
float kCa_plus_ref,
float kCa_minus_ref,
float percent_drug,
float k_force_drug,
float k_force_baseline,
float k_plus_SR_drug,
float k_plus_SR_baseline,
float k_minus_SR,
float * M1,
float * M2,
float * M3,
float * C,
float * B,
float * SR,
float * ATPase,
int cc,
float protocol, 
float Calc_conc_exp
)

{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    /* initialize random number generation per thread */
    float randNum[N_RU];
    float rand_drug[N_RU];
    int RU[N_RU];
    bool caRU[N_RU];
    curandState_t state;
    float f;
    int flag = 0;
    float f_prev;
    float SR_prev;
    float kCa_plus;
    float kCa_minus;
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
    float ATPcounter;
    for (int n = 0; n < MAX_TSTEPS; ++n)  // time marching
    {
        // begin n-loop for time marching
        int count_M1_state  = 0;
        int count_M2_state  = 0;
        int count_M3_state  = 0;                // used to find how many M3-state in each iteration
        int count_C_state  = 0;
        int count_B_state   = 0;
        int count_SR_state = 0;
        ATPcounter = 0;
        genrand(randNum, N_RU, &state); // fills array with random numbers
        genrand(rand_drug, N_RU, &state); // fills array with random numbers


        if (protocol == 1){
            kCa_plus = Calc_conc_exp*kCa_plus_ref; 
        }
        else if  (protocol == 0){
            current_time = n*DT;
            calcium = lin_interp_ca(current_time);
            kCa_plus = kCa_plus_ref * calcium;
        }

        kCa_minus = kCa_minus_ref;
    // float kCa_plus     = Cal_conc*kCa_plus_ref;
        // current_time = n*DT;
        // calcium = lin_interp_ca(current_time);
        // kCa_plus = kCa_plus_ref * calcium;
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
        
        // float k_plus_SR_ATP = k_plus_SR_baseline; //*(1+k_force_ATP*f); DELETE THIS
        // float k_plus_SR_drug = k_plus_SR_ref_drug; //*(1+k_force_drug*f); DELETE THIS LINE
        // float k_minus_SR = k_minus_SR_ref; // DELETE THIS LINE 
        //printf("%f\n",k_plus_SR);
        //printf("%f\n",k_minus_SR);
        update_RUs(lambda, DT, kCa_plus, kCa_minus, randNum, rand_drug, RU, caRU, kB_plus, kB_minus, k1_plus_drug, k1_plus_baseline, k1_minus, k2_plus_drug, k2_plus_baseline, k2_minus, k3_plus_drug, k3_plus_baseline, k3_minus, k4_plus_drug, k4_plus_baseline, k4_minus, percent_drug, k_force_drug, k_force_baseline, k_plus_SR_drug, k_plus_SR_baseline, k_minus_SR,f, &ATPcounter);
        // Print out all the kinetic variables that start with the letter k 
        // printf("k1_plus_drug = %f\n", k1_plus_drug[0]);
        // printf("k1_plus_baseline = %f\n", k1_plus_baseline[0]);
        // printf("k1_minus = %f\n", k1_minus[0]);
        // printf("k2_plus_drug = %f\n", k2_plus_drug);
        // printf("k2_plus_baseline = %f\n", k2_plus_baseline);
        // printf("k2_minus = %f\n", k2_minus);
        // printf("k3_plus_drug = %f\n", k3_plus_drug);
        // printf("k3_plus_baseline = %f\n", k3_plus_baseline);
        // printf("k3_minus = %f\n", k3_minus);
        // printf("k4_plus_drug = %f\n", k4_plus_drug[0]);
        // printf("k4_plus_baseline = %f\n", k4_plus_baseline[0]);
        // printf("k4_minus = %f\n", k4_minus[0]);
        // printf("kB_plus = %f\n", kB_plus[0]);


        //--------------------------------------------
        // Obtain Force estimate based on the M-state
        //--------------------------------------------
        for(int i = 0; i < N_RU; ++i)
        {
            if (RU[i]==0) // this represents B* (SRX)
            {
                ++count_SR_state;
            }
            else if (RU[i]==1) // this represents C* (SRX)
            {
                ++count_SR_state;
            }
            else if (RU[i]==2) // this represents B
            {
                ++count_B_state;
            }
            else if(RU[i]==3) // this represents C
            {
                ++count_C_state;
            }
            else if(RU[i]==4) // this represents M1
            {
                ++count_M1_state;
            }
            else if (RU[i]==5){ // this represents M2
                ++count_M2_state;
            }   
            else if(RU[i]==6) // this represents M3
            {
                ++count_M3_state;
            }
        }
        float forceValue = (float)count_M3_state / (N_RU); // Type casting because count_M3_state is defined as an int 
        float M1Value = (float)count_M1_state / (N_RU);
        float M2Value = (float)count_M2_state / (N_RU);
        float CValue = (float)count_C_state / (N_RU);
        float BValue = (float)count_B_state / (N_RU);
        float SRValue = (float)count_SR_state / (N_RU);
        
        f =  (float)count_M3_state + (float)count_M2_state; // Could also include some function of the M2 value here 
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

        // if ( n % 10000 == 0){
        //     printf("current ca = %f, %f\n",calcium, current_time);
        // }

        


        atomicAdd(&(M1[n]), M1Value); // add results every repeat
        atomicAdd(&(M2[n]), M2Value); // add results every repeat
        atomicAdd(&(M3[n]), forceValue); // add results every repeat
        atomicAdd(&(C[n]), CValue); // add results every repeat
        atomicAdd(&(B[n]), BValue); // add results every repeat
        atomicAdd(&(SR[n]), SRValue); // add results every repeat
        atomicAdd(&(ATPase[n]), ATPcounter); // add results every repeat
    } // end the (n-loop) of the time marching
}
