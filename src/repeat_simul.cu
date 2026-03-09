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
float k_plus_SS,
float k_minus_SS,
float k_plus_alt, // Backdoor pathway
float k_minus_alt, // Backdoor pathway
float K_D,
float coop_N,
float * M1,
float * M2,
float * M3,
float * C,
float * B,
float * SR,
float * SS,
float * ATPase,
int cc,
float protocol, 
float Calc_conc_exp
)

{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    /* initialize random number generation per thread */
    float randNum[N_RU];
    // float rand_drug[N_RU];
    int RU[N_RU];
    bool caRU[N_RU];
    bool drugboundRU[N_RU];
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
    memset(drugboundRU, 0, sizeof(bool)*N_RU);
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
        int count_SS_state = 0;
        ATPcounter = 0;
        genrand(randNum, N_RU, &state); // fills array with random numbers
        // genrand(rand_drug, N_RU, &state); // fills array with random numbers


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
 

        update_RUs(lambda, DT, kCa_plus, kCa_minus, randNum,
            //  rand_drug, 
             RU, caRU, drugboundRU, kB_plus, kB_minus, k1_plus_drug, k1_plus_baseline, k1_minus, k2_plus_drug, k2_plus_baseline, k2_minus, k3_plus_drug, k3_plus_baseline, k3_minus, k4_plus_drug, k4_plus_baseline, k4_minus, percent_drug, k_force_drug, k_force_baseline, k_plus_SR_drug, k_plus_SR_baseline, k_minus_SR, k_plus_SS, k_minus_SS, k_plus_alt, k_minus_alt, K_D, coop_N, f, &ATPcounter);



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
            else if(RU[i]==7) // this represents B**
            {
                ++count_SS_state;
            }
            else if(RU[i]==8) // this represents C**
            {
                ++count_SS_state;
            }
        }
        float forceValue = (float)count_M3_state / (N_RU); // Type casting because count_M3_state is defined as an int 
        float M1Value = (float)count_M1_state / (N_RU);
        float M2Value = (float)count_M2_state / (N_RU);
        float CValue = (float)count_C_state / (N_RU);
        float BValue = (float)count_B_state / (N_RU);
        float SRValue = (float)count_SR_state / (N_RU);
        float SSValue = (float)count_SS_state / (N_RU);

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
        atomicAdd(&(SS[n]), SSValue); // add results every repeat
        atomicAdd(&(ATPase[n]), ATPcounter); // add results every repeat
    } // end the (n-loop) of the time marching
}
