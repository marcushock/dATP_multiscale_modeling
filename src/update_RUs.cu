//--------------------------------------------------------------------------------------//
//                    |                                       |                         //
//                    |          Function Name                |                         //
//                    |           update_RUs()                |                         //
//                    |                                       |                         //
//--------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------
// This function is used to update the states of each RUs based on the Markov step
//---------------------------------------------------------------------------------
// Input|
//-------
// S :     Cell array of RU's states
// X :     Before this current state
// Y :     After  this current state
// Ca:     Calcium Concentration = constant value in this simulation
//--------
// Output
//---------
// RUs:     updated RUs according to Neighboring states (X,Y)
// Generally, these can be considered the 5 macrostates. 
// B* - 0 Good
// C* - 1 Good 
// B  - 2 Good
// C  - 3 Good 
// M1 - 4 
// M2 - 5
//--------------------------------------------------------------------------%
#include "update_RUs.h"
#include <stdio.h>



// arr[x][y] == arr[x * row_len + y]
__device__ void update_RUs(float lambda,
                            float dt,
                            float kCa_plus,
                            float kCa_minus,
                            float randNum[N_RU],
                            float rand_dATP[N_RU],
                            int   RU[N_RU],
                            bool caRU[N_RU],
                            float * kB_plus,
                            float * kB_minus,
                            float * k2_plus_dATP,
                            float * k2_plus_ATP,
                            float * k2_minus,
                            float k3_plus_dATP,
                            float k3_plus_ATP,
                            float k3_minus,
                            float * k4_plus_dATP,
                            float * k4_plus_ATP,
                            float * k4_minus,
                            float percent_dATP,
                            float k_force_dATP,
                            float k_force_ATP, 
                            float k_plus_SR_dATP,
                            float k_plus_SR_ATP, 
                            float k_minus_SR,
                            float f
                            )

{

    int state, x, y;
    bool caState;
    float p1, p2, p3, p4, p5;

    for (int i=1; i < N_RU-1; i++)   // only the interior RUs
    //for (int i = 0; i < N_RU; ++i) 
    {
        state = RU[i];     // get the current RU  state ID number i.e., B* = 0, C* = 1,  B = 2, C = 3, M1 = 4, M2 = 5
        caState = caRU[i];      // get current calcium status (1 = calcium present, 0 = calcium not present)
        x     = RU[i-1];
        y     = RU[i+1];

        // printf("%d %f %f\n", i, rand_dATP[i], randNum[i]);
        //-----------------------------------------------------------------
        // if (state = [B* = 0]): Then   [B1*]           else stay as [B0*]
        //     & caState = 0             ^  
        //                               | 
        //                             [B0*]---->[C0*]
        //				                 |
        //				                 ---->[B0]
        // So this is starting in state B0*
        // We have no calcium bound caState = 0
        //----------------------------------------------------------------
        if ((state == 0) && (caState == 0))
        {
            p1 = kCa_plus*dt; // This calculates the transition transition probability into calcium bound state...
            // kCa_plus appears to be calculated by multiplying the calcium concentration by the kCa_plus_ref value in the repeat_simul.cu function. 
            
            p2 = p1 + lambda*kB_plus[x*N_S+y]*dt;  // Unclear as to why we multiply by lambda, which I beleive is set to 0 
            // Per Abby, is due to a prevention of calcium from unbinding (in a general case). 
            // Here, it appears that you cannot move into a C0 state ever... 

            // This chunk of code is used to calculate kinetics based on either ATP parameters or dATP parameters 
            if (rand_dATP[i] <= percent_dATP) // percent dATP is somewhere between 0 and 1, which then helps to identify if we have ATP kinetics or dATP kinetics 
            {
                p3 = p2 + k_plus_SR_dATP*(1+k_force_dATP*f)*dt; // We calculate a new probability p3 using dATP kinetic parameters 
            }
            else
            {
                p3 = p2 + k_plus_SR_ATP*(1+k_force_ATP*f)*dt; // Otherwise use the basal ATP kinetic parameters to calcualte p3
            }
            // P3 is used for calculating whether we transition into or out of the SRX. 
            
            // Determine stochastically whether after this timestep if Ca will be bound to the thin finalment. 
            // Based on the fact that lambda is 0, it appears that we can only transition from B0* to B1*
            // In summary, the blocked (calcium free) and SRX (off) state can only transition to the 
            // Calcium bound confirmation of the SRX. 
            // if ((p3 < 1) || (p2 > 1) || (p1 > 1))
            // {
            //     // printf("Step i: %d. dATP rand %f, p1: %f, p2: %f, p3: %f\n", i, rand_dATP[i], p1, p2, p3);
            //     //asm("trap;"); // Force the kernel to terminate immediately
            // }
            if (randNum[i] < p1)
            {
                caRU[i] = 1;   // switch [B0*---->B1*]'
            }
            else if(randNum[i] < p2) // This checks if the 
            {
                RU[i] = 1; //switch [B0*---->C0*]
            }
            else if(randNum[i] < p3)
            {
            	RU[i] = 2; //switch [B0*---->B0]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [B* = 0]): Then	  ---->[B1]      else stay as [B1*]
        //     & caState = 1	          |	
        // 		                        [B1*]---->[C1*]
        //     			                  |
        //			                      v
        //			                    [B0*]
        //-----------------------------------------------------------------
        else if ((state == 0) && (caState == 1))
        {
            p1 = kCa_minus*dt; // Calculate unbinding probability of calcium 
            p2 = p1 + kB_plus[x*N_S+y]*dt; // Calculate the transition probability from B1* to C1* which is the unblocking of the thin filament. 
            if (rand_dATP[i] <= percent_dATP)
            {
                p3 = p2 + k_plus_SR_dATP*(1+k_force_dATP*f)*dt; // Calculate the probability of transitioning out of the SRX/OFF state
            }
            else
            {
                p3 = p2 + k_plus_SR_ATP*(1+k_force_ATP*f)*dt; //Calculate rate of out OFF state but instead assuming ATP kinetics. 
            }
            
            // Check if we have a state change in one of the 3 possible transitions based on above calculated probabilities. 
            if (randNum[i] < p1)
            {
                caRU[i] = 0; // switch [B1*---->B0*]
            }
            else if(randNum[i] < p2)
            {
                RU[i] = 1; //switch [B1*---->C1*]
            }
            else if(randNum[i] < p3)
            {
                RU[i] = 2; //switch [B1*---->B1]
            }

        }
        //-----------------------------------------------------------------
        // if (state = [C* = 1]): Then  [C1*]            else stay as [C0*]
        //     & caState = 0             ^  
        //		             	         | 
        //	                 [B0*]<----[C0*]
        //			                   	 |
        //				                 ---->[C0]
        //-----------------------------------------------------------------
        // Starting in the C* state (without calcium bound )
        else if ((state == 1) && (caState == 0))
        {
            p1 = kCa_plus*dt; // Calculate the probability for calcium binding. 
            p2 = p1 + kB_minus[x*N_S+y]*dt; // Calculate probability of moving back into the blocked state 
            if (rand_dATP[i] <= percent_dATP)
            {
                p3 = p2 + k_plus_SR_dATP*(1+k_force_dATP*f)*dt; // Calcualte probability out of the SRX/OFF state (dATP)
            }
            else
            {
                p3 = p2 + k_plus_SR_ATP*(1+k_force_ATP*f)*dt; // Calcualte probability out of the SRX/OFF state (dATP)
            }

            if  (randNum[i] < p1) // Check to see if moving into calcium bound state 
            {
                caRU[i] = 1;   // switch [C0*---->C1*]
            }
            else if (randNum[i] < p2) // Check to see if moving back into the blocked state 
            {
                RU[i] = 0; // switch [C0*---->B0*]
            }
            else if (randNum[i] < p3) // Check to see if moving out of the SRX state 
            {
                RU[i] = 3; // switch [C0*---->C0]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [C* = 1]):  Then     ---->[C1]   else stay as [C1*]
        //     & caState = 1               /
        //                    [B1*]<----[C1*]          
        //                               |
        //		                         v
        //                             [C0*]
        //-----------------------------------------------------------------
        // Starting with calcium bound, in the C* (SRX) state 
        // P1 is forced to be 0 because Lambda is also zero. Therefore we never move from the 
        // C* state with calcium bound back to the C* state without calcium bound. This appears to assume that ca cannot unbind when 
        // in the closed state. It must first transition back to the blocked state before calcium can unbind. 
        else if ((state == 1) && (caState ==1))
        {
            p1 = lambda*kCa_minus*dt; // Once again Lambda is included and still set to 0. 
            p2 = p1 + kB_minus[x*N_S+y]*dt; // Calculate prob of going back to the blocked state 
            if (rand_dATP[i] <= percent_dATP) // Calculate prob of going out of the SRX (again dATP dependent)
            {
                p3 = p2 + k_plus_SR_dATP*(1+k_force_dATP*f)*dt;
            }
            else
            {
                p3 = p2 + k_plus_SR_ATP*(1+k_force_ATP*f)*dt;
            }
  

            if  (randNum[i] < p1)
            {
                caRU[i] = 0; // switch [C1*---->C0*]
            }
            else if (randNum[i] < p2)
            {
                RU[i] = 0; // switch [C1*---->B1*]
            }
            else if (randNum[i] < p3)
            {
                RU[i] = 3; // switch [C1*---->C1]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [B = 2]):  Then    [B1]            else stay as [B0]
        //     & caState = 0               ^  
        //                                 | 
        //                               [B0]---->[C0]
        //		      	                  /
        //			            [B0*]<----
        //----------------------------------------------------------------
        // In the ON state (not SRX/OFF) and calcium unbound. 
        // Once again lambda is included in these calculations which appears to prevent certain transitions. 
        // Specifically, this appears to prevent the transition from B0 to C0. 
        if ((state == 2) && (caState == 0))
        {
            p1 = kCa_plus*dt;
            p2 = p1 + lambda*kB_plus[x*N_S+y]*dt;
            p3 = p2 + k_minus_SR*dt; // Calculate the probility of moving back into the SRX/OFF state 

            if (randNum[i] < p1)
            {
                caRU[i] = 1;   // switch [B0---->B1]
            }
            else if(randNum[i] < p2)
            {
                RU[i] = 3; //switch [B0---->C0]
            }
            else if(randNum[i] < p3)
            {
            	RU[i] = 0; //switch [B0---->B0*]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [B = 2]): Then  [B1*]<----	  else stay as [B1]
        //     & caState = 1		            |	
        // 		                               [B1]---->[C1]
        //     			                        |
        //			                            v
        //			                          [B0]
        //-----------------------------------------------------------------
        // Starting in the B state, ON state, with calcium bound 
        else if ((state == 2) && (caState == 1))
        {
            p1 = kCa_minus*dt; // Calculate the probility of calcium unbinding. 
            p2 = p1 + kB_plus[x*N_S+y]*dt; // Calculate probability of moving into the close state from blocked state 
            p3 = p2 + k_minus_SR*dt; // Calculate the probability of moving 

            if (randNum[i] < p1)
            {
                caRU[i] = 0; // switch [B1---->B0] calcium unbinding 
            }
            else if(randNum[i] < p2)
            {
                RU[i] = 3; //switch [B1---->C1] moving into the closed state
            }
            else if(randNum[i] < p3)
            {
                RU[i] = 0; //switch [B1---->B1*] moving back into the blocked SRX state
            }

        }
        //-----------------------------------------------------------------
        // if (state = [C = 3]):  Then     [C1]           else stay as [C0]
        //     & caState = 0                ^  ---->[M1,0]
        //			                        | /
        //	                               [B0]<----[C0]
        //			                       / \
        //	                     [C0*]<----   ---->[M2,0]
        //-----------------------------------------------------------------
        // With Lambda set to 0, I don't beleive that we can ever get to this state
        // Closed state of myosin, but without calcium bound... 
        // However, it should be noted, that with the cooperative parameters, maybe there is a way to 
        // Move into one of these states with the cooperativity alone.
        // Above diagrams seems to be wrong...

        else if ((state == 3) && (caState == 0))
        {
            p1 = kCa_plus*dt; // Calculate probability of calcium binding. 
            p2 = p1 + kB_minus[x*N_S+y]*dt;
	        p3 = p2 + k_minus_SR*dt;
	        if (rand_dATP[i] <= percent_dATP)
            {
                p4 = p3 + k2_plus_dATP[x*N_S+y]*dt;
            }
            else
            {
                p4 = p3 + k2_plus_ATP[x*N_S+y]*dt;
            }
            p5 = p4 + k4_minus[x*N_S+y]*dt;

            if  (randNum[i] < p1)
            {
                caRU[i] = 1;   // switch [C0---->C1]
            }
            else if (randNum[i] < p2)
            {
                RU[i] = 2; // switch [C0---->B0]
            }
            else if (randNum[i] < p3)
            {
                RU[i] = 1; // switch [C0---->C0*]
            }
            else if (randNum[i] < p4)
            {
                RU[i] = 4; // switch [C0---->M1,0]
            }
            else if (randNum[i] < p5)
            {
                RU[i] = 5; // switch [C0---->M2,0]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [C = 3]):   Then                   else stay as [C1]
        //     & caState = 1	         [C1*]<----    ---->[M1,1]
        //                                         \  /  
        //                               [B1]<----[C1]          
        //                                         | \
        //				                           v  ---->[M2,1]
        //                                       [C0]
        //-----------------------------------------------------------------
        else if ((state == 3) && (caState ==1))
        {
            p1 = lambda*kCa_minus*dt;
            p2 = p1 + kB_minus[x*N_S+y]*dt;
	        p3 = p2 + k_minus_SR*dt;
	        if (rand_dATP[i] <= percent_dATP)
            {
                p4 = p3 + k2_plus_dATP[x*N_S+y]*dt;
            }
            else
            {
                p4 = p3 + k2_plus_ATP[x*N_S+y]*dt;
            }
            p5 = p4 + k4_minus[x*N_S+y]*dt;
  

            if  (randNum[i] < p1)
            {
                caRU[i] = 0; // switch [C1---->C0]
            }
            else if (randNum[i] < p2)
            {
                RU[i] = 2; // switch [C1---->B1]
            }
            else if (randNum[i] < p3)
            {
                RU[i] = 1; // switch [C1---->C1*]
            }
            else if (randNum[i] < p4)
            {
                RU[i] = 4; // switch [C1---->M1,1]
            }
            else if (randNum[i] < p5)
            {
                RU[i] = 5; // switch [C1---->M2,1]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [M1 = 4]): Then     [M1,1]       else stay as [M1,0]
        //     & caState = 0		       ^
        //     		                       |	
        // 		                        [M1,0]
        //     			                 / |
        //		                [C0]<----  v
        //			                     [M2,0]
        //-----------------------------------------------------------------
        // This is in state M1, which is the weakly bound state. 
        // There is no calcium bound, so it's still surprising to be in this state to be honest. 
        else if ((state == 4) && (caState == 0))
        {
            p1 = kCa_plus*dt;
            if (rand_dATP[i] <= percent_dATP)
            {
                p2 = p3 + k3_plus_dATP*dt;
            }
            else
            {
                p2 = p3 + k3_plus_ATP*dt;
            }
            p3 = p2 + k2_minus[x*N_S+y]*dt;

            if  (randNum[i] < p1)
            {
                caRU[i] = 1;   // switch [M1,0---->M1,1]
            }
            else if (randNum[i] < p2)
            {
                RU[i] = 5; // switch [M1,0---->M2,0]
            }
            else if (randNum[i] < p3)
            {
                caRU[i] = 3; // switch [M1,0---->C0]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [M1 = 4]): Then     [M1,1]       else stay as [M1,1]
        //     & caState = 1	           / |   \
        // 		                  [C1]<----  v    ---->[M1,0]  
        //     			                   [M2,1]          
        //-----------------------------------------------------------------
        // This is in state M1, which is the weakly bound state. 
        // Now we do have calcium bound, so it's more likely that we do visit this state. 
        else if ((state == 4) && (caState == 1))
        {
            p1 = lambda*kCa_minus*dt;
            if (rand_dATP[i] <= percent_dATP)
            {
                p2 = p1 + k3_plus_dATP*dt;
            }
            else
            {
                p2 = p1 + k3_plus_ATP*dt;
            }
            p3 = p2 + k2_minus[x*N_S+y]*dt;

            if  (randNum[i] < p1)
            {
                caRU[i] = 0;   // switch [M1,1---->M1,0]
            }
            else if (randNum[i] < p2)
            {
                RU[i] = 5; // switch [M1,1---->M2,1]
            }
            else if (randNum[i] < p3)
            {
                caRU[i] = 3; // switch [M1,1---->C1]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [M2 = 5]): Then     [M1,0]    	else stay as [M2,0]
        //     & caState = 0      [C0]<----   ^   ---->[M2,1]       	
        // 		                   \  |  / 
        //			           [M2,0]
        //-----------------------------------------------------------------
        // Now in the strongly bound state. 
        // Again no calcium is bound. Not sure how we could ever visit this state. 
        else if ((state == 5) && (caState == 0))
        {
            p1 = kCa_plus*dt;
            if (rand_dATP[i] <= percent_dATP)
            {
                p2 = p1 + k4_plus_dATP[x*N_S+y]*dt;
            }
            else
            {
                p2 = p1 + k4_plus_ATP[x*N_S+y]*dt;
            }
            p3 = p2 + k3_minus*dt;

            if  (randNum[i] < p1)
            {
                caRU[i] = 1;   // switch [M2,0---->M2,1]
            }
            else if (randNum[i] < p2)
            {
                RU[i] = 3; // switch [M2,0---->C0]
            }
            else if (randNum[i] < p3)
            {
                caRU[i] = 4; // switch [M2,0---->M1,0]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [M2 = 5]): Then    [M1,1]  	else stay as [M2,1]
        //     & caState = 1     [C1]<----   ^       	
        // 		                  \  |      
        //			          [M2,1]
        //     			             \  
        //			             v
        //			          [M2,0]
        //-----------------------------------------------------------------
        // In the M2 force producing state. 
        // Calcium is actually bound, but is not able to unbind. 
        else if ((state == 5) && (caState == 1))
        {
            p1 = lambda*kCa_minus*dt;
            if (rand_dATP[i] <= percent_dATP)
            {
                p2 = p1 + k4_plus_dATP[x*N_S+y]*dt;
            }
            else
            {
                p2 = p1 + k4_plus_ATP[x*N_S+y]*dt;
            }
            p3 = p2 + k3_minus*dt;

            if  (randNum[i] < p1)
            {
                caRU[i] = 0;   // [M2,1---->M2,0]
            }
            else if (randNum[i] < p2)
            {
                RU[i] = 3; // switch [M2,1---->C1]
            }
            else if (randNum[i] < p3)
            {
                caRU[i] = 4; // switch [M2,1---->M1,1]
            }
        }
        ///-----------------------------------------------------------------

    } // close the for loop

} // close function
