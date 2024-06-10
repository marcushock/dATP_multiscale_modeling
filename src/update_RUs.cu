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
//--------------------------------------------------------------------------%
#include "update_RUs.h"
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

        //-----------------------------------------------------------------
        // if (state = [B* = 0]): Then   [B1*]           else stay as [B0*]
        //     & caState = 0             ^  
        //                               | 
        //                             [B0*]---->[C0*]
        //				 \
        //				  ---->[B0]
        //----------------------------------------------------------------
        if ((state == 0) && (caState == 0))
        {
            p1 = kCa_plus*dt;
            p2 = p1 + lambda*kB_plus[x*N_S+y]*dt;
            if (rand_dATP[i] <= percent_dATP)
            {
                p3 = p2 + k_plus_SR_dATP*(1+k_force_dATP*f)*dt;
            }
            else
            {
                p3 = p2 + k_plus_SR_ATP*(1+k_force_ATP*f)*dt;
            }
            
            if (randNum[i] < p1)
            {
                caRU[i] = 1;   // switch [B0*---->B1*]
            }
            else if(randNum[i] < p2)
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
        //     & caState = 1		 |	
        // 		               [B1*]---->[C1*]
        //     			        |
        //			        v
        //			      [B0*]
        //-----------------------------------------------------------------
        else if ((state == 0) && (caState == 1))
        {
            p1 = kCa_minus*dt;
            p2 = p1 + kB_plus[x*N_S+y]*dt;
            if (rand_dATP[i] <= percent_dATP)
            {
                p3 = p2 + k_plus_SR_dATP*(1+k_force_dATP*f)*dt;
            }
            else
            {
                p3 = p2 + k_plus_SR_ATP*(1+k_force_ATP*f)*dt;
            }
            
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
        //			         | 
        //	             [B0*]<----[C0*]
        //				\
        //				 ---->[C0]
        //-----------------------------------------------------------------
        else if ((state == 1) && (caState == 0))
        {
            p1 = kCa_plus*dt;
            p2 = p1 + kB_minus[x*N_S+y]*dt;
            if (rand_dATP[i] <= percent_dATP)
            {
                p3 = p2 + k_plus_SR_dATP*(1+k_force_dATP*f)*dt;
            }
            else
            {
                p3 = p2 + k_plus_SR_ATP*(1+k_force_ATP*f)*dt;
            }

            if  (randNum[i] < p1)
            {
                caRU[i] = 1;   // switch [C0*---->C1*]
            }
            else if (randNum[i] < p2)
            {
                RU[i] = 0; // switch [C0*---->B0*]
            }
            else if (randNum[i] < p3)
            {
                RU[i] = 3; // switch [C0*---->C0]
            }
        }
        //-----------------------------------------------------------------
        // if (state = [C* = 1]):  Then     ---->[C1]   else stay as [C1*]
        //     & caState = 1               /
        //                    [B1*]<----[C1*]          
        //                               |
        //		                 v
        //                             [C0*]
        //-----------------------------------------------------------------
        else if ((state == 1) && (caState ==1))
        {
            p1 = lambda*kCa_minus*dt;
            p2 = p1 + kB_minus[x*N_S+y]*dt;
            if (rand_dATP[i] <= percent_dATP)
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
        //		      	          /
        //			[B0*]<----
        //----------------------------------------------------------------
        if ((state == 2) && (caState == 0))
        {
            p1 = kCa_plus*dt;
            p2 = p1 + lambda*kB_plus[x*N_S+y]*dt;
            p3 = p2 + k_minus_SR*dt;

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
        //     & caState = 1		        |	
        // 		                      [B1]---->[C1]
        //     			               |
        //			               v
        //			             [B0]
        //-----------------------------------------------------------------
        else if ((state == 2) && (caState == 1))
        {
            p1 = kCa_minus*dt;
            p2 = p1 + kB_plus[x*N_S+y]*dt;
            p3 = p2 + k_minus_SR*dt;

            if (randNum[i] < p1)
            {
                caRU[i] = 0; // switch [B1---->B0]
            }
            else if(randNum[i] < p2)
            {
                RU[i] = 3; //switch [B1---->C1]
            }
            else if(randNum[i] < p3)
            {
                RU[i] = 0; //switch [B1---->B1*]
            }

        }
        //-----------------------------------------------------------------
        // if (state = [C = 3]):  Then     [C1]           else stay as [C0]
        //     & caState = 0                ^  ---->[M1,0]
        //			            | /
        //	                 [B0]<----[C0]
        //			           / \
        //	                 [C0*]<----   ---->[M2,0]
        //-----------------------------------------------------------------
        else if ((state == 3) && (caState == 0))
        {
            p1 = kCa_plus*dt;
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
        //				           v  ---->[M2,1]
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
        //     & caState = 0		     ^
        //     		                     |	
        // 		                   [M1,0]
        //     			          / |
        //		         [C0]<----  v
        //			          [M2,0]
        //-----------------------------------------------------------------
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
        //     & caState = 1		  / |   \
        // 		         [C1]<----  v    ---->[M1,0]  
        //     			          [M2,1]          
        //-----------------------------------------------------------------
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
