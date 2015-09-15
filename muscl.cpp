//
//  muscl.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/1/15.
//  
//

// standard includes
#include <vector>
#include <array>
#include <cmath>

// my includes
#include "muscl.h"
#include "Euler_Helper_Functions.h"

using namespace std;

void muscl_interpolate(array<double,4>& primVar_L,          \
                       array<double,4>& primVar_R,          \
                       const int& order,                    \
                       const double& kappa,                 \
                       const array<double,4>& pos_lim_L2L1, \
                       const array<double,4>& neg_lim_L1R1, \
                       const array<double,4>& pos_lim_L1R1, \
                       const array<double,4>& neg_lim_R1R2, \
                       const array<double,4>& primVar_L2,   \
                       const array<double,4>& primVar_L1,   \
                       const array<double,4>& primVar_R1,   \
                       const array<double,4>& primVar_R2)
{
    // declare variables for brevity
    double central;
    double upwind_l;
    double upwind_r;
	
	// declare conserved variables
	array<double,4> consVar_L2;
	array<double,4> consVar_L1;
	array<double,4> consVar_R1;
	array<double,4> consVar_R2;
	array<double,4> consVar_L;
	array<double,4> consVar_R;
    
    // loop over equations
    for (int n = 0; n < 4; ++n)
    {
		// convert primitive to conserved
		primitive_to_conserved( consVar_L2, primVar_L2);
		primitive_to_conserved( consVar_L1, primVar_L1);
		primitive_to_conserved( consVar_R1, primVar_R1);
		primitive_to_conserved( consVar_R2, primVar_R2);
		
		// define central and upwind components
        central  = consVar_R1[n] - consVar_L1[n];
        upwind_l = consVar_L1[n] - consVar_L2[n];
		upwind_r = consVar_R2[n] - consVar_R1[n];
        consVar_L[n] = consVar_L1[n] + \
                        0.25*(order-1)*((1.0-kappa)*pos_lim_L2L1[n]*upwind_l + \
                                    (1.0+kappa)*neg_lim_L1R1[n]*central);
		if (consVar_L[n] <= 0.0)
		{
			consVar_L[n] = consVar_L1[n];
		}
        consVar_R[n] = consVar_R1[n] - \
                        0.25*(order-1)*((1.0-kappa)*neg_lim_R1R2[n]*upwind_r + \
                                    (1.0+kappa)*neg_lim_L1R1[n]*central);
		if (consVar_R[n] <= 0.0)
		{
			consVar_R[n] = consVar_R1[n];
		}
		
		// convert conserved back to primitive
		conserved_to_primitive( primVar_L, consVar_L);
		conserved_to_primitive( primVar_R, consVar_R);
    }
}


void calculate_limiters(const char& pos_neg,               \
						array<double,4>& limiters,         \
						const array<double,4>& primVar_L1, \
						const array<double,4>& primVar_R1, \
						const array<double,4>& primVar_L2R2)
{
    // declare variables
    double numerator;
    double denomenator;
    double epsilon = 1.0e-12;
	double signed_epsilon;
    double R;
    double min;
 
    // loop over equations
    for (int n = 0; n < 4; ++n)
    {
        if (pos_neg == '+')
        {
            numerator = primVar_L2R2[n] - primVar_R1[n];
        }
        else
        {
            numerator = primVar_L1[n] - primVar_L2R2[n];
        }
        
        // ensure that denomenator is not zero:
        denomenator = primVar_R1[n] - primVar_L1[n];
        signed_epsilon = epsilon * denomenator/abs(denomenator);
        denomenator = (abs(denomenator) < epsilon) ? signed_epsilon : denomenator;
		
		R = numerator/denomenator;
        // calculate min-mod flux limiter.
        // Syntax is: limiter = max( 0, min( R, 1 ) )
        min = (R < 1.0) ? R : 1.0;
        limiters[n] = (0.0 < min) ? min : 0.0;
    }
}










