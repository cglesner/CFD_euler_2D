//
//  muscl.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/1/15.
//
//

#ifndef __Euler_2D_FVM__muscl__
#define __Euler_2D_FVM__muscl__

#include <stdio.h>
#include <vector>
#include <array>

using namespace std;

/* This header file contains the functions needed to calculate the flux
 * limiters and the muscl interpolation scheme needed to create a second
 * order accurate scheme. These functions will have to be used twice, 
 * once for the xi-fluxes and again for the eta-fluxes.
 */
void muscl_interpolate(array<double,4>& primVar_L,          \
					   array<double,4>& primVar_R,          \
					   const int& order,                 \
					   const double& kappa,                 \
					   const array<double,4>& pos_lim_L2L1, \
					   const array<double,4>& neg_lim_L1R1, \
					   const array<double,4>& pos_lim_L1R1, \
					   const array<double,4>& neg_lim_R1R2, \
					   const array<double,4>& primVar_L2,   \
					   const array<double,4>& primVar_L1,   \
					   const array<double,4>& primVar_R1,   \
					   const array<double,4>& primVar_R2);
/* This function will be called everytime flux needs to be calculated, and will
 * produce the appropriate left and right states for calculating the flux. It 
 * will take the appropriate flux limiters an primitive variable values at each
 * cell, not the entire vector of vectors of these quantities.
 *
 * Indexing scheme:
 *
 * ... | L2 | L1 |U_l,U_r| R1 | R2 | ...
 *
 */


void calculate_limiters(const char& pos_neg,               \
						array<double,4>& limiters,         \
						const array<double,4>& primVar_L1, \
						const array<double,4>& primVar_R1, \
						const array<double,4>& primVar_L2R2);
/* This function calcualtes the limiters at a given location, taking the flag
 * 'pos_neg' to indicate wheather the positive or negative flux limiters should
 * be calculated at that location. Implementation of the min-mod flux limiter. 
 * It takes the right and left value, the sign that the flux limiter should 
 * have, and the third value, which should correspond to the R2 value in the
 * case of the positive flux_limiter, or L2 in the case of the negative flux 
 * limiter.
 */


#endif /* defined(__Euler_2D_FVM__muscl__) */
