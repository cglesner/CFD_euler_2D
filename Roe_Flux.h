//
//  Roe_Flux.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/1/15.
//
//
//

#ifndef __Euler_2D_FVM__Roe_Flux__
#define __Euler_2D_FVM__Roe_Flux__

#include <stdio.h>
#include <array>

using namespace std;

/*
 * This header file includes all funcitons needed for the calculation of
 * flux using Roe's scheeme. These functions are needed for Roe
 * Flux calculation only.
 *
 *
 *  Primitive Variables:      Flux Variables:
 *  V[0] = density            F[0] = mass flux
 *  V[1] = x-velocity         F[1] = x-momentum flux
 *  V[2] = y-velocity         F[2] = y-momentum flux
 *  V[3] = Pressure           F[3] = energy flux
 *
 *  n[0] = x-normal
 *  n[1] = y-normal
 */

void Roe_Flux(array<double,4>& Fn,        \
			  const array<double,4>& V_L, \
			  const array<double,4>& V_R, \
			  const array<double,2>& n);
/* The actual function that calculates flux.
 *  Fn  -> normal flux
 *  V_L -> left primitive values
 *  V_R -> right primitive values
 *  n   -> normals
 */

void centered_flux(array<double,4>& c_flux,    \
				   const array<double,4>& V_L, \
				   const array<double,4>& V_R, \
				   const array<double,2>& n);
/* this function will calculate the centered flux at a cell interface based
 * on the value of primitive variables provided.
 */

void wave_flux(array<double,4>& w_flux,    \
			   const array<double,4>& V_L, \
			   const array<double,4>& V_R, \
			   const array<double,2>& n);
/* This function will take the primitive variables to the left and right of
 * a cell interface and return the flux due to the waves in the system, which is
 * how Roe's method determines the upwind contribution to the flux.
 */

void roe_avg(array<double,5>& roe_vars,  \
			 const array<double,4>& V_L, \
			 const array<double,4>& V_R);
/* This function will take the primitive variables to the left and right of
 * a cell interface and return the roe averaged variables for the interface.
 * this function will assume that no negative pressures or densities occur as 
 * this should be flagged in the top level function which calls this funciton.
 */

void roe_eigval(array<double,4>& lambda,         \
				const array<double,5>& roe_vars, \
				const array<double,2>& n);
/* This function will return the value of the eigenvalues as calculated by 
 * roe's method.
 */

void roe_eigvec(array<array<double,4>,4>& rightEig, \
				const array<double,5>& roe_vars,          \
				const array<double,2>& n);
/* This function will return the right eigen vectors at a given matrix,
 * taking the roe averaged values of the flow variables at the interface.
 * The three outputs will each be a [1X3] array.
 */

double harten(double lamdba, double roe_ss);
/* This function will modify the eigenvalues calculated in Roe's method to
 * prevent eigenvalues from approaching zero, in the method developed by
 * Harten.
 */

void wave_amp(array<double,4>& dw,             \
			  const array<double,4>& V_L,      \
			  const array<double,4>& V_R,      \
			  const array<double,5>& roe_vars, \
			  const array<double,2>& n);
/* This function will calculate the amplitude of the waves present in the
 * 2-d euler system.
 */


#endif /* defined(__Euler_2D_FVM__Roe_Flux__) */

