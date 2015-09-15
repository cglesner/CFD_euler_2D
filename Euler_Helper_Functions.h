//
//  Euler_Helper_Functions.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 4/24/15.
//

#ifndef __Euler_2D_FVM__Euler_Helper_Functions__
#define __Euler_2D_FVM__Euler_Helper_Functions__

#include <stdio.h>
#include <vector> 
#include <array>

using namespace std;

void input_to_prim(array<double,4>& primVar, double mach, double pres, \
				   double temp, double AoA);
/* Function that will take the flow input variables and convert them into
 * the coresponding primitive variables.
 */


double total_energy(array<double,4> V);
/* Function that takes the vector of primitive variables and
 * returns the total_energy at that location.
 */


double total_enthalpy(array<double,4> V);
/* Function that takes the vector of primitive variables and
 * returns the total_enthalpy at that location.
 */


double sound_speed(array<double,4> V);
/* Function that takes the vector of primitive variables and 
 * returns the speed of sound at that location.
 */


double total_pressure(array<double,4> V);
/* Function that takes the vector of primitive variables and
 * returns the total pressure at this location.
 */


double Cp(array<double,4> Vinf, double pres);
/* Function that takes the vector of freestream primitive variables
 * and the pressure at the current cell and returns the pressure coefficient.
 */


double Cl(array<double,4> Vinf, double L, double c);
/* Function that takes the vector of freestream primitive variables,
 * the lift per unit span and the chord length and returns the lift coefficient.
 */


double Cd(array<double,4> Vinf, double D, double c);
/* Function that takes the vector of freestream primitive variables,
 * the drag per unit span and the chord length and returns the drag coefficient.
 */


void primitive_to_conserved(array<double,4>& U, const array<double,4>& V);
/* This function will take the primitive variables at a location
 * and calculate the value of the equivalent conserved variables.
 */


void conserved_to_primitive(array<double,4>& V, const array<double,4>& U);
/* This function will take the conserved variables at a location
 * and calculate the value of the equivalent primitive variables.
 */


void L1_norm(array<double,4>& norm, const vector< vector< array<double,4> > >& variable_array);
/* This function will take an array and return the L1 norm of that array.
 */


void L2_norm(array<double,4>& norm, const vector< vector< array<double,4> > >& variable_array);
/* This function will take an array and return the L2 norm of that array.
 */


void Linf_norm(array<double,4>& norm, const vector< vector< array<double,4> > >& variable_array);
/* This function will take an array and return the Linf norm of that array.
 */


#endif /* defined(__Euler_2D_FVM__Euler_Helper_Functions__) */
