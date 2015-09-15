//
//  Euler_Helper_Functions.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 4/24/15.
//  Copyright (c) 2015 Colin Christopher Glesner. All rights reserved.
//

#include "Euler_Helper_Functions.h"
#include <cmath>
#include <array>
#include <vector>

using namespace std;

/*  Primitive Variables:
 *  V[0] = density
 *  V[1] = x-velocity
 *  V[2] = y-velocity
 *  V[3] = Pressure
 */


void input_to_prim(array<double,4>& primVar, double mach, double pres, \
				   double temp, double AoA)
{
	// define gas constants
	double GAMMA = 1.4;
	double R = 287.058;
	
	// Convert angle of attack into radians
	AoA = AoA*3.14159265359/180.0;
	
	// calculate density
	primVar[0] = pres/( R * temp );
	
	// determine speed of sound
	double a = sqrt(GAMMA * R * temp);
	
	// calculate Vx and Vy
	primVar[1] = cos(AoA)*mach*a;
	primVar[2] = sin(AoA)*mach*a;
	
	// assign pressure
	primVar[3] = pres;
}


double total_energy(array<double,4> V)
{
	// set gas constant
	const double GAMMA = 1.4;
	
	// read in vector for readability
	double dens = V[0];
	double xvel = V[1];
	double yvel = V[2];
	double pres = V[3];
	
	return pres/(dens*(GAMMA-1.0)) + 0.5*(xvel*xvel + yvel*yvel);
}


double total_enthalpy(array<double,4> V)
{
	// set gas constant
	const double GAMMA = 1.4;
	
	// read in vector for readability
	double dens = V[0];
	double xvel = V[1];
	double yvel = V[2];
	double pres = V[3];
	
	return (GAMMA/(GAMMA-1.0))*(pres/dens) + 0.5*(xvel*xvel + yvel*yvel);
}


double sound_speed(array<double,4> V)
{
	// set gas constant
	const double GAMMA = 1.4;
	
	// read in vector for readability
	double dens = V[0];
	double pres = V[3];
	
	return sqrt(GAMMA*pres/dens);
}


double total_pressure(array<double,4> V)
{
	// set gas constant
	const double GAMMA = 1.4;
	
	// read in vector for readability
	double dens = V[0];
	double xvel = V[1];
	double yvel = V[2];
	double pres = V[3];
	
	// define additional variables:
	double a = sqrt(GAMMA*pres/dens);
	double M = sqrt(xvel*xvel + yvel*yvel)/a;
	
	// return total pressure
	return pres/( pow((1.0 + 0.2*M*M),-3.5) );
}


double Cp(array<double,4> Vinf, double pres)
{
	// read in vector for readability
	double dens_inf = Vinf[0];
	double xvel_inf = Vinf[1];
	double yvel_inf = Vinf[2];
	double pres_inf = Vinf[3];
	
	// velocity magnitude squared
	double Vinf2 = xvel_inf*xvel_inf + yvel_inf*yvel_inf;
	
	// return pressure coefficient
	return (pres - pres_inf)/(0.5*dens_inf*Vinf2);
}


double Cl(array<double,4> Vinf, double L, double c)
{
	// read in vector for readability
	double dens_inf = Vinf[0];
	double xvel_inf = Vinf[1];
	double yvel_inf = Vinf[2];
	
	// velocity magnitude squared
	double Vinf2 = xvel_inf*xvel_inf + yvel_inf*yvel_inf;
	
	// return lift coefficient
	return L/(0.5*dens_inf*Vinf2*c);
}


double Cd(array<double,4> Vinf, double D, double c)
{
	// read in vector for readability
	double dens_inf = Vinf[0];
	double xvel_inf = Vinf[1];
	double yvel_inf = Vinf[2];
	
	// velocity magnitude squared
	double Vinf2 = xvel_inf*xvel_inf + yvel_inf*yvel_inf;
	
	// return drag coefficient
	return D/(0.5*dens_inf*Vinf2*c);
}


void primitive_to_conserved(array<double,4>& U, const array<double,4>& V)
{
	// read in vector for readability
	double dens = V[0];
	double xvel = V[1];
	double yvel = V[2];
	
	U[0] = dens;                 // mass
	U[1] = dens*xvel;            // x-momentum
	U[2] = dens*yvel;            // y-momentum
	U[3] = dens*total_energy(V); // energy
}


void conserved_to_primitive(array<double,4>& V, const array<double,4>& U)
{
	// set gas constant
	const double GAMMA = 1.4;
	
	// read in vector for readability
	double mass = U[0];
	double xmmt = U[1];
	double ymmt = U[2];
	double enrg = U[3];
	
	V[0] = mass;                                                   // density
	V[1] = xmmt/mass;                                              // x-velocity
	V[2] = ymmt/mass;                                              // y-velocity
	V[3] = (GAMMA -1.0)*(enrg - 0.5*(xmmt*xmmt + ymmt*ymmt)/mass); // pressure
}


void L1_norm(array<double,4>& norm, const vector< vector< array<double,4> > >& variable_array)
{
	// determine size of array.
	int CI = variable_array.size();
	int CJ = variable_array[0].size();
	array<double,4> array_sum = {{0.0, 0.0, 0.0, 0.0}};
	
	// sum over each element in the array
	for (int i = 0; i < CI; ++i)
	{
		for (int j = 0; j < CJ; ++j)
		{
			for (int n = 0; n < 4; ++n)
			{
				array_sum[n] += abs(variable_array[i][j][n]);
			}
		}
	}
	
	// normalize sum by the number of elements
	double num_elem = CI*CJ;

	for (int n = 0; n < 4; ++n)
	{
		norm[n] = array_sum[n]/num_elem;
	}
}


void L2_norm(array<double,4>& norm, const vector< vector< array<double,4> > >& variable_array)
{
	// determine size of the array
	int CI = variable_array.size();
	int CJ = variable_array[0].size();
	
	// claculate l2_norm of the residuals
	array<double,4> sumsquares = {{0.0, 0.0, 0.0, 0.0}};
	
	for (int i = 0; i < CI; ++i)
	{
		for (int j = 0; j < CJ; ++j)
		{
			for (int n = 0; n < 4; ++n)
			{
				sumsquares[n] += variable_array[i][j][n]*variable_array[i][j][n];
			}
		}
	}
	
	// determine the number of elements in the array
	double num_elem = CI*CJ;
	
	for (int n = 0; n < 4; ++n)
	{
		norm[n] = sqrt(sumsquares[n]/num_elem);
	}
}


void Linf_norm(array<double,4>& norm, const vector< vector< array<double,4> > >& variable_array)
{
	// determine size of the array
	int CI = variable_array.size();
	int CJ = variable_array[0].size();

	// define placeholder array
	array<double,4> array_max = {{0.0, 0.0, 0.0, 0.0}};
	
	
	// run through each element in the array
	for (int i = 0; i < CI; ++i)
	{
		for (int j = 0; j < CJ; ++j)
		{
			for (int n = 0; n < 4; ++n)
			{
				if (abs(variable_array[i][j][n]) > array_max[n])
				{
					array_max[n] = abs(variable_array[i][j][n]);
				}
			}
		}
	}
	
	// load array max into norm
	for (int n = 0; n < 4; ++n)
	{
		norm[n] = array_max[n];
	}
}

