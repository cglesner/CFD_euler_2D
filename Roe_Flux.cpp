//
//  Roe_Flux.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/1/15.
//


// Standard Includes
#include <cmath>
#include <array>


// My Includes
#include "Roe_Flux.h"
#include "Euler_Helper_Functions.h"

using namespace std;


void Roe_Flux(array<double,4>& Fn,        \
              const array<double,4>& V_L, \
              const array<double,4>& V_R, \
              const array<double,2>& n)
{
	// Calculate centered flux from left and right states
    array<double,4> c_flux;
	centered_flux(c_flux, V_L, V_R, n);

	// Calculate wave flux term
    array<double,4> w_flux;
	wave_flux(w_flux, V_L, V_R, n);
	
    // Combine flux components.
	for (int i = 0; i < 4; ++i)
	{
		Fn[i] = c_flux[i] - w_flux[i];
	}
}


void centered_flux(array<double,4>& c_flux,    \
                   const array<double,4>& V_L, \
                   const array<double,4>& V_R, \
                   const array<double,2>& n)
{
	// read in variables (for readability)
	double dens_L = V_L[0];
	double xvel_L = V_L[1];
    double yvel_L = V_L[2];
	double pres_L = V_L[3];
	
	double dens_R = V_R[0];
	double xvel_R = V_R[1];
    double yvel_R = V_R[2];
	double pres_R = V_R[3];
    
    double nx = n[0];
    double ny = n[1];
    
    // Calculate normal velocity
    double Un_L = xvel_L*nx + yvel_L*ny;
    double Un_R = xvel_R*nx + yvel_R*ny;
	
	// Calculate enthalpy
	double ht_L = total_enthalpy(V_L);
	double ht_R = total_enthalpy(V_R);
	
	// Calculate output
	c_flux[0] = 0.5*(dens_L*Un_L + dens_R*Un_R);
	c_flux[1] = 0.5*(dens_L*xvel_L*Un_L + nx*pres_L + dens_R*xvel_R*Un_R + nx*pres_R);
    c_flux[2] = 0.5*(dens_L*yvel_L*Un_L + ny*pres_L + dens_R*yvel_R*Un_R + ny*pres_R);
	c_flux[3] = 0.5*(dens_L*Un_L*ht_L + dens_R*Un_R*ht_R);
}

void wave_flux(array<double,4>& w_flux,    \
               const array<double,4>& V_L, \
               const array<double,4>& V_R, \
               const array<double,2>& n)
{
	// Calculate roe averaged variables
	array<double,5> roe_vars;
	roe_avg(roe_vars,V_L,V_R);
	
	// Calculate eigenvalues and vectors
	array<double,4> lambda;
    array<array<double,4>,4> rightEig;
	roe_eigval(lambda,roe_vars,n);
	roe_eigvec(rightEig,roe_vars,n);
	
	// Calculate wave amplitudes
	array<double,4> dw;
	wave_amp(dw,V_L,V_R,roe_vars,n);
	
	// calculate wave flux
	for (int i=0; i < 4; ++i)
	{
        w_flux[i] = 0.0; // reset value at memory address.
        for (int j=0; j < 4; ++j)
        {
            w_flux[i] += 0.5*abs(lambda[j])*dw[j]*rightEig[i][j];
        }
	}
}

void roe_avg(array<double,5>& roe_vars,  \
             const array<double,4>& V_L, \
             const array<double,4>& V_R)
{
    const double GAMMA = 1.4;
    
	// Calculate enthalpy
	double ht_L = total_enthalpy(V_L);
	double ht_R = total_enthalpy(V_R);
	
	// calculate R for efficient calculation of roe variables.
    double R = sqrt(V_R[0]/V_L[0]);
    double vel_mag;
    double a2;
	
	// calculate outputs
	roe_vars[0]  = R*V_L[0];                           // density
	roe_vars[1]  = (R*V_R[1] + V_L[1])/(R+1.0);        // x-velocity
    roe_vars[2]  = (R*V_R[2] + V_L[2])/(R+1.0);        // y-velocity
	roe_vars[3]  = (R*ht_R + ht_L)/(R+1.0);            // enthalpy

    vel_mag      = (roe_vars[1]*roe_vars[1] + roe_vars[2]*roe_vars[2]);
    a2           = (GAMMA - 1.0)*(roe_vars[3] - 0.5*vel_mag);
    roe_vars[4]  = sqrt(a2);                     // speed sound
}


void roe_eigval(array<double,4>& lambda,         \
                const array<double,5>& roe_vars, \
                const array<double,2>& n)
{
    // read in vectors for readability
    double roe_xvel = roe_vars[1];
    double roe_yvel = roe_vars[2];
    double roe_spsd = roe_vars[4];
    
    double nx = n[0];
    double ny = n[1];
    
    // calculate normal velocity
    double Un = roe_xvel*nx + roe_yvel*ny;
    
	// use Harten's method to control for vanishing eigenvalues.
	lambda[0] = harten(Un,roe_spsd);
    lambda[1] = lambda[0];
	lambda[2] = harten( (Un+roe_spsd) ,roe_spsd);
	lambda[3] = harten( (Un-roe_spsd) ,roe_spsd);
}


void roe_eigvec(array<array<double,4>,4>& rightEig, \
                const array<double,5>& roe_vars,    \
                const array<double,2>& n)
{
    // read in vectors for readability
    double roe_dens = roe_vars[0];
    double roe_xvel = roe_vars[1];
    double roe_yvel = roe_vars[2];
    double roe_ht   = roe_vars[3];
    double roe_spsd = roe_vars[4];

    double nx = n[0];
    double ny = n[1];

    // calculate normal velocity
    double Un = roe_xvel*nx + roe_yvel*ny;
    
    // pre-calculate
	double F = (roe_dens/(2.0*roe_spsd));
	
	// Calculate outputs
    rightEig[0][0] = 1.0;
    rightEig[1][0] = roe_xvel;
    rightEig[2][0] = roe_yvel;
    rightEig[3][0] = 0.5*(roe_xvel*roe_xvel + roe_yvel*roe_yvel);
    
    rightEig[0][1] = 0.0;
    rightEig[1][1] = ny*roe_dens;
    rightEig[2][1] = -nx*roe_dens;
    rightEig[3][1] = roe_dens*(ny*roe_xvel - nx*roe_yvel);
    
    rightEig[0][2] = F;
    rightEig[1][2] = F*(roe_xvel + nx*roe_spsd);
    rightEig[2][2] = F*(roe_yvel + ny*roe_spsd);
    rightEig[3][2] = F*(roe_ht + Un*roe_spsd);
    
    rightEig[0][3] = -F;
    rightEig[1][3] = -F*(roe_xvel - nx*roe_spsd);
    rightEig[2][3] = -F*(roe_yvel - ny*roe_spsd);
    rightEig[3][3] = -F*(roe_ht - Un*roe_spsd);
}


double harten(double lambda, double roe_ss)
{
	// set parameter used by Harten's method
	double e = 0.1;
	
	// calculate output
	if (abs(lambda) < 2.0*e*roe_ss)
	{
		double lambda_mag = lambda*lambda/(4.0*e*roe_ss) + e*roe_ss;
		return lambda*lambda_mag/abs(lambda);
	}
	else
	{
		return lambda;
	}
}


void wave_amp(array<double,4>& dw,             \
              const array<double,4>& V_L,      \
              const array<double,4>& V_R,      \
              const array<double,5>& roe_vars, \
              const array<double,2>& n)
{
	// Read in vectors for readability
    double roe_dens = roe_vars[0];
    double roe_spsd = roe_vars[4];
    
    double nx = n[0];
    double ny = n[1];
    
    // define the difference variables
	double d_dens = V_R[0] - V_L[0];
	double d_xvel = V_R[1] - V_L[1];
    double d_yvel = V_R[2] - V_L[2];
	double d_pres = V_R[3] - V_L[3];
	
	// Calculate outputs
	dw[0] = d_dens - d_pres/(roe_spsd*roe_spsd);
    dw[1] = ny*d_xvel - nx*d_yvel;
    dw[2] = nx*d_xvel + ny*d_yvel + d_pres/(roe_dens*roe_spsd);
    dw[3] = nx*d_xvel + ny*d_yvel - d_pres/(roe_dens*roe_spsd);
}



