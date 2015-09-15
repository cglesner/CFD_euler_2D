//
//  Van_Leer_Flux.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 4/24/15.
//


// Standard Includes
#include <cmath>
#include <array>

// my includes
#include "Van_Leer_Flux.h"
#include "Euler_Helper_Functions.h"


using namespace std;


void Van_Leer_Flux(array<double,4>& Fn,        \
				   const array<double,4>& V_L, \
				   const array<double,4>& V_R, \
				   const array<double,2>& n)
{

    // read in vectors to variable names for readabillity
    double dens_L = V_L[0];
    double dens_R = V_R[0];
    
    double xvel_L = V_L[1];
    double xvel_R = V_R[1];
    
    double yvel_L = V_L[2];
    double yvel_R = V_R[2];
    
    double pres_L = V_L[3];
    double pres_R = V_R[3];
    
    double nx = n[0];
    double ny = n[1];
    
    // Calculate normal velocity
    double Un_L = xvel_L*nx + yvel_L*ny;
    double Un_R = xvel_R*nx + yvel_R*ny;
    
    // calculate the speed of sound, mach, and total enthalpy for left and right states
    double a_L = sound_speed(V_L);
    double a_R = sound_speed(V_R);
    
    double mach_L = Un_L/a_L;
    double mach_R = Un_R/a_R;
    
    double ht_L = total_enthalpy(V_L);
    double ht_R = total_enthalpy(V_R);
    
    // Calculate the multipliers for the convective and pressure flux terms
    double c_pos = dens_L*a_L*C(mach_L,'+');
    double c_neg = dens_R*a_R*C(mach_R,'-');
    
    double p_pos = D(mach_L,'+');
    double p_neg = D(mach_R,'-');
    
    //Calcuate total flux
    Fn[0] = c_pos + c_neg;
    Fn[1] = (c_pos*xvel_L + c_neg*xvel_R) + (p_pos*nx*pres_L + p_neg*nx*pres_R);
    Fn[2] = (c_pos*yvel_L + c_neg*yvel_R) + (p_pos*ny*pres_L + p_neg*ny*pres_R);
    Fn[3] = c_pos*ht_L + c_neg*ht_R;
};


double C(double mach, char pos_neg)
{
    return ALPHA(mach,pos_neg)*(1.0 + BETA(mach))*mach - BETA(mach)*M_pm(mach,pos_neg);
}

double ALPHA(double mach, char pos_neg)
{
    double flip;

    if (pos_neg == '+')
    {
        flip = 1.0;
    }
    else
    {
        flip = -1.0;
    }
    
    return 0.5*(1.0 + flip * SIGN(1.0,mach));
}

double BETA(double mach)
{
	double x = 1.0 - INT(abs(mach));
    
    if (x < 0)
    {
        return 0.0;
    }
    else
    {
        return -1.0*x;
    }
}

double M_pm(double mach, char pos_neg)
{
    double flip = 0.0;
    
    if (pos_neg == '+')
    {
        flip = 1.0;
    }
    else
    {
        flip = -1.0;
    }
    
    double mpf = mach + flip;
    
    return 0.25*flip*mpf*mpf;
}

double D(double mach, char pos_neg)
{
    return ALPHA(mach,pos_neg)*(1.0+BETA(mach))-BETA(mach)*P_bar(mach,pos_neg);
}

double P_bar(double mach, char pos_neg)
{
    double flip;
    if (pos_neg == '+')
    {
        flip = 1.0;
    }
    else
    {
        flip = -1.0;
    }
    return M_pm(mach,pos_neg)*(-1.0*mach + 2.0*flip);
}

double SIGN(double A, double B)
{
    double signB;
    if (B == 0.0)
    {
        signB = 1.0;
    }
    else
    {
		signB = B/abs(B);
    }
	return signB * abs(A);
}

double INT(double A)
{
    if (A >= 0.0)
    {
		return floor(A);
    }
    else
    {
		return ceil(A);
    }
}





