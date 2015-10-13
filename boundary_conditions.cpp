//
//  boundary_conditions.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/5/15.
//
//


#include <array>

#include "boundary_conditions.h"

void freeStream_bc(array<double,4>& V_GC1,       \
				   array<double,4>& V_GC2,       \
				   const array<double,4>& V_FS )
{
	for (int n = 0; n < 4; ++n)
	{
		V_GC1[n] = V_FS[n];
		V_GC2[n] = V_FS[n];
	}
}


void slipWall_bc(array<double,4>& V_GC1,         \
					 array<double,4>& V_GC2,         \
				    const array<double,4>& V_IC1,   \
				    const array<double,4>& V_IC2,   \
				    const array<double,2>& wall_n )
{
	// read in wall normal variables
	double nx = wall_n[0];
	double ny = wall_n[1];
	
	// calculate terms used in velocity extapolation:
	double A = nx*nx - ny*ny;
	double B = 2.0*nx*ny;

	// calculate extrapolated velocities
	V_GC1[1] = -A*V_IC1[1] - B*V_IC1[2];
	V_GC1[2] = -B*V_IC1[1] + A*V_IC1[2];
	
	V_GC2[1] = -A*V_IC1[1] - B*V_IC1[2];
	V_GC2[2] = -B*V_IC1[1] + A*V_IC1[2];
	
	// calculate extrapolated pressure
	V_GC1[3] = V_IC1[3];
	V_GC2[3] = V_IC1[3];
	
	// calculate extrapolated density
	V_GC1[0] = V_IC1[0];
	V_GC2[0] = V_IC1[0];
}

void periodic_bc(array<double,4>& V_GC1,          \
				 array<double,4>& V_GC2,          \
				 const array<double,4>& V_per1,   \
				 const array<double,4>& V_per2)
{
	for (int n = 0; n < 4; ++n)
	{
		V_GC1[n] = V_per1[n];
		V_GC2[n] = V_per2[n];
	}
}


void super_outflow_bc(array<double,4>& V_GC1,         \
					  array<double,4>& V_GC2,         \
					  const array<double,4>& V_IC1,   \
					  const array<double,4>& V_IC2)
{
	for (int n = 0; n < 4; ++n)
	{
		V_GC1[n] = V_IC1[n];
		V_GC2[n] = V_IC1[n];
	}
}

void super_inflow_bc(array<double,4>& V_GC1,         \
					 array<double,4>& V_GC2,         \
					 const array<double,4>& V_IF)
{
	for (int n = 0; n < 4; ++n)
	{
		V_GC1[n] = V_IF[n];
		V_GC2[n] = V_IF[n];
	}
}

void new_bc()
{
	// new BC
}