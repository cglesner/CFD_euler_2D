//
//  Time_Step.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/4/15.
//  
//

#ifndef __Euler_2D_FVM__Time_Step__
#define __Euler_2D_FVM__Time_Step__

#include <stdio.h>
#include <vector>

using namespace std;

/* This header contains all of the functions used by the Time_Step funciton.
 * Time_Step takes the following inputs:
 * -> the array of primitive variables (passed by reference for efficiency)
 * -> the volumes and face areas of the mesh
 * -> the proscribed CFL number
 *
 * Using these, Time_Step calculates the maximum allowable timestep at each cell.
 * In this way, time-step makes all timesteps availiable to main if they are 
 * desired. another simple funciton will be able to choose a global timestep if needed.
 * Time_Step does not need to be told the size of the arrays passed to it
 * as long as the arrays are passed using the 'vector' container.
 */

void Time_Step(vector<vector<double> >& dt,                            \
			   const double& CFL,                                       \
			   const vector<vector<double> >& cell_Volume,              \
			   const vector<vector<double> >& xi_Area,                  \
			   const vector<vector<double> >& eta_Area,                 \
			   const vector<vector<array<double,2> > >& xi_Normal,      \
			   const vector<vector<array<double,2> > >& eta_Normal,     \
			   const vector<vector<array<double,4> > >& primitive_vars);
/* Top level function which will operate on the entire mesh at once, calculating the
 * maximum allowable local time step for each cell.
 */


double max_lambda(const array<double,4>& primitive_vars, \
				  const array<double,2>& Normal_avg);
/* This function will calcuate the maximum eigenvalue at a given cell, taking the vector
 * of primitive variables and the averaged normal vector as arguments.
 */

#endif /* defined(__Euler_2D_FVM__Time_Step__) */


