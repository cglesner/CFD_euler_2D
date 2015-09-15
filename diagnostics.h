//
//  diagnostics.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 7/18/15.
//

#ifndef __Euler_2D_FVM__diagnostics__
#define __Euler_2D_FVM__diagnostics__

#include <stdio.h>
#include <array>
#include <vector>

using namespace std;

/* This headerfile contains all of the functions that will be used to perform
 * the diagnostics required by this project, including inlet total pressure
 * loss, airfoil pressure, lift, and drag coefficients, as well as discretiation
 * error.
 */


void mms_diagnostics(const string sim_name, const int iter, const int& CI, const int& CJ, const int& sub_super, \
					 const vector<vector<array<double,2> > >& cell_ctr, const vector<vector<array<double,2> > >& node_loc, \
					 const vector<vector<array<double,4> > >& primVar);
/* This function will perform the diagnostics called for in the MMS case.
 * Once the simulation has converged this diagnostic function will be called, 
 * and will take the converged solution, calculate the L1 and Linf norms of the
 * discretization error and output it to the terminal, as well as outputting
 * a tecplot file showing the discretization error in the converged solution.
 */


void inlet_diagnostics(const int& CI, const int& CJ,                    \
						const array<double,4>& V_FS,                     \
						const vector<vector<array<double,4> > >& primVar);
/* This function will perform the diagnostics needed for the inlet case.
 * It will take the primitive variables by reference, as well as the number
 * of cells in the X and Y directions, and the free-stream primitive
 * variables. This fuction will not return anything by refference. It will
 * call the global output and terminal output files internally, and no other
 * function in the simulation needs to operate on the results of the 
 * diagnostics so it is not necissary to delcare variables in main to hold 
 * them, or return any variables from this function.
 */


void airfoil_diagnostics(const int& CI, const int& CJ, const int& edge,    \
						 const double& AoA, const string& sim_name,        \
						 const vector<vector<array<double,2> > >& node_loc,\
						 const vector<vector<array<double,2> > >& eta_norm, \
						 const vector<vector<double> >& eta_A,              \
						 const array<double,4>& V_FS,                      \
						 const vector<vector<array<double,4> > >& primVar);
/* This function will perform the diagnostics needed for the airfoil case.
 * Once the simulation has converged this diagnostic function will determine
 * the pressure coefficient along the surface of the airfoil, as well as the
 * lift coefficient and drag coefficients, outputting them to a diagnostics
 * file.
 */


#endif /* defined(__Euler_2D_FVM__diagnostics__) */
