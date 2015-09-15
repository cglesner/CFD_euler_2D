//
//  main.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 4/30/15.
//
//

#ifndef __Euler_2D_FVM__main__
#define __Euler_2D_FVM__main__

// standard includes
#include <stdio.h>
#include <array>
#include <vector>

using namespace std;

/* This header file contains all of the functions called just in main.
 * grouping of functions in other header files for specific types of tasks
 * or in this one is admittedly somewhat arbitrary.
 */
 
 
int main();
/* the main funciton that will be the core of our executable.
 */


void resize_vars(const int& NI, const int& NJ,                           \
				 vector<vector<array<double,4> > >& xi_pos_limiters,     \
				 vector<vector<array<double,4> > >& xi_neg_limiters,     \
				 vector<vector<array<double,4> > >& eta_pos_limiters,    \
				 vector<vector<array<double,4> > >& eta_neg_limiters,    \
				 vector<vector<array<double,4> > >& xi_Flux,             \
				 vector<vector<array<double,4> > >& eta_Flux,            \
				 vector<vector<array<double,4> > >& primVar,             \
				 vector<vector<array<double,4> > >& source,              \
				 vector<vector<array<double,4> > >& residual,            \
				 vector<vector<double> >& dt);
/* This funciton will take care of properly sizing all of the arrays of
 * pysical variables in our system. based on the number of nodes in the 
 * i and j direction, NI, and NJ.
 */


void mms_and_source_init(const char& sim_case, const int& sub_super,    \
                    const int& CI, const int& CJ,                       \
                    const vector<vector<array<double,2> > >& cell_ctr,  \
                    vector<vector<array<double,4> > >& primVar,         \
                    vector<vector<array<double,4> > >& source);
/* This function will set the boundary condition for the mms case, if that
 * case is being run, and will set the mms source terms if the mms case is
 * being run. if the mms case is NOT being run, it will set all source 
 * to zero.
 */


void update_BC(const char& sim_case, const int& edge,          \
               const int& CI, const int& CJ,                   \
               const vector<vector<array<double,2> > >& xi_N,  \
               const vector<vector<array<double,2> > >& eta_N, \
               const array<double,4>& inflow_condition,        \
               vector<vector<array<double,4> > >& primVar);
/* This funciton will perform the boundary condition update each iteration
 * of the simulation. It will take a simulation flag to set which of three
 * simulation cases is being performed: 
 * 'm' = method of manufactured soltions, 
 * 'i' = supersonic-inlet, 
 * 'a' = NACA-airfoil.
 */


void update_limiters(vector<vector<array<double,4> > >& xi_pos_limiters,     \
						vector<vector<array<double,4> > >& xi_neg_limiters,  \
						vector<vector<array<double,4> > >& eta_pos_limiters, \
						vector<vector<array<double,4> > >& eta_neg_limiters, \
						const vector<vector<array<double,4> > >& primVar,    \
						const int& CI, const int& CJ);
/* This function will calculate all the limiters that will be needed for the 
 * muscl extrapolation to be carried out in the 'update_fluxes' function.
 */


void update_fluxes(vector<vector<array<double,4> > >& xi_Flux,                \
                   vector<vector<array<double,4> > >& eta_Flux,               \
                   const char& method, const int& order, const double& kappa, \
                   const vector<vector<array<double,4> > >& xi_pos_limiters,  \
                   const vector<vector<array<double,4> > >& xi_neg_limiters,  \
                   const vector<vector<array<double,4> > >& eta_pos_limiters, \
                   const vector<vector<array<double,4> > >& eta_neg_limiters, \
                   const vector<vector<array<double,4> > >& primVar,          \
                   const vector<vector<array<double,2> > >& xi_N,             \
                   const vector<vector<array<double,2> > >& eta_N,            \
				   const int& CI, const int& CJ);
/* This function will update the fluxes at all cell interfaces in the domain,
 * it will overwrite the previous saved fluxes when it executes.
 */


void update_residuals(vector<vector<array<double,4> > >& residual,       \
                      const vector<vector<array<double,4> > >& xi_Flux,  \
					  const vector<vector<array<double,4> > >& eta_FLux, \
					  const vector<vector<array<double,4> > >& source,   \
					  const vector<vector<double> >& cell_V,             \
				      const vector<vector<double> >&  xi_A,              \
				      const vector<vector<double> >& eta_A,              \
				      const int& CI, const int& CJ);
/* This funciton will update the residual at all cells in the domain,
 * it will overwrite the previously saved residuals when it executes.
 *
 * To fit with the convension of the normal directions we have chosen, 
 * we will say tha the lower index faces will be fluxes in
 * while the higher index faces will be fluxes out. We will also
 * consider ther residual to be all of the spacial terms collected
 * on the oposite side of the equals sign as the time derivative,
 * i.e.
 * dU/dt*VOL = RESIDUAL
 * therefore:
 * RESIDUAL = SOURCE + (flux-in) - (flux-out)
 */


void update_solution(vector<vector<array<double,4> > >& primVar,         \
					 const vector<vector<double> >& dt,                  \
					 const vector<vector<double> >& cell_V,              \
					 const vector<vector<array<double,4> > >& residual,  \
					 const int& alpha, const int& CI, const int& CJ);
/* This funciton will update the solution at all cells in the domain.
 * it will overwrite the previously saved primitive variables when it executes.
 */


void residual_norms(array<double,4>& l2_norm, const int& CI, const int& CJ,   \
					const vector<vector<array<double,4> > >& residual);
/* This funciton will calculate l1, l2, and linf of all residual variables. */


#endif /* defined(__Euler_2D_FVM__main__) */
