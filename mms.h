//
//  mms.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/12/15.
//
//

#ifndef __Euler_2D_FVM__mms__
#define __Euler_2D_FVM__mms__

#include <array>
#include <vector>
#include <stdio.h>

using namespace std;

void MMS_initial_and_boundary(int sub_super, vector<vector<array<double,4> > >& primVar,     \
						   const vector<vector<array<double,2> > >& cell_ctr, \
						   const int& CI, const int& CJ);
/* set the initial condition for the MMS case.
 */

void MMS_exactsolution(int sub_super, array<double,4>& primVar,         \
					            const array<double,2>& cell_ctr);
/* This function returns the exact solution for the Euler MMS function used
 * in this code verification.
 */

void MMS_source(int sub_super, array<double,4>& source,         \
						 const array<double,2>& cell_ctr);
/* This function returns the analytical source terms for MMS when provided 
 * with a given location.
 */

void MMS_domain_and_ghost(vector<vector<array<double,2> > >& cell_and_ghost_ctr,\
						 const vector<vector<array<double,2> > >& cell_ctr,    \
						 const int& CI, const int& CJ);
/* This function takes the cell centers and extrapolates out ghost cell 
 * center locations.
 */



#endif /* defined(__Euler_2D_FVM__mms__) */
