//
//  boundary_conditions.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/5/15.
//
//

#ifndef __Euler_2D_FVM__boundary_conditions__
#define __Euler_2D_FVM__boundary_conditions__

#include <stdio.h>
#include <array>

using namespace std;

/* This header file contains all the functions used to calculate
 * boundary conditions in this simulation. Boundary conditions are
 * imposed in this simulation through the use of ghost cells. These
 * functions are all defined to create two cells to support second
 * order accuracy. The types of boundary conditions defined here are:
 * (1.) Free Stream (2.) Slip Wall (3.) Periodic (4.) Mach (5.) MMS
 */

void freeStream_bc(array<double,4>& V_GC1,       \
                   array<double,4>& V_GC2,       \
                   const array<double,4>& V_FS );
/* This boundary condition is used when we are assuming that at the
 * edge of our domain the perturbation of the flow is sufficiently
 * small to assume that the ghost cells remain at a constant,
 * predetermined set of primitive values.
 */

void slipWall_bc(array<double,4>& V_GC1,         \
                 array<double,4>& V_GC2,         \
                 const array<double,4>& V_IC1,   \
                 const array<double,4>& V_IC2,   \
                 const array<double,2>& wall_n);
/* This boundary condition is used for flow next to a frictionless wall.
 * it sets the values of the ghost cells s.t. the velocity of the flow
 * at the boundary will be only tangent to the wall.
 */

void periodic_bc(array<double,4>& V_GC1,          \
                 array<double,4>& V_GC2,          \
                 const array<double,4>& V_per1,   \
                 const array<double,4>& V_per2);
/* This boundary condition is used for enforcing a periodic boundary 
 * condition. It will coppy the values of V_per1 and V_per2 into 
 * V_GC1 and V_GC2 respectively.
 */


void super_outflow_bc(array<double,4>& V_GC1,         \
                      array<double,4>& V_GC2,         \
                      const array<double,4>& V_IC1,   \
                      const array<double,4>& V_IC2);
/* This boundary condition is used for calculating a supersoninc
 * outflow boundary. It will simply extrapolate the interior values
 * out to the ghost cells.
 */


void super_inflow_bc(array<double,4>& V_GC1,         \
                     array<double,4>& V_GC2,         \
                     const array<double,4>& V_IF);
/* This boundary condition is used for supersonic inflow conditions. 
 * it will simply set inflow conditions.
 */


#endif /* defined(__Euler_2D_FVM__boundary_conditions__) */
