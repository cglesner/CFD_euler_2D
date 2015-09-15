//
//  Geometry.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/4/15.
//
//

#ifndef __Euler_2D_FVM__Geometry__
#define __Euler_2D_FVM__Geometry__

#include <stdio.h>
#include <vector>
#include <array>

using namespace std;

/* This header contains all the functions that will be called by the
 * Geometry funciton. Calculate_Geometry will take the array 
 * of grid node coordinates and determine all geometric
 * quantities needed by the simulation:
 * -> Face areas
 * -> Face normals
 * -> Cell volumes
 * -> Cell centers
 *
 * These functions are needed for Calcuate_Geometry only. All functions
 * assume width into the page = 1.
 */

void Set_Geometry(vector< vector<double> >& xi_Area,               \
			  vector< vector<double> >& eta_Area,                  \
			  vector< vector< array<double,2> > >& xi_Normal,      \
			  vector< vector< array<double,2> > >& eta_Normal,     \
			  vector< vector<double> >& cell_Volume,               \
			  vector< vector< array<double,2> > >& cell_Center,    \
			  const vector< vector< array<double,2> > >& node_loc, \
			  const int& NI, const int& NJ);
/* Top level funciton which will take locations of the nodes of the grid in the form
 * of an [ni][nj][2] array. it will calculate the value of all relevant geometric
 * quantities. The vectors which contain these quantities must be initialized outside
 * of the Geometry funciton.
 */


void area(double& A,                     \
		  const array<double,2>& node_1, \
		  const array<double,2>& node_2);
/* This funciton will return the area of an arbitrary cell
 * face bracketed by node_1 and node_2. The order of the nodes
 * is arbitrary.
 */


void normal(array<double,2>& n,            \
			const array<double,2>& node_1, \
	        const array<double,2>& node_2);
/* This funciton will return the x and y normal components
 * of an arbitrary cell face. 
 * NOTE:------------------------------------------------------
 * It is assumed that node_2 is at a higher index than node_1.
 * -----------------------------------------------------------
 */


void volume(double& V,                     \
			const array<double,2>& node_A, \
			const array<double,2>& node_B, \
			const array<double,2>& node_C, \
			const array<double,2>& node_D);
/* This funciton will return the volume of an arbitrary cell.
 * The ordering of the nodes should be: 
 * A -> South West
 * B -> South East
 * C -> North East
 * D -> North West
 */


void center(array<double,2>& cell_center,  \
			const array<double,2>& node_A, \
			const array<double,2>& node_B, \
			const array<double,2>& node_C, \
			const array<double,2>& node_D);
/* This funciton will return the location of the center of an 
 * arbitrary cell. The ordering of the nodes should be:
 * A -> South West
 * B -> South East
 * C -> North East
 * D -> North West
 */



#endif /* defined(__Euler_2D_FVM__Geometry__) */
