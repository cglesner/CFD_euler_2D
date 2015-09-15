//
//  Geometry.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/4/15.
//  
//


// standard includes
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <array>

// my includes
#include "Geometry.h"


using namespace std;


void Set_Geometry(vector< vector<double> >& xi_Area,               \
			  vector< vector<double> >& eta_Area,                  \
			  vector< vector< array<double,2> > >& xi_Normal,      \
			  vector< vector< array<double,2> > >& eta_Normal,     \
			  vector< vector<double> >& cell_Volume,               \
			  vector< vector< array<double,2> > >& cell_Center,    \
			  const vector< vector< array<double,2> > >& node_loc, \
			  const int& NI, const int& NJ)
{
	// resize output vectors:
	/*********************************/
	xi_Area.resize(NI);
	xi_Normal.resize(NI);
	for	( int i = 0; i < NI; ++i)
	{
		xi_Area[i].resize(NJ-1);
		xi_Normal[i].resize(NJ-1);
	}
	
	eta_Area.resize(NI-1);
	eta_Normal.resize(NI-1);
	for	( int i = 0; i < NI-1; ++i)
	{
		eta_Area[i].resize(NJ);
		eta_Normal[i].resize(NJ);
	}
	
	cell_Volume.resize(NI-1);
	cell_Center.resize(NI-1);
	for	( int i = 0; i < NI-1; ++i)
	{
		cell_Volume[i].resize(NJ-1);
		cell_Center[i].resize(NJ-1);
	}
	/**********************************/
	
	// calculate xi-face areas and normals
	for(int i = 0; i < NI; ++i)
	{
		for (int j = 0; j < (NJ-1); ++j)
		{
			area( xi_Area[i][j],       \
				 node_loc[i][j],       \
				 node_loc[i][j+1]);
			normal(xi_Normal[i][j],    \
				   node_loc[i][j],     \
				   node_loc[i][j+1]);
		}
	}
	
	// calculate eta-face areas and normals
	for (int i = 0; i < (NI-1); ++i)
	{
		for (int j = 0; j < NJ; ++j)
		{
		    // node_Loc ordered for positive normal.
			area(eta_Area[i][j],       \
				 node_loc[i][j],       \
				 node_loc[i+1][j] );
			
			// node_Loc ordered for positive normal.
			normal(eta_Normal[i][j],   \
				   node_loc[i+1][j],   \
				   node_loc[i][j] );
		}
	}
	
	// calculate cell volumes and centers
	for (int j = 0; j < (NJ-1); ++j)
	{
		for (int i = 0; i < (NI-1); ++i)
		{
			volume(cell_Volume[i][j],  \
				   node_loc[i][j],     \
				   node_loc[i+1][j],   \
				   node_loc[i+1][j+1], \
				   node_loc[i][j+1]);

			center(cell_Center[i][j],  \
				   node_loc[i][j],     \
				   node_loc[i+1][j],   \
				   node_loc[i+1][j+1], \
				   node_loc[i][j+1]);
		}
	}
}


void area(double& A,                     \
		  const array<double,2>& node_1, \
		  const array<double,2>& node_2)
{
	// read in values from vectors for readability
	double x1 = node_1[0];
    double x2 = node_2[0];
	
    double y1 = node_1[1];
    double y2 = node_2[1];
	
	// calculate x-difference, y difference
    double dx = x2 - x1;
	double dy = y2 - y1;
	
	// calculate output
	A = sqrt(dx*dx + dy*dy);
}


void normal(array<double,2>& n,            \
			const array<double,2>& node_1, \
			const array<double,2>& node_2)
{
	// read in values from vectors for readability
	double x1 = node_1[0];
	double x2 = node_2[0];
	
	double y1 = node_1[1];
	double y2 = node_2[1];
	
	// calculate x-difference, y difference
	double dx = x2 - x1;
	double dy = y2 - y1;

	// Calculate face area
	double A;
	area(A,node_1,node_2);
	
	// calculate x-component of normal
	double nx = dy/A;
	double ny = -dx/A;
	
	// read calculated values into output vector.
	n[0] = nx;
	n[1] = ny;
}


void volume(double& V,                     \
			const array<double,2>& node_A, \
			const array<double,2>& node_B, \
			const array<double,2>& node_C, \
			const array<double,2>& node_D)
{
	// read in values from vectors for readability
	double xa = node_A[0];
	double xb = node_B[0];
	double xc = node_C[0];
	double xd = node_D[0];
	
	double ya = node_A[1];
	double yb = node_B[1];
	double yc = node_C[1];
	double yd = node_D[1];
	
	// calculate cross product
	double cross = (xc - xa)*(yd - yb) - (yc - ya)*(xd - xb);
	
	// calculate output
	V = 0.5*abs(cross);
}


void center(array<double,2>& cell_center,  \
			const array<double,2>& node_A, \
			const array<double,2>& node_B, \
			const array<double,2>& node_C, \
			const array<double,2>& node_D)
{
	// read in values from vectors for readability
	double xa = node_A[0];
	double xb = node_B[0];
	double xc = node_C[0];
	double xd = node_D[0];
	
	double ya = node_A[1];
	double yb = node_B[1];
	double yc = node_C[1];
	double yd = node_D[1];
	
	// calculate average x
	double xbar = 0.25*(xa + xb + xc + xd);
	double ybar = 0.25*(ya + yb + yc + yd);
	
	// read calcualted values into output vector.
	cell_center[0] = xbar;
	cell_center[1] = ybar;
}


