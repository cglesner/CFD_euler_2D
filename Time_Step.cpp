//
//  Time_Step.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/4/15.
//
//

// standard includes
#include <cmath>
#include <vector>
#include <array>


// my includes
#include "Time_Step.h"
#include "Euler_Helper_Functions.h"

const int NG = 2;

using namespace std;

void Time_Step(vector<vector<double> >& dt,                            \
			   const double& CFL,                                      \
			   const vector<vector<double> >& cell_Volume,             \
			   const vector<vector<double> >& xi_Area,                 \
			   const vector<vector<double> >& eta_Area,                \
			   const vector<vector<array<double,2> > >& xi_Normal,     \
			   const vector<vector<array<double,2> > >& eta_Normal,    \
			   const vector<vector<array<double,4> > >& primitive_vars)
{
	/* determine the size of the mesh in the i and j directions. take the size 
	 * of the cell_Volume array and add one, so as to recover the same NI and
	 * NJ used in the 'Geometry' function, so that indexing across functions
	 * will be consistent.
	 */
	int NI = cell_Volume.size() + 1;
	int NJ = cell_Volume[0].size() + 1;
	
	// declare all neded variables
	array<double,2> avg_xi_Normal;
	array<double,2> avg_eta_Normal;
	
	double avg_xi_Area;
	double avg_eta_Area;
	
	double xi_lambda_max;
	double eta_lambda_max;
	
	// placeholders
	double denom;
    double numer;
	
	// nested loop over cells
	for (int i = 0; i < (NI-1); ++i)
	{
		for (int j = 0; j < (NJ-1); ++j)
		{
			// calculate averaged normals and face areas
			avg_xi_Area  = 0.5*( xi_Area[i][j] +  xi_Area[i+1][j]);
			avg_eta_Area = 0.5*(eta_Area[i][j] + eta_Area[i][j+1]);
			
			avg_xi_Normal[0]  = 0.5*( xi_Normal[i][j][0] +  xi_Normal[i+1][j][0]);
			avg_xi_Normal[1]  = 0.5*( xi_Normal[i][j][1] +  xi_Normal[i+1][j][1]);
			
		    avg_eta_Normal[0] = 0.5*(eta_Normal[i][j][0] + eta_Normal[i][j+1][0]);
			avg_eta_Normal[1] = 0.5*(eta_Normal[i][j][1] + eta_Normal[i][j+1][1]);
			
			// calculate max eigenvalues
			xi_lambda_max  = max_lambda(primitive_vars[i+NG][j+NG], avg_xi_Normal );
			eta_lambda_max = max_lambda(primitive_vars[i+NG][j+NG], avg_eta_Normal);
			
			// calculate local max timestep
			denom = xi_lambda_max*avg_xi_Area + eta_lambda_max*avg_eta_Area;
           numer = CFL*cell_Volume[i][j];
            
			dt[i][j] = numer/denom;
		}
	}
}


double max_lambda(const array<double,4>& primitive_vars, \
				  const array<double,2>& Normal_avg )
{
	// read in values from vectors for readability
	double vx = primitive_vars[1];
	double vy = primitive_vars[2];
	
	double nx = Normal_avg[0];
	double ny = Normal_avg[1];
	
	// calculate speed of sound
	double a = sound_speed(primitive_vars);
	
	// calculate output
	return abs(nx*vx + ny*vy) + a;
}




