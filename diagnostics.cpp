//
//  diagnostics.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 7/18/15.
//

// std. includes
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>

// my includes
#include "diagnostics.h"
#include "Euler_Helper_Functions.h"
#include "mms.h"

const int NG = 2;

using namespace std;

void mms_diagnostics(const string sim_name, const int iter, const int& CI, const int& CJ, const int& sub_super, \
					 const vector<vector<array<double,2> > >& cell_ctr,  const vector<vector<array<double,2> > >& node_loc,\
					 const vector<vector<array<double,4> > >& primVar)
{
	// Initialize an array to contain the discretization error.
    vector< vector< array<double,4> > >  DE(CI,vector< array<double,4> >(CJ));
	
	// Initialize an array to contain the exact solution.
	vector< vector< array<double,4> > > exact_solution(CI,vector< array<double,4> >(CJ));
	
	// loop over the cell center locations to calculate exact solution at
	// each cell.
	for (int i = 0; i < CI; ++i)
	{
		for (int j = 0; j < CJ; ++j)
		{
			MMS_exactsolution(sub_super, exact_solution[i][j], cell_ctr[i][j]);
		}
	}
	
	// loop over the solution in each cell to determine the discretization
	// error for each of the primitive variables in each cell.
	for (int i = 0; i < CI; ++i)
	{
		for (int j = 0; j < CJ; ++j)
		{
			for (int n = 0; n < 4; ++n)
			{
				DE[i][j][n] = primVar[i+NG][j+NG][n] - exact_solution[i][j][n];
			}
		}
	}
	
	// Initialize variable to store L1, L2, and Linf norms.
	array<double,4> DE_L1_norm;
	array<double,4> DE_L2_norm;
	array<double,4> DE_Linf_norm;
	
	// Call norm functions on DE
	L1_norm(DE_L1_norm, DE);
	L2_norm(DE_L2_norm, DE);
	Linf_norm(DE_Linf_norm, DE);
	
	// write out norms to terminal
	cout << "\n Discretization error of the conserved variables of the MMS simulation: \n";
	cout << "             Density     X-vel      Y-vel      Energy \n";
	cout << "  L1-norm: " << DE_L1_norm[0] << "    " << DE_L1_norm[1] << "    " << DE_L1_norm[2] << "    " << DE_L1_norm[3] << endl;
	cout << "  L2-norm: " << DE_L2_norm[0] << "    " << DE_L2_norm[1] << "    " << DE_L2_norm[2] << "    " << DE_L2_norm[3] << endl;
	cout << "  Linf-norm: " << DE_Linf_norm[0] << "    " << DE_Linf_norm[1] << "    " << DE_Linf_norm[2] << "    " << DE_Linf_norm[3] << endl;

	// write out discretization error of converged solutions to tecplot file.
	
	// Define file name for DE_tecplot_output
	string DE_output = sim_name+"DE_output.dat";
	
	// create DE tecplot file
	ofstream DE_tecplot(DE_output);
	
	// Throw exception if file could not be created
	if (!DE_tecplot)
	{
		cerr << "unable to create file: " << DE_output << endl;
		exit(1);
	}
	
	else
	{
		// write header for tecplot file
		DE_tecplot << "TITLE = \"" << sim_name << "\" " << endl;
		DE_tecplot << "variables = \"x\",\"y\",\"Dens-DE\",";
		DE_tecplot << "\"x-vel-DE\",\"y-vel-DE\",\"pressure-DE\" " << endl;
		
		// close tecplot file
		DE_tecplot.close();
	}

	// open tecplot file in append mode
	ofstream output(DE_output, ios::app);
	
	// Throw exception if file could not be created
	if (!output)
	{
		cerr << "unable to open file: " << DE_output << endl;
		exit(1);
	}
	
	else
	{
		// get dimentions of node_loc
		int NI = node_loc.size();
		int NJ = node_loc[0].size();
		
		// write current iteration sub-header for tecplot file
		output << "ZONE \n";
		output << "T = \"" << iter << "\" \n";
		output << "I = " << NI << ", J = " << NJ << "\n";
		output << "DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, ";
		output << "DOUBLE) \n";
		output << "DATAPACKING = BLOCK \n";
		output << "VARLOCATION = ([3-6]=CELLCENTERED) \n";
		output << "SOLUTIONTIME =" << iter << "\n";
		
		// set the precision of the output.
		output.precision(14);
		
		// write all output:
		// x-cooridnates of nodes
		for (int j = 0; j < NJ; ++j)
		{
			for (int i = 0; i < NI; ++i)
			{
				output << node_loc[i][j][0] << "\n";
			}
		}
		
		// y-coordinates of nodes
		for (int j = 0; j < NJ; ++j)
		{
			for (int i = 0; i < NI; ++i)
			{
				output << node_loc[i][j][1] << "\n";
			}
		}
		
		// density DE
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << DE[i][j][0] << "\n";
			}
		}
		
		// x-vel DE
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << DE[i][j][1] << "\n";
			}
		}
		
		// y-vel DE
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << DE[i][j][2] << "\n";
			}
		}
		
		// pressure DE
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << DE[i][j][3] << "\n";
			}
		}
		
		// close file
		output.close();
	}
}

void inlet_diagnostics(const int& CI, const int& CJ,                    \
						const array<double,4>& V_FS,                     \
						const vector<vector<array<double,4> > >& primVar)
{
	// initialize variables that will be used here.
	double Pt_FS = 0.0;  // total pressure of the free-stream condition.
	double Pt_avg = 0.0; // average total pressure at the outlet.
	
	// calculate the total pressure of the free-stream condition
	Pt_FS = total_pressure(V_FS);
	
	// loop over the j indicies for the last row of cells in the simulation
	for (int j = 0; j < CJ; ++j)
	{
		Pt_avg += total_pressure(primVar[2+CI][j+2]);
	}
	
	// divide by the number of cells averaged over
	Pt_avg = Pt_avg/(1.0*CJ);
	
	// output results to the terminal.
	cout << "\n Total pressure of the free-stream condition: " << Pt_FS << " Pa \n";
	cout << "Average total pressure at the outlet: " << Pt_avg << " Pa" << endl;
}


void airfoil_diagnostics(const int& CI, const int& CJ, const int& edge,     \
						 const double& AoA, const string& sim_name,         \
						 const vector<vector<array<double,2> > >& node_loc, \
						 const vector<vector<array<double,2> > >& eta_norm,  \
						 const vector<vector<double> >& eta_A,               \
						 const array<double,4>& V_FS,                       \
						 const vector<vector<array<double,4> > >& primVar)
{
	// define variables for force components in (x,y). and (L,D)
	double Fx = 0.0;
	double Fy = 0.0;
	double L = 0.0;
	double D = 0.0;
	
	// convert AoA into radians.
	double alpha = AoA*3.14159265359/180.0;
	
	// determine chordlength from x location of trailing edge.
	double chord = node_loc[edge-NG][0][0];
	
	// loop over cells to sum up the force on each cell boundary in the (x,y) frame.
	for (int i = edge; i < (CI - edge); ++i)
	{
		// force is givein by -1*pressure*area*normalvector.
		Fx += -1.0*primVar[i+NG][NG][3]*eta_A[i][0]*eta_norm[i][0][0];
		Fy += -1.0*primVar[i+NG][NG][3]*eta_A[i][0]*eta_norm[i][0][1];
	}
	
	// rotate total force into the free-stream direction to determine lift and drag.
	D = cos(alpha)*Fx + sin(alpha)*Fy;       // ROTATION OF COORDINATE FRAME:
	L = -1.0*sin(alpha)*Fx + cos(alpha)*Fy;  // EQUIVALENT TO VECTOR ROTATION OF -ALPHA
	
	// Determine lift and drag coeficients.
	double Cl = 0.0;
	double Cd = 0.0;
	double Vsqrd = V_FS[1]*V_FS[1] + V_FS[2]*V_FS[2];
	
	Cl = L/(0.5*V_FS[0]*Vsqrd*chord);
	Cd = D/(0.5*V_FS[0]*Vsqrd*chord);
	
	// define array to hold all Cp values and x/c values.
	int panels = CI - 2*edge;

	vector<double> Cp;
	Cp.resize(panels);
	
	vector<double> x_c;
	x_c.resize(panels);
	
	// calclulate Cp and x/c at each cell
	for (int i = edge; i < (CI - edge); ++i)
	{
		Cp[i-edge] = (primVar[i+NG][NG][3] - V_FS[3])/(0.5*V_FS[0]*Vsqrd);
		x_c[i-edge] = 0.5*(node_loc[i][0][0] + node_loc[i+1][0][0])/chord;
	}
	
	// create file where Cp values will be stored.
	string Cp_output = sim_name+"_Cp_output.txt";
	
	// create Cp text file
	ofstream Cp_txt(Cp_output);
	
	// Throw exception if file could not be created
	if (!Cp_txt)
	{
		cerr << "unable to create file: " << Cp_output << endl;
		exit(1);
	}
	
	else
	{
		// write out header for data
		Cp_txt << "  xc  " << "  Cp  " << endl;
		
		// write out all values of x/c and Cp in two columns
		for (int i = 0; i < panels; ++i)
		{
			Cp_txt << x_c[i] << "  " << Cp[i] << endl;
		}
	}
	
	// finaly, output Cl and Cd to terminal
	cout << "\n angle of attack: " << AoA << " degrees" << endl;
	cout << " Coefficient of lift: " << Cl << endl;
	cout << " Coefficient of drag: " << Cd << endl;
	
}