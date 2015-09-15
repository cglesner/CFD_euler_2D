//
//  simulation_io.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/10/15.
//
//

// standard includes
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>

// my includes
#include "simulation_io.h"
const int NG = 2;

using namespace std;


void read_config(string& sim_name, string& mesh_file_path, char& sim_case,    \
				 char& flux_method, int& sub_super, double& angle_of_attack,  \
				 double& mach, double& temperature, double& pressure,         \
				 double& CFL, int& mstage, int& order, double& kappa,         \
				 int& first_till, int& freeze_at, double& converge_till,      \
				 double& blow_up, int& max_iter, int& file_out,               \
				 int& terminal_out, string& tecplot_file, string& go_file)
{
	// open file
	ifstream input("sim_config.txt");

	// set the precission of the input:
	input.precision(14);
	
	// throw exception if file cannot be opened.
	if (!input)
	{
		cerr << "could not open sim_config.txt" << endl;
		exit(1);
	}
	
	else
	{
		/* file is formatted in such a way that it will always have:
		 * paramteter_name <whitespace> parameter <whitespace>. This code 
		 * is therefore written to ignore odd inputs and store even inputs
		 * for the first set of inputs before the boundary conditions.
		 */
		
		//declare a place holder string, simulation name.
		string dummy;

		
		input >> dummy >> sim_name;
		input >> dummy >> mesh_file_path;
		input >> dummy >> sim_case;
		input >> dummy >> flux_method;
		input >> dummy >> sub_super;
		input >> dummy >> angle_of_attack;
		input >> dummy >> mach;
		input >> dummy >> temperature;
		input >> dummy >> pressure;
		input >> dummy >> CFL;
		input >> dummy >> mstage;
		input >> dummy >> order;
		input >> dummy >> kappa;
		input >> dummy >> first_till;
		input >> dummy >> freeze_at;
		input >> dummy >> converge_till;
		input >> dummy >> blow_up;
		input >> dummy >> max_iter;
		input >> dummy >> file_out;
		input >> dummy >> terminal_out;
		
		// close the input file
		input.close();
	}
	
	// create file names for tecplot_file, go_file.
	tecplot_file = sim_name+"_tecplot_output.dat";
	go_file = sim_name+"_global_output.dat";
}


void load_mesh(string mesh_file_path,                           \
				  vector<vector<array<double,2> > >& node_loc,     \
				  int& NI, int& NJ)
{
	// open file
	ifstream input(mesh_file_path);
	
	// set the precission of the input
	input.precision(14);
	
	// throw exception if file cannot be opened.
	if (!input)
	{
		cerr << "could not open mesh file at path: " << mesh_file_path << endl;
		exit(1);
	}
	
	else
	{
		// read in header
		int nzones;
		input >> nzones;  //number of zones
		
		input >> NI;      // number of Nodes in x direction
		input >> NJ;      // number of Nodes in y direction
		
		int NK;
		input >> NK;      // number of Nodes in z direction
		
		// resize node_loc container appropriately
		node_loc.resize(NI);
		for	( int i = 0; i < NI; ++i)
		{
			node_loc[i].resize(NJ);
		}
		
		/* because this is a code where we don't care about the z-coordinates,
		 * we will not be reading them in. Also, if NK = 2 is given, as is the 
		 * case for the MMS curvilinear grids it is assumed that the x and y 
		 * values are duplcated, and so they are simply read in twice; 
		 * overwriting the first write. They would be skipped if the data did
		 * not need to be read through systematically.
		 */
		
		// read in all x-locations
		for (int k = 0; k < NK; ++k)
		{
			for (int j = 0; j < NJ; ++j)
			{
				for (int i = 0; i < NI; ++i)
				{
					input >> node_loc[i][j][0];
				}
			}
		}
		
		// read in all y-locations
		for (int k = 0; k < NK; ++k)
		{
			for (int j = 0; j < NJ; ++j)
			{
				for (int i = 0; i < NI; ++i)
				{
					input >> node_loc[i][j][1];
				}
			}
		}
		
		// close file
		input.close();
	}
}


void setup_output(string sim_name, string tecplot_file, string go_file)
{
	// create tecplot file
	ofstream tecplot(tecplot_file);
	
	// Throw exception if file could not be created
	if (!tecplot)
   {
		cerr << "unable to create file: " << tecplot_file << endl;
		exit(1);
	}
	
	else
	{
		// write header for tecplot file
		tecplot << "TITLE = \"" << sim_name << "\" " << endl;
		tecplot << "variables = \"x\",\"y\",\"Density\",";
		tecplot << "\"x-vel\",\"y-vel\",\"pressure\",";
		tecplot << "\"R1\",\"R2\",\"R3\",\"R4\"" << endl;
		
		// close tecplot file
		tecplot.close();
	}
	
	// create global output file
	ofstream global(go_file);
	
	// Throw exception if file could not be created
	if (!global)
	{
		cerr << "unable to create file: " << go_file << endl;
		exit(1);
	}
	
	else
	{
		// write header for global file
		global << "TITLE = " << sim_name << endl;
		global << " iter " << " l1_norm " << " l2_norm " << " linf_norm " << endl;
		
		// close global file
		global.close();
	}
}


void tecplot_output(string tecplot_file, int iter,
					    const vector<vector<array<double,2> > >& node_loc,    \
					    const vector<vector<array<double,4> > >& primVar,     \
					    const vector<vector<array<double,4> > >& residual)
{
	// open tecplot file in append mode
	ofstream output(tecplot_file, ios::app);
	
	// Throw exception if file could not be created
	if (!output)
	{
		cerr << "unable to open file: " << tecplot_file << endl;
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
		output << "DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE) \n";
		output << "DATAPACKING = BLOCK \n";
		output << "VARLOCATION = ([3-10]=CELLCENTERED) \n";
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
		
		// density
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << primVar[i+NG][j+NG][0] << "\n";
			}
		}
		
		// x-vel
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << primVar[i+NG][j+NG][1] << "\n";
			}
		}
		
		// y-vel
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << primVar[i+NG][j+NG][2] << "\n";
			}
		}
		
		// pressure
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << primVar[i+NG][j+NG][3] << "\n";
			}
		}
		
		// residual mass
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << residual[i][j][0] << "\n";
			}
		}
		
		// residual x-mmt
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << residual[i][j][1] << "\n";
			}
		}
		
		// residual y-mmt
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << residual[i][j][2] << "\n";
			}
		}

		// residual energy
		for (int j = 0; j < NJ-1; ++j)
		{
			for (int i = 0; i < NI-1; ++i)
			{
				output << residual[i][j][3] << "\n";
			}
		}

      // close file
		output.close();
	}
}



void global_output(string go_file, int iter,      \
				   const array<double,4> l1_norm,    \
				   const array<double,4> l2_norm,    \
				   const array<double,4> linf_norm)
{
	// open global output file in append mode
	ofstream output(go_file, ios::app);
	
	// Throw exception if file could not be created
	if (!output)
	{
		cerr << "unable to open file: " << go_file << endl;
		exit(1);
	}
	
	else
	{
		
		// set the precision of the output.
		output.precision(14);

		// write out norms
		output << " " << iter << " " << " ";
		output << l1_norm[0] << " " << l1_norm[1] << " " << l1_norm[2] << " " << l1_norm[3] << " ";
		output << l2_norm[0] << " " << l2_norm[1] << " " << l2_norm[2] << " " << l2_norm[3] << " ";
		output << linf_norm[0] << " " << linf_norm[1] << " " << linf_norm[2] << " " << linf_norm[3] << "\n";
		// close file
		output.close();
	}
}

void terminal_output(const string& sim_name,                         \
					 const clock_t& start_time, const int& iter,         \
					 const double& error_norm)
{
	// calculate the current runtime
	double duration = ( clock() - start_time ) / (double) CLOCKS_PER_SEC;
	// write file name, current iteration, and total simulaiton time elapsed.
	cout << "Simulation case currently running: \n";
   cout << sim_name << endl;
	cout.width(20);
	cout.fill(' ');
	cout << left << "iterations completed: " << iter;
	cout.width(20);
	cout.fill(' ');
	cout << left << "  current runtime: " << duration << "\n \n";
	cout << left << "l2_norms of the iterative residuals: \n";
	cout.width(10);
	cout.fill(' ');
	cout << left << "Residual Mass  " << error_norm << "\n" << endl;
	/*cout.width(20);
	cout.fill(' ');
	cout << left << "Residual X-momentum";
	cout.width(20);
	cout.fill(' ');
	cout << left << "Residual Y-momentum";
	cout.width(20);
	cout.fill(' ');
	cout << left << "Residual Energy" << endl;
	cout.width(20);
	cout.fill(' ');
	 */

	/*cout << left << l2_norm[1];
	cout.width(20);
	cout.fill(' ');
	cout << left << l2_norm[2];
	cout.width(20);
	cout.fill(' ');
	cout << left << l2_norm[3];
	 */
}

void sim_end_output(const string& sim_name, const string& tecplot_file, \
					    const string& go_file,  const clock_t& start_time,  \
					    const int& iter)
{
	// calculate total runtime
	double total_duration;
	total_duration = ( clock() - start_time ) / (double) CLOCKS_PER_SEC;
	// notify user of completed simulation
	cout << "\n Simulation case: " << sim_name << " complete. \n";
	cout << " Solution converged in " << iter << " iterations. \n";
	cout << " Total run time: " << total_duration << " seconds \n";
	cout << " generated output files: " << tecplot_file << "  " << go_file << endl;
}


