//
//  simulation_io.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/10/15.
//
//

#ifndef __Euler_2D_FVM__simulation_io__
#define __Euler_2D_FVM__simulation_io__

// standard includes
#include <vector>
#include <array>
#include <stdio.h>
#include <string>

using namespace std;

/* This header and accompanying source file contains all the funcitons used for 
 * input output specific to this Euler_2D_FVM simulation code. This includes a 
 * function for reading in a configuration file and a set of grid mesh 
 * coordinates and a function for outputing data for display using TecPlot.
 */

void read_config(string& sim_name, string& mesh_file_path, char& sim_case,    \
                 char& flux_method, int& sub_super, double& angle_of_attack,  \
				 double& mach, double& temperature, double& pressure,         \
				 double& CFL, int& mstage, int& order, double& kappa,         \
				 int& first_till, int& freeze_at, double& converge_till,      \
				 double& blow_up, int& max_iter, int& file_out,               \
				 int& terminal_out, string& tecplot_file, string& go_file);
/* This function will read in the contents of a file which is always named
 * 'sim_config.txt' and contained in the same directory as the source
 * code. It reads in the following information:
 * 
 * Flux method
 * Order of accuracy
 * name of file containing mesh
 * Boundary conditions
 * ouput file name
 * first order till
 * freeze limiters at
 * converge till
 * maximum iteration
 * iter out
 */

void load_mesh(string mesh_file_path,                          \
				  vector<vector<array<double,2> > >& node_loc,    \
				  int& NI, int& NJ);
/* This function will open the file specified by the 'sim_config.dat' file, and 
 * load the grid coordinates located there. The first line of the grid file
 * contains the number of zones in the mesh, which can be ignored for the
 * purposes of this simulation. The next line contains the number of mesh nodes
 * in the x, y and z directions respectively. after that the file contains all
 * of the x coordinates in terms of increasing i then j then k index, followed
 * by all y coordinates in terms of increasing i then j then k index, followed 
 * by all z coordinates in terms of increasing i then j then k index. This 
 * simulation, since it is only 2-Dimentional, will only read in the x and y 
 * coordinates, and does not store the z-coordinates in the program.
 */

void setup_output(string sim_name, string tecplot_file, string go_file);
/* This function will create the files where the tecplot output and
 * global output data will be written. It will also create the appropriate
 * header for both files.
 */

void tecplot_output(string tecplot_file, int iter,
					    const vector<vector<array<double,2> > >& node_loc,    \
					    const vector<vector<array<double,4> > >& primVar,     \
					    const vector<vector<array<double,4> > >& residual);
/* This funciton will write the desired simulation output at the current 
 * itteration formated to be read in by TecPlot. It will generate one file 
 * for the entire simulation, with each iteration of output being placed in
 * a new zone. THE 'setup_output' FUNCTION MUST BE CALLED FIRST. This 
 * function will open the previously created file in append mode.
 */

void global_output(string go_file, int iter,         \
				      const array<double,4> l1_norm,    \
				      const array<double,4> l2_norm,    \
				      const array<double,4> linf_norm);
/* this function is for outputing any global quantities of interest; i.e. 
 * quantities that are not tied to a specific grid location. THE 'setup_output'
 * FUNCTION MUST BE CALLED FIRST. This function will open the previously 
 * created file in append mode. 
 * ---------NOTE:----------
 * CURRETLY THIS FUNCTION WILL ONLY RESIDUAL NORMS FOR ONE SOLUTION VARIABLE
 */

void terminal_output(const string& sim_name,                             \
						  const clock_t& start_time, const int& iter,         \
						  const double& error_norm);
/* this funciton outputs data to the terminal wich allows us to track the 
 * progress of an ongoing simulation. it will output the total current
 * runtime, how many iterations have been performed, and the current value 
 * of the iterative residuals.
 */

void sim_end_output(const string& sim_name, const string& tecplot_file, \
					const string& go_file,  const clock_t& start_time,      \
					const int& iter);
/* this function gives the final output for the simulation, to let the user
 * know that the simulation is finished, how many itterations it took, and 
 * how long, and what the name of the output files generated are.
 */

#endif /* defined(__Euler_2D_FVM__simulation_io__) */






