//
//  Van_Leer_Flux_Test.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/2/15.
//
//

#define CATCH_CONFIG_MAIN

#include <stdio.h>
#include <array>
#include <vector>
#include <string>

#include "catch.hpp"
#include "Van_Leer_Flux.h"
#include "Roe_Flux.h"
#include "Geometry.h"
#include "Time_Step.h"
#include "Flux_Determiner.h"
#include "Roe_Flux_Determiner.h"
#include "simulation_io.h"

using namespace std;
/*
TEST_CASE( "Check that Roe and Van Leer produce the same result", "[Flux]" )
{
	// define input vectors
	array<double,4> V_L = {{1, -1000, 100, 100000}};
	
	array<double,4> V_R = {{1, -1000, 100, 10000}};

	
	array<double,2> n = {{-0.5*sqrt(2.0), 0.5*sqrt(2.0)}};

	// define output vectors
	array<double,4> F_VL;
	array<double,4> F_ROE;

	
	// call Van Leer
	Van_Leer_Flux(F_VL, V_L, V_R, n);
	
	// call Roe
	Roe_Flux(F_ROE, V_L, V_R, n);
	
	// check that the outputs are the same
	CHECK( F_VL[0] == Approx( F_ROE[0] ) );
	CHECK( F_VL[1] == Approx( F_ROE[1] ) );
	CHECK( F_VL[2] == Approx( F_ROE[2] ) );
	CHECK( F_VL[3] == Approx( F_ROE[3] ) );
}

TEST_CASE( "Hunt Down SegFault in Roe", "[Roe]")
{
	// define input vectors for use
	array<double,4> V_L = {{1, -1000, 100, 100000}};
	array<double,4> V_R = {{1.1, -950, 200, 110000}};
	array<double,2> n = {{1, 0}};
	
	// check vector sizes
	REQUIRE( V_L.size() == 4 );
	REQUIRE( V_R.size() == 4 );
	REQUIRE(   n.size() == 2 );
	
	// calculate roe averaged variables
	array<double,5> roe_vars;
	roe_avg(roe_vars,V_L,V_R);
	
	REQUIRE( roe_vars.size() == 5 );

	// calculate eigenvalues
	array<double,4> lambda;
	roe_eigval(lambda,roe_vars,n);
		
	CHECK( lambda[0] == lambda[1] );
	CHECK( (lambda[2] - lambda[3]) == Approx( 2.0*roe_vars[4] ) );
	
	// calculate eigenvectors
    array<array<double,4>,4> rightEig;
	
	CHECK( rightEig.size()    == 4 );
	CHECK( rightEig[0].size() == 4 );
	
	roe_eigvec(rightEig,roe_vars,n);
	
	// calculate wave amplitudes
	array<double,4> dw;
	wave_amp(dw,V_L,V_R,roe_vars,n);
	
	CHECK( dw.size() == 4 );
	
	// calculate centered flux
	array<double,4> w_flux;
	
	// calculate wave flux
	for (int i=0; i < 4; ++i)
	{
		w_flux[i] = 0.0; // reset value at memory address.
		for (int j=0; j < 4; ++j)
		{
			w_flux[i] += 0.5*abs(lambda[j])*dw[j]*rightEig[i][j];
		}
	}
}
 */


TEST_CASE( "Test Geometry", "[Geometry]" )
{
	// create a simple cartesian array
	vector< vector< array<double,2> > > node_Loc(10, vector<array<double,2> >(10));
															   
	REQUIRE( node_Loc.size() == 10);
	REQUIRE( node_Loc[0].size() == 10);
	REQUIRE( node_Loc[0][0].size() == 2);
	
	// load the array
	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			node_Loc[i][j][0] = 0.1*i;
			node_Loc[i][j][1] = 0.1*j;
		}
	}
	REQUIRE( node_Loc[2][2][0] == node_Loc[2][2][1] );
	
	int NI = node_Loc.size();
	int NJ = node_Loc[0].size();
	
	REQUIRE( NI == 100 );
	REQUIRE( NJ == 200 );
	
	// initialize geometry variables
	vector< vector<double> > xi_A(1,vector<double>(1));
	vector< vector<double> > eta_A(1,vector<double>(1));
	vector< vector< array<double,2> > > xi_N(1,vector< array<double,2> >(1));
	vector< vector< array<double,2> > > eta_N(1,vector< array<double,2> >(1));
	vector< vector<double> > cell_V(1,vector<double>(1));
	vector< vector< array<double,2> > > cell_Loc(1,vector< array<double,2> >(1));

	// check my understanding of vector initialization
	CHECK( xi_N.size() == 1 );
	CHECK( xi_N[0].size() == 1 );
	CHECK( xi_N[0][0].size() == 2 );
	
    // call set geometry function
	Geometry(xi_A, eta_A, xi_N, eta_N, cell_V, cell_Loc, node_Loc);
	
	// check output vector size

	REQUIRE( xi_A.size() == NI );
	REQUIRE( xi_A[NI-1].size() == NJ-1 );
	REQUIRE( xi_N.size() == NI );
	REQUIRE( xi_N[0].size() == NJ-1 );
	REQUIRE( xi_N[2][1].size() == 2 );
	
	REQUIRE( eta_A.size() == NI-1 );
	REQUIRE( eta_A[0].size() == NJ );
	REQUIRE( eta_N.size() == NI-1 );
	REQUIRE( eta_N[0].size() == NJ );
	REQUIRE( eta_N[0][0].size() == 2 );
	
	REQUIRE( cell_V.size() == NI-1 );
	REQUIRE( cell_V[NI-2].size() == NJ-1 );
	REQUIRE( cell_Loc.size() == NI-1 );
	REQUIRE( cell_Loc[NI-2].size() == NJ-1 );
	REQUIRE( cell_Loc[NI-2][NJ-2].size() == 2);

	
	// Check that calculations are being performed correctly.
	REQUIRE( cell_V[79][85] == cell_V[20][176] );
	REQUIRE( cell_V[NI-2][NJ-2] == 1.0 );
	
	REQUIRE( xi_A[0][0]  == 1.0 );
	REQUIRE( eta_A[1][1] == 1.0 );
	
	REQUIRE( cell_Loc[0][1][0] == 0.5 );
	REQUIRE( cell_Loc[0][1][1] == 1.5 );
	
	REQUIRE( xi_N[NI-1][NJ-2][0] == 1.0 );
	REQUIRE( xi_N[34][23][1] == 0.0 );
	
	REQUIRE( eta_N[89][54][0] == 0.0 );
	REQUIRE( eta_N[57][34][1] == 1.0 );
	

	double CFL = .1;
	array<double,4> vars = {1,100,100,100000};
	vector< vector< array<double,4> > > primVar(NI-1, vector< array<double,4> >(NJ-1, vars));
	
	vector< vector<double> > dt(NI-1, vector<double>(NJ-1, 0.0));
	
	Time_Step(dt, CFL, cell_V, xi_A, eta_A, xi_N, eta_N, primVar);
	
	CHECK( dt[0][0] == Approx(1.05448361e-4) );

	
}
/*

int random(Flux_Determiner& fd) {
	vector<double> ting(4, 1);
	cout << ting.size() << endl;
	fd.flux(ting);
	cout << ting.size() << endl;
	
}

TEST_CASE( " class constructor test")
{
	vector<double> dummy(4);

	Roe_Flux_Determiner RFD = Roe_Flux_Determiner(dummy, dummy, dummy);
	Flux_Determiner FD = RFD;
	REQUIRE( (*FD).leftSize() == 4 );
	
	//random(FD);
	FD->flux(dummy);
	
	REQUIRE( dummy.size() == 5 );
}


*/

TEST_CASE( " simulation I/O " )
{
	// initialize parameter variables
	string sim_name;
	string mesh_file_path;
	char flux_method;
	int order;
	int first_till;
	int freeze_at;
	double converge_till;
	int max_iter;
	int iter_out;
	string tecplot_file;
	string go_file;
	
	// read configfile
	read_config(sim_name, mesh_file_path, flux_method, order, first_till, freeze_at, \
				converge_till, max_iter, iter_out, tecplot_file, go_file);
	
	// check that input file was read correctly
	REQUIRE( sim_name == "mms_ROE_o1" );
	REQUIRE( mesh_file_path == "curv2d9.grd" );
	REQUIRE( flux_method == 'r' );
	REQUIRE( order == 1 );
	REQUIRE( first_till == 1000 );
	REQUIRE( freeze_at == 3000 );
	REQUIRE( converge_till == 1.0e-1 );
	REQUIRE( max_iter == 1000 );
	REQUIRE( iter_out == 500 );
	REQUIRE( tecplot_file == "mms_ROE_o1_tecplot_output.dat");
	REQUIRE( go_file == "mms_ROE_o1_global_output.dat");
	
	
	// now try to load the grid.
	vector<vector<array<double,2> > > node_loc(1,vector<array<double,2> >(1));

	load_mesh(mesh_file_path, node_loc);
	
	// check that mesh was successfully loaded!
	REQUIRE( node_loc.size()       == 9 );
	REQUIRE( node_loc[0].size()    == 9 );
	REQUIRE( node_loc[5][7].size() == 2 );
	REQUIRE( node_loc[3][0][0] == 5.45319914800e-01 );
	
	// make dummy data
	int iter = 145;
	
	array<double,4> vars = {{1,100,100,100000}};
	vector< vector< array<double,4> > > primVar(8, vector< array<double,4> >(8, vars));

	array<double,4> res = {{3.14, 31.4, 31.4, 314}};
	vector< vector< array<double,4> > > residual(8, vector< array<double,4> >(8, res));
	
	double l1_norm = 0.2357;
	double l2_norm = 0.235711;
	double linf_norm = 0.23571113;
	
	// call setup_output, tecplot_output, global_output
	setup_output(sim_name, tecplot_file, go_file);
	tecplot_output(tecplot_file, iter, node_loc, primVar, residual);
	global_output(go_file, iter, l1_norm, l2_norm, linf_norm);
	
	// open files to see how it went!
}







