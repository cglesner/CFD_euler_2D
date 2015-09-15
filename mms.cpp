//
//  mms.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 5/12/15.


// standard includes
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <cstdlib>

// my includes
#include "mms.h"
const int NG = 2;

using namespace std;

void MMS_initial_and_boundary(int sub_super, vector<vector<array<double,4> > >& primVar,               \
						   const vector<vector<array<double,2> > >& cell_and_ghost_ctr, \
						   const int& CI, const int& CJ)
{
	double rho0;
	double uvel0;
	double vvel0;
	double press0;
	
	if (sub_super == 0)
	{
		// constants used in this exact MMS solution.
		rho0 = 1.0;
		uvel0 = 70.0;
		vvel0 = 90.0;
		press0 = 1.0e+5;
	}
	else
	{
		if (sub_super == 1)
		{
			rho0 = 1.0;
			uvel0 = 800.0;
			vvel0 = 800.0;
			press0 = 1.0e+5;
		}
		else
		{
			cerr << " invalid choice for sub_super. " << endl;
			exit(1);
		}
	}

	/*
	// initialize all cells in primVar exact solution.
	for (int i = NG; i < CI+NG; ++i)
	{
		for (int j = NG; j < CJ+NG; ++j)
		{
			MMS_exactsolution(sub_super,primVar[i][j],cell_and_ghost_ctr[i][j]);
		}
	}
	*/
	
    // set xi-lower bc.
	for (int i = NG; i < CI+NG; ++i)
	{
		for (int j =0; j < NG; ++j)
		{
			MMS_exactsolution(sub_super,primVar[i][j],cell_and_ghost_ctr[i][j]);
		}
	}
	
	// set xi-upper bc.
	for (int i = NG; i < CI+NG; ++i)
	{
		for (int j = CJ+NG; j < CJ+2*NG; ++j)
		{
			MMS_exactsolution(sub_super,primVar[i][j],cell_and_ghost_ctr[i][j]);
		}
	}
	
	// set eta-lower bc.
	for (int i = 0; i < NG; ++i)
	{
		for (int j = NG; j < CJ+NG; ++j)
		{
			MMS_exactsolution(sub_super,primVar[i][j],cell_and_ghost_ctr[i][j]);
		}
	}
	
	// set eta-upper bc.
	for (int i = NG+CI; i < CI+2*NG; ++i)
	{
		for (int j = NG; j < CJ+NG; ++j)
		{
			MMS_exactsolution(sub_super,primVar[i][j],cell_and_ghost_ctr[i][j]);
		}
	}

	
	// set interior values
	for (int i = NG; i < CI+NG; ++i)
	{
		for (int j = NG; j < CJ+NG; ++j)
		{
           //MMS_exactsolution_subsonic(primVar[i][j],cell_and_ghost_ctr[i][j]);
		   primVar[i][j][0] = rho0;
		   primVar[i][j][1] = uvel0;
		   primVar[i][j][2] = vvel0;
		   primVar[i][j][3] = press0;
		}
	}
	
	
}
void MMS_exactsolution(int sub_super, array<double,4>& primVar,          \
					            const array<double,2>& cell_ctr)
{
	double rho0;
	double rhox;
	double rhoy;
	
	double uvel0;
	double uvelx;
	double uvely;
	
	double vvel0;
	double vvelx;
	double vvely;
	
	double press0;
	double pressx;
	double pressy;
	
	double L = 1.0;
	double Pi = 3.141592654;
	
	if (sub_super == 0)
	{
		// constants used in this exact MMS solution (subsonic case).
		rho0 = 1.0;
		rhox = 0.15;
		rhoy = -0.1;
		
		uvel0 = 70.0;
		uvelx = 5.0;
		uvely = -7.0;
		
		vvel0 = 90.0;
		vvelx = -15.0;
		vvely = 8.5;
		
		press0 = 1.0e+5;
		pressx = 0.2e+5;
		pressy = 0.5e+5;

	}
	else
	{
		if (sub_super == 1)
		{
			rho0 = 1.0;
			rhox = 0.15;
			rhoy = -0.1;
			
			uvel0 = 800.0;
			uvelx = 50.0;
			uvely = -30.0;
			
			vvel0 = 800.0;
			vvelx = -75.0;
			vvely = 40.0;
			
			press0 = 1.0e+5;
			pressx = 0.2e+5;
			pressy = 0.5e+5;
		}
		else
		{
			cerr << " invalid choice for sub_super. " << endl;
			exit(1);
		}
	}
	
	
	// read in coordinates of cell_ctr
	double x = cell_ctr[0];
	double y = cell_ctr[1];
	
	// density exact solution
	primVar[0] = rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L);
	
	// u-velocity exact solution
	primVar[1] = uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L));
	
	// v-velocity exact solution
	primVar[2] = vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L));
	
	// pressure exact solution
	primVar[3] = press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L);
}


void MMS_source(int sub_super, array<double,4>& source,         \
						 const array<double,2>& cell_ctr)
{
	double rho0;
	double rhox;
	double rhoy;
	
	double uvel0;
	double uvelx;
	double uvely;
	
	double vvel0;
	double vvelx;
	double vvely;
	
	double press0;
	double pressx;
	double pressy;
	
	// constants used in this exact MMS solution.
	if (sub_super == 0)
	{
		// constants used in this exact MMS solution (subsonic case).
		rho0 = 1.0;
		rhox = 0.15;
		rhoy = -0.1;
		
		uvel0 = 70.0;
		uvelx = 5.0;
		uvely = -7.0;
		
		vvel0 = 90.0;
		vvelx = -15.0;
		vvely = 8.5;
		
		press0 = 1.0e+5;
		pressx = 0.2e+5;
		pressy = 0.5e+5;
		
	}
	else
	{
		if (sub_super == 1)
		{
			rho0 = 1.0;
			rhox = 0.15;
			rhoy = -0.1;
			
			uvel0 = 800.0;
			uvelx = 50.0;
			uvely = -30.0;
			
			vvel0 = 800.0;
			vvelx = -75.0;
			vvely = 40.0;
			
			press0 = 1.0e+5;
			pressx = 0.2e+5;
			pressy = 0.5e+5;
		}
		else
		{
			cerr << " invalid choice for sub_super. " << endl;
			exit(1);
		}
	}

	double L = 1.0;
	double Pi = 3.141592654;
	double gamma = 1.4;
	
	// read in coordinates of cell_ctr
	double x = cell_ctr[0];
	double y = cell_ctr[1];
	
	// mass source
	source[0] = (3*Pi*uvelx*cos((3*Pi*x*pow(L,-1))/2.)*
				 pow(L,-1)*(rho0 +
							rhoy*cos((Pi*y*pow(L,-1))/2.) +
							rhox*sin(Pi*x*pow(L,-1))))/2. +
	(2*Pi*vvely*cos((2*Pi*y*pow(L,-1))/
					3.)*pow(L,-1)*
	 (rho0 + rhoy*
	  cos((Pi*y*pow(L,-1))/2.) +
	  rhox*sin(Pi*x*pow(L,-1))))/3. +
	Pi*rhox*cos(Pi*x*pow(L,-1))*
	pow(L,-1)*(uvel0 +
      uvely*cos((3*Pi*y*pow(L,-1))/
				5.) + uvelx*
			   sin((3*Pi*x*pow(L,-1))/2.)) -
	(Pi*rhoy*pow(L,-1)*
	 sin((Pi*y*pow(L,-1))/2.)*
	 (vvel0 + vvelx*
	  cos((Pi*x*pow(L,-1))/2.) +
	  vvely*sin((2*Pi*y*pow(L,-1))/3.)
	  ))/2.;
	
	// x-momentum source
	source[1] = Pi*rhox*cos(Pi*x*pow(L,-1))*pow(L,-1)*
	pow(uvel0 + uvely*
		cos((3*Pi*y*pow(L,-1))/5.) +
		uvelx*sin((3*Pi*x*pow(L,-1))/2.),
		2) + 3*Pi*uvelx*
	cos((3*Pi*x*pow(L,-1))/2.)*
	pow(L,-1)*(rho0 +
      rhoy*cos((Pi*y*pow(L,-1))/2.) +
      rhox*sin(Pi*x*pow(L,-1)))*
	(uvel0 + uvely*
	 cos((3*Pi*y*pow(L,-1))/5.) +
	 uvelx*sin((3*Pi*x*pow(L,-1))/2.))\
	+ (2*Pi*vvely*
	   cos((2*Pi*y*pow(L,-1))/3.)*
	   pow(L,-1)*(rho0 +
				  rhoy*cos((Pi*y*pow(L,-1))/2.) +
				  rhox*sin(Pi*x*pow(L,-1)))*
	   (uvel0 + uvely*
		cos((3*Pi*y*pow(L,-1))/5.) +
		uvelx*sin((3*Pi*x*pow(L,-1))/2.)
		))/3. - 2*Pi*pressx*pow(L,-1)*
	sin(2*Pi*x*pow(L,-1)) -
	(Pi*rhoy*pow(L,-1)*
	 (uvel0 + uvely*
	  cos((3*Pi*y*pow(L,-1))/5.) +
	  uvelx*sin((3*Pi*x*pow(L,-1))/2.)
	  )*sin((Pi*y*pow(L,-1))/2.)*
	 (vvel0 + vvelx*
	  cos((Pi*x*pow(L,-1))/2.) +
	  vvely*sin((2*Pi*y*pow(L,-1))/3.)
	  ))/2. - (3*Pi*uvely*pow(L,-1)*
      (rho0 + rhoy*
	   cos((Pi*y*pow(L,-1))/2.) +
	   rhox*sin(Pi*x*pow(L,-1)))*
      sin((3*Pi*y*pow(L,-1))/5.)*
      (vvel0 + vvelx*
	   cos((Pi*x*pow(L,-1))/2.) +
	   vvely*sin((2*Pi*y*pow(L,-1))/3.)
	   ))/5.;
	
	// y-momentum source
	source[2] = Pi*pressy*cos(Pi*y*pow(L,-1))*
	pow(L,-1) - (Pi*vvelx*pow(L,-1)*
				 sin((Pi*x*pow(L,-1))/2.)*
				 (rho0 + rhoy*
				  cos((Pi*y*pow(L,-1))/2.) +
				  rhox*sin(Pi*x*pow(L,-1)))*
				 (uvel0 + uvely*
				  cos((3*Pi*y*pow(L,-1))/5.) +
				  uvelx*sin((3*Pi*x*pow(L,-1))/2.)
				  ))/2. - (Pi*rhoy*pow(L,-1)*
						   pow(vvel0 +
							   vvelx*cos((Pi*x*pow(L,-1))/
           2.) +
							   vvely*sin((2*Pi*y*pow(L,-1))/3.)
							   ,2)*sin((Pi*y*pow(L,-1))/2.))/2.
	+ (3*Pi*uvelx*
	   cos((3*Pi*x*pow(L,-1))/2.)*
	   pow(L,-1)*(rho0 +
				  rhoy*cos((Pi*y*pow(L,-1))/2.) +
				  rhox*sin(Pi*x*pow(L,-1)))*
	   (vvel0 + vvelx*
		cos((Pi*x*pow(L,-1))/2.) +
		vvely*sin((2*Pi*y*pow(L,-1))/3.)
		))/2. + (4*Pi*vvely*
				 cos((2*Pi*y*pow(L,-1))/3.)*
				 pow(L,-1)*(rho0 +
							rhoy*cos((Pi*y*pow(L,-1))/2.) +
							rhox*sin(Pi*x*pow(L,-1)))*
				 (vvel0 + vvelx*
				  cos((Pi*x*pow(L,-1))/2.) +
				  vvely*sin((2*Pi*y*pow(L,-1))/3.)
				  ))/3. + Pi*rhox*
	cos(Pi*x*pow(L,-1))*pow(L,-1)*
	(uvel0 + uvely*
	 cos((3*Pi*y*pow(L,-1))/5.) +
	 uvelx*sin((3*Pi*x*pow(L,-1))/2.))*
	(vvel0 + vvelx*
	 cos((Pi*x*pow(L,-1))/2.) +
	 vvely*sin((2*Pi*y*pow(L,-1))/3.));
	
	// energy source
	source[3] = (uvel0 + uvely*cos((3*Pi*y*pow(L,-1))/5.) +
				 uvelx*sin((3*Pi*x*pow(L,-1))/2.))*
	(-2*Pi*pressx*pow(L,-1)*sin(2*Pi*x*pow(L,-1)) +
	 (rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
	  rhox*sin(Pi*x*pow(L,-1)))*
	 (-2*Pi*pressx*pow(-1 + gamma,-1)*pow(L,-1)*
	  pow(rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
		  rhox*sin(Pi*x*pow(L,-1)),-1)*sin(2*Pi*x*pow(L,-1)) +
	  (3*Pi*uvelx*cos((3*Pi*x*pow(L,-1))/2.)*pow(L,-1)*
	   (uvel0 + uvely*cos((3*Pi*y*pow(L,-1))/5.) +
		uvelx*sin((3*Pi*x*pow(L,-1))/2.)) -
	   Pi*vvelx*pow(L,-1)*sin((Pi*x*pow(L,-1))/2.)*
	   (vvel0 + vvelx*cos((Pi*x*pow(L,-1))/2.) +
		vvely*sin((2*Pi*y*pow(L,-1))/3.)))/2. -
	  Pi*rhox*cos(Pi*x*pow(L,-1))*pow(-1 + gamma,-1)*pow(L,-1)*
	  pow(rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
		  rhox*sin(Pi*x*pow(L,-1)),-2)*
	  (press0 + pressx*cos(2*Pi*x*pow(L,-1)) +
	   pressy*sin(Pi*y*pow(L,-1)))) +
	 Pi*rhox*cos(Pi*x*pow(L,-1))*pow(L,-1)*
	 ((pow(uvel0 + uvely*cos((3*Pi*y*pow(L,-1))/5.) +
		   uvelx*sin((3*Pi*x*pow(L,-1))/2.),2) +
	   pow(vvel0 + vvelx*cos((Pi*x*pow(L,-1))/2.) +
		   vvely*sin((2*Pi*y*pow(L,-1))/3.),2))/2. +
	  pow(-1 + gamma,-1)*
	  pow(rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
		  rhox*sin(Pi*x*pow(L,-1)),-1)*
	  (press0 + pressx*cos(2*Pi*x*pow(L,-1)) +
	   pressy*sin(Pi*y*pow(L,-1))))) +
	(3*Pi*uvelx*cos((3*Pi*x*pow(L,-1))/2.)*pow(L,-1)*
	 (press0 + pressx*cos(2*Pi*x*pow(L,-1)) +
	  pressy*sin(Pi*y*pow(L,-1)) +
	  (rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
	   rhox*sin(Pi*x*pow(L,-1)))*
	  ((pow(uvel0 + uvely*cos((3*Pi*y*pow(L,-1))/5.) +
			uvelx*sin((3*Pi*x*pow(L,-1))/2.),2) +
		pow(vvel0 + vvelx*cos((Pi*x*pow(L,-1))/2.) +
			vvely*sin((2*Pi*y*pow(L,-1))/3.),2))/2. +
	   pow(-1 + gamma,-1)*
	   pow(rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
		   rhox*sin(Pi*x*pow(L,-1)),-1)*
	   (press0 + pressx*cos(2*Pi*x*pow(L,-1)) +
		pressy*sin(Pi*y*pow(L,-1))))))/2. +
	(2*Pi*vvely*cos((2*Pi*y*pow(L,-1))/3.)*pow(L,-1)*
	 (press0 + pressx*cos(2*Pi*x*pow(L,-1)) +
	  pressy*sin(Pi*y*pow(L,-1)) +
	  (rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
	   rhox*sin(Pi*x*pow(L,-1)))*
	  ((pow(uvel0 + uvely*cos((3*Pi*y*pow(L,-1))/5.) +
			uvelx*sin((3*Pi*x*pow(L,-1))/2.),2) +
		pow(vvel0 + vvelx*cos((Pi*x*pow(L,-1))/2.) +
			vvely*sin((2*Pi*y*pow(L,-1))/3.),2))/2. +
	   pow(-1 + gamma,-1)*
	   pow(rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
		   rhox*sin(Pi*x*pow(L,-1)),-1)*
	   (press0 + pressx*cos(2*Pi*x*pow(L,-1)) +
		pressy*sin(Pi*y*pow(L,-1))))))/3. +
	(vvel0 + vvelx*cos((Pi*x*pow(L,-1))/2.) +
	 vvely*sin((2*Pi*y*pow(L,-1))/3.))*
	(Pi*pressy*cos(Pi*y*pow(L,-1))*pow(L,-1) -
	 (Pi*rhoy*pow(L,-1)*sin((Pi*y*pow(L,-1))/2.)*
	  ((pow(uvel0 + uvely*cos((3*Pi*y*pow(L,-1))/5.) +
			uvelx*sin((3*Pi*x*pow(L,-1))/2.),2) +
		pow(vvel0 + vvelx*cos((Pi*x*pow(L,-1))/2.) +
			vvely*sin((2*Pi*y*pow(L,-1))/3.),2))/2. +
	   pow(-1 + gamma,-1)*
	   pow(rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
		   rhox*sin(Pi*x*pow(L,-1)),-1)*
	   (press0 + pressx*cos(2*Pi*x*pow(L,-1)) +
		pressy*sin(Pi*y*pow(L,-1)))))/2. +
	 (rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
	  rhox*sin(Pi*x*pow(L,-1)))*
	 (Pi*pressy*cos(Pi*y*pow(L,-1))*pow(-1 + gamma,-1)*
	  pow(L,-1)*pow(rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
					rhox*sin(Pi*x*pow(L,-1)),-1) +
	  ((-6*Pi*uvely*pow(L,-1)*
		(uvel0 + uvely*cos((3*Pi*y*pow(L,-1))/5.) +
		 uvelx*sin((3*Pi*x*pow(L,-1))/2.))*
		sin((3*Pi*y*pow(L,-1))/5.))/5. +
	   (4*Pi*vvely*cos((2*Pi*y*pow(L,-1))/3.)*pow(L,-1)*
		(vvel0 + vvelx*cos((Pi*x*pow(L,-1))/2.) +
		 vvely*sin((2*Pi*y*pow(L,-1))/3.)))/3.)/2. +
	  (Pi*rhoy*pow(-1 + gamma,-1)*pow(L,-1)*
	   pow(rho0 + rhoy*cos((Pi*y*pow(L,-1))/2.) +
		   rhox*sin(Pi*x*pow(L,-1)),-2)*
	   sin((Pi*y*pow(L,-1))/2.)*
	   (press0 + pressx*cos(2*Pi*x*pow(L,-1)) +
		pressy*sin(Pi*y*pow(L,-1))))/2.));
}


void MMS_domain_and_ghost(vector<vector<array<double,2> > >& cell_and_ghost_ctr,\
						 const vector<vector<array<double,2> > >& cell_ctr,    \
						 const int& CI, const int& CJ)
{
	// load central cell locations
	for (int i = NG; i < (NG+CI); ++i)
	{
		for (int j = NG; j < (NG+CJ); ++j)
		{
			cell_and_ghost_ctr[i][j] = cell_ctr[i-NG][j-NG];
		}
	}
	
	// extrapolate xi-lower positions
	for (int i = (NG-1); i >= 0; --i)
	{
		for (int j = NG; j < (NG+CJ); ++j)
		{
			cell_and_ghost_ctr[i][j][0] = 2.0*cell_and_ghost_ctr[i+1][j][0] - cell_and_ghost_ctr[i+2][j][0];
		    cell_and_ghost_ctr[i][j][1] = 2.0*cell_and_ghost_ctr[i+1][j][1] - cell_and_ghost_ctr[i+2][j][1];
		}
	}
	
	// extrapolate xi-upper positions
	for (int i = (NG+CI); i < (2*NG+CI); ++i)
	{
		for (int j = NG; j < (NG+CJ); ++j)
		{
			cell_and_ghost_ctr[i][j][0] = 2.0*cell_and_ghost_ctr[i-1][j][0] - cell_and_ghost_ctr[i-2][j][0];
			cell_and_ghost_ctr[i][j][1] = 2.0*cell_and_ghost_ctr[i-1][j][1] - cell_and_ghost_ctr[i-2][j][1];
		}
	}

	// extrapolate eta-lower positions
	for (int i = NG; i < (NG+CI); ++i)
	{
		for (int j = (NG-1); j >= 0; --j)
		{
			cell_and_ghost_ctr[i][j][0] = 2.0*cell_and_ghost_ctr[i][j+1][0] - cell_and_ghost_ctr[i][j+2][0];
			cell_and_ghost_ctr[i][j][1] = 2.0*cell_and_ghost_ctr[i][j+1][1] - cell_and_ghost_ctr[i][j+2][1];
		}
	}

	// extrapolate eta-upper positions
	for (int i = NG; i < (NG+CI); ++i)
	{
		for (int j = (NG+CJ); j < (2*NG+CJ); ++j)
		{
			cell_and_ghost_ctr[i][j][0] = 2.0*cell_and_ghost_ctr[i][j-1][0] - cell_and_ghost_ctr[i][j-2][0];
			cell_and_ghost_ctr[i][j][1] = 2.0*cell_and_ghost_ctr[i][j-1][1] - cell_and_ghost_ctr[i][j-2][1];
		}
	}
}















