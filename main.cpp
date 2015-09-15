//
//  main.cpp
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 4/30/15.
//  
//


// standard includes
#include <iostream>
#include <cmath>
#include <fstream>
#include <array>
#include <string>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>


// my includes
#include "main.h"
#include "Geometry.h"
#include "Time_Step.h"
#include "Roe_Flux.h"
#include "Van_Leer_Flux.h"
#include "simulation_io.h"
#include "Euler_Helper_Functions.h"
#include "mms.h"
#include "muscl.h"
#include "boundary_conditions.h"
#include "diagnostics.h"

// Will always have two ghost cells.
const int NG = 2;

using namespace std;

int main()
{

    
    // start runtime clock
    clock_t start_time;
    start_time = clock();
    
    /*--------------------------DECLARE-PARAMETERS----------------------------*/

    string sim_name;       // name given to simulation case executed
    string mesh_file_path; // name of file containing mesh
    char sim_case;         // flag to set what type of simulation is being peformed.
    char flux_method;      // flux method: 'v' == Van Leer, 'r' == Roe
    int sub_super;         // sub sonic or super sonic for mms case.
    double angle_of_attack;// angle of attack with respect to free stream.
    double mach;           // mach of free stream flow.
    double temperature;    // temperature of free stream flow.
    double pressure;       // pressure of free stream flow.
    double CFL;            // stability criteria to be used
    int mstage;            // multi-stage time advance to be used. (1, 2, 4)
    int order;             // spacial order of accuracy to be used (1 or 2)
    double kappa;          // degree of upwinding.
    int first_till;        // number of iterations to simulate at first order
    int freeze_at;         // iteration after which to stop updating flux limiters
    double converge_till;  // tolerance that solution should be converged to.
    double blow_up;        // halt solution if residual grows beyond this point.
    int max_iter;          // maximum number of iterations to execute
    int file_out;          // output desired data every (#) iterations.
    int terminal_out;      // output status to the terminal every (#) iterations.
    string tecplot_file;   // name of file in which to write tecplot data.
    string go_file;        // name of file in which to write global data.
 
    /*--------------------------READ-IN-CONFIG-FILE---------------------------*/

    read_config(sim_name, mesh_file_path, sim_case, flux_method, sub_super,  \
                angle_of_attack, mach, temperature, pressure, CFL, mstage,   \
                order, kappa, first_till, freeze_at, converge_till, blow_up, \
                max_iter, file_out, terminal_out, tecplot_file, go_file);

    /*----------------------DECLARE-SIMULAITON-VARIABLES----------------------*/
    
    // Node location
    vector< vector< array<double,2> > > node_loc(1,vector< array<double,2> >(1));
    int NI;  // nodes in the i-direction
    int NJ;  // nodes in the j-direction
    
    // gemometric variables
    vector< vector<double> >  xi_A(1,vector<double>(1));
    vector< vector<double> > eta_A(1,vector<double>(1));
    
    vector< vector< array<double,2> > >  xi_N(1,vector< array<double,2> >(1));
    vector< vector< array<double,2> > > eta_N(1,vector< array<double,2> >(1));
    
    vector< vector<double> > cell_V(1,vector<double>(1));
    vector< vector< array<double,2> > > cell_ctr(1,vector< array<double,2> >(1));
    
    // physical variables
    vector< vector< array<double,4> > >  xi_pos_limiters(1,vector< array<double,4> >(1));
    vector< vector< array<double,4> > >  xi_neg_limiters(1,vector< array<double,4> >(1));
    vector< vector< array<double,4> > > eta_pos_limiters(1,vector< array<double,4> >(1));
    vector< vector< array<double,4> > > eta_neg_limiters(1,vector< array<double,4> >(1));
    
    vector< vector< array<double,4> > >  xi_Flux(1,vector< array<double,4> >(1));
    vector< vector< array<double,4> > > eta_Flux(1,vector< array<double,4> >(1));
    
    vector< vector< array<double,4> > >  primVar(1,vector< array<double,4> >(1));
    vector< vector< array<double,4> > >   source(1,vector< array<double,4> >(1));
    vector< vector< array<double,4> > > residual(1,vector< array<double,4> >(1));
    
    vector< vector<double> > dt(1, vector<double>(1));
    
    array<double,4> inflow_condition;
    
    // Iterative residual storage variables
    //array<double,4> l1_norm;
    array<double,4> l2_norm;
    //array<double,4> linf_norm;
    
    /*---------------------------LOAD-MESH------------------------------------*/

    load_mesh(mesh_file_path, node_loc, NI, NJ);
 
    /*--------------------------CALCULATE-GEOMETRY----------------------------*/

    Set_Geometry(xi_A, eta_A, xi_N, eta_N, cell_V, cell_ctr, node_loc, NI, NJ);
    int CI = cell_V.size();
    int CJ = cell_V[0].size();

    /*--------------------------RESIZE-VARIABLES------------------------------*/

    resize_vars(NI, NJ,                             \
                xi_pos_limiters, xi_neg_limiters,   \
                eta_pos_limiters, eta_neg_limiters, \
                xi_Flux,    eta_Flux,               \
                primVar, source, residual, dt);

    /*------------------INITALIZE-SOURCE-AND-MMS-AS-NEEDED--------------------*/
    
    mms_and_source_init(sim_case, sub_super, CI, CJ, cell_ctr, primVar, source);

    /*--------------------------SETUP-OUTPUT-FILES----------------------------*/
    
    setup_output(sim_name, tecplot_file, go_file);
    tecplot_output(tecplot_file, 0, node_loc, primVar, residual);

    /*-------------------------MAIN-LOOP-OVER-TIME----------------------------*/

    // initialize misc. variables.
    bool loop_ctrl = true;
    int iter = 1;
    double error_norm = 1.0;
    array<double,4> l2_initial;
    int current_order = 1;
    
    // initialize xi limiters
    for (int i = 0; i < CI+2; ++i)
    {
        for (int j = 0; j < CJ; ++j)
        {
            for (int n = 0; n < 4; ++n)
            {
                xi_pos_limiters[i][j][n]  = 1.0;
                xi_neg_limiters[i][j][n]  = 1.0;
            }
        }
    }
   
    // initialize eta limiters
    for (int i = 0; i < CI; ++i)
    {
        for (int j = 0; j < CJ+2; ++j)
        {
            for (int n = 0; n < 4; ++n)
            {
                eta_pos_limiters[i][j][n]  = 1.0;
                eta_neg_limiters[i][j][n]  = 1.0;
                
            }
        }
    }
    
    int edge = 0;

    // set initial condition (INLET CASE)
    if (sim_case == 'i')
    {
        
        // calculate the primitive form of the input variables,
        // store in inflow_condition.
        input_to_prim(inflow_condition, mach, pressure, temperature, angle_of_attack);
        
        // initialize domain with the inflow conditon.
        for (int i = 0; i < CI+2*NG; ++i)
        {
            for (int j = 0; j < CJ+2*NG; ++j)
            {
                for (int n = 0; n < 4; ++n)
                {
                    primVar[i][j][n] = inflow_condition[n];
                }
            }
        }
        
        // find the proper location for the edge in this case.
        while (eta_N[edge][CJ][0] != 0.0)
        {
            ++edge;
        }
    }
    
    // set initial condition (AIRFOIL CASE)
    if (sim_case == 'a')
    {
        
        // calculate the primitive form of the input variables,
        // store in inflow_condition.
        input_to_prim(inflow_condition, mach, pressure, temperature, angle_of_attack);
        
        // initialize domain with the inflow conditon.
        for (int i = 0; i < CI+2*NG; ++i)
        {
            for (int j = 0; j < CJ+2*NG; ++j)
            {
                for (int n = 0; n < 4; ++n)
                {
                    primVar[i][j][n] = inflow_condition[n];
                }
            }
        }
 
        while (eta_N[edge][0][0] == 0.0)
        {
            ++edge;
        }
    }
    
    // continue to loop while loop_controller == TRUE
	while ( loop_ctrl )
	{
        // calculate local timesteps
        Time_Step(dt, CFL, cell_V, xi_A, eta_A, xi_N, eta_N, primVar);
        
        // apply boundary conditions
        update_BC(sim_case, edge, CI, CJ, xi_N, eta_N, inflow_condition, primVar);
        
        
        /*--------------------perform-explicit-time-advance-------------------*/
        // decide if first order
        if (iter < first_till)
        {
            current_order = 1;
        }
        else
        {
            current_order = order;
            // calculate flux limiters
            if (iter < freeze_at)
            {
                update_limiters(xi_pos_limiters,  xi_neg_limiters,  \
                                eta_pos_limiters, eta_neg_limiters, \
                                primVar, CI, CJ);
            }
        }
        
        for (int alpha = mstage; alpha > 0; --alpha)
        {
            // calculate fluxes
            update_fluxes(xi_Flux, eta_Flux,                  \
                          flux_method, current_order, kappa,  \
                          xi_pos_limiters, xi_neg_limiters,   \
                          eta_pos_limiters, eta_neg_limiters, \
                          primVar, xi_N, eta_N, CI, CJ);
        
            // calculate residuals
            update_residuals(residual, xi_Flux, eta_Flux,     \
                             source, cell_V, xi_A, eta_A,     \
                             CI, CJ);
        
            // update solution
            update_solution(primVar,dt,cell_V,residual,alpha,CI,CJ);
        }
        /*------------------------output-data---------------------------------*/
        if (iter % file_out == 0)
        {
            if ( isnan(error_norm) )
            {
                cerr << "error_norm is not a number. ending simulation" << endl;
                exit(1);
            }
            tecplot_output(tecplot_file, iter, node_loc, primVar, residual);
            //global_output(go_file, iter, l1_norm, l2_norm, linf_norm);
        }
        
        /*----------------MONITOR-SOLUTION-CONVERGENCE------------------------*/
        if (iter == 1)
        {
            // store first residual for normalization
            residual_norms(l2_initial,CI,CJ,residual);
        }
        
        if (iter % terminal_out == 0)
        {
            // clalculate residual norms
            residual_norms(l2_norm,CI,CJ,residual);
            // read in chosen diagnostic to the loop controller
            error_norm = l2_norm[0]/l2_initial[0];
            if ( isnan(error_norm) )
            {
                cerr << "error_norm is not a number. ending simulation" << endl;
                exit(1);
            }
            
            terminal_output(sim_name,start_time,iter,error_norm);
            
            // update the loop controller
            loop_ctrl=(iter<max_iter)&&(error_norm>converge_till)&&(error_norm<blow_up);
        }
        // iterate
        iter += 1;
	}
    
    /*------------------------END-LOOP-OVER-TIME------------------------------*/
    
    // final output
    sim_end_output(sim_name,tecplot_file,go_file,start_time,iter);
    
    switch (sim_case)
    {
        case 'm':
        {
            mms_diagnostics(sim_name, iter, CI, CJ, sub_super, cell_ctr, \
                            node_loc, primVar);
            break;
        }
        
        case 'i':
        {
            inlet_diagnostics(CI, CJ, inflow_condition, primVar);
            break;
        }
            
        case 'a':
        {
            airfoil_diagnostics(CI, CJ, edge, angle_of_attack, sim_name, node_loc, \
                                eta_N, eta_A, inflow_condition, primVar);
            break;
        }
        default:
        {
            cerr << "  " << sim_case << " is not a valid simulation case flag" << endl;
            break;
        }
    }
    
    // end main
	return 0;
}

void resize_vars(const int& NI, const int& NJ,                           \
                 vector<vector<array<double,4> > >& xi_pos_limiters,     \
                 vector<vector<array<double,4> > >& xi_neg_limiters,     \
                 vector<vector<array<double,4> > >& eta_pos_limiters,    \
                 vector<vector<array<double,4> > >& eta_neg_limiters,    \
                 vector<vector<array<double,4> > >& xi_Flux,             \
                 vector<vector<array<double,4> > >& eta_Flux,            \
                 vector<vector<array<double,4> > >& primVar,             \
                 vector<vector<array<double,4> > >& source,              \
                 vector<vector<array<double,4> > >& residual,            \
                 vector<vector<double> >& dt)
{
    // resize limiters
    xi_pos_limiters.resize(NI+1);
    xi_neg_limiters.resize(NI+1);
    for (int i = 0; i < (NI+1); ++i)
    {
        xi_pos_limiters[i].resize(NJ-1);
        xi_neg_limiters[i].resize(NJ-1);
    }
    
    eta_pos_limiters.resize(NI-1);
    eta_neg_limiters.resize(NI-1);
    for (int i = 0; i < (NI-1); ++i)
    {
        eta_pos_limiters[i].resize(NJ+1);
        eta_neg_limiters[i].resize(NJ+1);
    }
    
    
    // resize fluxes
    xi_Flux.resize(NI);
    for (int i = 0; i < NI; ++i)
    {
        xi_Flux[i].resize(NJ-1);
    }

    eta_Flux.resize(NI-1);
    for (int i = 0; i < (NI-1); ++i)
    {
        eta_Flux[i].resize(NJ);
    }
    
    // resize variables located at cell centers
    source.resize(NI-1);
    residual.resize(NI-1);
    dt.resize(NI-1);
    for (int i = 0; i < (NI-1); ++i)
    {
        source[i].resize(NJ-1);
        residual[i].resize(NJ-1);
        dt[i].resize(NJ-1);
    }
    
    // resize primitive variables (must include ghost cells
    primVar.resize(NI-1+2*NG);
    for (int i = 0; i < (NI-1+2*NG); ++i)
    {
        primVar[i].resize(NJ-1+2*NG);
    }
}


void mms_and_source_init(const char& sim_case, const int& sub_super,    \
                    const int& CI, const int& CJ,                       \
                    const vector<vector<array<double,2> > >& cell_ctr,  \
                    vector<vector<array<double,4> > >& primVar,         \
                    vector<vector<array<double,4> > >& source)
{
    // SPECIFIC TO MMS initialize primitive_variables
    if (sim_case == 'm')
    {
        vector<vector<array<double,2> > > cell_and_ghost_ctr((2*NG+CI),vector<array<double,2> >(2*NG+CJ));
        
        
        MMS_domain_and_ghost(cell_and_ghost_ctr, cell_ctr, CI, CJ);
        
        MMS_initial_and_boundary(sub_super, primVar, cell_and_ghost_ctr, CI, CJ);
        
        // SPECIFIC TO MMS calculate source terms
        for (int i = 0; i < CI; ++i)
        {
            for (int j = 0; j < CJ; ++j)
            {
                MMS_source(sub_super, source[i][j], cell_ctr[i][j]);
            }
        }
    }
    
    else
    {
        for (int i = 0; i < CI; ++i)
        {
            for (int j = 0; j < CJ; ++j)
            {
                for (int n = 0; n < 4; ++n)
                {
                    source[i][j][n] = 0.0;
                }
            }
        }
    }
}


void update_BC(const char& sim_case, const int& edge,          \
               const int& CI, const int& CJ,                   \
               const vector<vector<array<double,2> > >& xi_N,  \
               const vector<vector<array<double,2> > >& eta_N, \
               const array<double,4>& inflow_condition,        \
               vector<vector<array<double,4> > >& primVar)
{
    // take simulation case flag
    switch (sim_case)
    {
        case 'm':
        {
            // no Boundary conditon update needed for MMS.
            break;
        }
            
        case 'i':
        {
            // eta-boundary conditons
            for (int i = NG; i < (NG+CI); ++i)
            {
                // eta-min-bc
                slipWall_bc(primVar[i][1],  \
                            primVar[i][0],  \
                            primVar[i][2],  \
                            primVar[i][3],  \
                            eta_N[i-NG][0]);
                
                // eta-max-bc
                if (i < (edge + NG) )
                {
                    super_inflow_bc(primVar[i][NG+CJ],   \
                                    primVar[i][NG+CJ+1], \
                                    inflow_condition);
                }
                else
                {
                    slipWall_bc(primVar[i][NG+CJ],   \
                                primVar[i][NG+CJ+1], \
                                primVar[i][NG+CJ-1], \
                                primVar[i][NG+CJ-2], \
                                eta_N[i-NG][CJ]);
                }
            }
            
            // xi-boundary condition
            for (int j = NG; j < (NG+CJ); ++j)
            {
                // xi-min-bc
                super_inflow_bc(primVar[1][j], \
                                primVar[0][j], \
                                inflow_condition);
                
                // xi-max-bc
                super_outflow_bc(primVar[NG+CI][j],   \
                                 primVar[NG+CI+1][j], \
                                 primVar[NG+CI-1][j], \
                                 primVar[NG+CI-2][j]);
            }
            break;
        }

        case 'a':
        {
            // eta boundary conditions
            for (int i = NG; i < (NG+CI); ++i)
            {
               // eta-min-bc
                if (i > (edge + NG) && i < (CI - edge + NG -1))
                {
                    slipWall_bc(primVar[i][1],  \
                                primVar[i][0],  \
                                primVar[i][2],  \
                                primVar[i][3],  \
                                eta_N[i-NG][0]);
                }
                else
                {
                    periodic_bc(primVar[i][1],          \
                                primVar[i][0],          \
                                primVar[CI+3-i][2],  \
                                primVar[CI+3-i][3]);
                }
                
                // eta-max-bc
                freeStream_bc(primVar[i][NG+CJ],   \
                              primVar[i][NG+CJ+1], \
                              inflow_condition);
            }
            
            // xi boundary conditions
            for (int j = NG; j < (NG+CJ); ++j)
            {
                // xi-min-bc
                freeStream_bc(primVar[1][j], \
                              primVar[0][j], \
                              inflow_condition);
                
                // xi-max-bc
                freeStream_bc(primVar[NG+CI][j],   \
                              primVar[NG+CI+1][j], \
                              inflow_condition);
            }
            break;
        }
            
        default:
        {
            cerr << "Invalid sim_case flag passed" << endl;
            exit(1);
            break;
        }

    }
}

void update_limiters(vector<vector<array<double,4> > >& xi_pos_limiters,     \
                        vector<vector<array<double,4> > >& xi_neg_limiters,  \
                        vector<vector<array<double,4> > >& eta_pos_limiters, \
                        vector<vector<array<double,4> > >& eta_neg_limiters, \
                        const vector<vector<array<double,4> > >& primVar,    \
                        const int& CI, const int& CJ)
{
    // calculate xi-face limiters
    for (int i = 0; i < (CI+2); ++i)
    {
        for (int j = 0; j < CJ; ++j)
        {
            calculate_limiters('+', xi_pos_limiters[i][j], \
                               primVar[i+NG-2][j+NG],      \
                               primVar[i+NG-1][j+NG],      \
                               primVar[i+NG][j+NG]);
            
            calculate_limiters('-', xi_neg_limiters[i][j], \
                               primVar[i+NG-1][j+NG],      \
                               primVar[i+NG][j+NG],        \
                               primVar[i+NG-2][j+NG]);
        }
    }
    
    // calculate eta-face limiters
    for (int i = 0; i < CI; ++i)
    {
        for (int j = 0; j < (CJ+2); ++j)
        {
            calculate_limiters('+', xi_pos_limiters[i][j], \
                               primVar[i+NG][j+NG-2],      \
                               primVar[i+NG][j+NG-1],      \
                               primVar[i+NG][j+NG]);
            
            calculate_limiters('-', xi_neg_limiters[i][j], \
                               primVar[i+NG][j+NG-1],      \
                               primVar[i+NG][j+NG],        \
                               primVar[i+NG][j+NG-2]);
        }
    }
}

void update_fluxes(vector<vector<array<double,4> > >& xi_Flux,                \
                   vector<vector<array<double,4> > >& eta_Flux,               \
                   const char& method, const int& order, const double& kappa, \
                   const vector<vector<array<double,4> > >& xi_pos_limiters,  \
                   const vector<vector<array<double,4> > >& xi_neg_limiters,  \
                   const vector<vector<array<double,4> > >& eta_pos_limiters, \
                   const vector<vector<array<double,4> > >& eta_neg_limiters, \
                   const vector<vector<array<double,4> > >& primVar,          \
                   const vector<vector<array<double,2> > >& xi_N,             \
                   const vector<vector<array<double,2> > >& eta_N,            \
                   const int& CI, const int& CJ)
{
    // initialize left and right variable arrays for use in flux calculation.
    array<double,4> V_L;
    array<double,4> V_R;
    switch (method)
    {
        case 'v':
        {
            // update xi-face fluxes
            for (int i = 0; i < (CI+1); ++i)
            {
                for (int j = 0; j < CJ; ++j)
                {
                    muscl_interpolate(V_L, V_R, order, kappa,  \
                                      xi_pos_limiters[i][j],   \
                                      xi_neg_limiters[i][j],   \
                                      xi_pos_limiters[i+1][j], \
                                      xi_neg_limiters[i+1][j], \
                                      primVar[i+NG-2][j+NG],   \
                                      primVar[i+NG-1][j+NG],   \
                                      primVar[i+NG][j+NG],     \
                                      primVar[i+NG+1][j+NG]);
                    
                    // call Van Leer to calculate flux
                    Van_Leer_Flux(xi_Flux[i][j],V_L,V_R,xi_N[i][j]);
                }
            }
            
            // update eta-face fluxes
            for (int i = 0; i < CI; ++i)
            {
                for (int j = 0; j < (CJ+1); ++j)
                {
                    muscl_interpolate(V_L, V_R, order, kappa,  \
                                      xi_pos_limiters[i][j],   \
                                      xi_neg_limiters[i][j],   \
                                      xi_pos_limiters[i][j+1], \
                                      xi_neg_limiters[i][j+1], \
                                      primVar[i+NG][j+NG-2],   \
                                      primVar[i+NG][j+NG-1],   \
                                      primVar[i+NG][j+NG],     \
                                      primVar[i+NG][j+NG+1]);
                    
                    // call Van Leer to calculate flux
                    Van_Leer_Flux(eta_Flux[i][j],V_L,V_R,eta_N[i][j]);
                }
            }
            
        }
            break;
            
        case 'r':
        {
            // update xi-face fluxes
            for (int i = 0; i < (CI+1); ++i)
            {
                for (int j = 0; j < CJ; ++j)
                {
                    
                    muscl_interpolate(V_L, V_R, order, kappa,  \
                                      xi_pos_limiters[i][j],   \
                                      xi_neg_limiters[i][j],   \
                                      xi_pos_limiters[i+1][j], \
                                      xi_neg_limiters[i+1][j], \
                                      primVar[i+NG-2][j+NG],   \
                                      primVar[i+NG-1][j+NG],   \
                                      primVar[i+NG][j+NG],     \
                                      primVar[i+NG+1][j+NG]);
                    
                    // call Roe to calculate flux
                    Roe_Flux(xi_Flux[i][j],V_L,V_R,xi_N[i][j]);
                }
            }
            
            // update eta-face fluxes
            for (int i = 0; i < CI; ++i)
            {
                for (int j = 0; j < (CJ+1); ++j)
                {
                    muscl_interpolate(V_L, V_R, order, kappa,  \
                                      xi_pos_limiters[i][j],   \
                                      xi_neg_limiters[i][j],   \
                                      xi_pos_limiters[i][j+1], \
                                      xi_neg_limiters[i][j+1], \
                                      primVar[i+NG][j+NG-2],   \
                                      primVar[i+NG][j+NG-1],   \
                                      primVar[i+NG][j+NG],     \
                                      primVar[i+NG][j+NG+1]);
                    
                    // call Roe to calculate flux
                    Roe_Flux(eta_Flux[i][j],V_L,V_R,eta_N[i][j]);
                }
            }
        break;
        }
        default:
        {
            cerr << "Invalid choice for flux option passed." << endl;
            exit(1);
        break;
        }
    }
}


void update_residuals(vector<vector<array<double,4> > >& residual,       \
                      const vector<vector<array<double,4> > >& xi_Flux,  \
                      const vector<vector<array<double,4> > >& eta_Flux, \
                      const vector<vector<array<double,4> > >& source,   \
                      const vector<vector<double> >& cell_V,             \
                      const vector<vector<double> >&  xi_A,              \
                      const vector<vector<double> >& eta_A,              \
                      const int& CI, const int& CJ)
{
    // setup placeholder variables.
    double xi_in;
    double xi_out;
    double eta_in;
    double eta_out;
    double src;
    
    // update all residuals
    for (int i = 0; i < CI; ++i)
    {
        for (int j = 0; j < CJ; ++j)
        {
            for (int n = 0; n < 4; ++n)
            {
                xi_in = xi_Flux[i][j][n] * xi_A[i][j];
                xi_out = xi_Flux[i+1][j][n] * xi_A[i+1][j];
        
                eta_in = eta_Flux[i][j][n] * eta_A[i][j];
                eta_out = eta_Flux[i][j+1][n] * eta_A[i][j+1];
                
                src = source[i][j][n]*cell_V[i][j];
                
                residual[i][j][n] = (xi_out - xi_in + eta_out - eta_in) - src;
            }

        }
    }
}


void update_solution(vector<vector<array<double,4> > >& primVar,         \
                     const vector<vector<double> >& dt,                  \
                     const vector<vector<double> >& cell_V,              \
                     const vector<vector<array<double,4> > >& residual,  \
                     const int& alpha, const int& CI, const int& CJ)
{
    // set up needed temp variables.
    array<double,4> temp_prim;
    array<double,4> temp_cons;
    
    // update solution at every cell location
    for (int i = 0; i < CI; ++i)
    {
        for (int j = 0; j < CJ; ++j)
        {
            // read in current primitive variable to temp_prim
            temp_prim = primVar[i+NG][j+NG];
            // convert to conserved, writting to temp_cons
            primitive_to_conserved(temp_cons,temp_prim);
            
            // update value of temp_conserved
            for (int n = 0; n < 4; ++n)
            {
                temp_cons[n] -= dt[i][j]*residual[i][j][n]/(cell_V[i][j]*alpha);
            }
            
            // convert new conserved value back to primitve
            conserved_to_primitive(temp_prim,temp_cons);
            
            // copy new value of temp_prim into primVar array.
            primVar[i+NG][j+NG] = temp_prim;
        }
    }
}


void residual_norms(array<double,4>& l2_norm, const int& CI, const int& CJ,  \
                    const vector<vector<array<double,4> > >& residual)
{
    // claculate l2_norm of the residuals
    array<double,4> sumsquares = {{0.0, 0.0, 0.0, 0.0}};
    
    for (int i = 0; i < CI; ++i)
    {
        for (int j = 0; j < CJ; ++j)
        {
            for (int n = 0; n < 4; ++n)
            {
                sumsquares[n] += residual[i][j][n]*residual[i][j][n];
            }
        }
    }
    
    for (int n = 0; n < 4; ++n)
    {
        l2_norm[n] = sqrt( sumsquares[n]/(CI*CJ) );
    }
}
