//
//  Van_Leer_Flux.h
//  Euler_2D_FVM
//
//  Created by Colin Christopher Glesner on 4/24/15.
//
//

#ifndef __Euler_2D_FVM__Van_Leer_Flux__
#define __Euler_2D_FVM__Van_Leer_Flux__

#include <stdio.h>
#include <array>

using namespace std;

/* 
 * This header file includes all funcitons needed for the calculation of
 * flux using Van Leer's scheeme. These functions are needed for Van Leer
 * Flux calculation only.
 *
 *
 *  As reference:
 *
 *  Primitive Variables:      Flux Variables:
 *  V[0] = density            F[0] = mass flux
 *  V[1] = x-velocity         F[1] = x-momentum flux
 *  V[2] = y-velocity         F[2] = y-momentum flux
 *  V[3] = Pressure           F[3] = energy flux
 *
 *  n[0] = x-normal
 *  n[1] = y-normal
 */

void Van_Leer_Flux(array<double,4>& Fn,          \
                   const array<double,4>& V_L,   \
                   const array<double,4>& V_R,   \
                   const array<double,2>& n);
/* The actual function that calculates flux.
 *  Fn  -> normal flux
 *  V_L -> left primitive values
 *  V_R -> right primitive values
 *  n   -> normals
 */

double C(double mach, char pos_neg);
/* This function is the top level function of the switching mechanism in the
 * Van Leer flux vector splitting scheme. It takes the mach number at a
 * given position and the sign that it shoud have, called 'pos_neg'.
 * 'pos_neg' must be given the value of '+' or '-' only. C's logic
 * diagram is:
 *
 *       | (+) | (-) |
 * ------+-----+-----+
 *  M<=-1|  0  |  M  |
 * ------+-----+-----+
 * |M|<1 |  M+ |  M- |
 * ------+-----+-----+
 *  M>=1 |  M  |  0  |
 * ------+-----+-----+
 */


double ALPHA(double mach, char pos_neg);
/* This function is part of the switching mechanism in the Van Leer flux
 * vector splitting scheme. It takes the mach number at a given position
 * and the sign that it shoud have, called 'pos_neg'. 'pos_neg' must be
 * given the value of '+' or '-' only. ALPHA's logic diagram is:
 *
 *    | M>0 | M<0 |
 * ---+-----+-----+
 * (+)|  1  |  0  |
 * ---+-----+-----+
 * (-)|  0  |  1  |
 * ---+-----+-----+
 */


double BETA(double mach);
/* This function determines the value of the switching function BETA as part
 * of the standard Van Leer flux vector splitting scheme. BETA takes the
 * mach number at a given location and returns the appropriate value.
 * BETA's function is:
 *
 * BETA == 0 when |Mach| >= 1, BETA == -1 when |Mach| < 1
 *
 */


double M_pm(double mach, char pos_neg);
/* This function is part of the switching mechanism in the Van Leer flux
 * vector splitting scheme. It takes the mach number at a given position and
 * the sign that it should have, called 'pos_neg'. 'pos_neg' must be given
 * the value of '+' or '-' only.
 */


double D(double mach, char pos_neg);
/* This function is one of the top level function of the switching mechanism
 * in the Van Leer flux vector splitting scheme. It takes the mach number at
 * a given position and the sign that it shoud have, called 'pos_neg'.
 * 'pos_neg' must be given the value of '+' or '-' only. D's logic diagram
 * is as follows:
 *
 *       |    (+)    |    (-)    |
 * ------+-----------+-----------+
 *  M<=-1|     0     |     1     |
 * ------+-----------+-----------+
 * |M|<1 |  M+(-M+2) |  M-(-M-2) |
 * ------+-----------+-----------+
 *  M>=1 |     M     |     0     |
 * ------+-----------+-----------+
 */


double P_bar(double mach, char pos_neg);
/* This function is part of the switching mechanism in the
 * Van Leer flux vector splitting scheme. It takes the mach number at a
 * given position and the sign that it shoud have, called 'pos_neg'.
 * 'pos_neg' must be given the value of '+' or '-' only.
 */


double SIGN(double A, double B);
/* This function will Change the sign of A to match the sign of B. This
 * funciton will interpret the sign of 0 as being positive.
 */


double INT(double A);
/* This function will truncate the decimal portion of A, returning the
 * integer value only.
 */

#endif /* defined(__Euler_2D_FVM__Van_Leer_Flux__) */
