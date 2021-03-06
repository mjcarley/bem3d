/* gradient.c
 * 
 * Copyright (C) 2012, 2013, 2018 by Michael Carley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <bem3d.h>

#include "polar.h"

static void subtriangle_gradient_polar_quad(gdouble r1, gdouble r2, 
					    gdouble theta,  
					    gdouble z, gint order,
					    gdouble *phi, gdouble *I)

{
  gdouble a, b, b2, z2, k, kd,  t1, t0, S0, S1, C0, C1 ;
  gdouble L1, L0, d1, d0 ;

  a = (r2*cos(theta) - r1)/r2/sin(theta) ; *phi = atan(a) ;

  b2 = (r1*r1 + z*z*(1+a*a))/(1+a*a) ; b = sqrt(b2) ;
  k = fabs(z)/b ; kd = sqrt(1.0-k*k) ;

  g_assert(z != 0.0) ;
  z2 = z*z ;

  t1 = theta+(*phi) ; t0 = (*phi) ;
  S0 = sin(t0) ; S1 = sin(t1) ; C0 = cos(t0) ; C1 = cos(t1) ;
  d1 = sqrt(1-k*k*S1*S1) ; d0 = sqrt(1-k*k*S0*S0) ;
  L1 = log((d1+kd)/(d1-kd)) ; L0 = log((d0+kd)/(d0-kd)) ;

  /*1/R^{3}*/
  I[0]  = -1.0/b*(quad_sin_cos_nmr(0, 1, -1, k, t1, S1, C1) -
		  quad_sin_cos_nmr(0, 1, -1, k, t0, S0, C0)) ;
  I[0] += theta/fabs(z) ;

  /*x_1/R^{3}*/
  I[1]  = -kd*(quad_sin_cos_nmr(0, 1, -1, k, t1, S1, C1) -
		quad_sin_cos_nmr(0, 1, -1, k, t0, S0, C0)) ;
  I[1] += 0.5*(S1*L1 - S0*L0) ;
  I[1] -= kd*(quad_sin_cos_nmr(2, -1, -1, k, t1, S1, C1) -
	      quad_sin_cos_nmr(2, -1, -1, k, t0, S0, C0)) ;

  /*y_1/R^{3}*/
  I[2]  = -0.5*(C1*L1 - C0*L0) ;

  /*x_1^{2}/R^{3}*/
  I[3]  = (quad_sin_cos_nmr(0, 1, 1, k, t1, S1, C1) -
	   quad_sin_cos_nmr(0, 1, 1, k, t0, S0, C0))*b ;
  I[3] += (quad_sin_cos_nmr(0, 3, -1, k, t1, S1, C1) -
	   quad_sin_cos_nmr(0, 3, -1, k, t0, S0, C0))/b*z2 ;
  I[3] -=  (quad_sin_cos(0, 2, t1, S1, C1)-
	    quad_sin_cos(0, 2, t0, S0, C0))*2.0*fabs(z) ;

  /*x_1 y_{1}/R^{3}*/
  I[4]  = (quad_sin_cos_nmr(1, 0, 1, k, t1, S1, C1) -
	   quad_sin_cos_nmr(1, 0, 1, k, t0, S0, C0))*b ;
  I[4] += (quad_sin_cos_nmr(1, 2, -1, k, t1, S1, C1) -
	   quad_sin_cos_nmr(1, 2, -1, k, t0, S0, C0))/b*z2 ;
  I[4] -=  (quad_sin_cos(1, 1, t1, S1, C1)-
	    quad_sin_cos(1, 1, t0, S0, C0))*2.0*fabs(z) ;

  /*y_1^{2}/R^{3}*/
  I[5]  = (quad_sin_cos_nmr(2, -1, 1, k, t1, S1, C1) -
	   quad_sin_cos_nmr(2, -1, 1, k, t0, S0, C0))*b ;
  I[5] += (quad_sin_cos_nmr(2, 1, -1, k, t1, S1, C1) -
	   quad_sin_cos_nmr(2, 1, -1, k, t0, S0, C0))/b*z2 ;
  I[5] -=  (quad_sin_cos(2, 0, t1, S1, C1)-
	    quad_sin_cos(2, 0, t0, S0, C0))*2.0*fabs(z) ;

  /*1/R^{5}*/
  I[6]  = -((quad_sin_cos_nmr(0, 3, -3, k, t1, S1, C1) -
	     quad_sin_cos_nmr(0, 3, -3, k, t0, S0, C0))/b2/b -
	    theta/fabs(z)/z2)/3.0 ;
  /*x_1/R^{5}*/
  I[7]  = (quad_sin_cos_nmr(0, 1, -3, k, t1, S1, C1) -
	   quad_sin_cos_nmr(0, 1, -3, k, t0, S0, C0))*kd*kd*kd/3.0/z2 ;
  /*y_1/R^{5}*/
  I[8]  = (quad_sin_cos_nmr(1, 0, -3, k, t1, S1, C1) -
	   quad_sin_cos_nmr(1, 0, -3, k, t0, S0, C0))*kd*kd*kd/3.0/z2 ;

  /*x_1^{2}/R^{5}*/
  I[9]  = -(quad_sin_cos_nmr(0, 3, -1, k, t1, S1, C1) -
	    quad_sin_cos_nmr(0, 3, -1, k, t0, S0, C0))/b ;
  I[9] += (quad_sin_cos_nmr(0, 5, -3, k, t1, S1, C1) -
	   quad_sin_cos_nmr(0, 5, -3, k, t0, S0, C0))/b/b2/3.*z2 ;
  I[9] +=  (quad_sin_cos(0, 2, t1, S1, C1)-
	    quad_sin_cos(0, 2, t0, S0, C0))*2.0/3.0/fabs(z) ;

  /*x_1 y_1/R^{5}*/
  I[10]  = -(quad_sin_cos_nmr(1, 2, -1, k, t1, S1, C1) -
	     quad_sin_cos_nmr(1, 2, -1, k, t0, S0, C0))/b ;
  I[10] += (quad_sin_cos_nmr(1, 4, -3, k, t1, S1, C1) -
	    quad_sin_cos_nmr(1, 4, -3, k, t0, S0, C0))/b/b2/3.*z2 ;
  I[10] +=  (quad_sin_cos(1, 1, t1, S1, C1)-
	     quad_sin_cos(1, 1, t0, S0, C0))*2.0/3.0/fabs(z) ;

  /*y_1^{2}/R^{5}*/
  I[11]  = -(quad_sin_cos_nmr(2, 1, -1, k, t1, S1, C1) -
	     quad_sin_cos_nmr(2, 1, -1, k, t0, S0, C0))/b ;
  I[11] += (quad_sin_cos_nmr(2, 3, -3, k, t1, S1, C1) -
	    quad_sin_cos_nmr(2, 3, -3, k, t0, S0, C0))/b/b2/3.*z2 ;
  I[11] +=  (quad_sin_cos(2, 0, t1, S1, C1)-
	     quad_sin_cos(2, 0, t0, S0, C0))*2.0/3.0/fabs(z) ;
    
  if ( order == 1 ) return ;

  g_assert_not_reached() ;

  return ;
}

static void subtriangle_gradient_polar_quad_ip(gdouble r1, gdouble r2, 
					       gdouble theta, 
					       gint order,
					       gdouble *phi, gdouble *I)

{
  gdouble a, b, b2, t1, t0, S0, S1, C0, C1 ;

  a = (r2*cos(theta) - r1)/r2/sin(theta) ; *phi = atan(a) ;

  b2 = (r1*r1)/(1+a*a) ; b = sqrt(b2) ;

  t1 = theta+(*phi) ; t0 = (*phi) ;
  S0 = sin(t0) ; S1 = sin(t1) ; C0 = cos(t0) ; C1 = cos(t1) ;

  /*1/R^{3}*/
  I[0] = -(quad_sin_cos(0, 1, t1, S1, C1) -
	   quad_sin_cos(0, 1, t0, S0, C0))/b ;
  /*x_{1}/R^{3}*/
  I[1] = -(S1*log(C1/b) - S0*log(C0/b)) ;
  I[1] -= (quad_sin_cos(2, -1, t1, S1, C1) -
	   quad_sin_cos(2, -1, t0, S0, C0)) ;
  /*y_{1}/R^{3}*/
  I[2]  = C1*log(C1/b) - C0*log(C0/b) ;
  I[2] += (quad_sin_cos(1, 0, t1, S1, C1) -
	   quad_sin_cos(1, 0, t0, S0, C0)) ;

  /*x_1^{2}/R^{3}*/
  I[3] = (quad_sin_cos(0, 1, t1, S1, C1) -
	  quad_sin_cos(0, 1, t0, S0, C0))*b ;
  /*x_1 y_1/R^{3}*/
  I[4] = (quad_sin_cos(1, 0, t1, S1, C1) -
	  quad_sin_cos(1, 0, t0, S0, C0))*b ;
  /*x_1^{2}/R^{3}*/
  I[5] = (quad_sin_cos(2, -1, t1, S1, C1) -
	  quad_sin_cos(2, -1, t0, S0, C0))*b ;
  
  /*1/R^{5}*/
  I[6] = -(quad_sin_cos(0, 3, t1, S1, C1) -
	   quad_sin_cos(0, 3, t0, S0, C0))/3.0/b/b2 ;
  /*x_{1}/R^{5}*/
  I[7] = -(quad_sin_cos(0, 3, t1, S1, C1) -
	   quad_sin_cos(0, 3, t0, S0, C0))/2.0/b2 ;
  /*y_{1}/R^{5}*/
  I[8] = -(quad_sin_cos(1, 2, t1, S1, C1) -
	   quad_sin_cos(1, 2, t0, S0, C0))/2.0/b2 ;

  /*x_1^{2}/R^{5}*/
  I[9] = -(quad_sin_cos(0, 3, t1, S1, C1) -
	   quad_sin_cos(0, 3, t0, S0, C0))/b ;
  /*x_1 y_1/R^{5}*/
  I[10] = -(quad_sin_cos(1, 2, t1, S1, C1) -
	    quad_sin_cos(1, 2, t0, S0, C0))/b ;
  /*x_1^{2}/R^{5}*/
  I[11] = -(quad_sin_cos(2, 1, t1, S1, C1) -
	    quad_sin_cos(2, 1, t0, S0, C0))/b ;
  
  /* if ( order == 1 ) return ; */

  /* g_assert_not_reached() ; */

  return ;
}

static void subtriangle_gradient_quad(gdouble *p,
				      gdouble *x1, gdouble *x2, 
				      gdouble r1, gdouble r2,
				      gdouble ont,
				      gint order,
				      gdouble *I)

{
  gdouble J[16], th, phi, psi, C1, S1, C0, S0 ;
  gdouble IC, IS, IC2, ICS, IS2 ;

  if ( r1 == 0.0 || r2 == 0.0 ) return ;

  if ( (th = ont*acos(((x2[0]-p[0])*(x1[0]-p[0]) + 
		       (x2[1]-p[1])*(x1[1]-p[1]))/r1/r2)) == 0.0 ) return ;

  if ( p[2] != 0.0 ) 
    subtriangle_gradient_polar_quad(r1, r2, th, p[2], order, &phi, J) ;
  else
    subtriangle_gradient_polar_quad_ip(r1, r2, th, order, &phi, J) ;

  C0 = cos(phi) ; S0 = sin(phi) ;
  psi = atan2(x1[1]-p[1], x1[0]-p[0]) ; C1 = cos(psi) ; S1 = sin(psi) ;  

  IC =   J[1]*C0 + J[2]*S0 ;
  IS =   J[2]*C0 - J[1]*S0 ;
  IC2 =  J[3]*C0*C0 + 2*J[4]*C0*S0 + J[5]*S0*S0 ;
  ICS = -J[3]*C0*S0 + J[4]*(C0*C0-S0*S0) + J[5]*S0*C0 ;
  IS2 =  J[5]*C0*C0 - 2*J[4]*C0*S0 + J[3]*S0*S0 ;

  I[0] += J[0] ;
  I[1] += IC*C1 - IS*S1 + p[0]*J[0] ;
  I[2] += IC*S1 + IS*C1 + p[1]*J[0] ;

  I[3] += p[0]*p[0]*J[0] + 2*p[0]*(IC*C1 - IS*S1) +
    IC2*C1*C1 - 2*ICS*C1*S1 + IS2*S1*S1 ;
  I[4] += p[0]*p[1]*J[0] + p[0]*(IS*C1 + IC*S1) + p[1]*(IC*C1 - IS*S1) +
    IC2*S1*C1 + ICS*(C1*C1-S1*S1) - IS2*C1*S1 ;  
  I[5] += p[1]*p[1]*J[0] + 2*p[1]*(IS*C1 + IC*S1) +
    IS2*C1*C1 + 2*ICS*C1*S1 + IC2*S1*S1 ;
  

  IC =   J[7]*C0 + J[8]*S0 ;
  IS =   J[8]*C0 - J[7]*S0 ;
  IC2 =  J[9]*C0*C0 + 2*J[10]*C0*S0 + J[11]*S0*S0 ;
  ICS = -J[9]*C0*S0 + J[10]*(C0*C0-S0*S0) + J[11]*S0*C0 ;
  IS2 =  J[11]*C0*C0 - 2*J[10]*C0*S0 + J[9]*S0*S0 ;

  I[6] += J[6] ; 
  I[7] += IC*C1 - IS*S1 + p[0]*J[6] ;
  I[8] += IC*S1 + IS*C1 + p[1]*J[6] ;

  I[9] += p[0]*p[0]*J[6] + 2*p[0]*(IC*C1 - IS*S1) +
    IC2*C1*C1 - 2*ICS*C1*S1 + IS2*S1*S1 ;
  I[10] += p[0]*p[1]*J[6] + p[0]*(IS*C1 + IC*S1) + p[1]*(IC*C1 - IS*S1) +
    IC2*S1*C1 + ICS*(C1*C1-S1*S1) - IS2*C1*S1 ;  
  I[11] += p[1]*p[1]*J[6] + 2*p[1]*(IS*C1 + IC*S1) +
    IS2*C1*C1 + 2*ICS*C1*S1 + IC2*S1*S1 ;

  if ( order == 1 ) return ;

  return ;
}

static gint triangle_gradient_quad(gdouble *x1, gdouble *x2, gdouble *x3,
				   gdouble *p, 
				   gdouble o12, gdouble o23, gdouble o31,
				   gint order, gboolean hs,
				   gdouble *I)

/*
  evaluate gradient of potential integrals over triangle

  x1, x2, x3: vertices of triangle, all in plane z=0;
  p:          field point for potential, arbitrary position;
  o12, o23, o13: orientation (sign of signed area) of subtriangles (p, x1, x2),
                 (p, x2, x3) and (p, x3, x1);

  I:          pre-allocated array to hold output

  Output: integrals over triangle of:

  [[1 x1 y1 x1^2 x1y1 y1^2]/R^3 
   [1 x1 y1 x1^2 x1y1 y1^2 x1^3 x1^2y1 x1y1^2 y1^3]/R^5]

  Returns 0 on success.
 */

{
  gdouble r1, r2, r3 ;  
  gint i ;

  g_assert(order >= 0 && order <= 1) ;

  for ( i = 0 ; i < 18 ; i ++ ) I[i] = 0.0 ;

  r1 = sqrt((p[0]-x1[0])*(p[0]-x1[0]) +  (p[1]-x1[1])*(p[1]-x1[1])) ;
  r2 = sqrt((p[0]-x2[0])*(p[0]-x2[0]) +  (p[1]-x2[1])*(p[1]-x2[1])) ;
  r3 = sqrt((p[0]-x3[0])*(p[0]-x3[0]) +  (p[1]-x3[1])*(p[1]-x3[1])) ;

  subtriangle_gradient_quad(p, x1, x2, r1, r2, o12, order, I) ;
  subtriangle_gradient_quad(p, x2, x3, r2, r3, o23, order, I) ;
  subtriangle_gradient_quad(p, x3, x1, r3, r1, o31, order, I) ;

  return 0 ;
}

static gint triangle_gradient_quad_shape0(gdouble *x1, gdouble *x2, gdouble *x3,
					  gdouble *p,
					  gdouble o12, gdouble o23, gdouble o31,
					  gdouble *G, gdouble *dG)

/*
  Integration of 1/R and normal derivatives over a triangle, weighted
  with linear shape functions.

  Input parameters as for triangle_potential_quad

  Output:

  G: [d/dx (\int L_{1,2,3}/R) d/dy (\int L_{1,2,3}/R) d/dz (\int L_{1,2,3}/R)]
  dG: dG1/dz (i.e. dG/dn)
*/

{
  gdouble I[18], Ai[9], A[9], J[3] ;

  triangle_gradient_quad(x1, x2, x3, p, o12, o23, o31, 1, TRUE, I) ;

  A[0] = 1.0   ; A[1] = 1.0   ; A[2] = 1.0   ;
  A[3] = x1[0] ; A[4] = x2[0] ; A[5] = x3[0] ;
  A[6] = x1[1] ; A[7] = x2[1] ; A[8] = x3[1] ;

  invert3x3(Ai, A) ;

  /*integrals of d/dx 1/R*/
  J[0] = p[0]*I[0] - I[1] ; J[1] = p[0]*I[1] - I[3] ; J[2] = p[0]*I[2] - I[4] ;
  matvecmul3x3(&(G[0]), Ai, J) ;

  /*integrals of d/dy 1/R*/
  J[0] = p[1]*I[0] - I[2] ; J[1] = p[1]*I[1] - I[4] ; J[2] = p[1]*I[2] - I[5] ;
  matvecmul3x3(&(G[3]), Ai, J) ;

  /*integrals of d/dz 1/R*/
  J[0] = p[2]*I[0] ; J[1] = p[2]*I[1] ; J[2] = p[2]*I[2] ;
  matvecmul3x3(&(G[6]), Ai, J) ;

  /*integrals of d^2/dxdz 1/R*/
  J[0] = p[0]*I[6] - I[7] ; J[1] = p[0]*I[7] - I[9] ; J[2] = p[0]*I[8] - I[10] ;
  matvecmul3x3(&(dG[0]), Ai, J) ;
  dG[0] *= 3*p[2] ; dG[1] *= 3*p[2] ; dG[2] *= 3*p[2] ;

  /*integrals of d^2/dydz 1/R*/
  J[0] = p[1]*I[6] - I[8] ; J[1] = p[1]*I[7] - I[10] ; 
  J[2] = p[1]*I[8] - I[11] ;
  matvecmul3x3(&(dG[3]), Ai, J) ;
  dG[3] *= 3*p[2] ; dG[4] *= 3*p[2] ; dG[5] *= 3*p[2] ;

  /*integrals of d^2/dz^2 1/R*/
  J[0] = I[0] - 3*p[2]*p[2]*I[6] ; 
  J[1] = I[1] - 3*p[2]*p[2]*I[7] ; 
  J[2] = I[2] - 3*p[2]*p[2]*I[8] ; 
  matvecmul3x3(&(dG[6]), Ai, J) ;

  return 0 ;
}

gint triangle_gradient_quad_shape(gdouble *x1, gdouble *x2, gdouble *x3,
				  gdouble *p,
				  gdouble o12, gdouble o23, gdouble o31,
				  gdouble *G, gdouble *dG)

{
  gdouble A[9], y[3], y1[3], y2[3], y3[3], g[9], dg[9] ;
  gdouble tt[3], uu[3] ;

  triangle_rotation_matrix(x1, x2, x3, A) ;

  y[0] = A[0]*(p[0]-x1[0]) + A[1]*(p[1]-x1[1]) + A[2]*(p[2]-x1[2]) ;
  y[1] = A[3]*(p[0]-x1[0]) + A[4]*(p[1]-x1[1]) + A[5]*(p[2]-x1[2]) ;
  y[2] = A[6]*(p[0]-x1[0]) + A[7]*(p[1]-x1[1]) + A[8]*(p[2]-x1[2]) ;

  y2[0] = A[0]*(x2[0]-x1[0]) + A[1]*(x2[1]-x1[1]) + A[2]*(x2[2]-x1[2]) ;
  y2[1] = A[3]*(x2[0]-x1[0]) + A[4]*(x2[1]-x1[1]) + A[5]*(x2[2]-x1[2]) ;
  y2[2] = A[6]*(x2[0]-x1[0]) + A[7]*(x2[1]-x1[1]) + A[8]*(x2[2]-x1[2]) ;

  y3[0] = A[0]*(x3[0]-x1[0]) + A[1]*(x3[1]-x1[1]) + A[2]*(x3[2]-x1[2]) ;
  y3[1] = A[3]*(x3[0]-x1[0]) + A[4]*(x3[1]-x1[1]) + A[5]*(x3[2]-x1[2]) ;
  y3[2] = A[6]*(x3[0]-x1[0]) + A[7]*(x3[1]-x1[1]) + A[8]*(x3[2]-x1[2]) ;

  y1[0] = y1[1] = y1[2] = 0.0 ;

  triangle_orientations(y1, y2, y3, y, &o12, &o23, &o31) ;

  triangle_gradient_quad_shape0(y1, y2, y3, y, o12, o23, o31, g, dg) ;
  
  tt[0] = g[0] ; tt[1] = g[3] ; tt[2] = g[6] ;
  matvecmul3x3trans(uu, A, tt) ;
  G[0] = -uu[0] ; G[3] = -uu[1] ; G[6] = -uu[2] ;

  tt[0] = g[1] ; tt[1] = g[4] ; tt[2] = g[7] ;
  matvecmul3x3trans(uu, A, tt) ;
  G[1] = -uu[0] ; G[4] = -uu[1] ; G[7] = -uu[2] ;

  tt[0] = g[2] ; tt[1] = g[5] ; tt[2] = g[8] ;
  matvecmul3x3trans(uu, A, tt) ;
  G[2] = -uu[0] ; G[5] = -uu[1] ; G[8] = -uu[2] ;

  tt[0] = dg[0] ; tt[1] = dg[3] ; tt[2] = -dg[6] ;
  matvecmul3x3trans(uu, A, tt) ;
  dG[0] = -uu[0] ; dG[3] = -uu[1] ; dG[6] = -uu[2] ;

  tt[0] = dg[1] ; tt[1] = dg[4] ; tt[2] = -dg[7] ;
  matvecmul3x3trans(uu, A, tt) ;
  dG[1] = -uu[0] ; dG[4] = -uu[1] ; dG[7] = -uu[2] ;

  tt[0] = dg[2] ; tt[1] = dg[5] ; tt[2] = -dg[8] ;
  matvecmul3x3trans(uu, A, tt) ;
  dG[2] = -uu[0] ; dG[5] = -uu[1] ; dG[8] = -uu[2] ;

  return 0 ;
}
