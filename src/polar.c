#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <bem3d.h>

#include "polar.h"

gint triangle_quad_shape(gdouble *x1, gdouble *x2, gdouble *x3,
			 gdouble *p,
			 gdouble o12, gdouble o23, gdouble o31,
			 gint order, gboolean hs,
			 gdouble *G1, gdouble *G3, gdouble *G5) ;

static void subtriangle_polar_quad(gdouble r1, gdouble r2, gdouble theta, 
				   gdouble z, gint order, gboolean hs,
				   gdouble *phi, 
				   gdouble *I1, gdouble *I3, gdouble *I5)

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

  I1[0] = quad_sin_cos_nmr(0, -1, 1, k, t1, S1, C1) -
    quad_sin_cos_nmr(0, -1, 1, k, t0, S0, C0) ;
  I1[0] *= b ;
    
  I1[0] -= fabs(z)*theta ;
 
  I3[0]  = -1.0/b*(quad_sin_cos_nmr(0, 1, -1, k, t1, S1, C1) -
		quad_sin_cos_nmr(0, 1, -1, k, t0, S0, C0)) ;
  I3[0] += theta/fabs(z) ;
  if ( hs ) {
    I5[0] = -((quad_sin_cos_nmr(0, 3, -3, k, t1, S1, C1) -
	       quad_sin_cos_nmr(0, 3, -3, k, t0, S0, C0))/b2/b -
	      theta/fabs(z)/z2)/3.0 ;
  }

  if ( order == 0 ) return ;
    
  d1 = sqrt(1-k*k*S1*S1) ; d0 = sqrt(1-k*k*S0*S0) ;
  L1 = log((d1+kd)/(d1-kd)) ; L0 = log((d0+kd)/(d0-kd)) ;

  I1[1]  = 0.5*b2*kd*(quad_sin_cos_nmr(0, -1, 1, k, t1, S1, C1) -
		      quad_sin_cos_nmr(0, -1, 1, k, t0, S0, C0)) ;
  I1[1] += 0.5*z2*kd*(quad_sin_cos_nmr(2, -1, -1, k, t1, S1, C1) -
		      quad_sin_cos_nmr(2, -1, -1, k, t0, S0, C0)) ;
  I1[1] -= 0.25*z2*(L1*S1 - L0*S0) ;

  I1[2]  = 0.5*b2*kd*(quad_sin_cos_nmr(1, -2, 1, k, t1, S1, C1) -
		      quad_sin_cos_nmr(1, -2, 1, k, t0, S0, C0)) ;
  I1[2] -= 0.5*z2*kd*(quad_sin_cos_nmr(1, 0, -1, k, t1, S1, C1) -
		      quad_sin_cos_nmr(1, 0, -1, k, t0, S0, C0)) ;
  I1[2] += 0.25*z2*(L1*C1 - L0*C0) ;  

  I3[1]  = -kd*(quad_sin_cos_nmr(0, 1, -1, k, t1, S1, C1) -
		quad_sin_cos_nmr(0, 1, -1, k, t0, S0, C0)) ;
  I3[1] += 0.5*(S1*L1 - S0*L0) ;
  I3[1] -= kd*(quad_sin_cos_nmr(2, -1, -1, k, t1, S1, C1) -
	       quad_sin_cos_nmr(2, -1, -1, k, t0, S0, C0)) ;

  I3[2] = -0.5*(C1*L1 - C0*L0) ;

  if ( hs ) {
    I5[1] = (quad_sin_cos_nmr(0, 1, -3, k, t1, S1, C1) -
	     quad_sin_cos_nmr(0, 1, -3, k, t0, S0, C0))*kd*kd*kd/3.0/z2 ;
    I5[2] = (quad_sin_cos_nmr(1, 0, -3, k, t1, S1, C1) -
	     quad_sin_cos_nmr(1, 0, -3, k, t0, S0, C0))*kd*kd*kd/3.0/z2 ;
  }

  if ( order == 1 ) return ;

  g_assert_not_reached() ;

  return ;
}

static void subtriangle_polar_quad_ip(gdouble r1, gdouble r2, gdouble theta, 
				      gint order, gboolean hs,
				      gdouble *phi, 
				      gdouble *I1, gdouble *I3, gdouble *I5)

{
  gdouble a, b, b2, t1, t0, S0, S1, C0, C1 ;

  a = (r2*cos(theta) - r1)/r2/sin(theta) ; *phi = atan(a) ;

  b2 = (r1*r1)/(1+a*a) ; b = sqrt(b2) ;

  t1 = theta+(*phi) ; t0 = (*phi) ;
  S0 = sin(t0) ; S1 = sin(t1) ; C0 = cos(t0) ; C1 = cos(t1) ;

  I1[0] = (quad_sin_cos(0, -1, t1, S1, C1) -
	   quad_sin_cos(0, -1, t0, S0, C0))*b ;

  I3[0] = -(quad_sin_cos(0, 1, t1, S1, C1) -
	    quad_sin_cos(0, 1, t0, S0, C0))/b ;
 
  if ( hs ) {
    I5[0] = -(quad_sin_cos(0, 3, t1, S1, C1) -
	      quad_sin_cos(0, 3, t0, S0, C0))/3.0/b/b2 ;
  }

  if ( order == 0 ) return ;

  I1[1] = 0.5*b2*(quad_sin_cos(0, -1, t1, S1, C1) -
		  quad_sin_cos(0, -1, t0, S0, C0)) ;
  I1[2] = 0.5*b2*(quad_sin_cos(1, -2, t1, S1, C1) -
		  quad_sin_cos(1, -2, t0, S0, C0)) ;
  I3[1] = -S1*log(C1/b) + S0*log(C0/b) -
    (quad_sin_cos(2, -1, t1, S1, C1) -
     quad_sin_cos(2, -1, t0, S0, C0)) ;
  I3[2] = C1*log(C1/b) - C0*log(C0/b) +
    (quad_sin_cos(1, 0, t1, S1, C1) -
     quad_sin_cos(1, 0, t0, S0, C0)) ;

  if ( hs ) {
    I5[1] =  -(quad_sin_cos(0, 3, t1, S1, C1) -
	       quad_sin_cos(0, 3, t0, S0, C0))/2.0/b2 ;
    I5[2] =  -(quad_sin_cos(1, 2, t1, S1, C1) -
	       quad_sin_cos(1, 2, t0, S0, C0))/2.0/b2 ;
  }

  if ( order == 1 ) return ;

  g_assert_not_reached() ;

  return ;
}

static void subtriangle_potential_quad(gdouble *p,
				       gdouble *x1, gdouble *x2, 
				       gdouble r1, gdouble r2,
				       gdouble ont,
				       gint order, gboolean hs,
				       gdouble *I)

{
  gdouble J1[16], J3[16], J5[16] = {0}, th, phi, psi, C1, S1, C0, S0 ;
  gdouble *I3 = NULL, *I5 = NULL ;

  if ( r1 == 0.0 || r2 == 0.0 ) return ;
  if ( fabs(th = ((x2[0]-p[0])*(x1[0]-p[0]) + 
		  (x2[1]-p[1])*(x1[1]-p[1]))/r1/r2) >= 1.0 ) return ;
  if ( (th = ont*acos(th)) == 0.0 ) return ;
  g_assert(!isnan(th)) ;

  switch (order) {
  default: g_assert_not_reached() ; break ;
  case 0: I3 = &(I[1]) ; I5 = &(I[2]) ; break ;
  case 1: I3 = &(I[3]) ; I5 = &(I[6]) ; break ;
  }

  if ( p[2] != 0.0 ) 
    subtriangle_polar_quad(r1, r2, th, p[2], order, hs, &phi, J1, J3, J5) ;
  else
    subtriangle_polar_quad_ip(r1, r2, th, order, hs, &phi, J1, J3, J5) ;
  C0 = cos(phi) ; S0 = sin(phi) ;

  I[0] += J1[0] ; I3[0] += J3[0] ; if ( hs ) I5[0] += J5[0] ;
  if ( order == 0 ) return ;
  psi = atan2(x1[1]-p[1], x1[0]-p[0]) ; C1 = cos(psi) ; S1 = sin(psi) ;

  I[1]  += (J1[1]*C0 + J1[2]*S0)*C1 - (J1[2]*C0 - J1[1]*S0)*S1 + p[0]*J1[0] ;
  I[2]  += (J1[1]*C0 + J1[2]*S0)*S1 + (J1[2]*C0 - J1[1]*S0)*C1 + p[1]*J1[0] ;
  I3[1] += (J3[1]*C0 + J3[2]*S0)*C1 - (J3[2]*C0 - J3[1]*S0)*S1 + p[0]*J3[0] ;
  I3[2] += (J3[1]*C0 + J3[2]*S0)*S1 + (J3[2]*C0 - J3[1]*S0)*C1 + p[1]*J3[0] ;
  if ( hs ) {
    I5[1] += (J5[1]*C0 + J5[2]*S0)*C1 - (J5[2]*C0 - J5[1]*S0)*S1 + p[0]*J5[0] ;
    I5[2] += (J5[1]*C0 + J5[2]*S0)*S1 + (J5[2]*C0 - J5[1]*S0)*C1 + p[1]*J5[0] ;
  }

  if ( order == 1 ) return ;

  return ;
}

static gint triangle_potential_quad(gdouble *x1, gdouble *x2, gdouble *x3,
				    gdouble *p, 
				    gdouble o12, gdouble o23, gdouble o31,
				    gint order, gboolean hs,
				    gdouble *I)

/*
  evaluate potential integrals over triangle

  x1, x2, x3: vertices of triangle, all in plane z=0;
  p:          field point for potential, arbitrary position;
  o12, o23, o13: orientation (sign of signed area) of subtriangles (p, x1, x2),
                 (p, x2, x3) and (p, x3, x1);
  order:      order of integral with:
                  0: constant source;
		  1: linear variation;

  hs:         if TRUE hypersingular integrals are included, otherwise not.
  I:          pre-allocated array to hold output

  Returns 0 on success.

  Output arrangement is integrals over the triangle in the following
  order in I[]:

  Order 0: 1/R 1/R^{3} [1/R^{5}]
  Order 1: 1/R x/R y/R 1/R^{3} x/R^{3} y/R^{3} [1/R^{5} x/R^{5} y/R^{5}]

  Terms in square brackets are included if hs == TRUE.
 */


{
  gdouble r1, r2, r3 ;  
  gint ni, nhs, i ;

  g_assert(order >= 0 && order <= 1) ;

  switch ( order ) {
  case 0: ni = 2 ; nhs = 1 ; break ;
  case 1: ni = 6 ; nhs = 3 ; break ;
  default: g_assert_not_reached() ; break ;
  }

  for ( i = 0 ; i < ni ; i ++ ) I[i] = 0.0 ;
  if ( hs ) for ( i = ni ; i < ni+nhs ; i ++ ) I[i] = 0.0 ;

  r1 = sqrt((p[0]-x1[0])*(p[0]-x1[0]) +  (p[1]-x1[1])*(p[1]-x1[1])) ;
  r2 = sqrt((p[0]-x2[0])*(p[0]-x2[0]) +  (p[1]-x2[1])*(p[1]-x2[1])) ;
  r3 = sqrt((p[0]-x3[0])*(p[0]-x3[0]) +  (p[1]-x3[1])*(p[1]-x3[1])) ;

  subtriangle_potential_quad(p, x1, x2, r1, r2, o12, order, hs, I) ;
  subtriangle_potential_quad(p, x2, x3, r2, r3, o23, order, hs, I) ;
  subtriangle_potential_quad(p, x3, x1, r3, r1, o31, order, hs, I) ;

  return 0 ;
}

static gint triangle_quad_shape0(gdouble *x1, gdouble *x2, gdouble *x3,
				 gdouble *p, 
				 gdouble o12, gdouble o23, gdouble o31,
				 gint order, gboolean hs,
				 gdouble *G1, gdouble *G3, gdouble *G5)

/*
  Integration of 1/R and normal derivatives over a triangle, weighted
  with linear shape functions.

  Input parameters as for triangle_potential_quad

  Output: 
  
  G1: \int L_{1,2,3}/R dS
  G3: \int L_{1,2,3}/R^{3} dS
  G5: \int L_{1,2,3}/R^{5} dS (if requested with hs == TRUE, otherwise 
      G5 can be NULL)
 */

{
  gdouble I[16], Ai[9], A[9] ;

  triangle_potential_quad(x1, x2, x3, p, o12, o23, o31, 1, hs, I) ;

  A[0] = 1.0   ; A[1] = 1.0   ; A[2] = 1.0   ;
  A[3] = x1[0] ; A[4] = x2[0] ; A[5] = x3[0] ;
  A[6] = x1[1] ; A[7] = x2[1] ; A[8] = x3[1] ;

  invert3x3(Ai, A) ;

  matvecmul3x3(G1, Ai, I) ;
  matvecmul3x3(G3, Ai, &(I[3])) ;

  G3[0] *= -p[2] ; G3[1] *= -p[2] ; G3[2] *= -p[2] ;

  if ( hs ) matvecmul3x3(G5, Ai, &(I[6])) ;

  return 0 ;
}

gint triangle_quad_shape(gdouble *x1, gdouble *x2, gdouble *x3,
			 gdouble *p,
			 gdouble o12, gdouble o23, gdouble o31,
			 gint order, gboolean hs,
			 gdouble *G1, gdouble *G3, gdouble *G5)
  
{
  gdouble A[9], y[3], y1[3], y2[3], y3[3] ;

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

  triangle_quad_shape0(y1, y2, y3, y, o12, o23, o31, order, hs, G1, G3, G5) ;

  return 0 ;
}
