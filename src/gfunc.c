/* gfunc.c
 * 
 * Copyright (C) 2006, 2008 Michael Carley
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


/**
 * @defgroup gfunc Green's functions
 * @{
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <glib.h>
#include <gts.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "bem3d.h"
#include "bem3d-private.h"

/** 
 * Green's function for Laplace equation:
 * \f$G=1/4\pi R\f$, \f$R=|\mathbf{x}-\mathbf{y}|\f$
 * 
 * @param x field point
 * @param y source point
 * @param n normal at source point
 * @param p a ::BEM3DParameters struct, ignored in this function;
 * @param G single element GArray containing Green's function
 * @param dGdn single element GArray containing normal derivative of
 * Green's function
 * 
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_laplace(GtsPoint *x, GtsPoint *y,
			       GtsVector n, BEM3DParameters *p,
			       GArray *G, GArray *dGdn)

{
  gdouble R2, R, g, dRdn ;
  GtsVector r ;

/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(y != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(n != NULL, BEM3D_NULL_ARGUMENT) ;
  /* g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ; */
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  g_array_set_size(G, 1) ; g_array_set_size(dGdn, 1) ;
  R2 = gts_point_distance2(x, y) ;
  g_assert(R2 != 0.0) ;
  R = sqrt(R2) ;

  gts_vector_init(r, y, x) ;
  
  g = 0.25*M_1_PI/R ;
  g_array_index(G, gdouble, 0) = g ;

  dRdn = -gts_vector_scalar(r,n)/R ;
  g_array_index(dGdn, gdouble, 0) = -dRdn*g/R ;

  return BEM3D_SUCCESS ;
}

/** 
 * Green's function for gradient of field in Laplace equation:
 * \f$G=-1/4\pi R^{2}\nabla R$, $R=|\mathbf{x}-\mathbf{y}|\f$
 * 
 * @param x field point
 * @param y source point
 * @param n normal at source point
 * @param p a ::BEM3DParameters struct, ignored in this function;
 * @param G four element GArray containing Green's function and its gradient
 * @param dGdn four element GArray containing value and gradient of normal 
 * derivative of Green's function
 * 
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_gradient_laplace(GtsPoint *x, GtsPoint *y,
					GtsVector n, BEM3DParameters *p,
					GArray *G, GArray *dGdn)

{
  gdouble R2, R, g, dRdn ;
  GtsVector r, nablaR ;

/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(y != NULL, BEM3D_NULL_ARGUMENT) ;
/*   g_return_val_if_fail(GTS_IS_POINT(y), BEM3D_ARGUMENT_WRONG_TYPE) ; */
  g_return_val_if_fail(n != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  g_array_set_size(G, 4) ; g_array_set_size(dGdn, 4) ;
  R2 = gts_point_distance2(x, y) ;
  if ( R2 == 0.0 ) 
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	  "%s: coincident points x=(%g,%g,%g); y=(%g,%g,%g)",
	  __FUNCTION__,
	  GTS_POINT(x)->x, GTS_POINT(x)->y, GTS_POINT(x)->z,
	  GTS_POINT(y)->x, GTS_POINT(y)->y, GTS_POINT(y)->z) ;
  R = sqrt(R2) ;

  gts_vector_init(r, y, x) ;
  nablaR[0] = r[0]/R ; nablaR[1] = r[1]/R ;
  nablaR[2] = r[2]/R ;

  g = 0.25*M_1_PI/R ;
  g_array_index(G, gdouble, 0) = g ;

  dRdn = -gts_vector_scalar(r,n)/R ;
  g_array_index(dGdn, gdouble, 0) = -dRdn*g/R ;

  g = -0.25*M_1_PI/R2 ;
  g_array_index(G,gdouble,1) = g*nablaR[0] ;
  g_array_index(G,gdouble,2) = g*nablaR[1] ;
  g_array_index(G,gdouble,3) = g*nablaR[2] ;

  g_array_index(dGdn,gdouble,1) = 3*dRdn*nablaR[0] + n[0] ;
  g_array_index(dGdn,gdouble,2) = 3*dRdn*nablaR[1] + n[1] ;
  g_array_index(dGdn,gdouble,3) = 3*dRdn*nablaR[2] + n[2] ;

  g /= R ;
  g_array_index(dGdn,gdouble,1) *= -g ;
  g_array_index(dGdn,gdouble,2) *= -g ;
  g_array_index(dGdn,gdouble,3) *= -g ;
  
  return BEM3D_SUCCESS ;
}

/** 
 * Green's function for gradient of Helmholtz equation:
 * \f$G=\exp(\mathrm{J} kR)/4\pi R\f$, \f$R=|\mathbf{x}-\mathbf{y}|\f$
 *
 * @param x field point
 * @param y source point
 * @param n normal at source point
 * @param p a ::BEM3DParameters struct;
 * @param G eight element GArray containing real and imaginary parts of
 * value and gradient of Green's function
 * @param dGdn eight element GArray containing real and imaginary parts
 * of value and gradient of normal derivative of Green's function
 *
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_gradient_helmholtz(GtsPoint *x, GtsPoint *y,
					  GtsVector n, BEM3DParameters *p,
					  GArray *G, GArray *dGdn)

{
  gdouble R2, R, k, dRdn ;
  GtsVector r, nablaR ;
  gsl_complex E, *g, *dgdn ;
  gint i ;
  
  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(y != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(n != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( G->len != 8 ) g_array_set_size(G, 8) ;
  if ( dGdn->len != 8 ) g_array_set_size(dGdn, 8) ;

  k = bem3d_parameters_wavenumber(p) ;

  g_assert((R2 = gts_point_distance2(x, y)) != 0.0) ;
  R = sqrt(R2) ;
  gts_vector_init(r, y, x) ;
  nablaR[0] = r[0]/R ; nablaR[1] = r[1]/R ; nablaR[2] = r[2]/R ;
  dRdn = -gts_vector_scalar(r,n)/R ;

  GSL_SET_COMPLEX(&E, 0.25*M_1_PI*cos(k*R), 0.25*M_1_PI*sin(k*R)) ;
  g = (gsl_complex *)(&(g_array_index(G,gdouble,0))) ; 
  dgdn = (gsl_complex *)(&(g_array_index(dGdn,gdouble,0))) ;

  GSL_SET_COMPLEX(g, 1.0/R, 0) ;
  GSL_SET_COMPLEX(dgdn, -1.0/R2*dRdn, k/R*dRdn) ;
  *g = gsl_complex_mul(*g, E) ;
  *dgdn = gsl_complex_mul(*dgdn, E) ;

  E = gsl_complex_div_real(E, R2) ;
  for ( i = 0 ; i < 3 ; i ++ ) {
    GSL_SET_COMPLEX(&(g[i+1]), (-1.0*nablaR[i]), (k*R*nablaR[i])) ;
    g[i+1] = gsl_complex_mul(g[i+1], E) ;
  }

  E = gsl_complex_div_real(E, R) ;
  
  for ( i = 0 ; i < 3 ; i ++ ) {
    GSL_SET_COMPLEX(&(dgdn[i+1]),
		    ((3.0-k*k*R2)*dRdn*nablaR[i]+n[i]),
		    ((-3.0*k*R)*dRdn*nablaR[i]-k*R*n[i])) ;
    dgdn[i+1] = gsl_complex_mul(dgdn[i+1], E) ;
  }
  
  return BEM3D_SUCCESS ;
}

/** 
 * Green's function for Helmholtz equation:
 * \f$G=\exp(\mathrm{J} kR)/4\pi R\f$, \f$R=|\mathbf{x}-\mathbf{y}|\f$
 *
 * @param x field point
 * @param y source point
 * @param n normal at source point
 * @param p a ::BEM3DParameters struct;
 * @param G two element GArray containing real and imaginary parts of
 * Green's function
 * @param dGdn two element GArray containing real and imaginary parts
 * of normal derivative of Green's function
 *
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_helmholtz(GtsPoint *x, GtsPoint *y,
				 GtsVector n, BEM3DParameters *p,
				 GArray *G, GArray *dGdn)

{
  gdouble R2, R, rny, k ;
  GtsVector r ;
  gsl_complex E, *g, *dgdn ;

/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(y != NULL, BEM3D_NULL_ARGUMENT) ;
/*   g_return_val_if_fail(GTS_IS_POINT(y), BEM3D_ARGUMENT_WRONG_TYPE) ; */
  g_return_val_if_fail(n != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( G->len != 2 ) g_array_set_size(G, 2) ; 
  if ( dGdn->len != 2) g_array_set_size(dGdn, 2) ;

  k = bem3d_parameters_wavenumber(p) ;

  g_assert((R2 = gts_point_distance2(x, y)) != 0.0) ;
  R = sqrt(R2) ;

  GSL_SET_COMPLEX(&E, 0.25*M_1_PI*cos(k*R), 0.25*M_1_PI*sin(k*R)) ;
  g = (gsl_complex *)(&(g_array_index(G,gdouble,0))) ; 
  dgdn = (gsl_complex *)(&(g_array_index(dGdn,gdouble,0))) ;

  gts_vector_init(r, y, x) ; rny = gts_vector_scalar(r,n) ;

  GSL_SET_COMPLEX(g, 1.0/R, 0) ;
  GSL_SET_COMPLEX(dgdn, 1.0/R/R2*rny, -k/R2*rny) ;
  
  *g = gsl_complex_mul(*g, E) ;
  *dgdn = gsl_complex_mul(*dgdn, E) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Green's function for Helmholtz equation, hypersingular formulation:
 * \f$G=\exp(\mathrm{J} kR)/4\pi R\f$, \f$R=|\mathbf{x}-\mathbf{y}|\f$
 *
 * @param x field point
 * @param y source point
 * @param ny normal at source point
 * @param p a ::BEM3DParameters struct;
 * @param G two element GArray containing real and imaginary parts of
 * Green's function
 * @param dGdn two element GArray containing real and imaginary parts
 * of normal derivative of Green's function
 *
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_helmholtz_hs(GtsPoint *x, GtsPoint *y,
				    GtsVector ny, BEM3DParameters *p,
				    GArray *G, GArray *dGdn)

{
  gdouble R2, R, rnx, rny, k, alpha, *nx ;
  GtsVector r ;
  gsl_complex E, *g, *dgdn, dgdnx, dgdny, d2gdnxdny, lambda, t1, t2 ;

/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(y != NULL, BEM3D_NULL_ARGUMENT) ;
/*   g_return_val_if_fail(GTS_IS_POINT(y), BEM3D_ARGUMENT_WRONG_TYPE) ; */
  g_return_val_if_fail(ny != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( G->len != 2 ) g_array_set_size(G, 2) ;
  if ( dGdn->len != 2) g_array_set_size(dGdn, 2) ;

  k = bem3d_parameters_wavenumber(p) ;
  nx = bem3d_parameters_normal(p) ;
  alpha = bem3d_parameters_conditioning(p) ;
  lambda = *((gsl_complex *)(&(bem3d_parameters_lambda_real(p)))) ;

  g_assert((R2 = gts_point_distance2(x, y)) != 0.0) ;
  R = sqrt(R2) ;

  GSL_SET_COMPLEX(&E, 0.25*M_1_PI*cos(k*R), 0.25*M_1_PI*sin(k*R)) ;
  g = (gsl_complex *)(&(g_array_index(G,gdouble,0))) ;
  dgdn = (gsl_complex *)(&(g_array_index(dGdn,gdouble,0))) ;

  gts_vector_init(r, y, x) ; 
  rny = gts_vector_scalar(r,ny) ; rnx = gts_vector_scalar(r,nx) ;

  /*(1-\alpha)G*/
  GSL_SET_COMPLEX(g, (1.0-alpha)/R, 0) ;

  GSL_SET_COMPLEX(&t1, -1.0, k*R) ;
  /*(1-\alpha)dG/dny*/
  /* GSL_SET_COMPLEX(&dgdny,  (1.0-alpha)/R/R2*rny, -k*(1.0-alpha)/R2*rny) ; */
  dgdny = gsl_complex_mul_real(t1, -rny/R/R2) ;
  dgdny = gsl_complex_mul_real(dgdny, (1.0-alpha)) ;
  /*\alpha\lambda dG/dnx*/
  /* GSL_SET_COMPLEX(&dgdnx, -alpha/R/R2*rnx,  k*alpha/R2*rnx) ; */
  dgdnx = gsl_complex_mul_real(t1, rnx/R/R2) ;
  dgdnx = gsl_complex_mul_real(dgdnx, alpha) ;
  dgdnx = gsl_complex_mul(dgdnx, lambda) ;

  /*-\alpha\lambda d^2G/dnx dny*/
  GSL_SET_COMPLEX(&t2, -3.0, k*R) ;
  t2 = gsl_complex_mul(t1, t2) ;
  t2 = gsl_complex_add_imag(t2, k*R) ;
  t2 = gsl_complex_mul_real(t2, rnx*rny/R2/R2/R) ;
  t1 = gsl_complex_mul_real(t1, gts_vector_scalar(nx,ny)/R/R2) ;
  d2gdnxdny = gsl_complex_add(t2, t1) ;
  d2gdnxdny = gsl_complex_mul(d2gdnxdny, lambda) ;
  d2gdnxdny = gsl_complex_mul_real(d2gdnxdny, alpha) ;

  *dgdn = dgdny ;
  *dgdn = gsl_complex_sub(*dgdn, d2gdnxdny) ;

  *g = gsl_complex_add(*g, dgdnx) ;

  *g = gsl_complex_mul(*g, E) ;
  *dgdn = gsl_complex_mul(*dgdn, E) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Green's function for convected Helmholtz equation:
 * \f$G=\exp(\mathrm{j} k\sigma)/4\pi S\f$
 *
 * @param x field point
 * @param y source point
 * @param n normal at source point
 * @param p a ::BEM3DParameters struct;
 * @param G two element GArray containing real and imaginary parts of
 * Green's function
 * @param dGdn two element GArray containing real and imaginary parts
 * of normal derivative of Green's function
 *
 * @return ::BEM3D_SUCCESS on success 
 */

/*written by Lydia Meli and Vera Castiglione*/

gint bem3d_greens_func_convected_helmholtz(GtsPoint *x, GtsPoint *y,
					   GtsVector n, BEM3DParameters *p,
					   GArray *G, GArray *dGdn)

{
  gdouble  C, S, Ai, Ar ;
  GtsVector r ;
  gdouble beta2, S2, S1, sigma , dS1dn, k, M ;

/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(y != NULL, BEM3D_NULL_ARGUMENT) ;
/*   g_return_val_if_fail(GTS_IS_POINT(y), BEM3D_ARGUMENT_WRONG_TYPE) ; */
  g_return_val_if_fail(n != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  g_array_set_size(G,2) ; g_array_set_size(dGdn, 2) ;
  g_assert_not_reached() ; /*unchecked code*/
  k = bem3d_parameters_wavenumber(p) ;
  M = bem3d_parameters_mach_number(p) ;

  gts_vector_init(r, y, x) ;

  beta2 = 1.0-M*M ;
  g_assert(beta2 > 0.0) ;
  S2 =(r[0]*r[0]) + beta2*((r[1]*r[1]) + (r[2]*r[2])) ;
  g_assert(S2 != 0.0) ;
  S1 = sqrt(S2) ;
  sigma = (S1+M*r[0])/beta2 ;
                                      
  C = cos(k*sigma) ; S = sin(k*sigma) ;
  
  g_array_index(G, gdouble, 0) = C*0.25*M_1_PI/S1 ;
  g_array_index(G, gdouble, 1) = S*0.25*M_1_PI/S1 ;

  dS1dn = -(r[0]*n[0]+ beta2*(r[1]*n[1] + r[2]*n[2]))/S1 ;

  Ai = k/beta2*S1 ; Ar = -1 ;

  g_array_index(dGdn, gdouble, 0) = 
    (C*Ar - S*Ai)*dS1dn - M*n[0]*S*Ai ;
  g_array_index(dGdn, gdouble, 0) *= -0.25*M_1_PI/S2 ;
  g_array_index(dGdn, gdouble, 1) = 
    (S*Ar + C*Ai)*dS1dn + M*n[0]*C*Ai ;
  g_array_index(dGdn, gdouble, 1) *= -0.25*M_1_PI/S2 ;
  
  return BEM3D_SUCCESS ;
}

BEM3DParameters *bem3d_parameters_new(void)

{
  BEM3DParameters *p ;

  p = (BEM3DParameters *)g_malloc(sizeof(BEM3DParameters)) ;

  memset(p->f,0,BEM3D_PARAMETERS_SIZE*sizeof(gdouble)) ;
  memset(p->n,0,BEM3D_PARAMETERS_SIZE*sizeof(gint)) ;

  return p ;
}

/**
 * @}
 * 
 */

