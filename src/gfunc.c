/* gfunc.c
 * 
 * Copyright (C) 2006, 2008, 2018 Michael Carley
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
 * @param G single element array containing Green's function
 * @param dGdn single element array containing normal derivative of
 * Green's function
 * 
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_laplace(GtsPoint *x, GtsPoint *y,
			       GtsVector n, BEM3DParameters *p,
			       gdouble *G, gdouble *dGdn)

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

  R2 = gts_point_distance2(x, y) ;
  g_assert(R2 != 0.0) ;
  R = sqrt(R2) ;

  gts_vector_init(r, y, x) ;
  
  G[0] = g = 0.25*M_1_PI/R ;

  dRdn = -gts_vector_scalar(r,n)/R ;
  dGdn[0] = -dRdn*g/R ;

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
 * @param G four element array containing Green's function and its gradient
 * @param dGdn four element array containing value and gradient of normal 
 * derivative of Green's function
 * 
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_gradient_laplace(GtsPoint *x, GtsPoint *y,
					GtsVector n, BEM3DParameters *p,
					gdouble *G, gdouble *dGdn)

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

  G[0] = g = 0.25*M_1_PI/R ;

  dRdn = -gts_vector_scalar(r,n)/R ;
  dGdn[0] = -dRdn*g/R ;

  g = -0.25*M_1_PI/R2 ;
  G[1] = g*nablaR[0] ;
  G[2] = g*nablaR[1] ;
  G[3] = g*nablaR[2] ;

  dGdn[1] = 3*dRdn*nablaR[0] + n[0] ;
  dGdn[2] = 3*dRdn*nablaR[1] + n[1] ;
  dGdn[3] = 3*dRdn*nablaR[2] + n[2] ;

  g /= R ;
  dGdn[1] *= -g ;
  dGdn[2] *= -g ;
  dGdn[3] *= -g ;
  
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
 * @param G eight element array containing real and imaginary parts of
 * value and gradient of Green's function
 * @param dGdn eight element array containing real and imaginary parts
 * of value and gradient of normal derivative of Green's function
 *
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_gradient_helmholtz(GtsPoint *x, GtsPoint *y,
					  GtsVector n, BEM3DParameters *p,
					  gdouble *G, gdouble *dGdn)

{
  /* g_assert_not_reached() ; /\*untested code*\/ */

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

  /* if ( G->len != 8 ) g_array_set_size(G, 8) ; */
  /* if ( dGdn->len != 8 ) g_array_set_size(dGdn, 8) ; */

  k = bem3d_parameters_wavenumber(p) ;

  g_assert((R2 = gts_point_distance2(x, y)) != 0.0) ;
  R = sqrt(R2) ;
  gts_vector_init(r, y, x) ;
  nablaR[0] = r[0]/R ; nablaR[1] = r[1]/R ; nablaR[2] = r[2]/R ;
  dRdn = -gts_vector_scalar(r,n)/R ;

  GSL_SET_COMPLEX(&E, 0.25*M_1_PI*cos(k*R), 0.25*M_1_PI*sin(k*R)) ;
  /* g = (gsl_complex *)(&(g_array_index(G,gdouble,0))) ;  */
  /* dgdn = (gsl_complex *)(&(g_array_index(dGdn,gdouble,0))) ; */

  g = (gsl_complex *)(&(G[0])) ;
  dgdn = (gsl_complex *)(&(dGdn[0])) ;

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

#if 0
#endif
  
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
 * @param G two element array containing real and imaginary parts of
 * Green's function
 * @param dGdn two element array containing real and imaginary parts
 * of normal derivative of Green's function
 *
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_helmholtz(GtsPoint *x, GtsPoint *y,
				 GtsVector n, BEM3DParameters *p,
				 gdouble *G, gdouble *dGdn)

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

  k = bem3d_parameters_wavenumber(p) ;

  g_assert((R2 = gts_point_distance2(x, y)) != 0.0) ;
  R = sqrt(R2) ;

  GSL_SET_COMPLEX(&E, 0.25*M_1_PI*cos(k*R), 0.25*M_1_PI*sin(k*R)) ;
  g = (gsl_complex *)(&(G[0])) ;
  dgdn = (gsl_complex *)(&(dGdn[0])) ;

  gts_vector_init(r, y, x) ; rny = gts_vector_scalar(r,n) ;

  GSL_SET_COMPLEX(g, 1.0/R, 0) ;
  GSL_SET_COMPLEX(dgdn, 1.0/R/R2*rny, -k/R2*rny) ;
  
  *g = gsl_complex_mul(*g, E) ;
  *dgdn = gsl_complex_mul(*dgdn, E) ;

  return BEM3D_SUCCESS ;
}

#if 0
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
#endif

/** 
 * Green's function for Helmholtz equation, Burton and Miller
 * hypersingular formulation: \f$G=\exp(\mathrm{J} kR)/4\pi R\f$,
 * \f$R=|\mathbf{x}-\mathbf{y}|\f$ . The field point surface normal
 * should be set in the parameters input \a p before calling the
 * Green's function.
 *
 * @param x field point
 * @param y source point
 * @param ny normal at source point
 * @param p a ::BEM3DParameters struct;
 * @param G two element array containing real and imaginary parts of
 * Green's function
 * @param dGdn two element array containing real and imaginary parts
 * of normal derivative of Green's function
 *
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_helmholtz_hs(GtsPoint *x, GtsPoint *y,
				    GtsVector ny, BEM3DParameters *p,
				    gdouble *G, gdouble *dGdn)

{
  gdouble R2, R, dRx, dRy, k, *nx, nxny ;
  GtsVector r ;
  gsl_complex E, *g, *dgdn, g0, dgy, dgx, dgxy, lambda, al ;
  
/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(y != NULL, BEM3D_NULL_ARGUMENT) ;
/*   g_return_val_if_fail(GTS_IS_POINT(y), BEM3D_ARGUMENT_WRONG_TYPE) ; */
  g_return_val_if_fail(ny != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  k = bem3d_parameters_wavenumber(p) ;
  nx = bem3d_parameters_normal(p) ;
  GSL_SET_COMPLEX(&lambda,
		  bem3d_parameters_lambda_real(p),
		  bem3d_parameters_lambda_imag(p)) ;
  GSL_SET_COMPLEX(&al,
		  bem3d_parameters_coupling_real(p),
		  bem3d_parameters_coupling_imag(p)) ;

  g_assert((R2 = gts_point_distance2(x, y)) != 0.0) ;
  R = sqrt(R2) ;

  GSL_SET_COMPLEX(&E, 0.25*M_1_PI*cos(k*R), 0.25*M_1_PI*sin(k*R)) ;
  g = (gsl_complex *)(&(G[0])) ;
  dgdn = (gsl_complex *)(&(dGdn[0])) ;

  nxny = gts_vector_scalar(nx,ny) ;
  gts_vector_init(r, y, x) ; 
  dRx =  gts_vector_scalar(r, nx)/R ;
  dRy = -gts_vector_scalar(r, ny)/R ;
			  
  /*build up different bits of the formulation*/
  /*basic Green's function*/
  g0 = gsl_complex_div_real(E, R) ;

  /*normal derivatives*/
  GSL_SET_COMPLEX(&dgy, -1.0, k*R) ;
  dgy = gsl_complex_mul(dgy, E) ;

  /*normal derivative w.r.t. field normal*/
  dgx = gsl_complex_mul_real(dgy, dRx/R2) ;

  /*normal derivative w.r.t. source normal*/
  dgy = gsl_complex_mul_real(dgy, dRy/R2) ;

  /*double derivative*/
  GSL_SET_COMPLEX(&dgxy,
		  ((3.0-k*k*R*R)*dRx*dRy + nxny)/R/R2,
		  (-3.0*k*R*dRx*dRy - k*R*nxny)/R/R2) ;
  dgxy = gsl_complex_mul(dgxy, E) ;

  /* fprintf(stderr, "x=[%lg %lg %lg] ;\n", */
  /* 	  GTS_POINT(x)->x, GTS_POINT(x)->y, GTS_POINT(x)->z) ; */
  /* fprintf(stderr, "y=[%lg %lg %lg] ;\n", */
  /* 	  GTS_POINT(y)->x, GTS_POINT(y)->y, GTS_POINT(y)->z) ; */
  
  /* fprintf(stderr, "nx=[%lg %lg %lg] ;\n", nx[0], nx[1], nx[2]) ; */
  /* fprintf(stderr, "ny=[%lg %lg %lg] ;\n", ny[0], ny[1], ny[2]) ; */

  /* fprintf(stderr, "G = %lg + j*%lg ; \n", GSL_REAL(g0), GSL_IMAG(g0)) ; */
  /* fprintf(stderr, "Gx = %lg + j*%lg ; \n", GSL_REAL(dgx), GSL_IMAG(dgx)) ; */
  /* fprintf(stderr, "Gy = %lg + j*%lg ; \n", GSL_REAL(dgy), GSL_IMAG(dgy)) ; */
  /* fprintf(stderr, "Gxy = %lg + j*%lg ; \n", GSL_REAL(dgxy), GSL_IMAG(dgxy)) ; */

  /* exit(0) ; */
  
  *g = gsl_complex_mul(lambda, dgx) ;
  g0 = gsl_complex_mul(g0, al) ;
  *g = gsl_complex_add(*g, g0) ;

  *dgdn = gsl_complex_mul(lambda, dgxy) ;
  dgy = gsl_complex_mul(dgy, al) ;
  *dgdn = gsl_complex_add(*dgdn, dgy) ;
  
  return BEM3D_SUCCESS ;
}

/** 
 * Green's function for Helmholtz equation, Burton and Miller
 * hypersingular formulation: \f$G=\exp(\mathrm{J} kR)/4\pi R\f$,
 * \f$R=|\mathbf{x}-\mathbf{y}|\f$ for the method given by Chen and
 * Harris, Applied Numerical Mathematics, 36:475-489, 2001. The field
 * point surface normal should be set in the parameters input \a p
 * before calling the Green's function.
 *
 * @param x field point
 * @param y source point
 * @param ny normal at source point
 * @param p a ::BEM3DParameters struct;
 * @param G two element array containing real and imaginary parts of
 * Green's function
 * @param dGdn six element array containing real and imaginary parts
 * of normal derivative of Green's function and self-terms required by
 * formulation
 *
 * @return ::BEM3D_SUCCESS on success 
 */

gint bem3d_greens_func_helmholtz_ch(GtsPoint *x, GtsPoint *y,
				    GtsVector ny, BEM3DParameters *p,
				    gdouble *G, gdouble *dGdn)

{
  gdouble R2, R, dRx, dRy, k, *nx, nxny ;
  GtsVector r ;
  gsl_complex E, *g, *dgdn, g0, dgy, dgx, dgxy, gkk, lambda, al ;
  
/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(y != NULL, BEM3D_NULL_ARGUMENT) ;
/*   g_return_val_if_fail(GTS_IS_POINT(y), BEM3D_ARGUMENT_WRONG_TYPE) ; */
  g_return_val_if_fail(ny != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  k = bem3d_parameters_wavenumber(p) ;
  nx = bem3d_parameters_normal(p) ;
  GSL_SET_COMPLEX(&lambda,
		  bem3d_parameters_lambda_real(p),
		  bem3d_parameters_lambda_imag(p)) ;
  GSL_SET_COMPLEX(&al,
		  bem3d_parameters_coupling_real(p),
		  bem3d_parameters_coupling_imag(p)) ;

  g_assert(GSL_REAL(al) == 1.0) ;
  g_assert(GSL_IMAG(al) == 0.0) ;
  
  g_assert((R2 = gts_point_distance2(x, y)) != 0.0) ;
  R = sqrt(R2) ;

  G[0] = G[1] = G[2] = G[3] = G[4] = G[5] = 0.0 ;
  dGdn[0] = dGdn[1] = dGdn[2] = dGdn[3] = dGdn[4] = dGdn[5] = 0.0 ;
  
  GSL_SET_COMPLEX(&E, 0.25*M_1_PI*cos(k*R), 0.25*M_1_PI*sin(k*R)) ;
  g = (gsl_complex *)(&(G[0])) ;
  dgdn = (gsl_complex *)(&(dGdn[0])) ;

  nxny = gts_vector_scalar(nx,ny) ;
  gts_vector_init(r, y, x) ; 
  dRx =  gts_vector_scalar(r, nx)/R ;
  dRy = -gts_vector_scalar(r, ny)/R ;
			  
  /*build up different bits of the formulation*/
  /*basic Green's function*/
  g0 = gsl_complex_div_real(E, R) ;

  /*normal derivatives*/
  GSL_SET_COMPLEX(&dgy, -1.0, k*R) ;
  dgy = gsl_complex_mul(dgy, E) ;

  /*normal derivative w.r.t. field normal*/
  dgx = gsl_complex_mul_real(dgy, dRx/R2) ;

  /*normal derivative w.r.t. source normal*/
  dgy = gsl_complex_mul_real(dgy, dRy/R2) ;

  /*double derivative*/
  GSL_SET_COMPLEX(&dgxy,
		  ((3.0-k*k*R*R)*dRx*dRy + nxny)/R/R2,
		  (-3.0*k*R*dRx*dRy - k*R*nxny)/R/R2) ;
  dgxy = gsl_complex_mul(dgxy, E) ;

  gkk = gsl_complex_mul(g0, lambda) ;
  /*in the Chen and Harris formulation G is G+\lambda dG/dn*/
  *g = gsl_complex_mul(lambda, dgx) ;
  g0 = gsl_complex_mul(g0, al) ;
  *g = gsl_complex_add(*g, g0) ;

  /* G[0] = GSL_REAL(g0) + */
  /*   GSL_REAL(lambda)*GSL_REAL(dgx) - */
  /*   GSL_IMAG(lambda)*GSL_IMAG(dgx) ; */
  /* G[1] = GSL_IMAG(g0) + */
  /*   GSL_REAL(lambda)*GSL_IMAG(dgx) + */
  /*   GSL_IMAG(lambda)*GSL_REAL(dgx) ; */
  
  G[4] = G[2] = G[0] ;
  G[5] = G[3] = G[1] ;
  
  /*dG is dG/dn1, with two extra elements for the self-terms in the
    formulation*/
  *dgdn = gsl_complex_mul(lambda, dgxy) ;
  dGdn[4] = GSL_REAL(*dgdn) ;
  dGdn[5] = GSL_IMAG(*dgdn) ;

  /* dgy = gsl_complex_mul(dgy, al) ; */
  /* *dgdn = gsl_complex_add(*dgdn, dgy) ; */
  *dgdn = gsl_complex_mul(dgy, al) ;

  dGdn[2] = GSL_REAL(gkk)*nxny*k*k ;
  dGdn[3] = GSL_IMAG(gkk)*nxny*k*k ;

  
  return BEM3D_SUCCESS ;
}


BEM3DParameters *bem3d_parameters_new(void)

{
  BEM3DParameters *p ;

  p = (BEM3DParameters *)g_malloc0(sizeof(BEM3DParameters)) ;

  memset(p->f,0,BEM3D_PARAMETERS_REAL_SIZE*sizeof(gdouble)) ;
  memset(p->n,0,BEM3D_PARAMETERS_INT_SIZE*sizeof(gint)) ;
  memset(p->p,0,BEM3D_PARAMETERS_POINTER_SIZE*sizeof(gpointer)) ;

  return p ;
}

/**
 * @}
 * 
 */

