/* hayami.c
 * 
 * Copyright (C) 2006--2013, 2018 Michael Carley
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gts.h>
#include <glib.h>

#ifdef HAVE_GQR
#include <gqr.h>
#else /*HAVE_GQR*/
#error "GQR is required"
#endif /*HAVE_GQR*/

#include "bem3d.h"
#include "bem3d-private.h"

/**
 * @ingroup quadrature
 * @{
 * 
 */

static gint _bem3d_hayami_triangle_rule(gdouble hj, gdouble aj,
					gdouble d, gdouble dtheta,
					gqr_rule_t *qn,
					gqr_rule_t *qm,
					gdouble *rule)

{
  gdouble rmin, rmax, r, drdR ;
  gdouble Rmin, Rmax, Rbar, dR, R ;
  gdouble tmin, tmax, tbar, dt, t ;
  gdouble theta, E, J ;
  gdouble r0, r1 ;
  gdouble Sdth ;
  gint n, m, M ;

  M = gqr_rule_length(qm) ;

  tmin = 0.5*hj*log((1+sin(0-aj))/(1-sin(0-aj))) ;
  tmax = 0.5*hj*log((1+sin(dtheta-aj))/(1-sin(dtheta-aj))) ;

  gqr_rule_scale(qn, tmin, tmax, &tbar, &dt) ;
  Sdth = sin(dtheta) ;
  r0 = hj/cos(aj) ; r1 = hj/cos(dtheta-aj) ;

  J = 1.0 ;
  for ( n = 0 ; n < gqr_rule_length(qn) ; n ++ ) {
    t = tbar + dt*gqr_rule_abscissa(qn,n) ;
    E = exp(2.0*t/hj) ;
    theta = asin((E-1)/(E+1)) + aj ;
    rmin = 0.0 ; rmax = hj/cos(theta-aj) ;
    Rmin = log(rmin+d) ; Rmax = log(rmax+d) ;
    gqr_rule_scale(qm, Rmin, Rmax, &Rbar, &dR) ;
    for ( m = 0 ; m < gqr_rule_length(qm) ; m ++ ) {
      R = Rbar + dR*gqr_rule_abscissa(qm, m) ;
      r = exp(R)-d ; drdR = r+d ;
      g_assert(r != 0.0) ;
      g_assert((rule[3*(n*M+m)+0] = r*sin(theta)/r1/Sdth) != 0.0) ;
      rule[3*(n*M+m)+1] = r*sin(dtheta-theta)/r0/Sdth ;
      rule[3*(n*M+m)+2] = J*r*drdR/rmax*dt*dR/r0/r1/Sdth*
	gqr_rule_weight(qn,n)*gqr_rule_weight(qm,m) ;
    }
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Quadrature rule for near-singular integrals, using the method of
 * Hayami, Ken, `Variable transformations for nearly singular
 * integrals in the boundary element method', Publ. RIMS, Kyoto Univ.,
 * 41:821--842, 2005 and Hayami, Ken & Matsumoto, Hideki, `A numerical
 * quadrature for nearly singular boundary element integrals',
 * Engineering Analysis with Boundary Elements, 13:143--154, 1994. The
 * integration point \a xs may not lie on the element.
 * 
 * @param xs field point;
 * @param e element to integrate over;
 * @param q quadrature rule to fill;
 * @param gfunc ignored;
 * @param param ignored;
 * @param data gint[2], first element number of quadrature points in radius;
 * second element number of points in angle;
 * @param work workspace.
 * 
 * @return ::BEM3D_SUCCESS on success. 
 */

gint bem3d_quadrature_rule_hayami(GtsPoint *xs, BEM3DElement *e,
				  BEM3DQuadratureRule *q, 
				  BEM3DGreensFunction *gfunc,
				  BEM3DParameters *param,
				  gpointer data,
				  BEM3DWorkspace *work)

{
  GtsPoint *xj, *xjp1 ;
  gdouble d, L[32] ;
  gdouble hj, aj, dtheta ;
  BEM3DShapeFunc shfunc = NULL, shfsub = NULL ;
  gint i, j, k, nc ;
  gint M, N ;
  gdouble xi0, eta0, xi1, eta1, xi2, eta2 ;
  gdouble xn, en, ww ;
  gint rot[32] ;
  GtsPoint *x = NULL, *xst = NULL, *fj = NULL ;
  gqr_rule_t *qn, *qm ;

  g_debug("%s: ", __FUNCTION__) ;

  /* g_error("%s: untested code", __FUNCTION__) ; */

  g_return_val_if_fail(xs != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(xs), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  x = bem3d_workspace_gts_point_get(work) ;
  xst = bem3d_workspace_gts_point_get(work) ;
  fj = bem3d_workspace_gts_point_get(work) ;
  qn = bem3d_workspace_gqr_rule_get(work) ;
  qm = bem3d_workspace_gqr_rule_get(work) ;

  /*number of corners on element*/
  nc = bem3d_element_corner_number(e) ; 
  switch ( nc ) {
  default: g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
		 "%s: element must have three or four sides",
		 __FUNCTION__) ;
    break ;
  case 3: shfsub = bem3d_shfunc_t1 ; break ;
  case 4: shfsub = bem3d_shfunc_q1 ; break ;
  }
  rotation_indices(nc, rot) ;

  /*number of radial and angular quadrature points*/
  N = ((gint *)data)[0] ; M = ((gint *)data)[1] ;
  if ( N <= 0 ) 
    g_error("%s: number of radial quadrature points (%d) must be greater "
	    "than zero", __FUNCTION__, N) ;
  if ( M <= 0 ) 
    g_error("%s: number of angular quadrature points (%d) must be greater "
	    "than zero", __FUNCTION__, M) ;

  bem3d_quadrature_rule_realloc(q, 3*M*N) ;
  bem3d_quadrature_clear(q) ;
  
  qm = gqr_rule_realloc(qm, M) ; qn = gqr_rule_realloc(qn, N) ;
  gqr_rule_select(qn, GQR_GAUSS_LEGENDRE, N, NULL) ;
  gqr_rule_select(qm, GQR_GAUSS_LEGENDRE, M, NULL) ;

  /*location of quadrature point*/  
  xi0 = eta0 = 0.0 ;
  bem3d_element_nearest_point(e, xs, &xi0, &eta0, TRUE) ;

  shfunc = bem3d_element_shape_func(e) ;
  shfunc(xi0, eta0, L, NULL, NULL, NULL) ;
  bem3d_element_position(e, L, x) ;
  d = gts_point_distance(x, xs) ;

  shfsub(xi0, eta0, L, NULL, NULL, NULL) ;
  gts_point_set(xst, 0, 0, 0) ;
  for ( i = 0 ; i < bem3d_element_corner_number(e) ; i ++ ) {
    GTS_POINT(xst)->x += L[i]*GTS_POINT(bem3d_element_corner(e,i))->x ;
    GTS_POINT(xst)->y += L[i]*GTS_POINT(bem3d_element_corner(e,i))->y ;
    GTS_POINT(xst)->z += L[i]*GTS_POINT(bem3d_element_corner(e,i))->z ;
  }

  for ( j = 0 ; j < nc ; j ++ ) {
    k = bem3d_element_corner_index(e, rot[j]) ;
    xj = bem3d_element_vertex(e, k) ;
    xi1 = bem3d_element_vertex_xi(e,k) ;
    eta1 = bem3d_element_vertex_eta(e,k) ;

    k = bem3d_element_corner_index(e, rot[j+1]) ;
    xjp1 = bem3d_element_vertex(e, k) ;
    xi2 = bem3d_element_vertex_xi(e,k) ;
    eta2 = bem3d_element_vertex_eta(e,k) ;

    _bem3d_edge_nearest_point(xj, xjp1, xst, fj) ;
    hj = gts_point_distance(xst, fj) ;
    if ( hj > 1e-9 ) {
      dtheta = _bem3d_point_internal_angle(xj, xst, xjp1) ;
      aj = _bem3d_point_internal_angle(xj, xst, fj) ;
      i = bem3d_quadrature_vertex_number(q) ;
      _bem3d_hayami_triangle_rule(hj, aj, d, dtheta, qn, qm, q->rule) ;
      for ( k = 0 ; k < N*M ; k ++ ) {
	_bem3d_quadrature_rule_remap(xi0, eta0, xi1, eta1,
				     xi2, eta2,
				     q->rule[3*(i+k)+0],
				     q->rule[3*(i+k)+1],
				     q->rule[3*(i+k)+2],		    
				     &xn, &en, &ww) ;
 	_bem3d_quadrature_add_point(q, xn, en, ww) ;
      }
    }
  }

  bem3d_workspace_gqr_rule_put(work, qn) ;
  bem3d_workspace_gqr_rule_put(work, qm) ;
  bem3d_workspace_gts_point_put(work, x) ;
  bem3d_workspace_gts_point_put(work, xst) ;
  bem3d_workspace_gts_point_put(work, fj) ;

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
