/* quadrature.c
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

/**
 * @defgroup quadrature Quadrature rules
 * @{
 * 
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

#include "polar.h"
#include "htriquad.h"

#include "trace.h"

gint triangle_axes(gdouble *x1, gdouble *x2, gdouble *x3,
		   gdouble *s, gdouble *t, gdouble *n) ;
gint triangle_angles(gdouble a, gdouble b, gdouble c,
		     gdouble *A, gdouble *B, gdouble *C) ;

gint triangle_gradient_quad_shape(gdouble *x1, gdouble *x2, gdouble *x3,
				  gdouble *p,
				  gdouble o12, gdouble o23, gdouble o31,
				  gdouble *G, gdouble *dG) ;
gint qtri_quad_rule(gdouble x1[], gdouble x2[], gdouble x3[],
		    gdouble x4[], gdouble x5[], gdouble x6[],
		    gint ntmin, gint ntmax, gdouble dt,
		    gint nrmin, gint nrmax, gdouble dr,
		    gqr_rule_t *gt, gqr_rule_t *gr,
		    gboolean hypersingular,
		    gdouble *xi, gint xistr, 
		    gdouble *eta, gint etastr, 
		    gdouble *wt, gint wtstr,
		    gint *ngp) ;
gint ctri_quad_rule(gdouble x1[], gdouble x2[], gdouble x3[],
		    gdouble x4[], gdouble x5[], gdouble x6[],
		    gdouble x7[], gdouble x8[], gdouble x9[],
		    gdouble x10[],
		    gint ntmin, gint ntmax, gdouble dt,
		    gint nrmin, gint nrmax, gdouble dr,
		    gqr_rule_t *gt, gqr_rule_t *gr,
		    gdouble *xi, gint xistr, 
		    gdouble *eta, gint etastr, 
		    gdouble *wt, gint wtstr,
		    gint *ngp) ;

#define GAUSS_ALPHA_1 0.0597158717
#define GAUSS_BETA_1  0.4701420641
#define GAUSS_ALPHA_2 0.7974269853
#define GAUSS_BETA_2  0.1012865073

const gdouble GAUSS_TRIANGLE_1[] = {
  1/3.0, 1/3.0, 1.0/2.0
} ;
const gdouble GAUSS_TRIANGLE_3[] = {
  0.5, 0.5, 1/3.0/2.0, 
  0.0, 0.5, 1/3.0/2.0, 
  0.5, 0.0, 1/3.0/2.0
} ;

const gdouble GAUSS_TRIANGLE_4[] = {
  1/3.0, 1/3.0, -27.0/48.0/2.0,
  0.6, 0.2, 25.0/48.0/2.0,
  0.2, 0.6, 25.0/48.0/2.0,
  0.2, 0.2, 25.0/48.0/2.0
} ;

const gdouble GAUSS_TRIANGLE_7[] = {
  1/3.0, 1/3.0, 0.225/2.0,
  GAUSS_ALPHA_1, GAUSS_BETA_1, 0.1323941527/2.0,
  GAUSS_BETA_1, GAUSS_ALPHA_1, 0.1323941527/2.0,
  GAUSS_BETA_1, GAUSS_BETA_1,  0.1323941527/2.0,
  GAUSS_ALPHA_2, GAUSS_BETA_2, 0.1259391805/2.0,
  GAUSS_BETA_2, GAUSS_ALPHA_2, 0.1259391805/2.0,
  GAUSS_BETA_2, GAUSS_BETA_2,  0.1259391805/2.0
} ;

extern gdouble WANDZURA_7[], WANDZURA_25[], WANDZURA_54[], 
  WANDZURA_85[], WANDZURA_126[], WANDZURA_175[] ;

static gboolean nodes_coincide(GtsPoint *x, GtsPoint *y)

{
  return (x->x == y->x && x->y == y->y && x->z == y->z) ;
}

static gint _quadrature_add_point(BEM3DQuadratureRule *q,
				  gdouble xi, gdouble eta,
				  gdouble w)

{
  gint n ;

  n = bem3d_quadrature_vertex_number(q) ;
  bem3d_quadrature_xi(q,n) = xi ;
  bem3d_quadrature_eta(q,n) = eta ;
  bem3d_quadrature_weight(q,n) = w ;
  q->n ++ ;

  return BEM3D_SUCCESS ;
}

static gint _quadrature_rule_remap(gdouble xi0, gdouble eta0,
				   gdouble xi1, gdouble eta1,
				   gdouble xi2, gdouble eta2,
				   gdouble s, gdouble t, gdouble w,
				   gdouble *sn, gdouble *tn,
				   gdouble *wn)
{
  gdouble L0, L1, L2 ;

  L0 = 1-s-t ; L1 = s ; L2 = t ;
  *sn = L0*xi0 + L1*xi1 + L2*xi2 ;
  *tn = L0*eta0 + L1*eta1 + L2*eta2 ;

  *wn = w*((xi1-xi0)*(eta2-eta0)-((xi2-xi0)*(eta1-eta0))) ;

  return BEM3D_SUCCESS ;
}
 
static gint edge_nearest_point(GtsPoint *x1, GtsPoint *x2,
			       GtsPoint *x, GtsPoint *c)

{
  gdouble t ;
  GtsVector n, s ;

  gts_vector_init(n, x1, x2) ;
  gts_vector_init(s, x1, x) ;

  t = gts_vector_scalar(s,n)/gts_vector_scalar(n,n) ;

  gts_point_set(c,
		(1-t)*GTS_POINT(x1)->x+t*GTS_POINT(x2)->x,
		(1-t)*GTS_POINT(x1)->y+t*GTS_POINT(x2)->y,
		(1-t)*GTS_POINT(x1)->z+t*GTS_POINT(x2)->z) ;

  return BEM3D_SUCCESS ;
}

static gdouble point_internal_angle(GtsPoint *x1, GtsPoint *x2,
				    GtsPoint *x3)

{
  GtsVector r1, r2 ;
  gdouble c, n1, n2 ;

  gts_vector_init(r1, x1, x2) ;
  gts_vector_init(r2, x3, x2) ;

  c = gts_vector_scalar(r1, r2) ;
  g_assert( (n1 = gts_vector_norm(r1)) != 0.0) ;
  g_assert( (n2 = gts_vector_norm(r2)) != 0.0) ;

  return (acos(c/n1/n2)) ;
}

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

static gint _bem3d_rule_fill_wx(gint n, 
				gdouble xi0, gdouble eta0,
				gdouble xi1, gdouble eta1,
				gdouble xi2, gdouble eta2,
				gdouble *rule)

     /*
       Fill arrays (s,t) with Wandzura-Xiao rule of length n remapped
       onto the triangle (xi0,eta0), (xi1,eta1), (xi2,eta2). Do not
       call this function directly.
      */

{
  gint i ;

  switch (n) {
  case 7:   memcpy(rule, WANDZURA_7,   3*n*sizeof(gdouble)) ; break ;
  case 25:  memcpy(rule, WANDZURA_25,  3*n*sizeof(gdouble)) ; break ;
  case 54:  memcpy(rule, WANDZURA_54,  3*n*sizeof(gdouble)) ; break ;
  case 85:  memcpy(rule, WANDZURA_85,  3*n*sizeof(gdouble)) ; break ;
  case 126: memcpy(rule, WANDZURA_126, 3*n*sizeof(gdouble)) ; break ;
  case 175: memcpy(rule, WANDZURA_175, 3*n*sizeof(gdouble)) ; break ;
  default:
    g_log(G_LOG_DOMAIN,
  	  G_LOG_LEVEL_ERROR,
  	  "%s: %d point Wandzura-Xiao quadrature is not available",
  	  __FUNCTION__, n) ;
    break ;
  }

  for ( i = 0 ; i < n ; i ++ ) 
    _quadrature_rule_remap(xi0, eta0, xi1, eta1, xi2, eta2,
			   rule[3*i+0], rule[3*i+1], rule[3*i+2], 
			   &(rule[3*i+0]), &(rule[3*i+1]), &(rule[3*i+2])) ;

  return BEM3D_SUCCESS ;
}

void rotation_indices(gint nc, gint rot[])

     /*
       Generate indices for cycling round the corners of an element.
      */

{
  gint i ;
  
  for ( i = 0 ; i < nc-1 ; i ++ ) rot[i] = i+1 ;
  rot[nc-1] = 0 ; rot[nc] = 1 ;

  return ;
}

/** 
 * Allocate a new quadrature rule
 * 
 * @param n maximum number of points in quadrature rule;
 * @param nc number of components in integral.
 * 
 * @return pointer to new quadrature rule
 */

BEM3DQuadratureRule *bem3d_quadrature_rule_new(gint n, gint nc)

{
  BEM3DQuadratureRule *q ;

  g_log(G_LOG_DOMAIN,
	G_LOG_LEVEL_DEBUG,
	"%s: allocating for %d quadrature points", 
	__FUNCTION__, n) ;

  g_return_val_if_fail(n >= 0, NULL) ;
  g_return_val_if_fail(nc > 0, NULL) ;

  q = (BEM3DQuadratureRule *)g_malloc(sizeof(BEM3DQuadratureRule)) ;
  q->nmax = n ; q->n = 0 ; q->nc = nc ;

  q->nfree = 0 ; q->nfree_max = 32 ; q->wfree = 1 ;

  if ( n != 0 )
    q->rule = (gdouble *)g_malloc(3*nc*n*sizeof(gdouble)) ;
  else
    q->rule = NULL ;

  return q ;
}

gint bem3d_quadrature_rule_free(BEM3DQuadratureRule *q)

{
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  g_free(q->rule) ;

  g_free(q) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Reallocate a quadrature rule to hold a new number of quadrature
 * points. Note that the previous data in the rule are lost.
 * 
 * @param q quadrature rule to reallocate
 * @param n number of points in reallocated rule
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_quadrature_rule_realloc(BEM3DQuadratureRule *q, gint n)

{
  g_return_val_if_fail(n >= 0, BEM3D_ARGUMENT_OUT_OF_RANGE) ;

  if ( q->nmax >= n ) return BEM3D_SUCCESS ;

  g_log(G_LOG_DOMAIN,
	G_LOG_LEVEL_DEBUG,
	"%s: reallocating for %d quadrature points", 
	__FUNCTION__, n) ;

  g_free(q->rule) ;
  q->nmax = n ; q->n = 0 ;
  q->rule = (gdouble *)g_malloc(3*n*(q->nc)*sizeof(gdouble)) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Add a new quadrature point to a rule
 * 
 * @param q quadrature rule
 * @param xi area coordinate of quadrature point
 * @param eta area coordinate of quadrature point
 * @param w weight of quadrature point
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_quadrature_add_point(BEM3DQuadratureRule *q,
				gdouble xi, gdouble eta,
				gdouble w)

{
  g_assert( bem3d_quadrature_vertex_number(q) < 
	    bem3d_quadrature_vertex_number_max(q) ) ;

  return _quadrature_add_point(q, xi, eta, w) ;

}

/** 
 * Quadrature selection parameter used in choosing quadrature rules
 * 
 * @param p field point
 * @param e element over which to integrate
 * 
 * @return sigma parameter.
 */

gdouble bem3d_quadrature_parameter(GtsPoint *p, BEM3DElement *e)
  
{
  gdouble Rs, sc = M_SQRT2, Rp, len, s[3], t[3], n[3] ;
  gint i ;
  GtsPoint y ;

/*   g_debug("%s:", __FUNCTION__) ; */

  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ )
    if ( nodes_coincide(p, bem3d_element_vertex(e,i)) ) return 0.0 ;

  y.x = y.y = y.z = 0.0 ;
  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) {
    y.x += GTS_POINT(bem3d_element_vertex(e,i))->x ;
    y.y += GTS_POINT(bem3d_element_vertex(e,i))->y ;
    y.z += GTS_POINT(bem3d_element_vertex(e,i))->z ;
  }
  y.x /= bem3d_element_vertex_number(e) ;
  y.y /= bem3d_element_vertex_number(e) ;
  y.z /= bem3d_element_vertex_number(e) ;

  Rs = 0.0 ;
  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ )
    Rs = MAX(gts_point_distance2(bem3d_element_vertex(e,i), &y), Rs) ;

  Rs = sc*sqrt(Rs) ;

  Rp = gts_point_distance(p, &y) ;

  if ( Rp > Rs ) return Rp/Rs ;

  triangle_axes(&(GTS_POINT(bem3d_element_corner(e,0))->x),
		&(GTS_POINT(bem3d_element_corner(e,1))->x),
		&(GTS_POINT(bem3d_element_corner(e,2))->x),
		s, t, n) ;
  gts_vector_init(s, p, GTS_POINT(bem3d_element_corner(e,0))) ;
  len = fabs(gts_vector_scalar(s,n)) ;

  return len/Rs ;
}

/** 
 * Set quadrature rule to default
 * 
 * @param p field point;
 * @param e element to integrate over;
 * @param q quadrature rule;
 * @param gfunc ignored;
 * @param param ignored;
 * @param data ignored.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_quadrature_rule_default(GtsPoint *p,
				   BEM3DElement *e,
				   BEM3DQuadratureRule *q,
				   BEM3DGreensFunction *gfunc,
				   BEM3DParameters *param,
				   gpointer data)

{
  gpointer qd ;
  BEM3DQuadratureRuleFunc f ;
  
/*   g_debug("%s: e=%p; q=%p; data = %p", __FUNCTION__, e, q, data) ; */

  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(p), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(data != NULL, BEM3D_NULL_ARGUMENT) ;

  bem3d_quadrature_select((BEM3DQuadratureSelector *)data,
			  bem3d_quadrature_parameter(p, e), &f, &qd) ;
  f(p, e, q, gfunc, param, qd) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Khayat-Wilton quadrature rule for singular and near-singular
 * integrals on triangular elements (Khayat, Michael A. and Wilton,
 * Donald R., `Numerical evaluation of singular and near-singular
 * potential integrals', IEEE Transactions on Antennas and
 * Propagation, 53(10):3180--3190, 2005, doi: 10.1109/TAP.2005.856342
 * 
 * @param p field point
 * @param e element to integrate over
 * @param q quadrature rule to fill
 * @param gfunc ignored;
 * @param param ignored;
 * @param data gint[2], first element number of quadrature points in radius;
 * second element number of points in angle
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_quadrature_rule_kw(GtsPoint *p, BEM3DElement *e,
			      BEM3DQuadratureRule *q, 
			      BEM3DGreensFunction *gfunc,
			      BEM3DParameters *param,
			      gpointer data)

{
  gdouble A, Ad, x, y, z, xd, len, h1d, xi2, xi3 ;
  GtsVector n, nd, h1, h1h, L1, L2, L3, L1h ;
  GtsVector vt1, vt2, vt3, vt4 ;
  GtsVector xi ;
  GtsPoint *r1, *r2, *r3 ;
  static GtsPoint *r0 = NULL, *xsc = NULL ;
  gdouble du, ubar, dy, ybar ;
  gdouble xL, xU, u, uL, uU ;
  gdouble R0, R ;
  gdouble w ;
  gint i, j, t, nc ;
  gint nr, na ;
  static gint rot[32] ;
  static gqr_rule_t *qu, *qy ;

  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(p), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(data != NULL, BEM3D_NULL_ARGUMENT) ;

  /*number of angular and radial quadrature points*/
  na = ((gint *)data)[0] ; nr = ((gint *)data)[1] ;

  /*number of corners on element*/
  nc = bem3d_element_corner_number(e) ; g_assert(nc == 3) ;
  rotation_indices(nc, rot) ;
  
  bem3d_quadrature_rule_realloc(q, nc*nr*na) ;
  bem3d_quadrature_clear(q) ;

  if ( r0 == NULL ) {
    r0 = gts_point_new(gts_point_class(), 0, 0, 0) ;
    xsc = gts_point_new(gts_point_class(), 0, 0, 0) ;  
    qu = gqr_rule_alloc(32) ; qy = gqr_rule_alloc(32) ;
  }

  /*one-dimensional quadrature selection*/
  qu = gqr_rule_realloc(qu, na) ; qy = gqr_rule_realloc(qy, nr) ;
  gqr_rule_select(qu, GQR_GAUSS_LEGENDRE, na, NULL) ;
  gqr_rule_select(qy, GQR_GAUSS_LEGENDRE, nr, NULL) ;
  gqr_rule_scale(qy, 0, 1, &ybar, &dy) ;

  /*basic, whole-element quantities*/
  gts_vector_init(L1, 
		  GTS_POINT(bem3d_element_corner(e,0)), 
		  GTS_POINT(bem3d_element_corner(e,1))) ;
  gts_vector_init(L2, 
		  GTS_POINT(bem3d_element_corner(e,1)), 
		  GTS_POINT(bem3d_element_corner(e,2))) ;
  gts_vector_cross(n, L1, L2) ;
  len = gts_vector_norm(n) ;
  A = gts_vector_scalar(n, n)/len*0.5 ;
  n[0] /= len ; n[1] /= len ; n[2] /= len ; 

  gts_vector_init(L1, GTS_POINT(bem3d_element_corner(e,1)), p) ;
  z = gts_vector_scalar(n, L1) ;

  /*projection of point onto element plane (thanks Eric)*/
  gts_point_set(r0, 
		GTS_POINT(p)->x - z*n[0],
		GTS_POINT(p)->y - z*n[1],
		GTS_POINT(p)->z - z*n[2]) ;

  for ( t = 0 ; t < 3 ; t ++ ) {
    r1 = r0 ; 
    r2 = GTS_POINT(bem3d_element_corner(e,rot[t])) ;
    r3 = GTS_POINT(bem3d_element_corner(e,rot[t+1])) ;
    gts_vector_init(L1, r2, r3) ;
    gts_vector_init(L2, r3, r1) ;
    gts_vector_init(L3, r1, r2) ;
    gts_vector_cross(nd, L1, L2) ;
    len = gts_vector_norm(nd) ;
    if ( len > 1e-16 ) {
      /*check for observation point directly over triangle node*/
      nd[0] /= len ; nd[1] /= len ; nd[2] /= len ; 
      gts_vector_cross(vt1, L1, L2) ;
      Ad = 0.5*gts_vector_scalar(nd, vt1) ;
      gts_vector_cross(h1, L1, nd) ;
      len = gts_vector_norm(L1) ; len *= len ;
      h1[0] *= 2.0*Ad/len ; h1[1] *= 2.0*Ad/len ;
      h1[2] *= 2.0*Ad/len ; 
      h1d = gts_vector_norm(h1) ;
      h1h[0] = h1[0]/h1d ; h1h[1] = h1[1]/h1d ;
      h1h[2] = h1[2]/h1d ;
      len = gts_vector_norm(L1) ;
      L1h[0] = L1[0]/len ; L1h[1] = L1[1]/len ;
      L1h[2] = L1[2]/len ;
      
      gts_vector_cross(vt1, h1h, L2) ;
      gts_vector_cross(vt2, h1h, L3) ;
      for ( j = 0 ; j < gqr_rule_length(qy) ; j ++ ) {
	y = h1d*(1.0-(ybar+dy*gqr_rule_abscissa(qy,j))) ;
	xd = 1.0 - y/h1d ;
	xL = gts_vector_scalar(nd,vt1)*(1.0-xd) ;
	xU = -gts_vector_scalar(nd,vt2)*(1.0-xd) ;
	R0 = sqrt(y*y+z*z) ;
	uL = asinh(xL/R0) ; uU = asinh(xU/R0) ;
	gqr_rule_scale(qu, uL, uU, &ubar, &du) ;
	for ( i = 0 ; i < gqr_rule_length(qu) ; i ++ ) {
	  u = ubar + du*gqr_rule_abscissa(qu, i) ;
	  x = R0*sinh(u) ;
	  R = sqrt(x*x+y*y+z*z) ;
	  w = gts_vector_scalar(nd,n)*h1d*du*dy*R*0.5/A*
	    gqr_rule_weight(qu, i)*gqr_rule_weight(qy,j) ;
	  vt3[0] = y*h1h[0] - x*L1h[0] ;
	  vt3[1] = y*h1h[1] - x*L1h[1] ;
	  vt3[2] = y*h1h[2] - x*L1h[2] ;
	  gts_vector_cross(vt4, L3, vt3) ;
	  xi3 = gts_vector_scalar(nd,vt4)*0.5/Ad ;
	  xi2 = 1.0 - xd - xi3 ;
	  gts_point_set(xsc,
			xd*GTS_POINT(r1)->x +
			xi2*GTS_POINT(r2)->x +
			xi3*GTS_POINT(r3)->x,
			xd*GTS_POINT(r1)->y +
			xi2*GTS_POINT(r2)->y +
			xi3*GTS_POINT(r3)->y,
			xd*GTS_POINT(r1)->z +
			xi2*GTS_POINT(r2)->z +
			xi3*GTS_POINT(r3)->z) ;
	  bem3d_area_coordinates_tri(bem3d_element_corner(e,0),
				     bem3d_element_corner(e,1),
				     bem3d_element_corner(e,2),
				     xsc, xi) ;
  	  _quadrature_add_point(q, xi[0], xi[1], w) ;
	}
      }
    }
  }

  return BEM3D_SUCCESS ;
}

static gint _bem3d_polar_triangle_rule(gqr_rule_t *gN, 
				       gqr_rule_t *gM, 
				       gdouble *rule)

     /*
       Nodes and weights for quadrature on a unit triangle with
       singularity at the (0,0) vertex.
      */

{
  gdouble th, a, r, S, C ;
  gint i, j, k ;
#if BEM3D_QUADRATURE_CACHING
  static GArray *cache_mn = NULL ;
  static GPtrArray *cache_rules = NULL ;
  gdouble *rule ;
#endif /*BEM3D_QUADRATURE_CACHING*/

#if BEM3D_QUADRATURE_CACHING
  if ( cache_mn == NULL ) {
    cache_mn = g_array_new(TRUE, TRUE, sizeof(gint)) ;
    cache_rules = g_ptr_array_new() ;
  }

  /*look up the rule in the cache*/
  for ( i = 0 ; (i < (cache_mn->len)/2) && 
	  (g_array_index(cache_mn,gint,2*i) != M) &&
	  (g_array_index(cache_mn,gint,2*i+1) != N) ;
	i ++ ) ;
  if ( i < cache_mn->len/2 ) {
    rule = g_ptr_array_index(cache_rules, i) ;
    memcpy(zeta, rule, M*N*sizeof(gdouble)) ;
    memcpy(eta, &(rule[M*N]), M*N*sizeof(gdouble)) ;
    memcpy(w, &(rule[2*M*N]), M*N*sizeof(gdouble)) ;
    return BEM3D_SUCCESS ;
  }
#endif /*BEM3D_QUADRATURE_CACHING*/

  for ( i = 0 ; i < gqr_rule_length(gN) ; i ++ ) {
    th = 0.25*M_PI*(1.0 + gqr_rule_abscissa(gN, i)) ;
    S = sin(th) ; C = cos(th) ;
    a = 1./(S+C) ;
    for ( j = 0 ; j < gqr_rule_length(gM) ; j ++ ) {
      k = i*gqr_rule_length(gM) + j ;
      r = 0.5*a*(1.0+gqr_rule_abscissa(gM, j)) ;
      rule[3*k+0] = C*r ; 
      rule[3*k+1] = S*r ;
      rule[3*k+2] = 0.25*M_PI*a*r*gqr_rule_weight(gN,i)*gqr_rule_weight(gM,j) ;
    }
  }

#if BEM3D_QUADRATURE_CACHING
  /*cache the new rule*/
  g_array_append_val(cache_mn, M) ;
  g_array_append_val(cache_mn, N) ;
  rule = (gdouble *)g_malloc(3*M*N*sizeof(gdouble)) ;
  memcpy(rule, zeta, M*N*sizeof(gdouble)) ;
  memcpy(&rule[M*N], eta, M*N*sizeof(gdouble)) ;
  memcpy(&rule[2*M*N], w, M*N*sizeof(gdouble)) ;
  g_ptr_array_add(cache_rules, rule) ;
#endif /*BEM3D_QUADRATURE_CACHING*/

  return BEM3D_SUCCESS ;
}

static void lock_to_origin(gdouble *x)

{
  if ( fabs(x[0]) < 1e-9 ) x[0] = 0 ;
  if ( fabs(x[1]) < 1e-9 ) x[1] = 0 ;
  return ;
}

static gdouble triangle_level(gdouble *x1, gdouble *x2, gdouble *x3,
			      gdouble *x4, gdouble *x5, gdouble *x6,
			      gdouble *p,
			      gdouble *y1, gdouble *y2, gdouble *y3,
			      gdouble *y4, gdouble *y5, gdouble *y6)

{
  gdouble o[3], n[3], t[3], s[3], len, A[9], Ai[9] ;

  triangle_axes(x1, x2, x3, s, t, n) ;

  o[0] = p[0] - x1[0] ; o[1] = p[1] - x1[1] ; o[2] = p[2] - x1[2] ;
  len = gts_vector_scalar(n, o) ;

  /*origin of coordinates in plane*/
  o[0] = p[0] - len*n[0] ; o[1] = p[1] - len*n[1] ; o[2] = p[2] - len*n[2] ;
  
  A[0] = s[0] ; A[3] = s[1] ; A[6] = s[2] ;
  A[1] = t[0] ; A[4] = t[1] ; A[7] = t[2] ;
  A[2] = n[0] ; A[5] = n[1] ; A[8] = n[2] ;

  invert3x3(Ai, A) ;

  matrix_transform(Ai, o, x1, y1) ;
  matrix_transform(Ai, o, x2, y2) ;
  matrix_transform(Ai, o, x3, y3) ;
  matrix_transform(Ai, o, x4, y4) ;
  matrix_transform(Ai, o, x5, y5) ;
  matrix_transform(Ai, o, x6, y6) ;

  lock_to_origin(y1) ; lock_to_origin(y2) ; lock_to_origin(y3) ;
  lock_to_origin(y4) ; lock_to_origin(y5) ; lock_to_origin(y6) ;

  return len ;
}

static gdouble triangle_level_cubic(gdouble *x1, gdouble *x2, gdouble *x3,
				    gdouble *x4, gdouble *x5, gdouble *x6,
				    gdouble *x7, gdouble *x8, gdouble *x9,
				    gdouble *x10,
				    gdouble *p,
				    gdouble *y1, gdouble *y2, gdouble *y3,
				    gdouble *y4, gdouble *y5, gdouble *y6,
				    gdouble *y7, gdouble *y8, gdouble *y9,
				    gdouble *y10)

{
  gdouble o[3], n[3], t[3], s[3], len, A[9], Ai[9] ;

  triangle_axes(x1, x2, x3, s, t, n) ;

  o[0] = p[0] - x1[0] ; o[1] = p[1] - x1[1] ; o[2] = p[2] - x1[2] ;
  len = gts_vector_scalar(n, o) ;

  /*origin of coordinates in plane*/
  o[0] = p[0] - len*n[0] ; o[1] = p[1] - len*n[1] ; o[2] = p[2] - len*n[2] ;
  
  A[0] = s[0] ; A[3] = s[1] ; A[6] = s[2] ;
  A[1] = t[0] ; A[4] = t[1] ; A[7] = t[2] ;
  A[2] = n[0] ; A[5] = n[1] ; A[8] = n[2] ;

  invert3x3(Ai, A) ;

  matrix_transform(Ai, o, x1, y1) ; matrix_transform(Ai, o, x2, y2) ;
  matrix_transform(Ai, o, x3, y3) ; matrix_transform(Ai, o, x4, y4) ;
  matrix_transform(Ai, o, x5, y5) ; matrix_transform(Ai, o, x6, y6) ;
  matrix_transform(Ai, o, x7, y7) ; matrix_transform(Ai, o, x8, y8) ;
  matrix_transform(Ai, o, x9, y9) ; matrix_transform(Ai, o, x10, y10) ;

  lock_to_origin(y1) ; lock_to_origin(y2) ; lock_to_origin(y3) ;
  lock_to_origin(y4) ; lock_to_origin(y5) ; lock_to_origin(y6) ;
  lock_to_origin(y7) ; lock_to_origin(y8) ; lock_to_origin(y9) ;
  lock_to_origin(y10) ;

  return len ;
}

static gint _quadrature_rule_quadratic_polar(GtsPoint *p, BEM3DElement *e,
					     BEM3DQuadratureRule *q, 
					     BEM3DGreensFunction *gfunc,
					     BEM3DParameters *param,
					     gpointer data,
					     gqr_rule_t *gt, 
					     gqr_rule_t *gr,
					     gboolean hypersingular)

{
  gdouble *x1, *x2, *x3, *x4, *x5, *x6 ;
  gdouble y1[3], y2[3], y3[3], y4[3], y5[3], y6[3], z ;
  gdouble dt, dr, t1, t2, t3 ;
  gint N, M ;

  N = ((gint *)data)[0] ; M = ((gint *)data)[1] ;

  x1 = &(GTS_POINT(bem3d_element_vertex(e,0))->x) ;
  x2 = &(GTS_POINT(bem3d_element_vertex(e,1))->x) ;
  x3 = &(GTS_POINT(bem3d_element_vertex(e,2))->x) ;
  x4 = &(GTS_POINT(bem3d_element_vertex(e,3))->x) ;
  x5 = &(GTS_POINT(bem3d_element_vertex(e,4))->x) ;
  x6 = &(GTS_POINT(bem3d_element_vertex(e,5))->x) ;

  z = triangle_level(x1, x2, x3, x4, x5, x6, &(p->x),
		     y1, y2, y3, y4, y5, y6) ;
  
  triangle_angles(gts_point_distance(bem3d_element_vertex(e,0),
				     bem3d_element_vertex(e,1)),
		  gts_point_distance(bem3d_element_vertex(e,1),
				     bem3d_element_vertex(e,2)),
		  gts_point_distance(bem3d_element_vertex(e,2),
				     bem3d_element_vertex(e,0)),
		  &t1, &t2, &t3) ;
  
  dt = MAX(t1,MAX(t2,t3))/N ;
  dr = MAX(gts_point_distance(bem3d_element_vertex(e,0),
			      bem3d_element_vertex(e,1)),
	   MAX(gts_point_distance(bem3d_element_vertex(e,1),
				  bem3d_element_vertex(e,2)),
	       gts_point_distance(bem3d_element_vertex(e,2),
				  bem3d_element_vertex(e,0))))/M ;

  if ( qtri_quad_rule(y1, y2, y3, y4, y5, y6,
		      4, N, dt, 4, M, dr, gt, gr, hypersingular,
		      &(q->rule[0]), 3, &(q->rule[1]), 3,
		      &(q->rule[2]), 3, &(q->n)) != 0 ) {
    g_error("%s: quadrature error,\n"
	    "  element:\n"
	    "  %e %e %e\n  %e %e %e\n  %e %e %e\n"
	    "  %e %e %e\n  %e %e %e\n  %e %e %e\n"
	    "  field point:\n"
	    "  %e %e %e\n"
	    "  transformed element:\n"
	    "  %e %e %e\n  %e %e %e\n  %e %e %e\n"
	    "  %e %e %e\n  %e %e %e\n  %e %e %e\n"
	    "  transformed field point:\n"
	    "  %e %e %e\n",
	    __FUNCTION__,
	    GTS_POINT(e->v[0])->x, GTS_POINT(e->v[0])->y, GTS_POINT(e->v[0])->z,
	    GTS_POINT(e->v[1])->x, GTS_POINT(e->v[1])->y, GTS_POINT(e->v[1])->z,
	    GTS_POINT(e->v[2])->x, GTS_POINT(e->v[2])->y, GTS_POINT(e->v[2])->z,
	    GTS_POINT(e->v[3])->x, GTS_POINT(e->v[3])->y, GTS_POINT(e->v[3])->z,
	    GTS_POINT(e->v[4])->x, GTS_POINT(e->v[4])->y, GTS_POINT(e->v[4])->z,
	    GTS_POINT(e->v[5])->x, GTS_POINT(e->v[5])->y, GTS_POINT(e->v[5])->z,
	    p->x, p->y, p->z,
	    y1[0], y1[1], y1[2], y2[0], y2[1], y2[2], y3[0], y3[1], y3[2],
	    y4[0], y4[1], y4[2], y5[0], y5[1], y5[2], y6[0], y6[1], y6[2],
	    0., 0., z) ;
  }

  return BEM3D_SUCCESS ;
}

static gint _quadrature_rule_cubic_polar(GtsPoint *p, BEM3DElement *e,
					 BEM3DQuadratureRule *q, 
					 BEM3DGreensFunction *gfunc,
					 BEM3DParameters *param,
					 gpointer data,
					 gqr_rule_t *gt, 
					 gqr_rule_t *gr)

{
  gdouble *x1, *x2, *x3, *x4, *x5, *x6, *x7, *x8, *x9, *x10 ;
  gdouble y1[3], y2[3], y3[3], y4[3], y5[3], y6[3], y7[3], y8[3], y9[3], y10[3];
  gdouble dt, dr, t1, t2, t3, z ;
  gint N, M ;

  N = ((gint *)data)[0] ; M = ((gint *)data)[1] ;

  x1 = &(GTS_POINT(bem3d_element_vertex(e,0))->x) ;
  x2 = &(GTS_POINT(bem3d_element_vertex(e,1))->x) ;
  x3 = &(GTS_POINT(bem3d_element_vertex(e,2))->x) ;
  x4 = &(GTS_POINT(bem3d_element_vertex(e,3))->x) ;
  x5 = &(GTS_POINT(bem3d_element_vertex(e,4))->x) ;
  x6 = &(GTS_POINT(bem3d_element_vertex(e,5))->x) ;
  x7 = &(GTS_POINT(bem3d_element_vertex(e,6))->x) ;
  x8 = &(GTS_POINT(bem3d_element_vertex(e,7))->x) ;
  x9 = &(GTS_POINT(bem3d_element_vertex(e,8))->x) ;
  x10 = &(GTS_POINT(bem3d_element_vertex(e,9))->x) ;

  z = triangle_level_cubic(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
			   &(p->x),
			   y1, y2, y3, y4, y5, y6, y7, y8, y9, y10) ;
  
  triangle_angles(gts_point_distance(bem3d_element_vertex(e,0),
				     bem3d_element_vertex(e,1)),
		  gts_point_distance(bem3d_element_vertex(e,1),
				     bem3d_element_vertex(e,2)),
		  gts_point_distance(bem3d_element_vertex(e,2),
				     bem3d_element_vertex(e,0)),
		  &t1, &t2, &t3) ;
  
  dt = 0.5*MIN(t1,MIN(t2,t3))/N ;
  dr = 0.5*MIN(gts_point_distance(bem3d_element_vertex(e,0),
				  bem3d_element_vertex(e,1)),
	       MIN(gts_point_distance(bem3d_element_vertex(e,1),
				      bem3d_element_vertex(e,2)),
		   gts_point_distance(bem3d_element_vertex(e,2),
				      bem3d_element_vertex(e,0))))/M ;
  /* fprintf(stderr, "dt=%lg; dr=%lg;\n", dt, dr) ; */
  if ( ctri_quad_rule(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10,
		      4, N, dt, 4, M, dr, gt, gr, 
		      &(q->rule[0]), 3, &(q->rule[1]), 3,
		      &(q->rule[2]), 3, &(q->n)) != 0 ) {
    g_error("%s: quadrature error,\n"
	    "  element:\n"
	    "  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n"
	    "  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n"
	    "  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n"
	    "  %1.16e %1.16e %1.16e\n"
	    "  field point:\n"
	    "  %1.16e %1.16e %1.16e\n"
	    "  transformed element:\n"
	    "  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n"
	    "  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n"
	    "  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n  %1.16e %1.16e %1.16e\n"
	    "  %1.16e %1.16e %1.16e\n"
	    "  transformed field point:\n"
	    "  %1.16e %1.16e %1.16e\n",
	    __FUNCTION__,
	    GTS_POINT(e->v[0])->x, GTS_POINT(e->v[0])->y, GTS_POINT(e->v[0])->z,
	    GTS_POINT(e->v[1])->x, GTS_POINT(e->v[1])->y, GTS_POINT(e->v[1])->z,
	    GTS_POINT(e->v[2])->x, GTS_POINT(e->v[2])->y, GTS_POINT(e->v[2])->z,
	    GTS_POINT(e->v[3])->x, GTS_POINT(e->v[3])->y, GTS_POINT(e->v[3])->z,
	    GTS_POINT(e->v[4])->x, GTS_POINT(e->v[4])->y, GTS_POINT(e->v[4])->z,
	    GTS_POINT(e->v[5])->x, GTS_POINT(e->v[5])->y, GTS_POINT(e->v[5])->z,
	    GTS_POINT(e->v[6])->x, GTS_POINT(e->v[6])->y, GTS_POINT(e->v[6])->z,
	    GTS_POINT(e->v[7])->x, GTS_POINT(e->v[7])->y, GTS_POINT(e->v[7])->z,
	    GTS_POINT(e->v[8])->x, GTS_POINT(e->v[8])->y, GTS_POINT(e->v[8])->z,
	    GTS_POINT(e->v[9])->x, GTS_POINT(e->v[9])->y, GTS_POINT(e->v[9])->z,
	    p->x, p->y, p->z,
	    y1[0], y1[1], y1[2], y2[0], y2[1], y2[2], y3[0], y3[1], y3[2],
	    y4[0], y4[1], y4[2], y5[0], y5[1], y5[2], y6[0], y6[1], y6[2],
	    y7[0], y7[1], y7[2], y8[0], y8[1], y8[2], y9[0], y9[1], y9[2],
	    y10[0], y10[1], y10[2],
	    0., 0., z) ;
  }

  return BEM3D_SUCCESS ;
}


/** 
 * Polar transformation quadrature rule for singular and near-singular
 * integrals on triangular elements.
 * 
 * @param p field point;
 * @param e element to integrate over;
 * @param q quadrature rule to fill;
 * @param gfunc ignored;
 * @param param ignored;
 * @param data gint[2], first element number of quadrature points in radius;
 * second element number of points in angle.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_quadrature_rule_polar(GtsPoint *p, BEM3DElement *e,
				 BEM3DQuadratureRule *q, 
				 BEM3DGreensFunction *gfunc,
				 BEM3DParameters *param,
				 gpointer data)

{
  gdouble A ;
  gdouble xi1, xi2, eta1, eta2, ww ;
  gint i, j, k, t, nc ;
  gint M, N ;
  gdouble xi0, eta0, xn, en ;
  gint rot[32] ;
  static gqr_rule_t *gN = NULL, *gM = NULL ;

/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(p), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(data != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( gN == NULL ) {
    gN = gqr_rule_alloc(256) ; gM = gqr_rule_alloc(256) ;
  }

  /*number of corners on element*/
  nc = bem3d_element_corner_number(e) ;

  rotation_indices(nc, rot) ;

  /*number of radial and angular quadrature points*/
  N = ((gint *)data)[0] ; M = ((gint *)data)[1] ;
  if ( N <= 0 ) 
    g_error("%s: number of angular quadrature points (%d) must be greater "
	    "than zero", __FUNCTION__, N) ;
  if ( M <= 0 ) 
    g_error("%s: number of radial quadrature points (%d) must be greater "
	    "than zero", __FUNCTION__, M) ;

  bem3d_quadrature_rule_realloc(q, 32*M*N) ;
  bem3d_quadrature_clear(q) ;
  
  gN = gqr_rule_realloc(gN, N) ; gM = gqr_rule_realloc(gM, M) ;
  gqr_rule_select(gN, GQR_GAUSS_LEGENDRE, N, NULL) ;
  gqr_rule_select(gM, GQR_GAUSS_LEGENDRE, M, NULL) ;

  /*polar rule for a second-order triangle*/
  if ( nc == 3 && bem3d_element_vertex_number(e) == 6 )
    return _quadrature_rule_quadratic_polar(p, e, q, gfunc, param, data,
					    gM, gN, FALSE) ;
  if ( nc == 3 && bem3d_element_vertex_number(e) == 10 )
    return _quadrature_rule_cubic_polar(p, e, q, gfunc, param, data,
					gM, gN) ;

  /*location of quadrature point*/
  xi0 = eta0 = 0.0 ;
  i = bem3d_element_nearest_point(e, p, &xi0, &eta0, FALSE) ;

  for ( t = 0 ; t < nc ; t ++ ) {
    k = bem3d_element_corner_index(e, rot[t]) ;
    xi1 = bem3d_element_vertex_xi(e,k) ;
    eta1 = bem3d_element_vertex_eta(e,k) ;

    k = bem3d_element_corner_index(e, rot[t+1]) ;
    xi2 = bem3d_element_vertex_xi(e,k) ;
    eta2 = bem3d_element_vertex_eta(e,k) ;
    A = (xi1-xi0)*(eta2-eta0) - (eta1-eta0)*(xi2-xi0) ;

    if ( fabs(A) > 1e-9 ) {
      i = bem3d_quadrature_vertex_number(q) ;
      _bem3d_polar_triangle_rule(gN, gM, &(q->rule[3*i])) ;
      for ( j = 0 ; j < M*N ; j ++ ) {
	_quadrature_rule_remap(xi0, eta0, xi1, eta1,
			       xi2, eta2, 
			       q->rule[3*(i+j)+0],
			       q->rule[3*(i+j)+1],
			       0.5*q->rule[3*(i+j)+2],
			       &xn, &en, &ww) ;
 	_quadrature_add_point(q, xn, en, ww) ;
      }
    }
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Polar transformation quadrature rule for hypersingular.
 * 
 * @param p field point;
 * @param e element to integrate over;
 * @param q quadrature rule to fill;
 * @param gfunc ignored;
 * @param param ignored;
 * @param data gint[2], first element number of quadrature points in radius;
 * second element number of points in angle.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_quadrature_rule_polar_hs(GtsPoint *p, BEM3DElement *e,
				    BEM3DQuadratureRule *q, 
				    BEM3DGreensFunction *gfunc,
				    BEM3DParameters *param,
				    gpointer data)

{
  gdouble A ;
  gdouble xi1, xi2, eta1, eta2, ww ;
  gint i, j, k, t, nc ;
  gint M, N ;
  gdouble xi0, eta0, xn, en ;
  gint rot[32] ;
  static gqr_rule_t *gN = NULL, *gM = NULL ;

/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(p), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(data != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( gN == NULL ) {
    gN = gqr_rule_alloc(256) ; gM = gqr_rule_alloc(256) ;
  }

  /*number of corners on element*/
  nc = bem3d_element_corner_number(e) ;

  rotation_indices(nc, rot) ;

  /*number of radial and angular quadrature points*/
  N = ((gint *)data)[0] ; M = ((gint *)data)[1] ;
  if ( N <= 0 ) 
    g_error("%s: number of angular quadrature points (%d) must be greater "
	    "than zero", __FUNCTION__, N) ;
  if ( M <= 0 ) 
    g_error("%s: number of radial quadrature points (%d) must be greater "
	    "than zero", __FUNCTION__, M) ;

  bem3d_quadrature_rule_realloc(q, 32*M*N) ;
  bem3d_quadrature_clear(q) ;
  
  gN = gqr_rule_realloc(gN, N) ; gM = gqr_rule_realloc(gM, M) ;
  gqr_rule_select(gN, GQR_GAUSS_LEGENDRE, 4*N, NULL) ;
  gqr_rule_select(gM, GQR_GAUSS_LEGENDRE, M, NULL) ;

  /*polar rule for a second-order triangle*/
  if ( nc == 3 && bem3d_element_vertex_number(e) == 6 )
    return _quadrature_rule_quadratic_polar(p, e, q, gfunc, param, data,
					    gM, gN, TRUE) ;
  if ( nc == 3 && bem3d_element_vertex_number(e) == 10 )
    return _quadrature_rule_cubic_polar(p, e, q, gfunc, param, data,
					gM, gN) ;
  g_assert_not_reached() ;

  /*location of quadrature point*/
  xi0 = eta0 = 0.0 ;
  i = bem3d_element_nearest_point(e, p, &xi0, &eta0, FALSE) ;

  for ( t = 0 ; t < nc ; t ++ ) {
    k = bem3d_element_corner_index(e, rot[t]) ;
    xi1 = bem3d_element_vertex_xi(e,k) ;
    eta1 = bem3d_element_vertex_eta(e,k) ;

    k = bem3d_element_corner_index(e, rot[t+1]) ;
    xi2 = bem3d_element_vertex_xi(e,k) ;
    eta2 = bem3d_element_vertex_eta(e,k) ;
    A = (xi1-xi0)*(eta2-eta0) - (eta1-eta0)*(xi2-xi0) ;

    if ( fabs(A) > 1e-9 ) {
      i = bem3d_quadrature_vertex_number(q) ;
      _bem3d_polar_triangle_rule(gN, gM, &(q->rule[3*i])) ;
      for ( j = 0 ; j < M*N ; j ++ ) {
	_quadrature_rule_remap(xi0, eta0, xi1, eta1,
			       xi2, eta2, 
			       q->rule[3*(i+j)+0],
			       q->rule[3*(i+j)+1],
			       0.5*q->rule[3*(i+j)+2],
			       &xn, &en, &ww) ;
 	_quadrature_add_point(q, xn, en, ww) ;
      }
    }
  }

  return BEM3D_SUCCESS ;
}


/** 
 * Symmetric Gaussian quadrature rules for triangles, taken from
 * Wandzura, S. and Xiao, H., `Symmetric quadrature rules on a
 * triangle', Computers and Mathematics with Applications,
 * 45:1829--1840. If the element \a e is triangular, the rules are
 * passed unmodified; if the element has more than three corners, the
 * rules are mapped on to subtriangles formed by the element corners
 * and the centroid of the element in intrinsic coordinates.
 * 
 * @param p field point (ignored);
 * @param e element to integrate over;
 * @param q quadrature rule to fill;
 * @param gfunc ignored;
 * @param param ignored;
 * @param data gint *, number of points in rule (7, 25, 54, 85, 
 * 126, 175).
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_quadrature_rule_wx(GtsPoint *p, BEM3DElement *e,
			      BEM3DQuadratureRule *q, 
			      BEM3DGreensFunction *gfunc,
			      BEM3DParameters *param,
			      gpointer data)

{
  gint n, i, k, nc, rot[32], nq ;
  gdouble xi0, eta0, xi1, eta1, xi2, eta2, A ;

/*   g_debug("%s: p=%p; e=%p", __FUNCTION__, p, e) ; */

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(data != NULL, BEM3D_NULL_ARGUMENT) ;

  nc = bem3d_element_corner_number(e) ;
  n = *(gint *)data ;

  if ( nc == 3 ) {
    rotation_indices(nc, rot) ;

    k = bem3d_element_corner_index(e, rot[0]) ;
    xi0 = bem3d_element_vertex_xi(e, k) ;
    eta0 = bem3d_element_vertex_eta(e, k) ;

    k = bem3d_element_corner_index(e, rot[1]) ;
    xi1 = bem3d_element_vertex_xi(e, k) ;
    eta1 = bem3d_element_vertex_eta(e, k) ;

    k = bem3d_element_corner_index(e, rot[2]) ;
    xi2 = bem3d_element_vertex_xi(e, k) ;
    eta2 = bem3d_element_vertex_eta(e, k) ;

    bem3d_quadrature_rule_realloc(q, n) ;
    bem3d_quadrature_clear(q) ;

    _bem3d_rule_fill_wx(n, xi0, eta0, xi1, eta1, xi2, eta2, q->rule) ;
    q->n = n ;
  } else {
    bem3d_quadrature_rule_realloc(q, n*nc) ;
    rotation_indices(nc, rot) ;
    for ( (xi0 = eta0 = 0.0), (i = 0) ; i < nc ; i ++ ) {
      xi0 += bem3d_element_vertex_xi(e,i) ;
      eta0 += bem3d_element_vertex_eta(e,i) ;
    }
    xi0 /= (gdouble)nc ; eta0 /= (gdouble)nc ;
    for ( (i = 0), (nq = 0) ; i < nc ; i ++ ) {

      k = bem3d_element_corner_index(e, rot[i]) ;
      xi1 = bem3d_element_vertex_xi(e, k) ;
      eta1 = bem3d_element_vertex_eta(e, k) ;

      k = bem3d_element_corner_index(e, rot[i+1]) ;
      xi2 = bem3d_element_vertex_xi(e, k) ;
      eta2 = bem3d_element_vertex_eta(e, k) ;

      A = (xi1-xi0)*(eta2-eta0) - (eta1-eta0)*(xi2-xi0) ;
      if ( fabs(A) > 1e-16 ) {
	_bem3d_rule_fill_wx(n, xi0, eta0, xi1, eta1, xi2, eta2,
			    &(q->rule[3*nq])) ;
	nq += n ;
      }
    }
    q->n = nq ;
  }
  
  return BEM3D_SUCCESS ;
}

/** 
 * Remap a node and weight of a quadrature rule to a new base
 * triangle. The main use for the function is in the decomposition of
 * elements into triangles and the remapping of those triangles to the
 * original element.
 * 
 * @param xi0 coordinate of first corner of triangle;
 * @param eta0 coordinate of first corner of triangle;
 * @param xi1 coordinate of second corner of triangle;
 * @param eta1 coordinate of second corner of triangle;
 * @param xi2 coordinate of third corner of triangle;
 * @param eta2 coordinate of third corner of triangle;
 * @param s coordinate to remap;
 * @param t coordinate to remap;
 * @param w weight of node;
 * @param sn on exit, remapped coordinate;
 * @param tn on exit, remapped coordinate;
 * @param wn on exit, new weight of node.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_quadrature_rule_remap(gdouble xi0, gdouble eta0,
				 gdouble xi1, gdouble eta1,
				 gdouble xi2, gdouble eta2,
				 gdouble s, gdouble t, gdouble w,
				 gdouble *sn, gdouble *tn,
				 gdouble *wn)
{
  g_return_val_if_fail(sn != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(tn != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(wn != NULL, BEM3D_NULL_ARGUMENT) ;

  return _quadrature_rule_remap(xi0, eta0, xi1, eta1, xi2, eta2, 
				s, t, w, sn, tn, wn) ;

}

/** 
 * Write the nodes and weights of a ::BEM3DQuadratureRule to a file,
 * mainly of use in debugging.
 * 
 * @param q ::BEM3DQuadratureRule to write;
 * @param f file pointer.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_quadrature_rule_write(BEM3DQuadratureRule *q, FILE *f)

{
  gint i ;

/*   g_debug("%s: ", __FUNCTION__) ; */

  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  for ( i = 0 ; i < bem3d_quadrature_vertex_number(q) ; i ++ ) 
    fprintf(f, "%1.16e %1.16e %1.16e\n", 
	    bem3d_quadrature_xi(q,i),
	    bem3d_quadrature_eta(q,i),
	    bem3d_quadrature_weight(q,i)) ;

  return BEM3D_SUCCESS ;
}


gdouble bem3d_quadrature_rule_sum_weights(BEM3DQuadratureRule *q)

{
  gint i ;
  gdouble s ;

  for ( (s = 0.0), (i = 0) ; i < bem3d_quadrature_vertex_number(q) ; 
	i ++ ) s += bem3d_quadrature_weight(q,i) ;
  return s ;
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
 * second element number of points in angle.
 * 
 * @return ::BEM3D_SUCCESS on success. 
 */

gint bem3d_quadrature_rule_hayami(GtsPoint *xs, BEM3DElement *e,
				  BEM3DQuadratureRule *q, 
				  BEM3DGreensFunction *gfunc,
				  BEM3DParameters *param,
				  gpointer data)

{
  static GtsPoint *x = NULL, *xst = NULL, *fj = NULL ;
  GtsPoint *xj, *xjp1 ;
  gdouble d, L[32] ;
  gdouble hj, aj, dtheta ;
  BEM3DShapeFunc shfunc = NULL, shfsub = NULL ;
  gint i, j, k, nc ;
  gint M, N ;
  gdouble xi0, eta0, xi1, eta1, xi2, eta2 ;
  gdouble xn, en, ww ;
  static gint rot[32] ;
  static gqr_rule_t *qn, *qm ;

  g_debug("%s: ", __FUNCTION__) ;

  /* g_error("%s: untested code", __FUNCTION__) ; */

  g_return_val_if_fail(xs != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(xs), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( x == NULL ) {
    x = gts_point_new(gts_point_class(), 0, 0, 0) ;
    xst = gts_point_new(gts_point_class(), 0, 0, 0) ;
    fj = gts_point_new(gts_point_class(), 0, 0, 0) ;
    qm = gqr_rule_alloc(32) ; qn = gqr_rule_alloc(32) ;
  }

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

    edge_nearest_point(xj, xjp1, xst, fj) ;
    hj = gts_point_distance(xst, fj) ;
    if ( hj > 1e-9 ) {
      dtheta = point_internal_angle(xj, xst, xjp1) ;
      aj = point_internal_angle(xj, xst, fj) ;
      i = bem3d_quadrature_vertex_number(q) ;
      _bem3d_hayami_triangle_rule(hj, aj, d, dtheta, qn, qm, q->rule) ;
      for ( k = 0 ; k < N*M ; k ++ ) {
	_quadrature_rule_remap(xi0, eta0, xi1, eta1,
			       xi2, eta2,
			       q->rule[3*(i+k)+0],
			       q->rule[3*(i+k)+1],
			       q->rule[3*(i+k)+2],		    
			       &xn, &en, &ww) ;
 	_quadrature_add_point(q, xn, en, ww) ;
      }
    }
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Generate a quadrature rule for linear elements (triangular for now)
 * using the method of Newman, J. N., `Distributions of sources and
 * normal dipoles over a quadrilateral panel', Journal of Engineering
 * Mathematics, 20:113--126, 1986. The quadrature gives the weights
 * for integration of the Laplace Green's function in the
 * ::BEM3DQuadratureRule free terms so that, for example, \f$\int
 * L_{i}/4\pi R\mathrm{d}S = f_{i}\f$, with \f$f_{i}\f$ the \f$i\f$th
 * quadrature free term.
 * 
 * @param xs field point for the integration;
 * @param e element for integration;
 * @param q quadrature rule;
 * @param gfunc this should contain ::bem3d_greens_func_laplace or NULL;
 * @param param ignored;
 * @param data ignored.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_quadrature_rule_newman(GtsPoint *xs, BEM3DElement *e,
				  BEM3DQuadratureRule *q, 
				  BEM3DGreensFunction *gfunc,
				  BEM3DParameters *param,
				  gpointer data)

{
  gdouble p[3], x1[3], x2[3], x3[3] ;
  GtsVector s, t, n, r ;
  gdouble orient ;

  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(xs != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(xs), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( bem3d_element_node_number(e) != 3 &&
       bem3d_element_node_number(e) != 4 )
    g_error("%s: only implemented for linear triangular elements", 
	    __FUNCTION__) ;
  g_assert(bem3d_element_node_number(e) == 3) ;

  if ( (bem3d_greens_function_func(gfunc) != bem3d_greens_func_laplace) &&
       (bem3d_greens_function_func(gfunc) != NULL ) )
    g_error("%s: Green's function should be "
	    "bem3d_greens_func_laplace or NULL", __FUNCTION__) ;

  if ( !bem3d_greens_function_is_real(gfunc) ) 
    g_error("%s: only implemented for real quadratures", __FUNCTION__) ;

  bem3d_quadrature_clear(q) ;

  gts_vector_init(s, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  GTS_POINT(bem3d_element_node(e,1))) ;
  x2[0] = gts_vector_norm(s) ; x2[1] = 0.0 ;
  gts_vector_init(t, 
		  GTS_POINT(bem3d_element_node(e,1)),
		  GTS_POINT(bem3d_element_node(e,2))) ;
  gts_vector_cross(n, s, t) ;
  gts_vector_normalize(n) ; gts_vector_normalize(s) ;
  gts_vector_cross(t, n, s) ;
  x1[0] = x1[1] = x1[2] = x2[2] = x3[2] = 0.0 ;
  
  gts_vector_init(r, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  GTS_POINT(bem3d_element_node(e,2))) ;
  x3[0] = gts_vector_scalar(r,s) ;
  x3[1] = gts_vector_scalar(r,t) ;
  gts_vector_init(r, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  xs) ;
  p[0] = gts_vector_scalar(r,s) ;
  p[1] = gts_vector_scalar(r,t) ;
  p[2] = gts_vector_scalar(r,n) ;

  if ( nodes_coincide(xs, bem3d_element_node(e,0)) ) {
    p[0] = p[1] = p[2] = 0.0 ;
  }
  if ( nodes_coincide(xs, bem3d_element_node(e,1)) ) {
    p[0] = x2[0] ; p[1] = x2[1] ; p[2] = 0.0 ;
  }
  if ( nodes_coincide(xs, bem3d_element_node(e,2)) ) {
    p[0] = x3[0] ; p[1] = x3[1] ; p[2] = 0.0 ;
  }

  orient = gts_point_orientation_3d(bem3d_element_vertex(e,0),
				    bem3d_element_vertex(e,1),
				    bem3d_element_vertex(e,2),
				    xs) ;
  if ( orient == 0 ) p[2] = 0.0 ;

  if ( e->Imn != NULL && bem3d_element_moment_order(e) > 0 )
    newman_tri_shape(p, x1, x2, x3, 
		     &(bem3d_element_moment(e,0,0)),
		     bem3d_element_moment_order(e),
		     q->free_g, q->free_dg) ;
  else
    newman_tri_shape(p, x1, x2, x3, NULL, 0, q->free_g, q->free_dg) ;

  bem3d_quadrature_free_number(q) = 3 ;
  q->free_g[0] *= -0.25*M_1_PI ; q->free_g[1] *= -0.25*M_1_PI ;
  q->free_g[2] *= -0.25*M_1_PI ;

  q->free_dg[0] *= -0.25*M_1_PI ; q->free_dg[1] *= -0.25*M_1_PI ;
  q->free_dg[2] *= -0.25*M_1_PI ;

  return BEM3D_SUCCESS ;
}


/** 
 * Generate a quadrature rule for gradient of potential field using
 * linear elements (triangular for now) using the method of Newman,
 * J. N., `Distributions of sources and normal dipoles over a
 * quadrilateral panel', Journal of Engineering Mathematics,
 * 20:113--126, 1986. The quadrature gives the weights for integration
 * of the Laplace Green's function in the ::BEM3DQuadratureRule free
 * terms so that, for example, \f$\nabla\int L_{i}/4\pi R\mathrm{d}S =
 * f_{i}\f$, with \f$f_{i}\f$ the \f$i\f$th quadrature free term. The
 * ordering of the free terms \f$g_{i}\f$ is such that \f$\int G
 * \phi\mathrm{d}S=\sum_{i=0}^{2}g_{i}\phi_{i}\f$,
 * \f$\frac{\partial}{\partial x}\int G
 * \phi\mathrm{d}S=\sum_{i=3}^{5}g_{i}\phi_{i}\f$,
 * \f$\frac{\partial}{\partial y}\int G
 * \phi\mathrm{d}S=\sum_{i=6}^{8}g_{i}\phi_{i}\f$, and
 * \f$\frac{\partial}{\partial z}\int G
 * \phi\mathrm{d}S=\sum_{i=9}^{11}g_{i}\phi_{i}\f$, and likewise for
 * the double layer term.
 * 
 * @param xs field point for the integration;
 * @param e element for integration;
 * @param q quadrature rule;
 * @param gfunc this should contain ::bem3d_greens_func_laplace or NULL;
 * @param param ignored;
 * @param data ignored.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_quadrature_rule_newman_gradient(GtsPoint *xs, BEM3DElement *e,
					   BEM3DQuadratureRule *q, 
					   BEM3DGreensFunction *gfunc,
					   BEM3DParameters *param,
					   gpointer data)

{
  gdouble p[3], x1[3], x2[3], x3[3] ;
  GtsVector s, t, n, r ;
  gdouble orient, G[12], dG[12] ;
  gint i ;

  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(xs != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(xs), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( bem3d_element_node_number(e) != 3 &&
       bem3d_element_node_number(e) != 4 )
    g_error("%s: only implemented for linear triangular elements", 
	    __FUNCTION__) ;

  g_assert(bem3d_element_node_number(e) == 3) ;

  if ( (bem3d_greens_function_func(gfunc) != bem3d_greens_func_laplace) &&
       (bem3d_greens_function_func(gfunc) != NULL ) )
    g_error("%s: Green's function should be "
	    "bem3d_greens_func_laplace or NULL", __FUNCTION__) ;

  if ( !bem3d_greens_function_is_real(gfunc) ) 
    g_error("%s: only implemented for real quadratures", __FUNCTION__) ;

  bem3d_quadrature_clear(q) ;

  gts_vector_init(s, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  GTS_POINT(bem3d_element_node(e,1))) ;
  x2[0] = gts_vector_norm(s) ; x2[1] = 0.0 ;
  gts_vector_init(t, 
		  GTS_POINT(bem3d_element_node(e,1)),
		  GTS_POINT(bem3d_element_node(e,2))) ;
  gts_vector_cross(n, s, t) ;
  gts_vector_normalize(n) ; gts_vector_normalize(s) ;
  gts_vector_cross(t, n, s) ;
  x1[0] = x1[1] = x1[2] = x2[2] = x3[2] = 0.0 ;
  
  gts_vector_init(r, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  GTS_POINT(bem3d_element_node(e,2))) ;
  x3[0] = gts_vector_scalar(r,s) ;
  x3[1] = gts_vector_scalar(r,t) ;
  gts_vector_init(r, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  xs) ;
  p[0] = gts_vector_scalar(r,s) ;
  p[1] = gts_vector_scalar(r,t) ;
  p[2] = gts_vector_scalar(r,n) ;

  if ( nodes_coincide(xs, bem3d_element_node(e,0)) ) {
    p[0] = p[1] = p[2] = 0.0 ;
  }
  if ( nodes_coincide(xs, bem3d_element_node(e,1)) ) {
    p[0] = x2[0] ; p[1] = x2[1] ; p[2] = 0.0 ;
  }
  if ( nodes_coincide(xs, bem3d_element_node(e,2)) ) {
    p[0] = x3[0] ; p[1] = x3[1] ; p[2] = 0.0 ;
  }

  orient = gts_point_orientation_3d(bem3d_element_vertex(e,0),
				   bem3d_element_vertex(e,1),
				   bem3d_element_vertex(e,2),
				   xs) ;
  if ( orient == 0 ) p[2] = 0.0 ;

  newman_tri_shape_gradient(p, x1, x2, x3, NULL, 0, G, dG) ;

  q->free_g [0] =  G[0] ; q->free_g [1] =  G[1] ; q->free_g [2] =  G[2] ;
  q->free_dg[0] = dG[0] ; q->free_dg[1] = dG[1] ; q->free_dg[2] = dG[2] ;

  q->free_g[3]  = s[0]*G[3] + t[0]*G[6] + n[0]*G[ 9] ;
  q->free_g[4]  = s[0]*G[4] + t[0]*G[7] + n[0]*G[10] ;
  q->free_g[5]  = s[0]*G[5] + t[0]*G[8] + n[0]*G[11] ;

  q->free_g[6]  = s[1]*G[3] + t[1]*G[6] + n[1]*G[ 9] ;
  q->free_g[7]  = s[1]*G[4] + t[1]*G[7] + n[1]*G[10] ;
  q->free_g[8]  = s[1]*G[5] + t[1]*G[8] + n[1]*G[11] ;

  q->free_g[ 9] = s[2]*G[3] + t[2]*G[6] + n[2]*G[ 9] ;
  q->free_g[10] = s[2]*G[4] + t[2]*G[7] + n[2]*G[10] ;
  q->free_g[11] = s[2]*G[5] + t[2]*G[8] + n[2]*G[11] ;


  q->free_dg[3]  = s[0]*dG[3] + t[0]*dG[6] + n[0]*dG[ 9] ;
  q->free_dg[4]  = s[0]*dG[4] + t[0]*dG[7] + n[0]*dG[10] ;
  q->free_dg[5]  = s[0]*dG[5] + t[0]*dG[8] + n[0]*dG[11] ;

  q->free_dg[6]  = s[1]*dG[3] + t[1]*dG[6] + n[1]*dG[ 9] ;
  q->free_dg[7]  = s[1]*dG[4] + t[1]*dG[7] + n[1]*dG[10] ;
  q->free_dg[8]  = s[1]*dG[5] + t[1]*dG[8] + n[1]*dG[11] ;

  q->free_dg[9 ] = s[2]*dG[3] + t[2]*dG[6] + n[2]*dG[ 9] ;
  q->free_dg[10] = s[2]*dG[4] + t[2]*dG[7] + n[2]*dG[10] ;
  q->free_dg[11] = s[2]*dG[5] + t[2]*dG[8] + n[2]*dG[11] ;

  bem3d_quadrature_free_number(q) = 12 ;
  for ( i = 0 ; i < 12 ; i ++ ) {
    q->free_g[i] *= -0.25*M_1_PI ; 
    q->free_dg[i] *= -0.25*M_1_PI ;
  }

  bem3d_quadrature_component_number(q) = 4 ;

  return BEM3D_SUCCESS ;
}

gint bem3d_quadrature_rule_decomp(GtsPoint *xs, BEM3DElement *e,
				  BEM3DQuadratureRule *q, 
				  BEM3DGreensFunction *gfunc,
				  BEM3DParameters *param,
				  gpointer data)

{
  gdouble p[3], x1[3], x2[3], x3[3], o12, o23, o31 ;
  GtsVector s, t, n, r ;
  gdouble orient ;

  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(xs != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(xs), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( bem3d_element_node_number(e) != 3 )
    g_error("%s: only implemented for linear triangular elements", 
	    __FUNCTION__) ;

  if ( (bem3d_greens_function_func(gfunc) != bem3d_greens_func_laplace) &&
       (bem3d_greens_function_func(gfunc) != NULL ) )
    g_error("%s: Green's function should be "
	    "bem3d_greens_func_laplace or NULL", __FUNCTION__) ;

  if ( !bem3d_greens_function_is_real(gfunc) ) 
    g_error("%s: only implemented for real quadratures", __FUNCTION__) ;

  bem3d_quadrature_clear(q) ;

  gts_vector_init(s, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  GTS_POINT(bem3d_element_node(e,1))) ;
  x2[0] = gts_vector_norm(s) ; x2[1] = 0.0 ;
  gts_vector_init(t, 
		  GTS_POINT(bem3d_element_node(e,1)),
		  GTS_POINT(bem3d_element_node(e,2))) ;
  gts_vector_cross(n, s, t) ;
  gts_vector_normalize(n) ; gts_vector_normalize(s) ;
  gts_vector_cross(t, n, s) ;
  x1[0] = x1[1] = x1[2] = x2[2] = x3[2] = 0.0 ;
  
  gts_vector_init(r, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  GTS_POINT(bem3d_element_node(e,2))) ;
  x3[0] = gts_vector_scalar(r,s) ;
  x3[1] = gts_vector_scalar(r,t) ;
  gts_vector_init(r, 
		  GTS_POINT(bem3d_element_node(e,0)),
		  xs) ;
  p[0] = gts_vector_scalar(r,s) ;
  p[1] = gts_vector_scalar(r,t) ;
  p[2] = gts_vector_scalar(r,n) ;

  orient = gts_point_orientation_3d(bem3d_element_vertex(e,0),
				    bem3d_element_vertex(e,1),
				    bem3d_element_vertex(e,2),
				    xs) ;
  if ( nodes_coincide(xs, bem3d_element_node(e,0)) ) {
    p[0] = p[1] = p[2] = 0.0 ;
  }
  if ( nodes_coincide(xs, bem3d_element_node(e,1)) ) {
    p[0] = x2[0] ; p[1] = x2[1] ; p[2] = 0.0 ;
  }
  if ( nodes_coincide(xs, bem3d_element_node(e,2)) ) {
    p[0] = x3[0] ; p[1] = x3[1] ; p[2] = 0.0 ;
  }

  if ( orient == 0 || fabs(p[2]) < 1e-3 ) p[2] = 0.0 ;

  o12 = orient2d(p, x1, x2) ;
  o12 = ( o12 != 0 ? (o12 < 0 ? -1 : 1) : 0 ) ;
  o23 = orient2d(p, x2, x3) ;
  o23 = ( o23 != 0 ? (o23 < 0 ? -1 : 1) : 0 ) ;
  o31 = orient2d(p, x3, x1) ;
  o31 = ( o31 != 0 ? (o31 < 0 ? -1 : 1) : 0 ) ;

  triangle_quad_shape(x1, x2, x3, p, o12, o23, o31, 1, FALSE, 
		      q->free_g, q->free_dg, NULL) ;

  bem3d_quadrature_free_number(q) = 3 ;
  q->free_g[0] *= 0.25*M_1_PI ; 
  q->free_g[1] *= 0.25*M_1_PI ;
  q->free_g[2] *= 0.25*M_1_PI ;

  q->free_dg[0] *= 0.25*M_1_PI*p[2] ; 
  q->free_dg[1] *= 0.25*M_1_PI*p[2] ;
  q->free_dg[2] *= 0.25*M_1_PI*p[2] ;

  return BEM3D_SUCCESS ;
}

gint bem3d_quadrature_rule_decomp_gradient(GtsPoint *xs, BEM3DElement *e,
					   BEM3DQuadratureRule *q, 
					   BEM3DGreensFunction *gfunc,
					   BEM3DParameters *param,
					   gpointer data)

{
  gdouble o12, o23, o31 ;
  GtsVertex *v1, *v2, *v3 ;
  gint i ;

  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(xs != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(xs), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( bem3d_element_node_number(e) != 3 )
    g_error("%s: only implemented for linear triangular elements", 
	    __FUNCTION__) ;

  if ( (bem3d_greens_function_func(gfunc) != bem3d_greens_func_laplace) &&
       (bem3d_greens_function_func(gfunc) != NULL ) )
    g_error("%s: Green's function should be "
	    "bem3d_greens_func_laplace or NULL", __FUNCTION__) ;

  if ( !bem3d_greens_function_is_real(gfunc) ) 
    g_error("%s: only implemented for real quadratures", __FUNCTION__) ;

  bem3d_quadrature_clear(q) ;

  gts_triangle_vertices(GTS_TRIANGLE(e->f[0]), &v1, &v2, &v3) ;

  triangle_orientations(&(GTS_POINT(v1)->x),
			&(GTS_POINT(v2)->x),
			&(GTS_POINT(v3)->x),
			&(GTS_POINT(xs)->x),
			&o12, &o23, &o31) ;

  triangle_gradient_quad_shape(&(GTS_POINT(v1)->x),
			       &(GTS_POINT(v2)->x),
			       &(GTS_POINT(v3)->x),
			       &(GTS_POINT(xs)->x),
			       o12, o23, o31,
			       q->free_g, q->free_dg) ;

  bem3d_quadrature_free_number(q) = 9 ;
  for ( i = 0 ; i < 9 ; i ++ ) {
    q->free_g[i]  *= 0.25*M_1_PI ; 
    q->free_dg[i] *= 0.25*M_1_PI ;
  }

  bem3d_quadrature_component_number(q) = 3 ;

  return BEM3D_SUCCESS ;
}

/** 
 * Standard Gaussian quadrature for triangular elements
 * 
 * @param p field point (ignored)
 * @param e element to integrate over
 * @param q quadrature rule to fill
 * @param gfunc ignored;
 * @param param ignored;
 * @param data pointer to gint containing the number of points in rule (1, 
 * 3, 4, or 7)
 * 
 * @return ::BEM3D_SUCCESS on success, exit with error if number of
 * points in rule is not 1, 3, 4, or 7.
 */

gint bem3d_quadrature_rule_gauss(GtsPoint *p, BEM3DElement *e,
				 BEM3DQuadratureRule *q, 
				 BEM3DGreensFunction *gfunc,
				 BEM3DParameters *param,
				 gpointer data)

{
  gint n, nc ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(data != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( (nc = bem3d_element_corner_number(e)) != 3 ) 
    g_error("%s: only triangular elements can use this rule", 
	    __FUNCTION__) ;
  n = *((gint *)data) ;
  bem3d_quadrature_rule_realloc(q, n) ;

  /* gauss_select_quadrature(n, &xk, &wt, GAUSS_TRIANGULAR) ; */

  /* g_error("%s: untested code", __FUNCTION__) ; */
  
  switch (n) {
  default: g_error("%s: no %d point rule available", __FUNCTION__, n) ;
  case 1: memcpy(q->rule, GAUSS_TRIANGLE_1, n*3*sizeof(gdouble)) ; break ;
  case 3: memcpy(q->rule, GAUSS_TRIANGLE_3, n*3*sizeof(gdouble)) ; break ;
  case 4: memcpy(q->rule, GAUSS_TRIANGLE_4, n*3*sizeof(gdouble)) ; break ;
  case 7: memcpy(q->rule, GAUSS_TRIANGLE_7, n*3*sizeof(gdouble)) ; break ;
  }

  bem3d_quadrature_vertex_number(q) = n ;

  return BEM3D_SUCCESS ;
}

/** 
 * Semi-analytical series expansion for Helmholtz potential
 * 
 * @param p field point (ignored)
 * @param e element to integrate over
 * @param q quadrature rule to fill
 * @param gfunc ignored;
 * @param param ignored;
 * @param data pointer ignored.
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_quadrature_rule_series(GtsPoint *xs, BEM3DElement *e,
				  BEM3DQuadratureRule *q, 
				  BEM3DGreensFunction *gfunc,
				  BEM3DParameters *param,
				  gpointer data)

{
  gdouble tol, buf[24], lmr, lmi ;
  gint qmax, Q = 25, i, polar[] = {8, 8} ;
  
  tol = 1e-9 ; qmax = 0 ;
  
  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(xs != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(xs), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( bem3d_greens_function_func(gfunc) == bem3d_greens_func_laplace )
    return bem3d_quadrature_rule_wx(xs, e, q, gfunc, param, &Q) ;
    /* return bem3d_quadrature_rule_polar(xs, e, q, gfunc, param, polar) ; */

  g_assert(bem3d_parameters_wavenumber(param) != 0.0) ;
  
  if ( bem3d_element_vertex_number(e) != 3 )
    g_error("%s: only implemented for triangular elements", 
	    __FUNCTION__) ;

  if ( bem3d_element_node_number(e) != 3
       &&
       bem3d_element_node_number(e) != 1 )
    g_error("%s: only implemented for linear or constant elements", 
	    __FUNCTION__) ;

  if ( (bem3d_greens_function_func(gfunc) != bem3d_greens_func_helmholtz) &&
       (bem3d_greens_function_func(gfunc) != bem3d_greens_func_helmholtz_hs) &&
       (bem3d_greens_function_func(gfunc) != NULL ) )
    g_error("%s: Green's function should be "
	    "bem3d_greens_func_helmholtz or NULL", __FUNCTION__) ;

  if ( bem3d_greens_function_is_real(gfunc) ) 
    g_error("%s: only implemented for complex quadratures", __FUNCTION__) ;

  bem3d_quadrature_clear(q) ;

  if ( bem3d_element_node_number(e) == 3 ) {
    /*linear triangular elements*/
    i = bem3d_element_find_node(e, GTS_VERTEX(xs)) ;

    Q = htri_quad_shape_1(&(GTS_POINT(xs)->x),
			  &(GTS_POINT(bem3d_element_vertex(e,0))->x),
			  &(GTS_POINT(bem3d_element_vertex(e,1))->x),
			  &(GTS_POINT(bem3d_element_vertex(e,2))->x),
			  bem3d_parameters_wavenumber(param),
			  tol, qmax, q->free_g, q->free_dg, buf) ;
  
    if ( Q != 0 )
      return bem3d_quadrature_rule_polar(xs, e, q, gfunc, param, polar) ;

    bem3d_quadrature_free_number(q) = 3 ;

    for ( Q = 0 ; Q < 6 ; Q ++ ) {
      q->free_g[Q] *= 0.25*M_1_PI ;
      g_assert(!isnan(q->free_g[Q])) ;
      q->free_dg[Q] *= -0.25*M_1_PI ;
      g_assert(!isnan(q->free_dg[Q])) ;
    }

    i = bem3d_element_find_node(e, GTS_VERTEX(xs)) ;
    if ( i != -1 ) {
      q->free_dg[0] = q->free_dg[1] =
    	q->free_dg[2] = q->free_dg[3] =
    	q->free_dg[4] = q->free_dg[5] = 0.0 ;
    }
    
    return BEM3D_SUCCESS ;
  }

  if ( bem3d_greens_function_func(gfunc) == bem3d_greens_func_helmholtz ) {
    /*regular Helmholtz problem*/
    /* if ( _bem3d_trace_set(0) && _bem3d_trace_set(1) ) */
    /*   fprintf(stderr, "Hello\n") ; */

    Q = htri_quad_shape_0(&(GTS_POINT(xs)->x),
    			  &(GTS_POINT(bem3d_element_vertex(e,0))->x),
    			  &(GTS_POINT(bem3d_element_vertex(e,1))->x),
    			  &(GTS_POINT(bem3d_element_vertex(e,2))->x),
    			  bem3d_parameters_wavenumber(param),
    			  tol, qmax, q->free_g, q->free_dg, buf) ;
    /* g_assert(Q != -1) ; */
    /* fprintf(stderr, "%d\n", Q) ; */
    if ( Q != 0 )
      return bem3d_quadrature_rule_polar(xs, e, q, gfunc, param, polar) ;

    bem3d_quadrature_free_number(q) = 1 ;

    for ( Q = 0 ; Q < 2 ; Q ++ ) {
      q->free_g[Q] *= 0.25*M_1_PI ;
      g_assert(!isnan(q->free_g[Q])) ;
      q->free_dg[Q] *= -0.25*M_1_PI ;
      g_assert(!isnan(q->free_dg[Q])) ;
    }

    i = bem3d_element_find_node(e, GTS_VERTEX(xs)) ;
    if ( i != -1 ) {
      q->free_dg[0] = q->free_dg[1] = 0.0 ;
    }
    
    return BEM3D_SUCCESS ;
  }

  if ( bem3d_greens_function_func(gfunc) != bem3d_greens_func_helmholtz_hs )
    g_error("%s: if you have arrived here you should be using a "
	    "hypersingular formulation", __FUNCTION__) ;

  /* polar[0] = 8 ; polar[1] = 8 ; */
  
  i = bem3d_element_find_node(e, GTS_VERTEX(xs)) ;
  if ( i == -1 ) {
    /*collocation point not on element, fall back to numerical quadrature*/
    return bem3d_quadrature_rule_polar(xs, e, q, gfunc, param, polar) ;
    /* Q = 25 ; */
    /* return bem3d_quadrature_rule_wx(xs, e, q, gfunc, param, &Q) ;     */
  }

  qmax = 0 ;
  /*buf will contain integrals of G, dGdz, d2Gdz2*/
  Q = htri_quad_shape_0(&(GTS_POINT(xs)->x),
			&(GTS_POINT(bem3d_element_vertex(e,0))->x),
			&(GTS_POINT(bem3d_element_vertex(e,1))->x),
			&(GTS_POINT(bem3d_element_vertex(e,2))->x),
			bem3d_parameters_wavenumber(param),
			tol, qmax, &(buf[0]), &(buf[2]), &(buf[4])) ;

  /* if ( Q != 0 ) */
  /*   return bem3d_quadrature_rule_polar(xs, e, q, gfunc, param, polar) ; */
  
  bem3d_quadrature_free_number(q) = 1 ;
  
  lmr = bem3d_parameters_lambda_real(param) ;
  lmi = bem3d_parameters_lambda_imag(param) ;

  /* q->free_g[0] = buf[0] + (lmr*buf[2] - lmi*buf[3]) ; */
  /* q->free_g[1] = buf[1] + (lmr*buf[3] + lmi*buf[2]) ; */
  /*zero here because at the moment we are only dealing with self term*/
  q->free_g[0] = buf[0] + 0*(lmr*buf[2] - lmi*buf[3]) ;
  q->free_g[1] = buf[1] + 0*(lmr*buf[3] + lmi*buf[2]) ;

  /*negative on d/dn because this should be d/dn1*/
  q->free_dg[0] = -0*buf[2] - (lmr*buf[4] - lmi*buf[5]) ;
  q->free_dg[1] = -0*buf[3] - (lmr*buf[5] + lmi*buf[4]) ;

  q->free_g[0]  *= 0.25*M_1_PI ; q->free_g[1]  *= 0.25*M_1_PI ;
  q->free_dg[0] *= 0.25*M_1_PI ; q->free_dg[1] *= 0.25*M_1_PI ;
  
  return BEM3D_SUCCESS ;
}

static gint quad_subtriangle_mzht(gdouble k, gdouble *y1, gdouble *y2,
				  gqr_rule_t *rule, gdouble *g, gdouble *dg)

{
  gdouble dth, thbar, th, S, C, R, l1, l2, l3, phi ;
  /* static gdouble A = 0.0 ; */
  gint i ;
  
  l1 = sqrt(y1[0]*y1[0] + y1[1]*y1[1]) ;
  l2 = sqrt(y2[0]*y2[0] + y2[1]*y2[1]) ;
  l3 = sqrt((y1[0]-y2[0])*(y1[0]-y2[0]) +
	    (y1[1]-y2[1])*(y1[1]-y2[1])) ;
  
  dth = acos((l1*l1 + l2*l2 - l3*l3)/l1/l2*0.5) ;
  phi = atan2(l1 - l2*cos(dth), l2*sin(dth)) ;

  g_assert(dth > 0.0) ;
  
  gqr_rule_scale(rule, 0, dth, &thbar, &dth) ;
  /* fprintf(stderr, "%lg %lg\n", thbar, dth) ; */

  for ( i = 0 ; i < gqr_rule_length(rule) ; i ++ ) {
    th = thbar + gqr_rule_abscissa(rule, i)*dth ;

    R = l1*cos(phi)/cos(th-phi) ;
    
    S = sin(k*R)*0.25*M_1_PI ; C = cos(k*R)*0.25*M_1_PI ;
    g[0] += S/k*dth*gqr_rule_weight(rule, i) ;
    g[1] -= C/k*dth*gqr_rule_weight(rule, i) ;

    dg[0] -= C/R*dth*gqr_rule_weight(rule, i) ;
    dg[1] -= S/R*dth*gqr_rule_weight(rule, i) ;

    /* A += dth*gqr_rule_weight(rule, i)*R*R*0.5 ; */
  }

  /* fprintf(stderr, "%lg\n", A) ; */
  
  g_assert(!isnan(g[0])) ; g_assert(!isnan(g[1])) ;
  g_assert(!isnan(dg[0])) ; g_assert(!isnan(dg[1])) ;

  return 0 ;
}


/*
  Quadrature rule based on Matsumoto, Zheng, Harada, Takahashi, 2010,
  doi:10.1299/jcst.4.194
*/

gint bem3d_quadrature_rule_mzht(GtsPoint *xs, BEM3DElement *e,
				BEM3DQuadratureRule *q, 
				BEM3DGreensFunction *gfunc,
				BEM3DParameters *param,
				gpointer data)

{
  gint i, ngp, polar[2] ;
  gdouble *x1, *x2, *x3, *xf, lmr, lmi, alr, ali ;
  gdouble y1[3], y2[3], y3[3], yf[3], s[3], t[3], n[3], og[3], k ;
  GtsPoint *v ;
  gqr_rule_t *rule ;
  
  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(xs != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(xs), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;

  /* bem3d_quadrature_clear(q) ; */

  if ( bem3d_element_vertex_number(e) != 3 )
    g_error("%s: only implemented for triangular elements", 
	    __FUNCTION__) ;
  if ( bem3d_element_node_number(e) != 1 )
    g_error("%s: only implemented for zero-order elements", 
	    __FUNCTION__) ;

  /*check for collocation point on element*/
  i = bem3d_element_find_node(e, GTS_VERTEX(xs)) ;
  ngp = 25 ; polar[0] = 16 ; polar[1] = 16 ;
  k = bem3d_parameters_wavenumber(param) ;
  if ( (i == -1) ) {
    /*collocation point not on element, fall back to numerical quadrature*/
    /* return bem3d_quadrature_rule_wx(xs, e, q, gfunc, param, &ngp) ; */
    v = GTS_POINT(bem3d_element_node(e,0)) ;
    if ( gts_point_distance(v, GTS_POINT(xs)) > 0.125 || k == 0.0 ) 
      return bem3d_quadrature_rule_wx(xs, e, q, gfunc, param, &ngp) ;
    return bem3d_quadrature_rule_polar(xs, e, q, gfunc, param, polar) ;
    return bem3d_quadrature_rule_hayami(xs, e, q, gfunc, param, polar) ;
  }

  if ( k == 0.0 ) {
    return bem3d_quadrature_rule_polar(xs, e, q, gfunc, param, polar) ;
  }
  
  bem3d_quadrature_clear(q) ;
  /* memset(q->free_g, 0, 32*sizeof(gdouble)) ; */
  /* memset(q->free_dg, 0, 32*sizeof(gdouble)) ; */
  /* memset(q->rule, 0, sizeof(gdouble)*(q->n)*3) ; */
  ngp = 8 ;
  
  x1 = &(GTS_POINT(bem3d_element_vertex(e,0))->x) ;
  x2 = &(GTS_POINT(bem3d_element_vertex(e,1))->x) ;
  x3 = &(GTS_POINT(bem3d_element_vertex(e,2))->x) ;
  xf = &(GTS_POINT(xs)->x) ;

  htri_triangle_axes(x1, x2, x3, xf, og, s, t, n) ;
  htri_triangle_project(og, s, t, n, x1, y1) ;
  htri_triangle_project(og, s, t, n, x2, y2) ;
  htri_triangle_project(og, s, t, n, x3, y3) ;
  htri_triangle_project(og, s, t, n, xf, yf) ;

  /* fprintf(stderr, "%lg %lg %lg %lg\n", y1[2], y2[2], y3[2], yf[2]) ; */
  /* fprintf(stderr, "%lg %lg %lg %lg %lg %lg\n", */
  /* 	  xf[0], xf[1], xf[2], og[0], og[1], og[2]) ; */
  
  rule = gqr_rule_alloc(ngp) ;
  gqr_rule_select(rule, GQR_GAUSS_LEGENDRE, ngp, NULL) ;

  g_assert(k != 0.0) ;
  
  bem3d_quadrature_free_number(q) = 1 ;
  s[0] = 0.0 ; s[1] = 0.5/k ;
  og[0] = 0.0 ; og[1] = 0.5*k ;

  /*accumulate the dg terms in og for post-multiplication by lambda*/
  quad_subtriangle_mzht(k, y1, y2, rule, s, og) ;
  quad_subtriangle_mzht(k, y2, y3, rule, s, og) ;
  quad_subtriangle_mzht(k, y3, y1, rule, s, og) ;

  gqr_rule_free(rule) ;

  /*multiply dg by lambda*/
  lmr = bem3d_parameters_lambda_real(param) ;
  lmi = bem3d_parameters_lambda_imag(param) ;
  q->free_dg[0] = lmr*og[0] - lmi*og[1] ;
  q->free_dg[1] = lmr*og[1] + lmi*og[0] ;

  alr = bem3d_parameters_coupling_real(param) ;
  ali = bem3d_parameters_coupling_imag(param) ;
  q->free_g[0] = alr*s[0] - ali*s[1] ;
  q->free_g[1] = alr*s[1] + ali*s[0] ;
  
  g_assert(!isnan(q->free_g[0])) ; g_assert(!isnan(q->free_g[1])) ;
  g_assert(!isnan(q->free_dg[0])) ; g_assert(!isnan(q->free_dg[1])) ;

  return 0 ;
}

/**
 * @}
 * 
 */
