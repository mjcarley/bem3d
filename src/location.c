/* location.c
 * 
 * Copyright (C) 2006 Michael Carley
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

#include <glib.h>
#include <gts.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>

#include "bem3d.h"
#include "bem3d-private.h"

/**
 * 
 * @addtogroup belement
 * @{
 *
 */

/** 
 * Find element vertex which is closest to a GtsPoint.
 * 
 * @param e ::BEM3DElement;
 * @param p GtsPoint;
 * @param i local index of closest vertex;
 * @param R distance to closest vertex.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_nearest_vertex(BEM3DElement *e, GtsPoint *p, gint *i, 
				  gdouble *R) 

{
  gint j ;
  gdouble Rv ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(p != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(p), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(i != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(R != NULL, BEM3D_NULL_ARGUMENT) ;

  *R = G_MAXDOUBLE ;
  for ( j = 0 ; j < bem3d_element_vertex_number(e) ; j ++ ) {
    Rv = gts_point_distance2(p, e->v[j]) ;
    if ( Rv < *R ) { *i = j ; *R = Rv ; }
  }
  
  *R = sqrt(*R) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Given the local coordinates of a point on an element, check if it
 * lies inside the element.
 * 
 * @param e ::BEM3DElement
 * @param xi first local coordinate of point to check
 * @param eta second local coordinate of point to check
 * 
 * @return TRUE if (\a xi, \a eta) lies inside the element, FALSE
 * otherwise.
 */

gboolean bem3d_element_point_inside(BEM3DElement *e, gdouble xi, gdouble eta)

{
  gint i, rot[32] ;
  gdouble xi1, eta1, xi2, eta2 ;

  g_return_val_if_fail(e != NULL, FALSE) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), FALSE) ;

  g_debug("%s: e=%p", __FUNCTION__, e) ;

  rotation_indices(bem3d_element_corner_number(e), rot) ;

  for ( i = 0 ; i < bem3d_element_corner_number(e) ; i ++ ) {
    g_debug("%s: r=%d", __FUNCTION__, rot[i]) ;
    xi1 = bem3d_element_vertex_xi(e,bem3d_element_corner_index(e,rot[i])) ;
  g_debug("%s:", __FUNCTION__) ;
    eta1 = bem3d_element_vertex_eta(e,bem3d_element_corner_index(e,rot[i])) ;
    xi2 = bem3d_element_vertex_xi(e,bem3d_element_corner_index(e,rot[i+1])) ;
    eta2 = bem3d_element_vertex_eta(e,bem3d_element_corner_index(e,rot[i+1])) ;
    if ( ((xi2-xi)*(eta1-eta) - (xi1-xi)*(eta2-eta)) < 0 )
      return FALSE ;
  }


  return TRUE ;
}

static gdouble _element_area_coordinates_f(const gsl_vector *y,
					   gpointer q)

{
  gpointer *p = (gpointer *)q ;
  BEM3DElement *e = (BEM3DElement *)p[0] ;
  GtsPoint *x = (GtsPoint *)p[1] ;
  gdouble s, t, L[32], dLds[32], dLdt[32] ;
  gdouble R ;
  GtsPoint x0 ;
  GtsVector r ;

  s = gsl_vector_get(y, 0) ; t = gsl_vector_get(y, 1) ;
  e->shf(s, t, L, dLds, dLdt, NULL) ;
  bem3d_element_position(e, L, &x0) ;
  gts_vector_init(r, &x0, x) ;
  R = gts_vector_norm(r) ;

  return (R*R) ;
}

static void _element_area_coordinates_df(const gsl_vector *y,
					 gpointer q, gsl_vector *df)

{
  gpointer *p = (gpointer *)q ;
  BEM3DElement *e = (BEM3DElement *)p[0] ;
  GtsPoint *x = (GtsPoint *)p[1] ;
  gdouble s, t, L[32], dLds[32], dLdt[32] ;
  gdouble R ;
  GtsPoint x0 ;
  GtsVector r, drds, drdt ;

  s = gsl_vector_get(y, 0) ; t = gsl_vector_get(y, 1) ;
  e->shf(s, t, L, dLds, dLdt, NULL) ;
  bem3d_element_position(e, L, &x0) ;
  gts_vector_init(r, &x0, x) ;
  R = gts_vector_norm(r) ; R*= R ;

  if ( R == 0 ) {
    gsl_vector_set(df, 0, 0) ; gsl_vector_set(df, 1, 0) ;
    return ;
  }

  bem3d_element_slopes(e, dLds, dLdt, drds, drdt) ;
  gsl_vector_set(df, 0, -gts_vector_scalar(r,drds)) ;
  gsl_vector_set(df, 1, -gts_vector_scalar(r,drdt)) ;

  return ; 
}

static void _element_area_coordinates_fdf(const gsl_vector *y,
					  gpointer q,
					  gdouble *f, gsl_vector *df)

{
  gpointer *p = (gpointer *)q ;
  BEM3DElement *e = (BEM3DElement *)p[0] ;
  GtsPoint *x = (GtsPoint *)p[1] ;
  gdouble s, t, L[32], dLds[32], dLdt[32] ;
  gdouble R ;
  GtsPoint x0 ;
  GtsVector r, drds, drdt ;

  s = gsl_vector_get(y, 0) ; t = gsl_vector_get(y, 1) ;

  e->shf(s, t, L, dLds, dLdt, NULL) ;
  bem3d_element_position(e, L, &x0) ;
  gts_vector_init(r, &x0, x) ;
  R = gts_vector_norm(r) ; *f = R*R ;

  if ( R == 0 ) {
    gsl_vector_set(df, 0, 0) ; gsl_vector_set(df, 1, 0) ;
    return ;
  }

  bem3d_element_slopes(e, dLds, dLdt, drds, drdt) ;
  gsl_vector_set(df, 0, -gts_vector_scalar(r,drds)) ;
  gsl_vector_set(df, 1, -gts_vector_scalar(r,drdt)) ;

  return ; 
}

/** 
 * Given a ::BEM3DElement, find the area coordinates of the nearest
 * point to a specified GtsPoint, using a rootfinding method. The
 * point found need not lie on the ::BEM3DElement proper. For a plane
 * element, this is equivalent to projection onto the plane.
 * 
 * @param e ::BEM3DElement;
 * @param x GtsPoint;
 * @param xi first area coordinate of nearest point on ::BEM3DElement;
 * @param eta second area coordinate of nearest point on ::BEM3DElement;
 * @param constrain if TRUE, limit the search to points which lie on the
 * element or its boundary; if FALSE, the returned point may lie outside
 * the element.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_nearest_point(BEM3DElement *e, GtsPoint *x,
				 gdouble *xi, gdouble *eta, 
				 gboolean constrain)

{
  gsl_multimin_function_fdf f ;
  const gsl_multimin_fdfminimizer_type *solver ;
  gsl_multimin_fdfminimizer *s ;
  gsl_vector *init ;
  gpointer p[2] ;
  gint i, iter_max, status ;

  g_debug("%s: x=(%lg,%lg,%lg)", __FUNCTION__, x->x, x->y, x->z) ;

  g_return_val_if_fail(e != NULL, FALSE) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), FALSE) ;
  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(xi != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(eta != NULL, BEM3D_NULL_ARGUMENT) ;

  i = bem3d_element_find_vertex(e, GTS_VERTEX(x)) ;
  if ( i != BEM3D_FAILURE ) {
    *xi = bem3d_element_vertex_xi(e,i) ;
    *eta = bem3d_element_vertex_eta(e,i) ;
    return BEM3D_SUCCESS ;
  }

  iter_max = 16 ;
  f.f = &_element_area_coordinates_f ;
  f.df = &_element_area_coordinates_df ;
  f.fdf = &_element_area_coordinates_fdf ;
  f.n = 2 ;
  f.params = p ;
  p[0] = e ; p[1] = x ;

  init = gsl_vector_alloc(2) ;
  gsl_vector_set(init, 0, *xi) ; gsl_vector_set(init, 1, *eta) ;
  solver = gsl_multimin_fdfminimizer_conjugate_fr ;
  s = gsl_multimin_fdfminimizer_alloc(solver, 2) ;
  gsl_multimin_fdfminimizer_set(s, &f, init, 1e-1, 1e-2) ;
  
  for ( (i = 0), (status = GSL_CONTINUE) ;
	( i < 2 ) ||
	  ((i < iter_max ) && ( status == GSL_CONTINUE)) ; i ++ ) {
    status = gsl_multimin_fdfminimizer_iterate(s) ;
    if ( status ) break ;
    status = gsl_multimin_test_gradient(s->gradient, 1e-12) ;
    if ( status != GSL_CONTINUE ) break ;
    status = gsl_multimin_test_size(s->f, 1e-12) ;
    if ( status != GSL_CONTINUE ) break ;
  }

  *xi = gsl_vector_get(s->x, 0) ; *eta = gsl_vector_get(s->x, 1) ;

  if ( constrain && ! bem3d_element_point_inside(e, *xi, *eta) ) {
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: point (%g,%g) lies outside element, constraining",
	  __FUNCTION__, *xi, *eta) ;
/*     fprintf(stderr, "Entering constraint\n") ; */
    bem3d_element_boundary_nearest_point(e, x, xi, eta) ;
/*     fprintf(stderr, "Exiting constraint\n") ; */
  }

  g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	"%s: %d iterations; (s,t)=(%g,%g); f: %lg; df: (%lg,%lg); %s",
	__FUNCTION__, i, *xi, *eta,
	s->f,
	gsl_vector_get(s->gradient, 0),
	gsl_vector_get(s->gradient, 1),
	gsl_strerror(status)) ;

  gsl_multimin_fdfminimizer_free(s) ;
  gsl_vector_free(init) ;

  g_debug("%s:", __FUNCTION__) ;

  if ( i == iter_max ) return BEM3D_ITERMAX ;

  return BEM3D_SUCCESS ;
}


static gdouble edge_nearest_point_parameter(GtsPoint *x1, GtsPoint *x2,
					    GtsPoint *x)

{
  gdouble t ;
  GtsVector n, s ;

  gts_vector_init(n, x1, x2) ;
  gts_vector_init(s, x1, x) ;

  t = gts_vector_scalar(s,n)/gts_vector_scalar(n,n) ;

  return t ;
}

static gdouble edge_nearest(gdouble t, gpointer p)

{
  gpointer *q = (gpointer *)p ;
  BEM3DElement *e = BEM3D_ELEMENT(q[0]) ;
  GtsPoint *x = GTS_POINT(q[1]) ;
  gint i = *((gint *)q[2]) ;
  gint j = *((gint *)q[3]) ;
  GtsPoint *y = GTS_POINT(q[4]) ;
  gdouble L[32], xi, eta ;
  BEM3DShapeFunc shfunc ;

  shfunc = bem3d_element_shape_func(e) ;

  xi = (1-t)*bem3d_element_vertex_xi(e,i) + t*bem3d_element_vertex_xi(e,j) ;
  eta = (1-t)*bem3d_element_vertex_eta(e,i) + t*bem3d_element_vertex_eta(e,j) ;
  shfunc(xi, eta, L, NULL, NULL, NULL) ;
  bem3d_element_position(e, L, y) ;

  return (gts_point_distance(x,y)) ;
}

/** 
 * Find the area coordinates of the point on the edge of an element
 * which is closest to a point x.
 * 
 * @param e ::BEM3DElement
 * @param i index of first corner of edge
 * @param j index of second corner of edge
 * @param x  ::GtsPoint 
 * @param xi first area coordinate of point on edge \a c0--\a c1 which 
 * is closest to \a x
 * @param eta second area coordinate of point on edge \a i--\a j which 
 * is closest to \a x
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_edge_nearest_point(BEM3DElement *e, gint i, gint j,
				      GtsPoint *x,
				      gdouble *xi, gdouble *eta)

{
  gint status;
  gint k = 0, kmax = 8 ;
  static const gsl_min_fminimizer_type *m;
  static gsl_min_fminimizer *s;
  gdouble a, b, t ;
  gsl_function f ;
  static GtsPoint *y = NULL ;
  gpointer p[8] ;

  g_debug("%s:", __FUNCTION__) ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(xi != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(eta != NULL, BEM3D_NULL_ARGUMENT) ;

  f.function = &edge_nearest ;
  f.params = p ;
  if ( y == NULL ) {
    m = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc(m);
    y = gts_point_new(gts_point_class(), 0, 0, 0) ;
  }
  p[0] = e ; p[1] = x ; p[2] = &i ; p[3] = &j ; p[4] = y ;
  
  t = edge_nearest_point_parameter(bem3d_element_vertex(e,i),
				   bem3d_element_vertex(e,j),
				   x) ;

  gsl_min_fminimizer_set(s, &f, t, t-1.0, t+1.0) ;

  k = 0 ;
  do {
    k ++ ;
    status = gsl_min_fminimizer_iterate(s) ;
    t = gsl_min_fminimizer_x_minimum(s) ;
    a = gsl_min_fminimizer_x_lower(s) ;
    b = gsl_min_fminimizer_x_upper(s) ;

    status = gsl_min_test_interval(a, b, 1e-6, 0.0);

  } while ( (k < kmax) && (status == GSL_CONTINUE) ) ;

  if ( t < 0 ) t = 0.0 ; if ( t > 1 ) t = 1.0 ;
  *xi = (1-t)*bem3d_element_vertex_xi(e,i) + t*bem3d_element_vertex_xi(e,j) ;
  *eta = (1-t)*bem3d_element_vertex_eta(e,i) + t*bem3d_element_vertex_eta(e,j) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Find the area coordinates of the point on any edge of an element
 * which is closest to a point x.
 * 
 * @param e ::BEM3DElement
 * @param x  ::GtsPoint 
 * @param xi first area coordinate of point on edge \a c0--\a c1 which 
 * is closest to \a x
 * @param eta second area coordinate of point on edge \a i--\a j which 
 * is closest to \a x
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_element_boundary_nearest_point(BEM3DElement *e, GtsPoint *x,
					  gdouble *xi, gdouble *eta)

{
  gint i, rot[32] ;
  gdouble L[32], s, t ;
  gdouble Rmin, R ;
  static GtsPoint *y  = NULL ;
  BEM3DShapeFunc shfunc ;

  g_debug("%s:", __FUNCTION__) ;

  g_return_val_if_fail(e != NULL, FALSE) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), FALSE) ;
  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_POINT(x), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(xi != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(eta != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( y == NULL ) {
    y = gts_point_new(gts_point_class(), 0, 0, 0) ;
  }

  shfunc = bem3d_element_shape_func(e) ;
  rotation_indices(bem3d_element_corner_number(e), rot) ;

  Rmin = G_MAXDOUBLE ;
  for ( i = 0 ; i < bem3d_element_corner_number(e) ; i ++ ) {
    bem3d_element_edge_nearest_point(e, 
				   bem3d_element_corner_index(e,rot[i]),
				   bem3d_element_corner_index(e,rot[i+1]),
				   x, &s, &t) ;
    shfunc(s, t, L, NULL, NULL, NULL) ;
    bem3d_element_position(e, L, y) ;
    if ( (R = gts_point_distance(x, y)) < Rmin ) {
      Rmin = R ; *xi = s ; *eta = t ;
    }
  }

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
