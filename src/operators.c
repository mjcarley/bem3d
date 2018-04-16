/* operators.c
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

/**
 * @defgroup operators Differential (and other) operators
 *
 * BEM3D contains a number of functions for computing differential
 * operators applied to surface quantities. These are returned in a
 * ::BEM3DOperator which can be applied to surface data using the
 * appropriate functions. For example:
 * @code
 * BEM3DOperator *op ;
 * GtsVector df ;
 * gdouble *opi ;
 *
 * op = bem3d_operator_new() ;
 * bem3d_operator_gradient(m, i, op, BEM3D_AVERAGE_MWSELR) ;
 * df[0] = df[1] = df[2] = 0.0 ;
 * for ( j = 0 ; j < bem3d_operator_length(op) ; j ++ ) {
 *    opi = bem3d_operator_weight(op,j) ;
 *    df[0] += f[bem3d_operator_index(op,j)]*opi[0] ;
 *    df[1] += f[bem3d_operator_index(op,j)]*opi[1] ;
 *    df[2] += f[bem3d_operator_index(op,j)]*opi[2] ;
 * }
 *
 * @endcode
 * @{
 * 
 */

/** 
 * Allocate a new ::BEM3DOperator.
 * 
 * 
 * @return a newly allocated ::BEM3DOperator.
 */

BEM3DOperator *bem3d_operator_new(void)

{
  BEM3DOperator *op ;

  op = (BEM3DOperator *)g_malloc(sizeof(BEM3DOperator)) ;

  op->nc = 0 ;
  op->w = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
  op->id = g_array_new(FALSE, TRUE, sizeof(int)) ;

  return op ;
}

static void insert_element_nodes(BEM3DElement *e, GArray *id)

{
  gint i, j ;
  gboolean insert ;

  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    for ( (j = 0), (insert = TRUE) ; j < id->len ; j ++ ) {
      if ( bem3d_element_global_index(e,i) == g_array_index(id,gint,j) ) {
	insert = FALSE ; break ;
      }
    }
    if ( insert ) g_array_append_val(id, bem3d_element_global_index(e,i)) ;
  }

  return ;
}

static void edge_lengths_angle(BEM3DElement *e, gint i,
			       gdouble *E1, gdouble *E2, gdouble *S)

{
  GtsPoint *v1, *v2, *v3 ;
  gint j1, j2, j3, nc ;
  GtsVector e1, e2, es ;

  nc = bem3d_element_corner_number(e) ;
  j1 = i-1 ; j2 = i ; j3 = i+1 ;
  if ( i == 0 ) { j1 = nc-1 ; j2 = i ; j3 = i+1 ; } 
  if ( i == nc-1 ) { j1 = nc-2 ; j2 = i ; j3 = 0 ; }

  v1 = GTS_POINT(bem3d_element_vertex(e, bem3d_element_corner_index(e,j1))) ;
  v2 = GTS_POINT(bem3d_element_vertex(e, bem3d_element_corner_index(e,j2))) ;
  v3 = GTS_POINT(bem3d_element_vertex(e, bem3d_element_corner_index(e,j3))) ;

  gts_vector_init(e1, v2, v1) ; gts_vector_init(e2, v2, v3) ;
  gts_vector_cross(es, e1, e2) ;
  
  *E1 = gts_vector_norm(e1) ; *E2 = gts_vector_norm(e2) ;
  *S = gts_vector_norm(es)/(*E1)/(*E2) ;

  if ( *S > 1.0 ) *S = 1.0 ;
  if ( *S < -1.0 ) *S = -1.0 ;

  return ;
}

static void average_normal_mwa(BEM3DElement *e, GtsVertex *v, 
			       GtsVector n) 

{
  gint i, j ;
  gdouble L[32], dLds[32], dLdt[32], J, E1, E2, S ;
  BEM3DShapeFunc shf ;
  GtsVector ne ;

  if ( (i = bem3d_element_find_vertex(e, v)) == -1 ) 
    g_error("%s: vertex %p is not on element %p", __FUNCTION__, v, e) ;

  shf = bem3d_element_shape_func(e) ;
  shf(bem3d_element_vertex_xi(e,i), bem3d_element_vertex_eta(e,i), 
      L, dLds, dLdt, NULL) ;

  bem3d_element_normal(e, dLds, dLdt, ne, &J) ;
  gts_vector_normalize(ne) ;

  if ( (j = bem3d_element_vertex_is_corner(e, v)) == -1 ) {
    n[0] += ne[0] ; n[1] += ne[1] ; n[2] += ne[2] ;
    return ;
  }

  edge_lengths_angle(e, j, &E1, &E2, &S) ;

  S = asin(S) ;

  n[0] += S*ne[0] ; n[1] += S*ne[1] ; n[2] += S*ne[2] ;

  return ;
}

static void average_normal_mwselr(BEM3DElement *e, GtsVertex *v, 
				  GtsVector n) 

{
  gint i, j ;
  gdouble L[32], dLds[32], dLdt[32], J, E1, E2, S ;
  BEM3DShapeFunc shf ;
  GtsVector ne ;

  if ( (i = bem3d_element_find_node(e, v)) == -1 ) 
    g_error("%s: node %p is not on element %p", __FUNCTION__, v, e) ;

  shf = bem3d_element_shape_func(e) ;
  shf(bem3d_element_node_xi(e,i), bem3d_element_node_eta(e,i), 
      L, dLds, dLdt, NULL) ;

  bem3d_element_normal(e, dLds, dLdt, ne, &J) ;
  gts_vector_normalize(ne) ;

  if ( (j = bem3d_element_vertex_is_corner(e, v)) == -1 ) {
    n[0] += ne[0] ; n[1] += ne[1] ; n[2] += ne[2] ;
    return ;
  }

  edge_lengths_angle(e, j, &E1, &E2, &S) ;

  n[0] += S*ne[0]/E1/E2 ; n[1] += S*ne[1]/E1/E2 ; n[2] += S*ne[2]/E1/E2 ;

  return ;
}

static void average_normal_mwaat(BEM3DElement *e, GtsVertex *v, 
				  GtsVector n) 

{
  gint i, j ;
  gdouble L[32], dLds[32], dLdt[32], J, E1, E2, S ;
  BEM3DShapeFunc shf ;
  GtsVector ne ;

  if ( (i = bem3d_element_find_vertex(e, v)) == -1 ) 
    g_error("%s: vertex %p is not on element %p", __FUNCTION__, v, e) ;

  shf = bem3d_element_shape_func(e) ;
  shf(bem3d_element_vertex_xi(e,i), bem3d_element_vertex_eta(e,i), 
      L, dLds, dLdt, NULL) ;

  bem3d_element_normal(e, dLds, dLdt, ne, &J) ;
  gts_vector_normalize(ne) ;

  if ( (j = bem3d_element_vertex_is_corner(e, v)) == -1 ) {
    n[0] += ne[0] ; n[1] += ne[1] ; n[2] += ne[2] ;
    return ;
  }

  edge_lengths_angle(e, j, &E1, &E2, &S) ;

  n[0] += S*ne[0]*E1*E2 ; n[1] += S*ne[1]*E1*E2 ; n[2] += S*ne[2]*E1*E2 ;

  return ;
}

static void average_normal_mwelr(BEM3DElement *e, GtsVertex *v, 
				  GtsVector n) 

{
  gint i, j ;
  gdouble L[32], dLds[32], dLdt[32], J, E1, E2, S ;
  BEM3DShapeFunc shf ;
  GtsVector ne ;

  if ( (i = bem3d_element_find_vertex(e, v)) == -1 ) 
    g_error("%s: vertex %p is not on element %p", __FUNCTION__, v, e) ;

  shf = bem3d_element_shape_func(e) ;
  shf(bem3d_element_vertex_xi(e,i), bem3d_element_vertex_eta(e,i), 
      L, dLds, dLdt, NULL) ;

  bem3d_element_normal(e, dLds, dLdt, ne, &J) ;
  gts_vector_normalize(ne) ;

  if ( (j = bem3d_element_vertex_is_corner(e, v)) == -1 ) {
    n[0] += ne[0] ; n[1] += ne[1] ; n[2] += ne[2] ;
    return ;
  }

  edge_lengths_angle(e, j, &E1, &E2, &S) ;

  n[0] += ne[0]/E1/E2 ; n[1] += ne[1]/E1/E2 ; n[2] += ne[2]/E1/E2 ;

  return ;
}

static void average_normal_mwrelr(BEM3DElement *e, GtsVertex *v, 
				  GtsVector n) 

{
  gint i, j ;
  gdouble L[32], dLds[32], dLdt[32], J, E1, E2, S ;
  BEM3DShapeFunc shf ;
  GtsVector ne ;

  if ( (i = bem3d_element_find_vertex(e, v)) == -1 ) 
    g_error("%s: vertex %p is not on element %p", __FUNCTION__, v, e) ;

  shf = bem3d_element_shape_func(e) ;
  shf(bem3d_element_vertex_xi(e,i), bem3d_element_vertex_eta(e,i), 
      L, dLds, dLdt, NULL) ;

  bem3d_element_normal(e, dLds, dLdt, ne, &J) ;
  gts_vector_normalize(ne) ;

  if ( (j = bem3d_element_vertex_is_corner(e, v)) == -1 ) {
    n[0] += ne[0] ; n[1] += ne[1] ; n[2] += ne[2] ;
    return ;
  }

  edge_lengths_angle(e, j, &E1, &E2, &S) ;

  n[0] += ne[0]/sqrt(E1*E2) ; n[1] += ne[1]/sqrt(E1*E2) ; 
  n[2] += ne[2]/sqrt(E1*E2) ;

  return ;
}

static void average_normal_mwe(BEM3DElement *e, GtsVertex *v, 
			       GtsVector n) 

{
  gint i ;
  gdouble L[32], dLds[32], dLdt[32], J ;
  BEM3DShapeFunc shf ;
  GtsVector ne ;

  g_assert(bem3d_element_shape_func(e) == bem3d_shfunc_t1 ) ;

  if ( (i = bem3d_element_find_vertex(e, v)) == -1 ) 
    g_error("%s: vertex %p is not on element %p", __FUNCTION__, v, e) ;

  shf = bem3d_element_shape_func(e) ;
  shf(bem3d_element_vertex_xi(e,i), bem3d_element_vertex_eta(e,i), 
      L, dLds, dLdt, NULL) ;

  bem3d_element_normal(e, dLds, dLdt, ne, &J) ;

  n[0] += ne[0] ; n[1] += ne[1] ; n[2] += ne[2] ;

  return ;
}

/** 
 * Compute the surface normal at a node as the weighted average, where
 * appropriate, of the normals on each of the elements of a mesh using
 * the node, using the methods in Jin, S., Lewis, R. and West, D., `A
 * comparison of algorithms for vertex normal computation', The Visual
 * Computer, 21(1-2):71--82,
 * http://dx.doi.org/10.1007/s00371-004-0271-1
 * 
 * @param m a ::BEM3DMesh;
 * @param i index of a node of \a m;
 * @param n the normal at node \a i;
 * @param mode a ::BEM3DAverage defining the averaging method to be used. 
 * 
 * @return 
 */

gint bem3d_node_normal(BEM3DMesh *m, gint i, GtsVector n, 
		       BEM3DAverage mode)

{
  GSList *e, *j ;
  GtsVertex *v ;

  g_return_val_if_fail(m != NULL, BEM3D_EINVAL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_EINVAL) ;
  g_return_val_if_fail(i >= 0, BEM3D_EINVAL) ;  
  g_return_val_if_fail(n != NULL, BEM3D_EINVAL) ;

  e = bem3d_mesh_node_elements(m, i) ;
  v = bem3d_mesh_node_from_index(m,i) ;

  n[0] = n[1] = n[2] = 0.0 ;

  switch ( mode ) {
  default:
    g_error("%s: unrecognized mode %d", __FUNCTION__, mode) ;
    break ;
  case BEM3D_AVERAGE_MWE:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_normal_mwe(BEM3D_ELEMENT(j->data), v, n) ;
    break ;
  case BEM3D_AVERAGE_MWA:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_normal_mwa(BEM3D_ELEMENT(j->data), v, n) ;
    break ;
  case BEM3D_AVERAGE_MWSELR:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_normal_mwselr(BEM3D_ELEMENT(j->data), v, n) ;
    break ;
  case BEM3D_AVERAGE_MWAAT:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_normal_mwaat(BEM3D_ELEMENT(j->data), v, n) ;
    break ;
  case BEM3D_AVERAGE_MWELR:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_normal_mwelr(BEM3D_ELEMENT(j->data), v, n) ;
    break ;
  case BEM3D_AVERAGE_MWRELR:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_normal_mwrelr(BEM3D_ELEMENT(j->data), v, n) ;
    break ;
  }

  gts_vector_normalize(n) ;

  return BEM3D_SUCCESS ;
}

static void vector_gamma(GtsVertex *p, GtsVertex *q, GtsVertex *r,
			 GtsVector v)

{
  GtsVector qp, rq, pr ;
  gdouble dp ;

  gts_vector_init(qp, GTS_POINT(q), GTS_POINT(p)) ;
  gts_vector_init(rq, GTS_POINT(r), GTS_POINT(q)) ;
  gts_vector_init(pr, GTS_POINT(p), GTS_POINT(r)) ;
  dp = gts_vector_scalar(qp,rq) ;

  v[0] = dp*pr[0] ; v[1] = dp*pr[1] ; v[2] = dp*pr[2] ;

  return ;
}

static gint search_int_array(gint *i, gint n, gint i0)

{
  gint j ;

  for ( j = 0 ; (j < n) && (i[j] != i0) ; j ++ ) ;
  if ( i[j] != i0 ) g_error("%s: indexed value is not correct i[%d] != %d",
			    __FUNCTION__, j, i0) ;

  return j ;
}

static void element_gradient(BEM3DElement *e, gdouble wt,
			     BEM3DOperator *op)

{
  gint j ;
  GtsVector df0, df1 ;
  gdouble *w ;

  j = search_int_array((gint *)(op->id->data), op->id->len,
		       bem3d_element_global_index(e,0)) ;
  vector_gamma(bem3d_element_node(e,0),bem3d_element_node(e,1),
	       bem3d_element_node(e,2), df0) ;
  vector_gamma(bem3d_element_node(e,0),bem3d_element_node(e,2),
	       bem3d_element_node(e,1), df1) ;
  w = bem3d_operator_weight(op,j) ;
  w[0] += wt*(df0[0]+df1[0]) ; w[1] += wt*(df0[1]+df1[1]) ;
  w[2] += wt*(df0[2]+df1[2]) ;

  j = search_int_array((gint *)(op->id->data), op->id->len,
		       bem3d_element_global_index(e,1)) ;
  vector_gamma(bem3d_element_node(e,1),bem3d_element_node(e,0),
	       bem3d_element_node(e,2), df0) ;
  vector_gamma(bem3d_element_node(e,1),bem3d_element_node(e,2),
	       bem3d_element_node(e,0), df1) ;
  w = bem3d_operator_weight(op,j) ;
  w[0] += wt*(df0[0]+df1[0]) ; w[1] += wt*(df0[1]+df1[1]) ;
  w[2] += wt*(df0[2]+df1[2]) ;
  
  j = search_int_array((gint *)(op->id->data), op->id->len,
		       bem3d_element_global_index(e,2)) ;
  vector_gamma(bem3d_element_node(e,2),bem3d_element_node(e,1),
	       bem3d_element_node(e,0), df0) ;
  vector_gamma(bem3d_element_node(e,2),bem3d_element_node(e,0),
	       bem3d_element_node(e,1), df1) ;
  w = bem3d_operator_weight(op,j) ;
  w[0] += wt*(df0[0]+df1[0]) ; w[1] += wt*(df0[1]+df1[1]) ;
  w[2] += wt*(df0[2]+df1[2]) ;  

  return ;
}

static void average_gradient_mwe(BEM3DElement *e, GtsVertex *v,
				 BEM3DOperator *op, gdouble *tw)

{
  gint i ;
  gdouble A, wt ;

  g_assert(bem3d_element_shape_func(e) == bem3d_shfunc_t1 ) ;

  if ( (i = bem3d_element_find_vertex(e, v)) == -1 ) 
    g_error("%s: vertex %p is not on element %p", __FUNCTION__, v, e) ;

  A = gts_triangle_area(bem3d_element_face(e,0)) ;
  wt = 1.0/(4.0*A*A) ;

  element_gradient(e, wt, op) ; (*tw) += 1.0 ;

  return ;
}

static void average_gradient_mwa(BEM3DElement *e, GtsVertex *v,
				 BEM3DOperator *op, gdouble *tw)

{
  gint i, j ;
  gdouble A, S, E1, E2, wt ;

  g_assert(bem3d_element_shape_func(e) == bem3d_shfunc_t1 ) ;

  if ( (i = bem3d_element_find_vertex(e, v)) == -1 ) 
    g_error("%s: vertex %p is not on element %p", __FUNCTION__, v, e) ;

  j = bem3d_element_vertex_is_corner(e, v) ;
  edge_lengths_angle(e, j, &E1, &E2, &S) ;

  wt = asin(S) ; (*tw) += wt ;

  A = gts_triangle_area(bem3d_element_face(e,0)) ;
  wt /= 4.0*A*A ; 

  element_gradient(e, wt, op) ;

  return ;
}

static void average_gradient_mwselr(BEM3DElement *e, GtsVertex *v,
				    BEM3DOperator *op, gdouble *tw)

{
  gint i, j ;
  gdouble A, S, E1, E2, wt ;

  g_assert(bem3d_element_shape_func(e) == bem3d_shfunc_t1 ) ;

  if ( (i = bem3d_element_find_vertex(e, v)) == -1 ) 
    g_error("%s: vertex %p is not on element %p", __FUNCTION__, v, e) ;

  j = bem3d_element_vertex_is_corner(e, v) ;
  edge_lengths_angle(e, j, &E1, &E2, &S) ;

  wt = S/E1/E2 ; (*tw) += wt ;

  A = gts_triangle_area(bem3d_element_face(e,0)) ;
  wt /= 4.0*A*A ; 

  element_gradient(e, wt, op) ;

  return ;
}

/** 
 * Compute the operator for computation of the surface gradient of a
 * quantity at a node of a BEM3DMesh. Currently, the function is only
 * implemented for linear triangles, using the gradient discretization
 * of Xu, G., `Convergent discrete Laplace-Beltrami operators over
 * triangular surfaces', Proceedings of the Geometric Modeling and
 * Processing 2004 (GMP 04),
 * http://dx.doi.org/10.1109/GMAP.2004.1290041
 * 
 * @param m a ::BEM3DMesh;
 * @param i index of a node of \a m;
 * @param op on exit, contains the operator for \f$\nabla_{T}f\f$;
 * @param mode a ::BEM3DAverage specifying the mode for averaging over
 * elements (not all of these are implemented).
 * 
 * @return BEM3D_SUCCESS on success, an error code otherwise.
 */

gint bem3d_operator_gradient(BEM3DMesh *m, gint i, BEM3DOperator *op,
			     BEM3DAverage mode)

{
  GSList *e, *j ;
  GtsVertex *v ;
  gdouble tw ;
  gint k ;

  g_return_val_if_fail(m != NULL, BEM3D_EINVAL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_EINVAL) ;
  g_return_val_if_fail(i >= 0, BEM3D_EINVAL) ;  
  g_return_val_if_fail(op != NULL, BEM3D_EINVAL) ;

  e = bem3d_mesh_node_elements(m, i) ;
  v = bem3d_mesh_node_from_index(m,i) ;

  g_array_set_size(op->id, 0) ; g_array_set_size(op->w, 0) ;
  for ( j = e ; j != NULL ; j = j->next ) 
    insert_element_nodes(BEM3D_ELEMENT(j->data), op->id) ;

  g_array_set_size(op->w, 3*(op->id->len)) ;
  op->nc = 3 ;

  tw = 0.0 ;
  switch ( mode ) {
  default:
    g_error("%s: unrecognized mode %d", __FUNCTION__, mode) ;
    break ;
  case BEM3D_AVERAGE_MWE:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_gradient_mwe(BEM3D_ELEMENT(j->data), v, op, &tw) ;
    break ;
  case BEM3D_AVERAGE_MWA:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_gradient_mwa(BEM3D_ELEMENT(j->data), v, op, &tw) ;
    break ;
  case BEM3D_AVERAGE_MWSELR:
    for ( j = e ; j != NULL ; j = j->next ) 
      average_gradient_mwselr(BEM3D_ELEMENT(j->data), v, op, &tw) ;
    break ;
  case BEM3D_AVERAGE_MWAAT:
    g_error("%s: mode %d not implemented", __FUNCTION__, mode) ;
    break ;
  case BEM3D_AVERAGE_MWELR:
    g_error("%s: mode %d not implemented", __FUNCTION__, mode) ;
    break ;
  case BEM3D_AVERAGE_MWRELR:
    g_error("%s: mode %d not implemented", __FUNCTION__, mode) ;
    break ;
  }

  for ( k = 0 ; k < 3*(op->id->len) ; k ++ )
    g_array_index(op->w,gdouble,k) /= tw ;

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
