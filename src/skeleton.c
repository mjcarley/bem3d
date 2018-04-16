/* skeleton.c
 * 
 * Copyright (C) 2017 Michael Carley
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
 * @defgroup skeleton Mesh skeletons
 *
 * BEM3D skeletons are representations of meshes as arrays of point
 * sources, mainly intended for use in interfacing to third-party
 * libraries such as Fast Multipole Methods which take input as a set
 * of elementary sources.
 *
 * @{
 * 
 */

/** 
 * Allocate a new ::BEM3DMeshSkeleton for a given mesh
 * 
 * @param m a ::BEM3DMesh
 * @param order_max maximum number of source points per element on 
 * the skeleton 
 * 
 * @return a new ::BEM3DMeshSkeleton on success.
 */

BEM3DMeshSkeleton *bem3d_mesh_skeleton_new(BEM3DMesh *m, gint order_max)

/*
  packing of data inside s:

  geometry: sources followed by targets, 3 entries per block
  x[3*i]: x_i y_i z_i 


  n[3*i]: nx_i ny_i nz_i

  weights and indices for interpolation:

  idx: ppe entries per block, one per source point

  w: ppe entries per block, one per source point
*/

{
  BEM3DMeshSkeleton *s ;

  s = (BEM3DMeshSkeleton *)g_malloc(sizeof(BEM3DMeshSkeleton)) ;

  s->m = m ;

  /*number of target (collocation) points*/
  s->nnodes = bem3d_mesh_node_number(m) ;
  /*number of mesh elements*/
  s->nelem = bem3d_mesh_element_number(m) ;

  s->order = s->ns = s->nt = 0 ;
  s->npts = s->nelem*order_max + s->nnodes ;
  s->ppe = bem3d_mesh_element_node_number_max(m) ;

  s->x = (gdouble *)g_malloc(((s->npts)*6 + 
			      s->nelem*order_max*s->ppe)*sizeof(gdouble)) ;
  s->n = &(s->x[3*(s->npts)]) ;
  s->w = &(s->n[3*(s->npts)]) ;

  s->idx = (gint *)g_malloc(s->nelem*order_max*s->ppe*sizeof(gint)) ;

  bem3d_mesh_index_range(m, &(s->imin), &(s->imax)) ;
  
  s->e = g_hash_table_new(NULL, NULL) ;

  if ( s->imax - s->imin + 1 != s->nnodes ) 
    g_error("%s: indices must be contiguous on mesh", 
	    __FUNCTION__) ;

  return s ;
}

static gint source_points(BEM3DElement *e, gpointer data[])

{
  gint *ne = data[0] ;
  BEM3DQuadratureRule *q = data[1] ;
  BEM3DMeshSkeleton *skel = data[2] ;
  BEM3DShapeFunc shfunc = bem3d_element_shape_func(e) ;
  BEM3DShapeFunc cpfunc = bem3d_element_node_func(e) ;
  gdouble L[32], dLds[32], dLdt[32], s, t, w, J ;
  gint i, j ;
  GtsPoint y ;
  GtsVector n ;

  g_hash_table_insert(skel->e, (gpointer)e, GINT_TO_POINTER(skel->ns)) ;
  for ( i = 0 ; i < bem3d_quadrature_vertex_number(q) ; i ++ ) {
    s = bem3d_quadrature_xi(q,i) ;
    t = bem3d_quadrature_eta(q,i) ;
    w = bem3d_quadrature_weight(q,i) ;
    shfunc(s, t, L, dLds, dLdt, NULL) ;
    bem3d_element_position(e, L, &y) ;
    bem3d_element_normal(e, dLds, dLdt, n, &J) ;

    skel->x[3*(skel->ns)+0] = y.x ;
    skel->x[3*(skel->ns)+1] = y.y ;
    skel->x[3*(skel->ns)+2] = y.z ;

    skel->n[3*(skel->ns)+0] = n[0] ;
    skel->n[3*(skel->ns)+1] = n[1] ;
    skel->n[3*(skel->ns)+2] = n[2] ;

    cpfunc(s, t, L, NULL, NULL, NULL) ;
    
    for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
      /* g_assert(w*L[j]*J >= 0.0) ; */
      skel->w[(skel->ppe)*(skel->ns)+j] = w*L[j]*J ;
      skel->idx[(skel->ppe)*(skel->ns)+j] = 
	bem3d_element_global_index(e, j) - skel->imin ;
    }

    for ( j = bem3d_element_node_number(e) ; j < skel->ppe ; j ++ ) {
      skel->w[(skel->ppe)*(skel->ns)+j] = 0.0 ;
      skel->idx[(skel->ppe)*(skel->ns)+j] = 0 ;
    }

    skel->ns ++ ;
  }

  return 0 ;
}

static gint target_points(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMeshSkeleton *skel = data[2] ;
  gint off, j ;
  GtsVector n ;

  g_assert(i <= skel->imax) ;

  off = 3*(skel->ns) ; j = i - skel->imin ;

  bem3d_node_normal(skel->m, i, n, skel->anorm) ;

  skel->x[off+3*j+0] = GTS_POINT(v)->x ;
  skel->x[off+3*j+1] = GTS_POINT(v)->y ;
  skel->x[off+3*j+2] = GTS_POINT(v)->z ;

  skel->n[off+3*j+0] = n[0] ;
  skel->n[off+3*j+1] = n[1] ;
  skel->n[off+3*j+2] = n[2] ;

  skel->nt ++ ;

  return 0 ;
}

/** 
 * Initialize a ::BEM3DMeshSkeleton which has been allocated with
 * ::bem3d_mesh_skeleton_new. You should be careful to use a
 * ::BEM3DQuadratureRule which guarantees source points which are
 * distinct from the collocation points. For example, a three-point
 * Gaussian quadrature with nodes at the midpoints of triangle edges
 * will generate source points which coincide with the mid-edge nodes
 * of a second-order triangular element.
 * 
 * @param s a ::BEM3DMeshSkeleton
 * @param q a ::BEM3DQuadratureRule used to find points on mesh elements
 * @param anorm a ::BEM3DAverage specification used in computing collocation
 * point normals
 * 
 * @return BEM3D_SUCCESS on success
 */

gint bem3d_mesh_skeleton_init(BEM3DMeshSkeleton *s, BEM3DQuadratureRule *q,
			      BEM3DAverage anorm)

{
  gint ne ;
  gpointer data[] = {&ne, q, s} ;

  g_assert(bem3d_mesh_element_number(s->m)*bem3d_quadrature_vertex_number(q)+
	   bem3d_mesh_node_number(s->m) <= s->npts) ;

  /*number of source and target points, will be incremented on the traverse*/
  s->ns = s->nt = 0 ; 

  s->order = bem3d_quadrature_vertex_number(q) ;
  s->anorm = anorm ;

  ne = 0 ;

  bem3d_mesh_foreach_element(s->m, (BEM3DElementFunc)source_points, data) ;

  bem3d_mesh_foreach_node(s->m, (BEM3DNodeFunc)target_points, data) ;

  return 0 ;
}

/** 
 * Write a ::BEM3DMeshSkeleton to file
 * 
 * @param s a ::BEM3DMeshSkeleton
 * @param f FILE pointer
 * 
 * @return BEM3D_SUCCESS on success
 */

gint bem3d_mesh_skeleton_write(BEM3DMeshSkeleton *s, FILE *f)

{
  gint i, j, *idx ;
  gdouble *w ;

  fprintf(f, "%d %d %d %d %d %d %d %d %d\n",
	  s->nnodes, s->nelem, s->npts, s->order,
	  s->ns, s->nt, s->imin,
	  s->imax, s->ppe) ;

  for ( i = 0 ; i < s->ns+s->nt ; i ++ ) {
    fprintf(f, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",
	    s->x[3*i+0], s->x[3*i+1], s->x[3*i+2],
	    s->n[3*i+0], s->n[3*i+1], s->n[3*i+2]) ;
  }

  for ( i = 0 ; i < s->ns ; i ++ ) {
    w = &(s->w[i*s->ppe]) ; idx = &(s->idx[i*s->ppe]) ;
    for ( j = 0 ; j < s->ppe ; j ++ ) {
      fprintf(f, "%d %1.16e ", idx[j], w[j]) ;
    }
    fprintf(f, "\n") ;
  }

  return 0 ;
}

gint bem3d_mesh_skeleton_read(BEM3DMeshSkeleton *s, FILE *f)

{
  gint i, j, *idx ;
  gdouble *w ;

  g_assert_not_reached() ; /*untested, uncompleted code*/

  fprintf(f, "%d %d %d %d %d %d %d %d %d\n",
	  s->nnodes, s->nelem, s->npts, s->order,
	  s->ns, s->nt, s->imin,
	  s->imax, s->ppe) ;

  for ( i = 0 ; i < s->ns+s->nt ; i ++ ) {
    fprintf(f, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",
	    s->x[3*i+0], s->x[3*i+1], s->x[3*i+2],
	    s->n[3*i+0], s->n[3*i+1], s->n[3*i+2]) ;
  }

  for ( i = 0 ; i < s->ns ; i ++ ) {
    w = &(s->w[i*s->ppe]) ; idx = &(s->idx[i*s->ppe]) ;
    for ( j = 0 ; j < s->ppe ; j ++ ) {
      fprintf(f, "%d %1.16e ", idx[j], w[j]) ;
    }
    fprintf(f, "\n") ;
  }

  return 0 ;
}

/**
 * @}
 * 
 */
