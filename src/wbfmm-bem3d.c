/* wbfmm-bem3d.c
 * 
 * Copyright (C) 2019 Michael Carley
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

#ifdef HAVE_WBFMM

#include <wbfmm.h>
#include "wbfmm-bem3d.h"

#else /*HAVE_WBFMM*/
#endif /*HAVE_WBFMM*/

gint _bem3d_fmm_helmholtz_wbfmm(BEM3DFastMultipole solver,
				BEM3DFastMultipoleProblem problem,
				BEM3DParameters *param,
				BEM3DMeshSkeleton *s, 
				gdouble tol,
				gdouble *q, gdouble *dq,
				gdouble *p, gdouble *dp,
				BEM3DFMMWorkspace *w)
  
{
  g_assert(solver == BEM3D_FMM_WBFMM) ;
  g_assert(problem == BEM3D_FMM_HELMHOLTZ) ;

#ifdef HAVE_WBFMM
  wbfmm_tree_t *tree ;
  wbfmm_shift_operators_t *shifts ;
  gdouble xtree[3], xtmax[3], D, del, k, *xt, *normals, phi[2] ;
  gint depth, order[64], order_max, i, nda, nt ;
  gsize pstr ;
  guint64 b ;
    
  del = 1e-2 ;
  k = bem3d_parameters_wavenumber(param) ;
  depth = w->i[0] ;

  /*initialize everything if this is the first time through*/
  if ( w->p[0] == NULL ) {
    pstr = 3*sizeof(gdouble) ;
    
    wbfmm_points_origin_width(s->x, 3, s->npts, xtree, xtmax, &D, TRUE) ;

    xtree[0] -= del ; xtree[1] -= del ; xtree[2] -= del ;
    D += 2.0*del ;

    tree = wbfmm_tree_new(xtree, D, s->ns) ;
    w->p[0] = tree ;

    order_max = 0 ;
    for ( i = 1 ; i <= depth ; i ++ ) {
      order[2*i+0] = order[2*i+1] =
	wbfmm_truncation_number(tree, k, i, tol) ;
      order_max = MAX(order_max, order[2*i+0]) ;
      order_max = MAX(order_max, order[2*i+1]) ;      
    }

    nda = wbfmm_element_number_rotation(2*order_max) ;
    nda = MAX(nda, 16*(wbfmm_coefficient_index_nm(order_max+1,0))) ;

    g_assert(w->nda == 0) ;
    w->nda = 2*nda ;
    w->d = (gdouble *)g_malloc0(w->nda*sizeof(gdouble)) ;
    
    wbfmm_shift_angle_table_init() ;

    shifts = wbfmm_shift_operators_new(order_max, w->d) ;
    w->p[1] = shifts ;

    for ( i = 1 ; i <= depth ; i ++ ) {
      wbfmm_shift_operators_coaxial_SR_init(shifts, D, i, order[2*i+0], 
					    k, w->d) ;
    }

    for ( i = 2 ; i <= depth ; i ++ ) {
      wbfmm_shift_operators_coaxial_SS_init(shifts, D, i, 
					    order[2*(i-1)+0], 
					    k, w->d) ;
    }

    /*attach the points*/
    wbfmm_tree_add_points(tree, (gpointer)(s->x), s->ns, pstr) ;
    
    for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine(tree) ;

    for ( i = 1 ; i <= depth ; i ++ ) {
      wbfmm_tree_coefficient_init(tree, i, order[2*i+1], order[2*i+0]) ;
    }

  } else {
    tree = w->p[0] ; shifts = w->p[1] ;
  }

  normals = s->n ;
  if ( q == NULL ) {
    /*dipole sources only*/
    /* fprintf(stderr, "%s: matrix A\n", __FUNCTION__) ; */
    wbfmm_tree_leaf_expansions(tree, k, NULL, 2, normals, 3, dq, 2,
			       TRUE, w->d) ;
  }

  if ( dq == NULL ) {
    /*monopole sources only*/
    /* fprintf(stderr, "%s: matrix B\n", __FUNCTION__) ; */
    wbfmm_tree_leaf_expansions(tree, k, q, 2, NULL, 3, NULL, 2, TRUE, w->d) ;
  }

  /* fprintf(stderr, "%s: clearing coefficients\n", __FUNCTION__) ; */
  for ( i = 1 ; i < tree->depth ; i ++ ) 
    wbfmm_tree_coefficient_clear(tree, i) ;
  
  /*leaf expansions are computed as required, now perform FMM operations*/
  /* fprintf(stderr, "%s: upward pass\n", __FUNCTION__) ; */
  for ( i = tree->depth ; i >= 3 ; i -- ) 
    wbfmm_upward_pass(tree, shifts, i, w->d) ;

  /* fprintf(stderr, "%s: downward pass\n", __FUNCTION__) ; */
  for ( i = 2 ; i <= tree->depth ; i ++ )
    wbfmm_downward_pass(tree, shifts, i, w->d) ;

  /*run over targets and compute field*/
  nt = s->nt ; xt = &(s->x[3*(s->ns)]) ;

  if ( dq == NULL ) normals = NULL ;
  /* fprintf(stderr, "%s: local terms\n", __FUNCTION__) ; */
  for ( i = 0 ; i < nt ; i ++ ) {
    b = wbfmm_point_box(tree, tree->depth, &(xt[i*3])) ;
    phi[0] = phi[1] = 0.0 ;
    wbfmm_tree_box_local_field(tree, tree->depth, b, k, 
			       &(xt[i*3]), phi,
			       q, 2, normals, 3, dq, 2,
			       TRUE, w->d) ;
    p[2*i+0] = phi[0] ; p[2*i+1] = phi[1] ;
  }
  
#else /*HAVE_WBFMM*/
  g_error("%s: WBFMM support not compiled in", __FUNCTION__) ;
#endif /*HAVE_WBFMM*/

  return 0 ;
}
