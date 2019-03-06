/* reduction.c
 * 
 * Copyright (C) 2018 Michael Carley
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

#include "bem3d.h"
#include "bem3d-private.h"

static gchar *_operations[] = {"max",
			      "min",
			      "sum",
			      "int",
			      ""} ;
static gchar *_descriptions[] = {"maximum value in each column",
				 "minimum value in each column",
				 "sum of values in each column",
				 "integral of data over BEM3D surface",
				 ""} ;
static BEM3DReductionFunc _funcs[] = {bem3d_reduction_func_max,
				      bem3d_reduction_func_min,
				      bem3d_reduction_func_sum,
				      bem3d_reduction_func_int,
				      NULL} ;

/**
 * @defgroup reduction Functions for data reduction
 *
 * @{
 */

/**
 * Maximum value of entries in each column of data
 *
 * @param m a ::BEM3DMesh
 * @param f a ::BEM3DMeshData containing data for \a m
 * @param idata on exit, index of maximum value in each column
 * @param ni on exit 1
 * @param ddata maximum value in each column of \a f
 * @param nd on exit 1
 * @param data ignored
 *
 * @return 0 on success
 */

gint bem3d_reduction_func_max(BEM3DMesh *m, BEM3DMeshData *f,
			      gint *idata, gint *ni,
			      gdouble *ddata, gint *nd,
			      gpointer data)

{
  gint i, j, nf ;
  gdouble *d ;

  nf = bem3d_mesh_data_element_number(f) ;

  for ( j = 0 ; j < nf ; j ++ ) {
    idata[j] = 0.0 ; ddata[j] = -G_MAXDOUBLE ;
  }

  *ni = *nd = 1 ;
  
  for ( i = 0 ; i < bem3d_mesh_node_number(m) ; i ++ ) {
    d = bem3d_mesh_data_get(f, i) ;
    for ( j = 0 ; j < nf ; j ++ ) {
      if ( d[j] > ddata[j] ) {
	ddata[j] = d[j] ; idata[j] = i ;
      }
    }
  }
  
  return 0 ;
}

/**
 * Minimum value of entries in each column of data
 *
 * @param m a ::BEM3DMesh
 * @param f a ::BEM3DMeshData containing data for \a m
 * @param idata on exit, index of minimum value in each column
 * @param ni on exit 1
 * @param ddata minimum value in each column of \a f
 * @param nd on exit 1
 * @param data ignored
 *
 * @return 0 on success
 */

gint bem3d_reduction_func_min(BEM3DMesh *m, BEM3DMeshData *f,
			      gint *idata, gint *ni,
			      gdouble *ddata, gint *nd,
			      gpointer data)

{
  gint i, j, nf ;
  gdouble *d ;

  nf = bem3d_mesh_data_element_number(f) ;

  for ( j = 0 ; j < nf ; j ++ ) {
    idata[j] = 0.0 ; ddata[j] = G_MAXDOUBLE ;
  }

  *ni = *nd = 1 ;
  
  for ( i = 0 ; i < bem3d_mesh_node_number(m) ; i ++ ) {
    d = bem3d_mesh_data_get(f, i) ;
    for ( j = 0 ; j < nf ; j ++ ) {
      if ( d[j] < ddata[j] ) {
	ddata[j] = d[j] ; idata[j] = i ;
      }
    }
  }
  
  return 0 ;
}

/**
 * Sum of entries in each column
 *
 * @param m a ::BEM3DMesh
 * @param f a ::BEM3DMeshData containing data for \a m
 * @param idata ignored
 * @param ni ignored
 * @param ddata sum of entries in each column of \a f
 * @param nd on exit 1
 * @param data ignored
 *
 * @return 0 on success
 */

gint bem3d_reduction_func_sum(BEM3DMesh *m, BEM3DMeshData *f,
			      gint *idata, gint *ni,
			      gdouble *ddata, gint *nd,
			      gpointer data)

{
  gint i, j, nf ;
  gdouble *d ;

  nf = bem3d_mesh_data_element_number(f) ;

  for ( j = 0 ; j < nf ; j ++ ) {
    ddata[j] = 0.0 ;
  }

  *ni = 0 ; *nd = 1 ;
  
  for ( i = 0 ; i < bem3d_mesh_node_number(m) ; i ++ ) {
    d = bem3d_mesh_data_get(f, i) ;
    for ( j = 0 ; j < nf ; j ++ ) {
      ddata[j] += d[j] ;
    }
  }
  
  return 0 ;
}

static gint integrate_func(BEM3DElement *e, gpointer *data)

{
  BEM3DMeshData *f = data[0] ;
  gdouble *ddata = data[1] ;
  gint nf = *((gint *)(data[2])) ;
  BEM3DQuadratureRule *q = data[3] ;
  gint ngp = *((gint *)(data[4])) ;
  BEM3DShapeFunc shfunc = bem3d_element_shape_func(e) ;
  BEM3DShapeFunc cpfunc = bem3d_element_node_func(e) ;
  gdouble L[32], dLds[32], dLdt[32], s, t, wt, J, n[3], *d ;
  gint i, j, k ;
  
  bem3d_quadrature_rule_wx(NULL, e, q, NULL, NULL, &ngp) ;

  for ( i = 0 ; i < bem3d_quadrature_vertex_number(q) ; i ++ ) {
    s = bem3d_quadrature_xi(q,i) ;
    t = bem3d_quadrature_eta(q,i) ;
    wt = bem3d_quadrature_weight(q,i) ;
    shfunc(s, t, L, dLds, dLdt, NULL) ;
    bem3d_element_normal(e, dLds, dLdt, n, &J) ;
    J *= wt ;
    cpfunc(s, t, L, NULL, NULL, NULL) ;
    for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
      d = bem3d_mesh_data_get(f, bem3d_element_global_index(e, j)) ;
      for ( k = 0 ; k < nf ; k ++ ) {
	ddata[k] += d[k]*J*L[j] ;
      }
    }
  }
  
  return 0 ;
}

/**
 * Surface integral of data
 *
 * @param m a ::BEM3DMesh
 * @param f a ::BEM3DMeshData containing data for \a m
 * @param idata ignored
 * @param ni ignored
 * @param ddata surface integral on \a m of data in each column of \a f
 * @param nd on exit 1
 * @param data ignored
 *
 * @return 0 on success
 */

gint bem3d_reduction_func_int(BEM3DMesh *m, BEM3DMeshData *f,
			      gint *idata, gint *ni,
			      gdouble *ddata, gint *nd,
			      gpointer data)

{
  gpointer mdata[8] ;
  gint j, nf, ngp ;
  BEM3DQuadratureRule *rule ;

  ngp = 7 ;
  nf = bem3d_mesh_data_element_number(f) ;

  for ( j = 0 ; j < nf ; j ++ ) ddata[j] = 0.0 ;
  *ni = 0 ; *nd = 1 ;

  rule = bem3d_quadrature_rule_new(ngp, 1) ;

  mdata[0] = f ;
  mdata[1] = ddata ;
  mdata[2] = &nf ;
  mdata[3] = rule ;
  mdata[4] = &ngp ;

  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)integrate_func,
			     (gpointer)mdata) ;
  
  bem3d_quadrature_rule_free(rule) ;
  
  return 0 ;
}

static BEM3DReductionFunc parse_reduction(gchar *str)

{
  gint i ;

  for ( i = 0 ; strlen(_operations[i]) != 0 ; i ++ ) {
    if ( strcmp(_operations[i], str) == 0 ) return _funcs[i] ;
  }
  
  return NULL ;
}

/**
 * Apply a named reduction operation to a ::BEM3DMeshData on a
 * ::BEM3DMesh. Reduction operations can be listed by calling
 * ::bem3d_reduction_func_list.
 *
 * @param func the name of reduction operation
 * @param m a ::BEM3DMesh
 * @param f a ::BEM3DMeshData containing data for \a m
 * @param idata array containing integer data on exit
 * @param ni number of entries in \a idata on exit
 * @param ddata array containing double precision data on exit
 * @param nd number of entries in \a ddata on exit
 * @param data ignored
 *
 * @return 0 on success
 */

gint bem3d_reduction_func_apply(gchar *func,
				BEM3DMesh *m, BEM3DMeshData *f,
				gint *idata, gint *ni,
				gdouble *ddata, gint *nd,
				gpointer data)

{
  BEM3DReductionFunc rfunc ;

  rfunc = parse_reduction(func) ;
  if ( rfunc == NULL )
    g_error("%s: cannot parse reduction operation \"%s\"\n",
	    __FUNCTION__, func) ;

  rfunc(m, f, idata, ni, ddata, nd, data) ;

  return 0 ;
}

/**
 * Short description of a reduction operation.
 *
 * @param str the name of the reduction operation, e.g. `max'
 *
 * @return a pointer to a string containing a description of the
 * operation, or NULL if not found.
 */

gchar *bem3d_reduction_func_description(gchar *str)

{
  gint i ;
  
  for ( i = 0 ; strlen(_operations[i]) != 0 ; i ++ ) {
    if ( strcmp(_operations[i], str) == 0 ) return _descriptions[i] ;
  }

  return NULL ;
}

/**
 * Write a list of available reduction operations.
 *
 * @param output a FILE stream to write information to;
 * @param format a `printf' format string for the data;
 * @param describe if TRUE, also write the short description of the 
 * reduction operation 
 *
 * The format string should be appropriate to the setting of \a
 * describe. If \a describe is FALSE, \a format should contain one
 * entry of `%s' (for the name of the operation); if TRUE, \a format
 * should contain `%s' twice (once for the name of the operation and
 * once for its description).
 *
 * @return 0 on success
 */

gint bem3d_reduction_func_list(FILE *output, gchar *format, gboolean describe)

{
  gint i ;

  if ( describe ) {
    for ( i = 0 ; strlen(_operations[i]) != 0 ; i ++ ) {
      fprintf(output, format, _operations[i], _descriptions[i]) ;
    }
    return 0 ;
  }

  for ( i = 0 ; strlen(_operations[i]) != 0 ; i ++ ) {
    fprintf(output, format, _operations[i]) ;
  }

  return 0 ;
}

/**
 * @}
 * 
 */
