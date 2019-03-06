/* radiation.c
 * 
 * Copyright (C) 2006, 2018 Michael Carley
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
 * @defgroup radiation Computing fields
 * @{
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdlib.h>
#include <string.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

/** 
 * Standard built-in radiation function for Laplace equation, which
 * can be used for any real problem
 * 
 * @param G integral of Green's function
 * @param dG integral of normal derivative of Green's function
 * @param phi solution at nodes
 * @param dphi normal derivative of solution at nodes
 * @param f radiated solution integral of \f$Gd\phi/dn+dG/dn\phi\f$
 * @param data user data, ignored in this function
 * 
 * @return 0 on success
*/

gint bem3d_radiation_func_laplace(GArray *G, GArray *dG,
				  GArray *phi, GArray *dphi,
				  GArray *f, gpointer data)

{
  gint i, j, stride ;

  g_debug("%s: ", __FUNCTION__) ;

  stride = G->len/phi->len ;
/*   g_array_set_size(f,stride) ; */
  if ( f->len < stride ) g_array_set_size(f,stride) ;
  for ( j = 0 ; j < stride ; j ++ ) {
    for ( i = 0 ; i < phi->len ; i ++ ) 
      g_array_index(f, gdouble, j) +=
	g_array_index(dG, gdouble, i*stride+j)*g_array_index(phi, gdouble, i) -
	g_array_index(G, gdouble, i*stride+j)*g_array_index(dphi, gdouble, i) ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Standard built-in radiation function for Helmholtz equation, which
 * can be used for any complex problem
 * 
 * @param G integral of Green's function
 * @param dG integral of normal derivative of Green's function
 * @param phi solution at nodes
 * @param dphi normal derivative of solution at nodes
 * @param f radiated solution integral of \f$Gd\phi/dn+dG/dn\phi\f$
 * @param data user data, ignored in this function
 * 
 * @return 0 on success
 */

gint bem3d_radiation_func_helmholtz(GArray *G, GArray *dG,
				    GArray *phi, GArray *dphi,
				    GArray *f, gpointer data)

{
  gint i ;

  g_debug("%s: ", __FUNCTION__) ;

  g_array_set_size(f,2) ;
  for ( i = 0 ; i < (G->len)/2 ; i ++ ) {
    g_array_index(f, gdouble, 0) +=
      g_array_index(dG, gdouble, 2*i)*g_array_index(phi, gdouble, 2*i) -
      g_array_index(dG, gdouble, 2*i+1)*g_array_index(phi, gdouble, 2*i+1) -
      g_array_index(G, gdouble, 2*i)*g_array_index(dphi, gdouble, 2*i) +
      g_array_index(G, gdouble, 2*i+1)*g_array_index(dphi, gdouble, 2*i+1) ;
    g_array_index(f, gdouble, 1) +=
      g_array_index(dG, gdouble, 2*i)*g_array_index(phi, gdouble, 2*i+1) +
      g_array_index(dG, gdouble, 2*i+1)*g_array_index(phi, gdouble, 2*i) -
      g_array_index(G, gdouble, 2*i)*g_array_index(dphi, gdouble, 2*i+1) -
      g_array_index(G, gdouble, 2*i+1)*g_array_index(dphi, gdouble, 2*i) ;
  }

  return BEM3D_SUCCESS ;
}

static void radiation_real(GArray *G, GArray *dG, gint nc,
			   GArray *phi, GArray *dphi,
			   GArray *f)
{
  gint i, j, idx ;

  /* for ( i = 0 ; i < G->len ; i ++ ) */
  /*   fprintf(stderr, "%lg ", g_array_index(G, gdouble, i)) ; */
  /* fprintf(stderr, "\n") ; */

  /* for ( i = 0 ; i < dG->len ; i ++ ) */
  /*   fprintf(stderr, "%lg ", g_array_index(dG, gdouble, i)) ; */
  /* fprintf(stderr, "\n") ; */

  /* for ( i = 0 ; i < phi->len ; i ++ ) */
  /*   fprintf(stderr, "%lg ", g_array_index(phi, gdouble, i)) ; */
  /* fprintf(stderr, "\n") ; */

  /* for ( i = 0 ; i < dphi->len ; i ++ ) */
  /*   fprintf(stderr, "%lg ", g_array_index(dphi, gdouble, i)) ; */
  /* fprintf(stderr, "\n") ; */
  
  for ( j = 0 ; j < nc ; j ++ ) {
    for ( i = 0 ; i < phi->len ; i ++ ) {
      idx = i*nc + j ;
      g_array_index(f, gdouble, j) +=
	g_array_index(dG, gdouble, idx)*
	g_array_index(phi, gdouble, i) -
	g_array_index(G, gdouble, idx)*
	g_array_index(dphi, gdouble, i) ;
    }
  }

  /* for ( i = 0 ; i < f->len ; i ++ ) */
  /*   fprintf(stderr, "%lg ", g_array_index(f, gdouble, i)) ; */
  /* fprintf(stderr, "\n") ; */

  /* exit(0) ; */
  
  return ;
}

static void radiation_complex(GArray *G, GArray *dG, gint nc,
			      GArray *phi, GArray *dphi,
			      GArray *f)
{
  gint i, j, idx ;

  for ( j = 0 ; j < nc ; j ++ ) {
    for ( i = 0 ; i < (phi->len)/2 ; i ++ ) {
      idx = nc*i + j ;
      g_array_index(f, gdouble, 2*j+0) +=
	g_array_index(dG, gdouble, 2*idx+0)*
	g_array_index(phi, gdouble, 2*i+0) -
	g_array_index(dG, gdouble, 2*idx+1)*
	g_array_index(phi, gdouble, 2*i+1) -
	g_array_index(G, gdouble, 2*idx+0)*
	g_array_index(dphi, gdouble, 2*i+0) +
	g_array_index(G, gdouble, 2*idx+1)*
	g_array_index(dphi, gdouble, 2*i+1) ;
      g_array_index(f, gdouble, 2*j+1) += 
	g_array_index(dG, gdouble, 2*idx+0)*
	g_array_index(phi, gdouble, 2*i+1) +
	g_array_index(dG, gdouble, 2*idx+1)*
	g_array_index(phi, gdouble, 2*i+0) -
	g_array_index(G, gdouble, 2*idx+0)*
	g_array_index(dphi, gdouble, 2*i+1) -
	g_array_index(G, gdouble, 2*idx+1)*
	g_array_index(dphi, gdouble, 2*i+0) ;
    }
  }
  return ;
}

static void element_radiation_point(BEM3DElement *e, 
				    BEM3DConfiguration *config,
				    BEM3DParameters *gdata,
				    BEM3DLookupFunc lfunc, gpointer ldata,
				    GtsPoint *x, 
				    GArray *G, GArray *dGdn,
				    GArray *phi, GArray *dphi,
				    GArray *f)
{
  BEM3DShapeFunc shfunc, cpfunc ;
  BEM3DGreensFunction gfunc ;
  static BEM3DQuadratureRule *q = NULL ;
  gdouble L[32], dLds[32], dLdt[32], g[16], dgdn[16] ;
  BEM3DQuadratureRuleFunc qfunc ;
  gpointer qdata ;
  gint i, j, k, stride = 0, nc ;
  gdouble s, t, J, w, wt ;
  GtsPoint y ;
  GtsVector n ;

  g_debug("%s: e=%p", __FUNCTION__, e) ;

  shfunc = bem3d_element_shape_func(e) ;
  cpfunc = bem3d_element_node_func(e) ;
  
  if ( q == NULL ) q = bem3d_quadrature_rule_new(1024, 1) ;

  g_debug("%s: %d point quadrature rule", __FUNCTION__,
	  bem3d_quadrature_vertex_number(q)) ;

  gfunc = config->gfunc ;

  qfunc = config->qrule ; qdata = config->qdata ;

  qfunc(x, e, q, &gfunc, gdata, qdata) ;

  /* get the size of the Green's function */
  if ( bem3d_greens_function_is_real(&gfunc) ) stride = 1 ;
  else stride = 2 ;
  nc = bem3d_greens_function_component_number(&gfunc) ;
  
  g_array_set_size(G,bem3d_element_node_number(e)*stride*nc) ;
  g_array_set_size(dGdn,bem3d_element_node_number(e)*stride*nc) ;
  g_array_set_size(phi,bem3d_element_node_number(e)*stride) ; 
  g_array_set_size(dphi,bem3d_element_node_number(e)*stride) ;  
  g_array_set_size(f,stride*nc) ; 

  /* g_array_set_size(G, stride*bem3d_element_node_number(e)*nc) ; */
  /* g_array_set_size(dGdn, stride*bem3d_element_node_number(e)*nc) ; */

  for ( i = 0 ; i < G->len ; i ++ ) 
    g_array_index(G,gdouble,i) = g_array_index(dGdn,gdouble,i) = 0.0 ;

  for ( i = 0 ; i < bem3d_quadrature_vertex_number(q) ; i ++ ) {
    g_assert(i < bem3d_quadrature_vertex_number_max(q)) ;
    s = bem3d_quadrature_xi(q,i) ;
    t = bem3d_quadrature_eta(q,i) ;
    wt = bem3d_quadrature_weight(q,i) ;
    shfunc(s, t, L, dLds, dLdt, NULL) ;
    bem3d_element_position(e, L, &y) ;
    bem3d_element_normal(e, dLds, dLdt, n, &J) ;
    bem3d_greens_function_func(&gfunc)(x, &y, n, gdata, g, dgdn) ;
    cpfunc(s, t, L, NULL, NULL, NULL) ;
    J *= wt ;
    for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
      w = J*L[j] ;
      for ( k = 0 ; k < stride*nc ; k ++ ) {
  	g_array_index(G,gdouble,j*stride*nc+k) += g[k]*w ;
 	  /* g_array_index(g,gdouble,k)*w ; */
	g_array_index(dGdn,gdouble,j*stride*nc+k) += dgdn[k]*w ;
	  /* g_array_index(dgdn,gdouble,k)*w ; */
      }
    }
    /* fprintf(stderr, "%lg\n", w) ; */
    
    /* fprintf(stderr, "G: ") ; */
    /* for ( j = 0 ; j < G->len ; j ++ ) */
    /*   fprintf(stderr, "%lg ", g_array_index(G,gdouble,j)) ; */
    /* fprintf(stderr, "\n") ; */
    /* fprintf(stderr, "g: ") ; */
    /* for ( j = 0 ; j < g->len ; j ++ ) */
    /*   fprintf(stderr, "%lg ", g_array_index(g,gdouble,j)) ; */
    /* fprintf(stderr, "\n") ; */
    /* fprintf(stderr, "dG: ") ; */
    /* for ( j = 0 ; j < dGdn->len ; j ++ ) */
    /*   fprintf(stderr, "%lg ", g_array_index(dGdn,gdouble,j)) ; */
    /* fprintf(stderr, "\n") ; */
    /* fprintf(stderr, "dgdn: ") ; */
    /* for ( j = 0 ; j < g->len ; j ++ ) */
    /*   fprintf(stderr, "%lg ", g_array_index(dgdn,gdouble,j)) ; */
    /* fprintf(stderr, "\n") ; */

    /* exit(0) ; */
  }

  if ( bem3d_quadrature_free_number(q) != 0 ) {
    if ( bem3d_quadrature_free_number(q) != bem3d_element_node_number(e)*nc )
      g_error("%s: number of free terms (%d) in quadrature rule "
	      "does not match number of nodes on element (%d)",
	      __FUNCTION__, bem3d_quadrature_free_number(q), 
	      bem3d_element_node_number(e)) ;
    /*we should not arrive here in complex problems*/
    g_assert(q->wfree == 1) ;
    /* g_assert(stride == 1) ; */
    k = 0 ;
    for ( j = 0 ; j < bem3d_element_node_number(e)*nc*stride ; j ++ ) {
      g_array_index(G,gdouble,j) += bem3d_quadrature_free_term_g(q,j) ;
      g_array_index(dGdn,gdouble,j) += bem3d_quadrature_free_term_dg(q,j) ;
      /* g_array_index(G,gdouble,j*stride+k) +=  */
      /* 	bem3d_quadrature_free_term_g(q,j) ; */
      /* g_array_index(dGdn,gdouble,j*stride+k) +=  */
      /* 	bem3d_quadrature_free_term_dg(q,j) ; */
    }
  }

  g_array_set_size(phi, stride*bem3d_element_node_number(e)) ;
  g_array_set_size(dphi, stride*bem3d_element_node_number(e)) ;
  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
    lfunc(bem3d_element_global_index(e, j), j, ldata, phi, dphi) ;
  }

  if ( bem3d_greens_function_is_real(&gfunc) ) {
    radiation_real(G, dGdn, nc, phi, dphi, f) ;
  } else {
    radiation_complex(G, dGdn, nc, phi, dphi, f) ;
  }
  
  return ;
}

static gint mesh_radiation_point(BEM3DElement *e, gpointer data[])

{
  GtsPoint *x = (GtsPoint *)data[BEM3D_BMESH_DATA_POINT] ;
  BEM3DConfiguration *config = data[BEM3D_BMESH_DATA_CONFIG] ; 
  BEM3DParameters *gdata = data[BEM3D_BMESH_DATA_GDATA] ;
  BEM3DLookupFunc lf = (BEM3DLookupFunc)data[BEM3D_BMESH_DATA_LFUNC] ;
  gpointer ldata = data[BEM3D_BMESH_DATA_LDATA] ;
  GArray *f = (GArray *)data[BEM3D_BMESH_DATA_F] ;
  static GArray *g = NULL, *dg = NULL, *phi = NULL, *dphi = NULL ;

  if ( g == NULL ) {
    g = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
    dg = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
    phi = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
    dphi = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
  }

  element_radiation_point(e, config, gdata, lf, ldata,
			  x, g, dg, phi, dphi, f) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Compute the field radiated to a point from a BEM3DMesh with known
 * solution
 * 
 * @param m BEM3DMesh for geometry;
 * @param config ::BEM3DConfiguration for problem;
 * @param gdata ::BEM3DParameters to pass to Green's function;
 * @param lf lookup function for solution
 * @param ldata user data to pass to lf
 * @param x field point
 * @param f radiated solution
 * 
 * @return 0 on success
 */

gint bem3d_mesh_radiation_point(BEM3DMesh *m,
				BEM3DConfiguration *config,
				BEM3DParameters *gdata,
				BEM3DLookupFunc lf, gpointer ldata,
				GtsPoint *x, GArray *f)

{
  gpointer data[BEM3D_BMESH_DATA_WIDTH] ;
  gint j, len ;
  BEM3DGreensFunction gfunc ;

  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(config != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(gdata != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(lf != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  data[BEM3D_BMESH_DATA_MESH] = m ; 
  data[BEM3D_BMESH_DATA_POINT] = x ; 
  data[BEM3D_BMESH_DATA_CONFIG] = config ; 
  data[BEM3D_BMESH_DATA_GDATA] = gdata ;
  data[BEM3D_BMESH_DATA_LFUNC] = lf ;
  data[BEM3D_BMESH_DATA_LDATA] = ldata ;
  data[BEM3D_BMESH_DATA_IMIN] = &(bem3d_mesh_node_index_min(m)) ;
  data[BEM3D_BMESH_DATA_IMAX] = &(bem3d_mesh_node_index_max(m)) ;
  data[BEM3D_BMESH_DATA_F] = f ;  

  gfunc = config->gfunc ;
  len = bem3d_greens_function_component_number(&gfunc) ;
  if ( !bem3d_greens_function_is_real(&gfunc) ) len *= 2 ;
  
  g_array_set_size(f, len) ;
  for ( j = 0 ; j < f->len ; j ++ ) g_array_index(f,gdouble,j) = 0.0 ;  
  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)mesh_radiation_point,
			     data) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Compute the field radiated to a point from a BEM3DElement with known
 * solution
 * 
 * @param e BEM3DElement
 * @param config ::BEM3DConfiguration for problem;
 * @param gdata ::BEM3DParameters to pass to Green's function;
 * @param lfunc lookup function for solution
 * @param ldata user data to pass to lf
 * @param x field point
 * @param G integral of Green's function
 * @param dGdn integral of Green's function normal derivative
 * @param phi solution at element nodes
 * @param dphi solution normal derivative at element nodes
 * 
 * @return 0 on success
 */

gint bem3d_element_radiation_point(BEM3DElement *e, 
				   BEM3DConfiguration *config,
				   BEM3DParameters *gdata,
				   BEM3DLookupFunc lfunc, gpointer ldata,
				   GtsPoint *x, 
				   GArray *G, GArray *dGdn,
				   GArray *phi, GArray *dphi)
{
  BEM3DGreensFunction gfunc ;
  BEM3DShapeFunc shfunc, cpfunc ;
  static BEM3DQuadratureRule *q = NULL ;
  gdouble L[32], dLds[32], dLdt[32], g[16], dgdn[16] ;
  BEM3DQuadratureRuleFunc qfunc ;
  gpointer qdata ;
  gint i, j, k, stride = 0, nc ;
  gdouble s, t, J, w, wt ;
  GtsPoint y ;
  GtsVector n ;

  g_debug("%s: e=%p", __FUNCTION__, e) ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(config != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(gdata != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(lfunc != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(phi != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dphi != NULL, BEM3D_NULL_ARGUMENT) ;

  shfunc = bem3d_element_shape_func(e) ;
  cpfunc = bem3d_element_node_func(e) ;

  if ( q == NULL ) q = bem3d_quadrature_rule_new(1024, 1) ;

  gfunc = config->gfunc ;
  qfunc = config->qrule ; qdata = config->qdata ;

  qfunc(x, e, q, &gfunc, gdata, qdata) ;

  g_debug("%s: %d point quadrature rule", __FUNCTION__,
	  bem3d_quadrature_vertex_number(q)) ;
  
  /* get the size of the Green's function */
  if ( bem3d_greens_function_is_real(&gfunc) ) stride = 1 ;
  else stride = 2 ;
  nc = bem3d_quadrature_component_number(q) ;

  g_array_set_size(G,bem3d_element_node_number(e)*stride*nc) ; 
  g_array_set_size(dGdn,bem3d_element_node_number(e)*stride*nc) ;  
  g_array_set_size(phi,bem3d_element_node_number(e)*stride) ; 
  g_array_set_size(dphi,bem3d_element_node_number(e)*stride) ;  

  for ( i = 0 ; i < G->len ; i ++ ) 
    g_array_index(G,gdouble,i) = g_array_index(dGdn,gdouble,i) = 0.0 ;

  if ( bem3d_greens_function_is_real(&gfunc) ) stride = 1 ;
  else stride = 2 ;
  stride *= bem3d_greens_function_component_number(&gfunc) ;
  
  for ( i = 0 ; i < bem3d_quadrature_vertex_number(q) ; i ++ ) {
    g_assert(i < bem3d_quadrature_vertex_number_max(q)) ;
    g_assert(nc == 1) ;
    s = bem3d_quadrature_xi(q,i) ;
    t = bem3d_quadrature_eta(q,i) ;
    wt = bem3d_quadrature_weight(q,i) ;
    shfunc(s, t, L, dLds, dLdt, NULL) ;
    bem3d_element_position(e, L, &y) ;
    bem3d_element_normal(e, dLds, dLdt, n, &J) ;
    bem3d_greens_function_func(&gfunc)(x, &y, n, gdata, g, dgdn) ;
    cpfunc(s, t, L, NULL, NULL, NULL) ;
    J *= wt ;

    for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
      w = J*L[j] ;
      for ( k = 0 ; k < stride ; k ++ ) {
  	g_array_index(G,gdouble,j*stride+k) += g[k]*w ;
	g_array_index(dGdn,gdouble,j*stride+k) += dgdn[k]*w ;
      }
    }

  }

  if ( bem3d_quadrature_free_number(q) != 0 ) {
    if ( bem3d_quadrature_free_number(q) != bem3d_element_node_number(e) )
      g_error("%s: number of free terms (%d) in quadrature rule "
	      "does not match number of nodes on element (%d)",
	      __FUNCTION__, bem3d_quadrature_free_number(q), 
	      bem3d_element_node_number(e)) ;
    /*we should not arrive here in complex problems*/
    g_assert(q->wfree == 1) ; g_assert(stride == 1) ;
    k = 0 ;
    for ( j = 0 ; j < bem3d_element_node_number(e)*nc ; j ++ ) {
      g_array_index(G,gdouble,j*stride+k) += 
	bem3d_quadrature_free_term_g(q,j) ;
      g_array_index(dGdn,gdouble,j*stride+k) += 
	bem3d_quadrature_free_term_dg(q,j) ;
    }
  }

  g_array_set_size(phi, bem3d_element_node_number(e)) ;
  g_array_set_size(dphi, bem3d_element_node_number(e)) ;
  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ )
    lfunc(bem3d_element_global_index(e, j), j, ldata, phi, dphi) ;

  return BEM3D_SUCCESS ;
}

static gint mesh_radiation_mesh(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMesh *m = (BEM3DMesh *)data[BEM3D_BMESH_DATA_MESH] ;
  BEM3DConfiguration *config = data[BEM3D_BMESH_DATA_CONFIG] ; 
  BEM3DParameters *gdata = data[BEM3D_BMESH_DATA_GDATA] ;
  BEM3DLookupFunc lf = (BEM3DLookupFunc)data[BEM3D_BMESH_DATA_LFUNC] ;
  gpointer ldata = data[BEM3D_BMESH_DATA_LDATA] ;
  BEM3DMeshData *fm = (BEM3DMeshData *)data[BEM3D_BMESH_DATA_F] ;
  static GArray *f = NULL ;

  if ( f == NULL ) {
    f = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
  }

  bem3d_mesh_radiation_point(m, config, gdata, lf, ldata, GTS_POINT(v), f) ;
  bem3d_mesh_data_add(fm, i, f) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Compute the field radiated to a BEM3DMesh from a BEM3DMesh with
 * known solution
 * 
 * @param m BEM3DMesh for geometry;
 * @param config ::BEM3DConfiguration for problem;
 * @param gdata ::BEM3DParameters to pass to Green's function;
 * @param lfunc lookup function for solution
 * @param ldata user data to pass to lfunc
 * @param s BEM3DMesh on which to compute field
 * @param f BEM3DMeshData block containing computed field
 * 
 * @return 0 on success 
 */

gint bem3d_mesh_radiation_mesh(BEM3DMesh *m,
			       BEM3DConfiguration *config,
			       BEM3DParameters *gdata,
			       BEM3DLookupFunc lfunc, gpointer ldata,
			       BEM3DMesh *s, BEM3DMeshData *f)

{
  gpointer data[BEM3D_BMESH_DATA_WIDTH] ;

  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(config != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(gdata != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(lfunc != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(s != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(s), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  data[BEM3D_BMESH_DATA_MESH] = m ; 
  data[BEM3D_BMESH_DATA_CONFIG] = config ; 
  data[BEM3D_BMESH_DATA_GDATA] = gdata ;
  data[BEM3D_BMESH_DATA_LFUNC] = lfunc ;
  data[BEM3D_BMESH_DATA_LDATA] = ldata ;
  data[BEM3D_BMESH_DATA_IMIN] = &(bem3d_mesh_node_index_min(m)) ;
  data[BEM3D_BMESH_DATA_IMAX] = &(bem3d_mesh_node_index_max(m)) ;
  data[BEM3D_BMESH_DATA_F] = f ;  

  bem3d_mesh_foreach_node(s, (BEM3DNodeFunc)mesh_radiation_mesh,
			  data) ;

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
