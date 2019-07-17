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

static void radiation_real(gdouble *G, gdouble *dG, gint nc,
			   gint len, gdouble *phi, gdouble *dphi,
			   gdouble *f)

{
  gint i, j, idx ;
  
  for ( j = 0 ; j < nc ; j ++ ) {
    for ( i = 0 ; i < len ; i ++ ) {
      idx = i*nc + j ;
      f[j] += dG[idx]*phi[i] - G[idx]*dphi[i] ;
    }
  }
  
  return ;
}

static void radiation_complex(gdouble *G, gdouble *dG, gint nc, gint len,
			      gdouble *phi, gdouble *dphi, gdouble *f)

{
  gint i, j, idx ;

  for ( j = 0 ; j < nc ; j ++ ) {
    for ( i = 0 ; i < len/2 ; i ++ ) {
      idx = nc*i + j ;
      f[2*j+0] +=
	dG[2*idx+0]*phi[ 2*i+0] - dG[2*idx+1]*phi[ 2*i+1] -
	G[ 2*idx+0]*dphi[2*i+0] +  G[2*idx+1]*dphi[2*i+1] ;
      f[2*j+1] += 
	dG[2*idx+0]*phi[2*i+1] + dG[2*idx+1]*phi[2*i+0] -
	G[2*idx+0]*dphi[2*i+1] - G[2*idx+1]*dphi[2*i+0] ;
    }
  }
  return ;
}

static void element_radiation_point(BEM3DElement *e, 
				    BEM3DConfiguration *config,
				    BEM3DParameters *gdata,
				    BEM3DLookupFunc lfunc, gpointer ldata,
				    GtsPoint *x, 
				    GArray *f,
				    BEM3DWorkspace *work)

{
  BEM3DShapeFunc shfunc, cpfunc ;
  BEM3DGreensFunction gfunc ;
  GArray *phi, *dphi ;
  BEM3DQuadratureRule *q = work->q ;
  gdouble L[64], dLds[64], dLdt[64], g[16], dgdn[16] ;
  gdouble G[64], dGdn[64] ;
  BEM3DQuadratureRuleFunc qfunc ;
  gpointer qdata ;
  gint i, j, k, stride = 0, nc, ne ;
  gdouble s, t, J, w, wt ;
  GtsPoint y ;
  GtsVector n ;

  g_debug("%s: e=%p", __FUNCTION__, e) ;

  g_assert(bem3d_element_node_number(e) < 17) ;

   phi = bem3d_workspace_double_array_get(work) ;
  dphi = bem3d_workspace_double_array_get(work) ;
  
  shfunc = bem3d_element_shape_func(e) ;
  cpfunc = bem3d_element_node_func(e) ;
  
  g_debug("%s: %d point quadrature rule", __FUNCTION__,
	  bem3d_quadrature_vertex_number(q)) ;

  gfunc = config->gfunc ;

  qfunc = config->qrule ; qdata = config->qdata ;

  qfunc(x, e, q, &gfunc, gdata, qdata, work) ;

  /* get the size of the Green's function */
  if ( bem3d_greens_function_is_real(&gfunc) ) stride = 1 ;
  else stride = 2 ;
  nc = bem3d_greens_function_component_number(&gfunc) ;
  
  g_array_set_size(f,stride*nc) ; 

  ne = bem3d_element_node_number(e)*stride*nc ;
  memset(G,    0, ne*sizeof(gdouble)) ;
  memset(dGdn, 0, ne*sizeof(gdouble)) ;
  
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
  	G   [j*stride*nc+k] += g[k]*w ;
	dGdn[j*stride*nc+k] += dgdn[k]*w ;
      }
    }
  }

  if ( bem3d_quadrature_free_number(q) != 0 ) {
    if ( bem3d_quadrature_free_number(q) != bem3d_element_node_number(e)*nc )
      g_error("%s: number of free terms (%d) in quadrature rule "
	      "does not match number of nodes on element (%d)",
	      __FUNCTION__, bem3d_quadrature_free_number(q), 
	      bem3d_element_node_number(e)) ;
    /*we should not arrive here in complex problems*/
    g_assert(q->wfree == 1) ;
    k = 0 ;
    for ( j = 0 ; j < bem3d_element_node_number(e)*nc*stride ; j ++ ) {
      G   [j] += bem3d_quadrature_free_term_g(q,j) ;
      dGdn[j] += bem3d_quadrature_free_term_dg(q,j) ;
    }
  }

  g_array_set_size(phi, stride*bem3d_element_node_number(e)) ;
  g_array_set_size(dphi, stride*bem3d_element_node_number(e)) ;
  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
    lfunc(bem3d_element_global_index(e, j), j, ldata, phi, dphi) ;
  }

  if ( bem3d_greens_function_is_real(&gfunc) ) {
    radiation_real(G, dGdn, nc, phi->len,
		   &(g_array_index(phi, gdouble, 0)),
		   &(g_array_index(dphi, gdouble, 0)),
		   &(g_array_index(f, gdouble, 0))) ;
  } else {
    radiation_complex(G, dGdn, nc, phi->len,
		   &(g_array_index(phi, gdouble, 0)),
		   &(g_array_index(dphi, gdouble, 0)),
		   &(g_array_index(f, gdouble, 0))) ;
  }
  
  bem3d_workspace_double_array_put(work,  phi) ;
  bem3d_workspace_double_array_put(work, dphi) ;

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
  BEM3DWorkspace *work = data[BEM3D_BMESH_DATA_WORK] ;

  element_radiation_point(e, config, gdata, lf, ldata, x, f, work) ;

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
 * @param work a ::BEM3DWorkspace
 * 
 * @return 0 on success
 */

gint bem3d_mesh_radiation_point(BEM3DMesh *m,
				BEM3DConfiguration *config,
				BEM3DParameters *gdata,
				BEM3DLookupFunc lf, gpointer ldata,
				GtsPoint *x, GArray *f,
				BEM3DWorkspace *work)

{
  gpointer data[BEM3D_BMESH_DATA_WIDTH] = {NULL} ;
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
  data[BEM3D_BMESH_DATA_WORK] = work ;

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
 * @param work a ::BEM3DWorkspace
 * 
 * @return 0 on success
 */

gint bem3d_element_radiation_point(BEM3DElement *e, 
				   BEM3DConfiguration *config,
				   BEM3DParameters *gdata,
				   BEM3DLookupFunc lfunc, gpointer ldata,
				   GtsPoint *x, 
				   GArray *G, GArray *dGdn,
				   GArray *phi, GArray *dphi,
				   BEM3DWorkspace *work)
{
  BEM3DGreensFunction gfunc ;
  BEM3DShapeFunc shfunc, cpfunc ;
  BEM3DQuadratureRule *q = work->q ;
  gdouble L[64], dLds[64], dLdt[64], g[16], dgdn[16] ;
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

  gfunc = config->gfunc ;
  qfunc = config->qrule ; qdata = config->qdata ;

  qfunc(x, e, q, &gfunc, gdata, qdata, work) ;

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
  BEM3DWorkspace *work = data[BEM3D_BMESH_DATA_WORK] ;
  GArray *f ;
  
  g_assert(work->used_doubles[2] == FALSE) ;
  work->used_doubles[2] = TRUE ;
  f = work->doubles[2] ;
  
  bem3d_mesh_radiation_point(m, config, gdata, lf, ldata, GTS_POINT(v),
			     f, work) ;
  bem3d_mesh_data_add(fm, i, f) ;

  work->used_doubles[2] = FALSE ;

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
 * @param work a ::BEM3DWorkspace
 * 
 * @return 0 on success 
 */

gint bem3d_mesh_radiation_mesh(BEM3DMesh *m,
			       BEM3DConfiguration *config,
			       BEM3DParameters *gdata,
			       BEM3DLookupFunc lfunc, gpointer ldata,
			       BEM3DMesh *s, BEM3DMeshData *f,
			       BEM3DWorkspace *work)

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
  data[BEM3D_BMESH_DATA_WORK] = work ;
  
  bem3d_mesh_foreach_node(s, (BEM3DNodeFunc)mesh_radiation_mesh,
			  data) ;

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
