/* belement.c
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

#include "bem3d.h"
#include "bem3d-private.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>

static void element_destroy(GtsObject *object)

{
  BEM3DElement *e = BEM3D_ELEMENT(object) ;

  g_free(e->i) ;  g_free(e->s) ; 
  g_free(e->f) ;  
  if ( e->c != e->v ) g_free(e->c) ;
  g_free(e->v) ;
  if ( e->xc != e->xs ) g_free(e->xc) ;
  g_free(e->xs) ; 
  

  (*GTS_OBJECT_CLASS(bem3d_element_class ())->parent_class->destroy)
    (object) ;  
  return ;
}

static void bem3d_element_class_init (BEM3DElementClass * klass)
{
  /* define new methods and overload inherited methods here */
  GTS_OBJECT_CLASS (klass)->destroy = element_destroy ;
}


static void bem3d_element_init (BEM3DElement * object)
{
  /* initialize object here */
}

/**
 * @defgroup belement BEM3D elements
 *
 * The definitions of the basic ::BEM3DElement and the functions for
 * handling it. To see how to define new types of element, consult
 * \link elements Building elements \endlink
 *
 * @{
 * 
 */

/** 
 * The basic class for a ::BEM3DElement. 
 * 
 * @return the ::BEM3DElementClass
 */

BEM3DElementClass * bem3d_element_class (void)
{
  static BEM3DElementClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo bem3d_element_info = {
      "BEM3DElement",
      sizeof (BEM3DElement),
      sizeof (BEM3DElementClass),
      (GtsObjectClassInitFunc) bem3d_element_class_init,
      (GtsObjectInitFunc) bem3d_element_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &bem3d_element_info);
  }

  return klass;
}

/** 
 * Allocate a new ::BEM3DElement
 * 
 * @param klass ::BEM3DElementClass (use ::bem3d_element_class());
 * @param nf number of GtsFace on the element;
 * @param nv number of vertices for geometric interpolation; 
 * @param nc number of collocation points (not necessarily the same as nv);
 * set nc to zero make the collocation points coincide with the geometric
 * vertices;
 * @param ns number of sides/corners (e.g. three for triangular elements);
 * @param shf pointer to shape function for geometry;
 * @param cpf pointer to shape function for collocation.
 * 
 * @return pointer to new element.
 */

BEM3DElement *bem3d_element_new(BEM3DElementClass * klass,
				gint nf, gint nv, gint nc, gint ns,
				BEM3DShapeFunc shf, BEM3DShapeFunc cpf)

{
  BEM3DElement *e ;
  gint i ;

  g_return_val_if_fail(klass != NULL, NULL) ;
  g_return_val_if_fail(nf > 0, NULL) ;
  g_return_val_if_fail(nv > 0, NULL) ;
  g_return_val_if_fail(nc >= 0, NULL) ;
  g_return_val_if_fail(ns > 0, NULL) ;
  g_return_val_if_fail(shf != NULL, NULL) ;
  g_return_val_if_fail(cpf != NULL, NULL) ;

  e = BEM3D_ELEMENT (gts_object_new (GTS_OBJECT_CLASS (klass)));

  e->nf = nf ;
  e->f = (gpointer *)g_malloc(nf*sizeof(gpointer)) ;
  e->nv = nv ;
  e->v = (gpointer *)g_malloc(nv*sizeof(gpointer)) ;
  e->ns = ns ;
  e->s = (gint *)g_malloc(ns*sizeof(gint)) ;
  e->xs = (gdouble *)g_malloc(2*nv*sizeof(gdouble)) ;
  bem3d_element_moment_order(e) = 0 ;
  e->Imn = NULL ;

  if ( nc > 0 ) {
    e->nc = nc ;
    e->c = (gpointer *)g_malloc(nc*sizeof(gpointer)) ;
    e->xc = (gdouble *)g_malloc(2*nc*sizeof(gdouble)) ;
  } else {
    e->nc = nv ; e->c = e->v ; e->xc = e->xs ;
  }

  e->i = (gint *)g_malloc((e->nc)*sizeof(gint)) ;
  for ( i = 0 ; i < e->nc ; i ++ ) e->i[i] = 0 ;
  e->shf = shf ; e->cpf = cpf ;

  e->reserved = NULL ;

  return e ;
}

/** 
 * Add a new face to a ::BEM3DElement. 
 * 
 * @param e pointer to element;
 * @param f pointer to face to add;
 * @param i local face index for new face.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_add_face(BEM3DElement *e, gpointer f, gint i)

{
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( i >= bem3d_element_face_number(e) || i < 0 )
    g_error("%s: cannot add face number %d to %d face element",
	    __FUNCTION__, i, bem3d_element_face_number(e)) ;

  e->f[i] = f ;

  return BEM3D_SUCCESS ;
}

/** 
 * Add a GtsVertex to a ::BEM3DElement.
 * 
 * @param e pointer to ::BEM3DElement;
 * @param v pointer to GtsVertex to add;
 * @param i local index of GtsVertex.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_add_vertex(BEM3DElement *e, gpointer v, gint i)

{
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( i >= bem3d_element_vertex_number(e) || i < 0 )
    g_error("%s: cannot add vertex number %d to %d node element",
	    __FUNCTION__, i, bem3d_element_vertex_number(e)) ;

  e->v[i] = v ;

  return BEM3D_SUCCESS ;
}

/** 
 * Add a collocation point to a ::BEM3DElement
 * 
 * @param e pointer to ::BEM3DElement;
 * @param v pointer to vertex to add;
 * @param i local index of collocation point.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_add_node(BEM3DElement *e, gpointer v, gint i)

{
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( i >= bem3d_element_node_number(e) )
    g_error("%s: cannot add collocation point number %d to %d "
	    "collocation point element",
	  __FUNCTION__, i, bem3d_element_node_number(e)) ;

  e->c[i] = v ;

  return BEM3D_SUCCESS ;
}

/** 
 * Set the global index of a collocation point.
 * 
 * @param e pointer to ::BEM3DElement;
 * @param i local index;
 * @param j global index.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_set_index(BEM3DElement *e, gint i, gint j)

{
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;

  if ( i >= bem3d_element_node_number(e) || i < 0 )
    g_error("%s: cannot set global index of %dth point of %d "
	    "collocation point element",
	    __FUNCTION__, i, bem3d_element_node_number(e)) ;
  
  e->i[i] = j ;

  return BEM3D_SUCCESS ;
}

/** 
 * Set a corner of a ::BEM3DElement.
 * 
 * @param e BEM3DElement to set;
 * @param i local index of corner vertex;
 * @param j index of corner to set.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_set_corner(BEM3DElement *e, gint i, gint j) 

{
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;

  if ( i >= bem3d_element_corner_number(e) || i < 0 )
    g_error("%s: cannot set %dth corner of %d corner element",
	    __FUNCTION__, i, bem3d_element_corner_number(e)) ;

  e->s[i] = j ;

  return BEM3D_SUCCESS ;
}

/** 
 * Write a ::BEM3DElement to file. Usually called from bem3d_mesh_write. 
 *
 * Format: nf nv nc ns shf cpf [faces] [vertices] [collocation points]
 *         [vertex area coordinate] [collocation point area coordinates]
 *         [collocation point global indices] [corner indices]
 * 
 * @param e element to write;
 * @param v hash table of vertices;
 * @param t hash table of faces;
 * @param f file pointer.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_write(BEM3DElement *e, 
			 GHashTable *v, GHashTable *t,
			 FILE *f)

{
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(t != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  fprintf(f, "%u %u", e->nf, e->nv) ;

  if ( e->v == e->c ) fprintf(f, " 0") ;
  else fprintf(f, " %u", e->nc) ;

  fprintf(f, " %u", e->ns) ;

  fprintf(f, " %s %s", 
	  bem3d_shapefunc_lookup_name(e->shf),
	  bem3d_shapefunc_lookup_name(e->cpf)) ;

  for ( i = 0 ; i < e->nf ; i ++ ) {
    fprintf(f, " %u",
	    GPOINTER_TO_UINT (g_hash_table_lookup (t, e->f[i]))) ;
  }

  for ( i = 0 ; i < e->nv ; i ++ ) {
    fprintf(f, " %u",
	    GPOINTER_TO_UINT (g_hash_table_lookup (v, e->v[i]))) ;
  }

  if ( e->v != e->c ) 
    for ( i = 0 ; i < e->nc ; i ++ ) {
      fprintf(f, " %u",
	      GPOINTER_TO_UINT (g_hash_table_lookup (v, e->c[i]))) ;
    }

  if ( e->v != e->c ) {
    for ( i = 0 ; i < 2*(e->nc) ; i ++ ) {
      fprintf(f, " %e", e->xc[i]) ;	    
    }
  }
  for ( i = 0 ; i < 2*(e->nv) ; i ++ ) {
    fprintf(f, " %e", e->xs[i]) ;	    
  }

  for ( i = 0 ; i < e->nc ; i ++ )  fprintf(f, " %u", e->i[i]) ;
  for ( i = 0 ; i < e->ns ; i ++ )  fprintf(f, " %u", e->s[i]) ;

  fprintf(f, "\n") ;

  return BEM3D_SUCCESS ;
}

/** 
 * Find the local index of a GtsVertex on a ::BEM3DElement. 
 * 
 * @param e ::BEM3DElement to search;
 * @param v GtsVertex to find.
 * 
 * @return local index of vertex if on element, -1 otherwise. 
 */

gint bem3d_element_find_vertex(BEM3DElement *e, GtsVertex *v)

{
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, BEM3D_NULL_ARGUMENT) ;

  for ( i = 0 ; i < e->nv ; i ++ ) if ( e->v[i] == v ) return i ;

  return BEM3D_FAILURE ;
}

/** 
 * Find the local index of a collocation point on a ::BEM3DElement.
 * 
 * @param e ::BEM3DElement to search;
 * @param v collocation point (GtsVertex) to find.
 * 
 * @return local index of collocation point if on element, -1 otherwise. 
 */

gint bem3d_element_find_node(BEM3DElement *e, GtsVertex *v)

{
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v), BEM3D_NULL_ARGUMENT) ;

  for ( i = 0 ; i < e->nc ; i ++ ) if ( e->c[i] == v ) return i ;

  return BEM3D_FAILURE ;
}

/** 
 * Surface area of a ::BEM3DElement, computed using a Gaussian
 * quadrature rule.
 * 
 * @param e ::BEM3DElement;
 * @param ngp number of quadrature points to use.
 * 
 * @return computed area of \a e.
 */

gdouble bem3d_element_area(BEM3DElement *e, gint ngp)

{
  gdouble J, dA ;
  static gdouble *L = NULL, *dLds = NULL, *dLdt = NULL ;
  BEM3DShapeFunc sf = bem3d_element_shape_func(e) ;
  BEM3DQuadratureRule *q ;
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(ngp > 0, BEM3D_NULL_ARGUMENT) ;

  if ( L == NULL ) {
    L = (gdouble *)g_malloc(bem3d_element_node_number(e)*4*sizeof(gdouble)) ;
    dLds = (gdouble *)g_malloc(bem3d_element_node_number(e)*4*sizeof(gdouble)) ;
    dLdt = (gdouble *)g_malloc(bem3d_element_node_number(e)*4*sizeof(gdouble)) ;
  }

  q = bem3d_quadrature_rule_new(ngp, 1) ;
  bem3d_quadrature_rule_wx(NULL, e, q, NULL, NULL, &ngp) ;

  dA = 0.0 ;
  for ( i = 0 ; i < bem3d_quadrature_vertex_number(q) ; i ++ ) {
    sf(bem3d_quadrature_xi(q,i),
       bem3d_quadrature_eta(q,i), L, dLds, dLdt, NULL) ;
    J = bem3d_element_jacobian(e, dLds, dLdt) ;
    dA += bem3d_quadrature_weight(q,i)*J ;
  }

  return dA ;
}

/** 
 * Interpolate position on a ::BEM3DElement given a set of shape
 * functions.
 * 
 * @param e ::BEM3DElement;
 * @param L array of shape functions for element;
 * @param q preallocated GtsPoint for interpolated position.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_position(BEM3DElement *e, gdouble *L,
			    GtsPoint *q)

{
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(L != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(q != NULL, BEM3D_NULL_ARGUMENT) ;
/*   g_return_val_if_fail(GTS_IS_POINT(q), BEM3D_ARGUMENT_WRONG_TYPE) ; */

  gts_point_set(GTS_POINT(q), 0.0, 0.0, 0.0) ;

  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) {
    GTS_POINT(q)->x += GTS_POINT(e->v[i])->x*L[i] ;
    GTS_POINT(q)->y += GTS_POINT(e->v[i])->y*L[i] ;
    GTS_POINT(q)->z += GTS_POINT(e->v[i])->z*L[i] ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Compute Jacobian on a ::BEM3DElement given derivatives of shape
 * functions.
 * 
 * @param e ::BEM3DElement;
 * @param dLds shape function derivatives \f$\partial L_{i}/\partial \xi\f$; 
 * @param dLdt shape function derivatives \f$\partial L_{i}/\partial \eta\f$.
 * 
 * @return Jacobian.
 */

gdouble bem3d_element_jacobian(BEM3DElement *e, gdouble *dLds, gdouble *dLdt)
  
{
  gint i ;
  GtsVector dxds, dxdt, normal ;
  gdouble J ;

  g_debug("%s: ", __FUNCTION__) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(dLds != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dLdt != NULL, BEM3D_NULL_ARGUMENT) ;

  dxdt[0] = dxdt[1] = dxdt[2] = 0.0 ;
  dxds[0] = dxds[1] = dxds[2] = 0.0 ;

  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) {
    dxds[0] += GTS_POINT(e->v[i])->x*dLds[i] ;
    dxds[1] += GTS_POINT(e->v[i])->y*dLds[i] ;
    dxds[2] += GTS_POINT(e->v[i])->z*dLds[i] ;
    dxdt[0] += GTS_POINT(e->v[i])->x*dLdt[i] ;
    dxdt[1] += GTS_POINT(e->v[i])->y*dLdt[i] ;
    dxdt[2] += GTS_POINT(e->v[i])->z*dLdt[i] ;
  }

  gts_vector_cross(normal,dxds,dxdt) ;

/*   g_assert(!isnan(normal[0])) ; g_assert(!isnan(normal[1])) ; */
/*   g_assert(!isnan(normal[2])) ; */

  J = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + 
	   normal[2]*normal[2]) ;

  return J ;
}

/** 
 * Compute normal to, and Jacobian on, a BEM3DElement given shape
 * function derivatives.
 * 
 * @param e ::BEM3DElement;
 * @param dLds shape function derivatives \f$\partial L_{i}/\partial \xi\f$; 
 * @param dLdt shape function derivatives \f$\partial L_{i}/\partial \eta\f$; 
 * @param normal GtsVector for normal;
 * @param J Jacobian.
 * 
 * @return ::BEM3D_SUCCESS on success;
 */

gint bem3d_element_normal(BEM3DElement *e, gdouble *dLds, gdouble *dLdt,
			  GtsVector normal, gdouble *J)

{
  gint i ;
  GtsVector dxds, dxdt ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(dLds != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dLdt != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(normal != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(J != NULL, BEM3D_NULL_ARGUMENT) ;

  dxdt[0] = dxdt[1] = dxdt[2] = 0.0 ;
  dxds[0] = dxds[1] = dxds[2] = 0.0 ;

  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) {
    dxds[0] += GTS_POINT(e->v[i])->x*dLds[i] ;
    dxds[1] += GTS_POINT(e->v[i])->y*dLds[i] ;
    dxds[2] += GTS_POINT(e->v[i])->z*dLds[i] ;
    dxdt[0] += GTS_POINT(e->v[i])->x*dLdt[i] ;
    dxdt[1] += GTS_POINT(e->v[i])->y*dLdt[i] ;
    dxdt[2] += GTS_POINT(e->v[i])->z*dLdt[i] ;
  }

  gts_vector_cross(normal,dxds,dxdt) ;

  if ( isnan(normal[0]) || isnan(normal[1]) || isnan(normal[2]) )
    g_error("%s: NaN in normal", __FUNCTION__) ;

  *J = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + 
	    normal[2]*normal[2]) ;

  g_assert(*J != 0.0) ;

  normal[0] /= (*J) ; normal[1] /= (*J) ; normal[2] /= (*J) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Assemble equations for ::BEM3DElement and GtsPoint.
 * 
 * The integration of the Green's function and its normal derivative
 * weighted by the element shape functions is performed for a
 * specified field point.
 *
 * @param e ::BEM3DElement;
 * @param x GtsPoint;
 * @param config ::BEM3DConfiguration configuration for this problem;
 * @param gdata ::BEM3DParameters to pass to Green's function;
 * @param G GArray of integral of \f$G L_{i}\f$ over element;
 * @param dGdn GArray of integral of \f$d G/d n L_{i}\f$ over element. 
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_assemble_equations(BEM3DElement *e, GtsPoint *x,
				      BEM3DConfiguration *config,
				      BEM3DParameters *gdata,
				      GArray *G, GArray *dGdn)
				    
{
  BEM3DShapeFunc shfunc = bem3d_element_shape_func(e) ;
  BEM3DShapeFunc cpfunc = bem3d_element_node_func(e) ;
  static BEM3DQuadratureRule *q = NULL ;
  static GArray *g = NULL ;
  static GArray *dgdn = NULL ;
  static gdouble *L = NULL, *dLds = NULL, *dLdt = NULL ;
  BEM3DQuadratureRuleFunc qf ;
  gpointer qd ;
  BEM3DGreensFunction gfunc ;
  gint i, j, k, stride ;
  gdouble s, t, J, w, wt ;
  GtsPoint y ;
  GtsVector n ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(x != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(G != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dGdn != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( g == NULL ) {
    g = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
    dgdn = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
    L = (gdouble *)g_malloc(4*bem3d_element_node_number(e)*sizeof(gdouble)) ;
    dLds = (gdouble *)g_malloc(4*bem3d_element_node_number(e)*sizeof(gdouble)) ;
    dLdt = (gdouble *)g_malloc(4*bem3d_element_node_number(e)*sizeof(gdouble)) ;
  }
  
  if ( q == NULL ) q = bem3d_quadrature_rule_new(0, 1) ;

  /* if ( qfunc == NULL ) qf = bem3d_quadrature_rule_default ; */
  /* else qf = qfunc ; */
  /* if ( qdata == NULL ) qd = bem3d_quadrature_selector_default() ; */
  /* else qd = qdata ; */

  gfunc = config->gfunc ;

  qf = config->qrule ; qd = config->qdata ;

  qf(x, e, q, &gfunc, gdata, qd) ;

  /* g_debug("%s: quadrature rule weights total: %g", */
  /* 	__FUNCTION__, bem3d_quadrature_rule_sum_weights(q)) ; */

  /*get the size of the Green's function*/
  if ( bem3d_greens_function_is_real(&gfunc) ) stride = 1 ;
  else stride = 2 ;
  g_array_set_size(G,bem3d_element_node_number(e)*stride) ; 
  g_array_set_size(dGdn,bem3d_element_node_number(e)*stride) ;  

  for ( i = 0 ; i < G->len ; i ++ )
    g_array_index(G,gdouble,i) = g_array_index(dGdn,gdouble,i) = 0.0 ;

  for ( i = 0 ; i < bem3d_quadrature_vertex_number(q) ; i ++ ) {
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
	g_array_index(G,gdouble,j*stride+k) +=
	  g_array_index(g,gdouble,k)*w ;
	g_array_index(dGdn,gdouble,j*stride+k) +=
	  g_array_index(dgdn,gdouble,k)*w ;
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
    for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
      g_array_index(G,gdouble,j*stride+k) += 
	bem3d_quadrature_free_term_g(q,j) ;
      g_array_index(dGdn,gdouble,j*stride+k) += 
	bem3d_quadrature_free_term_dg(q,j) ;
    }
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Check if a GtsVertex is on a ::BEM3DElement.
 * 
 * @param e ::BEM3DElement;
 * @param v GtsVertex.
 * 
 * @return TRUE if vertex is on element, FALSE otherwise. 
 */

gboolean bem3d_element_has_vertex(BEM3DElement *e, GtsVertex *v)

{
  gint i ;

  g_return_val_if_fail(e != NULL, FALSE) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, FALSE) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v), BEM3D_ARGUMENT_WRONG_TYPE) ;

  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) {
    if ( e->v[i] == v ) return TRUE ;
  }

  return FALSE ;
}

gint bem3d_element_replace_vertex(BEM3DElement *e, GtsVertex *v, GtsVertex *w)

{
  gint i ;
  
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(w != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v), BEM3D_ARGUMENT_WRONG_TYPE) ;

  i = bem3d_element_find_vertex(e, v) ;

  if ( i == BEM3D_FAILURE ) return BEM3D_SUCCESS ;

  e->c[i] = w ;

  return BEM3D_SUCCESS ;
}

gint bem3d_element_find_face(BEM3DElement *e, GtsFace *f)

{
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_FACE(f), BEM3D_ARGUMENT_WRONG_TYPE) ;

  for ( i = 0 ; i < e->nf ; i ++ ) if ( e->f[i] == f ) return i ;

  return BEM3D_FAILURE ;
}

gint bem3d_element_reset_index(BEM3DMesh *m, BEM3DElement *e, gint i, gint j)

{
#ifdef _FOREACH_USE_HASH_TABLE_
  gpointer k ;
#endif /*_FOREACH_USE_HASH_TABLE_*/

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;

  if ( i >= bem3d_element_node_number(e) || i < 0 )
    g_error("%s: cannot set global index of %dth point of %d "
	    "collocation point element",
	    __FUNCTION__, i, bem3d_element_node_number(e)) ;
  
  e->i[i] = j ;

#ifdef _FOREACH_USE_HASH_TABLE_
  k = g_hash_table_lookup(m->c,
			  GINT_TO_POINTER(bem3d_element_global_index(e,i)+1)
			  ) ;
  if ( k == NULL ) {
    g_hash_table_insert(m->c,
			GINT_TO_POINTER(bem3d_element_global_index(e,i)+1),
			e->c[i]) ;
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: inserting: %d; hash table size: %d", 
	  __FUNCTION__, bem3d_element_global_index(e,i),
	  g_hash_table_size(m->c)) ;
  } else {
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: not inserting: %d; vertex %p",
	  __FUNCTION__, bem3d_element_global_index(e,i),
	  e->c[i]) ;
  }
#else
  g_assert_not_reached() ;
#endif /*_FOREACH_USE_HASH_TABLE_*/


  return BEM3D_SUCCESS ;
}

gint bem3d_element_slopes(BEM3DElement *e, gdouble *dLds, gdouble *dLdt,
			  GtsVector dxds, GtsVector dxdt)

{
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(dLds != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dLdt != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dxds != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dxdt != NULL, BEM3D_NULL_ARGUMENT) ;

  dxdt[0] = dxdt[1] = dxdt[2] = 0.0 ;
  dxds[0] = dxds[1] = dxds[2] = 0.0 ;

  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) {
    dxds[0] += GTS_POINT(e->v[i])->x*dLds[i] ;
    dxds[1] += GTS_POINT(e->v[i])->y*dLds[i] ;
    dxds[2] += GTS_POINT(e->v[i])->z*dLds[i] ;
    dxdt[0] += GTS_POINT(e->v[i])->x*dLdt[i] ;
    dxdt[1] += GTS_POINT(e->v[i])->y*dLdt[i] ;
    dxdt[2] += GTS_POINT(e->v[i])->z*dLdt[i] ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Find collocation points shared between two ::BEM3DElement's. 
 * 
 * @param e1 ::BEM3DElement;
 * @param e2 ::BEM3DElement.
 * 
 * @return GSList of common collocation points.
 */

GSList *bem3d_element_common_nodes(BEM3DElement *e1, BEM3DElement *e2)

{
  GSList *c = NULL ;
  gint i ;

  g_return_val_if_fail(e1 != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e1), NULL) ;
  g_return_val_if_fail(e2 != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e2), NULL) ;

  for ( i = 0 ; i < bem3d_element_node_number(e1) ; i ++ ) {
    if ( bem3d_element_find_node(e2, bem3d_element_node(e1,i)) != -1 ) {
      c = g_slist_prepend(c, bem3d_element_node(e1,i)) ;
    }
  }

  return c ;
}

/** 
 * Find the elements of a ::BEM3DMesh which share an edge with a
 * specified element.
 * 
 * @param el ::BEM3DElement whose neighbours are sought;
 * @param m ::BEM3DMesh.
 * 
 * @return a GSList of unique elements which border the specified element. 
 */

GSList *bem3d_element_neighbours(BEM3DElement *el, BEM3DMesh *m)

{
  GSList *e = NULL, *f, *edges = NULL, *j ;
  BEM3DElement *g ;
  gint i ;

  g_return_val_if_fail(el != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(el), NULL) ;
  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;

  for ( i = 0 ; i < bem3d_element_face_number(el) ; i ++ ) {
    if ( !g_slist_find(edges, GTS_TRIANGLE(bem3d_element_face(el,i))->e1) )
      edges = g_slist_prepend(edges, 
			      GTS_TRIANGLE(bem3d_element_face(el,i))->e1) ;
    if ( !g_slist_find(edges, GTS_TRIANGLE(bem3d_element_face(el,i))->e2) )
      edges = g_slist_prepend(edges, 
			      GTS_TRIANGLE(bem3d_element_face(el,i))->e2) ;
    if ( !g_slist_find(edges, GTS_TRIANGLE(bem3d_element_face(el,i))->e3) )
      edges = g_slist_prepend(edges, 
			      GTS_TRIANGLE(bem3d_element_face(el,i))->e3) ;
  }

  for ( j = edges ; j != NULL ; j = j->next ) {
    for ( f = GTS_EDGE(j->data)->triangles ; f != NULL ; f = f->next ) {
      g = bem3d_mesh_face_element(m, GTS_FACE(f->data)) ;
      if ( !g_slist_find(e, g) && (g != el) ) {
	e = g_slist_prepend(e, g) ;
      }
    }
  }

  g_slist_free(edges) ;

  return e ;
}

/** 
 * Find an element containing a given collocation point, identified by
 * its GtsVertex and its global index.
 * 
 * @param m ::BEM3DMesh
 * @param v GtsVertex of collocation point
 * @param i global index of collocation point
 * 
 * @return ::BEM3DElement containing \a v with index \a i.
 */

BEM3DElement *bem3d_element_from_node(BEM3DMesh *m, GtsVertex *v, gint i)

{
  GSList *el ;
  BEM3DElement *e ;
  gint k ;
  
  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;
  g_return_val_if_fail(v != NULL, NULL) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v), NULL) ;

  el = bem3d_mesh_vertex_elements(m, v) ;
  e = NULL ;

  for ( ; el != NULL ; el = el->next ) {
    k = bem3d_element_find_node(BEM3D_ELEMENT(el->data), v) ;
    if ( bem3d_element_global_index(BEM3D_ELEMENT(el->data),k) == i )
      return (el->data) ;
  }

  return e ;
}

/** 
 * Index the collocation points of a ::BEM3DElement. If a point is
 * present in the ::BEM3DMesh \a m, the collocation point there is
 * used. If not, new indices are generated, starting from \a n, which
 * is incremented as required.
 * 
 * @param e ::BEM3DElement to index
 * @param m ::BEM3DMesh for index checking
 * @param n first new index to assign, incremented as new indices are
 * generated.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_index_nodes(BEM3DElement *e, BEM3DMesh *m, gint *n)

{
  gint i, j ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;  
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(n != NULL, BEM3D_NULL_ARGUMENT) ;  

  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    j = bem3d_mesh_index_from_node(m, bem3d_element_node(e,i)) ;
    if ( j == -1 ) {
      bem3d_element_set_index(e, i, (*n)) ; (*n) ++ ;
    } else {
      bem3d_element_set_index(e, i, j) ;
    }
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Generate two unit vectors, not necessarily orthogonal, which lie in
 * the surface at a point of a ::BEM3DElement.
 * 
 * @param e ::BEM3DElement
 * @param dLds derivative of geometric shape function at required
 * point of \a e
 * @param dLdt derivative of geometric shape function at required
 * point of \a e
 * @param u1 first local vector \f$d x/d\xi/h_{1}\f$
 * @param h1 scale factor \f$h_{1}=|d x/d\xi|\f$
 * @param u2 first local vector \f$d x/d\xi/h_{2}\f$
 * @param h2 scale factor \f$h_{2}=|d x/d\xi|\f$
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_element_local_vectors(BEM3DElement *e,
				 gdouble dLds[], gdouble dLdt[],
				 GtsVector u1, gdouble *h1,
				 GtsVector u2, gdouble *h2)
{
  gint i ;
  
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(dLds != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(dLdt != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(u1 != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(h1 != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(u2 != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(h2 != NULL, BEM3D_NULL_ARGUMENT) ;

  u1[0] = u1[1] = u1[2] = u2[0] = u2[1] = u2[2] = 0.0 ;

  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) {
    u1[0] += dLds[i]*GTS_POINT(bem3d_element_vertex(e,i))->x ;
    u1[1] += dLds[i]*GTS_POINT(bem3d_element_vertex(e,i))->y ;
    u1[2] += dLds[i]*GTS_POINT(bem3d_element_vertex(e,i))->z ;
    u2[0] += dLdt[i]*GTS_POINT(bem3d_element_vertex(e,i))->x ;
    u2[1] += dLdt[i]*GTS_POINT(bem3d_element_vertex(e,i))->y ;
    u2[2] += dLdt[i]*GTS_POINT(bem3d_element_vertex(e,i))->z ;
  }

  *h1 = gts_vector_norm(u1) ; *h2 = gts_vector_norm(u2) ;

  u1[0] /= (*h1) ; u1[1] /= (*h1) ; u1[2] /= (*h1) ;
  u2[0] /= (*h2) ; u2[1] /= (*h2) ; u2[2] /= (*h2) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Find the GtsEdge's common to two ::BEM3DElement's. 
 * 
 * @param e1 ::BEM3DElement
 * @param e2 ::BEM3DElement
 * 
 * @return a GSList of GtsEdges common to \a e1 and \a e2, or NULL if
 * the elements do not touch.
 */

GSList *bem3d_elements_common_edges(BEM3DElement *e1, BEM3DElement *e2)

{
  GSList *e = NULL ;
  GtsEdge *edge ;
  gint i, j ;

  g_return_val_if_fail(e1 != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e1), NULL) ;
  g_return_val_if_fail(e2 != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e2), NULL) ;

  for ( i = 0 ; i < bem3d_element_face_number(e1) ; i ++ ) {
    for ( j = 0 ; j < bem3d_element_face_number(e2) ; j ++ ) {
      if ( (edge = 
	    gts_triangles_common_edge(GTS_TRIANGLE(bem3d_element_face(e1,i)),
				      GTS_TRIANGLE(bem3d_element_face(e2,j)))
	    ) != NULL ) {
	e = g_slist_prepend(e, edge) ;
      }
    }
  }

  return e ;
}

/** 
 * Calculate the weights for computing the gradient of a function on a
 * ::BEM3DElement. The gradient \f$\partial\phi/\partial
 * x=\sum_{i=0}^{N}\phi_{i}w_{3i}\f$, \f$\partial\phi/\partial
 * y=\sum_{i=0}^{N}\phi_{i}w_{3i+1}\f$, \f$\partial\phi/\partial
 * z=\sum_{i=0}^{N}\phi_{i}w_{3i+2}\f$, where there are \f$N\f$
 * collocation points on the element \a e. The weights are computed by
 * differentiation of the element shape functions.
 * 
 * @param e ::BEM3DElement
 * @param xi local coordinate
 * @param eta local coordinate
 * @param w array of weights, which should be presized to at least
 * \f$3N\f$ entries.
 * 
 * @return ::BEM3D_SUCCESS on success. 
 */

gint bem3d_element_gradient_weights(BEM3DElement *e, gdouble xi, gdouble eta,
				    gdouble *w)

{
  gint i ;
  BEM3DShapeFunc shfunc ;
  BEM3DMatrix J, iJ ;
  gdouble L[32], dLds[32], dLdt[32] ;
  GtsVector u1, u2, u3 ;
  gdouble h1, h2 ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(w != NULL, BEM3D_NULL_ARGUMENT) ;

  shfunc = bem3d_element_shape_func(e) ;

  shfunc(xi, eta, L, dLds, dLdt, NULL) ;
  bem3d_element_local_vectors(e, dLds, dLdt, u1, &h1, u2, &h2) ;
  gts_vector_cross(u3,u1,u2) ; 
  u1[0] *= h1 ; u1[1] *= h1 ; u1[2] *= h1 ;
  u2[0] *= h2 ; u2[1] *= h2 ; u2[2] *= h2 ;
  u3[0] *= h1*h2 ; u3[1] *= h1*h2 ; u3[2] *= h1*h2 ;

  J[0] = u1[0] ; J[1] = u1[1] ; J[2] = u1[2] ; 
  J[3] = u2[0] ; J[4] = u2[1] ; J[5] = u2[2] ; 
  J[6] = u3[0] ; J[7] = u3[1] ; J[8] = u3[2] ; 
  bem3d_matrix_inverse(J, iJ) ;
  
  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    u1[0] = dLds[i] ; u1[1] = dLdt[i] ; u1[2] = 0.0 ;
    bem3d_matrix_vector_mul(iJ, u1, u2) ;
    w[i*3+0] = u2[0] ; w[i*3+1] = u2[1] ; w[i*3+2] = u2[2] ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Check if a vertex is a corner of an element. 
 * 
 * @param e a ::BEM3DElement;
 * @param v a GtsVertex of \a e.
 * 
 * @return the index of a corner of \a e if \a v is a corner, -1
 * otherwise.
 */

gint bem3d_element_vertex_is_corner(BEM3DElement *e, GtsVertex *v)

{
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v), BEM3D_ARGUMENT_WRONG_TYPE) ;

  for ( i = 0 ; i < bem3d_element_corner_number(e) ; i ++ ) 
    if ( bem3d_element_corner(e,i) == v ) return i ;

  return -1 ;
}

/** 
 * Locate a node of given global index on a ::BEM3DElement
 * 
 * @param e a ::BEM3DElement;
 * @param i global index of a node.
 * 
 * @return local index of a node of \a e with global index \a i, if
 * found, -1 otherwise.
 */

gint bem3d_element_find_index(BEM3DElement *e, gint i)

{
  gint j ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(i >= 0, BEM3D_ARGUMENT_OUT_OF_RANGE) ;

  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) 
    if ( bem3d_element_global_index(e, j) == i ) return j ;

  return -1 ;
}

/** 
 * Generate a ::GtsBBox containing the vertices of a ::BEM3DElement. 
 * 
 * @param klass a GtsBBoxClass;
 * @param e a ::BEM3DElement;
 * 
 * @return a ::GtsBBox containing the vertices of \a e or NULL on
 * error.
 */

GtsBBox *bem3d_element_bounding_box(GtsBBoxClass *klass,
				    BEM3DElement *e)

{
  GtsBBox *b ;
  gint i ;
  gdouble x1, y1, z1, x2, y2, z2 ;

  g_return_val_if_fail(klass != NULL, NULL) ;
  g_return_val_if_fail(e != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), NULL) ;

  x1 = y1 = z1 =  G_MAXDOUBLE ;
  x2 = y2 = z2 = -G_MAXDOUBLE ;
 
  for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) {
    x1 = MIN(x1, GTS_POINT(bem3d_element_vertex(e,i))->x) ;
    y1 = MIN(y1, GTS_POINT(bem3d_element_vertex(e,i))->y) ;
    z1 = MIN(z1, GTS_POINT(bem3d_element_vertex(e,i))->z) ;
    x2 = MAX(x2, GTS_POINT(bem3d_element_vertex(e,i))->x) ;
    y2 = MAX(y2, GTS_POINT(bem3d_element_vertex(e,i))->y) ;
    z2 = MAX(z2, GTS_POINT(bem3d_element_vertex(e,i))->z) ;
  }

  b = gts_bbox_new(klass, e, x1, y1, z1, x2, y2, z2) ;

  return b ;
}

gint bem3d_element_moments_make(BEM3DElement *e, gint H)

{
  gdouble x1[3], x2[3], x3[3] ;
  GtsVector s, t, n, r ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(H >= 0, BEM3D_ARGUMENT_OUT_OF_RANGE) ;
  
  if ( e->Imn == NULL ) e->Imn = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  g_array_set_size(e->Imn,(H+1)*(H+2)/2) ;

  if ( bem3d_element_shape_func(e) == bem3d_shfunc_t1 ) {
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
      newman_tri_moments(x1, x2, x3, H, 
			 (gdouble *)(&bem3d_element_moment(e,0,0))) ;
  } else 
    g_error("%s: only implemented for linear triangular elements",
	    __FUNCTION__) ;

  bem3d_element_moment_order(e) = H ;

  return BEM3D_SUCCESS ;
}

/** 
 * Find the ::BEM3DElement s shared by two GtsVertex s. 
 * 
 * @param m a ::BEM3DMesh;
 * @param v1 a vertex of \a m;
 * @param v2 another vertex of \a m;
 * 
 * @return a GSList of ::BEM3DElement s of \a m which contain both \a
 * v1 and \a v2.
 */

GSList *bem3d_elements_from_vertices(BEM3DMesh *m, GtsVertex *v1, 
				     GtsVertex *v2)

{
  GSList *e = NULL, *e1, *i ;

  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;
  g_return_val_if_fail(v1 != NULL, NULL) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v1), NULL) ;
  g_return_val_if_fail(v2 != NULL, NULL) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v2), NULL) ;  
  g_return_val_if_fail(v1 != v2, NULL) ;  

  e1 = bem3d_mesh_vertex_elements(m, v1) ;
  for ( i = e1 ; i != NULL ; i = i->next ) 
    if ( bem3d_element_find_vertex(BEM3D_ELEMENT(i->data), v2) != -1 )
      e = g_slist_prepend(e, i->data) ;
  
  return e ;
}

/**
 * @}
 * 
 */

