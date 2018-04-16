/* bmesh.c
 * 
 * Copyright (C) 2006, 2009 Michael Carley
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
 * @defgroup mesh BEM3D meshes
 * @{
 * 
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

#define BMESH_FOREACH_DATA_WIDTH 8
#define BMESH_FOREACH_FUNC       0
#define BMESH_FOREACH_DATA       1
#define BMESH_FOREACH_HASH       2
#define BMESH_FOREACH_MESH       3

static void destroy_element(gpointer key, BEM3DElement *e, gpointer data)

{
  gts_object_destroy(GTS_OBJECT(e)) ;

  return ;
}

static void mesh_destroy(GtsObject *object)

{
  BEM3DMesh *m = BEM3D_MESH(object) ;

  if ( m->e != NULL ) {
    g_hash_table_foreach(m->e, (GHFunc)destroy_element, NULL) ;
    g_hash_table_destroy(m->e) ;
  }
  if ( m->f != NULL ) g_hash_table_destroy(m->f) ;
  if ( m->c != NULL ) g_hash_table_destroy(m->c) ;

  (*GTS_OBJECT_CLASS(bem3d_mesh_class ())->parent_class->destroy)
    (object) ;

  return ;
}

static void bem3d_mesh_class_init (BEM3DMeshClass * klass)
{
  /* define new methods and overload inherited methods here */
  GTS_OBJECT_CLASS (klass)->destroy = mesh_destroy ;
}

static void bem3d_mesh_init (BEM3DMesh * object)
{
  /* initialize object here */
}

/** 
 * Definition of the ::BEM3DMeshClass. 
 * 
 * @return ::BEM3DMeshClass
 */

BEM3DMeshClass * bem3d_mesh_class (void)
{
  static BEM3DMeshClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo bem3d_mesh_info = {
      "BEM3DMesh",
      sizeof (BEM3DMesh),
      sizeof (BEM3DMeshClass),
      (GtsObjectClassInitFunc) bem3d_mesh_class_init,
      (GtsObjectInitFunc) bem3d_mesh_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_surface_class ()),
				  &bem3d_mesh_info);
  }

  return klass;
}

/** 
 * Create a new BEM3DMesh.
 * 
 * @param klass the BEM3D mesh class for the mesh;
 * @param face_class the GTS face class;
 * @param edge_class the GTS edge class;
 * @param vertex_class the GTS vertex class. 
 * 
 * @return a pointer to the newly created mesh
 */

BEM3DMesh *bem3d_mesh_new(BEM3DMeshClass *klass,
			  GtsFaceClass *face_class,
			  GtsEdgeClass *edge_class,
			  GtsVertexClass *vertex_class)

{
  BEM3DMesh *m ;

  g_return_val_if_fail(klass != NULL, NULL) ;
  g_return_val_if_fail(face_class != NULL, NULL) ;
  g_return_val_if_fail(edge_class != NULL, NULL) ;
  g_return_val_if_fail(vertex_class != NULL, NULL) ;
  
  m = BEM3D_MESH (gts_object_new (GTS_OBJECT_CLASS (klass)));

  GTS_SURFACE(m)->face_class = face_class ;
  GTS_SURFACE(m)->edge_class = edge_class ;
  GTS_SURFACE(m)->vertex_class = vertex_class ;

  m->e = g_hash_table_new(g_direct_hash, g_direct_equal) ;
  m->f = g_hash_table_new(g_direct_hash, g_direct_equal) ;
  m->c = g_hash_table_new(g_direct_hash, g_direct_equal) ;
#ifndef _FOREACH_USE_HASH_TABLE_
  m->cpt = g_byte_array_new() ;
#endif /*_FOREACH_USE_HASH_TABLE_*/

  m->i0 = 0 ; m->i1 = G_MAXINT ;

  return m ;
}

/** 
 * Add an element to a BEM3DMesh, checking to see if it is already
 * present.
 * 
 * @param m BEM3DMesh to add element to
 * @param e element to add
 * @param force if TRUE force addition of the element vertices, even 
 * if \a e is present in \a m
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_mesh_add_element(BEM3DMesh *m, BEM3DElement *e, gboolean force)

{
  gint i ;
  gpointer k ;
#ifndef _FOREACH_USE_HASH_TABLE_
  guint8 f = 0 ;
#endif /*_FOREACH_USE_HASH_TABLE_*/

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;

  g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	"%s: inserting element: %p", __FUNCTION__, e) ;

  if ( g_hash_table_lookup(m->e, e) != NULL && !force ) {
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: element %p already in mesh %p", __FUNCTION__, e, m) ;
    return BEM3D_SUCCESS ;
  }

  if ( g_hash_table_lookup(m->e, e) == NULL )
    g_hash_table_insert(m->e, e, e) ;

  for ( i = 0 ; i < e->nf ; i ++ ) {
    g_hash_table_insert(m->f, e->f[i], e) ;
    gts_surface_add_face(GTS_SURFACE(m), e->f[i]) ;
  }

#ifdef _FOREACH_USE_HASH_TABLE_
  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
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
  }  
#else
  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    k = g_hash_table_lookup(m->c,
			    GINT_TO_POINTER(bem3d_element_global_index(e,i)+1)
			    ) ;
    if ( k == NULL ) {
      g_hash_table_insert(m->c,
			  GINT_TO_POINTER(bem3d_element_global_index(e,i)+1),
			  GINT_TO_POINTER(m->cpt->len)) ;
      m->cpt = g_byte_array_append(m->cpt,&f,1) ;
    }
  }
#endif /*_FOREACH_USE_HASH_TABLE_*/
  return BEM3D_SUCCESS ;
}

/** 
 * Extract a linked list of the elements of a mesh containing 
 * a given vertex.
 * 
 * @param m BEM3DMesh containing vertex;
 * @param v vertex to check.
 * 
 * @return GSList of elements connected to v.
 */

GSList *bem3d_mesh_vertex_elements(BEM3DMesh *m, GtsVertex *v)

{
  GSList *f, *e, *i ;
  BEM3DElement *el ;

  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;
  g_return_val_if_fail(v != NULL, NULL) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v), NULL) ;

  f = gts_vertex_faces(v, GTS_SURFACE(m), NULL) ;
  
  e = NULL ;
  for ( i = f ; i != NULL ; i = i->next ) {
    el = g_hash_table_lookup(m->f, i->data) ;
    if ( !g_slist_find(e, el) ) e = g_slist_prepend(e, el) ;
  }

  g_slist_free(f) ;

  return e ;
}

static void reset_indices(BEM3DElement *e, gpointer data)

{
  gint i ;

  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ )
    bem3d_element_global_index(e,i) = -1 ;

  return ;
}

static void element_normal(BEM3DMesh *m, BEM3DElement *e, gint i, GtsVector n)

{
  gdouble L[32], dLds[32], dLdt[32], J ;
  
  (e->shf)(bem3d_element_node_xi(e,i), bem3d_element_node_eta(e,i),
	   L, dLds, dLdt, NULL) ;
  bem3d_element_normal(e, dLds, dLdt, n, &J) ;

  return ;
}

static gdouble vector_angle(GtsVector n, GtsVector m)

{
  gdouble a = gts_vector_scalar(n,m)/gts_vector_norm(n)/gts_vector_norm(m) ;
  if ( a >  1.0 ) return 0.0 ;
  if ( a < -1.0 ) return M_PI ;
  return acos(a) ;
}

static gint node_is_indexed(GtsVertex *v, BEM3DMesh *m, GtsVector n,
			    gdouble C)

{
  GSList *e, *i ;
  gint j ;
  BEM3DElement *el ;
  GtsVector nn ;
  gdouble va ;

  e = bem3d_mesh_vertex_elements(m, v) ;

  for ( i = e ; i != NULL ; i = i->next ) {
    el = BEM3D_ELEMENT(i->data) ;
    j = bem3d_element_find_node(el, v) ;
    element_normal(m, el, j, nn) ;
    va = vector_angle(n,nn) ;
    if ( (bem3d_element_global_index(el,j) != -1) && va <= C )
/* 	 fabs(vector_angle(n,nn)) <= fabs(C) ) */
/*     if ( GTS_POINT(v)->x == 2.0 ) { */
/*       fprintf(stderr, "%lg %lg %lg %lg %lg %lg %lg %lg\n", */
/* 	      GTS_POINT(v)->x, GTS_POINT(v)->y, */
/* 	      GTS_POINT(v)->z, nn[0], nn[1], nn[2], */
/* 	      gts_vector_scalar(n,nn), C */
/* 	      ) ; */
/*     } */
/*     if ( (bem3d_element_global_index(el,j) != -1) && */
/* 	 gts_vector_scalar(n,nn) >= C ) */
      return
	bem3d_element_global_index(el,j) ;
  }
  
  return -1 ;
}

static void set_indices(BEM3DElement *e, gpointer data[])

{
  BEM3DMesh *m = data[BEM3D_BMESH_DATA_MESH] ;
  gdouble angle = *(gdouble *)data[BEM3D_BMESH_DATA_ANGLE] ;
  gint *n = data[BEM3D_BMESH_DATA_N] ;
  GtsVector norm ;
  gint i, j ;

  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    element_normal(m, e, i, norm) ;
    if ( (j = node_is_indexed(bem3d_element_node(e,i), 
			      m, norm, angle)) != -1 ) {
/* 			      m, norm, cos(angle))) != -1 ) { */
      bem3d_element_global_index(e,i) = j ;
    } else {
      bem3d_element_global_index(e,i) = (*n) ; (*n) ++ ;
    }
  }

  return ;
}

static void add_element(BEM3DElement *e, gpointer data[])

{
  BEM3DMesh *m = data[BEM3D_BMESH_DATA_MESH] ;

  bem3d_mesh_add_element(m, e, TRUE) ;

  return ;
}

/** 
 * Number the collocation points of a mesh, with sharp corners double
 * indexed if the angle between two element normals at an edge is
 * greater than \a angle.
 * 
 * @param m ::BEM3DMesh to be numbered;
 * @param angle angle to be used in determining edge sharpness for 
 * double indexing;
 * @param ni initial index.
 * 
 * @return number of indices on success, which can be used as initial
 * index for a subsequent mesh.
 */

gint bem3d_mesh_index_nodes(BEM3DMesh *m, gdouble angle, gint ni)

{
  gint n ;
  gpointer data[BEM3D_BMESH_DATA_WIDTH] ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

  n = ni ;
  data[BEM3D_BMESH_DATA_MESH] = m ;
  data[BEM3D_BMESH_DATA_N] = &n ;
  data[BEM3D_BMESH_DATA_ANGLE] = &angle ;

  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)reset_indices, NULL) ;
  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)set_indices, data) ;
  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)add_element, data) ;
  
  return n ;
}

static void bem3d_edge_discretize(GtsEdge *e, gpointer *data)

{
  gint *n, nne, i ;
  gdouble t ;
  GtsEdge **es ;
  GHashTable *eh ;
  GtsVertex *v[8] ;

  n = (gint *)data[BEM3D_BMESH_DATA_N] ;
  es = (GtsEdge **)data[BEM3D_BMESH_DATA_EDGES] ;
  nne = *(gint *)data[BEM3D_BMESH_DATA_N_EDGES] ;
  eh = (GHashTable *)data[BEM3D_BMESH_DATA_HASH] ;

  v[0] = GTS_SEGMENT(e)->v1 ;
  v[nne] = GTS_SEGMENT(e)->v2 ;

  for ( i = 1 ; i < nne ; i ++ ) {
    t = (gdouble)i/(gdouble)(nne) ;
    v[i] = gts_vertex_new(gts_vertex_class(),
			  GTS_POINT(v[0])->x*(1.0-t)+
			  GTS_POINT(v[nne])->x*t,
			  GTS_POINT(v[0])->y*(1.0-t)+
			  GTS_POINT(v[nne])->y*t,
			  GTS_POINT(v[0])->z*(1.0-t)+
			  GTS_POINT(v[nne])->z*t) ;
  }

  for ( i = 0 ; i < nne ; i ++ )
    es[(*n)*nne+i] = gts_edge_new(gts_edge_class(), v[i], v[i+1]) ;

  g_hash_table_insert(eh, e, 
		      GUINT_TO_POINTER (*((guint *) n))) ;

  (*n)++ ;

  return ;
}

static void edge_split_order(GtsVertex *v1, GtsEdge **es,
			     gint nne, GtsEdge **enew)

     /*
       extract the components of a split edge in the correct order for
       generation of an element
     */

{
  gint i, j ;

  for ( i = 0 ; 
	(i < nne)
 	  && (GTS_SEGMENT(es[i])->v1 != v1) 
	  && (GTS_SEGMENT(es[i])->v2 != v1) ;
	i ++ ) ;

  if ( i == 0 ) {
    for ( j = 0 ; j < nne ; j ++ ) enew[j] = es[j] ;
  } else {
    for ( j = 0 ; j < nne ; j ++ ) enew[j] = es[nne-1-j] ;
  }

  return ;
}

static void bem3d_face_discretize(GtsFace *f, gpointer *data)

{
  gint n, nne, i, j ;
  GtsEdge **es, *et[3], *enew[16] ;
  GHashTable *eh ;
  GtsVertex *v[3], *vnew[16] ;
  BEM3DMesh *m ; 
  BEM3DElement *el ;
  BEM3DElementBuildFunc bfunc ;

  es = (GtsEdge **)data[BEM3D_BMESH_DATA_EDGES] ;
  nne = *(gint *)data[BEM3D_BMESH_DATA_N_EDGES] ;
  eh = (GHashTable *)data[BEM3D_BMESH_DATA_HASH] ;
  m = (BEM3DMesh *)data[BEM3D_BMESH_DATA_MESH] ;
  bfunc = (BEM3DElementBuildFunc)data[BEM3D_BMESH_DATA_ELEMENT] ;

  gts_triangle_vertices_edges(GTS_TRIANGLE(f), NULL,
			      &v[0], &v[1], &v[2], 
			      &et[0], &et[1], &et[2]) ;

  /*extract the new edges in order*/
  for ( i = 0 ; i < 3 ; i ++ ) {
    n = GPOINTER_TO_UINT(g_hash_table_lookup(eh, et[i])) ;
    edge_split_order(v[i], &(es[n*nne]), nne, &(enew[i*nne])) ;
  }

  /*extract the new vertices in order*/
  for ( i = 0 ; i < 16 ; i ++ ) vnew[i] = NULL ;
  vnew[0] = v[0] ;
  for ( i = 0 ; i < 3 ; i ++ ) {
    for ( j = 0 ; j < nne ; j ++ ) {
      if ( GTS_SEGMENT(enew[i*nne+j])->v1 != vnew[i*nne+j] )
	vnew[i*nne+j+1] = GTS_SEGMENT(enew[i*nne+j])->v1 ;
      else
	vnew[i*nne+j+1] = GTS_SEGMENT(enew[i*nne+j])->v2 ;
    }
  }
  /*for build functions that depend on having the element vertices
    once and once only*/
  vnew[3*nne] = NULL ; 

  el = bfunc(enew, vnew) ;
  bem3d_mesh_add_element(m, el, FALSE) ;
  
  return ;
}

/** 
 * Discretize a surface by splitting its triangles into higher order
 * elements.
 * 
 * @param s GtsSurface to be discretized;
 * @param nne number of nodes per edge of each triangular face;
 * @param bfunc ::BEM3DElementBuildFunc for the required element;
 * @param m ::BEM3DMesh to be generated.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_mesh_discretize(GtsSurface *s, gint nne,
			   BEM3DElementBuildFunc bfunc,
			   BEM3DMesh *m)

{
  gpointer data[BEM3D_BMESH_DATA_WIDTH] ;
  gint n, ne ;
  GtsEdge **es ;
  GHashTable *e ;

  g_return_val_if_fail(s != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_SURFACE(s), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(bfunc != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

  ne = gts_surface_edge_number(s) ;
  es = (GtsEdge **)g_malloc(ne*nne*sizeof(GtsEdge *)) ;
  e = g_hash_table_new(g_direct_hash, g_direct_equal) ;

  data[BEM3D_BMESH_DATA_MESH] = m ;
  data[BEM3D_BMESH_DATA_N] = &n ;
  data[BEM3D_BMESH_DATA_EDGES] = es ;
  data[BEM3D_BMESH_DATA_N_EDGES] = &nne ;
  data[BEM3D_BMESH_DATA_HASH] = e ;
  data[BEM3D_BMESH_DATA_ELEMENT] = (gpointer)bfunc ;

  n = 0 ;
  gts_surface_foreach_edge(s, (GtsFunc)bem3d_edge_discretize, data) ;
  g_assert(n == ne) ;

  gts_surface_foreach_face(s, (GtsFunc)bem3d_face_discretize, data) ;

  g_hash_table_destroy(e) ;

  return BEM3D_SUCCESS ;
}

static void count_nodes(gint i, GtsVertex *v, gint *n)

{
  (*n)++ ;

  return ;
}

/** 
 * Count the number of collocation points on a mesh
 * 
 * @param m ::BEM3DMesh 
 * 
 * @return number of nodes on \a m.
 */

gint bem3d_mesh_node_number(BEM3DMesh *m)

{  
  gint n = 0 ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)count_nodes, &n) ;

  return n ;
}

static void bem3d_element_subassemble(BEM3DElement *e, gpointer *data)

{
  BEM3DConfiguration *config = data[BEM3D_BMESH_DATA_CONFIG] ; 
  BEM3DParameters *gdata = data[BEM3D_BMESH_DATA_GDATA] ;
  gpointer edata = data[BEM3D_BMESH_DATA_EDATA] ;
  BEM3DEquationFunc efunc = (BEM3DEquationFunc)data[BEM3D_BMESH_DATA_EFUNC] ;
  GtsVertex *x = (GtsVertex *)data[BEM3D_BMESH_DATA_NODE] ;
  gint row = *(gint *)data[BEM3D_BMESH_DATA_ROW] ;
  static GArray *G = NULL ;
  static GArray *dGdn = NULL ;
  gint i, stride ;

  if ( G == NULL ) {
    G = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
    dGdn = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  }

  g_array_set_size(G,bem3d_element_node_number(e)) ; 
  g_array_set_size(dGdn,bem3d_element_node_number(e)) ;
  bem3d_element_assemble_equations(e, GTS_POINT(x), 
				   config, gdata,
				   G, dGdn) ;
  stride = (G->len)/bem3d_element_node_number(e) ;
  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ )
    efunc(row, bem3d_element_global_index(e, i),
	  &(g_array_index(G,gdouble,stride*i)),
	  &(g_array_index(dGdn,gdouble,stride*i)),
	  stride, edata) ;
	  
  return ;
}

static void bem3d_node_subassemble(gint i, GtsVertex *v, gpointer *data)

{
  BEM3DMesh *m = (BEM3DMesh *)data[BEM3D_BMESH_DATA_MESH] ;
  gint imin = *(gint *)data[BEM3D_BMESH_DATA_IMIN] ;
  gint imax = *(gint *)data[BEM3D_BMESH_DATA_IMAX] ;
  data[BEM3D_BMESH_DATA_ROW] = &i ;
  data[BEM3D_BMESH_DATA_NODE] = v ;

  if ( (i >= imin) && ( i < imax ) ) {
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, "%s: i=%d", __FUNCTION__, i) ;
    bem3d_mesh_foreach_element(m, 
			     (BEM3DElementFunc)bem3d_element_subassemble, 
			     data) ;
  }

  return ;
}

/** 
 * Assemble the BEM3D equations for two meshes (which may be the same 
 * mesh)
 * 
 * @param m ::BEM3DMesh whose elements are to be processed;
 * @param n ::BEM3DMesh whose nodes are to be processed;
 * @param config ::BEM3DConfiguration for the problem;
 * @param gdata ::BEM3DParameters for the Green's function;
 * @param efunc ::BEM3DEquationFunc for the problem
 * @param edata data to pass to efunc
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_mesh_assemble_equations(BEM3DMesh *m, BEM3DMesh *n,
				   BEM3DConfiguration *config,
				   BEM3DParameters *gdata,
				   BEM3DEquationFunc efunc, gpointer edata)
{
  gpointer data[BEM3D_BMESH_DATA_WIDTH] ;
  gint np ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(n != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(n), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(gdata != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(efunc != NULL, BEM3D_NULL_ARGUMENT) ;

  np = bem3d_mesh_node_number(m) ;

  data[BEM3D_BMESH_DATA_CONFIG] = config ;
  data[BEM3D_BMESH_DATA_MESH] = m ; 
  data[BEM3D_BMESH_DATA_GDATA] = gdata ;
  data[BEM3D_BMESH_DATA_EFUNC] = (gpointer)efunc ;
  data[BEM3D_BMESH_DATA_EDATA] = edata ;
  data[BEM3D_BMESH_DATA_N] = &np ;
  data[BEM3D_BMESH_DATA_IMIN] = &(bem3d_mesh_node_index_min(m)) ;
  data[BEM3D_BMESH_DATA_IMAX] = &(bem3d_mesh_node_index_max(m)) ;

  g_log(G_LOG_DOMAIN, G_LOG_LEVEL_INFO, "%s", __FUNCTION__) ;

  bem3d_mesh_foreach_node(n, (BEM3DNodeFunc)bem3d_node_subassemble, data) ;

  return BEM3D_SUCCESS ;
}

static void bem3d_foreach_element(gpointer key, BEM3DElement *e, gpointer *fd)

{
  BEM3DElementFunc f = (BEM3DElementFunc)fd[BMESH_FOREACH_FUNC] ;
  gpointer data = fd[BMESH_FOREACH_DATA] ;

  f(e, data) ;

  return ;
}

/** 
 * Visit each element of a mesh and execute a function
 * 
 * @param m mesh to process
 * @param f BEM3DElementFunc to call
 * @param data user data to pass to f
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_mesh_foreach_element(BEM3DMesh *m, BEM3DElementFunc f, gpointer data)

{
  gpointer fd[BMESH_FOREACH_DATA_WIDTH] ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  fd[BMESH_FOREACH_FUNC] = (gpointer)f ;
  fd[BMESH_FOREACH_DATA] = data ;

  g_hash_table_foreach(m->e, (GHFunc)bem3d_foreach_element, fd) ;

  return BEM3D_SUCCESS ;
}

static void bem3d_foreach_node(gpointer key, BEM3DElement *e, gpointer *fd)

{
  BEM3DMesh *m = (BEM3DMesh *)fd[BMESH_FOREACH_MESH] ;
  BEM3DNodeFunc func = (BEM3DNodeFunc)fd[BMESH_FOREACH_FUNC] ;
  gpointer data = fd[BMESH_FOREACH_DATA] ;
#ifdef _FOREACH_USE_HASH_TABLE_
  GHashTable *h = (GHashTable *)fd[BMESH_FOREACH_HASH] ;
#else
  gint j ;
#endif /*_FOREACH_USE_HASH_TABLE_*/
  gint i ;

  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    if ( (e->i[i] >= bem3d_mesh_node_index_min(m)) &&
	 (e->i[i] < bem3d_mesh_node_index_max(m)) ) {
#ifdef _FOREACH_USE_HASH_TABLE_
/*       fprintf(stderr, "index: %d\n", e->i[i]) ; */
      if ( g_hash_table_lookup(h, GINT_TO_POINTER((e->i[i]+1))) == NULL ) {
	func(e->i[i], e->c[i], data) ;
	g_hash_table_insert(h, 
			    GINT_TO_POINTER(e->i[i]+1), 
			    GINT_TO_POINTER(e->i[i]+1)) ;
      }
#else
      g_assert_not_reached() ;
#endif /*_FOREACH_USE_HASH_TABLE_*/
    }
  }

  return ;
}

/** 
 * Visit each collocation point of a mesh and execute a function
 * 
 * @param m mesh to process
 * @param f BEM3DNodeFunc to call
 * @param data user data to pass to f
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_mesh_foreach_node(BEM3DMesh *m, BEM3DNodeFunc f, 
			     gpointer data)

{
  gpointer fd[BMESH_FOREACH_DATA_WIDTH] ;
#ifdef _FOREACH_USE_HASH_TABLE_
  GHashTable *h ;
  h = g_hash_table_new(NULL, NULL) ;
  fd[BMESH_FOREACH_HASH] = h ;
#endif /*_FOREACH_USE_HASH_TABLE_*/

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  fd[BMESH_FOREACH_FUNC] = (gpointer)f ;
  fd[BMESH_FOREACH_DATA] = data ;
  fd[BMESH_FOREACH_MESH] = m ;

#ifndef _FOREACH_USE_HASH_TABLE_
  memset(m->cpt->data, 0, m->cpt->len) ;
#endif /*_FOREACH_USE_HASH_TABLE_*/
  g_hash_table_foreach(m->e, (GHFunc)bem3d_foreach_node, fd) ;

#ifdef _FOREACH_USE_HASH_TABLE_
  g_hash_table_destroy(h) ;
#endif /*_FOREACH_USE_HASH_TABLE_*/

  return BEM3D_SUCCESS ;
}

gint bem3d_equation_func_simple(gint i, gint j,
				gdouble *G, gdouble *dGdn, 
				gint n, gpointer *e)
{
  gint np = *(gint *)e[0] ;
  gdouble *A = (gdouble *)e[1] ;
  gdouble *B = (gdouble *)e[2] ;

  g_assert(n == 1) ;

  B[i*np+j] += G[0] ; A[i*np+j] += dGdn[0] ;

  return BEM3D_SUCCESS ;
}

static void element_set_bc(BEM3DElement *e, gpointer *data)

{
  GHashTable *h = (GHashTable *)data[BEM3D_BMESH_DATA_HASH] ;
  BEM3DBCFunc bcf = (BEM3DBCFunc)data[BEM3D_BMESH_DATA_BFUNC] ;
  gpointer bdata = data[BEM3D_BMESH_DATA_BDATA] ;
  static gdouble *dLds, *dLdt ;
  gdouble s, t, J ;
  GtsVector n ;
  gint i ;

  if ( dLds == NULL ) {
    dLds = (gdouble *)g_malloc(4*bem3d_element_node_number(e)*sizeof(gdouble)) ;
    dLdt = (gdouble *)g_malloc(4*bem3d_element_node_number(e)*sizeof(gdouble)) ;
  }

  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    if ( !g_hash_table_lookup(h, GINT_TO_POINTER((e->i[i]+1))) ) {
      g_hash_table_insert(h, GINT_TO_POINTER(e->i[i]+1), 
			  GINT_TO_POINTER(e->i[i]+1)) ;
      s = bem3d_element_vertex_xi(e,i) ;
      t = bem3d_element_vertex_eta(e,i) ;
      e->shf(s, t, NULL, dLds, dLdt, NULL) ;
      bem3d_element_normal(e, dLds, dLdt, n, &J) ;
      bcf(GTS_VERTEX(e->c[i]), n, e->i[i], bdata) ;
    }
  }

  return ;
}

/** 
 * Set the boundary conditions on a mesh
 * 
 * @param m mesh for boundary conditions
 * @param bcf BEM3DBCFunc
 * @param bdata user data to pass to bcf 
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_mesh_set_bc(BEM3DMesh *m, BEM3DBCFunc bcf, gpointer bdata)

{
  gpointer data[BEM3D_BMESH_DATA_WIDTH] ;
  GHashTable *h ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(bcf != NULL, BEM3D_NULL_ARGUMENT) ;

  h = g_hash_table_new(NULL, NULL) ;
  data[BEM3D_BMESH_DATA_HASH] = h ;
  data[BEM3D_BMESH_DATA_BFUNC] = (gpointer)bcf ;
  data[BEM3D_BMESH_DATA_BDATA] = bdata ;

  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)element_set_bc, data) ;
  
  g_hash_table_destroy(h) ;

  return BEM3D_SUCCESS ;
}

static void mesh_quad_gfunc(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMesh *m = (BEM3DMesh *)data[BEM3D_BMESH_DATA_MESH] ;
  BEM3DConfiguration *config = data[BEM3D_BMESH_DATA_CONFIG] ; 
  BEM3DParameters *gdata = data[BEM3D_BMESH_DATA_GDATA] ;
  BEM3DLookupFunc lf = (BEM3DLookupFunc)data[BEM3D_BMESH_DATA_LFUNC] ;
  gpointer ldata = data[BEM3D_BMESH_DATA_LDATA] ;
  BEM3DEquationFunc efunc = (BEM3DEquationFunc)data[BEM3D_BMESH_DATA_EFUNC] ;
  gpointer edata = data[BEM3D_BMESH_DATA_EDATA] ;
  gint j ;
  static GArray *g = NULL, *f = NULL, *zero = NULL ;

  if ( g == NULL ) {
    f = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
    zero = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
    g = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
  }

  for ( j = 0 ; j < f->len ; j ++ ) g_array_index(f,gdouble,j)=0.0 ;
  bem3d_mesh_radiation_point(m, config, gdata,
			     lf, ldata, GTS_POINT(v), f) ;
  g_array_set_size(zero, f->len) ;
  efunc(i, i, 
	&(g_array_index(zero,gdouble,0)),
	&(g_array_index(f,gdouble,0)),
	f->len, edata) ;

  return ;
}

/** 
 * Integrate the normal derivative of the Green's function over a
 * mesh.
 * 
 * @param m BEM3DMEsh 
 * @param config ::BEM3DConfiguration for the problem;
 * @param gdata ::BEM3DParameters to pass to Green's function;
 * @param lfunc BEM3DLookupFunc; use ::bem3d_lookup_func_unit 
 * (for real problems) 
 * or bem3d_lookup_func_unit_c (for complex) 
 * @param ldata data to pass to the lookup function 
 * @param efunc equation function for problem. This will be called with 
 * the row and column numbers equal.
 * @param edata user data to pass to efunc
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_mesh_quad_dgdn(BEM3DMesh *m,
			  BEM3DConfiguration *config,
			  BEM3DParameters *gdata,
			  BEM3DLookupFunc lfunc, gpointer ldata,
			  BEM3DEquationFunc efunc, gpointer edata)

{
  gpointer data[BEM3D_BMESH_DATA_WIDTH] ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(gdata != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(efunc != NULL, BEM3D_NULL_ARGUMENT) ;

  data[BEM3D_BMESH_DATA_MESH] = m ; 
  data[BEM3D_BMESH_DATA_CONFIG] = config ; 
  data[BEM3D_BMESH_DATA_GDATA] = gdata ;
  data[BEM3D_BMESH_DATA_LFUNC] = (gpointer)lfunc ;
  data[BEM3D_BMESH_DATA_LDATA] = ldata ;
  data[BEM3D_BMESH_DATA_IMIN] = &(bem3d_mesh_node_index_min(m)) ;
  data[BEM3D_BMESH_DATA_IMAX] = &(bem3d_mesh_node_index_max(m)) ;
  data[BEM3D_BMESH_DATA_EFUNC] = efunc ;
  data[BEM3D_BMESH_DATA_EDATA] = edata ;

  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)mesh_quad_gfunc,
			  (gpointer )(&data[0])) ;

  return BEM3D_SUCCESS ;
}

static void add_mesh(BEM3DElement *e, BEM3DMesh *m)

{
  bem3d_mesh_add_element(m, e, FALSE) ;

  return ;
}

gint bem3d_mesh_merge(BEM3DMesh *m, BEM3DMesh *n)

{  
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(n != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(n), BEM3D_ARGUMENT_WRONG_TYPE) ;
  
  bem3d_mesh_foreach_element(n, (BEM3DElementFunc)add_mesh, m) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Find the ::BEM3DElement of a ::BEM3DMesh containing a given GtsFace. 
 * 
 * @param m BEM3DMesh to search
 * @param f GtsFace to locate
 * 
 * @return ::BEM3DElement of \a m containing \a f or NULL if \a f is not
 * in \a m. 
 */

BEM3DElement *bem3d_mesh_face_element(BEM3DMesh *m, GtsFace *f)

{
  BEM3DElement *e ;

  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;
  g_return_val_if_fail(f != NULL, NULL) ;
  g_return_val_if_fail(GTS_IS_FACE(f), NULL) ;

  e = g_hash_table_lookup(m->f, f) ;

  return e ;
}

/** 
 * Find the GtsVertex on a ::BEM3DMesh containing a given collocation
 * point.
 * 
 * @param m BEM3DMesh to search;
 * @param i index of collocation point.
 * 
 * @return GtsVertex of \a m with collocation point index \a i or NULL
 * if point \a i is not on \a m.
 */

GtsVertex *bem3d_mesh_node_from_index(BEM3DMesh *m, gint i)

{
  GtsVertex *v = NULL ;

  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;

  v = g_hash_table_lookup(m->c, GINT_TO_POINTER(i+1)) ;

  return v ;
}

static void bem3d_element_clear_reserved(BEM3DElement *e, gpointer data)

{
  e->reserved = NULL ;

  return ;
}

gint bem3d_mesh_element_clear_reserved(BEM3DMesh *m)

{
  bem3d_mesh_foreach_element(m, 
			     (BEM3DElementFunc)bem3d_element_clear_reserved, 
			     NULL) ;
  return BEM3D_SUCCESS ;
}

/** 
 * Find an index for a collocation point in a ::BEM3DMesh. Note that
 * this index may not be unique, for example if the mesh has sharp edges.
 * 
 * @param m ::BEM3DMesh
 * @param v GtsVertex for collocation point to find
 * 
 * @return index of \a v in \a m, or -1 if not found.
 */

gint bem3d_mesh_index_from_node(BEM3DMesh *m, GtsVertex *v)

{
  gint i = -1, j ;
  GSList *e ;
  
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(v != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_VERTEX(v), BEM3D_ARGUMENT_WRONG_TYPE) ;

  e = bem3d_mesh_vertex_elements(m, v) ;

  if ( e == NULL ) return i ;

  j = bem3d_element_find_node(BEM3D_ELEMENT(e->data), v) ;
  i = bem3d_element_global_index(BEM3D_ELEMENT(e->data), j) ;

  return i ;
}

gint bem3d_mesh_remove_element(BEM3DMesh *m, BEM3DElement *e)

{
  gint i ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  
  if ( !bem3d_element_has_parent_mesh(m, e) ) return 0 ;

  g_hash_table_remove(m->e, e) ;
  for ( i = 0 ; i < bem3d_element_face_number(e) ; i ++ ) {
    g_hash_table_remove(m->f, e->f[i]) ;
    gts_surface_remove_face(GTS_SURFACE(m), e->f[i]) ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Extract a linked list of the elements of a mesh containing a given
 * node (collocation point).
 * 
 * @param m BEM3DMesh;
 * @param i index of node on \a m.
 * 
 * @return linked list of elements containing node \a i, NULL otherwise.
 */

GSList *bem3d_mesh_node_elements(BEM3DMesh *m, gint i)

{
  GSList *f, *e, *j ;
  GtsVertex *v ;
  BEM3DElement *el ;

  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;

  if ( (v = bem3d_mesh_node_from_index(m, i)) == NULL ) return NULL ;
  
  f = gts_vertex_faces(v, GTS_SURFACE(m), NULL) ;
    
  for ( (e = NULL), (j = f) ; j != NULL ; j = j->next ) {
    el = g_hash_table_lookup(m->f, j->data) ;
    if ( !g_slist_find(e, el) && 
	 (bem3d_element_find_index(el, i) != -1 ) ) e = g_slist_prepend(e, el) ;
  }

  g_slist_free(f) ;

  return e ;
}

static gint element_moments(BEM3DElement *e, gint *H)

{  
  bem3d_element_moments_make(e, (*H)) ;

  return 0 ;
}

gint bem3d_mesh_element_moments(BEM3DMesh *m, gint H)

{
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)element_moments, &H) ;

  return BEM3D_SUCCESS ;
}

static gint biggest_element(BEM3DElement *e, gint *nv)

{
  *nv = MAX(*nv, bem3d_element_node_number(e)) ;

  return 0 ;
}

/** 
 * Find the largest number of nodes on an element of a ::BEM3DMesh
 * 
 * @param m a ::BEM3DMesh
 * 
 * @return the maximum number of nodes on an element of \a m 
 */

gint bem3d_mesh_element_node_number_max(BEM3DMesh *m)

{
  gint nv = 0 ;

  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)biggest_element, &nv) ;

  return nv ;
}

static gint index_range(gint i, GtsVertex *v, gpointer *data)

{
  gint *imin = data[0] ;
  gint *imax = data[1] ;

  *imin = MIN(i, *imin) ;
  *imax = MAX(i, *imax) ;

  return 0 ;
}

/** 
 * Find the maximum and minimum indices of a ::BEM3DMesh
 * 
 * @param m a ::BEM3DMesh
 * @param imin minimum index
 * @param imax maximum index
 * 
 * @return BEM3D_SUCCESS on success
 */

gint bem3d_mesh_index_range(BEM3DMesh *m, gint *imin, gint *imax)

{
  gpointer data[2] = {imin, imax} ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(imin, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(imax, BEM3D_NULL_ARGUMENT) ;

  *imin = G_MAXINT ; *imax = -1 ;
  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)index_range, data) ;

  return 0 ;
}

static gint element_sample(BEM3DElement *e, BEM3DElement **s)

{
  *s = e ;

  return 1 ;
}

/** 
 * Select a sample element from a ::BEM3DMesh
 * 
 * @param m a ::BEM3DMesh
 * 
 * @return a randomly chosen element from \a m
 */

BEM3DElement *bem3d_mesh_element_sample(BEM3DMesh *m)

{
  BEM3DElement *e ;

  g_return_val_if_fail(m != NULL, NULL) ;

  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)element_sample, &e) ;
  
  return e ;
}

static gint element_area(BEM3DElement *e, gpointer data[])

{
  gint *ngp = data[0] ;
  gdouble *S = data[1] ;

  *S += bem3d_element_area(e, *ngp) ;

  return 0 ;
}

/** 
 * Surface area of a ::BEM3DMesh by integrating over elements. This
 * gives an accurate estimate of surface area for meshes with
 * curvilinear elements.
 * 
 * @param m a ::BEM3DMesh
 * @param ngp number of points in quadrature rule, for call to 
 * ::bem3d_element_area 
 * 
 * @return surface area of \a m
 */

gdouble bem3d_mesh_surface_area(BEM3DMesh *m, gint ngp)

{
  gdouble S ;
  gpointer data[] = {&ngp, &S} ;

  S = 0 ;

  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)element_area, data) ;

  return S ;
}

/**
 * @}
 * 
 */
