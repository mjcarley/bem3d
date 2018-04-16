/* edge.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

#include <bem3d.h>

/* BEM3DEdge: Object */

static void bem3d_edge_class_init (BEM3DEdgeClass * klass)
{
  /* define new methods and overload inherited methods here */

}

static void bem3d_edge_init (BEM3DEdge * object)
{
  /* initialize object here */
  object->i = g_array_new(FALSE, FALSE, sizeof(gint)) ;
  object->e = g_ptr_array_new() ;
  object->v = g_ptr_array_new() ;
}

/**
 * @defgroup edge BEM3D edges
 *
 * A data type and set of functions which allow the handling of edges,
 * oriented lists of sharp nodes on a mesh.
 *
 * @{
 * 
 */


BEM3DEdgeClass * bem3d_edge_class (void)
{
  static BEM3DEdgeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo bem3d_edge_info = {
      "BEM3DEdge",
      sizeof (BEM3DEdge),
      sizeof (BEM3DEdgeClass),
      (GtsObjectClassInitFunc) bem3d_edge_class_init,
      (GtsObjectInitFunc) bem3d_edge_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &bem3d_edge_info);
  }

  return klass;
}

/** 
 * Allocate a new ::BEM3DEdge.
 * 
 * @param klass a ::BEM3DEdgeClass.
 * 
 * @return the newly allocated object.
 */

BEM3DEdge * bem3d_edge_new (BEM3DEdgeClass * klass)
{
  BEM3DEdge * object;

  object = BEM3D_EDGE (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

static gboolean int_in_list(gint id[], gint n, gint i)

{
  gint j ;
  for ( j = 0 ; j < n ; j ++ ) if ( id[j] == i ) return TRUE ;
  return FALSE ;
}

/** 
 * Find the number of indices assigned to a mesh vertex.
 * 
 * @param m a BEM3DMesh;
 * @param v a GtsVertex forming part of \a m.
 * 
 * @return the number of indices used by \a v on \a m.
 */

gint bem3d_mesh_vertex_index_number(BEM3DMesh *m, GtsVertex *v)

{
  GSList *i, *el = NULL ;
  gint idx[32], j, k ;
  
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

  el = bem3d_mesh_vertex_elements(m, v) ;

  for ( (j = 0), (i = el) ; i != NULL ; i = i->next ) {
    k = bem3d_element_find_vertex(BEM3D_ELEMENT(i->data), v) ;
    k = bem3d_element_global_index(BEM3D_ELEMENT(i->data), k) ;
    if ( !int_in_list(idx, j, k) ) {
      idx[j] = k ; j ++ ;
    }
  }
  return j ;
}

static void sharp_nodes(GtsVertex *v, gpointer data[])

{
  BEM3DMesh *m = data[0] ;
  GSList **list = data[1] ;
  
  if ( bem3d_mesh_vertex_index_number(m, v) <= 1 ) return ;

  *list = g_slist_prepend(*list, v) ;

  return ;
}

/** 
 * Find the sharp nodes of a ::BEM3DMesh, i.e. those with multiple
 * indices, as found by ::bem3d_mesh_vertex_index_number.
 * 
 * @param m a ::BEM3DMesh.
 * 
 * @return a GSList of the GtsVertex's of \a m which have multiple indices. 
 */

GSList *bem3d_mesh_sharp_vertices(BEM3DMesh *m)

{
  GSList *v  = NULL ;
  gpointer data[2] ;

  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;

  data[0] = m ; data[1] = &v ;

  gts_surface_foreach_vertex(GTS_SURFACE(m), (GtsFunc)sharp_nodes, data) ;

  return v ;
}

gint bem3d_edge_add_node(BEM3DEdge *e, GtsVertex *v, gint i, gint j)

{  
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e), BEM3D_ARGUMENT_WRONG_TYPE) ;

  g_array_append_val(e->i,i) ; g_array_append_val(e->i,j) ;
  g_ptr_array_add(e->v, v) ;

  return BEM3D_SUCCESS ;
}

static gboolean vertices_on_sharp_edge(GtsVertex *v1, GtsVertex *v2,
				       BEM3DMesh *m)

{
  BEM3DElement *e1, *e2 ;
  GSList *el ;
  gint i, i1, i2 ;

  if ( gts_vertices_are_connected(v1, v2) == NULL ||
       v1 == v2 ) return FALSE ;

  el = bem3d_elements_from_vertices(m, v1, v2) ;
  if ( g_slist_length(el) != 2 ) return FALSE ;

  if ( g_slist_length(el) != 2 ) {
    fprintf(stderr, "Non-conforming sharp edge at: ") ;
    for ( ; el != NULL ; el = el->next ) {
      fprintf(stderr, "(%lg,%lg,%lg) ", 
	      GTS_POINT(el->data)->x,
	      GTS_POINT(el->data)->y,
	      GTS_POINT(el->data)->z) ;
    }
    g_error("%s: ", __FUNCTION__) ;
  }

  e1 = el->data ; e2 = el->next->data ;

  i = bem3d_element_find_vertex(e1, v1) ;
  i1 = bem3d_element_global_index(e1, i) ;
  i = bem3d_element_find_vertex(e2, v1) ;
  i2 = bem3d_element_global_index(e2, i) ;

  if ( i1 == i2 ) return FALSE ;
  return TRUE ;
}

static GSList *connect_edges(GSList **v, BEM3DMesh *m)

{
  GSList *i, *ee = NULL, *en ;
  gboolean add_start, add_end ;
  gint added ;

  ee = *v ;
  while ( (ee != NULL) && 
	  (bem3d_mesh_vertex_index_number(m, GTS_VERTEX(ee->data)) != 2) ) 
    ee = ee->next ;
  if ( ee == NULL ) return NULL ;

  en = g_slist_prepend(NULL, ee->data) ;
  ee = en ;
  add_start = TRUE ; add_end = TRUE ;

  added = 1 ;

  *v = g_slist_remove(*v, ee->data) ;

  while ( added != 0 ) {
    added = 0 ;
    for ( i = *v ; i != NULL ; i = i->next ) {
      if ( bem3d_mesh_vertex_index_number(m, GTS_VERTEX(ee->data)) != 2 )
	add_start = FALSE ;
      if ( bem3d_mesh_vertex_index_number(m, GTS_VERTEX(en->data)) != 2 )
	add_end = FALSE ;
      if ( vertices_on_sharp_edge(GTS_VERTEX(i->data),
				  GTS_VERTEX(ee->data), m)
	   && add_start ) {
	ee = g_slist_prepend(ee, i->data) ;
	added ++ ;
      } else
	if ( vertices_on_sharp_edge(GTS_VERTEX(i->data),
				    GTS_VERTEX(en->data), m) 
	     && add_end ) {
	  ee = g_slist_append(ee, i->data) ;
	  en = g_slist_last(ee) ;
	  added ++ ;
	}

/*       if ( gts_vertices_are_connected(GTS_VERTEX(i->data), */
/* 				      GTS_VERTEX(ee->data))  */
/* 	   && add_start && */
/* 	   g_slist_find(ee, i->data) == NULL ) { */
/* 	ee = g_slist_prepend(ee, i->data) ; */
/* 	added ++ ; */
/*       } else */
/* 	if ( gts_vertices_are_connected(GTS_VERTEX(i->data), */
/* 				      GTS_VERTEX(en->data))  */
/* 	     && add_end && */
/* 	   g_slist_find(ee, i->data) == NULL ) { */
/* 	  ee = g_slist_append(ee, i->data) ; */
/* 	  en = g_slist_last(ee) ; */
/* 	  added ++ ; */
/* 	} */
    }

    for ( i = ee ; i != NULL ; i = i->next ) 
      if ( bem3d_mesh_vertex_index_number(m, GTS_VERTEX(i->data)) == 2 ) 
	*v = g_slist_remove(*v, i->data) ;
  }

  return ee ;
}

static void orient_edges(BEM3DMesh *m, GSList *e, BEM3DEdge *edge)

{
  GSList *el, *i ;
  GtsVector n1, n2 ;
  gint i1, i2 ;
  BEM3DElement *el1, *el2, *eu ;
  GtsPoint p1, p2 ;
  gdouble ee = 1e-3 ;

  el1 = el2 = eu = NULL ;

  for ( i = e ; i != NULL && i->next != NULL ; i = i->next ) {
    el = bem3d_elements_from_vertices(m, i->data, i->next->data) ;
    g_assert(g_slist_length(el) == 2) ;
    el1 = el->data ; el2 = el->next->data ;

    i1 = bem3d_element_find_vertex(el1, i->data) ;
    i1 = bem3d_element_global_index(el1, i1) ;

/*     i2 = i1 ; */
/*     while ( i2 == i1 ) { */
      el = el->next ;
/*       g_assert(el != NULL) ; */
      el2 = el->data ;
      i2 = bem3d_element_find_vertex(el2, i->data) ;
      i2 = bem3d_element_global_index(el2, i2) ;
/*     } */

    bem3d_node_normal(m, i1, n1, BEM3D_AVERAGE_MWE) ;
    bem3d_node_normal(m, i2, n2, BEM3D_AVERAGE_MWE) ;
    
    p1.x = GTS_POINT(i->data)->x + ee*n1[0] ;
    p1.y = GTS_POINT(i->data)->y + ee*n1[1] ;
    p1.z = GTS_POINT(i->data)->z + ee*n1[2] ;
    p2.x = GTS_POINT(i->data)->x + ee*n2[0] ;
    p2.y = GTS_POINT(i->data)->y + ee*n2[1] ;
    p2.z = GTS_POINT(i->data)->z + ee*n2[2] ;

    if ( gts_point_orientation_3d(GTS_POINT(i->data),
				  GTS_POINT(i->next->data),
				  &p1, &p2) > 0.0 ) {
      bem3d_edge_add_node(edge, i->data, i1, i2) ;
      eu = el1 ;
    } else {
      bem3d_edge_add_node(edge, i->data, i2, i1) ;
      eu = el2 ;
    }
  }

    if ( el1 == NULL || el2 == NULL ) 
      g_error("%s: vertices %p and %p "
	      "not on common element", __FUNCTION__, 
	      i->data, i->next->data) ;

  i1 = bem3d_element_find_vertex(el1, i->data) ;
  i2 = bem3d_element_find_vertex(el2, i->data) ;
  i1 = bem3d_element_global_index(el1, i1) ;
  i2 = bem3d_element_global_index(el2, i2) ;
  if ( eu == el1 ) bem3d_edge_add_node(edge, i->data, i1, i2) ;
  else bem3d_edge_add_node(edge, i->data, i2, i1) ;

  return ;
}

/** 
 * Extract the edge connectivity data from a list of sharp nodes.
 * 
 * @param m a ::BEM3DMesh; 
 * @param e GSList of vertices of \a m which are sharp (from 
 * ::bem3d_mesh_sharp_vertices).
 * 
 * @return a GSList of ::BEM3DEdges, each of which is a connected,
 * oriented edge on \a m.
 */

GSList *bem3d_mesh_extract_edges(BEM3DMesh *m, GSList **e)

{
  GSList *ee, *edges = NULL ;
  gboolean stop ;

  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;

  stop = FALSE ;
  do {
    if ( (ee = connect_edges(e, m)) != NULL && g_slist_length(ee) > 1 ) {
      edges = g_slist_prepend(edges, bem3d_edge_new(bem3d_edge_class())) ;
      orient_edges(m, ee, edges->data) ;
    }
    if ( ee == NULL ) stop = TRUE ;
    if ( g_slist_length(ee) == 1 ) 
      *e = g_slist_remove(*e, ee->data) ;
  } while ( !stop ) ;
/* while ( ee != NULL && g_slist_length(ee) > 1) ; */

  return edges ;
}

/** 
 * Clear a ::BEM3DEdge, i.e. set its number of nodes to zero.
 * 
 * @param e the ::BEM3DEdge to clear.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_edge_clear(BEM3DEdge *e)

{
  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e), BEM3D_ARGUMENT_WRONG_TYPE) ;

  g_array_set_size(e->i, 0) ;
  g_ptr_array_set_size(e->e, 0) ;
  g_ptr_array_set_size(e->v, 0) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Copy a ::BEM3DEdge. This is a `deep' copy, i.e. the GtsVertex at
 * each node is copied into a new instance, rather than being linked
 * to the same address, so that the copy can be swept without
 * affecting the original points.
 * 
 * @param e the destination ::BEM3DEdge;
 * @param f the source ::BEM3DEdge.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_edge_copy(BEM3DEdge *e, BEM3DEdge *f)

{
  gint i ;
  GtsVertex *v ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(f), BEM3D_ARGUMENT_WRONG_TYPE) ;

  for ( i = 0 ; i < bem3d_edge_node_number(f) ; i ++ ) {
    v = gts_vertex_new((GtsVertexClass *)
		       (GTS_OBJECT(bem3d_edge_vertex(f,i))->klass),
		       GTS_POINT(bem3d_edge_vertex(f,i))->x,
		       GTS_POINT(bem3d_edge_vertex(f,i))->y,
		       GTS_POINT(bem3d_edge_vertex(f,i))->z) ;
    bem3d_edge_add_node(e, v,
			bem3d_edge_node_index_upper(f,i),
			bem3d_edge_node_index_lower(f,i)) ;
  }

  for ( i = 0 ; i < bem3d_edge_node_number(e)-1 ; i ++ )
    gts_edge_new(gts_edge_class(), 
		 bem3d_edge_vertex(e,i), bem3d_edge_vertex(e,i+1)) ;

  return BEM3D_SUCCESS ;
}

static void index_element(BEM3DMesh *m, BEM3DElement *e, gint *i)

{
  gint j, k ;

  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
    if ( (k = bem3d_mesh_index_from_node(m, bem3d_element_node(e,j))) == -1 ) {
      bem3d_element_global_index(e,j) = (*i) ;
      (*i) ++ ;
    } else bem3d_element_global_index(e,j) = k ;
      
  }

  return ;
}

gint bem3d_link_edges(BEM3DMesh *m, BEM3DEdge *e1, BEM3DEdge *e2,
		      BEM3DElementBuildFunc build, gint nv, 
		      gint nc, gint *idx)
{
  gint i ;
  GtsEdge *e[32] ;
  GtsVertex *v[16] ;
  BEM3DElement *el ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e1 != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e1), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(e2 != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e2), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(build != NULL, BEM3D_NULL_ARGUMENT) ;

  g_assert(bem3d_edge_node_number(e1) == bem3d_edge_node_number(e2)) ;

  switch ( nc ) {
  default: 
    g_error("%s: nc=%d-sided element undefined", __FUNCTION__, nc) ;
    break ;
  case 3:
    for ( i = 0 ; i <  bem3d_edge_node_number(e1)-nv+1 ; i += nv-1 ) {
      v[0] = bem3d_edge_vertex(e1,i) ;
      v[1] = bem3d_edge_vertex(e1,i+nv-1) ;
      v[2] = bem3d_edge_vertex(e2,i) ;

      memset(e, 0, 32*sizeof(gpointer)) ;
      el = build(e, v) ;
      index_element(m, el, idx) ;
      bem3d_mesh_add_element(m, el, FALSE) ;

      v[0] = bem3d_edge_vertex(e2,i+nv-1) ;
      v[1] = bem3d_edge_vertex(e2,i) ;
      v[2] = bem3d_edge_vertex(e1,i+nv-1) ;

      memset(e, 0, 32*sizeof(gpointer)) ;
      el = build(e, v) ;
      index_element(m, el, idx) ;
      bem3d_mesh_add_element(m, el, FALSE) ;
    }
    break ;
  case 4:
    for ( i = 0 ; i <  bem3d_edge_node_number(e1)-nv+1 ; i += nv-1 ) {
      v[0] = bem3d_edge_vertex(e1,i) ;
      v[1] = bem3d_edge_vertex(e1,i+nv-1) ;
      v[2] = bem3d_edge_vertex(e2,i+nv-1) ;
      v[3] = bem3d_edge_vertex(e2,i) ;

      memset(e, 0, 32*sizeof(gpointer)) ;
      el = build(e, v) ;
      index_element(m, el, idx) ;
      bem3d_mesh_add_element(m, el, FALSE) ;
    }
    break ;
  }

  return 0 ;
}

/** 
 * Write a ::BEM3DEdge to a file. The output format is one line:
 *
 * @e N BEM3DEdge
 *
 * where @e N is the number of edge nodes followed by @e N pairs of
 * indices, one pair per line. 
 * 
 * @param e the ::BEM3DEdge to write;
 * @param f the file stream to which @a e should be written.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_edge_write(BEM3DEdge *e, FILE *f)

{
  gint i ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  fprintf(f, "%d BEM3DEdge\n", bem3d_edge_node_number(e)) ;
  for ( i = 0 ; i < bem3d_edge_node_number(e) ; i ++ ) {
    fprintf(f, "%d\t %d\n", 
	    bem3d_edge_node_index_upper(e,i),
	    bem3d_edge_node_index_lower(e,i)) ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Read a ::BEM3DEdge from a file. The format is that used by
 * ::bem3d_edge_write.
 * 
 * @param e a ::BEM3DEdge to which the edge in \a f will be added;
 * @param f a ::GtsFile.
 * 
 * @return BEM3D_SUCCESS on success.
 */

guint bem3d_edge_read(BEM3DEdge *e, GtsFile *f)

{
  gint n, nn, upper, lower ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( f->type != GTS_INT ) {
    gts_file_error (f, "expecting an integer (number of edge nodes)");
    return f->line;
  }
  nn = atoi (f->token->str) ;
  
  gts_file_next_token (f);
  if (f->type != GTS_STRING) {
    gts_file_error (f, "expecting a string (BEM3DEdge)");
    return f->line;
  }  
  gts_file_first_token_after(f, '\n') ;

  if ( nn <= 0 ) return 0 ;

  n = 0;
  while (n < nn && f->type != GTS_ERROR) {
    if ( f->type != GTS_INT ) {
      gts_file_error (f, "expecting an integer (upper node index)");
      return f->line;
    }
    upper = atoi(f->token->str) ;
    gts_file_next_token(f) ;
    if ( f->type != GTS_INT ) {
      gts_file_error (f, "expecting an integer (lower node index)");
      return f->line;
    }
    lower = atoi(f->token->str) ;

    bem3d_edge_add_node(e, NULL, upper, lower) ;

    gts_file_first_token_after(f, '\n') ;
  }

  return BEM3D_SUCCESS ;
}

gint bem3d_edge_link_to_mesh(BEM3DEdge *e, BEM3DMesh *m)

{
  gint i ;
  GtsVertex *v ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

  for ( i = 0 ; i < bem3d_edge_node_number(e) ; i ++ ) {
    v = bem3d_mesh_node_from_index(m, bem3d_edge_node_index_upper(e,i)) ;
    if ( v == NULL ) 
      g_error("%s: node %d not found on mesh m", __FUNCTION__,
	      bem3d_edge_node_index_upper(e,i)) ;
    if ( v != bem3d_mesh_node_from_index(m, bem3d_edge_node_index_lower(e,i)) )
      g_error("%s: nodes %d and %d do not refer to the same vertex",
	      __FUNCTION__, bem3d_edge_node_index_upper(e,i),
	      bem3d_edge_node_index_lower(e,i)) ;
      
    bem3d_edge_vertex(e,i) = v ;
  }

  bem3d_edge_mesh(e) = m ;

  return BEM3D_SUCCESS ;
}

static gboolean triangle_has_vertices(GtsTriangle *t,
				      GtsVertex *v,
				      GtsVertex *w)

{
  if ( GTS_SEGMENT(t->e1)->v1 == v && GTS_SEGMENT(t->e1)->v2 == w )
    return TRUE ;
  if ( GTS_SEGMENT(t->e1)->v1 == w && GTS_SEGMENT(t->e1)->v2 == v )
    return TRUE ;
  if ( GTS_SEGMENT(t->e2)->v1 == v && GTS_SEGMENT(t->e2)->v2 == w )
    return TRUE ;
  if ( GTS_SEGMENT(t->e2)->v1 == w && GTS_SEGMENT(t->e2)->v2 == v )
    return TRUE ;
  if ( GTS_SEGMENT(t->e3)->v1 == v && GTS_SEGMENT(t->e3)->v2 == w )
    return TRUE ;
  if ( GTS_SEGMENT(t->e3)->v1 == w && GTS_SEGMENT(t->e3)->v2 == v )
    return TRUE ;

  return FALSE ;
}

static GtsFace *face_from_element_and_vertices(BEM3DElement *e,
					       GtsVertex *v,
					       GtsVertex *w)

{
  gint i ;
  GtsFace *f ;

  for ( i = 0 ; i < bem3d_element_face_number(e) ; i ++ ) {
    f = bem3d_element_face(e,i) ;
    if ( triangle_has_vertices(GTS_TRIANGLE(f), v, w) ) return f ;
  }

  return NULL ;
}

/** 
 * Check if a ::BEM3DEdge is aligned with its parent ::BEM3DMesh, by
 * comparing the direction of the edge with the corresponding edges on
 * the mesh.
 * 
 * @param e a ::BEM3DEdge;
 * @param m the parent ::BEM3DMesh of \a e.
 * 
 * @return TRUE if \a e is compatibly aligned with the sharp edges of \a m.
 */

gboolean bem3d_edge_is_oriented(BEM3DEdge *e, BEM3DMesh *m)

{
  gint i ;
  GtsVertex *v, *w, *v1, *v2, *v3 ;
  GSList *el ;
  BEM3DElement *eu ;
  GtsFace *f ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

  v = bem3d_mesh_node_from_index(m, bem3d_edge_node_index_upper(e,0)) ;
  g_assert(v != NULL) ;
  for ( i = 1 ; i < bem3d_edge_node_number(e) ; i ++ ) {
    w = bem3d_mesh_node_from_index(m, bem3d_edge_node_index_upper(e,i)) ;    
    g_assert(w != NULL) ;
    el = bem3d_elements_from_vertices(m, v, w) ;
    g_assert(g_slist_length(el) == 2) ;
    
    if ( bem3d_element_find_index(BEM3D_ELEMENT(el->data), 
				  bem3d_edge_node_index_upper(e,i)) != -1 )
      eu = el->data ;
    else
      eu = el->next->data ;

    f = face_from_element_and_vertices(eu, v, w) ;

    gts_triangle_vertices(GTS_TRIANGLE(f), &v1, &v2, &v3) ;

    if ( w == v1 && v == v2 ) return FALSE ;
    if ( w == v2 && v == v3 ) return FALSE ;
    if ( w == v3 && v == v1 ) return FALSE ;

    v = w ;
  }
  
  return TRUE ;
}

/** 
 * Invert a ::BEM3DEdge by swapping the upper and lower indices to
 * properly align it with its parent ::BEM3DMesh.
 * 
 * @param e a ::BEM3DEdge.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_invert_edge(BEM3DEdge *e)

{
  gint i, t ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(e), BEM3D_ARGUMENT_WRONG_TYPE) ;

  for ( i = 0 ; i < bem3d_edge_node_number(e) ; i ++ ) {
    t = bem3d_edge_node_index_upper(e,i) ;
    bem3d_edge_node_index_upper(e,i)  = bem3d_edge_node_index_lower(e,i) ;
    bem3d_edge_node_index_lower(e,i) = t ;
  }

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */

