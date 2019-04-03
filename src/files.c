/* files.c
 * 
 * Copyright (C) 2006, 2009, 2018 Michael Carley
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
 * @addtogroup mesh
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

#define GMSH_ELEMENT_LINE          1
#define GMSH_ELEMENT_TRIANGLE      2
#define GMSH_ELEMENT_QUADRANGLE    3
#define GMSH_ELEMENT_TETRAHEDRON   4
#define GMSH_ELEMENT_HEXAHEDRON    5
#define GMSH_ELEMENT_PRISM         6
#define GMSH_ELEMENT_PYRAMID       7
#define GMSH_ELEMENT_LINE_2        8
#define GMSH_ELEMENT_TRIANGLE_2    9
#define GMSH_ELEMENT_QUADRANGLE_2  10
#define GMSH_ELEMENT_TETRAHEDRON_2 11
#define GMSH_ELEMENT_HEXAHEDRON_2  12
#define GMSH_ELEMENT_PRISM_2       13
#define GMSH_ELEMENT_PYRAMID_2     14
#define GMSH_ELEMENT_POINT         15
#define GMSH_ELEMENT_TRIANGLE_3    21
#define GMSH_ELEMENT_LINE_3        26

static void write_vertex (GtsPoint * p, gpointer * data)
{
  (*GTS_OBJECT (p)->klass->write) (GTS_OBJECT (p), (FILE *) data[0]);
  if (!GTS_POINT_CLASS (GTS_OBJECT (p)->klass)->binary)
    fputc ('\n', (FILE *) data[0]);
  g_hash_table_insert (data[2], p, 
		       GUINT_TO_POINTER (++(*((guint *) data[1]))));
}

static void write_edge (GtsSegment * s, gpointer * data) 
{
  fprintf ((FILE *) data[0], "%u %u",
	   GPOINTER_TO_UINT (g_hash_table_lookup (data[2], s->v1)),
	   GPOINTER_TO_UINT (g_hash_table_lookup (data[2], s->v2)));
  if (GTS_OBJECT (s)->klass->write)
    (*GTS_OBJECT (s)->klass->write) (GTS_OBJECT (s), (FILE *) data[0]);
  fputc ('\n', (FILE *) data[0]);
  g_hash_table_insert (data[3], s, 
		       GUINT_TO_POINTER (++(*((guint *) data[1]))));
}

static void write_face (GtsTriangle * t, gpointer * data)

{
  fprintf (data[0], "%u %u %u",
	   GPOINTER_TO_UINT (g_hash_table_lookup (data[3], t->e1)),
	   GPOINTER_TO_UINT (g_hash_table_lookup (data[3], t->e2)),
	   GPOINTER_TO_UINT (g_hash_table_lookup (data[3], t->e3)));
  if (GTS_OBJECT (t)->klass->write)
    (*GTS_OBJECT (t)->klass->write) (GTS_OBJECT (t), data[0]);
  fputc ('\n', data[0]);
  g_hash_table_insert (data[4], t, 
		       GUINT_TO_POINTER (++(*((guint *) data[1]))));
}

static void write_element(gpointer key, gpointer value, 
			  gpointer *data)

{
  bem3d_element_write(key, data[2], data[4], data[0]) ;
  
  return ;
}

/** 
 * Write a BEM3DMesh to a file using part of Stephane Popinet's GTS
 * function for writing a GTS surface.
 * 
 * @param m mesh to write;
 * @param f file pointer to write to.
 * 
 * @return BEM3D_SUCCESS on success.
 *
 * The output file is a GtsSurface file followed by the element data 
 * as follows:
 *
 * number of elements
 *
 * element data as specified in ::bem3d_element_write, one element per line
 */

gint bem3d_mesh_write(BEM3DMesh *m, FILE *f)

{
  guint n ;
  gpointer data[8];
  GHashTable * vindex, * eindex, * findex ;
  GtsSurfaceStats stats;
  GtsSurface *s ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  s = GTS_SURFACE(m) ;
  data[0] = f ;
  data[1] = &n ;
  data[2] = vindex = g_hash_table_new (NULL, NULL) ;
  data[3] = eindex = g_hash_table_new (NULL, NULL) ;
  data[4] = findex = g_hash_table_new (NULL, NULL) ;

  gts_surface_stats (s, &stats) ;
  fprintf (f, "%u %u %u", 
	   stats.edges_per_vertex.n, 
	   stats.faces_per_edge.n, 
	   stats.n_faces);
  if (GTS_OBJECT (s)->klass->write)
    (*GTS_OBJECT (s)->klass->write) (GTS_OBJECT (s), f) ;
  fputc ('\n', f) ;
  n = 0 ;
  gts_surface_foreach_vertex (s, (GtsFunc) write_vertex, data) ;
  n = 0 ;
  if (GTS_POINT_CLASS (s->vertex_class)->binary)
    fputc ('\n', f) ;
  gts_surface_foreach_edge (s, (GtsFunc) write_edge, data) ;
  n = 0 ;
  gts_surface_foreach_face (s, (GtsFunc) write_face, data) ;

  /*added code for writing elements*/
  fprintf(f, "%u\n", bem3d_mesh_element_number(m)) ;
  g_hash_table_foreach(m->e, (GHFunc) write_element, data) ;

  g_hash_table_destroy (vindex);
  g_hash_table_destroy (eindex);
  g_hash_table_destroy (findex);

  return BEM3D_SUCCESS ;
}

/** 
 * Read a BEM3DMesh from input file. This is a modified version of
 * Stephane Popinet's original gts_surface_read
 * 
 * @param m BEM3DMesh to read
 * @param f GtsFile to read from
 * 
 * @return 0 on success, line number where error was encountered
 * otherwise
 */

guint bem3d_mesh_read(BEM3DMesh *m, GtsFile * f)

{
  GtsVertex ** vertices;
  GtsEdge ** edges;
  GtsFace ** faces ;
  GtsSurface *surface ;  
  guint n, nv, ne, nf, ns, nlm ;
  gdouble x ;
  BEM3DShapeFunc shf, cpf ;
  BEM3DElement *e ;
  
  /*everything from here to `new code' is Stephane Popinet's*/

  g_return_val_if_fail (m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;  
  g_return_val_if_fail (f != NULL, BEM3D_NULL_ARGUMENT) ;

  surface = GTS_SURFACE(m) ;

  if (f->type != GTS_INT) {
    gts_file_error (f, "expecting an integer (number of vertices)");
    return f->line;
  }
  nv = atoi (f->token->str);

  gts_file_next_token (f);
  if (f->type != GTS_INT) {
    gts_file_error (f, "expecting an integer (number of edges)");
    return f->line;
  }
  ne = atoi (f->token->str);

  gts_file_next_token (f);
  if (f->type != GTS_INT) {
    gts_file_error (f, "expecting an integer (number of faces)");
    return f->line;
  }
  nf = atoi (f->token->str);
  
  gts_file_next_token (f);
  if (f->type == GTS_STRING) {
    if (f->type != GTS_STRING) {
      gts_file_error (f, "expecting a string (GtsSurfaceClass)");
      return f->line;
    }
    gts_file_next_token (f);
    if (f->type != GTS_STRING) {
      gts_file_error (f, "expecting a string (GtsFaceClass)");
      return f->line;
    }
    gts_file_next_token (f);
    if (f->type != GTS_STRING) {
      gts_file_error (f, "expecting a string (GtsEdgeClass)");
      return f->line;
    }
    gts_file_next_token (f);
    if (f->type != GTS_STRING) {
      gts_file_error (f, "expecting a string (GtsVertexClass)");
      return f->line;
    }
    if (!strcmp (f->token->str, "GtsVertexBinary"))
      GTS_POINT_CLASS (surface->vertex_class)->binary = TRUE;
    else
      gts_file_first_token_after (f, '\n');
  }
  else
    gts_file_first_token_after (f, '\n');

  if (nf <= 0) return BEM3D_SUCCESS ;

  /* allocate nv + 1 just in case nv == 0 */
  vertices = g_malloc ((nv + 1)*sizeof (GtsVertex *));
  edges = g_malloc ((ne + 1)*sizeof (GtsEdge *));
  faces = g_malloc ((nf + 1)*sizeof (GtsFace *));
  
  n = 0;
  while (n < nv && f->type != GTS_ERROR) {
    GtsObject * new_vertex =
      gts_object_new (GTS_OBJECT_CLASS (surface->vertex_class));

    (* GTS_OBJECT_CLASS (surface->vertex_class)->read) (&new_vertex, f);
    if (f->type != GTS_ERROR) {
      if (!GTS_POINT_CLASS (surface->vertex_class)->binary)
	gts_file_first_token_after (f, '\n');
      vertices[n++] = GTS_VERTEX (new_vertex);
    }
    else
      gts_object_destroy (new_vertex);
  }
  if (f->type == GTS_ERROR)
    nv = n;
  if (GTS_POINT_CLASS (surface->vertex_class)->binary)
    gts_file_first_token_after (f, '\n');

  n = 0;
  while (n < ne && f->type != GTS_ERROR) {
    guint p1, p2;

    if (f->type != GTS_INT)
      gts_file_error (f, "expecting an integer (first vertex index)");
    else {
      p1 = atoi (f->token->str);
      if (p1 == 0 || p1 > nv)
	gts_file_error (f, "vertex index `%d' is out of range `[1,%d]'", 
			p1, nv);
      else {
	gts_file_next_token (f);
	if (f->type != GTS_INT)
	  gts_file_error (f, "expecting an integer (second vertex index)");
	else {
	  p2 = atoi (f->token->str);
	  if (p2 == 0 || p2 > nv)
	    gts_file_error (f, "vertex index `%d' is out of range `[1,%d]'", 
			    p2, nv);
	  else {
	    GtsEdge * new_edge =
	      gts_edge_new (surface->edge_class,
			    vertices[p1 - 1], vertices[p2 - 1]);

	    gts_file_next_token (f);
	    if (f->type != '\n')
	      if (GTS_OBJECT_CLASS (surface->edge_class)->read)
		(*GTS_OBJECT_CLASS (surface->edge_class)->read)
		  ((GtsObject **)(&new_edge), f);
	    gts_file_first_token_after (f, '\n');
	    edges[n++] = new_edge;
	  }
	}
      }
    }
  }
  if (f->type == GTS_ERROR)
    ne = n;

  n = 0;
  while (n < nf && f->type != GTS_ERROR) {
    guint s1, s2, s3;

    if (f->type != GTS_INT)
      gts_file_error (f, "expecting an integer (first edge index)");
    else {
      s1 = atoi (f->token->str);
      if (s1 == 0 || s1 > ne)
	gts_file_error (f, "edge index `%d' is out of range `[1,%d]'", 
			s1, ne);
      else {
	gts_file_next_token (f);
	if (f->type != GTS_INT)
	  gts_file_error (f, "expecting an integer (second edge index)");
	else {
	  s2 = atoi (f->token->str);
	  if (s2 == 0 || s2 > ne)
	    gts_file_error (f, "edge index `%d' is out of range `[1,%d]'", 
			    s2, ne);
	  else {
	    gts_file_next_token (f);
	    if (f->type != GTS_INT)
	      gts_file_error (f, "expecting an integer (third edge index)");
	    else {
	      s3 = atoi (f->token->str);
	      if (s3 == 0 || s3 > ne)
		gts_file_error (f, "edge index `%d' is out of range `[1,%d]'", 
				s3, ne);
	      else {
		GtsFace * new_face = gts_face_new (surface->face_class,
						   edges[s1 - 1],
						   edges[s2 - 1],
						   edges[s3 - 1]);

		gts_file_next_token (f);
		if (f->type != '\n')
		  if (GTS_OBJECT_CLASS (surface->face_class)->read)
		    (*GTS_OBJECT_CLASS (surface->face_class)->read)
		      ((GtsObject **)(&new_face), f);
		gts_file_first_token_after (f, '\n');
		gts_surface_add_face (surface, new_face);
		faces[n++] = new_face ;
	      }
	    }
	  }
	}
      }
    }
  }

  if (f->type == GTS_ERROR) {
    gts_allow_floating_vertices = TRUE;
    while (nv)
      gts_object_destroy (GTS_OBJECT (vertices[nv-- - 1]));
    gts_allow_floating_vertices = FALSE;
  }

  /*new code*/
  if (f->type != GTS_INT)
    gts_file_error (f, "expecting an integer (number of elements)");
  nlm = atoi(f->token->str) ;
  gts_file_first_token_after (f, '\n');
  n = 0;
  while (n < nlm && f->type != GTS_ERROR) {
    guint i, j, k ;
    if (f->type != GTS_INT)
      gts_file_error (f, 
		      "expecting an integer (number of faces on element)");
    i = atoi(f->token->str) ;
    gts_file_next_token(f) ;

    if (f->type != GTS_INT)
      gts_file_error (f, 
		      "expecting an integer (number of vertices on element)");
    j = atoi(f->token->str) ;    
    gts_file_next_token(f) ;

    if (f->type != GTS_INT)
      gts_file_error (f, 
		      "expecting an integer (number of collocation points)");
    k = atoi(f->token->str) ;
    gts_file_next_token(f) ;

    if (f->type != GTS_INT)
      gts_file_error (f, 
		      "expecting an integer (number of sides on element)");
    ns = atoi(f->token->str) ;    
    gts_file_next_token(f) ;

    if ( f->type != GTS_STRING)
      gts_file_error(f, 
		     "expecting a string (name of shape function)") ;
    shf = bem3d_shapefunc_lookup_func(f->token->str) ;
    gts_file_next_token(f) ;

    if ( f->type != GTS_STRING)
      gts_file_error(f, 
		     "expecting a string (name of collocation point function)"
		     ) ;
    cpf = bem3d_shapefunc_lookup_func(f->token->str) ;
    gts_file_next_token(f) ;

    e = bem3d_element_new(bem3d_element_class(), i, j, k, ns, shf, cpf) ;
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: new element generated address %p", __FUNCTION__, e) ;

    for ( i = 0 ; i < e->nf ; i ++ ) {
      if (f->type != GTS_INT)
	gts_file_error (f, "expecting an integer (%dth face index)", i);
      else {
	j = atoi(f->token->str) ;
	bem3d_element_add_face(e, faces[j-1], i) ;
      }
      gts_file_next_token(f) ;
    }
    for ( i = 0 ; i < e->nv ; i ++ ) {
      if (f->type != GTS_INT)
	gts_file_error (f, "expecting an integer (%dth vertex index)", i);
      else {
	j = atoi(f->token->str) ;
	bem3d_element_add_vertex(e, vertices[j-1], i) ;
      }
      gts_file_next_token(f) ;
    }

    if ( k > 0 ) {
      for ( i = 0 ; i < e->nc ; i ++ ) {      
	if (f->type != GTS_INT)
	  gts_file_error (f, 
			  "expecting an integer (%dth collocation point)", i);
	else {
	  j = atoi(f->token->str) ;
	  bem3d_element_add_node(e, vertices[j-1], i) ;
	}
	gts_file_next_token(f) ;
      }
    }

    for ( i = 0 ; i < e->nc ; i ++ ) {
      if ((f->type != GTS_FLOAT) && (f->type != GTS_INT))
	gts_file_error (f, 
			"expecting a float (%dth shape coordinate)", 
			i);
      else {
	x = atof(f->token->str) ;
	bem3d_element_node_xi(e,i) = x ;
      }
      gts_file_next_token(f) ;
      if ((f->type != GTS_FLOAT) && (f->type != GTS_INT))
	gts_file_error (f, 
			"expecting a float (%dth shape coordinate)", 
			i);
      else {
	x = atof(f->token->str) ;
	bem3d_element_node_eta(e,i) = x ;
      }
      gts_file_next_token(f) ;
    }
    
    if ( k > 0 ) {
      for ( i = 0 ; i < e->nv ; i ++ ) {
	if ((f->type != GTS_FLOAT) && (f->type != GTS_INT))
	  gts_file_error (f, 
			  "expecting a float (%dth collocation coordinate)", 
			  i);
	else {
	  x = atof(f->token->str) ;
	  bem3d_element_vertex_xi(e,i) = x ;
	}
	gts_file_next_token(f) ;
	if ((f->type != GTS_FLOAT) && (f->type != GTS_INT))
	  gts_file_error (f, 
			  "expecting a float (%dth collocation coordinate)", 
			  i);
	else {
	  x = atof(f->token->str) ;
	  bem3d_element_vertex_eta(e,i) = x ;
	}
	gts_file_next_token(f) ;
    }
    }

    for ( i = 0 ; i < e->nc ; i ++ ) {
      if (f->type != GTS_INT)
	gts_file_error (f, 
			"expecting an integer (%dth vertex global index)", 
			i);
      else {
	j = atoi(f->token->str) ;
	bem3d_element_set_index(e, i, j) ;
      }
      gts_file_next_token(f) ;
    }

    for ( i = 0 ; i < e->ns ; i ++ ) {
      if (f->type != GTS_INT)
	gts_file_error (f, 
			"expecting an integer (%dth vertex global index)", 
			i);
      else {
	j = atoi(f->token->str) ;
	bem3d_element_set_corner(e, i, j) ;
      }
      gts_file_next_token(f) ;
    }

    bem3d_mesh_add_element(m, e, FALSE) ;
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: new element inserted address %p", __FUNCTION__, e) ;
    n ++ ;
    gts_file_next_token(f) ;
  }

  g_free (vertices);
  g_free (edges);
  g_free (faces);

  if (f->type == GTS_ERROR) return f->line;

  return BEM3D_SUCCESS ;
}

static void write_node(gint i, GtsVertex *v, gpointer *data)

{
  BEM3DMesh *m = data[0] ;
  FILE *f = data[1] ;
  GtsVector n ;

  bem3d_node_normal(m, i, n, BEM3D_AVERAGE_MWSELR) ;

  fprintf(f, "%d %lg %lg %lg %lg %lg %lg\n",
	  i, 
	  GTS_POINT(v)->x, GTS_POINT(v)->y, GTS_POINT(v)->z,
	  n[0], n[1], n[2]) ;

  return ;
}

/** 
 * Write the nodes and normals of a BEM3DMesh, one per line, to output
 * in the format:
 *
 * [node index] x y z nx ny nz
 * 
 * @param m BEM3DMesh to write;
 * @param f file pointer.
 * 
 * @return 0 on success.
 */

gint bem3d_mesh_write_nodes(BEM3DMesh *m, FILE *f)

{
  gpointer data[2] ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  data[0] = m ; data[1] = f ;
  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)write_node, data) ;

  return BEM3D_SUCCESS ;
}

static gint gmsh_element_n_nodes(gint type)

{
  switch (type) {
  default: g_assert_not_reached() ; break ;
  case GMSH_ELEMENT_LINE: return 2 ; break ;
  case GMSH_ELEMENT_TRIANGLE: return 3 ; break ;
  case GMSH_ELEMENT_QUADRANGLE: return 4 ; break ;
  case GMSH_ELEMENT_TETRAHEDRON: return 4 ; break ;
  case GMSH_ELEMENT_HEXAHEDRON: return 8 ; break ;
  case GMSH_ELEMENT_PRISM: return 6 ; break ;
  case GMSH_ELEMENT_PYRAMID: return 5 ; break ;
  case GMSH_ELEMENT_LINE_2: return 2+1 ; break ;
  case GMSH_ELEMENT_TRIANGLE_2: return 3+3 ; break ;
  case GMSH_ELEMENT_QUADRANGLE_2: return 4+4+1 ; break ;
  case GMSH_ELEMENT_TETRAHEDRON_2: return 4+6 ; break ;
  case GMSH_ELEMENT_HEXAHEDRON_2: return 8+12+6+1 ; break ;
  case GMSH_ELEMENT_PRISM_2: return 6+9+3 ; break ;
  case GMSH_ELEMENT_PYRAMID_2: return 5+8+1 ; break ;
  case GMSH_ELEMENT_POINT: return 1 ; break ;
  case GMSH_ELEMENT_TRIANGLE_3: return 3+6+1 ; break ;
  case GMSH_ELEMENT_LINE_3: return 2+2 ; break ;
  }
  return 0 ;
}

static gint gmsh_add_element(BEM3DMesh *m, GHashTable *h,
			     gint elem, gint *nodes, gint nn) 

{
  GtsEdge *e[32] ;
  GtsVertex *w[32] ;
  BEM3DElement *el ;
  BEM3DElementBuildFunc build_func ;
  gint i, nb ;

  switch (elem) {
  case GMSH_ELEMENT_LINE:
    return 0 ;
    break ;
  case GMSH_ELEMENT_LINE_2:
    return 0 ;
    break ;
  case GMSH_ELEMENT_TRIANGLE: 
    g_assert(nn == 3) ;
    w[0] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[0])) ;
    w[1] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[1])) ;
    w[2] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[2])) ;
    g_assert(w[0] != NULL && w[1] != NULL && w[2] != NULL) ;
    build_func = bem3d_element_build_t1 ;
    nb = 3 ;
    break ;
  case GMSH_ELEMENT_TRIANGLE_2:
    g_assert(nn == 6) ;
    w[0] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[0])) ;
    w[1] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[3])) ;
    w[2] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[1])) ;
    w[3] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[4])) ;
    w[4] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[2])) ;
    w[5] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[5])) ;
    build_func = bem3d_element_build_t2 ;
    nb = 6 ;
    break ;
  case GMSH_ELEMENT_TRIANGLE_3:
    g_assert(nn == 10) ;
    w[0] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[0])) ;
    w[1] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[3])) ;
    w[2] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[4])) ;
    w[3] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[1])) ;
    w[4] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[5])) ;
    w[5] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[6])) ;
    w[6] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[2])) ;
    w[7] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[7])) ;
    w[8] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[8])) ;
    w[9] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[9])) ;
    build_func = bem3d_element_build_t3 ;
    nb = 9 ;
    break ;
  case GMSH_ELEMENT_QUADRANGLE:
    g_assert(nn == 4) ;
    w[0] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[0])) ;
    w[1] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[1])) ;
    w[2] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[2])) ;
    w[3] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[3])) ;
    build_func = bem3d_element_build_q1 ;
    nb = 4 ;
    break ;
  case GMSH_ELEMENT_QUADRANGLE_2:
    g_assert(nn == 9) ;
    w[0] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[0])) ;
    w[1] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[4])) ;
    w[2] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[1])) ;
    w[3] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[5])) ;
    w[4] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[2])) ;
    w[5] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[6])) ;
    w[6] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[3])) ;
    w[7] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[7])) ;
    w[8] = g_hash_table_lookup(h,GINT_TO_POINTER(nodes[8])) ;
    build_func = bem3d_element_build_q2 ;
    nb = 8 ;
    break ;
  default: 
    g_debug("%s: unimplemented GMSH element type %d", 
	    __FUNCTION__, elem) ;
    return 0 ; break ;
  }

  /*link boundary nodes*/
  for ( i = 0 ; i < nb-1 ; i ++ )
    if ( (e[i] = GTS_EDGE(gts_vertices_are_connected(w[i], 
						     w[i+1]))) == NULL )
      e[i] = gts_edge_new(gts_edge_class(), w[i], w[i+1]) ;
  if ( (e[nb-1] = GTS_EDGE(gts_vertices_are_connected(w[nb-1], 
						      w[0]))) == NULL )
    e[nb-1] = gts_edge_new(gts_edge_class(), w[nb-1], w[0]) ;

  el = build_func(e, w) ;
  bem3d_mesh_add_element(m, el, FALSE) ;

  return 0 ;
}

static gint gmsh_read_element1(FILE *f, gint *i, gint *e,
			      gint *rp, gint *re,
			      gint *nodes, gint *nn)

{
  gint j ;
  
  if ( fscanf(f, "%d", i) != 1) return -1 ;
  if ( fscanf(f, "%d", e) != 1) return -1 ;

  *nn = gmsh_element_n_nodes(*e) ;

  if ( fscanf(f, "%d", rp) != 1) return -1 ;
  if ( fscanf(f, "%d", re) != 1) return -1 ;

  if ( fscanf(f, "%d", &j) != 1 ) return -1 ;

  for ( j = 0 ; j < (*nn) ; j ++ )
    if ( fscanf(f, "%d", &nodes[j]) != 1 ) return -1 ;

  return 0 ;
}

static guint gmsh_read_file1(FILE *f, BEM3DMesh *m)

{
  gchar line[256] ;
  guint lineno = 2 ;
  gint i, j, np, ne, nv, elem, ntags, data[32], rp, re ;
  GHashTable *h ;
  gdouble x, y, z ;

  nv = fscanf(f, "%d", &np) ;
  if ( nv == 0 ) {
    g_debug("%s: error reading number of vertices", __FUNCTION__) ;
    return lineno ;
  }

  lineno ++ ;

  h = g_hash_table_new(NULL, NULL) ;
  for ( i = 0 ; i < np ; i ++ ) {
    nv = fscanf(f, "%d %lg %lg %lg", &j, &x, &y, &z) ;
    if ( nv != 4 ) {
      g_debug("%s: error reading vertex at line %u", __FUNCTION__, lineno) ;
      return lineno ;    
    }
    if ( j <= 0 ) {
      g_debug("%s: vertex index %d out of range at line %u", 
	      __FUNCTION__, j, lineno) ;
      return lineno ;    
    }
    g_hash_table_insert(h, GINT_TO_POINTER(j), 
			gts_vertex_new(GTS_SURFACE(m)->vertex_class, 
				       x, y, z)) ;
    lineno ++ ;
  }

  g_debug("%s: %d vertices read", __FUNCTION__, np) ;

  nv = fscanf(f, "%s", line) ; 
  if ( strcmp(line, "$ENDNOD") ) {
    g_debug("%s: no $ENDNOD marker found", __FUNCTION__) ;
    return lineno ;        
  }
  lineno ++ ;

  nv = fscanf(f, "%s", line) ; 
  if ( strcmp(line, "$ELM") ) {
    g_debug("%s: no $ELM marker found", __FUNCTION__) ;
    return lineno ;        
  }
  lineno ++ ;

  nv = fscanf(f, "%d", &ne) ;
  if ( nv == 0 ) {
    g_debug("%s: error reading number of elements", __FUNCTION__) ;
    return lineno ;
  }
  lineno ++ ;

  for ( i = 0 ; i < ne ; i ++ ) {
    if ( gmsh_read_element1(f, &j, &elem, &rp, &re, data, &ntags) != 0 ) {
      g_debug("%s: error reading element at line %u", __FUNCTION__, lineno) ;
      return lineno ;
    }
    if ( gmsh_add_element(m, h, elem, data, ntags) != 0 ) {
      g_debug("%s: error adding element at line %u", __FUNCTION__, lineno) ;
      return lineno ;
    }
    lineno ++ ;
  }

  g_hash_table_destroy(h) ;  

  return 0 ;
}

static gint gmsh_read_element2(FILE *f, gint *i, gint *e,
			      gint *tags,
			      gint *nodes, gint *nn)

{
  gint j, nt ;
  
  if ( fscanf(f, "%d", i) != 1) return -1 ;
  if ( fscanf(f, "%d", e) != 1) return -1 ;

  g_debug("%s: tag %d; element %d", __FUNCTION__, *i, *e) ;
  
  *nn = gmsh_element_n_nodes(*e) ;
  g_debug("%s: %d nodes", __FUNCTION__, *nn) ;

  if ( fscanf(f, "%d", &nt) != 1) return -1 ;

  for ( j = 0 ; j < nt ; j ++ ) 
    if ( fscanf(f, "%d", &(tags[j])) != 1) return -1 ;

  for ( j = 0 ; j < (*nn) ; j ++ )
    if ( fscanf(f, "%d", &nodes[j]) != 1 ) return -1 ;

  return 0 ;
}

static gint gmsh_read_element4(FILE *f, gint *i, gint elem,
			       gint *nodes, gint *nn)

{
  gint j ;
  
  if ( fscanf(f, "%d", i) != 1) return -1 ;

  g_debug("%s: tag %d; element %d", __FUNCTION__, *i, elem) ;
  
  *nn = gmsh_element_n_nodes(elem) ;
  g_debug("%s: %d nodes", __FUNCTION__, *nn) ;

  for ( j = 0 ; j < (*nn) ; j ++ )
    if ( fscanf(f, "%d", &nodes[j]) != 1 ) return -1 ;
  g_debug("%s: nodes: %d, ..., %d", __FUNCTION__, nodes[0], nodes[*nn-1]) ;

  return 0 ;
}

static gint search_line(FILE *f, gchar *string, guint lineno, gint *nv)

{
  gchar line[256] ;

  *nv = fscanf(f, "%[^\n]s", line) ;
  while ( strncmp(string, line, strlen(string)) != 0 ) {
    lineno ++ ;
    fscanf(f, "%*c") ;
    if ( ( *nv = fscanf(f, "%[^\n]s", line) ) == EOF ) return lineno ;
  }
  fscanf(f, "%*c") ;
  
  return lineno ;
}

static guint gmsh_read_file2_0(FILE *f, BEM3DMesh *m, guint lineno)

{
  gchar line[256] ;
  gint i, j, np, ne, nv, elem, ntags, data[32], tags[32] ;
  GHashTable *h ;
  gdouble x, y, z ;
  
  nv = fscanf(f, "%s", line) ;
  g_debug("%s: %s", __FUNCTION__, line) ;
  if ( strcmp(line, "$EndMeshFormat") || nv == 0) {
    g_debug("%s: no $EndMeshFormat marker found", __FUNCTION__) ;
    return lineno ;        
  }
  lineno ++ ;

  nv = fscanf(f, "%s", line) ;
  g_debug("%s: %s", __FUNCTION__, line) ;
  if ( strcmp(line, "$Nodes") || nv == 0) {
    g_debug("%s: no $Nodes marker found", __FUNCTION__) ;
    return lineno ;        
  }
  lineno ++ ;

  nv = fscanf(f, "%d", &np) ;
  if ( nv == 0 ) {
    g_debug("%s: error reading number of vertices", __FUNCTION__) ;
    return lineno ;
  }

  lineno ++ ;

  h = g_hash_table_new(NULL, NULL) ;
  for ( i = 0 ; i < np ; i ++ ) {
    nv = fscanf(f, "%d %lg %lg %lg", &j, &x, &y, &z) ;
    if ( nv != 4 ) {
      g_debug("%s: error reading vertex at line %u", __FUNCTION__, lineno) ;
      return lineno ;    
    }
    if ( j <= 0 ) {
      g_debug("%s: vertex index %d out of range at line %u", 
	      __FUNCTION__, j, lineno) ;
      return lineno ;    
    }
    g_hash_table_insert(h, GINT_TO_POINTER(j), 
			gts_vertex_new(GTS_SURFACE(m)->vertex_class, 
				       x, y, z)) ;
    lineno ++ ;
  }

  g_debug("%s: %d vertices read", __FUNCTION__, np) ;

  nv = fscanf(f, "%s", line) ; 
  g_debug("%s: %s", __FUNCTION__, line) ;
  if ( strcmp(line, "$EndNodes") || nv == 0) {
    g_debug("%s: no $EndNodes marker found", __FUNCTION__) ;
    return lineno ;        
  }
  lineno ++ ;

  nv = fscanf(f, "%s", line) ;
  g_debug("%s: %s", __FUNCTION__, line) ;
  if ( strcmp(line, "$Elements") ) {
    g_debug("%s: no $Elements marker found", __FUNCTION__) ;
    return lineno ;        
  }
  lineno ++ ;

  nv = fscanf(f, "%d", &ne) ;
  if ( nv == 0 ) {
    g_debug("%s: error reading number of elements", __FUNCTION__) ;
    return lineno ;
  }
  lineno ++ ;

  for ( i = 0 ; i < ne ; i ++ ) {
    if ( gmsh_read_element2(f, &j, &elem, tags, data, &ntags) != 0 ) {
      g_debug("%s: error reading element at line %u", __FUNCTION__, lineno) ;
      return lineno ;
    }
    if ( gmsh_add_element(m, h, elem, data, ntags) != 0 ) {
      g_debug("%s: error adding element at line %u", __FUNCTION__, lineno) ;
      return lineno ;
    }
    lineno ++ ;
  }

  g_hash_table_destroy(h) ;  

  return 0 ;
}

static guint gmsh_read_file4_0(FILE *f, BEM3DMesh *m, guint lineno)

{
  gint i, j, k, ne, np, nv, elem, ntags, data[32] ;
  gint n_ent, n_nodes, n_elem, imin, imax, dim, etag, nnblock, *tags ;
  gboolean parametric ;
  GHashTable *h ;
  gdouble x, y, z ;
  
  lineno = search_line(f, "$EndMeshFormat", lineno, &nv) ;
  if ( nv == EOF ) {
    g_debug("%s: no $EndMeshFormat marker found", __FUNCTION__) ;
    return lineno ;
  }
  
  lineno = search_line(f, "$Nodes", lineno, &nv) ;
  if ( nv == EOF ) {
    g_debug("%s: no $Nodes marker found", __FUNCTION__) ;
    return lineno ;
  }

  nv = fscanf(f, "%d %d %d %d", &n_ent, &n_nodes, &imin, &imax) ;
  lineno ++ ;
  g_debug("%s: line %u; %d entities; %d nodes; min %d; max %d\n",
	  __FUNCTION__, lineno, n_ent, n_nodes, imin, imax) ;
  if ( nv == 0 ) {
    g_debug("%s: error reading number of vertices", __FUNCTION__) ;
    return lineno ;
  }

  h = g_hash_table_new(NULL, NULL) ;
  tags = (gint *)g_malloc(MAX(n_nodes, 32)*sizeof(gint)) ;
  for ( i = np = 0 ; i < n_ent ; i ++ ) {
    nv = fscanf(f, "%d %d %d %d", &dim, &etag, &parametric, &nnblock) ;

    lineno ++ ;
    g_debug("%s: line %u; dimension %d; tag %d; parametric %d; "
	    "nnodes block %d\n",
	    __FUNCTION__, lineno, dim, etag, parametric, nnblock) ;

    for ( k = 0 ; k < nnblock ; k ++ ) {
      nv = fscanf(f, "%d", &j) ;
      if ( nv == 0 ) {
	g_debug("%s: line %u; error reading vertex tag",
		__FUNCTION__, lineno) ;
	return lineno ;
      }
      lineno ++ ;
      tags[k] = j ;
    }
    
    for ( k = 0 ; k < nnblock ; k ++ ) {
      nv = fscanf(f, "%lg %lg %lg", &x, &y, &z) ;
      lineno ++ ;
      if ( nv != 3 ) {
	g_debug("%s: error reading vertex at line %u", __FUNCTION__, lineno) ;
      }
      g_hash_table_insert(h, GINT_TO_POINTER(tags[k]), 
			  gts_vertex_new(GTS_SURFACE(m)->vertex_class, 
					 x, y, z)) ;
      np ++ ;
    }
  }
  
  g_debug("%s: %d vertices read", __FUNCTION__, np) ;

  lineno = search_line(f, "$EndNodes", lineno, &nv) ;
  if ( nv == EOF ) {
    g_debug("%s: line %u; no $EndNodes marker found", __FUNCTION__, lineno) ;
    return lineno ;
  }

  lineno = search_line(f, "$Elements", lineno, &nv) ;
  if ( nv == EOF ) {
    g_debug("%s: no $Elements marker found", __FUNCTION__) ;
    return lineno ;
  }

  nv = fscanf(f, "%d %d %d %d", &n_ent, &n_elem, &imin, &imax) ;
  lineno ++ ;
  g_debug("%s: line %u; %d entities; %d elements; min %d; max %d\n",
	  __FUNCTION__, lineno, n_ent, n_elem, imin, imax) ;

  for ( i = ne = 0 ; i < n_ent ; i ++ ) {
    nv = fscanf(f, "%d %d %d %d", &dim, &etag, &elem, &nnblock) ;
    if ( nv != 4 )
      g_debug("%s: error reading element data at line %u",
	      __FUNCTION__, lineno) ;
    lineno ++ ;
    g_debug("%s: line %u; dimension %d; tag %d; type %d; "
	    "%d elements\n",
	    __FUNCTION__, lineno, dim, etag, elem, nnblock) ;
    for ( j = 0 ; j < nnblock ; j ++ ) {
      if ( gmsh_read_element4(f, &k, elem, data, &ntags) != 0 ) {
	g_debug("%s: error reading element at line %u", __FUNCTION__, lineno) ;
	return lineno ;
      }
      lineno ++ ;
      if ( gmsh_add_element(m, h, elem, data, ntags) != 0 ) {
	g_debug("%s: error adding element at line %u", __FUNCTION__, lineno) ;
	return lineno ;
      }
    }
    ne ++ ;
  }

  lineno ++ ;
  
  lineno = search_line(f, "$EndElements", lineno, &nv) ;
  if ( nv == EOF ) {
    g_debug("%s: no $EndElements marker found", __FUNCTION__) ;
    return lineno ;
  }

  g_debug("%s: line %u; %d elements read", __FUNCTION__, lineno, ne) ;

  g_hash_table_destroy(h) ;  

  return 0 ;
}

static guint gmsh_read_file2(FILE *f, BEM3DMesh *m)

{
  guint lineno = 2 ;
  gint ft, ds, nv ;
  gdouble version ;

  if ( (nv = fscanf(f, "%lg %d %d", &version, &ft, &ds) ) != 3 ) {
    g_debug("%s: error reading mesh format information", 
	    __FUNCTION__) ;
    return lineno ;
  }
  if ( version < 2.0 || version > 4.1 ) {
    g_debug("%s: file version (%lg) should be 2.0",  
	    __FUNCTION__, version) ;
    return lineno ;
  }
  lineno ++ ;

  g_debug("%s: file format %lg, file type %d, data size %d",  
	  __FUNCTION__, version, ft, ds) ;

  if ( version >= 2.0 && version < 4.0 )
    return gmsh_read_file2_0(f, m, lineno) ;
  
  if ( version >= 4.0 )
    return gmsh_read_file4_0(f, m, lineno) ;

  g_assert_not_reached() ;
  
  return 0 ;
}

/** 
 * Read a GMSH 1.0 or 2.0 mesh file to a ::BEM3DMesh.
 * 
 * @param m a ::BEM3DMesh;
 * @param f input file pointer.
 * 
 * @return BEM3D_SUCCESS on success.
 */

guint bem3d_gmsh_read(BEM3DMesh *m, FILE *f)

{
  gchar line[256] ;
  gint gmsh_format = 0, nf = 1 ;
  
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  while ( gmsh_format == 0 && nf != 0 ) {
    if ( (nf = fscanf(f, "%s", line)) != 0 ) {
      if ( strncmp(line, "$NOD", 4) == 0 ) gmsh_format = 1 ;
      if ( strncmp(line, "$MeshFormat", 11) == 0 ) gmsh_format = 2 ;
    }
  }

  if ( gmsh_format == 0 ) {
    g_warning("%s: unrecognized gmsh file format", __FUNCTION__) ;
    return BEM3D_UNKNOWN_FORMAT ;
  }

  g_debug("%s: GMSH file format: %d", __FUNCTION__, gmsh_format) ;

  switch ( gmsh_format ) {
  default: g_assert_not_reached() ; break ;
  case 1: return gmsh_read_file1(f, m) ; break ;
  case 2: return gmsh_read_file2(f, m) ; break ;
  }

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */

