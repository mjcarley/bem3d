/* gmsh.c
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
 * @defgroup gmsh GMSH meshing and visualization
 *
 * Functions for reading and writing GMSH files. GMSH is a meshing and
 * post-processing program available from http://geuz.org/gmsh/ and
 * described in Geuzaine, C. and Remacle, J.-F., `Gmsh: a
 * three-dimensional finite element mesh generator with built-in pre-
 * and post-processing facilities', International Journal for Numerical
 * Methods in Engineering, 2009, http://dx.doi.org/10.1002/nme.2579
 *
 * @{
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

#define _GMSH_DATA_WIDTH    16
#define _GMSH_DATA_FILE      0
#define _GMSH_DATA_DATA      1
#define _GMSH_DATA_INDEX     2
#define _GMSH_DATA_INDICES   3
#define _GMSH_DATA_N_NODES   4
#define _GMSH_DATA_MODE      5
#define _GMSH_DATA_HASH      6
#define _GMSH_DATA_ETYPE     7
#define _GMSH_DATA_ETYPES    8
#define _GMSH_DATA_N_COM     9
#define _GMSH_DATA_OFFSET   10

static gint _gmsh_element_tag(BEM3DElement *e)

{
  switch ( bem3d_element_node_number(e) ) {
  default: g_assert_not_reached() ; break ;
  case 3: return 2 ; break ; /*linear triangle*/
  }
  

  return -1 ;
}

static gint _write_element_pos(BEM3DElement *e, gpointer data[])

{
  FILE *fp = data[_GMSH_DATA_FILE] ;
  BEM3DMeshData *f = data[_GMSH_DATA_DATA] ;
  gint k = *(gint *)data[_GMSH_DATA_INDEX] ;
  bem3d_gmsh_mode_t mode = *(bem3d_gmsh_mode_t *)data[_GMSH_DATA_MODE] ;
  gint nf, nn, gt, *indices_g, indices_d[64], i ;
  gint indices_st[] = {0, 1, 2} ;
  gint indices_st2[] = {0, 1, 2, 3, 4, 5} ;
  /* gint indices_st3[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} ; */
  gint indices_sq[] = {0, 1, 2, 3} ;
  gint indices_sq2[] = {0, 1, 2, 3, 4, 5, 6, 7, 8} ;
  gchar *gmsh_types[] = {
    "SP", "VP", "TP", 
    "SL", "VL", "TL", 
    "ST", "VT", "TT", 
    "SQ", "VQ", "TQ", 
    "SS", "VS", "TS", 
    "SH", "VH", "TH", 
    "SI", "VI", "TI",
    "SY", "VY", "TY",
    "SL2", "VL2", "TL2", 
    "ST2", "VT2", "TT2", 
    "SQ2", "VQ2", "TQ2", 
    "SS2", "VS2", "TS2", 
    "SH2", "VH2", "TH2", 
    "SI2", "VI2", "TI2",
    "SY2", "VY2", "TY2"} ;

  gt = nf = nn = 0 ; indices_g = NULL ;
  switch ( mode ) {
  default: g_assert_not_reached() ; break ;
  case BEM3D_GMSH_SCALAR: nf = 1 ; gt = 0 ; break ;
  case BEM3D_GMSH_VECTOR: nf = 3 ; gt = 1 ; break ;
  case BEM3D_GMSH_TENSOR: nf = 9 ; gt = 2 ; break ;
  }

  switch ( (nn = bem3d_element_vertex_number(e)) ) {
  default: g_assert_not_reached() ; break ;
  case 3: gt += 6 ; indices_g = indices_st ; break ;
  case 4: gt += 9 ; indices_g = indices_sq ; break ;
  case 6: gt += 27 ; indices_g = indices_st2 ; break ;
  case 9: gt += 30 ; indices_g = indices_sq2 ; break ;
  case 10: 
    g_error("%s: cubic triangles not implemented in .pos format",
	    __FUNCTION__) ; 
    break ;
  }

  if ( bem3d_element_node_number(e) == bem3d_element_vertex_number(e) ) {
    for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ )
      indices_d[i] = bem3d_element_global_index(e,i) ;
  } else {
    switch ( bem3d_element_node_number(e) ) {
    default: g_assert_not_reached() ; break ;
    case 1:
      for ( i = 0 ; i < bem3d_element_vertex_number(e) ; i ++ ) 
	indices_d[i] = bem3d_element_global_index(e,0) ;
      break ;
    }
  }

  bem3d_element_write_pos(e, f, nn, indices_g, indices_d,
			  gmsh_types[gt], k, nf, fp) ; 

  return BEM3D_SUCCESS ;
}

/** 
 * Write an element in the GMSH .pos visualization format, including
 * the corresponding nodal data.
 * 
 * @param e element;
 * @param f data block for the mesh containing e;
 * @param nv number of vertices in GMSH element;
 * @param indices_g array of indices specifying the order in which to 
 *                write the element vertices; 
 * @param indices_d array of indices for the data for each vertex;
 * @param ename element name for GMSH parser;
 * @param k field of f to write;
 * @param nf number of data points to write per node of element;
 * @param fp file pointer for output.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_element_write_pos(BEM3DElement *e, BEM3DMeshData *f,
			     gint nv, 
			     gint *indices_g, gint *indices_d,
			     gchar *ename,
			     gint k, gint nf, FILE *fp)

{
  GtsVertex *v ;
  gdouble *x ;
  gint i, j, n ;

  g_return_val_if_fail(e != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_ELEMENT(e), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(indices_g != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(indices_d != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(ename != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(k+nf <= bem3d_mesh_data_element_number(f), 
		       BEM3D_ARGUMENT_OUT_OF_RANGE) ;

  fprintf(fp, "%s(", ename) ;
  i = indices_g[0] ; v = e->v[i] ;
  g_debug("%s: node 0 index %d", __FUNCTION__, i) ;
  fprintf(fp, "%lg,%lg,%lg", 
	  GTS_POINT(v)->x, GTS_POINT(v)->y, GTS_POINT(v)->z) ;

  for ( j = 1 ; j < nv ; j ++ ) {
    i = indices_g[j] ; v = e->v[i] ;
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: node %d index %d", __FUNCTION__, j, i) ;
    fprintf(fp, ",%lg,%lg,%lg", 
	    GTS_POINT(v)->x, GTS_POINT(v)->y, GTS_POINT(v)->z) ;
  }

  fprintf(fp, "){") ;
  i = indices_g[0] ; 
  if ( (x = bem3d_mesh_data_get(f, indices_d[0])) == NULL )
    g_error("%s: data for node %d not found in data block", 
	    __FUNCTION__, indices_d[0]) ;
  fprintf(fp, "%lg", x[k]) ;
  for ( n = 1 ; n < nf ; n ++ ) fprintf(fp, ",%lg", x[k+n]) ;

  for ( j = 1 ; j < nv ; j ++ ) {
    if ( (x = bem3d_mesh_data_get(f, indices_d[j])) == NULL )
    g_error("%s: data for node %d not found in data block", 
	    __FUNCTION__, indices_d[j]) ;
    for ( n = 0 ; n < nf ; n ++ ) fprintf(fp, ",%lg", x[k+n]) ;
  }
  fprintf(fp, "};\n") ;

  return BEM3D_SUCCESS ;
}

/** 
 * Write a BEM3DMesh and corresponding data to a GMSH .pos file for
 * visualization.
 * 
 * @param m mesh to write;
 * @param f data block;
 * @param k index of data block field to write;
 * @param view title of view or NULL;
 * @param mode BEM3D_GMSH_SCALAR | BEM3D_GMSH_VECTOR | BEM3D_GMSH_TENSOR, if
 * a vector or tensor mode is chosen, the 3 or 9 elements are written
 * starting from the indicated field index \a k. The data block must
 * be wide enough to accomodate this.
 * @param fp file pointer for output.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_mesh_write_pos(BEM3DMesh *m, BEM3DMeshData *f, gint k,
			  gchar *view, bem3d_gmsh_mode_t mode,
			  FILE *fp)

{
  gpointer data[_GMSH_DATA_WIDTH] ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(fp != NULL, BEM3D_NULL_ARGUMENT) ;

  data[_GMSH_DATA_FILE] = fp ; 
  data[_GMSH_DATA_DATA] = f ; 
  data[_GMSH_DATA_INDEX] = &k ;
  data[_GMSH_DATA_MODE] = &mode ;
  if ( view != NULL ) 
    fprintf(fp, "View \"%s\" {\n", view) ;
  else
    fprintf(fp, "View \"BEM3D\" {\n") ;

  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)_write_element_pos, data) ;

  fprintf(fp, "} ;\n") ;

  return BEM3D_SUCCESS ;
}

gint bem3d_edge_write_pos(BEM3DEdge *edge, gchar *view, FILE *f)

{
  gint i ;

  g_return_val_if_fail(edge != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_EDGE(edge), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( view == NULL )
    fprintf(f, "View \"BEM3D Edge\" {\n") ;
  else
    fprintf(f, "View \"%s\" {\n", view) ;

  for ( i = 0 ; i < bem3d_edge_node_number(edge)-1 ; i ++ ) {
    fprintf(f, "SL(%g,%g,%g,%g,%g,%g){0,0,0};\n",
	    GTS_POINT(bem3d_edge_vertex(edge,i))->x,
	    GTS_POINT(bem3d_edge_vertex(edge,i))->y,
	    GTS_POINT(bem3d_edge_vertex(edge,i))->z,
	    GTS_POINT(bem3d_edge_vertex(edge,i+1))->x,
	    GTS_POINT(bem3d_edge_vertex(edge,i+1))->y,
	    GTS_POINT(bem3d_edge_vertex(edge,i+1))->z) ;
  }

  fprintf(f, "} ;\n") ;

  return BEM3D_SUCCESS ;
}

static void write_node_msh_4_1(gint i, GtsVertex *v, gpointer *data)

{
  FILE *fp = data[_GMSH_DATA_FILE] ; 
  gint *idx = data[_GMSH_DATA_INDEX] ;
  GHashTable *h = data[_GMSH_DATA_HASH] ;

  g_hash_table_insert(h, GINT_TO_POINTER(i), GINT_TO_POINTER(*idx)) ;

  fprintf(fp, "%lg %lg %lg\n",
	  GTS_POINT(v)->x, GTS_POINT(v)->y, GTS_POINT(v)->z) ;

  (*idx) ++ ;
  
  return ;
}

static void write_element_msh_4_1(BEM3DElement *e, gpointer *data)

{
  FILE *fp = data[_GMSH_DATA_FILE] ; 
  gint *idx = data[_GMSH_DATA_INDEX] ;
  GHashTable *h = data[_GMSH_DATA_HASH] ;
  gint etype = *((gint *)data[_GMSH_DATA_ETYPE]) ;
  gint offset = *((gint *)data[_GMSH_DATA_OFFSET]) ;
  gint i, j ;
  
  (*idx) ++ ;
  /* 2D entity (surface) 1, element type, one element*/
  if ( _gmsh_element_tag(e) != etype ) return ;
  fprintf(fp, "%d", *idx) ;
  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    j = bem3d_element_global_index(e,i) ;
    j = GPOINTER_TO_INT(g_hash_table_lookup(h, GINT_TO_POINTER(j))) ;
    if ( j != 0 ) fprintf(fp, " %d", offset+j) ;
  }
  fprintf(fp, "\n") ;
  

  return ;
}

static void count_element_types(BEM3DElement *e, gpointer *data)

{
  gint *etypes = data[_GMSH_DATA_ETYPES] ;
  gint et ;

  et = _gmsh_element_tag(e) ;
  etypes[et] ++ ;
  
  return ;
}

static gint write_data_msh(gint i, gdouble *d, gint nf, gpointer *data)

{
  gint k = *((gint *)data[_GMSH_DATA_INDEX]) ;
  FILE *fp = data[_GMSH_DATA_FILE] ;
  GHashTable *h = data[_GMSH_DATA_HASH] ;
  gint ncom = *((gint *)data[_GMSH_DATA_N_COM]) ;
  gint offset = *((gint *)data[_GMSH_DATA_OFFSET]) ;
  gint j ;

  j = GPOINTER_TO_INT(g_hash_table_lookup(h, GINT_TO_POINTER(i))) ;

  if ( j == 0 ) return 0 ;

  fprintf(fp, "%d", offset+j) ;
  for ( j = 0 ; j < ncom ; j ++ ) {
    fprintf(fp, " %lg", d[k+j]) ;
  }
  fprintf(fp, "\n") ;
  
  return 0 ;
}

static gint bem3d_mesh_write_msh_4_1(BEM3DMesh *m, BEM3DMeshData *f,
				     gint k, gchar *view, gdouble t,
				     bem3d_gmsh_mode_t mode,
				     gint offset, FILE *fp)

{
  gpointer data[_GMSH_DATA_WIDTH] ;
  gint nnodes, nelem, i, n_etypes[16] = {0}, etype, ncom ;
  GHashTable *h ;
  
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(fp != NULL, BEM3D_NULL_ARGUMENT) ;

  data[_GMSH_DATA_FILE] = fp ; 
  data[_GMSH_DATA_DATA] = f ; 
  data[_GMSH_DATA_INDEX] = &i ;
  data[_GMSH_DATA_MODE] = &mode ;
  data[_GMSH_DATA_HASH] = h = g_hash_table_new(NULL, NULL) ;
  data[_GMSH_DATA_ETYPE] = &etype ;
  data[_GMSH_DATA_ETYPES] = n_etypes ;
  data[_GMSH_DATA_N_COM] = &ncom ;
  data[_GMSH_DATA_OFFSET] = &offset ;
  
  /*this really only works for isoparametric elements*/
  nnodes = bem3d_mesh_node_number(m) ;
  nelem = bem3d_mesh_element_number(m) ;

  /*we need to know how many elements of each type are to be written*/
  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)count_element_types, data) ;
  
  fprintf(fp, "$MeshFormat\n") ;
  fprintf(fp, "4.1 0 %lu\n", sizeof(gdouble)) ;
  fprintf(fp, "$EndMeshFormat\n") ;
  fprintf(fp, "$Nodes\n") ;
  /*one entity block, number of nodes, tags from 1 to nnodes*/
  fprintf(fp, "1 %d %d %d\n", nnodes, offset+1, offset+nnodes) ;
  /*2D entity (surface) 1, no parametric coordinates, nnodes */
  fprintf(fp, "2 1 0 %d\n", nnodes) ;
  /*node tags*/
  for ( i = 1 ; i <= nnodes ; i ++ ) fprintf(fp, "%d\n", offset+i) ;
  i = 1 ;
  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)write_node_msh_4_1, data) ;
  fprintf(fp, "$EndNodes\n") ;

  fprintf(fp, "$Elements\n") ;
  /*one entity block, nelem elements, tags from 1 to nelem*/
  fprintf(fp, "1 %d 1 %d\n", nelem, nelem) ;
  i = 0 ;
  for ( etype = 0 ; etype < 16 ; etype ++ ) {
    if ( n_etypes[etype] != 0 ) {
      fprintf(fp, "2 1 %d %d\n", etype, n_etypes[etype]) ;
      
      bem3d_mesh_foreach_element(m, (BEM3DElementFunc)write_element_msh_4_1,
				 data) ;
    }
  }
  
  fprintf(fp, "$EndElements\n") ;

  if ( f == NULL ) return BEM3D_SUCCESS ; /*no data to write*/

  ncom = ( (mode == BEM3D_GMSH_SCALAR) ? 1 :
	   ( mode == BEM3D_GMSH_VECTOR ? 3 : 9)) ;
  
  fprintf(fp, "$NodeData\n") ;
  fprintf(fp, "1\n") ;
  if ( view != NULL ) 
    fprintf(fp, "\"%s\"\n", view) ;
  else
    fprintf(fp, "\"BEM3D\"\n") ;
  fprintf(fp,
	  "1\n"
	  "%lg\n"
	  "3\n"
	  "0\n"
	  "%d\n"
	  "%d\n",
	  t, ncom, nnodes) ;

  data[_GMSH_DATA_INDEX] = &k ;
  bem3d_mesh_data_foreach(f, (BEM3DMeshDataEntryFunc)write_data_msh, data) ;
  
  fprintf(fp, "$EndNodeData\n") ;

  return BEM3D_SUCCESS ;
}


/** 
 * Write a BEM3DMesh and corresponding data to a GMSH .msh file for
 * visualization.
 * 
 * @param m mesh to write;
 * @param f data block;
 * @param k index of data block field to write;
 * @param view title of view or NULL;
 * @param t time for GMSH view;
 * @param mode BEM3D_GMSH_SCALAR | BEM3D_GMSH_VECTOR | BEM3D_GMSH_TENSOR, if
 * a vector or tensor mode is chosen, the 3 or 9 elements are written
 * starting from the indicated field index \a k. The data block must
 * be wide enough to accomodate this.
 * @param fp file pointer for output.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_mesh_write_msh(BEM3DMesh *m, BEM3DMeshData *f, gint k,
			  gchar *view, gdouble t,
			  bem3d_gmsh_mode_t mode,
			  gint offset,
			  FILE *fp)

{
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(fp != NULL, BEM3D_NULL_ARGUMENT) ;

  return bem3d_mesh_write_msh_4_1(m, f, k, view, t, mode, offset, fp) ;

  g_assert_not_reached() ;
  
  return BEM3D_SUCCESS ;
}


/**
 * @}
 * 
 */
