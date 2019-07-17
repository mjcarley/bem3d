/* data.c
 * 
 * Copyright (C) 2006, 2010, 2018 Michael Carley
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
 * @defgroup data Data on BEM3D surfaces
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

#include <wmpi.h>

#include "bem3d.h"
#include "bem3d-private.h"

#define BEM3D_EFDATA_WIDTH  8
#define BEM3D_EFDATA_DATA   0
#define BEM3D_EFDATA_FUNC   1
#define BEM3D_EFDATA_MDATA  2

#define DATA_FOREACH_DATA_WIDTH 8
#define DATA_FOREACH_FUNC       0
#define DATA_FOREACH_BLOCK      1
#define DATA_FOREACH_HASH       2
#define DATA_FOREACH_DATA       3

gdouble *_multiproc_buffer = NULL ;
gint _multiproc_nb ;

static gint _bem3d_init_data(gint i, GtsVertex *v, BEM3DMeshData *m)

{
  bem3d_mesh_data_add_node(m, i) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Allocate a ::BEM3DMeshData block, sized to the number of collocation
 * points on a ::BEM3DMesh.
 * 
 * @param m ::BEM3DMesh for which to allocate data
 * @param n number of data points per mesh node
 * 
 * @return pointer to new data block
 */

BEM3DMeshData *bem3d_mesh_data_new(BEM3DMesh *m, gint n)

{
  BEM3DMeshData *d = NULL ;
  gint i, imin, imax ;

  g_return_val_if_fail (m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MESH(m), NULL) ;
  g_return_val_if_fail (n > 0, NULL) ;
  
  imin = bem3d_mesh_node_index_min(m) ;
  imax = bem3d_mesh_node_index_max(m) ;  

  bem3d_mesh_node_index_min(m) = 0 ;
  bem3d_mesh_node_index_max(m) = G_MAXINT ;

  i = bem3d_mesh_node_number(m) ;
  d = bem3d_mesh_data_sized_new(n, i) ;

  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)_bem3d_init_data, d) ;

  bem3d_mesh_node_index_min(m) = imin ;
  bem3d_mesh_node_index_max(m) = imax ;
  
  return d ;  
}

/** 
 * Free a ::BEM3DMeshData.
 * 
 * @param d ::BEM3DMeshData to be freed.
 * 
 * @return BEM3D_SUCCESS on success ;
 */

gint bem3d_mesh_data_free(BEM3DMeshData *d)

{
  g_return_val_if_fail(d != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( d->t != NULL ) g_hash_table_destroy(d->t) ;
  if ( d->d != NULL ) g_array_free(d->d, TRUE) ;

  g_free(d) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Allocate ::BEM3DMeshData block of a given size, without linking to a
 * ::BEM3DMesh.
 * 
 * @param n number of data points per mesh node;
 * @param m number of points for which to allocate the ::BEM3DMeshData.
 * 
 * @return ::BEM3DMeshData block.
 */

BEM3DMeshData *bem3d_mesh_data_sized_new(gint n, gint m)

{
  BEM3DMeshData *d = NULL ;

  g_return_val_if_fail (n > 0, NULL) ;
  g_return_val_if_fail (m > 0, NULL) ;

  d = (BEM3DMeshData *)g_malloc(sizeof(BEM3DMeshData)) ;
  d->nd = n ;

  d->d = g_array_sized_new(TRUE, TRUE, sizeof(gdouble), n*m) ;
  d->t = g_hash_table_new(NULL, NULL) ;

  return d ;  
}

/** 
 * Set all entries of a ::BEM3DMeshData block to zero.
 * 
 * @param m ::BEM3DMeshData to clear.
 * 
 * @return 0 on success
 */

gint bem3d_mesh_data_clear(BEM3DMeshData *m)

{
  gint i ;
  
  g_return_val_if_fail (m != NULL, BEM3D_NULL_ARGUMENT) ;

  for ( i = 0 ; i < m->d->len ; i ++ ) g_array_index(m->d,gdouble,i) = 0.0 ;

  return BEM3D_SUCCESS ;
}

/** 
 * Look up the data for a given mesh node. 
 * 
 * @param m ::BEM3DMeshData block;
 * @param i global index of node to look up.
 * 
 * @return a pointer to the data for node i, NULL if the index is not
 * in the block.
 */

gdouble *bem3d_mesh_data_get(BEM3DMeshData *m, gint i)

{
  gdouble *x = NULL ;
  gint j ;

  g_return_val_if_fail (m != NULL, NULL) ;

  if ( (j = GPOINTER_TO_INT(g_hash_table_lookup(m->t, 
						GINT_TO_POINTER(i+1)))-1)
       == -1 ) return NULL ;
  x = &g_array_index(m->d,gdouble,(j*(m->nd))) ;

  return x ;
}


/** 
 * Add an array of data to a mesh data entry.
 * 
 * @param f ::BEM3DMeshData block;
 * @param i global node index;
 * @param g array of entries to add to field i.
 * 
 * @return 0 on success
 */

gint bem3d_mesh_data_add(BEM3DMeshData *f, gint i, GArray *g)

{
  gint j ;
  gdouble *x ;

  g_return_val_if_fail (f != NULL, BEM3D_NULL_ARGUMENT) ;
  
  g_return_val_if_fail(bem3d_mesh_data_element_number(f) >= g->len,
		       BEM3D_ARGUMENT_OUT_OF_RANGE) ;

  x = bem3d_mesh_data_get(f, i) ;
  if ( x == NULL ) return BEM3D_ARGUMENT_OUT_OF_RANGE ;
  for ( j = 0 ; j < g->len ; j ++ ) {
    x[j] += g_array_index(g,gdouble,j) ;
  }
  
  return BEM3D_SUCCESS ;
}

static void _data_write(gpointer key, gpointer value, gpointer data[])

{
  BEM3DMeshData *d = data[0] ;
  FILE *fp = data[1] ;
  gint i, j, k ;

  i = GPOINTER_TO_INT(key)-1 ;
  k = (GPOINTER_TO_INT(value)-1)*bem3d_mesh_data_element_number(d) ;

  fprintf(fp, "%d", i) ;

  for ( j = 0 ; j < bem3d_mesh_data_element_number(d) ; j ++ ) 
    fprintf(fp, " %1.16e", g_array_index(d->d,gdouble,k+j)) ;
  fprintf(fp, "\n") ;
  
  return ;
}

/** 
 * Write a ::BEM3DMeshData block to a file. The output format is one
 * line containing:
 *
 * @e NP @e NE BEM3DMeshData
 *
 * where @e NP is the number of nodes in the block and @e NE is the
 * number of numerical entries per point, followed by @e NP lines each
 * containing @e NE numerical entries. If \a header is not NULL, it is
 * used in place of the standard first line. If it contains C
 * formatting strings %d, these generate @e NP and @e NE. 
 * 
 * @param f ::BEM3DMeshData block;
 * @param fp file pointer;
 * @param header first line of output file (may be NULL)
 * 
 * @return 0 on success.
 */

gint bem3d_mesh_data_write(BEM3DMeshData *f, FILE *fp, gchar *header)

{
  gpointer data[2] ;

  g_return_val_if_fail (f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail (fp != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( header == NULL ) 
    fprintf(fp, "%d %d BEM3DMeshData\n", 
	    g_hash_table_size(f->t), bem3d_mesh_data_element_number(f)) ;
  else
    fprintf(fp, header, 
	    g_hash_table_size(f->t), bem3d_mesh_data_element_number(f)) ;
    
  data[0] = f ; data[1] = fp ;
  g_hash_table_foreach(f->t, (GHFunc)_data_write, data) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Read a ::BEM3DMeshData block from an input file, allocating it as
 * necessary.
 * 
 * @param f ::BEM3DMEshData block to allocate;
 * @param fp input file stream;
 * @param width data is set to the maximum of \a width and the number of
 * elements specified in the input file; set to 0 to use the number specified
 * in the file.
 * 
 * @return 0 on success.
 */

gint bem3d_mesh_data_read(BEM3DMeshData **f, FILE *fp, gint width)

{
  gint w, ne, np, i, j, k ;
  gchar line[1024] ;

  g_return_val_if_fail (f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail (fp != NULL, BEM3D_NULL_ARGUMENT) ;

  /*np: number of points in data block*/
  /*ne: number of data elements per point*/
  fscanf(fp, "%d %d", &np, &ne) ;
  fscanf(fp, "%[^\n]s", line) ;
  fscanf(fp, "%*c") ;

  w = MAX(ne, width) ;
  *f = bem3d_mesh_data_sized_new(w, np) ;
  bem3d_mesh_data_element_number(*f) = w ;

  g_array_set_size((*f)->d, np*w) ;
  for ( i = 0 ; i < np ; i ++ ) {
    fscanf(fp, "%d", &j) ;
    g_hash_table_insert((*f)->t, 
			GINT_TO_POINTER(j+1),
			GINT_TO_POINTER(i+1)) ;
    for ( k = 0 ; k < ne ; k ++ ) 
      fscanf(fp, "%lg", &(g_array_index((*f)->d,gdouble,w*i+k))) ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Add a node to a ::BEM3DMeshData. 
 * 
 * @param m mesh data block to add information to;
 * @param i index of node to add. A check is performed to ensure that 
 * it is not already present. 
 * 
 * @return 0 on success.
 */

gint bem3d_mesh_data_add_node(BEM3DMeshData *m, gint i)

{
  gdouble *x ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( (x = bem3d_mesh_data_get(m, i)) != NULL) {
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, 
	  "%s: node %d already in BEM3DMeshData", 
	  __FUNCTION__, i) ;
    return BEM3D_FAILURE ;
  }

  g_hash_table_insert(m->t, 
		      GINT_TO_POINTER(i+1),
		      GINT_TO_POINTER((m->d->len)/m->nd)+1) ;
  g_array_set_size(m->d,m->d->len+m->nd) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Sum the data stored in a ::BEM3DMeshData block across all processes
 * so that each process holds the same data.
 * 
 * @param m ::BEM3DMeshData block to sum
 * 
 * @return 0 on success.
 */

gint bem3d_mesh_data_multiproc_sum(BEM3DMeshData *m)

{
  /* static gdouble *buffer = NULL ; */
  /* static gint nb = 0 ; */

  if ( wmpi_process_number() == 1 ) return 0 ;

  if ( _multiproc_buffer == NULL ) {
    _multiproc_nb = m->d->len ;
    _multiproc_buffer = (gdouble *)g_malloc(_multiproc_nb*sizeof(gdouble)) ;
  }

  if ( _multiproc_nb < m->d->len) {
    _multiproc_nb = m->d->len ;
    _multiproc_buffer = (gdouble *)g_realloc(_multiproc_buffer,
					     _multiproc_nb*sizeof(gdouble)) ;
  }

  g_memmove(_multiproc_buffer, m->d->data, m->d->len*sizeof(gdouble)) ;

  wmpi_sum_all_double((gdouble *)(m->d->data),
		      _multiproc_buffer, m->d->len) ;

  wmpi_pause() ;

  return 0 ;
}

/** 
 * Find the number of nodes in a ::BEM3DMeshData.
 * 
 * @param d a ::BEM3DMeshData.
 * 
 * @return the number of nodes in \a d, or zero if d is NULL.
 */

gint bem3d_mesh_data_node_number(BEM3DMeshData *d)

{
  g_return_val_if_fail(d != NULL, 0) ;

  return g_hash_table_size(d->t) ;
}

gint bem3d_mesh_data_expand(BEM3DMeshData *d, gint ne)

{
  GArray *f ;
  gint i, j, n ;

  g_return_val_if_fail(d != NULL, 0) ;

  if ( ne < bem3d_mesh_data_element_number(d) ) return BEM3D_SUCCESS ;

  n = bem3d_mesh_data_node_number(d) ;
  
  f = g_array_sized_new(TRUE, TRUE, sizeof(gdouble), n*ne) ;
  g_array_set_size(f, n*ne) ;

  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j < bem3d_mesh_data_element_number(d) ; j ++ ) {
      g_array_index(f, gdouble, i*ne+j) = 
	g_array_index(d->d, gdouble, i*bem3d_mesh_data_element_number(d)+j) ;
    }
  }

  g_array_free(d->d, TRUE) ;
  d->d = f ; bem3d_mesh_data_element_number(d) = ne ;

  return BEM3D_SUCCESS ;
}

static void _data_merge(gpointer key, gpointer value, gpointer data[])

{
  BEM3DMeshData *f1 = data[0] ;
  BEM3DMeshData *fm = data[1] ;
  gdouble *d1, *dm ;
  gint i, j ;
  
  i = GPOINTER_TO_INT(key)-1 ;

  if ( bem3d_mesh_data_get(fm, i) != NULL ) return ;
  bem3d_mesh_data_add_node(fm, i) ;
  d1 = bem3d_mesh_data_get(f1, i) ;
  dm = bem3d_mesh_data_get(fm, i) ;
  
  for ( j = 0 ; j < bem3d_mesh_data_element_number(f1) ; j ++ ) dm[j] = d1[j] ;
  
  return ;
}

/** 
 * Merge two data blocks into a single block. This is intended for
 * merging data generated for two disjoint meshes, with different
 * indices. If an index appears in both blocks the behaviour is
 * undefined. The output block will be sized to accomodate all data
 * from both inputs, with zero padding where required.
 * 
 * @param f1 an input ::BEM3DMeshData;
 * @param f2 another input ::BEM3DMeshData;
 * @param strict (not currently used).
 * 
 * @return a new ::BEM3DMeshData containing the data from \a f1 and \a
 * f2.
 */

BEM3DMeshData *bem3d_mesh_data_merge(BEM3DMeshData *f1,
				     BEM3DMeshData *f2,
				     gboolean strict)

{
  BEM3DMeshData *d = NULL ;
  gint nn, ne ;
  gpointer data[2] ;
    
  /*number of elements in new data block*/
  ne = MAX(bem3d_mesh_data_element_number(f1),
	   bem3d_mesh_data_element_number(f2)) ;
  /*number of nodes*/
  nn = bem3d_mesh_data_node_number(f1) + bem3d_mesh_data_node_number(f2) ;

  d = bem3d_mesh_data_sized_new(ne, nn) ;

  data[1] = d ;
  if ( !strict ) {
    data[0] = f1 ;
    g_hash_table_foreach(f1->t, (GHFunc)_data_merge, data) ;
    data[0] = f2 ;
    g_hash_table_foreach(f2->t, (GHFunc)_data_merge, data) ;

    return d ;
  }

  /*yet to decide how strict "strict" is*/
  g_assert_not_reached() ; 
  
  return d ;  
}

/** 
 * Expand a ::BEM3DMeshData to include entries for another mesh, so
 * that the data block contains data for more than one mesh.
 * 
 * @param d an existing ::BEM3DMeshData
 * @param m a ::BEM3DMesh whose nodes are to be added to \a d
 * 
 * @return 0 on success
 */

gint bem3d_mesh_data_add_mesh(BEM3DMeshData *d, BEM3DMesh *m)

{
  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)_bem3d_init_data, d) ;  
  
  return 0 ;
}

static void _data_write_weights(gpointer key, gpointer value, gpointer data[])

{
  BEM3DMeshData *d = data[0] ;
  FILE *fp = data[1] ;
  gint *fields = data[2] ;
  gint nfields = *((gint *)data[3]) ;
  gint i, j, k ;
  gboolean zeros ;
  
  i = GPOINTER_TO_INT(key)-1 ;
  k = (GPOINTER_TO_INT(value)-1)*bem3d_mesh_data_element_number(d) ;

  zeros = TRUE ;
  for ( j = 0 ; j < nfields ; j ++ ) {
    g_assert(fields[j] < bem3d_mesh_data_element_number(d)) ;
    if ( g_array_index(d->d,gdouble,k+fields[j]) != 0.0 )
      zeros = FALSE ;
  }
  if ( zeros ) return ;
  
  fprintf(fp, "%d", i) ;

  for ( j = 0 ; j < nfields ; j ++ ) 
    fprintf(fp, " %1.16e", g_array_index(d->d,gdouble,k+fields[j])) ;
  
  return ;
}

/** 
 * Write a block of data in the form of a single row suitable for use
 * in a weighting matrix, e.g. of a type generated by
 * ::bem3d_function_integral_weights. Data are written on one row in
 * the form
 *
 * index f0 f1 f2 index f0 f1 f2 ...
 * 
 * with any entries with all data equal to zero neglected, effectively
 * treating the data as a sparse matrix. 
 *
 * @param f data to write;
 * @param fp file pointer for output;
 * @param fields indices of fields to write to output (must all be less
 * than the number of fields per node in \a f); 
 * @param nfields number of fields to write.
 * 
 * @return 0 on success
 *
 */

gint bem3d_mesh_data_write_weights(BEM3DMeshData *f, FILE *fp,
				   gint *fields, gint nfields)

{
  gpointer data[4] ;

  g_return_val_if_fail (f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail (fp != NULL, BEM3D_NULL_ARGUMENT) ;
    
  data[0] = f ;
  data[1] = fp ;
  data[2] = fields ;
  data[3] = &nfields ;
  g_hash_table_foreach(f->t, (GHFunc)_data_write_weights, data) ;

  return BEM3D_SUCCESS ;
}

static void _entry_func(gpointer key, gpointer value, gpointer fd[])

{
  BEM3DMeshDataEntryFunc func = fd[DATA_FOREACH_FUNC] ;
  BEM3DMeshData *f = fd[DATA_FOREACH_BLOCK] ;
  gpointer data = fd[DATA_FOREACH_DATA] ;
  gdouble *d ;
  gint i ;

  i = GPOINTER_TO_INT(key)-1 ;
  d = bem3d_mesh_data_get(f, i) ;
  func(i, d, bem3d_mesh_data_element_number(f), data) ;
  
  return ;
}

gint bem3d_mesh_data_foreach(BEM3DMeshData *f, BEM3DMeshDataEntryFunc func,
			     gpointer data)

{
  gpointer fd[DATA_FOREACH_DATA_WIDTH] ;

  fd[DATA_FOREACH_FUNC] = func ;
  fd[DATA_FOREACH_BLOCK] = f ;
  fd[DATA_FOREACH_DATA] = data ;

  g_hash_table_foreach(f->t, (GHFunc)_entry_func, fd) ;
  
  return 0 ;
}

/**
 * @}
 * 
 */

