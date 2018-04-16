/* data.c
 * 
 * Copyright (C) 2006, 2010 Michael Carley
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
 * Extract maximum and minimum values of a function from a BEM3DMeshData
 * block
 * 
 * @param f ::BEM3DMeshData block;
 * @param i field to check;
 * @param xmin minimum value of field;
 * @param xmax maximum value of field.
 * 
 * @return 0 on success
*/

gint bem3d_mesh_function_limits(BEM3DMeshData *f, gint i,
				gdouble *xmin, gdouble *xmax) 
  
{
  gint j ;

  g_return_val_if_fail (f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail (i >= 0 && i < bem3d_mesh_data_element_number(f), 
			BEM3D_ARGUMENT_OUT_OF_RANGE) ;  
  g_return_val_if_fail (xmin != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail (xmax != NULL, BEM3D_NULL_ARGUMENT) ;

  *xmin = G_MAXDOUBLE ; *xmax = -G_MAXDOUBLE ;
  
  for ( j = i ; j < f->d->len ; j += f->nd ) {
    *xmin = MIN(*xmin, g_array_index(f->d,gdouble,j)) ;
    *xmax = MAX(*xmax, g_array_index(f->d,gdouble,j)) ;
  }

  return BEM3D_SUCCESS ;
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
 * number of elements per point, followed by @e NP lines each
 * containing @e NE numerical entries.
 * 
 * @param f ::BEM3DMeshData block;
 * @param fp file pointer.
 * 
 * @return 0 on success.
 */

gint bem3d_mesh_data_write(BEM3DMeshData *f, FILE *fp)

{
  gpointer data[2] ;

  g_return_val_if_fail (f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail (fp != NULL, BEM3D_NULL_ARGUMENT) ;

  fprintf(fp, "%d %d BEM3DMeshData\n", 
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
  static gdouble *buffer = NULL ;
  static gint nb = 0 ;

  if ( wmpi_process_number() == 1 ) return 0 ;

  if ( buffer == NULL ) {
    nb = m->d->len ;
    buffer = (gdouble *)g_malloc(nb*sizeof(gdouble)) ;
  }

  if ( nb < m->d->len) {
    nb = m->d->len ;
    buffer = (gdouble *)g_realloc(buffer, nb*sizeof(gdouble)) ;
  }

  g_memmove(buffer, m->d->data, m->d->len*sizeof(gdouble)) ;

  wmpi_sum_all_double((gdouble *)(m->d->data), buffer, m->d->len) ;

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

/**
 * @}
 * 
 */

