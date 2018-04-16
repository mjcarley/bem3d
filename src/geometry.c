/* geometry.c
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
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <glib.h>
#include <gts.h>

#include "bem3d-private.h"

/**
 * @defgroup geometry Geometric utilities
 * @{
 * 
 */

/** 
 * Generate a GtsSurface of points on a regular grid, adding triangles
 * regularly spaced in x and y.
 * 
 * @param s a GtsSurface;
 * @param ni number of points in x grid direction;
 * @param nj number of points in y grid direction.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_geometry_plane(GtsSurface *s, gint ni, gint nj)

{
  GtsVertex *v1, *v2, *v3, *v ;
  GtsEdge *e1, *e2, *e3 ;
  GtsFace *f ;
  gdouble ee, xi, eta ;
  gint i, j ;

  g_return_val_if_fail(s != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(GTS_IS_SURFACE(s), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(ni > 0, BEM3D_ARGUMENT_OUT_OF_RANGE) ;  
  g_return_val_if_fail(nj > 0, BEM3D_ARGUMENT_OUT_OF_RANGE) ;  

  ee = 1000.0 ;
  v1 = gts_vertex_new(gts_vertex_class(), -ee, -ee, 0) ;
  v2 = gts_vertex_new(gts_vertex_class(), ee, -ee, 0) ;
  v3 = gts_vertex_new(gts_vertex_class(), 0, ee, 0) ;

  e1 = gts_edge_new(gts_edge_class(), v1, v2) ;
  e2 = gts_edge_new(gts_edge_class(), v2, v3) ;
  e3 = gts_edge_new(gts_edge_class(), v3, v1) ;
  f = gts_face_new(gts_face_class(), e1, e2, e3) ;
  gts_surface_add_face(s, f) ;

  for ( i = 0 ; i < ni ; i ++ ) {
    xi = (gdouble)i/(gdouble)(ni-1) ;
    for ( j = 0 ; j < nj ; j ++ ) {
      eta = (gdouble)j/(gdouble)(nj-1) ;
      v = gts_vertex_new(gts_vertex_class(), xi, eta, 0.0) ;
      gts_delaunay_add_vertex(s, v, NULL) ;
    }
  }
  
  gts_allow_floating_vertices = TRUE;
  gts_object_destroy (GTS_OBJECT (v1));
  gts_object_destroy (GTS_OBJECT (v2));
  gts_object_destroy (GTS_OBJECT (v3));
  gts_allow_floating_vertices = FALSE;

  return 0 ;
}

/**
 * @}
 * 
 */
