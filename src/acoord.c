/* acoord.c
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

/**
 * @defgroup acoord Area/Local coordinates
 * @{
 * 
 */

#include <math.h>
#include <stdlib.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

/** 
 * Calculate the area coordinates of a point on a triangle; taken from 
 * http://huizen.dto.tudelft.nl/deBruijn/programs/suna02.htm
 * 
 * @param t1 first vertex of triangle
 * @param t2 second vertex of triangle
 * @param t3 third vertex of triangle
 * @param x  point
 * @param xi GtsVector containing area coordintes
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_area_coordinates_tri(GtsPoint *t1, GtsPoint *t2, GtsPoint *t3,
				GtsPoint *x, GtsVector xi)

{
  GtsVector L1, L2, n, t ;
  gdouble A, len ;

  gts_vector_init(L1, t2, t1) ;
  gts_vector_init(L2, t3, t2) ;
  gts_vector_cross(n, L1, L2) ;
  len = gts_vector_norm(n) ;
  A = gts_vector_scalar(n,n)/len*0.5 ;
  n[0] /= len ; n[1] /= len ; n[2] /= len ;
  
  gts_vector_init(L1, t2, x) ;
  gts_vector_init(L2, t3, x) ;
  gts_vector_cross(t, L1, L2) ;
  xi[0] = gts_vector_scalar(t,n)*0.5/A ;

  gts_vector_init(L1, t3, x) ;
  gts_vector_init(L2, t1, x) ;
  gts_vector_cross(t, L1, L2) ;
  xi[1] = gts_vector_scalar(t,n)*0.5/A ;

  gts_vector_init(L1, t3, x) ;
  gts_vector_init(L2, t1, x) ;
  gts_vector_cross(t, L1, L2) ;
  xi[2] = gts_vector_scalar(t,n)*0.5/A ;

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
