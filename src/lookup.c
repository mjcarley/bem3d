/* lookup.c
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
 * @defgroup lookup Looking up data
 *
 * Various ::BEM3DLookupFunc's for dummy operations, such as calculating
 * local surface coefficients.
 *
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
 * A ::BEM3DLookupFunc which returns 0.0 for the normal derivative and
 * 1.0 for the surface term, to be used in integrating a real Green's
 * function over a surface.
 * 
 * @param i ignored
 * @param j local index
 * @param data ignored
 * @param s surface term
 * @param ds normal derivative
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_lookup_func_unit(gint i, gint j, 
			    gpointer data,
			    GArray *s, GArray *ds)

{
  g_array_index(s, gdouble, j) = 1.0 ;
  g_array_index(ds, gdouble, j) = 0.0 ;

  return BEM3D_SUCCESS ;
}

/** 
 * A BEM3DLookupFunc which returns 0.0+j0.0 for the normal derivative
 * and 1.0+j0.0 for the surface term, to be used in integrating a
 * complex Green's function over a surface.
 * 
 * @param i ignored
 * @param j local index
 * @param data ignored
 * @param s surface term
 * @param ds normal derivative
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_lookup_func_unit_c(gint i, gint j, 
			      gpointer data,
			      GArray *s, GArray *ds)

{
  g_array_index(s, gdouble, 2*j) = 1.0 ;
  g_array_index(s, gdouble, 2*j+1) = 0.0 ;
  g_array_index(ds, gdouble, 2*j) = 0.0 ;
  g_array_index(ds, gdouble, 2*j+1) = 0.0 ;

  return BEM3D_SUCCESS ;
}

/** 
 * A ::BEM3DLookupFunc which returns 1.0 for the normal derivative and
 * 1.0 for the surface term, to be used in testing quadrature rules. 
 * 
 * @param i ignored
 * @param j local index
 * @param data ignored
 * @param s surface term
 * @param ds normal derivative
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_lookup_func_both_unit(gint i, gint j, 
				 gpointer data,
				 GArray *s, GArray *ds)

{
  g_array_index(s, gdouble, j) = 1.0 ;
  g_array_index(ds, gdouble, j) = 1.0 ;

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
