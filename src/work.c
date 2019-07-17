/* work.c
 * 
 * Copyright (C) 2019 Michael Carley
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

#include <gqr.h>

#include "bem3d.h"
#include "bem3d-private.h"

/**
 * @defgroup work BEM3D workspaces
 *
 * Workspaces used for memory management and thread safety in various
 * calculations. Currently under development.
 *
 * @{
 * 
 */

/** 
 * Allocate a new ::BEM3DWorkspace. This should be called to allocate
 * space before any function using the workspace is called. 
 * 
 * @return a newly allocated ::BEM3DWorkspace
 */

BEM3DWorkspace *bem3d_workspace_new(void)

{
  BEM3DWorkspace *w ;
  gint i ;
  
  w = (BEM3DWorkspace *)g_malloc0(sizeof(BEM3DWorkspace)) ;

  w->q = bem3d_quadrature_rule_new(0, 1) ;

  for ( i = 0 ; i < BEM3D_WORKSPACE_DOUBLE_ARRAY_NUMBER ; i ++ ) {
    w->doubles[i] = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
    w->used_doubles[i] = FALSE ;
  }

  for ( i = 0 ; i < BEM3D_WORKSPACE_GTS_POINT_NUMBER ; i ++ ) {
    w->points[i] = gts_point_new(gts_point_class(), 0, 0, 0) ;
    w->used_points[i] = FALSE ;
  }

  for ( i = 0 ; i < BEM3D_WORKSPACE_GQR_RULE_NUMBER ; i ++ ) {
    w->gqr[i] = gqr_rule_alloc(8) ;
    w->used_gqr[i] = FALSE ;
  }

  return w ;
}

GArray *bem3d_workspace_double_array_get(BEM3DWorkspace *w)

{
  gint i ;

  for ( i = 0 ; i < BEM3D_WORKSPACE_DOUBLE_ARRAY_NUMBER ; i ++ ) {
    if ( w->used_doubles[i] == FALSE ) {
      w->used_doubles[i] = TRUE ;
  
      return w->doubles[i] ;
    }
  }

  g_error("%s: workspace ran out of double arrays", __FUNCTION__) ;
  
  return NULL ;
}
  
gint bem3d_workspace_double_array_put(BEM3DWorkspace *w, GArray *g)

{
  gint i ;

  for ( i = 0 ; i < BEM3D_WORKSPACE_DOUBLE_ARRAY_NUMBER ; i ++ ) {
    if ( w->doubles[i] == g ) {
      w->used_doubles[i] = FALSE ;
  
      return 0 ;
    }
  }

  g_error("%s: could not return double array to workspace", __FUNCTION__) ;
  
  return -1 ;
}

GtsPoint *bem3d_workspace_gts_point_get(BEM3DWorkspace *w)

{
  gint i ;

  for ( i = 0 ; i < BEM3D_WORKSPACE_GTS_POINT_NUMBER ; i ++ ) {
    if ( w->used_points[i] == FALSE ) {
      w->used_points[i] = TRUE ;
  
      return w->points[i] ;
    }
  }

  g_error("%s: workspace ran out of Gtspoints", __FUNCTION__) ;

  return NULL ;
}

gint bem3d_workspace_gts_point_put(BEM3DWorkspace *w, GtsPoint *p)

{
  gint i ;

  for ( i = 0 ; i < BEM3D_WORKSPACE_GTS_POINT_NUMBER ; i ++ ) {
    if ( w->points[i] == p ) {
      w->used_points[i] = FALSE ;
  
      return 0 ;
    }
  }  
  
  g_error("%s: could not return GtsPoint to workspace", __FUNCTION__) ;

  return -1 ;
}


gqr_rule_t *bem3d_workspace_gqr_rule_get(BEM3DWorkspace *w)

{
  gint i ;

  for ( i = 0 ; i < BEM3D_WORKSPACE_GQR_RULE_NUMBER ; i ++ ) {
    if ( w->used_gqr[i] == FALSE ) {
      w->used_gqr[i] = TRUE ;
  
      return w->gqr[i] ;
    }
  }

  g_error("%s: workspace ran out of GQR rules", __FUNCTION__) ;

  return NULL ;
}

gint bem3d_workspace_gqr_rule_put(BEM3DWorkspace *w, gqr_rule_t *g)

{
  gint i ;

  for ( i = 0 ; i < BEM3D_WORKSPACE_GQR_RULE_NUMBER ; i ++ ) {
    if ( w->gqr[i] == g ) {
      w->used_gqr[i] = FALSE ;
  
      return 0 ;
    }
  }  
  
  g_error("%s: could not return GQR rule to workspace", __FUNCTION__) ;

  return -1 ;
}


/**
 * @}
 * 
 */
