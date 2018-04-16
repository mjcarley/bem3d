/* logging.c
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
 * @defgroup logging Logging functions and function status codes
 * @{
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>

#include <glib.h>

#include "bem3d.h"
#include "bem3d-private.h"

#define BEM3D_LOGGING_DATA_WIDTH     4
#define BEM3D_LOGGING_DATA_FID       0
#define BEM3D_LOGGING_DATA_PREFIX    1
#define BEM3D_LOGGING_DATA_LEVEL     2
#define BEM3D_LOGGING_DATA_EXIT_FUNC 3

static const gchar *bem3d_logging_string(GLogLevelFlags level)

{
  const gchar *strings[] = {"RECURSION", 
			    "FATAL",
			    "ERROR",
			    "CRITICAL",
			    "WARNING",
			    "MESSAGE",
			    "INFO",
			    "DEBUG"} ;

  if ( G_LOG_LEVEL_ERROR & level) return strings[2] ; 
  if ( G_LOG_LEVEL_CRITICAL & level) return strings[3] ; 
  if ( G_LOG_LEVEL_WARNING & level) return strings[4] ; 
  if ( G_LOG_LEVEL_MESSAGE & level) return strings[5] ; 
  if ( G_LOG_LEVEL_INFO & level) return strings[6] ; 
  if ( G_LOG_LEVEL_DEBUG & level) return strings[7] ; 

  g_assert_not_reached() ;

  return NULL ;
}

static void bem3d_logging_func(const gchar *log_domain,
			       GLogLevelFlags log_level,
			       const gchar *message,
			       gpointer data[])

{
  FILE *f = (FILE *)data[BEM3D_LOGGING_DATA_FID] ;
  gchar *p = (gchar *)data[BEM3D_LOGGING_DATA_PREFIX] ;
  GLogLevelFlags level = *(GLogLevelFlags *)data[BEM3D_LOGGING_DATA_LEVEL] ;
  gint (*exit_func)(void) = data[BEM3D_LOGGING_DATA_EXIT_FUNC] ;

  if ( log_level > level ) return ;

  fprintf(f, "%s%s-%s: %s\n", p, 
	  G_LOG_DOMAIN, bem3d_logging_string(log_level),
	  message) ;
  fflush(f) ;

  if ( log_level <= G_LOG_LEVEL_ERROR ) {
    if ( exit_func != NULL ) exit_func() ;
  }

  return ;
}

/** 
 * Initialize BEM3D logging
 * 
 * @param f file stream for messages, if NULL, stderr is used;
 * @param p string to prepend to messages, or NULL;
 * @param log_level maximum logging level to handle (see g_log);
 * @param exit_func function to call if exiting on an error, or NULL.
 * 
 * @return BEM3D_SUCCESS on success. 
 */

gint bem3d_logging_init(FILE *f, gchar *p, 
			GLogLevelFlags log_level,
			gpointer exit_func)

{
  static gpointer data[BEM3D_LOGGING_DATA_WIDTH] ;
  static GLogLevelFlags level ;

  if ( f != NULL ) 
    data[BEM3D_LOGGING_DATA_FID] = f ;
  else
    data[BEM3D_LOGGING_DATA_FID] = stderr ;    
  if ( p != NULL ) 
    data[BEM3D_LOGGING_DATA_PREFIX] = g_strdup(p) ;
  else
    data[BEM3D_LOGGING_DATA_PREFIX] = g_strdup("") ;

  level = log_level ;
  data[BEM3D_LOGGING_DATA_LEVEL] = &level ;    
    
  g_log_set_handler (G_LOG_DOMAIN, 
		     G_LOG_FLAG_RECURSION |
		     G_LOG_FLAG_FATAL |   
		     G_LOG_LEVEL_ERROR |
		     G_LOG_LEVEL_CRITICAL |
		     G_LOG_LEVEL_WARNING |
		     G_LOG_LEVEL_MESSAGE |
		     G_LOG_LEVEL_INFO |
		     G_LOG_LEVEL_DEBUG,
		     (GLogFunc)bem3d_logging_func, data);

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
