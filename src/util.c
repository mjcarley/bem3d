/* bem3d
 * 
 * Copyright (C) 2006, 2009, 2017, 2018 Michael Carley
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

#include "bem3d-private.h"

gint range_double(gchar *arg, GArray *range) ;
gint range_int(gchar *arg, GArray *range) ;

static gchar *file_mode_string(gchar *mode)

{
  gchar *mode_strings[] = {"reading", "writing", "appending",
			   "unknown operation"} ;

  switch ( mode[0] ) {
  default: return mode_strings[3] ; break ;
  case 'r': return mode_strings[0] ; break ;
  case 'w': return mode_strings[1] ; break ;
  case 'a': return mode_strings[2] ; break ;
  }

  return NULL ;
}

FILE *file_open(gchar *fname, gchar *namedefault, gchar *mode, 
		FILE *fdefault)

{
  FILE *f ;

  if ( !strcmp(fname, namedefault) ) return fdefault ;

  f = fopen(fname, mode) ;

  if ( f == NULL )
    g_error("%s: cannot open file \"%s\" for %s (mode \"%s\")", 
	    __FUNCTION__, fname, file_mode_string(mode), mode) ;

  return f ;
}

gint file_close(FILE *f)

{
  if ( (f != stdin) && (f != stdout) && (f != stderr) ) fclose(f) ;

  return 0 ;
}

static gboolean is_imaginary_unit(gchar c)

{
  return ( (c == 'i') || (c == 'I') || (c == 'j') || (c == 'J') ) ;
}

gint parse_complex(gchar *v, gdouble z[])

{
  gint i ;
  gchar *nptr, *eptr ;
  gdouble val ;
  gboolean imag ;

  g_strchug(g_strchomp(v)) ;
  
  for ( i = 0 ; i < strlen(v) ; i ++ ) {
    if ( v[i] == ' ' ) g_strchug(&(v[i])) ;
  }

  z[0] = z[1] = 0.0 ;

  nptr = v ; imag = FALSE ;
  if ( is_imaginary_unit(nptr[0]) ) {
    val = strtod(&(nptr[1]), &eptr) ;
    if ( eptr == &(nptr[1]) ) return -1 ;
    z[1] = val ;
    if ( *eptr == '\0' ) return 0 ;
    imag = TRUE ;
  } else {  
    val = strtod(nptr, &eptr) ;
    if ( nptr == eptr ) return -1 ;
  }

  if ( is_imaginary_unit(eptr[0]) ) {
    z[1] = val ; nptr = &(eptr[1]) ;
  } else {
    if ( !imag ) z[0] = val ; 
    nptr = eptr ;
  }

  if ( *nptr == '\0' ) return 0 ;

  imag = FALSE ;
  if ( is_imaginary_unit(nptr[1]) ) {
    val = strtod(&(nptr[2]), &eptr) ;
    if ( eptr == &(nptr[2]) ) return -1 ;
    z[1] = val ;
    if ( nptr[0] == '-' ) z[1] = -z[1] ;
    if ( *eptr == '\0' ) return 0 ;
    imag = TRUE ;
  } else {  
    val = strtod(nptr, &eptr) ;
    if ( nptr == eptr ) return -1 ;
  }

  if ( is_imaginary_unit(eptr[0]) ) z[1] = val ;
  else { if ( !imag ) z[0] = val ; }

  return 0 ;
}

/* Functions which parse a Matlab/Octave style colon-delimited range
   specification and return a filled GArray of the appropriate
   values */

/* gint range_double(gchar *arg, GArray *range) ; */
/* gint range_int(gchar *arg, GArray *range) ; */

static gint range_colon_double(gchar *arg, GArray *range)

{
  gchar **tokens ;
  guint nt ;
  gint i, nst ;
  gdouble min, max, step, x ;

  tokens = g_strsplit(arg, ":",0) ;

  for ( nt = 0 ; tokens[nt] != NULL ; nt ++ ) ;

  min = max = step = 0.0 ;
  switch ( nt ) {
  case 1: min = max = g_strtod(tokens[0], NULL) ; step = 1.0 ; break ;
  case 2: 
    min = g_strtod(tokens[0], NULL) ;
    max = g_strtod(tokens[1], NULL) ;
    step = 1.0 ;
    break ;
  case 3:
    min = g_strtod(tokens[0], NULL) ;
    step = g_strtod(tokens[1], NULL) ;
    max = g_strtod(tokens[2], NULL) ;
    break ;
  default: g_assert_not_reached() ; break ;
  }

  /*needed to avoid rounding problems*/
  nst = (max-min+1e-9)/step ;

  g_assert(min <= max) ;

  for ( i = 0 ; i <= nst ; i ++ ) {
    x = min + i*step ;
    g_array_append_val(range, x) ;
  }

  return 0 ;
}

gint range_double(gchar *arg, GArray *range)

{
  gchar **tokens ;
  guint i ;

  tokens = g_strsplit(arg, ",",0) ;

  g_array_set_size(range, 0) ;
  for ( i = 0 ; tokens[i] != NULL ; i ++ )
    range_colon_double(tokens[i], range) ;

  return 0 ;
}

static gint range_colon_int(gchar *arg, GArray *range)

{
  gchar **tokens ;
  guint nt ;
  gint min, max, step, x ;

  tokens = g_strsplit(arg, ":",0) ;

  for ( nt = 0 ; tokens[nt] != NULL ; nt ++ ) ;

  min = max = step = 0.0 ;
  switch ( nt ) {
  case 1: min = max = strtol(tokens[0], NULL, 10) ; step = 1 ; break ;
  case 2: 
    min = strtol(tokens[0], NULL, 10) ;
    max = strtol(tokens[1], NULL, 10) ;
    step = 1 ;
    break ;
  case 3:
    min = strtol(tokens[0], NULL, 10) ; 
    step = strtol(tokens[1], NULL, 10) ;
    max = strtol(tokens[2], NULL, 10) ;
    break ;
  default: g_assert_not_reached() ; break ;
  }

  g_assert(min <= max) ;

  for ( x = min ; x <= max ; x += step )
    g_array_append_val(range, x) ;

  return 0 ;
}

gint range_int(gchar *arg, GArray *range)

{
  gchar **tokens ;
  guint i ;

  tokens = g_strsplit(arg, ",",0) ;

  g_array_set_size(range, 0) ;
  for ( i = 0 ; tokens[i] != NULL ; i ++ )
    range_colon_int(tokens[i], range) ;

  return 0 ;
}

gint printf_fixed_width(gchar *str, gint width, gchar *prefix, FILE *f)

{
  gint i, stride ;
  gchar *line, buf[1024] ;
  
  stride = width - strlen(prefix) ;

  line = str ;
  do { 
    if ( strlen(line) <= stride ) {
      fprintf(f, "%s%s\n", prefix, line) ;
      return 0 ;
    }

    /*there must be a better way of doing this than copying to a buffer
     but playing with pointers in-situ didn't work*/
    for ( i = stride ; line[i] != ' ' ; i -- ) ;
    strncpy(buf, line, i+1) ;
    buf[i] = '\0' ;
    fprintf(f, "%s%s\n", prefix, buf) ;
    line = &(line[i+1]) ;
  } while (1) ;
  
  return 0 ;
}
