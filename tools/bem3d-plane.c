/* bem3d-plane.c
 * 
 * Copyright (C) 2006, 2018 Michael Carley
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
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

static gint skew_vertices(GtsVertex *v, gpointer data[])

{
  gdouble xi, eta, L[4] ;
  gint i ;

  xi = GTS_POINT(v)->x ;
  eta = GTS_POINT(v)->y ;
  L[0] = (1.0-xi)*(1.0-eta) ; L[1] = xi*(1.0-eta) ;
  L[2] = xi*eta ; L[3] = (1.0-xi)*eta ;

  gts_point_set(GTS_POINT(v), 0, 0, 0) ;
  for ( i = 0 ; i < 4 ; i ++ ) {
    GTS_POINT(v)->x += L[i]*GTS_POINT(data[i])->x ;
    GTS_POINT(v)->y += L[i]*GTS_POINT(data[i])->y ;
    GTS_POINT(v)->z += L[i]*GTS_POINT(data[i])->z ;
  }  

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  GtsSurface *s ;
  gdouble x, y, z ;
  gchar ch ;
  gpointer data[8] ;
  gint i, ni, nj ;
  gchar *progname ;
  extern char *optarg;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  ni = nj = 16 ;

  data[0] = gts_point_new(gts_point_class(), 0, 0, 0) ;
  data[1] = gts_point_new(gts_point_class(), 1, 0, 0) ;
  data[2] = gts_point_new(gts_point_class(), 1, 1, 0) ;
  data[3] = gts_point_new(gts_point_class(), 0, 1, 0) ;
  
  i = 0 ;
  while ( (ch = getopt(argc, argv, "hi:j:v:")) != EOF ) {
    switch (ch) {
    default:
    case 'h':
      fprintf(stderr, 
	      "%s: generate a GTS triangulation of a quadrilateral patch\n\n", 
	      progname) ;
      fprintf(stderr, "Usage %s <options> > output.gts\n", progname) ;
      fprintf(stderr, "Options:\n") ;
      fprintf(stderr, "        -h (print this message and exit)\n") ;
      fprintf(stderr, "        -i # (number of nodes in u direction)\n") ;
      fprintf(stderr, "        -j # (number of nodes in v direction)\n") ;
      fprintf(stderr, 
	      "        -v #,#,# (position of corner vertex, repeat for subsequent corners)\n") ;
      exit(0) ;
      break ;
    case 'i': sscanf(optarg, "%d", &ni) ; break ;
    case 'j': sscanf(optarg, "%d", &nj) ; break ;
    case 'v': 
      sscanf(optarg, "%lg,%lg,%lg", &x, &y, &z) ; 
      gts_point_set(GTS_POINT(data[i]), x, y, z) ;
      i ++ ;
      break ;
    }
  }

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  s = gts_surface_new(gts_surface_class(),
		      gts_face_class(),
		      gts_edge_class(),
		      gts_vertex_class()) ;
  bem3d_geometry_plane(s, ni, nj) ;

  gts_surface_foreach_vertex(s, (GtsFunc)skew_vertices, data) ;

  gts_surface_write(s, stdout) ;

  return 0 ;
}
