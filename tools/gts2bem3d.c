/* gts2bem3d.c
 * 
 * Copyright (C) 2006, 2009, 2018 Michael Carley
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

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

gint main(gint argc, gchar **argv)

{
  gchar ch ;
  BEM3DMesh *m ;
  BEM3DElementBuildFunc bfunc = NULL ;
  GtsSurface *s ;
  GtsFile *fid ;
  gint order, nne = 0, n ;
  gdouble sharp_edge_angle ;
  extern char *optarg;
  gchar *progname ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  bem3d_logging_init(stderr, "", G_LOG_LEVEL_MESSAGE, NULL) ;

  bem3d_shapefunc_lookup_init() ;
  order = 1 ; n = 0 ; 
  sharp_edge_angle = 15.0 ;
  while ( (ch = getopt(argc, argv, "a:e:hn:")) != EOF ) {
    switch (ch) {
    case 'a': 
      sscanf(optarg, "%lg", &sharp_edge_angle) ;
      break ;
    case 'e': sscanf(optarg, "%d", &order) ; break ;
    default:
    case 'h':
      fprintf(stderr, 
	      "%s: convert GTS to BEM3D format\n\n", progname) ;
      fprintf(stderr, 
	      "Usage %s <options> < input.gts > output.bem3d\n", 
	      progname) ;
      fprintf(stderr, "Options:\n") ;
      fprintf(stderr, "   -h print this message and exit\n") ;
      fprintf(stderr, "   -a # (angle for determining sharp edges)\n") ;
      fprintf(stderr, "   -e # (order of element)\n") ;
      fprintf(stderr, "   -n # (index of first node)\n") ;
      return 0 ;
      break ;
    case 'n': sscanf(optarg, "%d", &n) ; break ;
    }
  }

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  sharp_edge_angle *= M_PI/180.0 ;

  switch ( order ) {
  default: g_assert_not_reached() ; break ;
  case 0:  bfunc = bem3d_element_build_t0 ; nne = 1 ; break ;
  case 1:  bfunc = bem3d_element_build_t1 ; nne = 1 ; break ;
  case 2:  bfunc = bem3d_element_build_t2 ; nne = 2 ; break ;
  case 3:  bfunc = bem3d_element_build_t3 ; nne = 3 ; break ;
  }

  s = gts_surface_new(gts_surface_class(),
		      gts_face_class(),
		      gts_edge_class(),
		      gts_vertex_class()) ;
  fid = gts_file_new(stdin)  ;
  gts_surface_read(s, fid) ;

  m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		   gts_edge_class(), gts_vertex_class()) ;
  bem3d_mesh_discretize(s, nne, bfunc, m) ;

  n = bem3d_mesh_index_nodes(m, sharp_edge_angle, n) ;

  fprintf(stderr, "%s: last index: %d\n", progname, n) ;

  bem3d_mesh_write(m, stdout) ;

  return 0 ;
}
