/* bem3d-sphere.c
 * 
 * Copyright (C) 2006, 2008 Michael Carley
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

static gint shift_vertex(GtsVertex *v, gpointer *data)

{
  gdouble sx = *(gdouble *)data[0] ;
  gdouble sy = *(gdouble *)data[1] ;
  gdouble sz = *(gdouble *)data[2] ;
  gdouble sc = *(gdouble *)data[3] ;
  gdouble du = *(gdouble *)data[4] ;
  gdouble dv = *(gdouble *)data[5] ;
  gdouble dw = *(gdouble *)data[6] ;
  gdouble phi, th ;

  phi = acos(GTS_POINT(v)->z) ;
  th = atan2(GTS_POINT(v)->y,GTS_POINT(v)->x) ;

  GTS_POINT(v)->x = sc*sx*cos(th)*sin(phi) ;
  GTS_POINT(v)->y = sc*sy*sin(th)*sin(phi) ;
  GTS_POINT(v)->z = sc*sz*cos(phi) ;
  
  GTS_POINT(v)->x += du ; GTS_POINT(v)->y += dv ;
  GTS_POINT(v)->z += dw ;

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *m ;
  GtsSurface *s ;
  BEM3DElementBuildFunc bfunc = NULL ;
  gint order, refine, nne = 0, n ;
  gboolean write_last_node ;
  gchar ch ;
  gchar *opfile, *progname ; 
  gdouble sx, sy, sz, R ;
  gdouble u, v, w ;
  gdouble sharp_edge_angle ;
  gpointer data[8] ;
  FILE *f ;
  extern char *optarg;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  bem3d_logging_init(stderr, "", G_LOG_LEVEL_MESSAGE, NULL) ;
  bem3d_shapefunc_lookup_init() ;
  order = 1 ; refine = 3 ; sx = sy = sz = 1.0 ; R = 1.0 ;
  u = v = w = 0.0 ; n = 0 ; write_last_node = FALSE ;
  opfile = NULL ; sharp_edge_angle = M_PI ;
  while ( (ch = getopt(argc, argv, "e:hln:o:R:r:u:v:w:x:y:z:")) != EOF ) {
    switch (ch) {
    case 'e': sscanf(optarg, "%d", &order) ; break ;
    default: 
    case 'h':
      fprintf(stderr, 
	      "%s: output a boundary element discretization of a "
	      "sphere or ellipsoid\n\n", progname) ;
      fprintf(stderr, "Usage: %s <options> > output.bem3d\n"
	      "Options:\n"
	      "        -e # (order of element)\n"
	      "        -h (print this message and exit)\n"
	      
	      "        -l (write the index of last node plus one to stderr)\n"
	      "        -n # (index of first node)\n"
	      "        -o <string> (name of output file)\n"
	      "        -R # (sphere radius)\n"
	      "        -r # (sphere refinement)\n"
	      "        -u # (sphere centre X)\n"
	      "        -v # (sphere centre Y)\n"
	      "        -w # (sphere centre Z)\n"
	      "        -x # (scaling factor in X)\n"
	      "        -y # (scaling factor in Y)\n"
	      "        -z # (scaling factor in Z)\n", progname) ;
      return 0 ;
      break ;
    case 'l': write_last_node = TRUE ; break ;
    case 'n': sscanf(optarg, "%d", &n) ; break ;
    case 'o': opfile = g_strdup(optarg) ; break ;
    case 'R': R = atof(optarg) ; break ;
    case 'r': refine = atoi(optarg) ; break ;
    case 'u': u = atof(optarg) ; break ;
    case 'v': v = atof(optarg) ; break ;
    case 'w': w = atof(optarg) ; break ;
    case 'x': sx = atof(optarg) ; break ;
    case 'y': sy = atof(optarg) ; break ;
    case 'z': sz = atof(optarg) ; break ;
    }
  }

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

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

  m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		     gts_edge_class(), gts_vertex_class()) ;
  s = gts_surface_generate_sphere(s, refine) ;  
  bem3d_mesh_discretize(s, nne, bfunc, m) ;
  data[0] = &sx ; data[1] = &sy ; data[2] = &sz ; data[3] = &R ;
  data[4] = &u ; data[5] = &v ; data[6] = &w ;
  gts_surface_foreach_vertex(GTS_SURFACE(m), (GtsFunc)shift_vertex, data) ;

  n = bem3d_mesh_index_nodes(m, sharp_edge_angle, 0) ;

  if ( opfile == NULL ) bem3d_mesh_write(m, stdout) ;
  else {
    f = file_open(opfile, "-", "w", stdout) ;
    bem3d_mesh_write(m, f) ;
    file_close(f) ;
  }

  if ( write_last_node ) 
    fprintf(stdout, "%s: last indexed point: %d\n", progname, n) ;

  return 0 ;
}
