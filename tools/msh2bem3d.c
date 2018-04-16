/* msh2bem3d.c
 * 
 * Copyright (C) 2006, 2009 Michael Carley
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

/** 
@page msh2bem3d msh2bem3d: convert GMSH mesh files to BEM3D format

@c msh2bem3d is used to convert geometries from the GMSH to the BEM3D
mesh format. This includes sharp edge detection and indexing the
nodes, including multiple indices for nodes on sharp edges. The most
basic invocation is:
@verbatim
msh2bem3d < input.msh > output.bem @endverbatim
which takes a GMSH geometry @c input.msh and converts it to a BEM3D
geometry @c output.bem, with indices starting at 0.

To check for and properly index nodes lying in sharp edges, an angle
must be specified as a criterion for deciding whether or not an edge
is sharp, like this:
@verbatim
msh2bem3d -a 89 -n 573 < input.msh > output.bem @endverbatim
This treats edges with an angle between elements greater than 89 degrees
as sharp and assigns multiple indices to the corresponding nodes as
necessary. As an example, this also shows the use of the option which
sets the first index, here 573. This is useful in problems with multiple 
geometries where each must be indexed starting from the last index of the
previous geometry. 

In scattering problems, it is sufficient merely to properly index
nodes on sharp edges. In aerodynamic problems where a wake is shed,
the trailing edge, from which the wake is generated, must be
identified. In this case: 
@verbatim
msh2bem3d -a 89 -e output-%d.edg < input.msh > output.bem 
@endverbatim will generate a set of @ref edge
"edge files", called @c output-0.edg, @c output-1.edg, etc. for the
connected sharp edges. Use @ref bem3d2pos "bem3d2pos" with the edge output
option to visualize the geometry and its edges and select the edge file(s)
for your problem.

**/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>

#include <glib.h>

#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

gint main(gint argc, gchar **argv)

{
  FILE *input, *output ;
  BEM3DMesh *m ;
  gint n ;
  gdouble sharp_edge_angle ;
  gchar *efile, *ufile ;
  GLogLevelFlags loglevel ;
  GString *fname ;
  GSList *s, *e ;
  gchar ch, *progname ;
  extern char *optarg;
  gboolean orient = TRUE ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  n = 0 ; efile = ufile = NULL ; 
  loglevel = G_LOG_LEVEL_MESSAGE ;
  sharp_edge_angle = 90.0 ;
  while ( (ch = getopt(argc, argv, "a:e:hl:n:u:")) != EOF ) {
    switch (ch) {
    default: 
    case 'h':
      fprintf(stderr, 
	      "%s: convert MSH to BEM3D format\n\n", progname) ;
      fprintf(stderr, 
	      "Usage %s <options> < input.msh > output.bem\n",
	      progname) ;
      fprintf(stderr, 
	      "Options:\n"
	      "        -a # (angle for determining sharp edges)\n"
	      "        -e <edge file name template>\n"
	      "        -h print this message and exit\n"
	      "        -l # (set logging level)\n"
	      "        -n # (index of first node)\n"
	      "        -u <file name for list of unlinked sharp nodes>\n") ;
      return 0 ;
      break ;
    case 'a': sharp_edge_angle = atof(optarg) ; break ;
    case 'e': efile = g_strdup(optarg) ; break ;
    case 'l': loglevel = 1 << atoi(optarg) ; break ;
    case 'n': n = atoi(optarg) ; break ;
    /* case 'o': orient = TRUE ; break ; */
    case 'u': ufile = g_strdup(optarg) ; break ;
    }
  }

  bem3d_logging_init(stderr, NULL, loglevel, NULL) ;
  bem3d_shapefunc_lookup_init() ;

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  sharp_edge_angle *= M_PI/180.0 ;
  input = stdin ;
  m = bem3d_mesh_new(bem3d_mesh_class(),
		   gts_face_class(),
		   gts_edge_class(),
		   gts_vertex_class()) ;

  bem3d_gmsh_read(m, input) ;

  n = bem3d_mesh_index_nodes(m, sharp_edge_angle, n) ;

  fprintf(stderr, "%s: last index: %d\n", progname, n) ;

  bem3d_mesh_write(m, stdout) ;

  s = bem3d_mesh_sharp_vertices(m) ;

  if ( efile != NULL ) {

    e = bem3d_mesh_extract_edges(m, &s) ;
    fname = g_string_new("") ;

    for ( n = 0 ; e != NULL ; (e = e->next), (n ++) ) {
      g_string_printf(fname, efile, n) ;
      output = file_open(fname->str, "-", "w", stdout) ;
      if ( orient && !bem3d_edge_is_oriented(e->data, m) ) {
	fprintf(stderr, "%s: inverting edge %s\n", progname, fname->str) ;
	bem3d_invert_edge(e->data) ;
      }
/*       if ( output == NULL ) { */
/* 	fprintf(stderr, "%s: cannot open file %s", progname, fname->str) ; */
/* 	return -1 ; */
/*       } */
      bem3d_edge_write(BEM3D_EDGE(e->data), output) ;
      file_close(output) ;
    }
  }

  if ( ufile != NULL ) {
    output = file_open(ufile, "-", "w", stdout) ;
/*     if ( output == NULL ) { */
/*       fprintf(stderr, "%s: cannot open file %s", progname, ufile) ; */
/*       return -1 ; */
/*     } */

    for ( ; s != NULL ; s = s->next ) {
      fprintf(output, "%lg %lg %lg\n", 
	      GTS_POINT(s->data)->x, GTS_POINT(s->data)->y,
	      GTS_POINT(s->data)->z) ;
    }
    
    file_close(output) ;
  }

  return 0 ;
}
