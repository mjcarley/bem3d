/* bem3d-skeleton.c
 * 
 * Copyright (C) 2017 Michael Carley
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
  BEM3DMesh *m ;
  BEM3DMeshSkeleton *s ;
  BEM3DQuadratureRule *q ;
  BEM3DAverage anorm ;
  GtsFile *fp ;
  gchar ch, *progname ;
  GLogLevelFlags loglevel ;
  gint order ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  bem3d_shapefunc_lookup_init() ;
  loglevel = G_LOG_LEVEL_MESSAGE ;

  order = 3 ; anorm = BEM3D_AVERAGE_MWE ;

  while ( (ch = getopt(argc, argv, "hl:q:")) != EOF ) {
    switch (ch) {
    default: 
    case 'h':
      fprintf(stderr, 
	      "%s: write the skeleton points and normals of a BEM3D "
	      "geometry file to stdout\n\n",
	      progname) ;
      fprintf(stderr, "Usage: %s [options] < file.bem > file.skel\n", 
	      progname) ;
      fprintf(stderr, "Points are written one per line:\n") ;
      fprintf(stderr, "  <index> <x y z> <nx ny nz>\n\n") ;
      fprintf(stderr, 
	      "Options:\n\n"
	      "  -q # order of quadrature rule for source points (%d)\n",
	      order) ;
      return 0 ;
      break ;
    case 'l': loglevel = 1 << atoi(optarg) ; break ;
    case 'q': order = atoi(optarg) ; break ;
    }
  }

  bem3d_logging_init(stderr, "", loglevel, NULL) ;
  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		   gts_edge_class(), gts_vertex_class()) ;
  fp = gts_file_new(stdin) ;
  bem3d_mesh_read(m, fp) ;

  q = bem3d_quadrature_rule_new(order, 1) ;

  if ( order <= 7 ) 
    bem3d_quadrature_rule_gauss(NULL, bem3d_mesh_element_sample(m), q, 
				NULL, NULL, &order) ;
  else
    bem3d_quadrature_rule_wx(NULL, bem3d_mesh_element_sample(m), q, 
			     NULL, NULL, &order) ;    

  s = bem3d_mesh_skeleton_new(m, order) ;
  bem3d_mesh_skeleton_init(s, q, anorm) ;

  bem3d_mesh_skeleton_write(s, stdout) ;

  /* bem3d_mesh_write_nodes(m, stdout) ; */

  return 0 ;
}
