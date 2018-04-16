/* bem3d-boundary.c
 * 
 * Copyright (C) 2006 Michael Carley
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
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

#include "tools.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

static gint boundary_condition_complex(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMesh *m = data[BEM3D_TOOLS_DATA_MESH] ;
  FILE *output = data[BEM3D_TOOLS_DATA_OUTPUT] ;
  BEM3DMeshData *f = data[BEM3D_TOOLS_DATA_DATA] ;
  GtsVector n ;
  gdouble *fi ;

  bem3d_node_normal(m, i, n, BEM3D_AVERAGE_MWAAT) ;
  /* bem3d_node_normal(m, i, n, BEM3D_AVERAGE_MWSELR) ; */
  /* bem3d_node_normal(m, i, n, BEM3D_AVERAGE_MWA) ; */

  g_assert( (fi = bem3d_mesh_data_get(f, i)) != NULL) ; 

  fprintf(output, "%d %1.16e %1.16e %1.16e %1.16e\n", i,
	  fi[0], fi[1],
	  n[0]*fi[2] + n[1]*fi[4] + n[2]*fi[6],
	  n[0]*fi[3] + n[1]*fi[5] + n[2]*fi[7]) ;

  return 0 ;
}

static gint boundary_condition_real(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMesh *m = data[BEM3D_TOOLS_DATA_MESH] ;
  FILE *output = data[BEM3D_TOOLS_DATA_OUTPUT] ;
  BEM3DMeshData *f = data[BEM3D_TOOLS_DATA_DATA] ;
  GtsVector n ;
  gdouble *fi ;

/*   bem3d_node_normal(m, i, n, BEM3D_AVERAGE_MWSELR) ; */
  bem3d_node_normal(m, i, n, BEM3D_AVERAGE_MWA) ;

  g_assert( (fi = bem3d_mesh_data_get(f, i)) != NULL) ; 

  fprintf(output, "%d %1.16e %1.16e\n", i,
	  fi[0], n[0]*fi[1] + n[1]*fi[2] + n[2]*fi[3]) ;

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *m ;
  GtsFile *fp ;
  gchar *ipfile, *opfile, *datfile ;
  gchar ch, p[32], *progname ;
  FILE *input, *output ;
  BEM3DMeshData *sf ;
  gpointer sdata[BEM3D_TOOLS_DATA_WIDTH] ;
  gint nf ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  p[0] = '\0' ;
  bem3d_logging_init(stderr, p, G_LOG_LEVEL_MESSAGE, NULL) ;
  bem3d_shapefunc_lookup_init() ;

  /*default is a Laplace equation (k == 0)*/
  ipfile = opfile = datfile = NULL ;
  m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		     gts_edge_class(), gts_vertex_class()) ;
  sdata[BEM3D_TOOLS_DATA_MESH] = m ;

  while ( (ch = getopt(argc, argv, "hd:i:o:")) != EOF ) {
    switch (ch) {
    default:
    case 'h':
	fprintf(stderr, 
		"%s: set boundary conditions for a mesh and output them\n"
		"as a data block\n\n",
		progname) ;
	fprintf(stderr, "Usage: %s <options>\n", progname) ;
	fprintf(stderr, 
		"Options:\n"
		"        -h (print this message and exit)\n"
		"        -d <surface data file>\n"
		"        -i <bem3d input file>\n"
		"        -o <output file name> (for mesh block data)\n") ;
      return 0 ;
      break ;
    case 'd': datfile = g_strdup(optarg) ; break ;
    case 'i': ipfile = g_strdup(optarg) ; break ;
    case 'o': opfile = g_strdup(optarg) ; break ;
    }
  }

  if ( ipfile == NULL ) ipfile = g_strdup("-") ;
  if ( opfile == NULL ) opfile = g_strdup("-") ;
  if ( datfile == NULL ) datfile = g_strdup("-") ;

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  input = file_open(ipfile, "-", "r", stdin) ;
  fp = gts_file_new(input) ;
  bem3d_mesh_read(m, fp) ;
  file_close(input) ;

  input = file_open(datfile, "-", "r", stdin) ;
  bem3d_mesh_data_read(&sf, input, 0) ;
  file_close(input) ;  

  nf = bem3d_mesh_data_element_number(sf)/2 ;

  output = file_open(opfile, "-", "w", stdout) ;
  sdata[BEM3D_TOOLS_DATA_OUTPUT] = output ;
  sdata[BEM3D_TOOLS_DATA_DATA] = sf ;

  fprintf(output, "%d %d BEM3DMeshData\n", bem3d_mesh_node_number(m), nf) ;

  if ( nf == 4 )
    bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)boundary_condition_complex, 
			    sdata) ;
  else
    bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)boundary_condition_real, sdata) ;

  file_close(output) ;  

  return 0 ;
}
