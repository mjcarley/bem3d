/* bem3d2msh.c
 * 
 * Copyright (C) 2019 Michael Carley
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
@page bem3d2msh bem3d2msh

bem3d2msh converts BEM3D meshes to GMSH .msh format, with the option
to include surface data.

The most basic invocation of @c bem3d2msh is
@verbatim
bem3d2msh -i input.bem -o output.msh @endverbatim
which generates a mesh file @c output.msh which can be viewed
with
@verbatim
gmsh output.msh @endverbatim

An extra option allows you to visualize the surface with solution
information added:
@verbatim
bem3d2msh -i input.bem -d input.dat -f 4 -o output.msh @endverbatim
will generate a file @c output.msh which contains the geometry of 
@c input.bem with the nodes coloured according to the value of field
4 of the data block in @c input.dat (the field option defaults to 0). 

The main purpose of the program is visualization of geometries and
data, rather than meshing.

*
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <glib.h>

#include <gts.h>
#include <gqr.h>

#include "bem3d.h"
#include "bem3d-private.h"

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *m ;
  BEM3DMeshData *f ;
  GtsFile *fid ;
  GPtrArray *edgefiles ;
  gchar *progname, *ipfile = NULL, *opfile = NULL, *datafile = NULL ;
  gchar *title = NULL ;
  GLogLevelFlags loglevel ;
  FILE *fs ;
  bem3d_gmsh_mode_t mode ;
  gdouble t ;
  gint field, idx0 ;
  gchar ch ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  t = 0.0 ;
  mode = BEM3D_GMSH_SCALAR ; field = 0 ; idx0 = 0 ; f = NULL ;
  loglevel = G_LOG_LEVEL_MESSAGE ;
  edgefiles = g_ptr_array_new() ;
  while ( (ch = getopt(argc, argv, "hd:e:f:gi:l:n:o:t:")) != EOF ) {
    switch (ch) {
    default: 
    case 'h':
      fprintf(stderr, 
	      "%s: convert BEM3D geometry and (optionally) data to "
	      "gmsh .msh format\n\n",
	      progname) ;
      fprintf(stderr, "Usage: %s <options>\n", progname) ;
      fprintf(stderr, 
	      "Options:\n"
	      "        -h print this message and exit\n"
	      "        -d <data file name>\n"
	      "        -e <edge file name (possibly multiple)>\n"
	      "        -f # field of data to plot (%d)\n"
	      "        -g output three components of gradient\n"
	      "        -l # set logging level\n"
	      "        -n # offset added to node indices (%d)\n"
	      "        -i <bem3d input file>\n"
	      "        -o <output file name>\n"
	      "        -t <title for GMSH view>\n",
	      field, idx0) ;
      return 0 ;
      break ;
    case 'd': datafile = g_strdup(optarg) ; break ;
    case 'e': g_ptr_array_add(edgefiles, g_strdup(optarg));
    case 'f': field = atoi(optarg) ; break ;
    case 'g': mode = BEM3D_GMSH_VECTOR ; break ;
    case 'i': ipfile = g_strdup(optarg) ; break ;
    case 'l': loglevel = 1 << atoi(optarg) ; break ;
    case 'n': idx0 = atoi(optarg) ; break ;
    case 'o': opfile = g_strdup(optarg) ; break ;
    case 't': title = g_strdup(optarg) ; break ;
    }
  }

  bem3d_logging_init(stderr, "", loglevel, NULL) ;
  gqr_logging_init(stderr, "", loglevel, NULL) ;
  bem3d_shapefunc_lookup_init() ;

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  if ( ipfile == NULL ) fs = stdin ;
  else {
    fs = file_open(ipfile, "-", "r", stdin) ;
    if ( fs == NULL ) {
      fprintf(stderr, "%s: cannot open input file %s.\n",
	      progname, ipfile) ;
      return 1 ;
    }
  }
    
  fid = gts_file_new(fs) ;
  m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		   gts_edge_class(), gts_vertex_class()) ;
  bem3d_mesh_read(m, fid) ;

  file_close(fs) ;

  if ( datafile != NULL ) {
    fs = file_open(datafile, "-", "r", stdin) ;
    if ( fs == NULL ) {
      fprintf(stderr, "%s: cannot open data file %s.\n",
	      progname, datafile) ;
      return 1 ;
    }
    bem3d_mesh_data_read(&f, fs, 0) ;
    file_close(fs) ;
  }

  if ( bem3d_mesh_data_node_number(f) < bem3d_mesh_node_number(m) ) {
    fprintf(stderr, "%s: number of mesh nodes (%d) is greater than number of "
	    "data points (%d)\n", progname, bem3d_mesh_node_number(m),
	    bem3d_mesh_data_node_number(f)) ;
    return 1 ;
  }

  if ( opfile == NULL ) fs = stdout ;
  else {
    fs = file_open(opfile, "-", "w", stdout) ;
    if ( fs == NULL ) {
      fprintf(stderr, "%s: cannot open output file %s.\n",
	      progname, opfile) ;
    }
  }
  
  bem3d_mesh_write_msh(m, f, field, title, t, mode, idx0, fs) ;

  /* if ( edgefiles->len != 0 ) { */
  /*   edge = bem3d_edge_new(bem3d_edge_class()) ; */
  /*   for ( i = 0 ; i < edgefiles->len ; i ++ ) { */
  /*     ip = file_open(g_ptr_array_index(edgefiles,i), "", "r", NULL) ; */
  /*     if ( ip == NULL ) { */
  /* 	fprintf(stderr, "%s: cannot open edge file %s.\n", */
  /* 		progname, opfile) ; */
  /*     } */
  /*     bem3d_edge_clear(edge) ; */
  /*     fid = gts_file_new(ip) ; */
  /*     bem3d_edge_read(edge, fid) ; */
  /*     file_close(ip) ; gts_file_destroy(fid) ; */
  /*     bem3d_edge_link_to_mesh(edge, m) ; */
  /*     bem3d_edge_write_(edge, g_ptr_array_index(edgefiles,i), fs) ; */
  /*   } */
  /* } */

  return 0 ;
}
