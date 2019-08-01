/* bem3d2pos.c
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

/** 
@page bem3d2pos bem3d2pos

bem3d2pos converts BEM3D meshes to GMSH .pos format. It is now
superseded by @ref bem3d2msh

The GMSH pos format is deprecated, since visualization is now done
with mesh files, but it is still sufficient for visualization of data
from BEM3D calculations.

The most basic invocation of @c bem3d2pos is
@verbatim
bem3d2pos -i input.bem -o output.pos @endverbatim
which generates a visualization file @c output.pos which can be viewed
with
@verbatim
gmsh output.pos @endverbatim

An extra option allows you to visualize the surface with solution
information added:
@verbatim
bem3d2pos -i input.bem -d input.dat -f 4 -o output.pos @endverbatim
will generate a file @c output.pos which contains the geometry of 
@c input.bem with the nodes coloured according to the value of field
4 of the data block in @c input.dat (the field option defaults to 0). 

The edge data for the input file can also be visualized, for example,
to select the trailing edge in a @ref bem3daero "bem3d-aero" 
calculation:
@verbatim
bem3d2pos -i input.bem -e input-0.edg -e input-1.edg -o output.pos @endverbatim
generates a .pos file of @c input.bem with additional entities for the
edge files @c input-0.edg and @c input-1.edg.

**/

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
  gchar *progname ;
  gchar *ipfile = NULL, *opfile = NULL, *datafile = NULL ;
  gchar *title = NULL ;
  GLogLevelFlags loglevel ;
  FILE *fs, *ip ;
  BEM3DEdge *edge ;
  bem3d_gmsh_mode_t mode ;
  gint field ;
  gchar ch ;
  gint i ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  mode = BEM3D_GMSH_SCALAR ; field = 0 ; f = NULL ;
  loglevel = G_LOG_LEVEL_MESSAGE ;
  edgefiles = g_ptr_array_new() ;
  while ( (ch = getopt(argc, argv, "d:e:f:ghl:i:o:t:")) != EOF ) {
    switch (ch) {
    default: 
    case 'h':
      fprintf(stderr, 
	      "%s: translate a BEM3D geometry and data to gmsh .pos format\n\n",
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
	      "        -i <bem3d input file>\n"
	      "        -o <output file name>\n"
	      "        -t <title for GMSH view>\n",
	      field) ;
      fprintf(stderr,
	      "%s is now superseded by bem3d2msh and is no longer "
	      "supported.\n", progname) ;
      return 0 ;
      break ;
    case 'd': datafile = g_strdup(optarg) ; break ;
    case 'e': g_ptr_array_add(edgefiles, g_strdup(optarg));
    case 'f': field = atoi(optarg) ; break ;
    case 'g': mode = BEM3D_GMSH_VECTOR ; break ;
    case 'l': loglevel = 1 << atoi(optarg) ; break ;
    case 'i': ipfile = g_strdup(optarg) ; break ;
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

  if ( f == NULL ) {
    f = bem3d_mesh_data_new(m, field+1) ;
    bem3d_mesh_data_clear(f) ;
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
  
  bem3d_mesh_write_pos(m, f, field, title, mode, fs) ;

  if ( edgefiles->len != 0 ) {
    edge = bem3d_edge_new(bem3d_edge_class()) ;
    for ( i = 0 ; i < edgefiles->len ; i ++ ) {
      ip = file_open(g_ptr_array_index(edgefiles,i), "", "r", NULL) ;
      if ( ip == NULL ) {
	fprintf(stderr, "%s: cannot open edge file %s.\n",
		progname, opfile) ;
      }
      bem3d_edge_clear(edge) ;
      fid = gts_file_new(ip) ;
      bem3d_edge_read(edge, fid) ;
      file_close(ip) ; gts_file_destroy(fid) ;
      bem3d_edge_link_to_mesh(edge, m) ;
      bem3d_edge_write_pos(edge, g_ptr_array_index(edgefiles,i), fs) ;
    }
  }

  return 0 ;
}
