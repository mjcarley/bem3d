/* bem3d-function.c
 * 
 * Copyright (C) 2010, 2017, 2018 Michael Carley
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
#include <unistd.h>
#include <string.h>

#include <glib.h>

#include <gts.h>

#ifdef HAVE_LIBMATHEVAL
#include <matheval.h>
#endif /*HAVE_LIBMATHEVAL*/

#include "bem3d.h"
#include "bem3d-private.h"

#define __USAGE_MESSAGE_1__						\
  "%s: evaluate functions on BEM3D meshes and data using\n"		\
  "analytically-specified functions\n\n"				\
  "Typical uses are generation of boundary conditions and post-processing\n" \
  "of results using functions. Functions are input as a BEM3DFunction\n" \
  "specified using the -F option (see below). A simple example would be" \
  "\n\n"								\
  "   BEM3DFunction\n"							\
  "   x0 = 0.2\n"							\
  "   y0 = 0.1\n"							\
  "   z0 = -0.3\n"							\
  "   k = 1.0\n"							\
  "   R = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)\n"			\
  "   dR = (nx*(x-x0) + ny*(y-y0) + nz*(z-z0))/R^3\n"			\
  "   C = cos(k*R)\n"							\
  "   S = sin(k*R)\n"							\
  "   f[0] = C/R/4/pi\n"						\
  "   f[1] = S/R/4/pi\n"						\
  "   f[2] = -(C + k*R*S)/4/pi*dR\n"					\
  "   f[3] =  (k*R*C - S)/4/pi*dR\n\n"					\
  "which generates the real and imaginary parts of the point source\n"	\
  "potential (cos(kR) + i sin(kR))/(4 pi R) and its normal derivative\n" \
  "for a point source of wavenumber k=1.0 located at (0.2,0.1,-0.3) and\n" \
  "inserts them in a data block suitably formatted to be passed to\n"	\
  "bem3d-solve. If the function is in a file called \"source.fn\", it can\n" \
  "be applied via\n\n"							\
  "   %s -E 4 -i (geometry file) -F source.fn > bc.dat\n\n"		\
  "The `-s' option can be used to set new values for the variables in the\n" \
  "function. For example,\n\n"						\
  "   %s -E 4 -i (geometry file) -F source.fn -s \"k=0.3\" > bc.dat\n\n" \
  "generates the potential with the wavenumber set to k=0.3. Quotes are\n" \
  "recommended to make sure the full expression is passed, especially\n" \
  "when it includes blank space.\n\n"

#define __USAGE_MESSAGE_2__						\
  "A set of reduction operations are also available via the `-R' option:\n\n"

#define __USAGE_MESSAGE_3__						\
  "\nand can be applied like this, for example:\n\n"			\
  "   %s -E 4 -i (geometry file) -F source.fn -R max > bc.dat\n\n"	\
  "which performs the same operation as in the previous example and then\n" \
  "outputs the maximum value in each column of the data block to stderr.\n"


static void detailed_usage(gchar *progname)

{
  fprintf(stderr, __USAGE_MESSAGE_1__, progname, progname, progname) ;

  fprintf(stderr, __USAGE_MESSAGE_2__) ;

  bem3d_reduction_func_list(stderr, "   %s: %s;\n", TRUE) ;

  fprintf(stderr, __USAGE_MESSAGE_3__, progname) ;
  
  return ;
}

#ifndef HAVE_LIBMATHEVAL
gint main(gint argc, gchar **argv)

{
  fprintf(stderr, 
	  "%s: requires libmatheval. Recompile with libmatheval enabled\n",
	  argv[0]) ; 

  return 0 ;
}

#else /*HAVE_LIBMATHEVAL*/

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *m, *m0 ;
  BEM3DMeshData *f, *g, *fmerge ;
  GtsFile *fid ;
  gchar *progname ;
  gchar *ipfile, *opfile, *datafile, *edatafile, *funcfile ;
  gchar **reductions ;
  GLogLevelFlags loglevel ;
  FILE *fs ;
  gint mesh_data_width, i, j, k ;
  gint idata[64], ni, nd ;
  gdouble ddata[64] ;
  gchar ch ;
  GPtrArray *overrides ;
  gboolean write_func, merge ;
  BEM3DFunction *efunc ;
  BEM3DMotion *motion ;
  gchar *header = NULL ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  f = g = NULL ;
  write_func = FALSE ;
  merge = FALSE ;
  overrides = g_ptr_array_new() ;
  loglevel = G_LOG_LEVEL_MESSAGE ;
  mesh_data_width = 1 ;
  ipfile = opfile = datafile = edatafile = funcfile = NULL ; fs = NULL ;
  reductions = NULL ;
  
  while ( (ch = getopt(argc, argv, "HhDd:E:e:F:l:mi:o:pr:s:")) != EOF ) {
    switch (ch) {
    default: 
    case 'h':
      fprintf(stderr, 
	      "%s: evaluate functions on BEM3D meshes and data using\n"
	      "analytically-specified functions\n\n",
	      progname) ;
      fprintf(stderr, "Usage: %s <options>\n", progname) ;
      fprintf(stderr, 
	      "Options:\n"
	      "        -h (print this message and exit)\n"
	      "        -D (write data as surface data file)\n"
	      "        -d <data file name>\n"
	      "        -e <extra data file name>\n"
	      "        -E # (expand output mesh data block to this size)\n"
	      "        -F <function definition file>\n"
	      "        -H (print longer help message and exit)\n"
	      "        -l # (set logging level)\n"
	      "        -m merge the input (-d) and extra (-e) data files\n"
	      "           and output as a single block\n"
	      "        -i <bem3d mesh input file>\n"
	      "        -o <output file name>\n"
	      "        -p (write expanded parsed function to stderr)\n"
	      "        -r <comma separated list of reduction operations>\n"
	      "        -s <expression> (set a function variable)\n") ;
      return 0 ;
      break ;
    case 'H':
      detailed_usage(progname) ;
      return 0 ;
      break ;
    case 'd': datafile = g_strdup(optarg) ; break ;
    case 'D':
      header = "%d %d BEM3DSurface\ndiagonal\n" ;
      break ;
    case 'e': edatafile = g_strdup(optarg) ; break ;
    case 'E': mesh_data_width = atoi(optarg) ; break ;
    case 'F': funcfile = g_strdup(optarg) ; break ;
    case 'l': loglevel = 1 << atoi(optarg) ; break ;
    case 'm': merge = TRUE ; break ;
    case 'i': ipfile = g_strdup(optarg) ; break ;
    case 'o': opfile = g_strdup(optarg) ; break ;
    case 'p': write_func = TRUE ; break ;
    case 'r': reductions = g_strsplit_set(optarg, " ,", 0) ; break ;
    case 's': g_ptr_array_add(overrides, g_strdup(optarg)) ; break ;
    }
  }

  bem3d_logging_init(stderr, "", loglevel, NULL) ;
  bem3d_shapefunc_lookup_init() ;

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  efunc = bem3d_function_new(bem3d_function_class()) ;
  
  if ( funcfile != NULL ) {
    fs = file_open(funcfile, "-", "r", stdin) ;

    fid = gts_file_new(fs) ;
    bem3d_function_read(efunc, fid) ;
    file_close(fs) ;

    for ( i = 0 ; i < overrides->len ; i ++ ) {
      bem3d_function_insert_string(efunc,
				   (gchar *)g_ptr_array_index(overrides,i)) ;
    }

    if ( write_func) bem3d_function_write(efunc, stderr) ;
  }

  if ( merge ) {
    if ( (datafile == NULL) || (edatafile == NULL) ) {
      fprintf(stderr,
	      "Data file and extra data file must be specified for data "
	      "merge\n") ;
      exit(1) ;
    }

    fs = file_open(datafile, "-", "r", stdin) ;
    if ( fs == NULL ) {
      fprintf(stderr, "%s: cannot open data file %s.\n",
	      progname, datafile) ;
      return 1 ;
    }
    bem3d_mesh_data_read(&f, fs, mesh_data_width) ;
    file_close(fs) ;
    
    fs = file_open(edatafile, "-", "r", stdin) ;
    if ( fs == NULL ) {
      fprintf(stderr, "%s: cannot open data file %s.\n",
	      progname, datafile) ;
      return 1 ;
    }
    bem3d_mesh_data_read(&g, fs, 0) ;
    file_close(fs) ;

    fmerge = bem3d_mesh_data_merge(f, g, FALSE) ;

    if ( opfile == NULL ) fs = stdout ;
    else
      fs = file_open(opfile, "-", "w", stdout) ;

    bem3d_mesh_data_write(fmerge, fs, header) ;
    
    return 0 ;
  }
  
  if ( ipfile == NULL ) fs = stdin ;
  else fs = file_open(ipfile, "-", "r", stdin) ;
    
  fid = gts_file_new(fs) ;
  m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		   gts_edge_class(), gts_vertex_class()) ;
  bem3d_mesh_read(m, fid) ;

  file_close(fs) ;

  if ( ipfile == NULL ) fs = stdin ;
  else fs = file_open(ipfile, "-", "r", stdin) ;
    
  fid = gts_file_new(fs) ;
  m0 = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		      gts_edge_class(), gts_vertex_class()) ;
  bem3d_mesh_read(m0, fid) ;

  file_close(fs) ;

  if ( datafile != NULL ) {
    fs = file_open(datafile, "-", "r", stdin) ;
    if ( fs == NULL ) {
      fprintf(stderr, "%s: cannot open data file %s.\n",
	      progname, datafile) ;
      return 1 ;
    }
    bem3d_mesh_data_read(&f, fs, mesh_data_width) ;
    file_close(fs) ;
  }

  if ( f == NULL ) {
    f = bem3d_mesh_data_new(m, mesh_data_width) ;
    bem3d_mesh_data_clear(f) ;
  }

  if ( edatafile != NULL ) {
    fs = file_open(edatafile, "-", "r", stdin) ;
    if ( fs == NULL ) {
      fprintf(stderr, "%s: cannot open data file %s.\n",
	      progname, datafile) ;
      return 1 ;
    }
    bem3d_mesh_data_read(&g, fs, 0) ;
    file_close(fs) ;
  }

  motion = bem3d_motion_new(bem3d_motion_class(), m, m0) ;

  bem3d_motion_expand_defs(motion) ;
  bem3d_motion_create_evaluators(motion) ;
  bem3d_motion_mesh_position(motion, 0.0) ;

  if ( funcfile != NULL ) {
    bem3d_function_expand_functions(efunc) ;
    bem3d_function_apply(efunc, motion, 0.0, f, g) ;
  }
  
  if ( opfile == NULL ) fs = stdout ;
  else
    fs = file_open(opfile, "-", "w", stdout) ;
  bem3d_mesh_data_write(f, fs, header) ;

  file_close(fs) ;

  if ( reductions != NULL ) {
    mesh_data_width = bem3d_mesh_data_element_number(f) ;
    for ( i = 0 ; reductions[i] != NULL ; i ++ ) {
      bem3d_reduction_func_apply(reductions[i],
				 m, f, idata, &ni, ddata, &nd, NULL) ;
      fprintf(stdout, "%s:", reductions[i]) ;
      for ( j = 0 ; j < mesh_data_width ; j ++ ) {
	for ( k = 0 ; k < ni ; k ++ )
	  fprintf(stdout, " %d", idata[j*ni+k]) ;
	for ( k = 0 ; k < nd ; k ++ )
	  fprintf(stdout, " %lg", ddata[j*nd+k]) ;
      }
      fprintf(stdout, "\n") ;
    }
  }
  
  return 0 ;
}

#endif /*HAVE_LIBMATHEVAL*/
