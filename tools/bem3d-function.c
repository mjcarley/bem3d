/* bem3d-function.c
 * 
 * Copyright (C) 2010, 2017, 2018, 2019 Michael Carley
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

#define SURFACE_DATA_WIDTH     16
#define SURFACE_DATA_INDEX      0
#define SURFACE_DATA_FUNCTION   1
#define SURFACE_DATA_QUADRATURE 2
#define SURFACE_DATA_DATA       3
#define SURFACE_DATA_MESHES     4
#define SURFACE_DATA_OUTPUT     5
#define SURFACE_DATA_FIELDS     6
#define SURFACE_DATA_

#define __USAGE_MESSAGE_1__						\
  "%s: evaluate functions on BEM3D meshes and data using\n"		\
  "analytically-defined functions\n\n"				\
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
  "A set of reduction operations are also available via the `-r' option:\n\n"

#define __USAGE_MESSAGE_3__						\
  "\nand can be applied like this, for example:\n\n"			\
  "   %s -E 4 -i (geometry file) -F source.fn -r max > bc.dat\n\n"	\
  "which performs the same operation as in the previous example and then\n" \
  "outputs the maximum value in each column of the data block and its node\n" \
  "index to stderr.\n\n" \
  "The -w option treats the function as an integrand and outputs a block of\n" \
  "data which can be used as a set of weights for estimating an integral on\n" \
  "a surface. For example, if a single element function is supplied, with\n" \
  "the -w option:\n\n"							\
  "  %s -w -F function.fn -i surface.bem -o weights.dat\n\n"		\
  "an integral of f(x)g(x) over the surface can be estimated using\n\n" \
  "  %s -F multiply.fn -d weights.dat -e data.dat -o output.dat\n\n"	\
  "where data.dat contains g(x) and the function multiply.fn is given by\n" \
  "f[0] = f[0]*g[0].\n\n"						\
  "Adding the -S option generates integral weighting output in the form of\n" \
  "a surface matrix which can be used to implement non-local boundary\n" \
  "conditions.\n\n"

static void detailed_usage(gchar *progname)

{
  fprintf(stderr, __USAGE_MESSAGE_1__, progname, progname, progname) ;

  fprintf(stderr, __USAGE_MESSAGE_2__) ;

  bem3d_reduction_func_list(stderr, "   %s: %s;\n", TRUE) ;

  fprintf(stderr, __USAGE_MESSAGE_3__, progname, progname, progname) ;
  
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

static void write_sparse(gint i, gdouble *d, gint nf, gpointer *sdata)

{
  FILE *fs = sdata[0] ;
  GArray *fields = sdata[1] ;
  gint row = *((gint *)sdata[2]) ;
  gint j, idx ;
  gboolean write_entry = FALSE ;

  for ( j = 0 ; j < fields->len ; j ++ ) {
    idx = g_array_index(fields, gint, j) ;
    if ( d[idx] != 0.0 ) write_entry = TRUE ;
  }
  if ( !write_entry ) return ;

  fprintf(fs, "%d %d", row, i) ;
  for ( j = 0 ; j < fields->len ; j ++ ) {
    idx = g_array_index(fields, gint, j) ;
    fprintf(fs, " %1.16e", d[idx]) ;
  }
  fprintf(fs, "\n") ;
  
  return ;
}

static gint write_surface_matrix(gint i, GtsVertex *v, gpointer *data)

{
  gint imesh = *((gint *)data[SURFACE_DATA_INDEX]) ;
  BEM3DFunction *efunc = data[SURFACE_DATA_FUNCTION] ;
  BEM3DQuadratureRule *q = data[SURFACE_DATA_QUADRATURE] ;
  BEM3DMeshData *f = data[SURFACE_DATA_DATA] ;
  GPtrArray *meshes = data[SURFACE_DATA_MESHES] ;
  FILE *fs = data[SURFACE_DATA_OUTPUT] ;
  GArray *fields = data[SURFACE_DATA_FIELDS] ;
  GtsVector n ;
  gint j, ilimits[64], ni, nd, idx ;
  gdouble limits[64] ;
  gboolean write_line ;
  gpointer sdata[] = {fs, fields, &i} ;
    
  bem3d_mesh_data_clear(f) ;
  bem3d_node_normal(g_ptr_array_index(meshes, imesh), i, n, BEM3D_AVERAGE_MWA) ;
  
  for ( j = 0 ; j < meshes->len ; j ++ ) {
    bem3d_function_integral_weights(efunc,
				    g_ptr_array_index(meshes, j), j,
				    v, n, i, q, f) ;
  }

  /*write data here*/
  write_line = FALSE ;
  bem3d_reduction_func_limits(NULL, f, ilimits, &ni, limits, &nd, NULL) ;
  for ( j = 0 ; j < fields->len ; j ++ ) {
    idx = g_array_index(fields, gint , j) ;
    if ( limits[2*idx+0] != 0.0 || limits[2*idx+1] != 0.0 )
      write_line = TRUE ;
  }

  if ( !write_line ) return 0 ; /*no non-zero entries for this vertex*/

  /* fprintf(fs, "%d ", i) ; */

  bem3d_mesh_data_foreach(f, (BEM3DMeshDataEntryFunc)write_sparse, sdata) ;

  /* fprintf(fs, "\n") ; */
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *m ;
  GPtrArray *meshes ;
  BEM3DMeshData *f, *g, *fmerge ;
  GtsFile *fid ;
  gchar *progname ;
  gchar *opfile, *datafile, *edatafile, *funcfile ;
  gchar **reductions ;
  GLogLevelFlags loglevel ;
  FILE *fs, *fo ;
  gint mesh_data_width, i, j, k, idata[64], nc, ni, nd, nnodes, lineno ;
  gdouble ddata[64] ;
  gchar ch, line[1024], *point_file ;
  GPtrArray *overrides ;
  GArray *isurf ;
  gboolean write_func, merge, write_data, integral_weights ;
  BEM3DFunction *efunc ;
  BEM3DQuadratureRule *q ;
  GtsPoint *x ;
  gchar *header = NULL ;
  gpointer surfdata[SURFACE_DATA_WIDTH] ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  f = g = NULL ;
  integral_weights = FALSE ;
  write_func = FALSE ;
  write_data = TRUE ;
  point_file = NULL ;
  merge = FALSE ;
  overrides = g_ptr_array_new() ;
  loglevel = G_LOG_LEVEL_MESSAGE ;
  mesh_data_width = 1 ;
  opfile = datafile = edatafile = funcfile = NULL ; fs = NULL ;
  reductions = NULL ;
  nnodes = 0 ;
  isurf = g_array_new(TRUE, TRUE, sizeof(gint)) ;
  
  meshes = g_ptr_array_new() ;
  m = NULL ;
  
  bem3d_shapefunc_lookup_init() ;

  while ( (ch = getopt(argc, argv, "HhDd:E:e:F:l:mi:o:pqr:S:s:wX:")) != EOF ) {
    switch (ch) {
    default: 
    case 'h':
      fprintf(stderr, 
	      "%s: evaluate functions on BEM3D meshes and data using\n"
	      "analytically-defined functions\n\n",
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
	      "        -m merge the input (-d) and extra (-e) data files and "
	      "output\n"
	      "           as a single block\n"
	      "        -i <bem3d mesh input file> (can be repeated for "
	      "multiple meshes)\n"
	      "        -o <output file name>\n"
	      "        -p (write expanded parsed function to stderr)\n"
	      "        -q do not write evaluated function data to output\n"
	      "        -r <comma separated list of reduction operations>\n"
	      "        -s <expression> (set a function variable)\n"
	      "        -S (list) write output as surface matrix using "
	      "components\n"
	      "           specified in (list)\n"
	      "        -w treat function as integrand of integral weights "
	      "and output\n"
	      "           as surface condition matrix\n"
	      "        -X <file name> list of points in bem3d-dump "
	      "format\n"
	      ) ;
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
    case 'i': append_mesh_from_file(meshes,  optarg) ; break ;
    case 'o': opfile = g_strdup(optarg) ; break ;
    case 'p': write_func = TRUE ; break ;
    case 'q': write_data = FALSE ; break ;
    case 'r': reductions = g_strsplit_set(optarg, " ,", 0) ; break ;
    case 'S': range_int(optarg, isurf) ; break ;
    case 's': g_ptr_array_add(overrides, g_strdup(optarg)) ; break ;
    case 'w': integral_weights = TRUE ; break ;      
    case 'X': point_file = g_strdup(optarg) ; break ;
    }
  }

  bem3d_logging_init(stderr, "", loglevel, NULL) ;
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

    bem3d_function_expand_functions(efunc) ;

    if ( write_func) bem3d_function_write(efunc, stderr) ;
  }

  if ( merge ) {
    if ( (datafile == NULL) || (edatafile == NULL) ) {
      fprintf(stderr,
	      "Data file (-f) and extra data file (-e) must be specified "
	      "for data merge\n") ;
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

  if ( point_file != NULL ) {
    lineno = 1 ;
    x = gts_point_new(gts_point_class(), 0, 0, 0) ;
    fs = file_open(point_file, "-", "r", stdin) ;

    if ( opfile == NULL ) fo = stdout ;
    else
      fo = file_open(opfile, "-", "w", stdout) ;
    
    while ( (nc = fscanf(fs, "%[^\n]c", line)) != EOF && ( nc != 0 ) ) {
      nc = sscanf(line, "%d %lg %lg %lg", &j, &(x->x), &(x->y), &(x->z)) ;
      if ( nc != 4 ) {
	fprintf(stderr, "%s: cannot parse line %d\n  %s\n",
		progname, lineno, line) ;
	exit(1) ;
      }
      fprintf(fo, "%d", j) ;
      j = bem3d_function_eval_point(efunc, x, NULL, 0, ddata, 8) ;
      for ( i = 0 ; i <= j ; i ++ ) fprintf(fo, " %lg", ddata[i]) ;
      fprintf(fo, "\n") ;
      lineno ++ ;
      if ( (nc = fscanf(fs, "%*c")) == EOF ) break ;
    }

    file_close(fo) ;
    
    return 0 ;
  }

  if ( f == NULL ) {
    if ( meshes->len == 0 ) {
      fprintf(stderr,
	      "%s: if no data file name is specified (-d), "
	      "a mesh must be specified\n", progname) ;
      exit(1) ;
    }
    m = g_ptr_array_index(meshes, 0) ;
    nnodes = bem3d_mesh_node_number(m) ;
    f = bem3d_mesh_data_new(m, mesh_data_width) ;
    for ( i = 1 ; i < meshes->len ; i ++ ) {
      m = g_ptr_array_index(meshes, i) ;
      nnodes += bem3d_mesh_node_number(m) ;
      bem3d_mesh_data_add_mesh(f, m) ;
    }
    bem3d_mesh_data_clear(f) ;
  }

  if ( integral_weights ) {
    if ( funcfile == NULL ) {
      fprintf(stderr,
	      "%s: need a function definition to generate integral weights\n",
	      progname) ;
      return 1 ;
    }
    if ( meshes->len == 0 ) {
      fprintf(stderr,
	      "%s: need meshes to generate integral weights\n",
	      progname) ;
      return 1 ;
    }
      
    if ( opfile == NULL ) fs = stdout ;
    else fs = file_open(opfile, "-", "w", stdout) ;

    q = bem3d_quadrature_rule_new(7, 1) ;

    surfdata[SURFACE_DATA_INDEX] = &i ;
    surfdata[SURFACE_DATA_FUNCTION] = efunc ;
    surfdata[SURFACE_DATA_QUADRATURE] = q ;
    surfdata[SURFACE_DATA_DATA] = f ;
    surfdata[SURFACE_DATA_MESHES] = meshes ;
    surfdata[SURFACE_DATA_OUTPUT] = fs ;
    surfdata[SURFACE_DATA_FIELDS] = isurf ;

    if ( isurf->len != 0 ) {
      fprintf(fs, "%d %d BEM3DSurface\n",
	      bem3d_mesh_data_node_number(f), isurf->len) ;
      fprintf(fs, "sparse\n") ;
    
      for ( i = 0 ; i < meshes->len ; i ++ ) {
	bem3d_mesh_foreach_node(g_ptr_array_index(meshes, i),
				(BEM3DNodeFunc)write_surface_matrix,
				surfdata) ;
      }
    } else {
      for ( i = 0 ; i < meshes->len ; i ++ ) {
  	bem3d_function_integral_weights(efunc,
  					g_ptr_array_index(meshes, i), i,
  					NULL, NULL, 0, q, f) ;
      }
      bem3d_mesh_data_write(f, fs, header) ;      
    }
    
    file_close(fs) ;

    return 0 ;
  }
  
  if ( funcfile != NULL && integral_weights == FALSE ) {
    /* basic evaluation of surface functions*/
    bem3d_function_apply_mesh_list(efunc, meshes, f, g) ;
  }
  
  if ( write_data ) {
    if ( opfile == NULL ) fs = stdout ;
    else fs = file_open(opfile, "-", "w", stdout) ;
    bem3d_mesh_data_write(f, fs, header) ;
    
    file_close(fs) ;
  }

  if ( reductions != NULL ) {
    mesh_data_width = bem3d_mesh_data_element_number(f) ;
    if ( meshes->len == 0 ) m = NULL ;
    else m = g_ptr_array_index(meshes, 0) ;
    
    for ( i = 0 ; reductions[i] != NULL ; i ++ ) {
      fprintf(stdout, "%s:", reductions[i]) ;
      bem3d_reduction_func_apply(reductions[i],
				 m, f, idata, &ni, ddata, &nd, NULL) ;
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
