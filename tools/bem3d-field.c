/* bem3d-field.c
 * 
 * Copyright (C) 2006, 2008, 2018, 2019, 2020 Michael Carley
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
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <mcheck.h>

#include <glib.h>
#include <gts.h>

#include <wmpi.h>

#include "bem3d.h"
#include "bem3d-private.h"

gint lookup_func_block(gint i, gint j, BEM3DMeshData *f,
		       GArray *s, GArray *ds) ;
gint lookup_func_block_c(gint i, gint j, 
			 BEM3DMeshData *f,
			 GArray *s, GArray *ds) ;

gint lookup_func_block(gint i, gint j, 
		       BEM3DMeshData *f,
		       GArray *s, GArray *ds)

{
  gdouble *x ;

  x = bem3d_mesh_data_get(f, i) ;

  if ( f != NULL ) {
    g_array_index(s, gdouble, j) = x[0] ;
    g_array_index(ds, gdouble, j) = x[1] ;
    return 0 ;
  }

  g_error("%s: vertex %d not found in data", __FUNCTION__, i) ;

  return 0 ;
}

gint lookup_func_block_c(gint i, gint j, 
			 BEM3DMeshData *f,
			 GArray *s, GArray *ds)

{
  gdouble *x ;

  x = bem3d_mesh_data_get(f, i) ;

  if ( f != NULL ) {
    g_array_index(s, gdouble, 2*j) = x[0] ;
    g_array_index(s, gdouble, 2*j+1) = x[1] ;
    g_array_index(ds, gdouble, 2*j) = x[2] ;
    g_array_index(ds, gdouble, 2*j+1) = x[3] ;
    return 0 ;
  }

  g_error("%s: vertex %d not found in data", __FUNCTION__, i) ;

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *s, *m ;
  BEM3DWorkspace *work ;
  GPtrArray *meshes, *mdata, *overrides ;
  GArray *field ;
  GtsVertex *x ;
  GtsFile *fp ;
  gint i, j, lineno, nc ;
  gchar ch, p[32], *progname, line[2048] ;
  gchar *ipfile, *sfile, *opfile ;
  BEM3DParameters gdata ;
  BEM3DLookupFunc lfunc ;
  BEM3DMeshData *sdata ;
  gint nmp, ndp ;
  gdouble k, M, result[32] = {0.0} ;
  gboolean helmholtz, point_input ;
  BEM3DFunction *efunc ;
  FILE *input, *output ;
  BEM3DConfiguration *config ;

  wmpi_initialize(&argc, &argv) ;
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  /*computational geometry*/
  meshes = g_ptr_array_new() ; mdata = g_ptr_array_new() ;
  /*surface for radiated field*/
  s = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		   gts_edge_class(), gts_vertex_class()) ;

  sprintf(p, "P%03d: ", wmpi_rank()) ;
  bem3d_logging_init(stderr, p, G_LOG_LEVEL_MESSAGE, wmpi_shutdown) ;
  bem3d_shapefunc_lookup_init() ;

  lfunc = (BEM3DLookupFunc)lookup_func_block ;
  k = 0.0 ; M = 0.0 ; helmholtz = FALSE ;
  ipfile = opfile = sfile = NULL ;
  
  point_input = FALSE ;
  efunc = NULL ;
  overrides = g_ptr_array_new() ;
  nmp = ndp = 0 ;

  /*default configuration*/
  config = bem3d_configuration_new() ;
  bem3d_parameters_init(&gdata) ;

  while ( (ch = getopt(argc, argv, "hC:d:F:i:k:M:o:s:v:X")) != EOF ) {
    switch (ch) {
    default:
          case 'h':
      fprintf(stderr, 
	      "%s: compute the radiated field from a BEM3D solution\n\n",
	      progname) ;
      fprintf(stderr, "Usage: %s <options>\n", progname) ;
      fprintf(stderr, 
	      "Options:\n"
	      "        -h (print this message and exit)\n"
	      "        -C <configuration file name>\n"
	      "        -d <data file name> (can be repeated, need one per \n"
	      "           input file)\n"
	      "        -F <function file> function to add to computed field\n"
	      "           when pointwise calculation is selected with -X\n"
	      "        -i <bem3d input file> (can be repeated)\n"
	      "        -k # (wave number for Helmholtz calculation)\n"
	      "        -M # (Mach number for convected Helmholtz equation)\n"
	      "        -o <output file name>\n"
	      "        -s <surface file name> (a BEM3D mesh file of points "
	      "where the\n"
	      "           field will be computed)\n"
	      "        -v <expression> set variables in the function "
	      "specified\n"
	      "           with -F option\n"
	      "        -X treat surface file as list of points in bem3d-dump\n"
	      "           format\n") ;
      return 0 ;
      break ;
    case 'C':
      bem3d_configuration_init() ;
      bem3d_configuration_read(config, optarg) ;
      break ;
    case 'd':
      input = file_open(optarg, "-", "r", stdin) ;
      bem3d_mesh_data_read(&sdata, input, 0) ;
      g_ptr_array_add(mdata, sdata) ;
      ndp += bem3d_mesh_data_node_number(sdata) ;
      file_close(input) ;
      break ;
    case 'F':
      efunc = bem3d_function_new(bem3d_function_class()) ;
      input = file_open(optarg, "-", "r", stdin) ;

      fp = gts_file_new(input) ;
      bem3d_function_read(efunc, fp) ;
      file_close(input) ;
      break ;
    case 'i': 
      ipfile = g_strdup(optarg) ;
      m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
			 gts_edge_class(), gts_vertex_class()) ;
      input = file_open(ipfile, "-", "r", stdin) ;
      fp = gts_file_new(input) ;
      bem3d_mesh_read(m, fp) ;
      file_close(input) ;
      g_ptr_array_add(meshes, m) ;
      nmp += bem3d_mesh_node_number(m) ;
      break ;
    case 'k':
      k = atof(optarg) ;
      /*use Helmholtz equation*/
      helmholtz = TRUE ;
      if ( wmpi_rank() == 0 ) 
	  fprintf(stderr, "%s: wave equation, k=%f\n", progname, k) ;
      lfunc = (BEM3DLookupFunc)lookup_func_block_c ;
      break ;
    case 'M': sscanf(optarg, "%lg", &M) ; break ;
    case 'o': opfile = g_strdup(optarg) ; break ;
    case 's': sfile = g_strdup(optarg) ; break ;
    case 'v': g_ptr_array_add(overrides, g_strdup(optarg)) ; break ;
    case 'X': point_input = TRUE ; break ;
    }
  }

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  if ( nmp > ndp ) {
    fprintf(stderr, "%s: number of data points (%d) does not match "
	    "number of nodes (%d)\n", progname, ndp, nmp) ;
    return 1 ;
  }

  if ( helmholtz && config->gfunc.func == bem3d_greens_func_laplace ) {
    fprintf(stderr, "%s: cannot use Laplace Greens function and specify"
	    "wavenumber (-k %lg)\n", progname, k) ;
    
    return 1 ;
  }

  if ( mdata->len != meshes->len ) {
    fprintf(stderr, "%s: require one data file per mesh (can be duplicates)\n",
	    progname) ;
    return 1 ;
  }
  
  bem3d_parameters_wavenumber(&gdata) = k ;
  bem3d_parameters_mach_number(&gdata) = M ;
  bem3d_parameters_quadrature_tol(&gdata) = config->quad_tol ;

  /* if ( M != 0.0 ) { */
  /*   config->gfunc.func = bem3d_greens_func_convected_helmholtz ; */
  /* } */

  if ( opfile != NULL )
    output = file_open(opfile, "-", "w", stdout) ;
  else 
    output = stdout ;
  
  if ( sfile == NULL )
    input = stdin ;
  else
    input = file_open(sfile, "-", "r", stdin) ;

  if ( efunc != NULL ) {
    for ( i = 0 ; i < overrides->len ; i ++ ) {
      bem3d_function_insert_string(efunc,
				   (gchar *)g_ptr_array_index(overrides,i)) ;
    }
    bem3d_function_expand_functions(efunc) ;
  }

  work = bem3d_workspace_new() ;
  
  if ( point_input ) {
    x = gts_vertex_new(gts_vertex_class(), 0, 0, 0) ;
    lineno = 1 ;
    field = g_array_new(FALSE, FALSE, sizeof(gdouble)) ;

    
    while ( (nc = fscanf(input, "%[^\n]c", line)) != EOF && ( nc != 0 ) ) {
      nc = sscanf(line, "%d %lg %lg %lg", &j,
		  &(GTS_POINT(x)->x),
		  &(GTS_POINT(x)->y),
		  &(GTS_POINT(x)->z)) ;
      if ( nc != 4 ) {
	fprintf(stderr, "%s: cannot parse line %d\n  %s\n",
		progname, lineno, line) ;
	exit(1) ;
      }
      for ( i = 0 ; i < meshes->len ; i ++ ) {
	bem3d_mesh_radiation_point(g_ptr_array_index(meshes,i),
				   config, &gdata,
				   lfunc, g_ptr_array_index(mdata,i),
				   GTS_POINT(x), field, work) ;
      }
      if ( efunc != NULL ) {
	memset(result, 0, 32*sizeof(gdouble)) ;
	i = bem3d_function_eval_point(efunc, GTS_POINT(x), NULL, 0, result, 8) ;
      }
      fprintf(output, "%d", j) ;
      for ( i = 0 ; i < field->len ; i ++ )
	fprintf(output, " %lg",
		g_array_index(field,gdouble,i) + result[i]) ;
      fprintf(output, "\n") ;
      lineno ++ ;
      if ( (nc = fscanf(input, "%*c")) == EOF ) break ;
    }
  } else {
    fp = gts_file_new(input) ;
    bem3d_mesh_read(s, fp) ;
    file_close(input) ;
    
    fprintf(stderr, "%s: surface file %s read\n", progname, sfile) ;
    fprintf(stderr, "%s: nodes: %d; elements: %d\n", progname,
	    bem3d_mesh_node_number(s), bem3d_mesh_element_number(s)) ;
  
    sdata = bem3d_mesh_data_new(s, 8) ;
    bem3d_mesh_data_clear(sdata) ;
    
    if ( efunc != NULL ) {
      bem3d_function_apply_mesh(efunc, s, sdata, NULL) ;
    }
    
    for ( i = 0 ; i < meshes->len ; i ++ ) {
      bem3d_mesh_radiation_mesh(g_ptr_array_index(meshes,i),
				config, &gdata,
				lfunc, g_ptr_array_index(mdata,i),
				s, sdata, work) ;
    }

    bem3d_mesh_data_write(sdata, output, NULL) ;
  }
  
  file_close(output) ;
  wmpi_shutdown() ;

  return 0 ;
}
