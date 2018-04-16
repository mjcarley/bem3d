/* bem3d-source.c
 * 
 * Copyright (C) 2006 Michael Carley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <mcheck.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

#include "tools.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

static void parse_amplitude(gchar *s, BEM3DParameters *p)

{
  gchar **tokens ;

  tokens = g_strsplit(s, " ", -1) ;
  if ( tokens[0] != NULL ) 
    bem3d_parameters_amplitude_real(p) = atof(tokens[0]) ;
  if ( tokens[1] != NULL ) 
    bem3d_parameters_amplitude_imag(p) = atof(tokens[1]) ;

  return ;
}

static gint source_func_helmholtz(gint i, GtsVertex *v, gpointer data[],
				  gdouble *f)

{
  GtsPoint *s = data[BEM3D_TOOLS_DATA_POSITION] ;
  BEM3DParameters *p = data[BEM3D_TOOLS_DATA_PARAMETERS] ;
  gdouble R, R3, k, C, S, fi, fr, Ar, Ai ;
  GtsVector r ;

  gts_vector_init(r, s, GTS_POINT(v)) ;
  R = gts_vector_norm(r) ;
  k = bem3d_parameters_wavenumber(p) ;
  Ar = bem3d_parameters_amplitude_real(p) ;
  Ai = bem3d_parameters_amplitude_imag(p) ;

  C = cos(k*R) ; S = sin(k*R) ;

  f[0] += Ar*C/R*0.25*M_1_PI - Ai*S/R*0.25*M_1_PI ; 
  f[1] += Ar*S/R*0.25*M_1_PI + Ai*C/R*0.25*M_1_PI ;

  fr = -(C + k*R*S) ; fi = k*R*C - S ;

  R3 = R*R*R ;
  f[2] += Ar*0.25*M_1_PI/R3*r[0]*fr - Ai*0.25*M_1_PI/R3*r[0]*fi ;
  f[3] += Ar*0.25*M_1_PI/R3*r[0]*fi + Ai*0.25*M_1_PI/R3*r[0]*fr ;
  f[4] += Ar*0.25*M_1_PI/R3*r[1]*fr - Ai*0.25*M_1_PI/R3*r[1]*fi ;
  f[5] += Ar*0.25*M_1_PI/R3*r[1]*fi + Ai*0.25*M_1_PI/R3*r[1]*fr ;
  f[6] += Ar*0.25*M_1_PI/R3*r[2]*fr - Ai*0.25*M_1_PI/R3*r[2]*fi ;
  f[7] += Ar*0.25*M_1_PI/R3*r[2]*fi + Ai*0.25*M_1_PI/R3*r[2]*fr ;

  return 8 ;
}

static gint source_func_plane(gint i, GtsVertex *v, gpointer data[],
			      gdouble *f)

{
  BEM3DParameters *p = data[BEM3D_TOOLS_DATA_PARAMETERS] ;
  gdouble x, k, C, S, Ar, Ai ;

  k = bem3d_parameters_wavenumber(p) ;
  Ar = bem3d_parameters_amplitude_real(p) ;
  Ai = bem3d_parameters_amplitude_imag(p) ;
  x = GTS_POINT(v)->x ;

  C = cos(k*x) ; S = sin(k*x) ;

  f[0] += Ar*C - Ai*S ;      f[1] += Ai*C + Ar*S ;

  f[2] += k*(-Ar*S - Ai*C) ; f[3] += k*(-Ai*S + Ar*C) ;

  return 8 ;
}

static gint source_func_laplace(gint i, GtsVertex *v, gpointer data[],
				gdouble *f)

{
  GtsPoint *s = data[BEM3D_TOOLS_DATA_POSITION] ;
  BEM3DParameters *p = data[BEM3D_TOOLS_DATA_PARAMETERS] ;
  gdouble R, R3, A ;
  GtsVector r ;

  gts_vector_init(r, s, GTS_POINT(v)) ;
  R = gts_vector_norm(r) ;
  
  A = bem3d_parameters_amplitude_real(p) ;

  f[0] += A*1.0/R*0.25*M_1_PI ; 
  R3 = R*R*R ;

  f[1] += -A*0.25*M_1_PI/R3*r[0] ; 
  f[2] += -A*0.25*M_1_PI/R3*r[1] ;
  f[3] += -A*0.25*M_1_PI/R3*r[2] ;

  return 4 ;
}

static gint source_func_velocity(gint i, GtsVertex *v, gpointer data[],
				 gdouble *f)

{
  GtsVector *U = data[BEM3D_TOOLS_DATA_VELOCITY] ;

  f[1] += (*U)[0] ; f[2] += (*U)[1] ; f[3] += (*U)[2] ;

  return 4 ;
}

static gint source_func(gint i, GtsVertex *v, gpointer data[])

{
  gint (*func)(gint, GtsVertex *, gpointer *, gdouble *) = 
    data[BEM3D_TOOLS_DATA_FUNC] ;
  FILE *output = data[BEM3D_TOOLS_DATA_OUTPUT] ;
  gdouble f[64] ;
  gint j, nf ;
  
  for ( j = 0 ; j < 64 ; j ++ ) f[j] = 0.0 ;

  nf = func(i, v, data, f) ;

  fprintf(output, "%d", i) ;

  for ( j = 0 ; j < nf ; j ++ )
    fprintf(output, " %1.16e", f[j]) ;
  fprintf(output, "\n") ;

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *m ;
  BEM3DParameters param ;
  GtsFile *fp ;
  gchar *ipfile, *opfile, *srcfile ;
  gchar ch, p[32], *progname ;
  gpointer sdata[BEM3D_TOOLS_DATA_WIDTH], sfunc ;
  FILE *input, *output ;
  GtsPoint *s ;
  GtsVector U ;
  BEM3DConfiguration *config ;
  gint nf ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  p[0] = '\0' ;
  bem3d_logging_init(stderr, p, G_LOG_LEVEL_MESSAGE, NULL) ;
  bem3d_shapefunc_lookup_init() ;

  /*default test case is a Laplace equation (k=0)*/
  bem3d_parameters_wavenumber(&param) = 0.0 ;
  bem3d_parameters_amplitude_real(&param) = 1.0 ;
  bem3d_parameters_amplitude_imag(&param) = 0.0 ;

  ipfile = opfile = srcfile = NULL ;
  m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
		     gts_edge_class(), gts_vertex_class()) ;

  nf = 4 ; sfunc = source_func_laplace ;
  sdata[BEM3D_TOOLS_DATA_PARAMETERS] = &param ;
  sdata[BEM3D_TOOLS_DATA_MESH] = m ;
  sdata[BEM3D_TOOLS_DATA_POSITION] = s = 
    gts_point_new(gts_point_class(), 0, 0, 0) ;
  sdata[BEM3D_TOOLS_DATA_VELOCITY] = &U ;

  U[0] = U[1] = U[2] = 0.0 ;

  /*default configuration*/
  bem3d_configuration_init() ;
  config = bem3d_configuration_new() ;

  while ( (ch = getopt(argc, argv, "A:C:hi:k:o:pS:x:y:z:u:v:w:")) != EOF ) {
    switch (ch) {
    default:
    case 'h':
	fprintf(stderr, 
		"%s: compute field on BEM3D surface for boundary conditions\n"
		"for specified source type (deprecated, now better to use\n"
		"bem3d-function with analytically specified field)\n\n",
		progname) ;
	fprintf(stderr, "Usage: %s <options>\n", progname) ;
	fprintf(stderr, 
		"Options:\n"
		"        -h (print this message and exit)\n"
		"        -A # (amplitude, possibly complex)\n"
		"        -C <configuration file>\n"
		"        -i <bem3d input file>\n"
		"        -k # (wave number for Helmholtz calculation)\n"
		"        -o <output file name> (for mesh block data)\n"
		"        -p (plane wave incident field)\n"
		"        -s <source strength (real or complex)>\n"
		"        -S <source file name>\n"
		"        -u # (x component of velocity)\n"
		"        -v # (y component of velocity)\n"
		"        -w # (z component of velocity)\n"
		"        -x # (coordinate of point source)\n"
		"        -y # (coordinate of point source)\n"
		"        -z # (coordinate of point source)\n") ;
      return 0 ;
      break ;
    case 'A': parse_amplitude(optarg, &param) ; break ;
    case 'C': bem3d_configuration_read(config, optarg) ; break ;
    case 'i': ipfile = g_strdup(optarg) ; break ;
    case 'k': 
      bem3d_parameters_wavenumber(&param) = atof(optarg) ; break ;
    /* case 'M':  */
    /*   sscanf(optarg, "%lg", &bem3d_parameters_mach_number(&param)) ; */
    /*   break ; */
    case 'o': opfile = g_strdup(optarg) ; break ;
    case 'p': /*plane wave incident field*/
      sfunc = source_func_plane ;
      break ;
    case 'S': srcfile = g_strdup(optarg) ; break ;
    case 'u': U[0] = atof(optarg) ; sfunc = source_func_velocity ;
      break ;
    case 'v': U[1] = atof(optarg) ; sfunc = source_func_velocity ;
      break ;
    case 'w': U[2] = atof(optarg) ; sfunc = source_func_velocity ;
      break ;
    case 'x': sscanf(optarg, "%lg", &(s->x)) ; break ;
    case 'y': sscanf(optarg, "%lg", &(s->y)) ; break ;
    case 'z': sscanf(optarg, "%lg", &(s->z)) ; break ;
    }
  }

  if ( ipfile == NULL ) ipfile = g_strdup("-") ;
  if ( opfile == NULL ) opfile = g_strdup("-") ;

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  if ( ( bem3d_parameters_wavenumber(&param) != 0.0 ) ) {    
    if ( config->gfunc.real == TRUE ) 
      g_error("%s: Green's function must be complex for "
	      "non-zero wavenumber (k=%lg)", __FUNCTION__, 
	      bem3d_parameters_wavenumber(&param)) ;
    nf = 8 ;
    sfunc = source_func_helmholtz ;
  }

  input = file_open(ipfile, "-", "r", stdin) ;
  fp = gts_file_new(input) ;
  bem3d_mesh_read(m, fp) ;
  file_close(input) ;

  output = file_open(opfile, "-", "w", stdout) ;
  sdata[BEM3D_TOOLS_DATA_OUTPUT] = output ;
  sdata[BEM3D_TOOLS_DATA_FUNC] = sfunc ;

  fprintf(output, "%d %d BEM3DMeshData\n", bem3d_mesh_node_number(m), nf) ;
  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)source_func, sdata) ;

  file_close(output) ;

  return 0 ;
}
