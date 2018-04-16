/* bem3d-assemble.c
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <wmpi.h>

#include "bem3d.h"
#include "bem3d-private.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#ifdef HAVE_GQR
#include <gqr.h>
#else /*HAVE_GQR*/
#error "GQR is required for quadrature rules"
#endif /*HAVE_GQR*/

/* #define FMM_DIRECT_CONSTANT */

static gint skeleton_set_unit_sources(BEM3DMeshSkeleton *s, gdouble *dq)

{
  gint i, j ;
  gdouble *w ;

  memset(dq, 0, (s->ns)*sizeof(gdouble)) ;

  for ( i = 0 ; i < s->ns ; i ++ ) {
    w = &(s->w[i*s->ppe]) ;

    for ( j = 0 ; j < s->ppe ; j ++ ) dq[i] +=  w[j]*0.25*M_1_PI ;
  }

  return 0 ;
}

static gint equation_func_C(gint i, gint j,
			    gdouble *G, gdouble *dGdn, 
			    gint n, gdouble *C)
{
  C[i] += dGdn[0] ;

  return 0 ;
}

static void _assemble_element(BEM3DElement *e, gpointer adata[])

{
  BEM3DConfiguration *config = adata[1] ;
  BEM3DParameters *gdata = adata[2] ;
  gdouble *a = adata[3] ;
  gdouble *b = adata[4] ;
  gint np = *(gint *)adata[5] ;
  gint nc = *(gint *)adata[6] ;
  GtsVertex *x = (GtsVertex *)adata[7] ;
  static GArray *G = NULL ;
  static GArray *dGdn = NULL ;
  gint i, j, k, stride ;

  if ( G == NULL ) {
    G = g_array_new(FALSE, FALSE, sizeof(gdouble)) ;
    dGdn = g_array_new(FALSE, FALSE, sizeof(gdouble)) ;
  }

  g_array_set_size(G,bem3d_element_node_number(e)*nc) ; 
  g_array_set_size(dGdn,bem3d_element_node_number(e)*nc) ;
  bem3d_element_assemble_equations(e, GTS_POINT(x), config, gdata,
				   G, dGdn) ;
  stride = nc ;
  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    k = bem3d_element_global_index(e, i) ;
    for ( j = 0 ; j < stride ; j ++ ) {
      g_assert(k*stride+j < np*nc) ;
      b[k*stride+j] += g_array_index(G,gdouble,(stride*i+j)) ;
      a[k*stride+j] += g_array_index(dGdn,gdouble,(stride*i+j)) ;
    }
  }

  return ;
}

static void _assemble_rows(gint i, GtsVertex *v, gpointer adata[])

{
  GPtrArray *meshes = adata[0] ;
  BEM3DParameters *gdata = adata[2] ;
  gdouble *a = adata[3] ;
  gdouble *b = adata[4] ;
  gint np = *(gint *)adata[5] ;
  gint nc = *(gint *)adata[6] ;
  FILE *output = (FILE *)adata[10] ;
  guint imin = *(guint *)adata[11] ;
  guint imax = *(guint *)adata[12] ;
  gdouble *C = (gdouble *)adata[13] ;
  gint im = *(gint *)adata[14] ;
  gboolean hypersingular = *(gboolean *)adata[15] ;
  gdouble alpha ;
  gsl_complex lambda ;
  gint j ;

  if ( (i < imin) || (i > imax ) ) return ;

  for ( j = 0 ; j < nc*np ; j ++ ) a[j] = b[j] = 0.0 ;
  adata[7] = v ;
  alpha = bem3d_parameters_conditioning(gdata) ;
  lambda = *((gsl_complex *)(&(bem3d_parameters_lambda_real(gdata)))) ;

  /*local normal for hypersingular formulations*/
  if ( hypersingular ) 
    bem3d_node_normal(g_ptr_array_index(meshes,im), i,
		      bem3d_parameters_normal(gdata), BEM3D_AVERAGE_MWA) ;

  for ( j = 0 ; j < meshes->len ; j ++ )
    bem3d_mesh_foreach_element(g_ptr_array_index(meshes,j),
			       (BEM3DElementFunc)_assemble_element, 
			       adata) ;
  if ( hypersingular ) {
    a[i*nc+0] -= C[i]*(1.0-alpha) ;
    b[i*nc+0] += alpha*GSL_REAL(lambda)*C[i] ;
    b[i*nc+1] += alpha*GSL_IMAG(lambda)*C[i] ;
  } else
    a[i*nc+0] -= C[i] ;

  fprintf(output, "%d ", i) ;
  for ( j = 0 ; j < np*nc ; j ++ ) fprintf(output, " %1.16e", a[j]) ;
  fprintf(output, "\n") ;

  fprintf(output, "%d ", i) ;
  for ( j = 0 ; j < np*nc ; j ++ ) fprintf(output, " %1.16e", b[j]) ;
  fprintf(output, "\n") ;

  fflush(output) ;

  return ;
}

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *m ;
  BEM3DMeshSkeleton *skel ;
  BEM3DFMMMatrix *mtx ;
  BEM3DAverage anorm ;
  BEM3DQuadratureRule *quad ;
  BEM3DFMMWorkspace *work ;
  GPtrArray *meshes ;
  GtsFile *fp ;
  GTimer *t ;
  gint i, j, idx, np, nc, w, order ;
  gchar ch, *ipfile, *opfile, p[32], *progname ;
  BEM3DLookupFunc dgfunc ;
  BEM3DParameters gdata ;
  gpointer adata[16] ;
  gdouble k ;
  gdouble *a, *b, *C, *unit, r_correct ;
  guint imin, imax ;
  GLogLevelFlags loglevel ;
  FILE *input, *output ;
  BEM3DConfiguration *config ;
  BEM3DGreensFunction gfunc ;
  gboolean hypersingular ;

  wmpi_initialize(&argc, &argv) ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  /*computational geometry*/
  meshes = g_ptr_array_new() ;

  dgfunc = (BEM3DLookupFunc)bem3d_lookup_func_unit ;

  order = 7 ; anorm = BEM3D_AVERAGE_MWA ;

  bem3d_parameters_wavenumber(&gdata) = 0.0 ;
  w = 1 ; hypersingular = FALSE ;
  ipfile = opfile = NULL ;
  output = stdout ;
  loglevel = G_LOG_LEVEL_MESSAGE ;
  bem3d_shapefunc_lookup_init() ;
  sprintf(p, "P%03d:", wmpi_rank()) ;
  bem3d_logging_init(stderr, p, loglevel, wmpi_shutdown) ;

  /*default configuration*/
  bem3d_configuration_init() ;
  config = bem3d_configuration_new() ;

  while ( (ch = getopt(argc, argv, "hl:C:i:k:o:")) != EOF ) {
    switch (ch) {
    default:
    case 'h':
      if ( wmpi_rank() == 0 ) {
	fprintf(stderr, 
		"%s: assemble equations for BEM3D calculation\n\n",
		progname) ;
	fprintf(stderr, "Usage: %s <options>\n", progname) ;
	fprintf(stderr, 
		"Options:\n"
		"        -h (print this message and exit)\n"
		"        -C <configuration file>\n"
		"        -i <bem3d input file>\n"
		"        -k # (wave number for Helmholtz calculation\n"
/* 		"        -M # (Mach number for convected Helmholtz\n" */
		"        -o <output file name>\n") ;
      }
      wmpi_pause() ;
      wmpi_shutdown() ;
      return 0 ;
      break ;
    case 'l': loglevel = 1 << atoi(optarg) ; break ;
    case 'C': bem3d_configuration_read(config, optarg) ; break ;
    case 'i': 
      ipfile = g_strdup(optarg) ;
      m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
			 gts_edge_class(), gts_vertex_class()) ;
      input = file_open(ipfile, "-", "r", stdin) ;
      fp = gts_file_new(input) ;
      bem3d_mesh_read(m, fp) ;
      file_close(input) ;
      g_free(ipfile) ;
      g_ptr_array_add(meshes, m) ;
      break ;
    case 'k': bem3d_parameters_wavenumber(&gdata) = atof(optarg) ;
      break ;
    /* case 'M':  */
    /*   sscanf(optarg, "%lg", &M) ; */
    /*   break ; */
    case 'o': 
      if ( wmpi_process_number() == 1 ) 
	opfile = g_strdup(optarg) ; 
      else
	opfile = g_strdup_printf("%s-%04d", optarg, wmpi_rank()) ; 
      break ;
    }
  }

  bem3d_logging_init(stderr, p, loglevel, wmpi_shutdown) ;
  gqr_logging_init(stderr, p, loglevel, wmpi_shutdown) ;

  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;

  if ( (bem3d_parameters_wavenumber(&gdata) != 0.0) &&  
       ( config->gfunc.real == TRUE ) ) 
      g_error("%s: Green's function must be complex for "
	      "non-zero wavenumber (k=%lg)", __FUNCTION__, 
	      bem3d_parameters_wavenumber(&gdata)) ;

  if ( config->gfunc.real == FALSE ) w = 2 ;
  if ( config->gfunc.func == bem3d_greens_func_helmholtz_hs ) {
    hypersingular = TRUE ;
    bem3d_parameters_lambda_real(&gdata) = 0.0 ;
    bem3d_parameters_lambda_imag(&gdata) = 
      1.0/bem3d_parameters_wavenumber(&gdata) ;
    bem3d_parameters_conditioning(&gdata) = 0.1 ;
  }

  for ( (i = 0), (np = 0) ; i < meshes->len ; i ++ ) {
    nc = bem3d_mesh_node_number(g_ptr_array_index(meshes,i)) ;
    if ( wmpi_rank() == 0 )
      fprintf(stderr, "%s: mesh %d, %d collocation points\n", progname, i, nc) ;
    np += nc ;
  }
  
  if ( np == 0 ) {
    if ( wmpi_rank() == 0 ) {
      fprintf(stderr, "%s: no meshes or empty meshes read: exiting\n", 
	      progname) ;
    }
    wmpi_pause() ;
    wmpi_shutdown() ;
    return 0 ;
  }

  if ( wmpi_rank() == 0 )
      fprintf(stderr, "%s: mesh collocation points: %d\n", progname, np) ;

  t = g_timer_new() ; g_timer_start(t) ;

  if ( config->solver == BEM3D_SOLVER_DIRECT ) {
    if ( opfile != NULL ) {
      output = file_open(opfile, "-", "w", stdout) ;
    }

    a = (gdouble *)g_malloc(np*w*sizeof(gdouble)) ;
    b = (gdouble *)g_malloc(np*w*sizeof(gdouble)) ;
    C = (gdouble *)g_malloc(np*sizeof(gdouble)) ;

    wmpi_split_range(0, np, &imin, &imax) ;
    fprintf(stderr, "%s: P%d: nodes: %d--%d\n", 
	    progname, wmpi_rank(), imin, imax) ;
  
    adata[0] = meshes ;
    adata[1] = config ;
    adata[2] = &gdata ;
    adata[3] = a ; adata[4] = b ;
    adata[5] = &np ; adata[6] = &w ;
    adata[10] = output ;
    adata[11] = &imin ; adata[12] = &imax ; 
    adata[13] = C ;
    adata[14] = &i ;
    adata[15] = &hypersingular ;

    fprintf(stderr, "%s: starting assembly: t=%f\n", 
	    progname, g_timer_elapsed(t, NULL)) ;

    k = bem3d_parameters_wavenumber(&gdata) ;
    bem3d_parameters_wavenumber(&gdata) = 0.0 ;  
    for ( i = 0 ; i < np ; i ++ ) C[i] = 1.0 ;
    gfunc = config->gfunc ;
    config->gfunc = greens_func_laplace ;
    for ( i = 0 ; i < meshes->len ; i ++ )
      bem3d_mesh_quad_dgdn(g_ptr_array_index(meshes,i),
			   config, &gdata, dgfunc, NULL, 
			   (BEM3DEquationFunc)equation_func_C, C) ;

    config->gfunc = gfunc ;
    bem3d_parameters_wavenumber(&gdata) = k ;

    fprintf(output, "%s %u %d %d %u %u %d\n", 
	    bem3d_solver_name(config->solver), 
	    wmpi_rank(), np, w, imin, imax, 0) ;

    for ( i = 0 ; i < meshes->len ; i ++ )
      bem3d_mesh_foreach_node(g_ptr_array_index(meshes,i),
			      (BEM3DNodeFunc)_assemble_rows, adata) ;
    file_close(output) ;

    fprintf(stderr, "%s: assembly completed: t=%f\n", 
	    progname, g_timer_elapsed(t, NULL)) ;

    g_ptr_array_free(meshes, FALSE) ; g_free(progname) ;

    wmpi_pause() ;
    wmpi_shutdown() ;

    return 0 ;
  }

  fprintf(stderr, "%s: assembling for FMM solver %d: t=%f\n", 
	  progname, config->fmm, g_timer_elapsed(t, NULL)) ;

  if ( opfile != NULL ) {
    output = file_open(opfile, "-", "w", stdout) ;
  }

  C = (gdouble *)g_malloc0(np*sizeof(gdouble)) ;

  bem3d_mesh_index_range(m, &imin, &imax) ;

  fprintf(output, "%s %u %d %d %u %u %d\n", 
	  bem3d_solver_name(config->solver), 
	  wmpi_rank(), np, w, imin, imax, config->fmm) ;

#ifdef FMM_DIRECT_CONSTANT
  for ( i = 0 ; i < np ; i ++ ) C[i] = 1.0 ;

  gfunc = config->gfunc ;
  config->gfunc = greens_func_laplace ;

  for ( i = 0 ; i < meshes->len ; i ++ )
    bem3d_mesh_quad_dgdn(g_ptr_array_index(meshes,i),
			 config, &gdata, dgfunc, NULL, 
			 (BEM3DEquationFunc)equation_func_C, C) ;

#else /*FMM_DIRECT_CONSTANT*/
  m = g_ptr_array_index(meshes, 0) ;
  quad = bem3d_quadrature_rule_new(order, 1) ;

  order = config->skel_order ;
  r_correct = config->fmm_radius ;

  if ( order <= 7 ) 
    bem3d_quadrature_rule_gauss(NULL, bem3d_mesh_element_sample(m), quad, 
				NULL, NULL, &order) ;
  else
    bem3d_quadrature_rule_wx(NULL, bem3d_mesh_element_sample(m), quad, 
			     NULL, NULL, &order) ;    

  skel = bem3d_mesh_skeleton_new(m, order) ;
  bem3d_mesh_skeleton_init(skel, quad, anorm) ;

  gfunc = config->gfunc ;
  config->gfunc = greens_func_laplace ;

  mtx = bem3d_fmm_matrix_new(config->fmm, BEM3D_FMM_LAPLACE,
			     skel, config, NULL, r_correct) ;

  mtx->tol = config->fmm_tol ;

  unit = (gdouble *)g_malloc(skel->ns*sizeof(gdouble)) ;

  work = bem3d_fmm_workspace_alloc(config->fmm, skel) ;

  skeleton_set_unit_sources(skel, unit) ;

  bem3d_fmm_calculate(mtx->solver, mtx->problem, NULL, mtx->skel, mtx->tol,
  		      NULL, unit, C, NULL, work) ;

  /*local correction terms in FMM integration*/
  for ( i = 0 ; i < np ; i ++ ) {
    C[i] += 1 ;
    for ( j = mtx->idxcorr[2*i+0] ; j < mtx->idxcorr[2*i+1] ; j ++ ) {
      C[i] += g_array_index(mtx->dgcorr, gdouble, j) ;
    }
  }
#endif /*FMM_DIRECT_CONSTANT*/

  fprintf(stderr, "%s: assembly completed: t=%f\n", 
	  progname, g_timer_elapsed(t, NULL)) ;

  for ( i = 0 ; i < np ; i ++ ) {
    fprintf(output, "%d %1.16e\n", i, C[i]) ;
  }

  file_close(output) ;

  wmpi_pause() ;
  wmpi_shutdown() ;

  return 0 ;
}
