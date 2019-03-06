/* bem3d-assemble.c
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

#include "trace.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#ifdef HAVE_GQR
#include <gqr.h>
#else /*HAVE_GQR*/
#error "GQR is required for quadrature rules"
#endif /*HAVE_GQR*/

/* #define FMM_DIRECT_CONSTANT */

#define ADATA_MESHES        0
#define ADATA_CONFIG        1
#define ADATA_GDATA         2
#define ADATA_ROW_A         3
#define ADATA_ROW_B         4
#define ADATA_N_NODES       5
#define ADATA_WIDTH         6
#define ADATA_POINT         7
#define ADATA_INDEX         8
#define ADATA_MESH          9
#define ADATA_OUTPUT       10
#define ADATA_IMIN         11
#define ADATA_IMAX         12
#define ADATA_C            13
#define ADATA_GRADIENT_W_S 14
#define ADATA_GRADIENT_W_T 15
#define ADATA_GRADIENT_W_I 16
#define ADATA_GRADIENT_N_W 17

#define ADATA_SIZE     32

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
  BEM3DConfiguration *config = adata[ADATA_CONFIG] ;
  BEM3DParameters *gdata = adata[ADATA_GDATA] ;
  gdouble *a = adata[ADATA_ROW_A] ;
  gdouble *b = adata[ADATA_ROW_B] ;
  GtsVertex *x = (GtsVertex *)adata[ADATA_POINT] ;
  gint i = *((gint *)(adata[ADATA_INDEX])) ;

  bem3d_element_assemble_equations_direct(e, GTS_POINT(x), i, config, gdata,
					  a, b) ;
  return ;
}

static void _assemble_element_ch(BEM3DElement *e, gpointer adata[])

{
  BEM3DConfiguration *config = adata[ADATA_CONFIG] ;
  BEM3DParameters *gdata = adata[ADATA_GDATA] ;
  gdouble *a = adata[ADATA_ROW_A] ;
  gdouble *b = adata[ADATA_ROW_B] ;
  GtsVertex *x = (GtsVertex *)adata[ADATA_POINT] ;
  gint i = *((gint *)(adata[ADATA_INDEX])) ;
  gint j, nc = 2, idx ;
  static GArray *G = NULL ;
  static GArray *dG = NULL ;

  if ( G == NULL ) {
    G = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
    dG = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  }
  
  bem3d_element_assemble_equations(e, GTS_POINT(x), config, gdata, G, dG) ;

  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
    idx = bem3d_element_global_index(e,j) ;
    a[idx*nc+0] += g_array_index(dG, gdouble, 6*j+0) ;
    a[idx*nc+1] += g_array_index(dG, gdouble, 6*j+1) ;
    b[idx*nc+0] += g_array_index(G, gdouble, 6*j+0) ;
    b[idx*nc+1] += g_array_index(G, gdouble, 6*j+1) ;
  }

  if ( !bem3d_element_has_node(e, x) ) {
    for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
      idx = bem3d_element_global_index(e,j) ;
      /*d2G/dndn1 ...*/
      a[idx*nc+0] += g_array_index(dG, gdouble, 6*j+4) ;
      a[idx*nc+1] += g_array_index(dG, gdouble, 6*j+5) ;
      a[i*nc+0] -= g_array_index(dG, gdouble, 6*j+4) ;
      a[i*nc+1] -= g_array_index(dG, gdouble, 6*j+5) ;

      a[idx*nc+0] += g_array_index(dG, gdouble, 6*j+2) ;
      a[idx*nc+1] += g_array_index(dG, gdouble, 6*j+3) ;
    }
    
    return ;
  }
  
  /*handle the hypersingular self terms for the element*/
  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
    /*k^2 G ...*/
    a[i*nc+0] += g_array_index(dG, gdouble, 6*j+2) ;
    a[i*nc+1] += g_array_index(dG, gdouble, 6*j+3) ;
  }
  
  return ;
}

static void _assemble_rows(gint i, GtsVertex *v, gpointer adata[])

{
  GPtrArray *meshes = adata[ADATA_MESHES] ;
  BEM3DConfiguration *config = adata[ADATA_CONFIG] ;
  BEM3DParameters *gdata = adata[ADATA_GDATA] ;
  gdouble *a = adata[ADATA_ROW_A] ;
  gdouble *b = adata[ADATA_ROW_B] ;
  gint np = *(gint *)adata[ADATA_N_NODES] ;
  gint nc = *(gint *)adata[ADATA_WIDTH] ;
  FILE *output = (FILE *)adata[ADATA_OUTPUT] ;
  guint imin = *(guint *)adata[ADATA_IMIN] ;
  guint imax = *(guint *)adata[ADATA_IMAX] ;
  gdouble *C = (gdouble *)adata[ADATA_C] ;
  gint im = *(gint *)adata[ADATA_MESH] ;
  gdouble *J, *Ji, L[32], dLds[32], dLdt[32] ;
  BEM3DElement *e = NULL ;
  BEM3DGreensFunction gfunc ;
  BEM3DShapeFunc shfunc ;
  gdouble *nx ;
  gsl_complex lambda, al ;
  gint j ;
  
  if ( (i < imin) || (i > imax ) ) return ;

  /* fprintf(stderr, "%lg ", C[i]) ; */
  
  nx = bem3d_parameters_normal(gdata) ;
  nx[0] = nx[1] = nx[2] = 0.0 ;
  memset(a, 0, nc*np*sizeof(gdouble)) ; memset(b, 0, nc*np*sizeof(gdouble)) ;

  adata[ADATA_POINT] = v ; adata[ADATA_INDEX] = &i ;
  GSL_SET_COMPLEX(&lambda,
		  bem3d_parameters_lambda_real(gdata),
		  bem3d_parameters_lambda_imag(gdata)) ;
  GSL_SET_COMPLEX(&al,
		  bem3d_parameters_coupling_real(gdata),
		  bem3d_parameters_coupling_imag(gdata)) ;
  
  gfunc = config->gfunc ;
  /*local normal for hypersingular formulations*/
  if ( bem3d_greens_function_compute_normal(&gfunc) ) {
    if ( bem3d_greens_function_compute_jacobian(&gfunc) ) {
      g_assert(meshes->len == 1) ;
      J = bem3d_parameters_jacobian(gdata) ;
      Ji = bem3d_parameters_jacobian_inverse(gdata) ;
      e = bem3d_element_from_node(g_ptr_array_index(meshes,0), v, i) ;
      j = bem3d_element_find_node(e, v) ;
      shfunc = bem3d_element_shape_func(e) ;
      shfunc(bem3d_element_node_xi(e,j), bem3d_element_node_eta(e,j),
	     L, dLds, dLdt, NULL) ;
      bem3d_element_jacobian_matrix_normal(e, dLds, dLdt, 
					   nx, J, Ji) ;      
      shfunc = bem3d_element_node_func(e) ;
      shfunc(bem3d_element_node_xi(e,j), bem3d_element_node_eta(e,j),
	     L, dLds, dLdt, NULL) ;
    }
    else {
      bem3d_node_normal(g_ptr_array_index(meshes,im), i, nx,
			BEM3D_AVERAGE_MWA) ;
    }
  }

  if ( bem3d_greens_function_compute_gradient(&gfunc) ) {
    if ( bem3d_parameters_gradient(gdata) == NULL ) {
      bem3d_parameters_gradient(gdata) = bem3d_operator_new() ;
    }
    bem3d_operator_gradient(g_ptr_array_index(meshes,0), i,
			    bem3d_parameters_gradient(gdata),
			    BEM3D_AVERAGE_MWA) ;
  }
  
  /* if ( i == 0 ) _bem3d_set_trace(0) ; */

  if ( gfunc.size == nc ) {
    /*"standard" Green's function with no special treatment for self
      terms*/
    for ( j = 0 ; j < meshes->len ; j ++ )
      bem3d_mesh_foreach_element(g_ptr_array_index(meshes,j),
				 (BEM3DElementFunc)_assemble_element, 
				 adata) ;
  } else {
    /*we need to do something different here*/
    /* fprintf(stderr, "Hello\n") ; */
    for ( j = 0 ; j < meshes->len ; j ++ )
      bem3d_mesh_foreach_element(g_ptr_array_index(meshes,j),
				 (BEM3DElementFunc)_assemble_element_ch, 
				 adata) ;
  }
  
  /* _bem3d_unset_trace(0) ; */

  /*this works for hypersingular and standard formulations as long as
    lambda is set to zero by default for standard methods*/
  /* fprintf(stderr, "%lg ", a[i*nc+0]) ; */
  /* C[i] = 0.5 ; */
  a[i*nc+0] -= C[i]*GSL_REAL(al) ;
  a[i*nc+1] -= C[i]*GSL_IMAG(al) ;
  b[i*nc+0] += GSL_REAL(lambda)*C[i] ;
  b[i*nc+1] += GSL_IMAG(lambda)*C[i] ;

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
  gboolean meshes_closed ;
  GtsFile *fp ;
  GTimer *t ;
  gint i, j, np, nc, w, order ;
  gchar ch, *ipfile, *opfile, p[32], *progname ;
  BEM3DLookupFunc dgfunc ;
  BEM3DParameters gdata ;
  gpointer adata[ADATA_SIZE] ;
  gdouble k ;
  gdouble *a, *b, *C, *unit, r_correct ;
  gint imin, imax ;
  GLogLevelFlags loglevel ;
  FILE *input, *output ;
  BEM3DConfiguration *config ;
  BEM3DGreensFunction gfunc ;

  wmpi_initialize(&argc, &argv) ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  /*computational geometry*/
  meshes = g_ptr_array_new() ;
  m = NULL ;
  
  dgfunc = (BEM3DLookupFunc)bem3d_lookup_func_unit ;

  order = 7 ; anorm = BEM3D_AVERAGE_MWA ;
  bem3d_parameters_init(&gdata) ;

  meshes_closed = TRUE ;
  
  bem3d_parameters_wavenumber(&gdata) = 0.0 ;
  /*this is required to ensure correct evaluation for
    non-hypersingular methods*/
  bem3d_parameters_lambda_real(&gdata) = 
    bem3d_parameters_lambda_imag(&gdata) = 0.0 ;
  bem3d_parameters_coupling_real(&gdata) = 1.0 ;
  bem3d_parameters_coupling_imag(&gdata) = 0.0 ;
  w = 1 ;
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
  if ( bem3d_greens_function_compute_lambda(&(config->gfunc)) ) {
    bem3d_parameters_lambda_real(&gdata) = 0.0 ;
    bem3d_parameters_lambda_imag(&gdata) = 0.0 ;
      /* -1.0/bem3d_parameters_wavenumber(&gdata) ; */
    bem3d_parameters_coupling_real(&gdata) = 1.0 ;
    bem3d_parameters_coupling_imag(&gdata) = 0.0 ;
      /* -1.0/bem3d_parameters_wavenumber(&gdata); */
  }

  for ( (i = 0), (np = 0) ; i < meshes->len ; i ++ ) {
    nc = bem3d_mesh_node_number(g_ptr_array_index(meshes,i)) ;

    if ( wmpi_rank() == 0 ) {
      fprintf(stderr, "%s: ", progname) ;
      if ( gts_surface_is_closed(GTS_SURFACE(g_ptr_array_index(meshes,i))) ) {
	fprintf(stderr, "closed ") ;
      } else {
	meshes_closed = FALSE ;
	fprintf(stderr, "open ") ;
      }

      fprintf(stderr, "mesh %d, %d collocation points\n", i, nc) ;
    }
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
  
    adata[ADATA_MESHES] = meshes ;
    adata[ADATA_CONFIG] = config ;
    adata[ADATA_GDATA] = &gdata ;
    adata[ADATA_ROW_A] = a ; adata[ADATA_ROW_B] = b ;
    adata[ADATA_N_NODES] = &np ; adata[ADATA_WIDTH] = &w ;
    adata[ADATA_OUTPUT] = output ;
    adata[ADATA_IMIN] = &imin ; adata[ADATA_IMAX] = &imax ; 
    adata[ADATA_C] = C ;
    adata[ADATA_MESH] = &i ;

    fprintf(stderr, "%s: starting assembly: t=%f\n", 
	    progname, g_timer_elapsed(t, NULL)) ;

    k = bem3d_parameters_wavenumber(&gdata) ;
    bem3d_parameters_wavenumber(&gdata) = 0.0 ;  
    for ( i = 0 ; i < np ; i ++ ) C[i] = 1.0 ;
    gfunc = config->gfunc ;
    config->gfunc = greens_func_laplace ;
    if ( meshes_closed ) {
      for ( i = 0 ; i < meshes->len ; i ++ )
	bem3d_mesh_quad_dgdn(g_ptr_array_index(meshes,i),
			     config, &gdata, dgfunc, NULL, 
			     (BEM3DEquationFunc)equation_func_C, C) ;
    } else {
      /*if we have any open meshes, assume they are disjoint parts of the
       same surface and treat accordingly*/
      for ( i = 0 ; i < meshes->len ; i ++ )
	for ( j = 0 ; j < meshes->len ; j ++ )
	  bem3d_mesh_disjoint_quad_dgdn(g_ptr_array_index(meshes,i),
					g_ptr_array_index(meshes,j),
					config, &gdata, dgfunc, NULL, 
					(BEM3DEquationFunc)equation_func_C, C) ;
    }

    /* for ( i = 0 ; i < np ; i ++ ) fprintf(stderr, "%d %lg\n", i, C[i]) ; */

    /* fprintf(stderr, "%lg %lg %lg\n", C[21], C[22], C[23]) ; */

    /* exit(1) ; */
    
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
