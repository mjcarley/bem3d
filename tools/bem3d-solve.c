/* bem3d-solve.c
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
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <glib.h>
#include <gts.h>

#include <wmpi.h>
#include <sisl.h>

#include "bem3d.h"
#include "bem3d-private.h"

#include "tools.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

/* #define SLOW_VECTOR_ACCESS 0 */

/* #define FMM_MATRIX_CHECK 0 */

static void read_real_matrices(FILE *f, gint n, 
			       sisl_matrix_t *A, sisl_matrix_t *B)

{
  gint i, j, k, nc ;
  gdouble x ;

  nc = sisl_matrix_column_number(A) ;
  g_assert(sisl_matrix_row_number(B) == nc) ;

  for ( i = 0 ; i < n ; i ++ ) {
    fscanf(f, "%d", &k) ;
    if ( sisl_matrix_has_local_row(A,k) ) {
      for ( j = 0 ; j < nc ; j ++ ) {
	fscanf(f, "%lg", &x) ;
	sisl_matrix_set(A, k, j, x) ;
      }
      fscanf(f, "%d", &k) ;
      for ( j = 0 ; j < nc ; j ++ ) {
	fscanf(f, "%lg", &x) ;
	sisl_matrix_set(B, k, j, x) ;
      }      
    } else {
      for ( j = 0 ; j < nc ; j ++ ) fscanf(f, "%*g") ;
      fscanf(f, "%*d") ;
      for ( j = 0 ; j < nc ; j ++ ) fscanf(f, "%*g") ;
    }
  }

  return ;
}

static void read_complex_matrices(FILE *f, gint n, 
				  sisl_matrix_t *A, sisl_matrix_t *B)

{
  gint i, j, k, nc ;
  gsl_complex x ;

  nc = sisl_matrix_column_number(A) ;
  g_assert(sisl_matrix_row_number(B) == nc) ;

  for ( i = 0 ; i < n ; i ++ ) {
    fscanf(f, "%d", &k) ;
    if ( k < 0 )
      g_error("P%d: A matrix row index (%d) out of range at line %d", 
	      wmpi_rank(), k, i) ;
    if ( sisl_matrix_has_local_row(A,k) ) {
      for ( j = 0 ; j < nc ; j ++ ) {
	fscanf(f, "%lg %lg", &GSL_REAL(x), &GSL_IMAG(x)) ;
	sisl_matrix_set_complex(A, k, j, x) ;
      }
      fscanf(f, "%d", &k) ;
      if ( k < 0 )
	g_error("P%d: B matrix row index (%d) out of range at line %d", 
		wmpi_rank(), k, i) ;
      for ( j = 0 ; j < nc ; j ++ ) {
	fscanf(f, "%lg %lg", &GSL_REAL(x), &GSL_IMAG(x)) ;
	sisl_matrix_set_complex(B, k, j, x) ;
      }      
    } else {
      for ( j = 0 ; j < 2*nc ; j ++ ) fscanf(f, "%*g") ;
      fscanf(f, "%*d") ;
      for ( j = 0 ; j < 2*nc ; j ++ ) fscanf(f, "%*g") ;
    }
  }

  return ;
}

static void set_complex_vectors(sisl_vector_t *phi, sisl_vector_t *dphi,
				BEM3DMeshData *data)

{
  gint i ;
  gdouble *f ;

  g_assert(sisl_vector_length(phi) == sisl_vector_length(dphi)) ;

  for ( i = 0 ; i < sisl_vector_length(phi) ; i ++ ) {
    f = bem3d_mesh_data_get(data, i) ;
    sisl_vector_set_complex(phi, i, *((gsl_complex *)(&(f[0])))) ;
    sisl_vector_set_complex(dphi, i, *((gsl_complex *)(&(f[2])))) ;
  }

  return ;
}

static void set_real_vectors(sisl_vector_t *phi, sisl_vector_t *dphi,
			     BEM3DMeshData *data)

{
  gint i ;
  gdouble *f ;

  g_assert(sisl_vector_length(phi) == sisl_vector_length(dphi)) ;

  for ( i = 0 ; i < sisl_vector_length(phi) ; i ++ ) {
    g_assert( (f = bem3d_mesh_data_get(data, i)) != NULL ) ;
    sisl_vector_set(phi, i, f[0]) ;
    sisl_vector_set(dphi, i, f[1]) ;
  }

  return ;
}

static gint sisl_matrix_vector_mul_A(sisl_matrix_t *A, sisl_vector_t *v, 
				     sisl_vector_t *w)

{
  gpointer *data = sisl_matrix_user_defined_data(A) ;
  BEM3DFMMMatrix *mtx = data[0] ;
  BEM3DFMMWorkspace *work = data[1] ;
  BEM3DParameters *param = data[2] ;
  BEM3DMeshSkeleton *s ;
  gdouble *scratch = data[3] ;
  gint i, j, k, *idx ;
  gdouble *wt ;

  g_assert(sisl_matrix_density(A) == SISL_MATRIX_USER_DEFINED) ;

  s = mtx->skel ;
  /*set the source strengths from the input vector*/
  if ( sisl_is_real(A) ) {
    g_assert_not_reached() ;
  } else {
    memset(scratch, 0, 2*s->ns*sizeof(gdouble)) ;
#ifndef SLOW_VECTOR_ACCESS
    gdouble *vp ;

    vp = sisl_vector_data(v) ;

    for ( i = 0 ; i < s->ns ; i ++ ) {
      wt = &(s->w[i*s->ppe]) ; idx = &(s->idx[i*s->ppe]) ;

      for ( j = 0 ; j < s->ppe ; j ++ ) {
	scratch[2*i+0] +=  wt[j]*vp[2*idx[j]+0]*0.25*M_1_PI ;
	scratch[2*i+1] +=  wt[j]*vp[2*idx[j]+1]*0.25*M_1_PI ;
      }
    }

#ifdef FMM_MATRIX_CHECK
    fprintf(stderr, "%lg %lg %lg ", vp[0], vp[1], mtx->C[0]) ;
#endif /*FMM_MATRIX_CHECK*/

#else /*SLOW_VECTOR_ACCESS*/
    gsl_complex vc, wc, *cc ;
    for ( i = 0 ; i < s->ns ; i ++ ) {
      wt = &(s->w[i*s->ppe]) ; idx = &(s->idx[i*s->ppe]) ;

      for ( j = 0 ; j < s->ppe ; j ++ ) {
	wc = sisl_vector_get_complex(v,idx[j]) ; 

	scratch[i*2+0] +=  wt[j]*GSL_REAL(wc)*0.25*M_1_PI ;
	scratch[i*2+1] +=  wt[j]*GSL_IMAG(wc)*0.25*M_1_PI ;
      }
    }
#endif /*SLOW_VECTOR_ACCESS*/
  }

  bem3d_fmm_calculate(mtx->solver, mtx->problem, param, mtx->skel,
  		      mtx->tol, NULL, scratch,
  		      sisl_vector_data(w), NULL, work) ;

#ifdef FMM_MATRIX_CHECK
  fprintf(stderr, "%lg %lg ", 
	  (sisl_vector_data(w))[0], (sisl_vector_data(w))[1]) ;
#endif /*FMM_MATRIX_CHECK*/

  if ( sisl_is_real(A) ) {
    g_assert_not_reached() ;
  } else {

#ifndef SLOW_VECTOR_ACCESS
    gdouble *wp, *vp, *cp ;

    cp = &(g_array_index(mtx->dgcorr, gdouble, 0)) ;

    vp = sisl_vector_data(v) ; wp = sisl_vector_data(w) ;

    for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
      wp[2*i+0] -= mtx->C[i]*vp[2*i+0] ;
      wp[2*i+1] -= mtx->C[i]*vp[2*i+1] ;
    }
#ifdef FMM_MATRIX_CHECK
  fprintf(stderr, "%lg %lg ", 
	  (sisl_vector_data(w))[0], (sisl_vector_data(w))[1]) ;
#endif /*FMM_MATRIX_CHECK*/

    for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
      for ( j = mtx->idxcorr[2*i+0] ; j < mtx->idxcorr[2*i+1] ; j ++ ) {
	k = g_array_index(mtx->icorr, gint, j) ;
	wp[2*i+0] += vp[2*k+0]*cp[2*j+0] - vp[2*k+1]*cp[2*j+1] ;
	wp[2*i+1] += vp[2*k+0]*cp[2*j+1] + vp[2*k+1]*cp[2*j+0] ;
      }
    }

#ifdef FMM_MATRIX_CHECK
  fprintf(stderr, "%lg %lg\n", 
	  (sisl_vector_data(w))[0], (sisl_vector_data(w))[1]) ;
#endif /*FMM_MATRIX_CHECK*/

#else /*SLOW_VECTOR_ACCESS*/
    for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
      vc = sisl_vector_get_complex(v, i) ;
      wc = sisl_vector_get_complex(w, i) ;
      wc = gsl_complex_add(wc, gsl_complex_mul_real(vc, mtx->C[i])) ;
      for ( j = mtx->idxcorr[2*i+0] ; j < mtx->idxcorr[2*i+1] ; j ++ ) {
    	k = g_array_index(mtx->icorr, gint, j) ;
    	vc = sisl_vector_get_complex(v, k) ;
    	cc = (gsl_complex *)(&(g_array_index(mtx->dgcorr, gdouble, 2*j))) ;
      
    	wc = gsl_complex_add(wc,
    			     gsl_complex_mul(vc, *cc)) ;
      }
      sisl_vector_set_complex(w, i, wc) ;
    }
#endif /*SLOW_VECTOR_ACCESS*/
  }

  return 0 ;
}

static gint sisl_matrix_vector_mul_B(sisl_matrix_t *B, sisl_vector_t *v, 
				     sisl_vector_t *w)

{
  gpointer *data = sisl_matrix_user_defined_data(B) ;
  BEM3DFMMMatrix *mtx = data[0] ;
  BEM3DFMMWorkspace *work = data[1] ;
  BEM3DParameters *param = data[2] ;
  BEM3DMeshSkeleton *s ;
  gdouble *scratch = data[3], *wp, *vp ;
  gsl_complex vc, wc, *cc ;
  gint i, j, k, *idx ;
  gdouble *wt ;

  g_assert(sisl_matrix_density(B) == SISL_MATRIX_USER_DEFINED) ;

  s = mtx->skel ;
  /*set the source strengths from the input vector*/
  if ( sisl_is_real(B) ) {
    g_assert_not_reached() ;
  } else {
    memset(scratch, 0, 2*s->ns*sizeof(gdouble)) ;

#ifndef SLOW_VECTOR_ACCESS
    gdouble *vp ;

    vp = sisl_vector_data(v) ;

    for ( i = 0 ; i < s->ns ; i ++ ) {
      wt = &(s->w[i*s->ppe]) ; idx = &(s->idx[i*s->ppe]) ;

      for ( j = 0 ; j < s->ppe ; j ++ ) {
	/* zc = sisl_vector_get_complex(v,idx[j]) ;  */

	scratch[2*i+0] +=  wt[j]*vp[2*idx[j]+0]*0.25*M_1_PI ;
	scratch[2*i+1] +=  wt[j]*vp[2*idx[j]+1]*0.25*M_1_PI ;
      }
    }

#else /*SLOW_VECTOR_ACCESS*/

    for ( i = 0 ; i < s->ns ; i ++ ) {
      wt = &(s->w[i*s->ppe]) ; idx = &(s->idx[i*s->ppe]) ;

      for ( j = 0 ; j < s->ppe ; j ++ ) {
	vc = sisl_vector_get_complex(v,idx[j]) ; 

	scratch[i*2+0] +=  wt[j]*GSL_REAL(vc)*0.25*M_1_PI ;
	scratch[i*2+1] +=  wt[j]*GSL_IMAG(vc)*0.25*M_1_PI ;
      }
    }

#endif /*SLOW_VECTOR_ACCESS*/

  }

  bem3d_fmm_calculate(mtx->solver, mtx->problem, param, mtx->skel,
  		      mtx->tol,
  		      scratch, NULL,
  		      &(g_array_index(w->x,gdouble,0)), NULL, work) ;

  if ( sisl_is_real(B) ) {
    g_assert_not_reached() ;
  } else {
    /*correction terms to FMM multiplication*/

#ifndef SLOW_VECTOR_ACCESS
    wp = sisl_vector_data(w) ;
    vp = sisl_vector_data(v) ;
    /* wp = &(g_array_index(w->x,gdouble,0)) ; */
    /* vp = &(g_array_index(v->x,gdouble,0)) ; */

    for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
      wc = *((gsl_complex *)(&(wp[2*i+0]))) ;
      for ( j = mtx->idxcorr[2*i+0] ; j < mtx->idxcorr[2*i+1] ; j ++ ) {
	k = g_array_index(mtx->icorr, gint, j) ;
	vc = *((gsl_complex *)(&(vp[2*k+0]))) ;
	cc = (gsl_complex *)(&(g_array_index(mtx->gcorr, gdouble, 2*j))) ;
      
	wc = gsl_complex_add(wc, 
			     gsl_complex_mul(vc, *cc)) ;
      }
      wp[2*i+0] = GSL_REAL(wc) ; wp[2*i+1] = GSL_IMAG(wc) ;
    }

#else /*SLOW_VECTOR_ACCESS*/
    for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
      wc = sisl_vector_get_complex(w, i) ;
      for ( j = mtx->idxcorr[2*i+0] ; j < mtx->idxcorr[2*i+1] ; j ++ ) {
	k = g_array_index(mtx->icorr, gint, j) ;
	vc = sisl_vector_get_complex(v, k) ;
	cc = (gsl_complex *)(&(g_array_index(mtx->gcorr, gdouble, 2*j))) ;
      
	wc = gsl_complex_add(wc, 
			     gsl_complex_mul(vc, *cc)) ;
      }
      sisl_vector_set_complex(w, i, wc) ;
    }
#endif /*SLOW_VECTOR_ACCESS*/

  }

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  GTimer *t ;
  BEM3DMesh *m ;
  BEM3DMeshData *data ;
  BEM3DParameters param ;
  GPtrArray *meshes ;
  sisl_matrix_t *A, *B ;
  sisl_vector_t *phi, *dphi, *rhs ;
  sisl_solver_workspace_t *w ;
  sisl_solver_performance_t perf ;
  sisl_complex_t rc ;
  gchar *ipfile, *opfile, *matfile, *datfile, solver_name[256] ;
  FILE *input, *output ;
  GtsFile *fp ;
  gchar *progname, ch, p[32] ;
  gint np, i, j, mstride, solver ;
  guint imin, imax, itmp0, itmp1 ;
  gsl_complex zc ;
  BEM3DConfiguration *config ;
  gdouble tol ;

  wmpi_initialize(&argc, &argv) ;
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  sprintf(p, "P%03d: ", wmpi_rank()) ;
  bem3d_logging_init(stderr, p, G_LOG_LEVEL_MESSAGE, wmpi_shutdown) ;
  bem3d_shapefunc_lookup_init() ;
  sisl_logging_init(stderr, p, G_LOG_LEVEL_MESSAGE, wmpi_shutdown) ;

  wmpi_log_status_set(wmpi_rank(), TRUE) ;

  ipfile = opfile = matfile = datfile = NULL ;
  tol = 1e-3 ;

  meshes = g_ptr_array_new() ;

  bem3d_configuration_init() ;
  config = bem3d_configuration_new() ;

  while ( (ch = getopt(argc, argv, "hC:d:i:k:m:t:o:")) != EOF ) {
    switch (ch) {
    default:
    case 'h':
      if ( wmpi_rank() == 0 ) {
	fprintf(stderr, 
		"%s: solve assembled BEM3D problem\n\n",
		progname) ;
	fprintf(stderr, "Usage: %s <options>\n", progname) ;
	fprintf(stderr, 
		"Options:\n"
		"        -h (print this message and exit)\n"
		"        -C <configuration file name>\n"
		"        -d <data file name> (for boundary conditions)\n"
		"        -k # (wave number for Helmholtz calculation\n"
		"        -m <matrix file name from bem3d-assemble>\n"
		"        -t # (iterative solver convergence tolerance: %lg)\n"
		"        -o <output file name> (for mesh block data)\n",
		tol) ;
      }
      wmpi_pause() ;
      wmpi_shutdown() ;
      return 0 ;
      break ;
    case 'C': bem3d_configuration_read(config, optarg) ; break ;
    case 'd': datfile = g_strdup(optarg) ; break ;
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
    case 'k': bem3d_parameters_wavenumber(&param) = atof(optarg) ;
      break ;
    case 'm': 
      if ( wmpi_process_number() == 1 ) {
	matfile = g_strdup(optarg) ; 
      } else {
	matfile = g_strdup_printf("%s-%04d", optarg, wmpi_rank()) ;
      }
      break ;
    case 't': tol = atof(optarg) ; break ;
    case 'o': opfile = g_strdup(optarg) ; break ;
    }
  }

  if ( opfile == NULL ) opfile = g_strdup("-") ;
  if ( datfile == NULL ) datfile = g_strdup("-") ;
  if ( matfile == NULL ) matfile = g_strdup("-") ;

  if ( wmpi_rank() == 0 )
    fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;
  
  t = g_timer_new() ; g_timer_start(t) ;

  input = file_open(matfile, "-", "r", stdin) ;

  if ( wmpi_rank() == 0 )
    fprintf(stderr, "%s: reading matrices: t=%f\n",
	    progname, g_timer_elapsed(t, NULL)) ;

  fscanf(input, "%[^ ]s", solver_name) ;
  fscanf(input, "%*u %d %d %u %u %d\n", 
	 &np, &mstride, &itmp0, &itmp1, &solver) ;

  if ( wmpi_rank() == 0 )
    fprintf(stderr, "%s: solver=%s, size=%d; t=%f\n",
	    progname, solver_name, np, g_timer_elapsed(t, NULL)) ;

  config->solver = bem3d_solver_type(solver_name) ;

  if ( config->solver == BEM3D_SOLVER_DIRECT ) {
    wmpi_split_range(0, np, &imin, &imax) ;
    if ( itmp0 != imin || itmp1 != imax ) {
      g_error("matrix split error, requested rows %u--%u not in input (%u--%u)",
	      imin, imax, itmp0, itmp1) ;
    }

    if ( mstride == 1 ) rc = SISL_REAL ;
    if ( mstride == 2 ) rc = SISL_COMPLEX ;

    A = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
    B = sisl_matrix_new(rc, SISL_MATRIX_DENSE) ;
    fprintf(stderr, "%s: P%d: matrices allocated %d rows\n", 
	    progname, wmpi_rank(), imax-imin) ;
    
    phi = sisl_vector_new(rc) ;
    dphi = sisl_vector_new(rc) ;
    rhs = sisl_vector_new(rc) ;

    sisl_matrix_set_block_size(A, imax-imin, np) ; 
    sisl_matrix_set_block_size(B, imax-imin, np) ;
    sisl_matrix_row_number(A) = sisl_matrix_row_number(B) = np ;
    sisl_matrix_column_number(A) = sisl_matrix_column_number(B) = np ;
    sisl_matrix_local_row_start(A) = imin ; 
    sisl_matrix_local_row_end(A) = imax ;
    sisl_matrix_local_row_start(B) = imin ; 
    sisl_matrix_local_row_end(B) = imax ;

    if ( wmpi_rank() == 0 )
      fprintf(stderr, "%s: matrix chunks split\n", progname) ;

    i = sisl_matrix_local_row_start(A) ;
    j = sisl_matrix_local_row_end(A) ;
    fprintf(stderr, "%s: P%d: local rows: %d--%d\n", 
	    progname, wmpi_rank(), i, j) ;
  
    if ( rc == SISL_REAL ) 
      read_real_matrices(input, imax-imin, A, B) ;
    else
      read_complex_matrices(input, imax-imin, A, B) ;
  } else {
    /*set up the matrices for fast multipole solver*/
    BEM3DMeshSkeleton *skel ;
    BEM3DQuadratureRule *quad ;
    BEM3DAverage anorm = BEM3D_AVERAGE_MWAAT; /* BEM3D_AVERAGE_MWA ; */
    BEM3DFMMMatrix *mtx ;
    gpointer Adata[4] ;
    gint order ;
    gdouble r_correct ;

    g_assert(config->fmm == solver) ;
    g_assert(meshes->len > 0) ;

    imin = itmp0 ; imax = itmp1 ;

    if ( mstride == 1 ) rc = SISL_REAL ;
    if ( mstride == 2 ) rc = SISL_COMPLEX ;

    order = config->skel_order ;
    r_correct = config->fmm_radius ;

    A = sisl_matrix_new(rc, SISL_MATRIX_USER_DEFINED) ;
    B = sisl_matrix_new(rc, SISL_MATRIX_USER_DEFINED) ;
    fprintf(stderr, "%s: P%d: matrices allocated\n", 
	    progname, wmpi_rank()) ;
    
    phi = sisl_vector_new(rc) ;
    dphi = sisl_vector_new(rc) ;
    rhs = sisl_vector_new(rc) ;

    m = g_ptr_array_index(meshes, 0) ;

    quad = bem3d_quadrature_rule_new(order, 1) ;
    if ( order < 7 ) 
      bem3d_quadrature_rule_gauss(NULL, bem3d_mesh_element_sample(m), quad, 
				NULL, NULL, &order) ;
    else
      bem3d_quadrature_rule_wx(NULL, bem3d_mesh_element_sample(m), quad, 
			       NULL, NULL, &order) ;    

    if ( wmpi_rank() == 0 ) 
      fprintf(stderr, "%s: generating mesh skeleton: t=%f\n",
	      progname, g_timer_elapsed(t, NULL)) ;

    skel = bem3d_mesh_skeleton_new(m, order) ;
    bem3d_mesh_skeleton_init(skel, quad, anorm) ;

    if ( wmpi_rank() == 0 ) 
      fprintf(stderr, "%s: initializing corrected FMM matrix: t=%f\n",
	      progname, g_timer_elapsed(t, NULL)) ;

    mtx = bem3d_fmm_matrix_new(solver, 
			       (rc == SISL_REAL ? 
				BEM3D_FMM_LAPLACE : BEM3D_FMM_HELMHOLTZ),
			       skel, config, &param, r_correct) ;
    mtx->tol = config->fmm_tol ;

    for ( i = 0 ; i < np ; i ++ ) fscanf(input, "%*d %lg", &(mtx->C[i])) ;
    
    Adata[0] = mtx ;
    Adata[1] = bem3d_fmm_workspace_alloc(config->fmm, skel) ;
    Adata[2] = &(param) ;
    Adata[3] = (gdouble *)g_malloc(skel->ns*2*sizeof(gdouble)) ;

    sisl_matrix_row_number(A) = np ;
    sisl_matrix_column_number(A) = np ;
    sisl_matrix_local_row_start(A) = imin ;
    sisl_matrix_local_row_end(A) = imax ;
    sisl_matrix_user_defined_multiply(A) = sisl_matrix_vector_mul_A ;
    sisl_matrix_user_defined_data(A) = Adata ;
    sisl_matrix_distribution(A) = SISL_SINGLE ;

    sisl_matrix_row_number(B) = np ;
    sisl_matrix_column_number(B) = np ;
    sisl_matrix_local_row_start(B) = imin ;
    sisl_matrix_local_row_end(B) = imax ;
    sisl_matrix_user_defined_multiply(B) = sisl_matrix_vector_mul_B ;
    sisl_matrix_user_defined_data(B) = Adata ;
    sisl_matrix_distribution(B) = SISL_SINGLE ;
  }

  file_close(input) ;

    if ( wmpi_rank() == 0 ) 
      fprintf(stderr, "%s: reading boundary conditions: t=%f\n",
	      progname, g_timer_elapsed(t, NULL)) ;

  input = file_open(datfile, "-", "r", stdin) ;

  bem3d_mesh_data_read(&data, input, 0) ;

  file_close(input) ;

  if ( bem3d_mesh_data_node_number(data) < np ) {
    fprintf(stderr, "%s: not enough data points (%d) for %dx%d matrix\n",
	    progname, bem3d_mesh_data_node_number(data), np, np) ;
    exit(1) ;
  }

  wmpi_pause() ;
  sisl_vector_set_length(phi, np) ; 
  sisl_vector_set_length(dphi, np) ;

  if ( mstride == 1 ) set_real_vectors(phi, dphi, data) ;
  if ( mstride == 2 ) set_complex_vectors(phi, dphi, data) ;
  
  if ( sisl_is_real(phi) ) sisl_vector_set_all(phi, 0.0) ;
  else sisl_vector_set_all_complex(phi, GSL_COMPLEX_ZERO) ;

  if ( wmpi_process_number() > 1 )
    sisl_matrix_distribution(A) = sisl_matrix_distribution(B) = 
      SISL_DISTRIBUTED ;

  if ( wmpi_rank() == 0 ) 
    fprintf(stderr, "%s: starting solution: t=%f\n",
	    progname, g_timer_elapsed(t, NULL)) ;

  w = sisl_solver_workspace_new() ;
  sisl_matrix_vector_mul(B, dphi, rhs) ;

#ifdef FMM_MATRIX_CHECK
  sisl_matrix_vector_mul(A, rhs, dphi) ;

  i = 0 ; 
  zc = sisl_vector_get_complex(dphi,i) ;      
  fprintf(stderr, "%lg %lg\n", GSL_REAL(zc), GSL_IMAG(zc)) ;

#endif /*FMM_MATRIX_CHECK*/

  sisl_solve(SISL_SOLVER_BICGSTAB, A, phi, rhs, tol, 128, w, &perf) ;
  wmpi_pause() ;
  if ( wmpi_rank() == 0 ) 
    fprintf(stderr, "%s: system solved: t=%f\n",
	    progname, g_timer_elapsed(t, NULL)) ;

  sisl_matrix_free(A) ; sisl_matrix_free(B) ;

  if ( wmpi_rank() == 0 ) {
    output = file_open(opfile, "-", "w", stdout) ;
    if ( !sisl_is_real(phi) ) {
      fprintf(output, "%d 4 BEM3DMeshData\n", sisl_vector_length(phi)) ;
      for ( i = 0 ; i < sisl_vector_length(phi) ; i ++ ) {
	zc = sisl_vector_get_complex(phi,i) ;      
	fprintf(output, "%d %1.16e %1.16e", i, GSL_REAL(zc), GSL_IMAG(zc)) ;
	zc = sisl_vector_get_complex(dphi,i) ;      
	fprintf(output, " %1.16e %1.16e\n", GSL_REAL(zc), GSL_IMAG(zc)) ;
      }
    } else {
      fprintf(output, "%d 2 BEM3DMeshData\n", sisl_vector_length(phi)) ;
      for ( i = 0 ; i < sisl_vector_length(phi) ; i ++ ) 
	fprintf(output, "%d %1.16e %1.16e\n", i,
		sisl_vector_get(phi,i),
		sisl_vector_get(dphi,i)) ;
    }
    file_close(output) ;
  }

  wmpi_pause() ;

  sisl_solver_workspace_free(w) ;

  wmpi_shutdown() ;

  return 0 ;
}
