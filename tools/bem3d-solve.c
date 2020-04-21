/* bem3d-solve.c
 * 
 * Copyright (C) 2006, 2018, 2019 Michael Carley
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

gchar *progname ;

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
  gdouble *wt, *wp, *vp ;

  g_assert(sisl_matrix_density(A) == SISL_MATRIX_USER_DEFINED) ;

  s = mtx->skel ;
  /*set the source strengths from the input vector*/
  if ( sisl_is_real(A) ) {
    g_assert_not_reached() ;
  } else {
    memset(scratch, 0, 2*s->ns*sizeof(gdouble)) ;
    vp = sisl_vector_data(v) ;

    for ( i = 0 ; i < s->ns ; i ++ ) {
      wt = &(s->w[i*s->ppe]) ; idx = &(s->idx[i*s->ppe]) ;

      for ( j = 0 ; j < s->ppe ; j ++ ) {
	scratch[2*i+0] +=  wt[j]*vp[2*idx[j]+0] ;
	scratch[2*i+1] +=  wt[j]*vp[2*idx[j]+1] ;
      }
    }
  }

  wp = sisl_vector_data(w) ;
  bem3d_fmm_calculate(mtx->solver, mtx->problem, param, mtx->skel,
  		      mtx->tol, NULL, scratch,
  		      wp, NULL, work) ;

  sisl_vector_scale_complex(w, *((gsl_complex *)(mtx->scaleA))) ;
  
  if ( sisl_is_real(A) ) {
    g_assert_not_reached() ;
  } else {

    gdouble *wp, *vp, *cp ;

    cp = &(g_array_index(mtx->dgcorr, gdouble, 0)) ;

    vp = sisl_vector_data(v) ; wp = sisl_vector_data(w) ;

    for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
      wp[2*i+0] -= mtx->C[i]*vp[2*i+0] ;
      wp[2*i+1] -= mtx->C[i]*vp[2*i+1] ;
    }

    for ( i = 0 ; i < sisl_vector_length(v) ; i ++ ) {
      for ( j = mtx->idxcorr[2*i+0] ; j < mtx->idxcorr[2*i+1] ; j ++ ) {
	k = g_array_index(mtx->icorr, gint, j) ;
	wp[2*i+0] += vp[2*k+0]*cp[2*j+0] - vp[2*k+1]*cp[2*j+1] ;
	wp[2*i+1] += vp[2*k+0]*cp[2*j+1] + vp[2*k+1]*cp[2*j+0] ;
      }
    }
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

    vp = sisl_vector_data(v) ;
    
    for ( i = 0 ; i < s->ns ; i ++ ) {
      wt = &(s->w[i*s->ppe]) ; idx = &(s->idx[i*s->ppe]) ;

      for ( j = 0 ; j < s->ppe ; j ++ ) {
	scratch[2*i+0] += wt[j]*vp[2*idx[j]+0] ;
	scratch[2*i+1] += wt[j]*vp[2*idx[j]+1] ;
      }
    }
  }
  wp = sisl_vector_data(w) ;
  bem3d_fmm_calculate(mtx->solver, mtx->problem, param, mtx->skel,
  		      mtx->tol,
  		      scratch, NULL, wp, NULL, work) ;
  sisl_vector_scale_complex(w, *((gsl_complex *)(mtx->scaleB))) ;
  
  if ( sisl_is_real(B) ) {
    g_assert_not_reached() ;
  } else {
    /*correction terms to FMM multiplication*/

    wp = sisl_vector_data(w) ;
    vp = sisl_vector_data(v) ;

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
  }

  return 0 ;
}

static sisl_matrix_t *read_surface_diagonal(FILE *f, gdouble *Ad, gint nc,
					    gint str, gint imin, gint imax,
					    sisl_matrix_t *Ai)

{
  sisl_matrix_t *A ;
  sisl_complex_t rc ;
  gint i, j, lineno ;
  gdouble x ;
  gsl_complex xc ;
  gchar line[1024] ;
  
  rc = (str == 1 ? SISL_REAL : SISL_COMPLEX) ;
  if ( Ai == NULL ) {
    A = sisl_matrix_new(rc, SISL_MATRIX_DIAGONAL) ;
    sisl_matrix_set_block_size(A, nc, nc) ;
    sisl_matrix_local_row_start(A) = imin ;
    sisl_matrix_local_row_end(A) = imax ;
    sisl_matrix_row_number(A) = sisl_matrix_column_number(A) = nc ;
  } else {
    A = Ai ;
  }

  lineno = 3 ;
  if ( sisl_is_real(A) ) {
    g_assert_not_reached() ; /*unchecked code*/
    sisl_matrix_set_all(A, Ad[0]) ;
    while ( (i = fscanf(f, "%[^\n]c", line)) != EOF && ( i != 0 ) ) {
      i = sscanf(line, "%d %lg", &j, &x) ;
      if ( i != 2 ) g_error("%s: cannot parse line %d\n  %s\n",
			    __FUNCTION__, lineno, line) ;
      sisl_matrix_set(A, j, j, x) ;
      lineno ++ ;
      if ( (i = fscanf(f, "%*c")) == EOF ) break ;
    }    
  } else {
    GSL_SET_COMPLEX(&xc, Ad[0], Ad[1]) ;
    sisl_matrix_set_all_complex(A, xc) ;
    while ( (i = fscanf(f, "%[^\n]c", line)) != EOF && ( i != 0 ) ) {
      i = sscanf(line, "%d %lg %lg", &j, &GSL_REAL(xc), &GSL_IMAG(xc)) ;
      if ( i != 3 ) g_error("%s: cannot parse line %d\n  %s\n",
			    __FUNCTION__, lineno, line) ;
      sisl_matrix_set_complex(A, j, j, xc) ;
      xc = sisl_matrix_get_complex(A, j, j) ;
      lineno ++ ;
      if ( (i = fscanf(f, "%*c")) == EOF ) break ;
    }
  }

  return A ;  
}

static sisl_matrix_t *read_surface_sparse(FILE *f, gdouble *Ad, gint nc,
					  gint str, gint imin, gint imax,
					  sisl_matrix_t *Ai)

{
  sisl_matrix_t *A ;
  sisl_complex_t rc ;
  gint i, j, n, lineno ;
  gsl_complex xc ;
  gchar line[1024] ;
  
  rc = (str == 1 ? SISL_REAL : SISL_COMPLEX) ;
  if ( Ai == NULL ) {
    A = sisl_matrix_new(rc, SISL_MATRIX_SPARSE) ;
    /* sisl_matrix_set_block_size(A, nc, nc) ; */
    sisl_matrix_local_row_start(A) = imin ;
    sisl_matrix_local_row_end(A) = imax ;
    sisl_matrix_row_number(A) = sisl_matrix_column_number(A) = nc ;
  } else {
    A = Ai ;
  }

  lineno = 3 ;
  if ( sisl_is_real(A) ) {
    g_assert_not_reached() ; /*unchecked code*/
    /* sisl_matrix_set_all(A, Ad[0]) ; */
    /* while ( (i = fscanf(f, "%[^\n]c", line)) != EOF && ( i != 0 ) ) { */
    /*   i = sscanf(line, "%d %lg", &j, &x) ; */
    /*   if ( i != 2 ) g_error("%s: cannot parse line %d\n  %s\n", */
    /* 			    __FUNCTION__, lineno, line) ; */
    /*   sisl_matrix_set(A, j, j, x) ; */
    /*   lineno ++ ; */
    /*   if ( (i = fscanf(f, "%*c")) == EOF ) break ; */
    /* }     */
  } else {
    /* GSL_SET_COMPLEX(&xc, Ad[0], Ad[1]) ; */
    /* sisl_matrix_set_all_complex(A, xc) ; */
    while ( (n = fscanf(f, "%[^\n]c", line)) != EOF && ( n != 0 ) ) {
      n = sscanf(line, "%d %d %lg %lg", &i, &j, &GSL_REAL(xc), &GSL_IMAG(xc)) ;
      if ( n != 4 ) g_error("%s: cannot parse line %d\n  %s\n",
			    __FUNCTION__, lineno, line) ;
      sisl_matrix_set_complex(A, i, j, xc) ;
      lineno ++ ;
      if ( (n = fscanf(f, "%*c")) == EOF ) break ;
    }
  }

  return A ;  
}

static sisl_matrix_t *read_surface_matrix(FILE *f, gdouble *Ad,
					  gint imin, gint imax,
					  gint nnodes, sisl_matrix_t *Ai)

/*
  read surface property data and return the corresponding surface
  matrix, with default diagonal entries Ad (real or complex, which is
  why it's a pointer); currently only implemented for diagonal
  matrices corresponding to locally-reacting boundary conditions, and
  for sparse non-local boundary conditions;

 returns NULL on error, which must be checked for.

  Input syntax, header:

  [size (number of nodes)] [stride=1, 2 for real or complex] BEM3DSurface 
  [matrix format]

  Currently implemented for diagonal matrix with second line (matrix format):

  diagonal [number of entries]

  For diagonal matrix not all entries need be specifed (Ad is used as
  default), and each line contains:

  [index] [entry]
  
*/
  
{
  sisl_matrix_t *A = Ai ;
  gchar line[1024], **tokens ;
  gint nc, str ;
    
  fscanf(f, "%[^\n]s", line) ;
  fscanf(f, "%*c") ;

  /*check the header is correct*/
  tokens = g_strsplit(line, " ",0) ;
  if ( tokens[0] == NULL ) return NULL ;
  if ( tokens[1] == NULL ) return NULL ;
  if ( tokens[2] == NULL ) return NULL ;

  if ( strcmp(tokens[2], "BEM3DSurface") != 0 )
    g_error("%s: cannot read surface file (incorrect header %s)",
	    __FUNCTION__, tokens[2]) ;

  nc = atoi(tokens[0]) ;
  if ( nc <= 0 )
    g_error("%s: cannot read surface file (invalid size %s)",
	    __FUNCTION__, tokens[0]) ;
  if ( (nc > nnodes) && (nnodes != -1) )
    g_error("%s: number of nodes in surface file (%d) is greater than "
	    "problem size (%d)", __FUNCTION__, nc, nnodes) ;
  nc = MAX(nc, nnodes) ;
  
  str = atoi(tokens[1]) ;
  if ( str <= 0 || str > 2 )
    g_error("%s: cannot read surface file (invalid stride %s)",
	    __FUNCTION__, tokens[1]) ;

  fscanf(f, "%[^\n]s", line) ;
  fscanf(f, "%*c") ;

  if ( strncmp(line, "diagonal", 8) == 0) {
    return read_surface_diagonal(f, Ad, nc, str, imin, imax, A) ;
  }

  if ( strncmp(line, "sparse", 6) == 0) {
    return read_surface_sparse(f, Ad, nc, str, imin, imax, A) ;
  }

  g_error("%s: unrecognized or unimplemented matrix type %s",
	  __FUNCTION__, line) ;
  
  return A ;
}

static gint parse_fmm_file(FILE *f, gchar **self_file, gchar **skel_file)

{
  gchar line[1024] ;

  *self_file = *skel_file = NULL ;
  while ( (fscanf(f, "%[^\n]c", line) != EOF) &&
	  ((*self_file == NULL) || (*skel_file == NULL)) ) {
    if ( strncmp(line, "self:", 5) == 0 ) {
      *self_file = g_strdup(&(line[6])) ;
    }
    if ( strncmp(line, "skeleton:", 9) == 0 ) {
      *skel_file = g_strdup(&(line[10])) ;
    }
    fscanf(f, "%*c") ;
  }
  
  if ( *skel_file == NULL ) {
    fprintf(stderr, "%s: could not find skeleton file\n", progname) ;
  }

  if ( *self_file == NULL ) {
    fprintf(stderr, "%s: could not find self term file\n", progname) ;
  }
  
  return 0 ;
}

static gint read_self_file(gchar *file, gdouble *C, gint n)

{
  FILE *f ;
  gint i, j ;
  gdouble tmp ;
  
  f = fopen(file, "r") ;
  if ( f == NULL ) {
    fprintf(stderr, "%s: cannot open file %s for self terms\n",
	    progname, file) ;
    exit(1) ;
  }

  fscanf(f, "self: %d\n", &i) ;
  if ( i != n ) {
    fprintf(stderr,
	    "%s: self-file (%s) header (%d) does not match problem size (%d)\n",
	    progname, file, i, n) ;
  }
  for ( i = 0 ; i < n ; i ++ ) {
    fscanf(f, "%d %lg\n", &j, &tmp) ;
    C[j] = tmp ;
  }
  
  fclose(f) ;

  return 0 ;
}

static gint read_skel_file(gchar *file, BEM3DMeshSkeleton *s)

{
  FILE *f ;
  
  f = fopen(file, "r") ;
  if ( f == NULL ) {
    fprintf(stderr, "%s: cannot open file %s for mesh skeleton\n",
	    progname, file) ;
    exit(1) ;
  }

  bem3d_mesh_skeleton_read(s, f) ;
  
  fclose(f) ;
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  GTimer *t ;
  BEM3DMesh *m ;
  BEM3DMeshData *data ;
  BEM3DParameters param ;
  BEM3DWorkspace *work ;
  GPtrArray *meshes ;
  sisl_matrix_t *A, *B, *surface_alpha, *surface_beta ;
  sisl_vector_t *phi, *dphi, *rhs, *v1, *v2 ;
  sisl_solver_workspace_t *w ;
  sisl_solver_performance_t perf ;
  sisl_complex_t rc ;
  gchar *ipfile, *opfile, *matfile, *datfile, solver_name[256] ;
  gchar *skel_file, *self_file ;
  FILE *input, *output ;
  GtsFile *fp ;
  gchar ch, p[32], *surface_alpha_file, *surface_beta_file ;
  gint np, i, j, mstride, solver ;
  guint imin, imax, itmp0, itmp1 ;
  gsl_complex zc, ac ;
  BEM3DConfiguration *config ;
  gdouble tol ;
  gboolean invertB ;
  
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
  bem3d_parameters_wavenumber(&param) = G_MAXDOUBLE ;

  surface_alpha_file = surface_beta_file = NULL ;
  surface_alpha = surface_beta = NULL ;
  invertB = FALSE ;
  
  while ( (ch = getopt(argc, argv, "hA:C:d:i:k:m:t:o:Z:")) != EOF ) {
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
		"        -A <surface admittance data file name>\n"
		"        -C <configuration file name>\n"
		"        -d <data file name> (for boundary conditions)\n"
		"        -i <bem3d input file> (for FMM calculations)\n"
		"        -k # (wave number for FMM Helmholtz calculation)\n"
		"        -m <matrix file name from bem3d-assemble>\n"
		"        -o <output file name> (for mesh block data)\n"
		"        -t # (iterative solver convergence tolerance: %lg)\n",
		/* "        -Z <surface admittance data file name>\n", */
		tol) ;
      }
      wmpi_pause() ;
      wmpi_shutdown() ;
      return 0 ;
      break ;
    case 'A': surface_alpha_file = g_strdup(optarg) ; break ;
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
    case 'Z':
      g_assert_not_reached() ; /*unchecked code*/
      surface_beta_file = g_strdup(optarg) ; break ;
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

  bem3d_parameters_quadrature_tol(&param) = config->quad_tol ;  
  
  work = bem3d_workspace_new() ;
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
    BEM3DFMMMatrix *mtx ;
    gpointer Adata[4] ;
    gint order ;
    gdouble r_correct ;

    g_assert(config->fmm == solver) ;
    g_assert(meshes->len > 0) ;

    imin = itmp0 ; imax = itmp1 ;

    if ( mstride == 1 ) rc = SISL_REAL ;
    if ( mstride == 2 ) rc = SISL_COMPLEX ;

    if ( rc == SISL_COMPLEX &&
	 bem3d_parameters_wavenumber(&param) == G_MAXDOUBLE ) {
      fprintf(stderr,
	      "%s: P%d: wavenumber not set for complex FMM calculation\n",
	      progname, wmpi_rank()) ;
      return 1 ;
    }
    
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

    parse_fmm_file(input, &self_file, &skel_file) ;

    /* fprintf(stderr, "self: %s\n", self_file) ; */
    /* fprintf(stderr, "skel: %s\n", skel_file) ; */
    
    skel = bem3d_mesh_skeleton_new(m, order) ;
    read_skel_file(skel_file, skel) ;
			     
    if ( wmpi_rank() == 0 ) 
      fprintf(stderr, "%s: initializing corrected FMM matrix: t=%f\n",
	      progname, g_timer_elapsed(t, NULL)) ;

    mtx = bem3d_fmm_matrix_new(solver, 
			       (rc == SISL_REAL ? 
				BEM3D_FMM_LAPLACE : BEM3D_FMM_HELMHOLTZ),
			       skel, config, &param, r_correct, work) ;
    mtx->tol = config->fmm_tol ;

    read_self_file(self_file, mtx->C, np) ;
    
    Adata[0] = mtx ;
    Adata[1] = bem3d_fmm_workspace_alloc(config->fmm, skel, config) ;
    Adata[2] = &(param) ;
    Adata[3] = (gdouble *)g_malloc0(skel->ns*2*sizeof(gdouble)) ;

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

  /*surface treatment (impedance) matrices*/
  /*limited to diagonal (locally-reacting) case for now*/
  if ( surface_beta_file == NULL ) {
    /*default setting for impedance term*/
    surface_beta = sisl_matrix_new(rc, SISL_MATRIX_DIAGONAL) ;    
    sisl_matrix_set_block_size(surface_beta, np, np) ;
    sisl_matrix_row_number(surface_beta) =
    sisl_matrix_column_number(surface_beta) = np ;
    sisl_matrix_local_row_start(surface_beta) = imin ;
    sisl_matrix_local_row_end(surface_beta) = imax ;

    /*this is the default for a locally-reacting surface*/
    sisl_matrix_set_all(surface_beta, 1.0) ;
  } else {
    g_assert_not_reached() ;
  }
  
  if ( surface_alpha_file == NULL ) {
    surface_alpha = sisl_matrix_new(rc, SISL_MATRIX_DIAGONAL) ;
    sisl_matrix_set_block_size(surface_alpha, np, np) ; 
    sisl_matrix_row_number(surface_alpha) =
      sisl_matrix_column_number(surface_alpha) = np ;
    sisl_matrix_local_row_start(surface_alpha) = imin ;
    sisl_matrix_local_row_end(surface_alpha) = imax ;

    /*take the local admittance from the configuration (usually zero)*/
    if ( rc == SISL_REAL ) {
      sisl_matrix_set_all(surface_alpha, config->bc_default_admittance[0]) ;
    } else {
      GSL_SET_COMPLEX(&zc,
		      config->bc_default_admittance[0],
		      config->bc_default_admittance[1]) ;
      sisl_matrix_set_all_complex(surface_alpha, zc) ;
    }
  } else {
    /*we need to read surface admittance data*/
    fprintf(stderr, "%s: reading surface admittance matrix %s: t=%f\n",
	    progname, surface_alpha_file, g_timer_elapsed(t, NULL)) ;
    GSL_SET_COMPLEX(&zc,
		    config->bc_default_admittance[0],
		    config->bc_default_admittance[1]) ;
    input = file_open(surface_alpha_file, "", "r", NULL) ;
    surface_alpha = read_surface_matrix(input, &(GSL_REAL(zc)), imin, imax,
				    np, surface_alpha) ;
    if ( surface_alpha == NULL ) {
      fprintf(stderr, "%s: could not read surface data matrix from %s",
	      progname, surface_alpha_file) ;
      exit(1) ;
    }
    
    file_close(input) ;
  }
  
  /*if required, invert the surface admittance matrix*/
  if ( invertB ) sisl_matrix_invert(surface_beta) ;

  if ( config->solver == BEM3D_SOLVER_DIRECT ) {
    /*apply the surface treatment to the A matrix: A-> A - B\beta^{-1}\alpha*/
    if ( rc == SISL_REAL ) {
      sisl_matrix_triple_mul_w(B, surface_beta, surface_alpha, A, -1.0, 1.0) ;
    } else {
      GSL_SET_COMPLEX(&zc, 1.0, 0.0) ; GSL_SET_COMPLEX(&ac, -1.0, 0.0) ;
      sisl_matrix_triple_mul_w_complex(B, surface_beta, surface_alpha,
				       A, ac, zc) ;
    }
  } else {
    fprintf(stderr, "%s: FMM solver, no surface conditions applied\n",
	    progname) ;
  }

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
  
  if ( rc == SISL_REAL ) set_real_vectors(phi, dphi, data) ;
  if ( rc == SISL_COMPLEX ) set_complex_vectors(phi, dphi, data) ;
  
  if ( wmpi_process_number() > 1 )
    sisl_matrix_distribution(A) = sisl_matrix_distribution(B) = 
      SISL_DISTRIBUTED ;

  if ( wmpi_rank() == 0 ) 
    fprintf(stderr, "%s: starting solution: t=%f\n",
	    progname, g_timer_elapsed(t, NULL)) ;

  v1 = sisl_vector_new(rc) ;
  v2 = sisl_vector_new(rc) ;
  sisl_vector_set_length(v1, np) ;
  sisl_vector_set_length(v2, np) ;

  /*initialize the right hand side*/
  sisl_matrix_vector_mul(surface_alpha, phi, v1) ;
  sisl_matrix_vector_mul(surface_beta, v1, v2) ;
  sisl_vector_sub(v2, dphi) ;
  sisl_matrix_vector_mul(B, v2, rhs) ;

  /*store the incident potential for later use*/
  sisl_vector_copy(v1, phi) ;
  
  w = sisl_solver_workspace_new() ;

  if ( sisl_is_real(phi) ) sisl_vector_set_all(phi, 0.0) ;
  else sisl_vector_set_all_complex(phi, GSL_COMPLEX_ZERO) ;

  sisl_solve(SISL_SOLVER_BICGSTAB, A, phi, rhs, tol, 128, w, &perf) ;

  /*phi now contains the scattered potential: generate the scattered
    dphi to make the solution complete*/
  sisl_vector_add(v1, phi) ;
  sisl_matrix_vector_mul(surface_alpha, v1, v2) ;
  sisl_matrix_vector_mul(surface_beta, v2, v1) ;
  sisl_vector_sub(v1, dphi) ;
  
  wmpi_pause() ;
  if ( wmpi_rank() == 0 ) 
    fprintf(stderr, "%s: system solved: t=%f\n",
	    progname, g_timer_elapsed(t, NULL)) ;

  sisl_matrix_free(A) ; sisl_matrix_free(B) ;

  /*output from here is the scattered potential (only)*/
  if ( wmpi_rank() == 0 ) {
    output = file_open(opfile, "-", "w", stdout) ;
    if ( !sisl_is_real(phi) ) {
      fprintf(output, "%d 4 BEM3DMeshData\n", sisl_vector_length(phi)) ;
      for ( i = 0 ; i < sisl_vector_length(phi) ; i ++ ) {
	zc = sisl_vector_get_complex(phi,i) ;      
	fprintf(output, "%d %1.16e %1.16e", i, GSL_REAL(zc), GSL_IMAG(zc)) ;
	/* zc = sisl_vector_get_complex(dphi,i) ;       */
	zc = sisl_vector_get_complex(v1,i) ;      
	fprintf(output, " %1.16e %1.16e\n", GSL_REAL(zc), GSL_IMAG(zc)) ;
      }
    } else {
      fprintf(output, "%d 2 BEM3DMeshData\n", sisl_vector_length(phi)) ;
      for ( i = 0 ; i < sisl_vector_length(phi) ; i ++ ) 
	fprintf(output, "%d %1.16e %1.16e\n", i,
		sisl_vector_get(phi,i),
		sisl_vector_get(v1,i)) ;
		/* sisl_vector_get(dphi,i)) ; */
    }
    file_close(output) ;
  }

  wmpi_pause() ;

  sisl_solver_workspace_free(w) ;

  wmpi_shutdown() ;

  return 0 ;
}
