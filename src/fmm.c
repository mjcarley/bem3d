/* fmm.c
 * 
 * Copyright (C) 2017, 2018 Michael Carley
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
#include <string.h>
#include <math.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

#ifdef HAVE_FMMLIB3D_1_2
#include "fmmlib3d_1_2.h"
#endif

#ifdef HAVE_WBFMM
#include "wbfmm-bem3d.h"
#endif

/**
 * @defgroup fmm Fast Multipole Method
 *
 * The Fast Multipole Method (FMM) is a standard technique for
 * accelerating BEM calculations. BEM3D includes an interface to FMM
 * solvers developed by other groups, identified by ::BEM3DFastMultipole
 *
 * @{
 * 
 */

/** 
 * Perform a Fast Multipole Method (FMM) calculation using a specified
 * solver. The source and dipole strengths, and the output arrays, may
 * be NULL, and are then ignored in the calculation, for example
 * depending on whether the single- or double-layer potential is being
 * computed.
 * 
 * @param solver a ::BEM3DFastMultipole solver (depending on what has
 * been compiled in to the library);
 * @param problem type of problem to be solved, for example,
 * ::BEM3D_FMM_LAPLACE, not all \a solver types can solve all problem types;
 * @param param parameters for the problem;
 * @param s an initialized ::BEM3DMeshSkeleton containing the point sources
 * to be used in the FMM evaluation;
 * @param tol FMM evaluation tolerance;
 * @param q packed array of source strengths for the points in \a s;
 * @param dq packed array of dipole strengths for the points in \a s;
 * @param p on output, potential evaluated at target points in \a s;
 * @param dp on output, normal derivative of potential evaluated at 
 * target points in \a s;
 * @param w a ::BEM3DFMMWorkspace allocated using
 * ::bem3d_fmm_workspace_alloc.
 * 
 * @return 0 on success.
 */

gint bem3d_fmm_calculate(BEM3DFastMultipole solver,
			 BEM3DFastMultipoleProblem problem,
			 BEM3DParameters *param,
			 BEM3DMeshSkeleton *s,
			 gdouble tol,
			 gdouble *q, gdouble *dq,
			 gdouble *p, gdouble *dp,
			 BEM3DFMMWorkspace *w)

{
  switch ( solver ) {
  default: g_error("%s: unrecognized solver type %d",
		   __FUNCTION__, solver) ;
    break ;
  case BEM3D_FMM_FMMLIB3D_1_2:
    switch ( problem ) {
    default: g_error("%s: unrecognized problem type %d for solver %d",
		     __FUNCTION__, problem, solver) ;
      break ;
    case BEM3D_FMM_LAPLACE: 
      _bem3d_fmm_laplace_fmmlib3d_1_2(solver, problem, param, s, tol,
				      q, dq, p, dp, w) ;
      break ;
    case BEM3D_FMM_HELMHOLTZ:
      _bem3d_fmm_helmholtz_fmmlib3d_1_2(solver, problem, param, s, tol,
					q, dq, p, dp, w) ;
      break ;      
    }
    break ;
  case BEM3D_FMM_WBFMM:
    switch (problem) {
    default: g_error("%s: unrecognized problem type %d for solver %d",
		     __FUNCTION__, problem, solver) ; break ;
    case BEM3D_FMM_LAPLACE: 
      _bem3d_fmm_laplace_wbfmm(solver, problem, param, s, tol,
			       q, dq, p, dp, w) ;
      break ;
    case BEM3D_FMM_HELMHOLTZ:
      _bem3d_fmm_helmholtz_wbfmm(solver, problem, param, s, tol,
				 q, dq, p, dp, w) ;
      break ;      
    }
    break ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Allocate a ::BEM3DFMMWorkspace for use in Fast Multipole Method
 * evaluations.
 * 
 * @param solver a ::BEM3DFastMultipole specifying the FMM method to be 
 * used;
 * @param skel the ::BEM3DMeshSkeleton which is to be passed to the 
 * solver;
 * 
 * @return the newly allocated ::BEM3DFMMWorkspace.
 */

BEM3DFMMWorkspace *bem3d_fmm_workspace_alloc(BEM3DFastMultipole solver,
					     BEM3DMeshSkeleton *skel,
					     BEM3DConfiguration *config)

{
  BEM3DFMMWorkspace *w ;
  gint nda ;

  w = (BEM3DFMMWorkspace *)g_malloc0(sizeof(BEM3DFMMWorkspace)) ;
  
  w->solver = solver ; nda = 0 ;
  
  switch ( solver ) {
  default: g_error("unrecognized solver type %d", solver) ; break ;
  case BEM3D_FMM_FMMLIB3D_1_2:
    nda = 
      skel->ns + /*number of charges*/
      skel->ns + /*number of dipole strengths*/
      skel->nt + /*number of field points for potential calculation*/  
      3*(skel->nt) ; /*points for field gradient calculation*/
    nda *= 2 ;   /*FMMLIB3D requires that everything be assumed complex*/
    break ;
  case BEM3D_FMM_WBFMM:
    w->i[0] = config->fmm_tree_depth ;
    /* g_assert_not_reached() ; /\*don't know what to do yet*\/ */
    break ;
  }

  w->nda = nda ;
  if ( w->nda != 0 ) 
    w->d = (gdouble *)g_malloc0(nda*sizeof(gdouble)) ;

  return w ;
}

static gint near_vertex(BEM3DElement *e, gpointer data[])

{
  gdouble r = *((gdouble *)data[1]) ;
  GSList **elements = data[4] ;
  GtsVertex *v = data[5] ;
  gdouble R ;
  gint i ;

  bem3d_element_nearest_vertex(e, GTS_POINT(v), &i, &R) ;

  if ( R < r ) *elements = g_slist_prepend(*elements, e) ;

  return 0 ;
}

static gint insert_correction_weights(GArray *gcorr, GArray *dgcorr,
				      GArray *icorr, gint i0, gint nc,
				      gint idx, gdouble *G, gdouble *dG)

{
  gint i, j ;

  g_assert(i0 <= icorr->len) ;

  for ( i = i0 ; 
	(i < icorr->len) && (idx != g_array_index(icorr,gint,i)) ; 
	i ++ ) ;    

  if ( i == icorr->len ) {
    g_array_append_val(icorr, idx) ;
    g_array_append_vals( gcorr,  G, nc) ;
    g_array_append_vals(dgcorr, dG, nc) ;

#ifdef TRACE_CORRECTION
    fprintf(stderr, "%lg %lg\n", G[0], dG[0]) ;
#endif /*TRACE_CORRECTION*/

    return 0 ;
  }

#ifdef TRACE_CORRECTION
  fprintf(stderr, "%lg %lg\n", G[0], dG[0]) ;
#endif /*TRACE_CORRECTION*/

  for ( j = 0 ; j < nc ; j ++ ) {
    g_array_index( gcorr, gdouble, nc*i+j) +=  G[j] ;
    g_array_index(dgcorr, gdouble, nc*i+j) += dG[j] ;
  }

  return 0 ;
}

static gint correction_gfunc(gdouble *x, gdouble *y, gdouble *n,
			     BEM3DParameters *p, gdouble *g, gdouble *dg,
			     BEM3DFastMultipoleProblem problem)

{
  gdouble R, R2, r[3], C, S, k ;

  r[0] = x[0] - y[0] ; r[1] = x[1] - y[1] ; r[2] = x[2] - y[2] ; 

  R2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] ;
  if ( R2 > 1e-12 ) {
    R = sqrt(R2) ;

    g[0] = -0.25*M_1_PI/R ; 
    dg[0] = -0.25*M_1_PI*(r[0]*n[0] + r[1]*n[1] + r[2]*n[2])/R2/R ;
  } else {
    R = 0.0 ; g[0] = dg[0] = 0.0 ;
  }

  if ( problem == BEM3D_FMM_LAPLACE ) return 0 ;

  k = bem3d_parameters_wavenumber(p) ;
  C = cos(k*R) ; S = sin(k*R) ;
  
  g[1] = g[0]*S ; g[0] *= C ;

  dg[1] = -dg[0]*(k*R*C - S) ; dg[0] *= (k*R*S + C) ;

  return 0 ;
}

static gint subtract_correction(gint *indices, gdouble *w, gint ppe,
				gdouble *g, gdouble *dg, 
				GArray *gcorr, GArray *dgcorr, GArray *icorr,
				gint i0, gint nc)

{
  gint i, j, k ;

  g_assert(i0 < icorr->len) ;

  for ( j = 0 ; j < ppe ; j ++ ) {
    for ( i = i0 ; 
	  (i < icorr->len) && (indices[j] != g_array_index(icorr,gint,i)) ; 
	  i ++ ) ;    

    g_assert(i < icorr->len) ;

#ifdef TRACE_CORRECTION
    fprintf(stderr, "%lg %lg\n", g[0]*w[j], dg[0]*w[j]) ;
#endif /*TRACE_CORRECTION*/

    for ( k = 0 ; k < nc ; k ++ ) {
      g_array_index( gcorr, gdouble, nc*i+k) +=  g[k]*w[j] ;
      g_array_index(dgcorr, gdouble, nc*i+k) += dg[k]*w[j] ;
    }
  }

  return 0 ;
}

static gint correction_terms(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DFMMMatrix *m = data[0] ;
  BEM3DConfiguration *config = data[2] ;
  BEM3DParameters *param = data[3] ;
  GArray *G  = data[6] ;
  GArray *dG = data[7] ;
  BEM3DWorkspace *work = data[8] ;
  GSList *elements = NULL ;
  BEM3DElement *e ;
  gdouble *w, *y, *n, g[8], dg[8] ;
  gint j, k, nc = 1, idx, *indices ;

  if ( m->problem == BEM3D_FMM_HELMHOLTZ ) nc = 2 ;

  data[4] = &elements ; data[5] = v ;
  
  bem3d_mesh_foreach_element(m->skel->m, (BEM3DElementFunc)near_vertex, data) ;

  g_assert(m->gcorr->len/nc == m->icorr->len) ;

  m->idxcorr[2*i+0] = m->icorr->len ;
  
  for ( ; elements != NULL ; elements = elements->next ) {
    e = BEM3D_ELEMENT(elements->data) ;

    /*add the integrated element contributions*/
    bem3d_element_assemble_equations(e, GTS_POINT(v), config, param,
				     G, dG, work) ;

    for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
      insert_correction_weights(m->gcorr, m->dgcorr, m->icorr,
				m->idxcorr[2*i+0], nc,
				bem3d_element_global_index(e,j),
				&(g_array_index(G,gdouble,nc*j)),
				&(g_array_index(dG,gdouble,nc*j))) ;
    }

/* #ifdef DIRECT_TEST */
    /*subtract the point source terms*/
    /*source locations and interpolation weights*/
    idx = GPOINTER_TO_INT(g_hash_table_lookup(m->skel->e, e)) ;
    y       = &(m->skel->x[3*idx]) ;
    n       = &(m->skel->n[3*idx]) ;
    indices = &(m->skel->idx[(m->skel->ppe)*idx]) ;
    w       = &(m->skel->w[(m->skel->ppe)*idx]) ;

    for ( k = 0 ; k < m->skel->order ; k ++ ) {
      correction_gfunc(&(GTS_POINT(v)->x), &(y[3*k]), &(n[3*k]), 
		       param, g, dg, m->problem) ;
      subtract_correction(indices, &(w[(m->skel->ppe)*k]), m->skel->ppe, 
			  g, dg, 
			  m->gcorr, m->dgcorr, m->icorr,
			  m->idxcorr[2*i+0], nc) ;
    }
/* #endif /\*DIRECT_TEST*\/ */
  }
  
  m->idxcorr[2*i+1] = m->icorr->len ;

#ifdef TRACE_CORRECTION
  exit(0) ;
#endif /*TRACE_CORRECTION*/

  return 0 ;
}

/** 
 * Generate a new ::BEM3DFMMMatrix which can be used in matrix-free
 * calculations. The underlying calculation method uses a set of point
 * sources to represent the mesh, in \a skel. This is then `corrected'
 * by subtracting sources which lie within distance \a r of a
 * collocation point and adding the integrated contribution from those
 * elements calculated using the quadrature rules set in the
 * configuration \a config.
 * 
 * @param solver the FMM solver to be used;
 * @param problem problem to be solved by FMM solver (e.g. Helmholtz or 
 * Laplace) 
 * @param skel a ::BEM3DMeshSkeleton for the source points which will be 
 * passed to the FMM solver;
 * @param config configuration for the problem (including quadrature rules
 * for correction calculation);
 * @param param parameters for the physics;
 * @param r radius for elements to included in correction calculation;
 * @param work a ::BEM3DWorkspace.
 * 
 * @return BEM3D_SUCCESS on success.
 */

BEM3DFMMMatrix *bem3d_fmm_matrix_new(BEM3DFastMultipole solver,
				     BEM3DFastMultipoleProblem problem,
				     BEM3DMeshSkeleton *skel,
				     BEM3DConfiguration *config,
				     BEM3DParameters *param,
				     gdouble r,
				     BEM3DWorkspace *work)

{
  BEM3DFMMMatrix *m ;
  gpointer data[16] = {NULL} ;
  GArray *G, *dG ;
  gint ns ;
  
  m = (BEM3DFMMMatrix *)g_malloc0(sizeof(BEM3DFMMMatrix)) ;

  m->skel = skel ;
  m->solver = solver ;
  m->problem = problem ;
  ns = skel->ns ;
  
  /*set the scaling factor for the output from different FMM solvers*/
  switch ( solver ) {
  default:
    g_error("%s: solver %d not recognized", __FUNCTION__, solver) ;
    break ;
  case BEM3D_FMM_FMMLIB3D_1_2:
    m->scaleA[0] = 0.25*M_1_PI ; m->scaleA[1] = 0.0 ;
    m->scaleB[0] = 0.25*M_1_PI ; m->scaleB[1] = 0.0 ;
    m->w = (gdouble *)g_malloc0(2*2*ns*sizeof(gdouble)) ;
    break ;
  case BEM3D_FMM_WBFMM:
    switch ( problem ) {
    default: g_assert_not_reached() ; break ;
    case BEM3D_FMM_LAPLACE:
      m->scaleA[0] = 1.0 ; m->scaleA[1] = 0.0 ;
      m->scaleB[0] = 1.0 ; m->scaleB[1] = 0.0 ;
      m->w = (gdouble *)g_malloc0(2*1*ns*sizeof(gdouble)) ;
      break ;
    case BEM3D_FMM_HELMHOLTZ:
      m->scaleA[0] = 0 ; m->scaleA[1] = -bem3d_parameters_wavenumber(param) ;
      m->scaleB[0] = 0 ; m->scaleB[1] =  bem3d_parameters_wavenumber(param) ;
      m->w = (gdouble *)g_malloc0(2*2*ns*sizeof(gdouble)) ;
      break ;
    }
    break ;
  }
  
  m->C = (gdouble *)g_malloc0((skel->nt)*sizeof(gdouble)) ;

  m->gcorr = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  m->dgcorr = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  m->icorr = g_array_new(TRUE, TRUE, sizeof(gint)) ;
  m->idxcorr = (gint *)g_malloc0(2*skel->nt*sizeof(gint)) ;

  data[0] = m ; data[1] = &r ; data[2] = config ; data[3] = param ;

  data[6] =  G = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  data[7] = dG = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  data[8] = work ;
  
  bem3d_mesh_foreach_node(skel->m, (BEM3DNodeFunc)correction_terms, data) ;

  g_array_free(G, TRUE) ; g_array_free(dG, TRUE) ;

  return m ;
}

gint bem3d_fmm_matrix_multiply(BEM3DFMMMatrix *m,
			       gdouble *q, gdouble *dq,
			       gdouble *p, gdouble *dp,
			       BEM3DFMMWorkspace *w)

{
  gint nc, ns, i, j, k ;
  
  g_assert(dp == NULL) ;
  g_assert(m->problem == BEM3D_FMM_LAPLACE) ; /*affects address for
						sources*/

  if ( q == NULL && dq == NULL ) return 0 ;
  
  nc = 1 ;
  ns = m->skel->ns ;
  /*generate skeleton source strengths*/
  if ( q  != NULL ) 
    bem3d_skeleton_set_sources(m->skel,  q, nc, &(m->w[ 0])) ;
  if ( dq != NULL ) 
    bem3d_skeleton_set_sources(m->skel, dq, nc, &(m->w[ns])) ;

  for ( i = 0 ; i < ns ; i ++ ) m->w[   i] *= m->scaleA[0] ;
  for ( i = 0 ; i < ns ; i ++ ) m->w[ns+i] *= m->scaleB[0] ;
  
  /*FMM calculation*/
  if ( dq == NULL ) {
    bem3d_fmm_calculate(m->solver, m->problem, NULL, m->skel, m->tol,
			NULL, m->w, p, NULL, w) ;
  
    /*add corrections*/
    for ( i = 0 ; i < m->skel->nnodes ; i ++ ) {
      for ( j = m->idxcorr[2*i+0] ; j < m->idxcorr[2*i+1] ; j ++ ) {
	k = g_array_index(m->icorr, gint, j) ;
	p[i] += g_array_index(m->dgcorr, gdouble, j)*q[k] ;
      }
    }
    return BEM3D_SUCCESS ;
  }

  if ( q == NULL ) {
    bem3d_fmm_calculate(m->solver, m->problem, NULL, m->skel, m->tol,
			&(m->w[ns]), NULL, p, NULL, w) ;
  
    /*add corrections*/
    for ( i = 0 ; i < m->skel->nnodes ; i ++ ) {
      for ( j = m->idxcorr[2*i+0] ; j < m->idxcorr[2*i+1] ; j ++ ) {
	k = g_array_index(m->icorr, gint, j) ;
	p[i] += g_array_index(m->gcorr, gdouble, j)*dq[k] ;
      }
    }
    return BEM3D_SUCCESS ;
  }
  
  bem3d_fmm_calculate(m->solver, m->problem, NULL, m->skel, m->tol,
		      &(m->w[ns]), m->w, p, NULL, w) ;
  
  /*add corrections*/
  for ( i = 0 ; i < m->skel->nnodes ; i ++ ) {
    for ( j = m->idxcorr[2*i+0] ; j < m->idxcorr[2*i+1] ; j ++ ) {
      k = g_array_index(m->icorr, gint, j) ;
      p[i] += g_array_index(m->gcorr , gdouble, j)*dq[k] ;
      p[i] += g_array_index(m->dgcorr, gdouble, j)* q[k] ;
    }
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Return a string giving a name to a ::BEM3DSolver
 * 
 * @param s a solver type
 * 
 * @return a string containing a human-readable name for \a s
 */

gchar *bem3d_solver_name(BEM3DSolver s)

{
  switch ( s ) {
  default: g_error("%s: unrecognized solver type %d", __FUNCTION__, s) ;
  case BEM3D_SOLVER_DIRECT: return "direct" ; break ;
  case BEM3D_SOLVER_FMM: return "fmm" ; break ;
 }

  return NULL ;
}

/** 
 * Inverse of ::bem3d_solver_name, given a human-readable name for a
 * solver, return the corresponding ::BEM3DSolver 
 * 
 * @param s a string containing a solver name returned by
 * ::bem3d_solver_name
 * 
 * @return the ::BEM3DSolver corresponding to \a s
 */

BEM3DSolver bem3d_solver_type(gchar *s)

{
  if ( !strcmp(s, "direct") ) return BEM3D_SOLVER_DIRECT ;
  if ( !strcmp(s, "fmm") ) return BEM3D_SOLVER_FMM ;
  
  g_error("%s: unrecognized solver type %s", __FUNCTION__, s) ;

  return 0 ;
}

gchar *bem3d_fmm_name(BEM3DFastMultipole solver)

{
  switch ( solver ) {
  default: return NULL ;
  case BEM3D_FMM_FMMLIB3D_1_2:
    return "fmmlib3d1.2" ;
    break ;
  case BEM3D_FMM_WBFMM:
    return "wbfmm" ;
    break ;
  }
  
  return NULL ;
}

/**
 * @}
 * 
 */
