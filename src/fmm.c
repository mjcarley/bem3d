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
  default: g_error("unrecognized solver type %d", solver) ; break ;
  case BEM3D_FMM_FMMLIB3D_1_2:
    switch (problem) {
    default: g_error("unrecognized problem type %d", problem) ; break ;
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
  }

  return BEM3D_SUCCESS ;
}

BEM3DFMMWorkspace *bem3d_fmm_workspace_alloc(BEM3DFastMultipole solver,
					     BEM3DMeshSkeleton *skel)

{
  BEM3DFMMWorkspace *w ;
  gint nda ;

  w = (BEM3DFMMWorkspace *)g_malloc(sizeof(BEM3DFMMWorkspace)) ;
  
  w->solver = solver ;

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
  }

  w->nda = nda ;
  w->d = (gdouble *)g_malloc(nda*sizeof(gdouble)) ;

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
  R = sqrt(R2) ;

  g[0] = -0.25*M_1_PI/R ; 
  dg[0] = -0.25*M_1_PI*(r[0]*n[0] + r[1]*n[1] + r[2]*n[2])/R2/R ;

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
    bem3d_element_assemble_equations(e, GTS_POINT(v), config, param, G, dG) ;

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
 * @param r radius for elements to included in correction calculation.
 * 
 * @return BEM3D_SUCCESS on success.
 */

BEM3DFMMMatrix *bem3d_fmm_matrix_new(BEM3DFastMultipole solver,
				     BEM3DFastMultipoleProblem problem,
				     BEM3DMeshSkeleton *skel,
				     BEM3DConfiguration *config,
				     BEM3DParameters *param,
				     gdouble r)

{
  BEM3DFMMMatrix *m ;
  gpointer data[8] ;
  GArray *G, *dG ;

  m = (BEM3DFMMMatrix *)g_malloc(sizeof(BEM3DFMMMatrix)) ;

  m->skel = skel ;
  m->solver = solver ;
  m->problem = problem ;

  m->C = (gdouble *)g_malloc((skel->nt)*sizeof(gdouble)) ;

  m->gcorr = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  m->dgcorr = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  m->icorr = g_array_new(TRUE, TRUE, sizeof(gint)) ;
  m->idxcorr = (gint *)g_malloc(2*skel->nt*sizeof(gint)) ;

  data[0] = m ; data[1] = &r ; data[2] = config ; data[3] = param ;

  data[6] =  G = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  data[7] = dG = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;

  bem3d_mesh_foreach_node(skel->m, (BEM3DNodeFunc)correction_terms, data) ;

  g_array_free(G, TRUE) ; g_array_free(dG, TRUE) ;

  return m ;
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

/**
 * @}
 * 
 */
