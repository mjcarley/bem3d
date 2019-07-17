/* bem3d-aero.c
 * 
 * Copyright (C) 2006, 2009, 2010, 2018 Michael Carley
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

/** 
@page bem3daero bem3d-aero: a potential code for unsteady aerodynamics

bem3d-aero is a Morino type code for quasi-potential flow
calculations. It is intended to handle higher order elements, although
at present it only deals with linear triangular panels (although this
is still better than zero order). The formulation, and in particular
the implementation of the trailing edge conditions, is described in
more detail in Tai, A, Boundary element methods for flow over
deforming bodies, PhD thesis, University of Bath, 2009.

To perform a calculation with bem3d-aero, the steps are the usual ones
of generating a geometry, assembling the matrix, applying the boundary
conditions and solving. There are some important differences, however:
the first is that the sharp edges must be identified (as in standard
problems) but the trailing edge(s) must be specified so that-the
second major difference-the wake(s) can be shed as required. This can
be done using the special options of @ref msh2bem3d "msh2bem3d" to
generate edge files and of @ref bem3d2pos "bem3d2pos" to visualize the
edges and select those which shed a wake. 

A typical invocation of bem3d-aero (given in the examples) would be 
@verbatim
bem3d-aero -i naca0012.bem -m naca0012.mtn \
    -e naca0012.edg \
    -o naca0012-%04d.bem -O naca0012-%04d.dat \
    -w wake-%04d.bem -W wake-%04d.dat \
    -S 0 -d 0.5 -f 10 -Z 0.01 @endverbatim
which uses a geometry file @c naca0012.bem and a @ref motion 
"motion file" @c naca0012.mtn. The wing has one wake-shedding edge, described
in the @ref edge "edge file" @c naca0012.edg. The results for each time
step are output as BEM files for the surfaces of the body and the wake
which is shed from the trailing edge, and as @ref data "data blocks" for
the potential and the potential gradient on the surfaces (actually potential
jump on the wake). These data are saved in files whose name is given by
the templates @c naca0012-%04d.*, etc, where the @c %04d is a C-style format
specifier used to generate an index from the time point. Any valid C format
can be used but this is a useful one since it generates a list of files in
time-point order.

The remaining options set the start time at zero, the time step to
0.5, the final time to 10 and the trailing edge condition stand off
distance to 0.01.

In the example directory, the edge file is already set up for you, but
for your own problems you should look at the information given in 
@ref msh2bem3d "msh2bem3d". 

A handy option to remember is @c -G, which generates `graphical
output', i.e. it does not run the solver but only computes the body
positions and generates the wake surface at each time step. You can
use this to check that you have set up the problem correctly,
including selecting the right edge for wake shedding. 

**/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#include <sisl.h>

#include <wmpi.h>

#include "bem3d.h"
#include "bem3d-private.h"

#define DATA_WIDTH         16
#define DATA_MESH_1         0
#define DATA_MESH_2         1
#define DATA_MATRIX_A       2
#define DATA_DPHI           3
#define DATA_RHS            4
#define DATA_EPSILON        5
#define DATA_POINT          6
#define DATA_INDEX          7
#define DATA_WAKE_HASH      8
#define DATA_WAKE_LINES     9
#define DATA_MESH_DATA_1   10
#define DATA_MESH_DATA_2   11
#define DATA_PHI           12
#define DATA_COUNTER       13
#define DATA_CONFIG        14
#define DATA_WORK          15

static gint element_assemble(BEM3DElement *e, gpointer data[])

{
  sisl_matrix_t *A = data[DATA_MATRIX_A] ;
  GArray *dphi = data[DATA_DPHI] ;
  sisl_vector_t *rhs = data[DATA_RHS] ;
  GtsPoint *p = data[DATA_POINT] ;
  gint i = *((gint *)(data[DATA_INDEX])) ;
  BEM3DConfiguration *config = data[DATA_CONFIG] ;
  BEM3DWorkspace *work = data[DATA_WORK] ;
  static GArray *G = NULL ;
  static GArray *dGdn = NULL ;
  gint j, k ;
  BEM3DParameters pm ;

  if ( G == NULL ) {
    G = g_array_new(FALSE, FALSE, sizeof(gdouble)) ;
    dGdn = g_array_new(FALSE, FALSE, sizeof(gdouble)) ;
  }
  
  g_array_set_size(G,bem3d_element_node_number(e)) ; 
  g_array_set_size(dGdn,bem3d_element_node_number(e)) ;
  bem3d_element_assemble_equations(e, p, config, &pm, G, dGdn, work) ;
  
  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
    k = bem3d_element_global_index(e,j) ;
    sisl_matrix_set(A, i, k,
		    sisl_matrix_get(A, i, k)-
		    g_array_index(dGdn,gdouble,j)) ;
    sisl_vector_set(rhs, i, sisl_vector_get(rhs,i)-
		    g_array_index(dphi,gdouble,k)*
		    g_array_index(G,gdouble,j)) ;
  }
  
  return 0 ;
}

static gint self_assemble(BEM3DElement *e, gpointer data[])

{
  sisl_matrix_t *A = data[DATA_MATRIX_A] ;
  GArray *dphi = data[DATA_DPHI] ;
  sisl_vector_t *rhs = data[DATA_RHS] ;
  GtsPoint *p = data[DATA_POINT] ;
  gint i = *((gint *)(data[DATA_INDEX])) ;
  BEM3DConfiguration *config = data[DATA_CONFIG] ;
  BEM3DWorkspace *work = data[DATA_WORK] ;  
  static GArray *G = NULL ;
  static GArray *dGdn = NULL ;
  gint j, k ;
  BEM3DParameters pm ;

  if ( G == NULL ) {
    G = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
    dGdn = g_array_new(TRUE, TRUE, sizeof(gdouble)) ;
  }
  
  g_array_set_size(G,bem3d_element_node_number(e)) ; 
  g_array_set_size(dGdn,bem3d_element_node_number(e)) ;
  bem3d_element_assemble_equations(e, p, config, &pm, G, dGdn, work) ;

  if ( isnan(g_array_index(G,gdouble,0)) ||
       isnan(g_array_index(G,gdouble,1)) ||
       isnan(g_array_index(G,gdouble,2)) ||
       isnan(g_array_index(dGdn,gdouble,0)) ||
       isnan(g_array_index(dGdn,gdouble,1)) ||
       isnan(g_array_index(dGdn,gdouble,2)) )
    fprintf(stderr, "NaN\n") ;
  
  for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
    k = bem3d_element_global_index(e,j) ;
    sisl_matrix_set(A, i, k,
		    sisl_matrix_get(A, i, k)-
		    g_array_index(dGdn,gdouble,j)) ;
    sisl_matrix_set(A, i, i,
		    sisl_matrix_get(A, i, i) +
		    g_array_index(dGdn,gdouble,j)) ;
    sisl_vector_set(rhs, i, sisl_vector_get(rhs,i)-
		    g_array_index(dphi,gdouble,k)*
		    g_array_index(G,gdouble,j)) ;
  }

  return 0 ;
}

static gint mesh_assemble(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMesh *m1 = data[DATA_MESH_1] ;
  BEM3DMesh *m2 = data[DATA_MESH_2] ;
  sisl_matrix_t *A = data[DATA_MATRIX_A] ;
  GArray *dphi = data[DATA_DPHI] ;
  sisl_vector_t *rhs = data[DATA_RHS] ;
  gdouble ee = *((gdouble *)(data[DATA_EPSILON])) ;
  GtsPoint *p = data[DATA_POINT] ;
  GHashTable *wakenodes = data[DATA_WAKE_HASH] ;
  gboolean wake_node = FALSE ;
  GtsVector n ;
  gint j ;

  for ( j = 0 ; j < sisl_matrix_column_number(A) ; j ++ ) 
    sisl_matrix_set(A, i, j, 0.0) ;

  if ( (g_hash_table_lookup(wakenodes, GINT_TO_POINTER(i+1))) != NULL )
    wake_node = TRUE ;

  sisl_vector_set(rhs, i, 0.0) ;
  sisl_matrix_set(A, i, i, 1.0) ;

  p->x = GTS_POINT(v)->x ; 
  p->y = GTS_POINT(v)->y ;
  p->z = GTS_POINT(v)->z ;

  if ( wake_node ) {
    bem3d_node_normal(m1, i, n, BEM3D_AVERAGE_MWA) ;
    p->x += n[0]*ee ; p->y += n[1]*ee ; p->z += n[2]*ee ;
  }

  data[DATA_INDEX] = &i ;
  if ( m1 != m2 || wake_node ) 
    bem3d_mesh_foreach_element(m2, (BEM3DElementFunc)element_assemble,
			       data) ;
  else
    bem3d_mesh_foreach_element(m1, (BEM3DElementFunc)self_assemble,
			       data) ;

  if ( wake_node ) {
    sisl_vector_set(rhs, i, sisl_vector_get(rhs,i)-
		    ee*g_array_index(dphi,gdouble,i)) ;
  }

  return 0 ;
}

static void mesh_mesh_assemble(BEM3DMesh *m1, BEM3DMesh *m2, 
			       BEM3DConfiguration *config,
			       GHashTable *wakenodes,
			       sisl_matrix_t *A, gdouble ee,
			       GArray *dphi, sisl_vector_t *rhs,
			       BEM3DWorkspace *work)

{
  gpointer data[DATA_WIDTH] = {NULL} ;
  gint count ;

  data[DATA_MESH_1] = m1 ; data[DATA_MESH_2] = m2 ;
  data[DATA_MATRIX_A] = A ; data[DATA_DPHI] = dphi ; data[DATA_RHS] = rhs ;
  data[DATA_EPSILON] = &ee ;

  data[DATA_POINT] = gts_point_new(gts_point_class(), 0, 0, 0) ;
  data[DATA_WAKE_HASH] = wakenodes ;

  data[DATA_COUNTER] = &count ; count = 0 ;
  data[DATA_CONFIG] = config ;

  data[DATA_WORK] = work ;
  
  bem3d_mesh_foreach_node(m1, (BEM3DNodeFunc)mesh_assemble, data) ;

  return ;
}

static gint wake_element_assemble(BEM3DElement *e, gpointer data[])

{
  GHashTable *wakelines = data[DATA_WAKE_LINES] ;
  GtsPoint *p = data[DATA_POINT] ;
  BEM3DMeshData *dp = data[DATA_MESH_DATA_2] ;
  sisl_matrix_t *A = data[DATA_MATRIX_A] ;
  sisl_vector_t *rhs = data[DATA_RHS] ;
  gint ind = *((gint *)(data[DATA_INDEX])) ;
  BEM3DConfiguration *config = data[DATA_CONFIG] ;
  BEM3DWorkspace *work = data[DATA_WORK] ;
  static GArray *G = NULL ;
  static GArray *dGdn = NULL ;
  gdouble *x ;
  gint i, j, *idx, iu, il ;
  BEM3DParameters pt ;

  if ( G == NULL ) {
    G = g_array_new(FALSE, FALSE, sizeof(gdouble)) ;
    dGdn = g_array_new(FALSE, FALSE, sizeof(gdouble)) ;
  }

  g_array_set_size(G,bem3d_element_node_number(e)) ; 
  g_array_set_size(dGdn,bem3d_element_node_number(e)) ;
  bem3d_element_assemble_equations(e, p, config, &pt, G, dGdn, work) ;
  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    if ( (idx = g_hash_table_lookup(wakelines, bem3d_element_node(e,i)))
	 == NULL ) {
      j = bem3d_element_global_index(e,i) ;
      x = bem3d_mesh_data_get(dp, j) ;
      sisl_vector_set(rhs, ind,
		      sisl_vector_get(rhs,ind) +
		      g_array_index(dGdn,gdouble,i)*x[0]) ;
/*       sisl_vector_set(rhs, ind, */
/* 		      sisl_vector_get(rhs,ind) -  */
/* 		      g_array_index(dGdn,gdouble,i)*x[0]) ; */
    } else {
      iu = idx[0] ; il = idx[1] ;
      sisl_matrix_set(A, ind, iu,
		      sisl_matrix_get(A, ind, iu) -
		      g_array_index(dGdn,gdouble,i)) ;
      sisl_matrix_set(A, ind, il,
		      sisl_matrix_get(A, ind, il) +
		      g_array_index(dGdn,gdouble,i)) ;
/*       sisl_matrix_set(A, ind, iu, */
/* 		      sisl_matrix_get(A, ind, iu) +  */
/* 		      g_array_index(dGdn,gdouble,i)) ; */
/*       sisl_matrix_set(A, ind, il, */
/* 		      sisl_matrix_get(A, ind, il) -  */
/* 		      g_array_index(dGdn,gdouble,i)) ; */
    }
  }

  return 0 ;
}

static gint wake_node_influence(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMesh *m = data[DATA_MESH_1] ;
  BEM3DMesh *w = data[DATA_MESH_2] ;
  GHashTable *wakenodes = data[DATA_WAKE_HASH] ;
  gdouble ee = *((gdouble *)(data[DATA_EPSILON])) ;
  GtsPoint *p = data[DATA_POINT] ;
  gboolean wake_node = FALSE ;
  GtsVector n ;

 if ( (g_hash_table_lookup(wakenodes, GINT_TO_POINTER(i+1))) != NULL )
   wake_node = TRUE ;
   
  p->x = GTS_POINT(v)->x ; 
  p->y = GTS_POINT(v)->y ;
  p->z = GTS_POINT(v)->z ;

  if ( wake_node ) {
    bem3d_node_normal(m, i, n, BEM3D_AVERAGE_MWA) ;
    p->x += n[0]*ee ; p->y += n[1]*ee ; p->z += n[2]*ee ;
  }

  data[DATA_INDEX] = &i ;

  bem3d_mesh_foreach_element(w, (BEM3DElementFunc)wake_element_assemble, 
			     data) ;

  return 0 ;
}

static void mesh_wake_assemble(BEM3DMesh *m,
			       BEM3DConfiguration *config,
			       GPtrArray *wakes,
			       GPtrArray *wakedata,
			       GHashTable *wakenodes,
			       GHashTable *wakelines,
			       gdouble ee,
			       sisl_matrix_t *A,
			       sisl_vector_t *rhs,
			       BEM3DWorkspace *work)
			       

{
  gpointer data[DATA_WIDTH] = {NULL} ;
  gint i ;

  data[DATA_MESH_1] = m ;
  data[DATA_MATRIX_A] = A ; data[DATA_RHS] = rhs ;
  data[DATA_EPSILON] = &ee ;
  data[DATA_WAKE_HASH] = wakenodes ;
  data[DATA_WAKE_LINES] = wakelines ;
  data[DATA_POINT] = gts_point_new(gts_point_class(), 0, 0, 0) ;
  data[DATA_CONFIG] = config ;
  data[DATA_WORK] = work ;
  
  for ( i = 0 ; i < wakes->len ; i ++ ) {
    data[DATA_MESH_2] = g_ptr_array_index(wakes,i) ;
    data[DATA_MESH_DATA_2] = g_ptr_array_index(wakedata,i) ;
    bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)wake_node_influence, data) ;
  }

  return ;
}

static gint init_node_data(BEM3DElement *e, gpointer data[])

{
  BEM3DMeshData *d = data[2] ;
  gdouble *f ;
  gint i ;

  for ( i = 0 ; i < bem3d_element_node_number(e) ; i ++ ) {
    if ( bem3d_mesh_data_get(d, bem3d_element_global_index(e,i))
	 == NULL ) {
      bem3d_mesh_data_add_node(d, bem3d_element_global_index(e,i)) ;
      f = bem3d_mesh_data_get(d, bem3d_element_global_index(e,i)) ;
      f[0] = 0.0 ;
    }
  }

  return 0 ;
}

static void init_wake_node_data(BEM3DMesh *w, BEM3DMeshData *d, gdouble t)

{
  gpointer data[4] ;

  data[0] = w ; data[2] = d ;

  bem3d_mesh_foreach_element(w, (BEM3DElementFunc)init_node_data, data) ;

  return ;
}

static gint set_bc(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMotion *m = data[0] ;
  GArray *dphi = data[1] ;
  gdouble t = *((gdouble *)(data[2])) ;
  BEM3DMeshData *d = data[3] ;
  GtsVector n, u ;
  gdouble *dd ;

  bem3d_node_normal(bem3d_motion_mesh(m), i, n, BEM3D_AVERAGE_MWA) ;
  bem3d_motion_node_velocity(m, i, t, u) ;
  dd = bem3d_mesh_data_get(d,i) ;
  g_array_index(dphi, gdouble, i) = dd[1] = gts_vector_scalar(u,n) ;

  fprintf(stderr, "%d %lg (%lg, %lg, %lg) (%lg, %lg, %lg)\n",
	  i, dd[1],
	  n[0], n[1], n[2],
	  u[0], u[1], u[2]) ;
  
  /* g_assert(!isnan(dd[1])) ; */
  
  return 0 ;
}

static void mesh_set_bc(BEM3DMesh *m, 
			BEM3DMeshData *d,
			BEM3DMotion *mt, 
			GHashTable *wakenodes,
			gdouble t, GArray *dphi)

{
  gpointer data[4] ;

  g_assert(m == bem3d_motion_mesh(mt)) ;

  data[0] = mt ; data[1] = dphi ; data[2] = &t ; data[3] = d ;

  bem3d_mesh_foreach_node(m, (BEM3DNodeFunc)set_bc, data) ;

  return ;
}

static gint node_potential_jump(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMeshData *d = data[DATA_MESH_DATA_1] ;
  GHashTable *wakelines = data[DATA_WAKE_LINES] ;
  sisl_vector_t *phi = data[DATA_PHI] ;
  gint *idx, iu, il ;
  gdouble *f ;

  if ( (idx = g_hash_table_lookup(wakelines, v)) == NULL ) return 0 ;

  iu = idx[0] ; il = idx[1] ;

  f = bem3d_mesh_data_get(d, i) ;

  f[0] = sisl_vector_get(phi,iu) - sisl_vector_get(phi,il) ;
/*   f[0] = -(sisl_vector_get(phi,iu) - sisl_vector_get(phi,il)) ; */

  return 0 ;
}

static void wake_potential_jump(BEM3DMesh *w, BEM3DMeshData *d,
				GHashTable *wakelines, 
				sisl_vector_t *phi)

{
  gpointer data[DATA_WIDTH] ;

  data[DATA_MESH_DATA_1] = d ;
  data[DATA_WAKE_LINES] = wakelines ;
  data[DATA_PHI] = phi ;

  bem3d_mesh_foreach_node(w, (BEM3DNodeFunc)node_potential_jump, data) ;

  return ;
}

gint main(gint argc, gchar **argv)

{
  BEM3DMesh *m ;
  BEM3DWorkspace *work ;
  GPtrArray *meshes, *refmeshes, *wakes, *wakedata, *aerodata ;
  GPtrArray *edges, *edges1, *edges2, *motions ;
  GPtrArray *aerofiles, *adatafiles, *wakefiles, *wdatafiles ;
  GArray *wake_indices ;
  GHashTable *wake_nodes, *wake_lines ;
  gint *wake_line_indices, iwl ;
  sisl_matrix_t *A ;
  sisl_vector_t *phi, *rhs ;
  sisl_solver_workspace_t *ws ;
  sisl_solver_performance_t perf ;
  GArray *dphi ;
  GtsFile *fp ;
  GTimer *timer ;
  gint i, j, np, nc, s, s0 ;
  gchar ch, *ipfile, p[32], *progname ;
  gdouble *xd = NULL, ee, tol ;
  GString *fname ;
  BEM3DParameters gdata ;
  BEM3DEdge *edge ;
  BEM3DElementBuildFunc ebuild ;
  BEM3DMotion *motion ;
  BEM3DMeshData *data ;
  BEM3DConfiguration *config ;
  gpointer pt ;
  gdouble t, dt, t0, tf ;
  guint imin, imax ;
  GLogLevelFlags loglevel ;
  FILE *input, *output ;
  gboolean calc_potential ;

  wmpi_initialize(&argc, &argv) ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  /*computational geometry*/
  meshes = g_ptr_array_new() ;
  refmeshes = g_ptr_array_new() ;
  wakes = g_ptr_array_new() ;
  edges = g_ptr_array_new() ;
  edges1 = g_ptr_array_new() ;
  edges2 = g_ptr_array_new() ;
  motions = g_ptr_array_new() ;
  aerodata = g_ptr_array_new() ;
  wakedata = g_ptr_array_new() ;
  wake_indices = g_array_new(TRUE, TRUE, sizeof(gint)) ;
  wake_nodes = g_hash_table_new(NULL, NULL) ;
  fname = g_string_new("") ;
  aerofiles = g_ptr_array_new() ;
  adatafiles = g_ptr_array_new() ;
  wakefiles = g_ptr_array_new() ;
  wdatafiles = g_ptr_array_new() ;

  ebuild = bem3d_element_build_t1 ;

  bem3d_parameters_wavenumber(&gdata) = 0.0 ;
  bem3d_parameters_mach_number(&gdata) = 0.0 ;

  t0 = 0.0 ; tf = 1.0 ; dt = 1/16.0 ; s0 = 0 ; ee = 2.5e-2 ;
  calc_potential = TRUE ; tol = 1e-3 ;

  phi = rhs = NULL ; A = NULL ; dphi = NULL ;

  output = stdout ;
  loglevel = G_LOG_LEVEL_MESSAGE ;
  sprintf(p, "P%03d:", wmpi_rank()) ;

  bem3d_shapefunc_lookup_init() ; bem3d_configuration_init() ;
  bem3d_logging_init(stderr, p, loglevel, wmpi_shutdown) ;
  sisl_logging_init(stderr, p, loglevel, wmpi_shutdown) ;

  /*default configuration*/
  config = bem3d_configuration_new() ;

  while ( (ch = getopt(argc, argv, "hC:d:e:f:Gl:m:i:o:O:s:S:t:w:W:Z:")) 
	  != EOF ) {
    switch (ch) {
    default:
    case 'h':
      if ( wmpi_rank() == 0 ) {
	fprintf(stderr, 
		"%s: BEM3D aerodynamic calculation\n\n",
		progname) ;
	fprintf(stderr, "Usage: %s <options>\n", progname) ;
	fprintf(stderr, 
		"Options:\n"
		"        -h (print this message and exit)\n"
		"        -C <configuration file name>\n"
		"        -d # (time step)\n"
		"        -e <edge file name>\n"
		"        -f # (final time)\n"
		"        -G (graphical output: do not compute potential)\n"
		"        -i <bem3d input file>\n"
		"        -m <motion file name>\n"
		"        -o <output file name template>\n"
		"        -O <data output file name template>\n"
		/* "        -R rigid body motion\n" */
		"        -S # (start index)\n"
		"        -s # (start time)\n"
		"        -t # (iterative solver convergence tolerance)\n"
		"        -w <wake output file name template>\n"
		"        -W <wake data output file name template>\n"
		"        -Z # (Zhoukowsky condition parameter)\n") ;
      }
      wmpi_pause() ;
      wmpi_shutdown() ;
      return 0 ;
      break ;
    case 'l': 
      loglevel = 1 << atoi(optarg) ; 
      break ;
    case 'C':
      bem3d_configuration_init() ;
      bem3d_configuration_read(config, optarg) ;
      break ;
    case 'd': dt = atof(optarg) ; break ;
    case 'e':
      ipfile = g_strdup(optarg) ;
      edge = bem3d_edge_new(bem3d_edge_class()) ;
      input = file_open(ipfile, "-", "r", stdin) ;
      fp = gts_file_new(input) ;
      bem3d_edge_read(edge, fp) ;
      file_close(input) ;
      g_free(ipfile) ;
      if ( meshes->len == 0 ) {
	fprintf(stderr,
		"%s: edge files must be loaded after their meshes\n",
		progname) ;
 	return -1 ;
      }
      bem3d_edge_link_to_mesh(edge,
			      g_ptr_array_index(meshes, meshes->len-1)) ;
      g_ptr_array_add(edges, edge) ;
      edge = bem3d_edge_new(bem3d_edge_class()) ;
      g_ptr_array_add(edges1, edge) ;
      bem3d_edge_copy(g_ptr_array_index(edges1, edges1->len-1),
		      g_ptr_array_index(edges, edges->len-1)) ;

      edge = bem3d_edge_new(bem3d_edge_class()) ;
      g_ptr_array_add(edges2, edge) ;
      data = bem3d_mesh_data_sized_new(2, 32) ;
      g_ptr_array_add(wakedata, data) ;
      i = 0 ; g_array_append_val(wake_indices, i) ;
      g_ptr_array_add(wakes, 
		      bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
				     gts_edge_class(), gts_vertex_class())) ;
      
      break ;
    case 'f': tf = atof(optarg) ; break ;
    case 'G': calc_potential = FALSE ; break ;
    case 'i': 
      ipfile = g_strdup(optarg) ;
      m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
			 gts_edge_class(), gts_vertex_class()) ;
      input = file_open(ipfile, "-", "r", stdin) ;
      fp = gts_file_new(input) ;
      bem3d_mesh_read(m, fp) ;
      file_close(input) ; 

      g_ptr_array_add(meshes, m) ;
      data = bem3d_mesh_data_new(m, 2) ;
      g_ptr_array_add(aerodata, data) ;      

      input = file_open(ipfile, "-", "r", stdin) ;
      fp = gts_file_new(input) ;
      m = bem3d_mesh_new(bem3d_mesh_class(), gts_face_class(),
			 gts_edge_class(), gts_vertex_class()) ;
      bem3d_mesh_read(m, fp) ;
      file_close(input) ;
      g_free(ipfile) ;
      g_ptr_array_add(refmeshes, m) ;
      break ;
    case 'm':
      ipfile = g_strdup(optarg) ;
      if ( motions->len >= meshes->len ) {
	fprintf(stderr, 
		"%s: motion files must be loaded after their meshes\n",
		progname) ;
	return -1 ;
      }
      motion = bem3d_motion_new(bem3d_motion_class(),
				g_ptr_array_index(meshes, meshes->len-1),
				g_ptr_array_index(refmeshes, 
						  refmeshes->len-1)) ;
      input = file_open(ipfile, "-", "r", stdin) ;
      fp = gts_file_new(input) ;
      i = bem3d_motion_read(motion, fp) ;
      if ( i != 0 ) {
	fprintf(stderr, "%s: error reading motion at line %d: %s",
		progname, i, fp->error) ;
      }
      file_close(input) ;
      g_free(ipfile) ;
      g_ptr_array_add(motions, motion) ;
      break ;
    case 'o': g_ptr_array_add(aerofiles, g_strdup(optarg)) ; break ;
    case 'O': g_ptr_array_add(adatafiles, g_strdup(optarg)) ; break ;
    /* case 'R': rigid_body = TRUE ; */
    case 'S': s0 = atoi(optarg) ; break ;
    case 's': t0 = atof(optarg) ; break ;
    case 't': tol = atof(optarg) ; break ;
    case 'w': g_ptr_array_add(wakefiles, g_strdup(optarg)) ; break ;
    case 'W': g_ptr_array_add(wdatafiles, g_strdup(optarg)) ; break ;
    case 'Z': ee = atof(optarg) ; break ;
    }
  }
  bem3d_logging_init(stderr, p, loglevel, wmpi_shutdown) ;

  /*sanity checking here*/
  fprintf(stderr, "%s", BEM3D_STARTUP_MESSAGE) ;
  if ( meshes->len == 0 ) {
    fprintf(stderr, 
	    "%s: at least one body geometry must be specified\n", 
	    progname) ;
    return -1 ;
  }

  if ( aerofiles->len != meshes->len ) {
    fprintf(stderr, 
	    "%s: all geometry file name templates must be specified\n",
	    progname) ;
    return -1 ;
  }

  if ( adatafiles->len != meshes->len ) {
    fprintf(stderr, 
	    "%s: all data file name templates must be specified\n",
	    progname) ;
    return -1 ;
  }

  if ( wakefiles->len != edges->len ) {
    fprintf(stderr, 
	    "%s: all wake geometry file name templates must be specified\n",
	    progname) ;
    return -1 ;
  }

  if ( wdatafiles->len != edges->len ) {
    fprintf(stderr, 
	    "%s: all wake data file name templates must be specified\n",
	    progname) ;
    return -1 ;
  }

  if ( motions->len == 0 ) {
    fprintf(stderr, 
	    "%s: the mesh motions must be specified\n", 
	    progname) ;
    return -1 ;    
  }

  /*set up the solver workspace*/
  ws = sisl_solver_workspace_new() ;
  work = bem3d_workspace_new() ;
  
  for ( (i = 0), (np = 0) ; i < meshes->len ; i ++ ) {
    nc = bem3d_mesh_node_number(g_ptr_array_index(meshes,i)) ;
    if ( wmpi_rank() == 0 )
      fprintf(stderr, "%s: mesh %d, %d collocation points\n", progname, i, nc) ;
    np += nc ;
  }
  
  if ( np == 0 ) {
    if ( wmpi_rank() == 0 ) {
      fprintf(stderr, "%s: empty meshes read (no nodes): exiting\n", 
	      progname) ;
    }
    wmpi_pause() ;
    wmpi_shutdown() ;
    return 0 ;
  }

  if ( wmpi_rank() == 0 )
    fprintf(stderr, "%s: mesh collocation points: %d\n", progname, np) ;

  timer = g_timer_new() ; g_timer_start(timer) ;

  wmpi_split_range(0, np, &imin, &imax) ;
  fprintf(stderr, "%s: P%d: nodes: %d--%d\n", 
	  progname, wmpi_rank(), imin, imax) ;

  if ( calc_potential ) {
    phi = sisl_vector_new(SISL_REAL) ;
    rhs = sisl_vector_new(SISL_REAL) ;
    sisl_vector_set_length(phi, np) ; 
    sisl_vector_set_length(rhs, np) ; 
    dphi = g_array_sized_new(FALSE, FALSE, sizeof(gdouble), np) ;
    g_array_set_size(dphi, np) ;

    A = sisl_matrix_new(SISL_REAL, SISL_MATRIX_DENSE) ; 
    sisl_matrix_set_block_size(A, imax-imin, np) ; 
    sisl_matrix_row_number(A) = 
      sisl_matrix_column_number(A) = np ;
    sisl_matrix_local_row_start(A) = imin ; 
    sisl_matrix_local_row_end(A) = imax ;

    if ( wmpi_rank() == 0 )
      fprintf(stderr, "%s: matrix chunks split\n", progname) ;
  }

  for ( i = 0 ; i < motions->len ; i ++ )
    bem3d_motion_expand_defs(g_ptr_array_index(motions,i)) ;

  /*hash the wake nodes*/
  for ( i = iwl = 0 ; i < edges->len ; i ++ ) {
    edge = g_ptr_array_index(edges,i) ;
    for ( j = 0 ; j < bem3d_edge_node_number(edge) ; j ++ ) {
      g_hash_table_insert(wake_nodes,
			  GINT_TO_POINTER(bem3d_edge_node_index_upper(edge,j)+1),
			  GINT_TO_POINTER(bem3d_edge_node_index_upper(edge,j)+1)) ;
      g_hash_table_insert(wake_nodes,
			  GINT_TO_POINTER(bem3d_edge_node_index_lower(edge,j)+1),
			  GINT_TO_POINTER(bem3d_edge_node_index_lower(edge,j)+1)) ;
    }
    iwl += bem3d_edge_node_number(edge) ;
  }

  wake_line_indices = (gint *)g_malloc(2*iwl*sizeof(gint)) ;

  for ( (s = s0), (t = t0+dt) ; t <= tf ; (s ++), (t += dt) ) {
    fprintf(stderr, "%s i=%d; t=%lg;\n", p, s, t) ;

    /*position the meshes*/
    for ( i = 0 ; i < meshes->len ; i ++ ) {
      bem3d_motion_create_evaluators(g_ptr_array_index(motions,i)) ;
      bem3d_motion_mesh_position(g_ptr_array_index(motions,i), t) ;
      g_string_printf(fname, g_ptr_array_index(aerofiles,i), s) ;
      output = file_open(fname->str, "-", "w", stdout) ;
      bem3d_mesh_write(g_ptr_array_index(meshes,i), output) ;
      file_close(output) ;
    }
    /*extend the wakes*/
    wake_lines = g_hash_table_new(NULL, NULL) ;
    iwl = 0 ;
    for ( i = 0 ; i < edges->len ; i ++ ) {
      bem3d_edge_copy(g_ptr_array_index(edges2, i),
		      g_ptr_array_index(edges, i)) ;
      bem3d_link_edges(g_ptr_array_index(wakes, i),
		       g_ptr_array_index(edges1, i),
		       g_ptr_array_index(edges2, i),
		       ebuild, 2, 3,
		       &(g_array_index(wake_indices, gint, i))) ;
      bem3d_edge_clear(g_ptr_array_index(edges1, i)) ;
      pt = g_ptr_array_index(edges2, i) ;
      g_ptr_array_index(edges2, i) = g_ptr_array_index(edges1, i) ;
      g_ptr_array_index(edges1, i) = pt ;
      init_wake_node_data(g_ptr_array_index(wakes, i),
			  g_ptr_array_index(wakedata, i), t) ;

      edge = g_ptr_array_index(edges1,i) ;
      for ( j = 0 ; j < bem3d_edge_node_number(edge) ; (j ++), (iwl += 2) ) {
	wake_line_indices[iwl+0] = bem3d_edge_node_index_upper(edge,j) ;
	wake_line_indices[iwl+1] = bem3d_edge_node_index_lower(edge,j) ;
	g_hash_table_insert(wake_lines,
			    bem3d_edge_vertex(edge,j),
			    &(wake_line_indices[iwl])) ;
      }
    }

    if ( calc_potential ) {
      sisl_vector_set_all(rhs, 0.0) ;
      for ( i = 0 ; i < motions->len ; i ++ ) {
	mesh_set_bc(g_ptr_array_index(meshes, i),
		    g_ptr_array_index(aerodata, i),
		    g_ptr_array_index(motions, i),
		    wake_nodes, t, dphi) ;
      }

      /*mesh-mesh interactions and boundary conditions*/
      for ( i = 0 ; i < meshes->len ; i ++ )
	for ( j = 0 ; j < meshes->len ; j ++ )
	  mesh_mesh_assemble(g_ptr_array_index(meshes,i),
			     g_ptr_array_index(meshes,j),
			     config, 
			     wake_nodes,
			     A, ee, dphi, rhs, work) ;

      /*mesh-wake interaction and Kutta condition*/
      for ( i = 0 ; i < meshes->len ; i ++ ) {
	mesh_wake_assemble(g_ptr_array_index(meshes,i),
			   config,
			   wakes, wakedata,
			   wake_nodes, wake_lines, ee, A, rhs, work) ;
      }
      sisl_solve(SISL_SOLVER_BICGSTAB, A, phi, rhs, tol, 128, ws, &perf) ;

      /*calculate the potential jump on the new wake vertices*/
      for ( i = 0 ; i < wakes->len ; i ++ ) 
	wake_potential_jump(g_ptr_array_index(wakes,i),
			    g_ptr_array_index(wakedata,i),
			    wake_lines, phi) ;

      for ( i = 0 ; i < np ; i ++ ) {
	for ( j = 0 ; j < aerodata->len ; j ++ ) {
	  if ( (xd = bem3d_mesh_data_get(g_ptr_array_index(aerodata, j), i))
	       != NULL ) {
	    xd[0] = sisl_vector_get(phi, i) ;
	    break ;
	  }
	}
	if ( xd == NULL ) {
	  g_error("cannot find point %d in mesh data", i) ;
	}
      }
    }

    for ( i = 0 ; i < edges->len ; i ++ ) {    
      g_string_printf(fname, g_ptr_array_index(wakefiles, i), s) ;
      output = file_open(fname->str, "-", "w", stdout) ;
      bem3d_mesh_write(g_ptr_array_index(wakes,i), output) ;
      file_close(output) ;
      g_string_printf(fname, g_ptr_array_index(wdatafiles, i), s) ;
      output = file_open(fname->str, "-", "w", stdout) ;
      bem3d_mesh_data_write(g_ptr_array_index(wakedata, i), output, NULL) ;
      file_close(output) ;
    }

    for ( i = 0 ; i < aerodata->len ; i ++ ) {    
      g_string_printf(fname, g_ptr_array_index(adatafiles, i), s) ;
      output = file_open(fname->str, "-", "w", stdout) ;
      bem3d_mesh_data_write(g_ptr_array_index(aerodata,i), output, NULL) ;
      file_close(output) ;
    }

    for ( i = 0 ; i < motions->len ; i ++ )
      bem3d_motion_free_evaluators(g_ptr_array_index(motions,i)) ;

    g_hash_table_destroy(wake_lines) ;
  }

  g_ptr_array_free(meshes, FALSE) ; 
  g_ptr_array_free(wakes, FALSE) ; 
  g_ptr_array_free(edges, FALSE) ; 
  g_free(progname) ;

  wmpi_pause() ;
  wmpi_shutdown() ;

  return 0 ;
}
