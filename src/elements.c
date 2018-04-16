/* elements.c
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

/**
 * @defgroup elements Building elements
 *
 * BEM3D can handle a wide range of elements and some of them are
 * built in to the main library. The BEMElement type contains the
 * geometric and other information needed to perform calculations; the
 * ::BEM3DElementBuildFunc type defines a function which takes as input
 * the vertices of an element and returns the new element. In each
 * case, the inputs are an array of GtsEdge's and an array of
 * GtsVertex's which define the boundary of the element, with the
 * exception of ::bem3d_element_build_t3 and ::bem3d_element_build_q2
 * which include an interior point in the centre of the element.
 *
 * @{
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

static GtsEdge *connect_vertices(GtsEdgeClass *klass, 
				 GtsVertex *v1, GtsVertex *v2)

{
  GtsEdge *e ;

  if ( (e = GTS_EDGE(gts_vertices_are_connected(v1, v2))) == NULL )
    e = gts_edge_new(klass, v1, v2) ;

  return e ;
}

/** 
 * Generate a zero order triangular element.  The shape function for
 * computation is ::bem3d_shfunc_t0 and for the geometry is
 * ::bem3d_shfunc_t1.
 *
 * @param e array of GtsEdge around the outside of the element
 * @param v array of GtsVertex around the outside of the element
 * 
 * @return a pointer to the new element
 */

BEM3DElement *bem3d_element_build_t0(GtsEdge **e, GtsVertex **v)

{
  BEM3DElement *el ;
  GtsFace *f ;
  GtsEdge *e1, *e2, *e3 ;
  GtsVertex *v1 ;

  g_return_val_if_fail(e != NULL, NULL) ;
  g_return_val_if_fail(v != NULL, NULL) ;

  el = bem3d_element_new(bem3d_element_class(),
		       3, 3, 1, 3,
		       bem3d_shfunc_t1, bem3d_shfunc_t0) ;

  bem3d_element_add_vertex(el, v[0], 0) ;
  bem3d_element_set_corner(el, 0, 0) ;
  bem3d_element_add_vertex(el, v[1], 1) ;
  bem3d_element_set_corner(el, 1, 1) ;
  bem3d_element_add_vertex(el, v[2], 2) ;
  bem3d_element_set_corner(el, 2, 2) ;

  v1 = gts_vertex_new(gts_vertex_class(), 
		      (GTS_POINT(v[0])->x+
		       GTS_POINT(v[1])->x+
		       GTS_POINT(v[2])->x)/3.0,
		      (GTS_POINT(v[0])->y+
		       GTS_POINT(v[1])->y+
		       GTS_POINT(v[2])->y)/3.0,
		      (GTS_POINT(v[0])->z+
		       GTS_POINT(v[1])->z+
		       GTS_POINT(v[2])->z)/3.0) ;

  e1 = gts_edge_new(gts_edge_class(), v[0], v1) ;
  e2 = gts_edge_new(gts_edge_class(), v[1], v1) ;
  e3 = gts_edge_new(gts_edge_class(), v[2], v1) ;

  f = gts_face_new(gts_face_class(), e[0], e2, e1) ;
  bem3d_element_add_face(el, f, 0) ;

  f = gts_face_new(gts_face_class(), e[1], e3, e2) ;
  bem3d_element_add_face(el, f, 1) ;

  f = gts_face_new(gts_face_class(), e[2], e1, e3) ;
  bem3d_element_add_face(el, f, 2) ;

  bem3d_element_add_node(el, v1, 0) ;

  bem3d_element_node_xi(el,0) = 1/3.0 ;
  bem3d_element_node_eta(el,0) = 1/3.0 ;

  bem3d_element_vertex_xi(el,0) = 0.0 ;
  bem3d_element_vertex_eta(el,0) = 0.0 ;
  bem3d_element_vertex_xi(el,1) = 1.0 ;
  bem3d_element_vertex_eta(el,1) = 0.0 ;
  bem3d_element_vertex_xi(el,2) = 0.0 ;
  bem3d_element_vertex_eta(el,2) = 1.0 ;
  
  return el ;
}

/** 
 * Generate a first order triangular element. The shape function for
 * geometry and computation is ::bem3d_shfunc_t1.
 * 
 * @param e array of GtsEdge around the outside of the element
 * @param v array of GtsVertex around the outside of the element
 * 
 * @return a pointer to the new ::BEM3DElement
 */

BEM3DElement *bem3d_element_build_t1(GtsEdge **e, GtsVertex **v)

{
  BEM3DElement *el ;
  GtsFace *f ;

  g_return_val_if_fail(e != NULL, NULL) ;
  g_return_val_if_fail(v != NULL, NULL) ;
  g_return_val_if_fail(v[0] != NULL && v[1] != NULL && v[2] != NULL, NULL) ;

  if ( e[0] == NULL ) e[0] = connect_vertices(gts_edge_class(), v[0], v[1]) ;
  if ( e[1] == NULL ) e[1] = connect_vertices(gts_edge_class(), v[1], v[2]) ;
  if ( e[2] == NULL ) e[2] = connect_vertices(gts_edge_class(), v[2], v[0]) ;

  if ( !gts_segment_connect(GTS_SEGMENT(e[0]), v[0], v[1]) )
    g_error("%s: e[0] does not connect v[0] and v[1];", __FUNCTION__) ;
  if ( !gts_segment_connect(GTS_SEGMENT(e[1]), v[1], v[2]) )
    g_error("%s: e[1] does not connect v[1] and v[2];", __FUNCTION__) ;
  if ( !gts_segment_connect(GTS_SEGMENT(e[2]), v[2], v[0]) )
    g_error("%s: e[2] does not connect v[2] and v[0];", __FUNCTION__) ;

  el = bem3d_element_new(bem3d_element_class(),
			 1, 3, 0, 3,
			 bem3d_shfunc_t1, bem3d_shfunc_t1) ;

  bem3d_element_add_vertex(el, v[0], 0) ;
  bem3d_element_set_corner(el, 0, 0) ;
  bem3d_element_add_vertex(el, v[1], 1) ;
  bem3d_element_set_corner(el, 1, 1) ;
  bem3d_element_add_vertex(el, v[2], 2) ;
  bem3d_element_set_corner(el, 2, 2) ;

  f = gts_face_new(gts_face_class(), e[0], e[1], e[2]) ;

  bem3d_element_add_face(el, f, 0) ;

  /*the shape and collocation point data are shared: we only need to
    deal with one of them*/
  bem3d_element_vertex_xi(el,0) = 0.0 ;
  bem3d_element_vertex_eta(el,0) = 0.0 ;
  bem3d_element_vertex_xi(el,1) = 1.0 ;
  bem3d_element_vertex_eta(el,1) = 0.0 ;
  bem3d_element_vertex_xi(el,2) = 0.0 ;
  bem3d_element_vertex_eta(el,2) = 1.0 ;

  return el ;
}

/** 
 * Generate a second order triangular element.  The shape function for
 * geometry and computation is ::bem3d_shfunc_t2.
 * 
 * @param e array of GtsEdge around the outside of the element
 * @param v array of GtsVertex around the outside of the element
 * 
 * @return a pointer to the new element
 */

BEM3DElement *bem3d_element_build_t2(GtsEdge **e, GtsVertex **v)

{
  BEM3DElement *el ;
  GtsEdge *e1, *e2, *e3 ;
  GtsFace *f ;

  g_return_val_if_fail(e != NULL, NULL) ;
  g_return_val_if_fail(v != NULL, NULL) ;

  el = bem3d_element_new(bem3d_element_class(),
		       4, 6, 0, 3,
		       bem3d_shfunc_t2, bem3d_shfunc_t2) ;

  e[0] = connect_vertices(gts_edge_class(), v[0], v[1]) ;
  e[1] = connect_vertices(gts_edge_class(), v[1], v[2]) ;
  e[2] = connect_vertices(gts_edge_class(), v[2], v[3]) ;
  e[3] = connect_vertices(gts_edge_class(), v[3], v[4]) ;
  e[4] = connect_vertices(gts_edge_class(), v[4], v[5]) ;
  e[5] = connect_vertices(gts_edge_class(), v[5], v[0]) ;

  e1 = connect_vertices(gts_edge_class(), v[1], v[5]) ;
  e2 = connect_vertices(gts_edge_class(), v[1], v[3]) ;
  e3 = connect_vertices(gts_edge_class(), v[3], v[5]) ;

  bem3d_element_add_vertex(el, v[0], 0) ;
  bem3d_element_set_corner(el, 0, 0) ;
  bem3d_element_add_vertex(el, v[2], 1) ;
  bem3d_element_set_corner(el, 1, 1) ;
  bem3d_element_add_vertex(el, v[4], 2) ;
  bem3d_element_set_corner(el, 2, 2) ;

  bem3d_element_add_vertex(el, v[1], 3) ;
  bem3d_element_add_vertex(el, v[3], 4) ;
  bem3d_element_add_vertex(el, v[5], 5) ;

  f = gts_face_new(gts_face_class(), e[0], e1, e[5]) ;
  bem3d_element_add_face(el, f, 0) ;

  f = gts_face_new(gts_face_class(), e[1], e[2], e2) ;
  bem3d_element_add_face(el, f, 1) ;

  f = gts_face_new(gts_face_class(), e1, e2, e3) ;
  bem3d_element_add_face(el, f, 2) ;

  f = gts_face_new(gts_face_class(), e3, e[3], e[4]) ;
  bem3d_element_add_face(el, f, 3) ;

  bem3d_element_vertex_xi(el,0) = 0.0 ;
  bem3d_element_vertex_eta(el,0) = 0.0 ;
  bem3d_element_vertex_xi(el,1) = 1.0 ;
  bem3d_element_vertex_eta(el,1) = 0.0 ;
  bem3d_element_vertex_xi(el,2) = 0.0 ;
  bem3d_element_vertex_eta(el,2) = 1.0 ;
  bem3d_element_vertex_xi(el,3) = 0.5 ;
  bem3d_element_vertex_eta(el,3) = 0.0 ;
  bem3d_element_vertex_xi(el,4) = 0.5 ;
  bem3d_element_vertex_eta(el,4) = 0.5 ;
  bem3d_element_vertex_xi(el,5) = 0.0 ;
  bem3d_element_vertex_eta(el,5) = 0.5 ;

  return el ;
}

/** 
 * Generate a third order triangular element.  The shape function for
 * geometry and computation is ::bem3d_shfunc_t3. If the tenth GtsVertex
 * in \a v is NULL, the centre vertex is the centroid of the supplied
 * vertices, otherwise \a v[9] is used. 
 * 
 * @param e array of GtsEdge around the outside of the element
 * @param v array of GtsVertex around the outside of the element
 * 
 * @return a pointer to the new element
 */

BEM3DElement *bem3d_element_build_t3(GtsEdge **e, GtsVertex **v)

{
  BEM3DElement *el ;
  GtsEdge *e1, *e2, *e3, *e4, *e5, *e6, *e7, *e8, *e9 ;
  GtsFace *f ;
  GtsVertex *v1 ;

  g_return_val_if_fail(e != NULL, NULL) ;
  g_return_val_if_fail(v != NULL, NULL) ;

  if ( v[9] == NULL ) 
    v1 = gts_vertex_new(gts_vertex_class(),
			(GTS_POINT(v[0])->x+
			 GTS_POINT(v[3])->x+
			 GTS_POINT(v[6])->x)/3.0,
			(GTS_POINT(v[0])->y+
			 GTS_POINT(v[3])->y+
			 GTS_POINT(v[6])->y)/3.0,
			(GTS_POINT(v[0])->z+
			 GTS_POINT(v[3])->z+
			 GTS_POINT(v[6])->z)/3.0) ;
  else 
    v1 = v[9] ;

  el = bem3d_element_new(bem3d_element_class(),
			 9, 10, 0, 3, bem3d_shfunc_t3, bem3d_shfunc_t3) ;

  /*element corners*/
  bem3d_element_add_vertex(el, v[0], 0) ;
  bem3d_element_set_corner(el, 0, 0) ;
  bem3d_element_add_vertex(el, v[3], 1) ;
  bem3d_element_set_corner(el, 1, 1) ;
  bem3d_element_add_vertex(el, v[6], 2) ;
  bem3d_element_set_corner(el, 2, 2) ;

  /*element edges*/
  bem3d_element_add_vertex(el, v[1], 3) ;
  bem3d_element_add_vertex(el, v[2], 4) ;
  bem3d_element_add_vertex(el, v[4], 5) ;
  bem3d_element_add_vertex(el, v[5], 6) ;
  bem3d_element_add_vertex(el, v[7], 7) ;
  bem3d_element_add_vertex(el, v[8], 8) ;
  bem3d_element_add_vertex(el, v1,   9) ;

  /*edges*/
  e1 = gts_edge_new(gts_edge_class(), el->v[3], el->v[8]) ;
  e2 = gts_edge_new(gts_edge_class(), el->v[3], el->v[9]) ;
  e3 = gts_edge_new(gts_edge_class(), el->v[4], el->v[9]) ;
  e4 = gts_edge_new(gts_edge_class(), el->v[4], el->v[5]) ;
  e5 = gts_edge_new(gts_edge_class(), el->v[8], el->v[9]) ;
  e6 = gts_edge_new(gts_edge_class(), el->v[9], el->v[5]) ;
  e7 = gts_edge_new(gts_edge_class(), el->v[7], el->v[9]) ;
  e8 = gts_edge_new(gts_edge_class(), el->v[6], el->v[9]) ;
  e9 = gts_edge_new(gts_edge_class(), el->v[7], el->v[6]) ;
  
  f = gts_face_new(gts_face_class(), e[0], e1, e[8]) ;
  bem3d_element_add_face(el, f, 0) ;

  f = gts_face_new(gts_face_class(), e1, e2, e5) ;
  bem3d_element_add_face(el, f, 1) ;

  f = gts_face_new(gts_face_class(), e[1], e3, e2) ;
  bem3d_element_add_face(el, f, 2) ;

  f = gts_face_new(gts_face_class(), e3, e4, e6) ;
  bem3d_element_add_face(el, f, 3) ;

  f = gts_face_new(gts_face_class(), e[2], e[3], e4) ;
  bem3d_element_add_face(el, f, 4) ;

  f = gts_face_new(gts_face_class(), e5, e7, e[7]) ;
  bem3d_element_add_face(el, f, 5) ;

  f = gts_face_new(gts_face_class(), e7, e8, e9) ;
  bem3d_element_add_face(el, f, 6) ;

  f = gts_face_new(gts_face_class(), e6, e[4], e8) ;
  bem3d_element_add_face(el, f, 7) ;

  f = gts_face_new(gts_face_class(), e[5], e[6], e9) ;
  bem3d_element_add_face(el, f, 8) ;

  bem3d_element_vertex_xi(el,0) = 0.0 ;
  bem3d_element_vertex_eta(el,0) = 0.0 ;

  bem3d_element_vertex_xi(el,1) = 1.0 ;
  bem3d_element_vertex_eta(el,1) = 0.0 ;

  bem3d_element_vertex_xi(el,2) = 0.0 ;
  bem3d_element_vertex_eta(el,2) = 1.0 ;

  bem3d_element_vertex_xi(el,3) = 1.0/3.0 ;
  bem3d_element_vertex_eta(el,3) = 0.0 ;

  bem3d_element_vertex_xi(el,4) = 2.0/3.0 ;
  bem3d_element_vertex_eta(el,4) = 0.0 ;

  bem3d_element_vertex_xi(el,5) = 2.0/3.0 ;
  bem3d_element_vertex_eta(el,5) = 1.0/3.0 ;

  bem3d_element_vertex_xi(el,6) = 1.0/3.0 ;
  bem3d_element_vertex_eta(el,6) = 2.0/3.0 ;

  bem3d_element_vertex_xi(el,7) = 0.0 ;
  bem3d_element_vertex_eta(el,7) = 2.0/3.0 ;

  bem3d_element_vertex_xi(el,8) = 0.0 ;
  bem3d_element_vertex_eta(el,8) = 1.0/3.0 ;

  bem3d_element_vertex_xi(el,9) = 1.0/3.0 ;
  bem3d_element_vertex_eta(el,9) = 1.0/3.0 ;

  return el ;
}

/** 
 * Generate a first order quadrilateral element.  The shape function
 * for geometry and computation is ::bem3d_shfunc_q1.
 * 
 * @param e array of GtsEdge around the outside of the element
 * @param v array of GtsVertex around the outside of the element
 * 
 * @return a pointer to the new element
 */

BEM3DElement *bem3d_element_build_q1(GtsEdge **e, GtsVertex **v)

{
  BEM3DElement *el ;
  GtsFace *f ;

  g_return_val_if_fail(e != NULL, NULL) ;
  g_return_val_if_fail(v != NULL, NULL) ;

  if ( e[0] == NULL ) e[0] = connect_vertices(gts_edge_class(), v[0], v[1]) ;
  if ( e[1] == NULL ) e[1] = connect_vertices(gts_edge_class(), v[1], v[2]) ;
  if ( e[2] == NULL ) e[2] = connect_vertices(gts_edge_class(), v[2], v[3]) ;
  if ( e[3] == NULL ) e[3] = connect_vertices(gts_edge_class(), v[3], v[0]) ;
  if ( e[4] == NULL ) e[4] = connect_vertices(gts_edge_class(), v[0], v[2]) ;

  if ( !gts_segment_connect(GTS_SEGMENT(e[0]), v[0], v[1]) )
    g_error("%s: e[0] does not connect v[0] and v[1];", __FUNCTION__) ;
  if ( !gts_segment_connect(GTS_SEGMENT(e[1]), v[1], v[2]) )
    g_error("%s: e[1] does not connect v[1] and v[2];", __FUNCTION__) ;
  if ( !gts_segment_connect(GTS_SEGMENT(e[2]), v[2], v[3]) )
    g_error("%s: e[2] does not connect v[2] and v[3];", __FUNCTION__) ;
  if ( !gts_segment_connect(GTS_SEGMENT(e[3]), v[3], v[0]) )
    g_error("%s: e[3] does not connect v[3] and v[0];", __FUNCTION__) ;
  if ( !gts_segment_connect(GTS_SEGMENT(e[4]), v[0], v[2]) )
    g_error("%s: e[4] does not connect v[0] and v[2];", __FUNCTION__) ;

  el = bem3d_element_new(bem3d_element_class(),
			 2, 4, 0, 4, bem3d_shfunc_q1, bem3d_shfunc_q1) ;  
  bem3d_element_add_vertex(el, v[0], 0) ;
  bem3d_element_add_vertex(el, v[1], 1) ;
  bem3d_element_add_vertex(el, v[2], 2) ;
  bem3d_element_add_vertex(el, v[3], 3) ;
  bem3d_element_set_corner(el, 0, 0) ;
  bem3d_element_set_corner(el, 1, 1) ;
  bem3d_element_set_corner(el, 2, 2) ;
  bem3d_element_set_corner(el, 3, 3) ;

  f = gts_face_new(gts_face_class(), e[0], e[1], e[4]) ;
  bem3d_element_add_face(el, f, 0) ;
  f = gts_face_new(gts_face_class(), e[4], e[2], e[3]) ;
  bem3d_element_add_face(el, f, 1) ;

  bem3d_element_vertex_xi(el,0) = 0.0 ;
  bem3d_element_vertex_eta(el,0) = 0.0 ;
  bem3d_element_vertex_xi(el,1) = 1.0 ;
  bem3d_element_vertex_eta(el,1) = 0.0 ;
  bem3d_element_vertex_xi(el,2) = 1.0 ;
  bem3d_element_vertex_eta(el,2) = 1.0 ;
  bem3d_element_vertex_xi(el,3) = 0.0 ;
  bem3d_element_vertex_eta(el,3) = 1.0 ;

  return el ;
}

/** 
 * Generate a second order quadrilateral element, based on the GMSH
 * nine point quad. The shape function for geometry and computation is
 * ::bem3d_shfunc_q2. 
 * 
 * @param e array of GtsEdge around the outside of the element
 * @param v array of GtsVertex around the outside of the element
 * 
 * @return a pointer to the new element
 */

BEM3DElement *bem3d_element_build_q2(GtsEdge **e, GtsVertex **v)

{
  BEM3DElement *el ;
  GtsEdge *e1, *e2, *e3, *e4, *e5, *e6, *e7, *e8 ;
  GtsFace *f ;

  g_return_val_if_fail(e != NULL, NULL) ;
  g_return_val_if_fail(v != NULL, NULL) ;

  el = bem3d_element_new(bem3d_element_class(),
		       8, 9, 0, 4, bem3d_shfunc_q2, bem3d_shfunc_q2) ;  

  bem3d_element_add_vertex(el, v[0], 0) ;
  bem3d_element_add_vertex(el, v[2], 1) ;
  bem3d_element_add_vertex(el, v[4], 2) ;
  bem3d_element_add_vertex(el, v[6], 3) ;

  bem3d_element_add_vertex(el, v[1], 4) ;
  bem3d_element_add_vertex(el, v[3], 5) ;
  bem3d_element_add_vertex(el, v[5], 6) ;
  bem3d_element_add_vertex(el, v[7], 7) ;

  if ( v[8] == NULL ) {
    v[8] = gts_vertex_new(gts_vertex_class(),
			  (GTS_POINT(v[0])->x+
			   GTS_POINT(v[2])->x+
			   GTS_POINT(v[4])->x+
			   GTS_POINT(v[6])->x)*0.25,
			  (GTS_POINT(v[0])->y+
			   GTS_POINT(v[2])->y+
			   GTS_POINT(v[4])->y+
			   GTS_POINT(v[6])->y)*0.25,
			  (GTS_POINT(v[0])->z+
			   GTS_POINT(v[2])->z+
			   GTS_POINT(v[4])->z+
			   GTS_POINT(v[6])->z)*0.25) ;
  }

  bem3d_element_add_vertex(el, v[8], 8) ;

  bem3d_element_set_corner(el, 0, 0) ;
  bem3d_element_set_corner(el, 1, 1) ;
  bem3d_element_set_corner(el, 2, 2) ;
  bem3d_element_set_corner(el, 3, 3) ;

  e[0] = connect_vertices(gts_edge_class(), v[0], v[1]) ;
  e[1] = connect_vertices(gts_edge_class(), v[1], v[2]) ;
  e[2] = connect_vertices(gts_edge_class(), v[2], v[3]) ;
  e[3] = connect_vertices(gts_edge_class(), v[3], v[4]) ;
  e[4] = connect_vertices(gts_edge_class(), v[4], v[5]) ;
  e[5] = connect_vertices(gts_edge_class(), v[5], v[6]) ;
  e[6] = connect_vertices(gts_edge_class(), v[6], v[7]) ;
  e[7] = connect_vertices(gts_edge_class(), v[7], v[0]) ;
  
  e1 = connect_vertices(gts_edge_class(), v[0], v[8]) ;
  e2 = connect_vertices(gts_edge_class(), v[1], v[8]) ;
  e3 = connect_vertices(gts_edge_class(), v[2], v[8]) ;
  e4 = connect_vertices(gts_edge_class(), v[3], v[8]) ;
  e5 = connect_vertices(gts_edge_class(), v[4], v[8]) ;
  e6 = connect_vertices(gts_edge_class(), v[5], v[8]) ;
  e7 = connect_vertices(gts_edge_class(), v[6], v[8]) ;
  e8 = connect_vertices(gts_edge_class(), v[7], v[8]) ;

  f = gts_face_new(gts_face_class(), e[0], e2, e1) ;
  bem3d_element_add_face(el, f, 0) ;
  f = gts_face_new(gts_face_class(), e[1], e3, e2) ;
  bem3d_element_add_face(el, f, 1) ;
  f = gts_face_new(gts_face_class(), e[2], e4, e3) ;
  bem3d_element_add_face(el, f, 2) ;
  f = gts_face_new(gts_face_class(), e[3], e5, e4) ;
  bem3d_element_add_face(el, f, 3) ;
  f = gts_face_new(gts_face_class(), e[4], e6, e5) ;
  bem3d_element_add_face(el, f, 4) ;
  f = gts_face_new(gts_face_class(), e[5], e7, e6) ;
  bem3d_element_add_face(el, f, 5) ;
  f = gts_face_new(gts_face_class(), e[6], e8, e7) ;
  bem3d_element_add_face(el, f, 6) ;
  f = gts_face_new(gts_face_class(), e[7], e1, e8) ;
  bem3d_element_add_face(el, f, 7) ;

  bem3d_element_vertex_xi(el,0) = -1.0 ;
  bem3d_element_vertex_eta(el,0) = -1.0 ;
  bem3d_element_vertex_xi(el,1) = 1.0 ;
  bem3d_element_vertex_eta(el,1) = -1.0 ;
  bem3d_element_vertex_xi(el,2) = 1.0 ;
  bem3d_element_vertex_eta(el,2) = 1.0 ;
  bem3d_element_vertex_xi(el,3) = -1.0 ;
  bem3d_element_vertex_eta(el,3) = 1.0 ;

  bem3d_element_vertex_xi(el,4) = 0.0 ;
  bem3d_element_vertex_eta(el,4) = -1.0 ;
  bem3d_element_vertex_xi(el,5) = 1.0 ;
  bem3d_element_vertex_eta(el,5) = 0.0 ;
  bem3d_element_vertex_xi(el,6) = 0.0 ;
  bem3d_element_vertex_eta(el,6) = 1.0 ;
  bem3d_element_vertex_xi(el,7) = -1.0 ;
  bem3d_element_vertex_eta(el,7) = 0.0 ;

  bem3d_element_vertex_xi(el,8) = 0.0 ;
  bem3d_element_vertex_eta(el,8) = 0.0 ;

  return el ;
}

/**
 * @}
 * 
 */
