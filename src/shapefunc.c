/* shapefunc.c
 * 
 * Copyright (C) 2006, 2018 Michael Carley
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
#include <math.h>

#include <gts.h>
#include <glib.h>

#include "bem3d.h"
#include "bem3d-private.h"

gdouble t3_shape_cffts[] = {
  1.0, -5.5, -5.5,  9.0,  18.0,  9.0,  -4.5,-13.5, -13.5,  -4.5,
  0.0,  1.0,  0.0, -4.5,   0.0,  0.0,   4.5,  0.0,   0.0,   0.0,
  0.0,  0.0,  1.0,  0.0,   0.0, -4.5,   0.0,  0.0,   0.0,   4.5,
  0.0,  9.0,  0.0,-22.5, -22.5,  0.0,  13.5, 27.0,  13.5,   0.0,
  0.0, -4.5,  0.0, 18.0,   4.5,  0.0, -13.5,-13.5,   0.0,   0.0,
  0.0,  0.0,  0.0,  0.0,  -4.5,  0.0,   0.0, 13.5,   0.0,   0.0,
  0.0,  0.0,  0.0,  0.0,  -4.5,  0.0,   0.0,  0.0,  13.5,   0.0,
  0.0,  0.0, -4.5,  0.0,   4.5, 18.0,   0.0,  0.0, -13.5, -13.5,
  0.0,  0.0,  9.0,  0.0, -22.5,-22.5,   0.0, 13.5,  27.0,  13.5,
  0.0,  0.0,  0.0,  0.0,  27.0,  0.0,   0.0,-27.0, -27.0,   0.0
} ;

/**
 * @defgroup shapefunc Shape functions
 *
 * Shape functions for triangular and quadrilateral elements, mainly
 * taken from the FEMPACK library of John Burkardt:
 *
 * http://www.csit.fsu.edu/~burkardt/m_src/fempack/fempack.html
 * http://orion.math.iastate.edu/burkardt/f_src/fempack/fempack.htm
 *
 * Shape functions are called with the arguments \a s and \a t and
 * return the shape functions and derivatives in arrays which must be
 * properly allocated by the user before the function
 * call. User-defined shape functions can be used but must be of the
 * BEM3DShapeFunc type. They must also be added to the internal lookup
 * table, using ::bem3d_shapefunc_lookup_add before any elements using
 * them can be read.
 *
 * For mathematical validity and in order for other functions in the
 * library to work, the shape functions should meet the requirements
 * of `Shape function magic', 
 * http://caswww.colorado.edu/courses.d/IFEM.d/IFEM.Ch18.d/IFEM.Ch18.pdf
 *
 * The definition of ::BEM3DShapeFunc has provision for passing data to
 * the shape function. Within bem3dlib, at present, all calls set the
 * data field to NULL. This behaviour will not change for the
 * isoparametric elements which are currently implemented but may be
 * modified in future to take data from the element. 
 *
 * @{
 */

GHashTable *_bem3d_shfunc_funcs = NULL, *_bem3d_shfunc_names = NULL ;

/** 
 * Zero order (constant) shape function for a triangle, \f$L=1\f$.
 * 
 * @param s coordinate
 * @param t coordinate
 * @param L \f$L_{i}\f$, shape functions;
 * @param dLds derivatives of \f$L_{i}\f$ with respect to \f$s\f$;
 * @param dLdt derivatives of \f$L_{i}\f$ with respect to \f$t\f$;
 * @param data data to be passed to the shape function (ignored)
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_shfunc_t0(gdouble s, gdouble t, 
		     gdouble *L, gdouble *dLds,
		     gdouble *dLdt, gpointer data)

{
  if ( L != NULL ) L[0] = 1.0 ; 
  if ( dLds != NULL ) dLds[0] = 0.0 ; 
  if ( dLdt != NULL ) dLdt[0] = 0.0 ; 

  return BEM3D_SUCCESS ;
}

/** 
 * Linear shape function for a triangle \f$(0,0)\f$, \f$(1,0)\f$,
 * \f$(0,1)\f$.
 * 
 * @param s local coordinate
 * @param t local coordinate
 * @param L \f$L_{i}\f$, shape functions;
 * @param dLds derivatives of \f$L_{i}\f$ with respect to \f$s\f$;
 * @param dLdt derivatives of \f$L_{i}\f$ with respect to \f$t\f$;
 * @param data data to be passed to the shape function (ignored)
 * 
 * @return ::BEM3D_SUCCESS on success. 
 */

gint bem3d_shfunc_t1(gdouble s, gdouble t, 
		     gdouble *L, gdouble *dLds,
		     gdouble *dLdt, gpointer data)
     /*
       |
       1  3
       |  |\
       |  | \
       T  |  \
       |  |   \
       |  |    \
       0  1-----2
       |
       +--0--S--1-->
     */

{
  if ( L != NULL ) {
    L[0] = 1.0 - s - t ; L[1] = s ; L[2] = t ;
  }

  if ( dLds != NULL ) {
    dLds[0] = -1.0 ; dLds[1] = 1.0 ; dLds[2] = 0.0 ;
  }

  if ( dLdt != NULL ) {
    dLdt[0] = -1.0 ; dLdt[1] = 0.0 ; dLdt[2] = 1.0 ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Second order shape function for a triangle \f$(0,0)\f$, \f$(1,0)\f$,
 * \f$(0,1)\f$.
 * 
 * @param s local coordinate
 * @param t local coordinate
 * @param L \f$L_{i}\f$, shape functions;
 * @param dLds derivatives of \f$L_{i}\f$ with respect to \f$s\f$;
 * @param dLdt derivatives of \f$L_{i}\f$ with respect to \f$t\f$;
 * @param data data to be passed to the shape function (ignored)
 * 
 * @return ::BEM3D_SUCCESS on success. 
 */

gint bem3d_shfunc_t2(gdouble s, gdouble t, 
		     gdouble *L, gdouble *dLds,
		     gdouble *dLdt, gpointer data)
     /*
       |
       1  3
       |  |\
       |  | \
       T  6  5
       |  |   \
       |  |    \
       0  1--4--2
       |
       +--0--S--1-->
     */

{
  if ( L != NULL ) {
    L[0] = 2.0*(1.0-s-t)*(0.5-s-t) ;
    L[1] = 2.0*s*(s-0.5);
    L[2] = 2.0*t*(t-0.5);
    L[3] = 4.0*s*(1.0-s-t);
    L[4] = 4.0*s*t;
    L[5] = 4.0*t*(1.0-s-t);
/*     if ( shapefunction_isnan(L, 6) )  */
/*       g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, */
/* 	    "%s: NaN in dLdt at (%g,%g)", __FUNCTION__, s, t) ; */
  }

  if ( dLds != NULL ) {
    dLds[0] = -3.0 + 4.0*s + 4.0*t ;
    dLds[1] = -1.0 + 4.0*s;
    dLds[2] =  0.0;
    dLds[3] =  4.0 - 8.0*s - 4.0*t;
    dLds[4] =  4.0*t;
    dLds[5] =  -4.0*t;
/*     if ( shapefunction_isnan(dLds, 6) )  */
/*       g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, */
/* 	    "%s: NaN in dLds at (%g,%g)", __FUNCTION__, s, t) ;   */
  }

  if ( dLdt != NULL ) {
    dLdt[0] = -3.0 + 4.0*s + 4.0*t ;
    dLdt[1] = 0.0 ;
    dLdt[2] = -1.0 + 4.0*t ;
    dLdt[3] = -4.0*s ;
    dLdt[4] = 4.0*s ;
    dLdt[5] = 4.0 - 4.0*s - 8.0*t ;
/*     if ( shapefunction_isnan(dLds, 6) )  */
/*       g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, */
/* 	    "%s: NaN in dLdt at (%g,%g)", __FUNCTION__, s, t) ; */

  }

  return BEM3D_SUCCESS ;
}

/** 
 * Third order shape function for a triangle \f$(0,0)\f$, \f$(1,0)\f$,
 * \f$(0,1)\f$.
 * 
 * @param s local coordinate
 * @param t local coordinate
 * @param L \f$L_{i}\f$, shape functions;
 * @param dLds derivatives of \f$L_{i}\f$ with respect to \f$s\f$;
 * @param dLdt derivatives of \f$L_{i}\f$ with respect to \f$t\f$;
 * @param data data to be passed to the shape function (ignored)
 * 
 * @return ::BEM3D_SUCCESS on success. 
 */

gint bem3d_shfunc_t3(gdouble s, gdouble t, 
		     gdouble *L, gdouble *dLds,
		     gdouble *dLdt, gpointer data)

     /*
       |
       1  3
       |  |\
       |  | \
       |  8  7
       |  |   \
       T  |    \
       |  9  10 6
       |  |      \
       |  |       \
       0  1--4--5--2
       |
       +--0----S---1-->
     */
{
  gdouble stm[10] ;
  gint i, j ;

  if ( L != NULL ) {
    stm[0] = 1.0 ; stm[1] = s ; stm[2] = t ;
    stm[3] = s*s ; stm[4] = s*t ; stm[5] = t*t ;
    stm[6] = s*s*s ; stm[7] = s*s*t ; stm[8] = s*t*t ; stm[9] = t*t*t ;
    for ( i = 0 ; i < 10 ; i ++ ) {
      L[i] = 0.0 ;
      for ( j = 0 ; j < 10 ; j ++ ) L[i] += t3_shape_cffts[10*i+j]*stm[j] ;
    }
  }

  if ( dLds != NULL ) {
    stm[0] = 0.0 ; stm[1] = 1.0 ; stm[2] = 0.0 ;
    stm[3] = 2*s ; stm[4] = t ; stm[5] = 0.0 ;
    stm[6] = 3*s*s ; stm[7] = 2*s*t ; stm[8] = t*t ; stm[9] = 0.0 ;
    for ( i = 0 ; i < 10 ; i ++ ) {
      dLds[i] = 0.0 ;
      for ( j = 0 ; j < 10 ; j ++ ) dLds[i] += t3_shape_cffts[10*i+j]*stm[j] ;
    }
  }

  if ( dLdt != NULL ) {
    stm[0] = 0.0 ; stm[1] = 0.0 ; stm[2] = 1.0 ;
    stm[3] = 0.0 ; stm[4] = s ; stm[5] = 2.0*t ;
    stm[6] = 0.0 ; stm[7] = s*s ; stm[8] = 2*s*t ; stm[9] = 3*t*t ;
    for ( i = 0 ; i < 10 ; i ++ ) {
      dLdt[i] = 0.0 ;
      for ( j = 0 ; j < 10 ; j ++ ) dLdt[i] += t3_shape_cffts[10*i+j]*stm[j] ;
    }
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Linear shape function for a quadrilateral \f$(0,0)\f$, \f$(1,0)\f$,
 * \f$(1,1)\f$, \f$(0,1)\f$.
 * 
 * @param s local coordinate
 * @param t local coordinate
 * @param L \f$L_{i}\f$, shape functions;
 * @param dLds derivatives of \f$L_{i}\f$ with respect to \f$s\f$;
 * @param dLdt derivatives of \f$L_{i}\f$ with respect to \f$t\f$;
 * @param data data to be passed to the shape function (ignored)
 * 
 * @return ::BEM3D_SUCCESS on success. 
 */

gint bem3d_shfunc_q1(gdouble s, gdouble t, 
		     gdouble *L, gdouble *dLds,
		     gdouble *dLdt, gpointer data)

/*
       |
       1  4-------3
       |  |       |
       T  |       |
       |  |       |
       0  1-------2
       |
       +--0----S---1-->
  
 */

{
  if ( L != NULL ) {
    L[0] = (1.0-s)*(1.0-t) ;
    L[1] = s*(1.0-t) ;
    L[2] = s*t ;
    L[3] = (1.0-s)*t ;
  }

  if ( dLds != NULL ) {
    dLds[0] = -1.0 + t ;
    dLds[1] = 1.0 - t ;
    dLds[2] = t ;
    dLds[3] = -t ;
  }

  if ( dLdt != NULL ) {
    dLdt[0] = -1.0 + s ;
    dLdt[1] = -s ;
    dLdt[2] = s ;
    dLdt[3] = 1.0-s ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Second order shape function for a quadrilateral \f$(0,0)\f$,
 * \f$(1,0)\f$, \f$(1,1)\f$, \f$(0,1)\f$.
 * 
 * @param s local coordinate
 * @param t local coordinate
 * @param L \f$L_{i}\f$, shape functions;
 * @param dLds derivatives of \f$L_{i}\f$ with respect to \f$s\f$;
 * @param dLdt derivatives of \f$L_{i}\f$ with respect to \f$t\f$;
 * @param data data to be passed to the shape function (ignored)
 * 
 * @return ::BEM3D_SUCCESS on success. 
 */


gint bem3d_shfunc_q2(gdouble s, gdouble t, 
		     gdouble *L, gdouble *dLds,
		     gdouble *dLdt, gpointer data)

     /*
    |
    1  7--8--9
    |  |  :  |
    |  |  :  |
    S  4..5..6
    |  |  :  |
    |  |  :  |
    0  1--2--3
    |
    +--0--R--1-->
     */

     /*
    |
    1  4--7--3
    |  |  :  |
    |  |  :  |
    S  8..9..6
    |  |  :  |
    |  |  :  |
    0  1--5--2
    |
    +--0--R--1-->
     */

     /*
       From http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/
      */

{
  if ( L != NULL ) {
    L[0] = 0.25*(1.0-s)*(1.0-t)*s*t ;
    L[1] = -0.25*(1.0+s)*(1.0-t)*s*t ;
    L[2] = 0.25*(1.0+s)*(1.0+t)*s*t ;
    L[3] = -0.25*(1.0-s)*(1.0+t)*s*t ;
    L[4] = -0.5*(1.0-s*s)*(1-t)*t ;
    L[5] = 0.5*(1.0+s)*(1-t*t)*s ;
    L[6] = 0.5*(1.0-s*s)*(1+t)*t ;
    L[7] = -0.5*(1.0-s)*(1-t*t)*s ;
    L[8] = (1.0-s*s)*(1.0-t*t) ;
  }

  if ( dLds != NULL ) {
    dLds[0] = 0.25*(1.0-2*s)*(1.0-t)*t ;
    dLds[1] = -0.25*(1.0+2*s)*(1.0-t)*t ;
    dLds[2] = 0.25*(1.0+2*s)*(1.0+t)*t ;
    dLds[3] = -0.25*(1.0-2*s)*(1.0+t)*t ;
    dLds[4] = -0.5*(-2*s)*(1-t)*t ;
    dLds[5] = 0.5*(1.0+2*s)*(1-t*t) ;
    dLds[6] = 0.5*(-2*s)*(1+t)*t ;
    dLds[7] = -0.5*(1.0-2*s)*(1-t*t) ;
    dLds[8] = (-2*s)*(1.0-t*t) ;
  }

  if ( dLdt != NULL ) {
    dLdt[0] = 0.25*(1.0-s)*(1.0-2*t)*s ;
    dLdt[1] = -0.25*(1.0+s)*(1.0-2*t)*s ;
    dLdt[2] = 0.25*(1.0+s)*(1.0+2*t)*s ;
    dLdt[3] = -0.25*(1.0-s)*(1.0+2*t)*s ;
    dLdt[4] = -0.5*(1.0-s*s)*(1-2*t) ;
    dLdt[5] = 0.5*(1.0+s)*(-2*t)*s ;
    dLdt[6] = 0.5*(1.0-s*s)*(1+2*t) ;
    dLdt[7] = -0.5*(1.0-s)*(-2*t)*s ;
    dLdt[8] = (1.0-s*s)*(-2*t) ;
  }

  return BEM3D_SUCCESS ;
}

/** 
 * Add a ::BEM3DShapeFunc to the lookup table used in mapping shape
 * functions to names. 
 * 
 * @param func ::BEM3DShapeFunc to add
 * @param name name of shape function, which must be unique.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_shapefunc_lookup_add(BEM3DShapeFunc func, gchar *name)

{
  g_return_val_if_fail(func != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(name != NULL, BEM3D_NULL_ARGUMENT) ;

  if (_bem3d_shfunc_names == NULL ) 
    g_error("%s: shape function lookup table has not been initialized; "
	    "call bem3d_shapefunc_lookup_init()",
	    __FUNCTION__) ;

  if ( g_hash_table_lookup(_bem3d_shfunc_funcs, name) != NULL ) 
    g_error("%s: shape function %s already added", __FUNCTION__, name) ;

  g_hash_table_insert(_bem3d_shfunc_funcs, name, func) ;
  g_hash_table_insert(_bem3d_shfunc_names, func, name) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Initialize the shape function lookup table. This function must be
 * called before any shape functions are used or elements read. As
 * well as allocating the lookup table, it also adds the basic
 * isoparametric shape functions of order 0--3 for triangles and
 * orders 1 and 2 for quadrilaterals.
 * 
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_shapefunc_lookup_init()

{
  _bem3d_shfunc_funcs = g_hash_table_new(g_str_hash, g_str_equal) ;
  _bem3d_shfunc_names = g_hash_table_new(g_direct_hash, g_direct_equal) ;

  bem3d_shapefunc_lookup_add(bem3d_shfunc_t0, "SF_ISO_T0") ;
  bem3d_shapefunc_lookup_add(bem3d_shfunc_t1, "SF_ISO_T1") ;
  bem3d_shapefunc_lookup_add(bem3d_shfunc_t2, "SF_ISO_T2") ;
  bem3d_shapefunc_lookup_add(bem3d_shfunc_t3, "SF_ISO_T3") ;
  bem3d_shapefunc_lookup_add(bem3d_shfunc_q1, "SF_ISO_Q1") ;
  bem3d_shapefunc_lookup_add(bem3d_shfunc_q2, "SF_ISO_Q2") ;

  return BEM3D_SUCCESS ;
}

/** 
 * Look up a shape function, given its name. 
 * 
 * @param name the name under which the shape function has been
 * registered.
 * 
 * @return the ::BEM3DShapeFunc of that name. Otherwise, the program
 * reports an error and exits.
 */

BEM3DShapeFunc bem3d_shapefunc_lookup_func(const gchar *name)

{
  BEM3DShapeFunc func ;

  if ( _bem3d_shfunc_names == NULL ) 
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	  "%s: shape function lookup table has not been initialized; "
	    "call bem3d_shapefunc_lookup_init()",
	  __FUNCTION__) ;

  func = g_hash_table_lookup(_bem3d_shfunc_funcs, name) ;

  if ( func == NULL ) 
    g_error("%s: shape function %s not found", __FUNCTION__, name) ;

  return func ;
}

/** 
 * Find the name under which a ::BEM3DShapeFunc has been registered. 
 * 
 * @param func a ::BEM3DShapeFunc
 * 
 * @return the name of \a func in the lookup table, or NULL if not found. 
 */

const gchar *bem3d_shapefunc_lookup_name(BEM3DShapeFunc func)

{
  const gchar *name ;

  if ( _bem3d_shfunc_names == NULL ) 
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	  "%s: shape function lookup table has not been initialized; "
	    "call bem3d_shapefunc_lookup_init()",
	  __FUNCTION__) ;

  name = g_hash_table_lookup(_bem3d_shfunc_names, func) ;
  
  return name ;
}

/**
 * @}
 * 
 */
