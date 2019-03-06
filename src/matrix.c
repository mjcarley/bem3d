/* matrix.c
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

#include <glib.h>
#include <gts.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>

#include "bem3d.h"
#include "bem3d-private.h"

/**
 * @defgroup matrix Matrix operations
 *
 * Some basic matrix definitions and operations.
 *
 * @{
 * 
 */

#define det2x2(a,b,c,d) ((a)*(d)-(b)*(c))

/** 
 * Multiply a ::BEM3DMatrix by a GtsVector.
 * 
 * @param m ::BEM3DMatrix
 * @param v GtsVector to multiply
 * @param w GtsVector for result of operation
 * 
 * @return 0 on success.
 */

gint bem3d_matrix_vector_mul(BEM3DMatrix m, GtsVector v, GtsVector w)

{
  w[0] = m[0]*v[0] + m[1]*v[1] + m[2]*v[2] ; 
  w[1] = m[3]*v[0] + m[4]*v[1] + m[5]*v[2] ; 
  w[2] = m[6]*v[0] + m[7]*v[1] + m[8]*v[2] ; 

  return BEM3D_SUCCESS ;
}

/** 
 * Evaluate the determinant of a ::BEM3DMatrix
 * 
 * @param m ::BEM3DMatrix
 * 
 * @return determinant of \a m
 */

gdouble bem3d_matrix_det(BEM3DMatrix m)

{
  return (m[0]*det2x2(m[4], m[5], m[7], m[8])
	  - m[1]*det2x2(m[3], m[5], m[6], m[8])
	  + m[2]*det2x2(m[3], m[4], m[6], m[7])) ;
}

/** 
 * Find the inverse of a ::BEM3DMatrix. 
 * 
 * @param m ::BEM3DMatrix
 * @param im ::BEM3DMatrix to hold inverse of \a m
 * 
 * @return ::BEM3D_SUCCESS on success or ::BEM3D_FAILURE if \a m is
 * singular.
 */

gint bem3d_matrix_inverse(BEM3DMatrix m, BEM3DMatrix im)

{
  gdouble det ;

  det = bem3d_matrix_det(m) ;

  if ( det == 0.0 ) return BEM3D_FAILURE ;

  im[0] = (m[4]*m[8] - m[5]*m[7])/det; 
  im[1] = (m[7]*m[2] - m[1]*m[8])/det;
  im[2] = (m[1]*m[5] - m[4]*m[2])/det; 
  im[3] = (m[5]*m[6] - m[3]*m[8])/det; 
  im[4] = (m[0]*m[8] - m[6]*m[2])/det; 
  im[5] = (m[3]*m[2] - m[0]*m[5])/det; 
  im[6] = (m[3]*m[7] - m[6]*m[4])/det; 
  im[7] = (m[6]*m[1] - m[0]*m[7])/det; 
  im[8] = (m[0]*m[4] - m[1]*m[3])/det; 

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */

