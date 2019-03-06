/* transforms.c
 * 
 * Copyright (C) 2012, 2013 by Michael Carley
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

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <gts.h>

#include "polar.h"

gint triangle_axes(gdouble *x1, gdouble *x2, gdouble *x3,
		   gdouble *s, gdouble *t, gdouble *n) ;
gint triangle_angles(gdouble a, gdouble b, gdouble c,
		     gdouble *A, gdouble *B, gdouble *C) ;

gint triangle_orientations(gdouble *x1, gdouble *x2, gdouble *x3,
			   gdouble *p, gdouble *o12, gdouble *o23, 
			   gdouble *o31)

{
  *o12 = orient2d(p, x1, x2) ;
  *o12 = ( *o12 != 0 ? (*o12 < 0 ? -1 : 1) : 0 ) ;
  *o23 = orient2d(p, x2, x3) ;
  *o23 = ( *o23 != 0 ? (*o23 < 0 ? -1 : 1) : 0 ) ;
  *o31 = orient2d(p, x3, x1) ;
  *o31 = ( *o31 != 0 ? (*o31 < 0 ? -1 : 1) : 0 ) ;

  return 0 ;
}

void matrix_transform(gdouble *A, gdouble *x, gdouble *p, gdouble *y)

/*
  y = [A](p-x)
 */

{
  y[0] = A[0]*(p[0]-x[0]) + A[1]*(p[1]-x[1]) + A[2]*(p[2]-x[2]) ;
  y[1] = A[3]*(p[0]-x[0]) + A[4]*(p[1]-x[1]) + A[5]*(p[2]-x[2]) ;
  y[2] = A[6]*(p[0]-x[0]) + A[7]*(p[1]-x[1]) + A[8]*(p[2]-x[2]) ;

  return ;
}

gint triangle_rotation_matrix(gdouble *x1, gdouble *x2, gdouble *x3,
			      gdouble *A)

/*
  Input: general triangle coordinates (x1, x2, x3);
  Output: 
          transformation matrix A such that a point p in the original system 
	  transforms to [A](p-x1)
	  in the transformed system, x1 is transformed to the origin
 */

{

  A[0] = x2[0] - x1[0] ; A[1] = x2[1] - x1[1] ; A[2] = x2[2] - x1[2] ; 
  
  A[6] = A[1]*(x3[2] - x2[2]) - A[2]*(x3[1] - x2[1]) ;
  A[7] = A[2]*(x3[0] - x2[0]) - A[0]*(x3[2] - x2[2]) ;
  A[8] = A[0]*(x3[1] - x2[1]) - A[1]*(x3[0] - x2[0]) ;

  {
    gdouble _len = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]) ;
    A[0] /= _len ; A[1] /= _len ; A[2] /= _len ; 
    _len = sqrt(A[6]*A[6] + A[7]*A[7] + A[8]*A[8]) ;
    A[6] /= _len ; A[7] /= _len ; A[8] /= _len ; 
  }

  A[3] = A[7]*A[2] - A[8]*A[1] ;
  A[4] = A[8]*A[0] - A[6]*A[2] ;
  A[5] = A[6]*A[1] - A[7]*A[0] ;

  return 0 ;
}

gint triangle_axes(gdouble *x1, gdouble *x2, gdouble *x3,
		   gdouble *s, gdouble *t, gdouble *n)

/*
  generate orthonormal coordinate system based on plane of triangle
  (x1,x2,x3), with (s,t) in plane, and n normal.
 */

{
  gdouble len ;

  s[0] = x2[0] - x1[0] ; s[1] = x2[1] - x1[1] ; s[2] = x2[2] - x1[2] ;
  t[0] = x3[0] - x2[0] ; t[1] = x3[1] - x2[1] ; t[2] = x3[2] - x2[2] ;

  gts_vector_cross(n, s, t) ; len = gts_vector_norm(n) ;
  n[0] /= len ; n[1] /= len ; n[2] /= len ; 
  
  len = gts_vector_norm(s) ;
  s[0] /= len ; s[1] /= len ; s[2] /= len ; 
  gts_vector_cross(t, n, s) ;

  return 0 ;
}

gint triangle_angles(gdouble a, gdouble b, gdouble c,
		     gdouble *A, gdouble *B, gdouble *C)

{
  *A = acos(0.5*(b*b + c*c - a*a)/b/c) ;
  *B = acos(0.5*(c*c + a*a - b*b)/c/a) ;

  *C = M_PI - *A - *B ;

  return 0 ;
}
