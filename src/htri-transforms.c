#ifdef  _USE_SINCOS_
#define _GNU_SOURCE
#endif /*_USE_SINCOS_*/

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <glib.h>

#include "htriquad.h"

/*this is Shewchuk's adaptive predicate, from
  https://www.cs.cmu.edu/~quake/robust.html
 */
gdouble orient2d(gdouble *pa, gdouble *pb, gdouble *pc) ;

gint htri_triangle_decompose(gdouble *x1, gdouble *x2, gdouble *x3,
			     gdouble *x,
			     gdouble psi[], gdouble th[],
			     gdouble r1[], gdouble r2[],
			     gint sgn[],
			     gint *ntri)

/*
  Decomposition of triangle (x1,x2,x3) in plane z=0 into up to three
  oriented triangles (r1,r2,th) with rotation angle psi and common
  vertex x

  on output ntri is set to the number of triangles in the
  decomposition

  return: 0 if x lies on an edge of the triangle, -1 if outside, 1 if
  inside;

*/

{
  gdouble s1, s2, s3, l2, D[3], d[3], y1[2], y2[2], y3[2] ;
  gint in ;

  *ntri = 0 ;

  D[0] = orient2d(x1, x2, x) ;
  D[1] = orient2d(x2, x3, x) ;
  D[2] = orient2d(x3, x1, x) ;

  if ( D[0] == 0 || D[1] == 0 || D[2] == 0 ) in = 0 ;
  else {
    if ( D[0] < 0 || D[1] < 0 || D[2] < 0 ) in = -1 ;
    else in = 1 ;
  }

  if ( D[0] > 0.0 ) {
    y1[0] = x1[0] ; y1[1] = x1[1] ; d[0] = D[0] ; 
    y2[0] = x2[0] ; y2[1] = x2[1] ; d[1] = D[1] ;
    y3[0] = x3[0] ; y3[1] = x3[1] ; d[2] = D[2] ;
  } else {
    if ( D[1] > 0.0 ) {
      y1[0] = x2[0] ; y1[1] = x2[1] ; d[0] = D[1] ; 
      y2[0] = x3[0] ; y2[1] = x3[1] ; d[1] = D[2] ;
      y3[0] = x1[0] ; y3[1] = x1[1] ; d[2] = D[0] ;      
    } else {
      if ( D[2] > 0.0 ) {
	y1[0] = x3[0] ; y1[1] = x3[1] ; d[0] = D[2] ; 
	y2[0] = x1[0] ; y2[1] = x1[1] ; d[1] = D[0] ;
	y3[0] = x2[0] ; y3[1] = x2[1] ; d[2] = D[1] ;      
      } else 
	g_assert_not_reached() ;
    }
  }

  g_assert(d[0] > 0.0) ;

  s1 = sqrt((y1[0]-x[0])*(y1[0]-x[0]) + (y1[1]-x[1])*(y1[1]-x[1])) ;
  s2 = sqrt((y2[0]-x[0])*(y2[0]-x[0]) + (y2[1]-x[1])*(y2[1]-x[1])) ;
  s3 = sqrt((y3[0]-x[0])*(y3[0]-x[0]) + (y3[1]-x[1])*(y3[1]-x[1])) ;

  l2 = (y2[0] - y1[0])*(y2[0] - y1[0]) + (y2[1] - y1[1])*(y2[1] - y1[1]) ;
  th[(*ntri)] = acos(0.5*(s1*s1 + s2*s2 - l2)/s1/s2) ;
  if ( isnan(th[*ntri]) ) th[*ntri] = 0.0 ;
  
  r1[(*ntri)] = s1 ; r2[(*ntri)] = s2 ; 
  psi[(*ntri)] = atan2(y1[1] - x[1], y1[0] - x[0]) ;
  sgn[(*ntri)] = SIGN(d[0]) ;
  (*ntri) ++ ;

  g_assert(!isnan(th[0])) ;
  
  if ( d[1] != 0.0 ) {
    l2 = (y3[0] - y2[0])*(y3[0] - y2[0]) + (y3[1] - y2[1])*(y3[1] - y2[1]) ;
    th[(*ntri)] = acos(0.5*(s2*s2 + s3*s3 - l2)/s2/s3) ;
    if ( isnan(th[*ntri]) ) th[*ntri] = 0.0 ;

    /* if ( !isnan(th[*ntri]) ) { */
    if ( d[1] > 0 ) {
      r1[(*ntri)] = s2 ; r2[(*ntri)] = s3 ;
      psi[(*ntri)] = atan2(y2[1] - x[1], y2[0] - x[0]) ;
    } else {
      r1[(*ntri)] = s3 ; r2[(*ntri)] = s2 ;
      psi[(*ntri)] = atan2(y3[1] - x[1], y3[0] - x[0]) ;      
    }
    sgn[(*ntri)] = SIGN(d[1]) ;
    (*ntri) ++ ;
    /* } */
  }

  if ( d[2] != 0.0 ) {
    l2 = (y1[0] - y3[0])*(y1[0] - y3[0]) + (y1[1] - y3[1])*(y1[1] - y3[1]) ;
    th[(*ntri)] = acos(0.5*(s3*s3 + s1*s1 - l2)/s3/s1) ;
    if ( isnan(th[*ntri]) ) th[*ntri] = 0.0 ;

    /* if ( !isnan(th[*ntri]) ) { */
    if ( d[2] > 0 ) {
      r1[(*ntri)] = s3 ; r2[(*ntri)] = s1 ; 
      psi[(*ntri)] = atan2(y3[1] - x[1], y3[0] - x[0]) ;
    } else {
      r1[(*ntri)] = s1 ; r2[(*ntri)] = s3 ; 
      psi[(*ntri)] = atan2(y1[1] - x[1], y1[0] - x[0]) ;
    }
    sgn[(*ntri)] = SIGN(d[2]) ;
    (*ntri) ++ ;
    /* } */
  }

  return in ;
}

gint htri_triangle_axes(gdouble *x1, gdouble *x2, gdouble *x3,
			gdouble *x,
			gdouble *og, gdouble *s, gdouble *t, gdouble *n)

/*
  projection quantities for in-plane coordinate system based on
  triangle (x1,x2,x3) and field point x

  on exit: og is projection of x onto triangle plane, (s,t,n) is
  right-handed unit coordinate system in plane, with n unit normal
*/

{
  gdouble len, z ;

  s[0] = x2[0] - x1[0] ; s[1] = x2[1] - x1[1] ; s[2] = x2[2] - x1[2] ; 
  t[0] = x3[0] - x2[0] ; t[1] = x3[1] - x2[1] ; t[2] = x3[2] - x2[2] ; 

  _htri_vector_cross(n, s, t) ;
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  n[0] /= len ; n[1] /= len ; n[2] /= len ; 

  z = n[0]*(x[0] - x1[0]) + n[1]*(x[1] - x1[1]) + n[2]*(x[2] - x1[2])  ;

  og[0] = x[0] - n[0]*z ; og[1] = x[1] - n[1]*z ; og[2] = x[2] - n[2]*z ;

  len = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]) ;
  s[0] /= len ; s[1] /= len ; s[2] /= len ; 

  _htri_vector_cross(t, n, s) ;
  /* len = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]) ; */
  /* t[0] /= len ; t[1] /= len ; t[2] /= len ;  */

  return 0 ;
}

gint htri_triangle_project(gdouble *og, gdouble *s, gdouble *t, gdouble *n,
			   gdouble *x, gdouble *y)

/*
  compute the coordinates of a point in the coordinate system found
  using htri_triangle_axes

  on output y is set such that x = og + y[0]*s + y[1]*t + y[2]*n
*/

{
  y[0] = (x[0] - og[0])*s[0] + (x[1] - og[1])*s[1] + (x[2] - og[2])*s[2] ;
  y[1] = (x[0] - og[0])*t[0] + (x[1] - og[1])*t[1] + (x[2] - og[2])*t[2] ;
  y[2] = (x[0] - og[0])*n[0] + (x[1] - og[1])*n[1] + (x[2] - og[2])*n[2] ;

  return 0 ;
}
