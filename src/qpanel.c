#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <glib.h>

#include <gqr.h>

#include <gsl/gsl_poly.h>

#ifndef REAL_FLOAT
#ifndef REAL_DOUBLE
#define REAL_DOUBLE
#endif
#endif

#ifdef REAL_FLOAT
/* #warning "Compiling in single precision" */
#define REAL gfloat
#define REAL_ATAN2(x,y)  atan2f((x),(y))
#define REAL_ATAN(x)     atanf((x))
#define REAL_SIN(x)      sinf((x))
#define REAL_COS(x)      cosf((x))
#define REAL_LOG(x)      logf((x))
#define REAL_SQRT(x)     sqrtf((x))
#define QTRIANGLE_QUAD_RULE  qtri_quad_rule_f
#else
/* #warning "Compiling in double precision" */
#define REAL gdouble
#define REAL_ATAN2(x,y)      atan2((x),(y))
#define REAL_ATAN(x)         atan((x))
#define REAL_SIN(x)          sin((x))
#define REAL_COS(x)          cos((x))
#define REAL_LOG(x)          log((x))
#define REAL_SQRT(x)         sqrt((x))
#define QTRIANGLE_QUAD_RULE  qtri_quad_rule
#endif /*REAL*/

gint QTRIANGLE_QUAD_RULE(REAL x1[], REAL x2[], REAL x3[],
			 REAL x4[], REAL x5[], REAL x6[],
			 gint ntmin, gint ntmax, REAL dt,
			 gint nrmin, gint nrmax, REAL dr,
			 gqr_rule_t *gt, gqr_rule_t *gr,
			 gboolean hypersingular,
			 REAL *xi, gint xistr, 
			 REAL *eta, gint etastr, 
			 REAL *wt, gint wtstr,
			 gint *ngp) ;

REAL EDGE_SHAPES_T2[] = {1, -3, 2, 0, -1, 2, 0, 4, -4} ;
REAL EDGE_DIFFS_T2[] = {-3, 4, -1, 4, 4, -8} ;

/* *INDENT-OFF* */

REAL EDGE_SHAPE_DIFFS_T2[] =				\
  { 
    -3,   13,  -18,    8,  /* L1.L1' */
    -1,    7,  -14,    8,  /* L1.L2' */
    4,  -20,   32,  -16,  /* L1.L3' */
    0,    3,  -10,    8,  /* L2.L1' */
    0,    1,   -6,    8,  /* L2.L2' */
    0,   -4,   16,  -16,  /* L2.L2' */
    0,  -12,   28,  -16,  /* L3.L1' */
    0,   -4,   20,  -16,  /* L3.L2' */
    0,   16,  -48,   32,  /* L3.L3' */
  } ;

/* *INDENT-ON* */

#define edge_shape_func_t2(_i,_g)					\
  (EDGE_SHAPES_T2[3*(_i)+0]+					\
   (_g)*(EDGE_SHAPES_T2[3*(_i)+1] + (_g)*EDGE_SHAPES_T2[3*(_i)+2]))

#define edge_shape_diff_t2(_i,_g)				\
  (EDGE_DIFFS_T2[2*(_i)+0] + (_g)*EDGE_DIFFS_T2[2*(_i)+1])

#define edge_shape_cfft_t2(_i, _j) EDGE_SHAPES_T2[3*(_i)+(_j)]

#define INTER_WIDTH 8
#define INTER_PHI   0
#define INTER_XI    1
#define INTER_ETA   2
#define INTER_X     3
#define INTER_Y     4
#define INTER_Z     5
#define INTER_TAG   6

enum {
  INTER_CORNER,
  INTER_TANGENT,
  INTER_ON_EDGE
} ;

static gint compare_real(const void *i1, const void *i2)

{
  if ( *((REAL *)i1) < *((REAL *)i2) ) return -1 ;
  if ( *((REAL *)i1) > *((REAL *)i2) ) return  1 ;

  return 0 ;
}

static void area_coordinates_from_edge(gint e, REAL g, REAL *xi, REAL *eta)

{
  switch ( e ) {
  default: g_assert_not_reached() ; break ;
  case 1: *xi = g ; *eta = 0.0 ; break ;
  case 2: *xi = 1.0-g ; *eta = g ; break ;
  case 3: *xi = 0.0 ; *eta = 1.0-g ; break ;
  }

  return ;
}

static void edge_shape_funcs(REAL g, REAL L[])

{
  L[0] = edge_shape_func_t2(0,g) ;
  L[1] = edge_shape_func_t2(1,g) ;
  L[2] = edge_shape_func_t2(2,g) ;

  return ;
}

static void area_shape_diffs(REAL xi, REAL eta, 
			     REAL L[], REAL dLx[], REAL dLe[])

{
  L[0] = 2*(1-xi-eta)*(0.5-xi-eta) ;
  L[1] = 2*xi*(xi-0.5) ;
  L[2] = 2*eta*(eta-0.5) ;
  L[3] = 4*xi*(1-xi-eta) ;
  L[4] = 4*xi*eta ;
  L[5] = 4*eta*(1-xi-eta) ;

  dLx[0] = -2*(1.5-2*xi-2*eta) ;
  dLx[1] = 4*xi-1 ;
  dLx[2] = 0.0 ;
  dLx[3] = 4-8*xi-4*eta ;
  dLx[4] = 4*eta ;
  dLx[5] = -4*eta ;

  dLe[0] = -2*(1.5-2*xi-2*eta) ;
  dLe[1] = 0.0 ;
  dLe[2] = 4*eta-1 ;
  dLe[3] = -4*xi ;
  dLe[4] = 4*xi ;
  dLe[5] = 4-8*eta-4*xi ;

  return ;
}

static gint solve_quadratic(REAL a, REAL b, REAL c, REAL *x0, REAL *x1)

{
  return gsl_poly_solve_quadratic(a, b, c, x0, x1) ;
}

static gint solve_cubic(REAL a, REAL b, REAL c, REAL d,
			REAL *x0, REAL *x1, REAL *x2)

{
  *x0 = *x1 = *x2 = -G_MAXDOUBLE ;
  if ( a == 0.0 ) return solve_quadratic(b, c, d, x0, x1) ;

  return gsl_poly_solve_cubic(b/a, c/a, d/a, x0, x1, x2) ;
}

static void area_interpolate(REAL x1[], REAL x2[], REAL x3[],
			     REAL x4[], REAL x5[], REAL x6[],
			     REAL L[], REAL x[])

{
  x[0] = L[0]*x1[0]+L[1]*x2[0]+L[2]*x3[0]+L[3]*x4[0]+L[4]*x5[0]+L[5]*x6[0] ;
  x[1] = L[0]*x1[1]+L[1]*x2[1]+L[2]*x3[1]+L[3]*x4[1]+L[4]*x5[1]+L[5]*x6[1] ;
  x[2] = L[0]*x1[2]+L[1]*x2[2]+L[2]*x3[2]+L[3]*x4[2]+L[4]*x5[2]+L[5]*x6[2] ;

  return ;
}

static void edge_interpolate(REAL x1[], REAL x2[], REAL x3[], REAL L[],
			     REAL x[])

{
  x[0] = L[0]*x1[0] + L[1]*x2[0] + L[2]*x3[0] ;
  x[1] = L[0]*x1[1] + L[1]*x2[1] + L[2]*x3[1] ;
  x[2] = L[0]*x1[2] + L[1]*x2[2] + L[2]*x3[2] ;

  return ;
}

static gint edge_intersection(REAL x1[], REAL x2[], REAL x3[], REAL th, 
			      REAL g[])

/*
  Compute intersections of ray of angle th with curved edge x(123),
  return number of intersections on edge.
*/

{
  gint n = 0, nr ;
  REAL a, b, c, r0, r1, S, C ;

  S = REAL_SIN(th) ; C = REAL_COS(th) ;

  a = 
    edge_shape_cfft_t2(0,2)*(S*x1[0] - C*x1[1]) +
    edge_shape_cfft_t2(1,2)*(S*x2[0] - C*x2[1]) +
    edge_shape_cfft_t2(2,2)*(S*x3[0] - C*x3[1]) ;
  b = 
    edge_shape_cfft_t2(0,1)*(S*x1[0] - C*x1[1]) +
    edge_shape_cfft_t2(1,1)*(S*x2[0] - C*x2[1]) +
    edge_shape_cfft_t2(2,1)*(S*x3[0] - C*x3[1]) ;
  c = 
    edge_shape_cfft_t2(0,0)*(S*x1[0] - C*x1[1]) +
    edge_shape_cfft_t2(1,0)*(S*x2[0] - C*x2[1]) +
    edge_shape_cfft_t2(2,0)*(S*x3[0] - C*x3[1]) ;
    
  nr = solve_quadratic(a, b, c, &r0, &r1) ;

  if ( nr == 0 ) return 0 ;
  if ( nr == 1 ) {
    if ( r0 < 0 || r0 > 1 ) return 0 ;
    g[0] = r0 ; return 1 ;
  }

  if ( r0 > 0 && r0 < 1 ) { g[0] = r0 ; n = 1 ; }
  if ( r1 > 0 && r1 < 1 && r1 != r0 ) { g[n] = r1 ; n ++ ; }  

  return n ;
}

static gint edge_tangents(REAL x1[], REAL x2[], REAL x3[], REAL g[])

/*
  Compute points on edge x(012) which have a tangent passing through
  the origin.  
*/

{
  gint n = 0, nr ;
  REAL r0, r1, r2, q[4] ;
  
  for ( nr = 0 ; nr < 4 ; nr ++ ) {
    q[nr] = 
      x1[0]*x2[1]*(EDGE_SHAPE_DIFFS_T2[1*4+nr] - EDGE_SHAPE_DIFFS_T2[3*4+nr]) +
      x1[0]*x3[1]*(EDGE_SHAPE_DIFFS_T2[2*4+nr] - EDGE_SHAPE_DIFFS_T2[6*4+nr]) +
      x2[0]*x1[1]*(EDGE_SHAPE_DIFFS_T2[3*4+nr] - EDGE_SHAPE_DIFFS_T2[1*4+nr]) +
      x2[0]*x3[1]*(EDGE_SHAPE_DIFFS_T2[5*4+nr] - EDGE_SHAPE_DIFFS_T2[7*4+nr]) +
      x3[0]*x1[1]*(EDGE_SHAPE_DIFFS_T2[6*4+nr] - EDGE_SHAPE_DIFFS_T2[2*4+nr]) +
      x3[0]*x2[1]*(EDGE_SHAPE_DIFFS_T2[7*4+nr] - EDGE_SHAPE_DIFFS_T2[5*4+nr]) ;
  }

  if ( (nr = solve_cubic(q[3], q[2], q[1], q[0], &r0, &r1, &r2)) == 0 )
    return 0 ;
  
  if ( r0 > 0 && r0 < 1 ) { g[0] = r0 ; n = 1 ; }
  if ( r1 > 0 && r1 < 1 && r1 != r0 ) { g[n] = r1 ; n ++ ; }  
  if ( r2 > 0 && r2 < 1 && r2 != r1 ) { g[n] = r2 ; n ++ ; }  

  return n ;
}

static void tangent_on_edge(gint e, REAL gm, REAL x1[], REAL x2[], REAL x3[],
			    gint *ne, REAL inter[])

{
  REAL dL[3], x[3] ;
  
  dL[0] = edge_shape_diff_t2(0,gm) ;
  dL[1] = edge_shape_diff_t2(1,gm) ;
  dL[2] = edge_shape_diff_t2(2,gm) ;

  edge_interpolate(x1, x2, x3, dL, x) ;

  inter[(*ne)*INTER_WIDTH + INTER_ETA] = e ;
  inter[(*ne)*INTER_WIDTH + INTER_XI] = gm ;
  inter[(*ne)*INTER_WIDTH + INTER_TAG] = INTER_ON_EDGE ;
  inter[(*ne)*INTER_WIDTH + INTER_PHI] = atan2(x[1], x[0]) ;

  inter[(*ne)*INTER_WIDTH + INTER_X] = 
    inter[(*ne)*INTER_WIDTH + INTER_Y] = 
    inter[(*ne)*INTER_WIDTH + INTER_Z] = 0.0 ;
  if ( gm == 1 ) {
    if ( inter[(*ne)*INTER_WIDTH + INTER_PHI] < 0 ) 
      inter[(*ne)*INTER_WIDTH + INTER_PHI] += M_PI ;
    else
      inter[(*ne)*INTER_WIDTH + INTER_PHI] -= M_PI ;
  }

  (*ne) ++ ;

  if ( gm != 0.5 ) return ;

  inter[(*ne)*INTER_WIDTH + INTER_ETA] = e ;
  inter[(*ne)*INTER_WIDTH + INTER_XI] = gm ;
  inter[(*ne)*INTER_WIDTH + INTER_TAG] = INTER_ON_EDGE ;
  if ( inter[((*ne)-1)*INTER_WIDTH + INTER_PHI] < M_PI ) 
    inter[(*ne)*INTER_WIDTH + INTER_PHI] = 
      inter[((*ne)-1)*INTER_WIDTH + INTER_PHI] + M_PI ;
  else
    inter[(*ne)*INTER_WIDTH + INTER_PHI] = 
      inter[((*ne)-1)*INTER_WIDTH + INTER_PHI] - M_PI ;    

  inter[(*ne)*INTER_WIDTH + INTER_X] = 
    inter[(*ne)*INTER_WIDTH + INTER_Y] = 
    inter[(*ne)*INTER_WIDTH + INTER_Z] = 0.0 ;

  (*ne) ++ ;
  
  return ;
}

static void edge_intersect(gint e, REAL x1[], REAL x2[], REAL x3[],
			   gint *ne, REAL inter[])

{
  gint i ;
  REAL g[4], L[3] ;

  if ( x1[0] == 0 && x1[1] == 0 ) {
    tangent_on_edge(e, 0, x1, x2, x3, ne, inter) ; return ;
  }
  if ( x2[0] == 0 && x2[1] == 0 ) {
    tangent_on_edge(e, 1, x1, x2, x3, ne, inter) ; return ;
  }
  if ( x3[0] == 0 && x3[1] == 0 ) {
    tangent_on_edge(e, 0.5, x1, x2, x3, ne, inter) ; return ;
  }

  i = edge_tangents(x1, x2, x3, g) - 1 ;
  for ( ; i >= 0 ; i -- ) {
    inter[(*ne)*INTER_WIDTH + INTER_ETA] = e ; 
    inter[(*ne)*INTER_WIDTH + INTER_XI]  = g[i] ;
    edge_shape_funcs(g[i], L) ;
    edge_interpolate(x1, x2, x3, L, &(inter[(*ne)*INTER_WIDTH + INTER_X])) ;
    inter[(*ne)*INTER_WIDTH + INTER_TAG]  = INTER_TANGENT ;
    (*ne) ++ ;
  }
  
  return ;
}

static void sub_intersections(REAL x1[], REAL x2[], REAL x3[], gint e,
			      REAL th, REAL r[], gint *nr)

{
  REAL x[3], g[4], L[3], sc ;
  gint i, n, si = 1 ;

  if ( (n = edge_intersection(x1, x2, x3, th, g)) == 0) return ;

  /*for checking radius is positive*/
  if ( fabs(sc = REAL_SIN(th)) < 0.5 ) { sc = REAL_COS(th) ; si = 0 ; }

  for ( i = 0 ; i < n ; i ++ ) {
    edge_shape_funcs(g[i], L) ;
    edge_interpolate(x1, x2, x3, L, x) ;
    if ( x[si]/sc > 0.0 ) {      
      r[4*(*nr)+0] = REAL_SQRT(x[0]*x[0] + x[1]*x[1]) ;
      r[4*(*nr)+1] = e ;
      area_coordinates_from_edge(e, g[i], &(r[4*(*nr)+2]), &(r[4*(*nr)+3])) ;
      (*nr) ++ ;
    }
  }

  return ;
}

static gint tri_intersections(REAL x1[], REAL x2[], REAL x3[],
			      REAL x4[], REAL x5[], REAL x6[],
			      REAL th, REAL *r)
{
  gint nr = 0 ;

  sub_intersections(x1, x2, x4, 1, th, r, &nr) ;
  sub_intersections(x2, x3, x5, 2, th, r, &nr) ;
  sub_intersections(x3, x1, x6, 3, th, r, &nr) ;

  qsort((gpointer)r, nr, 4*sizeof(REAL), compare_real) ;

  return nr ;
}

static void shape_coordinates(REAL x1[], REAL x2[], REAL x3[],
			      REAL x4[], REAL x5[], REAL x6[],
			      REAL x, REAL y, REAL *xi, REAL *eta)

{
  REAL tol = 1e-7, L[6], dLx[6], dLe[6], R2, x0[3], dxx[3], dxe[3] ;
  REAL det ;

  area_shape_diffs(*xi, *eta, L, dLx, dLe) ;
  area_interpolate(x1, x2, x3, x4, x5, x6, L, x0) ;

  do {
    area_interpolate(x1, x2, x3, x4, x5, x6, dLx, dxx) ;
    area_interpolate(x1, x2, x3, x4, x5, x6, dLe, dxe) ;
    /*inverse of Jacobian matrix*/
    det  =  dxx[0]*dxe[1] - dxx[1]*dxe[0] ;

    *xi  += (dxe[1]*(x-x0[0]) - dxe[0]*(y-x0[1]))/det ;
    *eta -= (dxx[1]*(x-x0[0]) - dxx[0]*(y-x0[1]))/det ;
    
    area_shape_diffs(*xi, *eta, L, dLx, dLe) ;
    area_interpolate(x1, x2, x3, x4, x5, x6, L, x0) ;

    R2 = (x0[0] - x)*(x0[0] - x) + (x0[1] - y)*(x0[1] - y) ;
  } while ( R2 > tol*tol ) ;
  
  return ;
}

static gint rule_select_length(REAL x0, REAL x1, REAL dx, 
			       gint nmin, gint nmax)

{
  REAL n = (x1-x0)/dx ;
  if ( n > nmax ) return nmax ;
  if ( n < nmin ) return nmin ;

  return (gint)(4*round(0.25*n)) ;  
}

static void clean_intersections(REAL *inter, gint *ne)

/* check for duplicated angles, mod 2\pi */

{
  gint i ;

  for ( i = 0 ; i < *ne ; i ++ ) {
    while ( inter[i*INTER_WIDTH + INTER_PHI] > M_PI ) 
      inter[i*INTER_WIDTH + INTER_PHI] -= 2*M_PI ;
    while ( inter[i*INTER_WIDTH + INTER_PHI] < -M_PI ) 
      inter[i*INTER_WIDTH + INTER_PHI] += 2*M_PI ;
  }

  return ;
}

static gint inside_triangle(REAL x1[], REAL x2[], REAL x3[],
			    REAL x4[], REAL x5[], REAL x6[])
/*
  check if the curved triangle encloses the origin. Return 1 if origin
  is inside triangle, or on the boundary, 0 otherwise.
 */


{
  REAL xi, eta ;

  if ( x1[0]<0 && x2[0]<0 && x3[0]<0 && x4[0]<0 && x5[0]<0 && x6[0]< 0 )
    return 0 ;
  if ( x1[0]>0 && x2[0]>0 && x3[0]>0 && x4[0]>0 && x5[0]>0 && x6[0]> 0 )
    return 0 ;
  if ( x1[1]<0 && x2[1]<0 && x3[1]<0 && x4[1]<0 && x5[1]<0 && x6[1]< 0 )
    return 0 ;
  if ( x1[1]>0 && x2[1]>0 && x3[1]>0 && x4[1]>0 && x5[1]>0 && x6[1]> 0 )
    return 0 ;

  /*maybe we are on a vertex*/
  if ( x1[0] == 0 && x1[1] == 0 ) return 1 ;
  if ( x2[0] == 0 && x2[1] == 0 ) return 1 ;
  if ( x3[0] == 0 && x3[1] == 0 ) return 1 ;
  if ( x4[0] == 0 && x4[1] == 0 ) return 1 ;
  if ( x5[0] == 0 && x5[1] == 0 ) return 1 ;
  if ( x6[0] == 0 && x6[1] == 0 ) return 1 ;

  xi = eta = 0.0 ;
  shape_coordinates(x1, x2, x3, x4, x5, x6, 0, 0, &xi, &eta) ;
  if ( fabs(xi) < 1e-9 ) return 1 ;
  if ( fabs(eta) < 1e-9 ) return 1 ;
  if ( fabs(xi + eta -1.0) < 1e-6 ) return 1 ;

  if ( xi < 0 || xi > 1 )   return 0 ;
  if ( eta < 0 || eta > 1 ) return 0 ;
  if ( eta > 1.0-xi )       return 0 ;
  
  return 1 ;
}

gint QTRIANGLE_QUAD_RULE(REAL x1[], REAL x2[], REAL x3[],
			 REAL x4[], REAL x5[], REAL x6[],
			 gint ntmin, gint ntmax, REAL dt,
			 gint nrmin, gint nrmax, REAL dr,
			 gqr_rule_t *gt, gqr_rule_t *gr,
			 gboolean hypersingular,			 
			 REAL *xi, gint xistr, 
			 REAL *eta, gint etastr, 
			 REAL *wt, gint wtstr,
			 gint *ngp)

/*
  Generate a quadrature rule for integration on a second-order
  triangle with nodes x1, x2, x3, x4, x5, x6, and field point at
  (0,0,z). Quadrature points in angle are spaced approximately dt
  apart, with a minimum of ntmin and a maximum of ntmax points, and
  likewise for points in radius (dr, nrmin, nrmax). Quadrature rules
  are generated by GQR in the rules gt and gr. The output quadrature
  rule is given as ngp coordinates (xi,eta) with weights wt, inserted
  into the corresponding arrays with stride xistr, etastr, wtstr.
 */

{
  gint ne, nr, i, j, k, d, ngt, ngr ;
  REAL inter[256], t0, t1, r[128], *intr = &(r[0]), y[3], dyx[3], dye[3] ;
  REAL th, thbar, dth, s, sbar, ds, wth, *xi1, *eta1, S, C ;
  REAL K, L[6], dLx[6], dLe[6], wsum ;
  gqr_t rule ;
  gqr_parameter_t pgr ;

  /*adjust calculation of radial limits of integration*/
  if ( ( d = inside_triangle(x1, x2, x3, x4, x5, x6) ) == 1 ) { 
    r[0] = r[1] = r[2] = r[3] = 0.0 ; intr = &(r[4]) ; 
  }

  /*
    inter holds all edge intersection data, see #define INTER_... at
    top of file
  */

  /*corners first*/
  inter[0*INTER_WIDTH + INTER_ETA] = 1 ; 
  inter[0*INTER_WIDTH + INTER_XI]  = 0.0 ;
  memcpy(&(inter[0*INTER_WIDTH+INTER_X]), x1, 3*sizeof(REAL)) ;
  inter[0*INTER_WIDTH + INTER_TAG] = INTER_CORNER ;
  inter[1*INTER_WIDTH + INTER_ETA] = 2 ; 
  inter[1*INTER_WIDTH + INTER_XI]  = 0.0 ;
  memcpy(&(inter[1*INTER_WIDTH+INTER_X]), x2, 3*sizeof(REAL)) ;
  inter[1*INTER_WIDTH + INTER_TAG] = INTER_CORNER ;
  inter[2*INTER_WIDTH + INTER_ETA] = 3 ; 
  inter[2*INTER_WIDTH + INTER_XI]  = 0.0 ;
  memcpy(&(inter[2*INTER_WIDTH+INTER_X]), x3, 3*sizeof(REAL)) ;
  inter[2*INTER_WIDTH + INTER_TAG] = INTER_CORNER ;
  ne = 3 ;

  /*tangents*/
  edge_intersect(1, x1, x2, x4, &ne, inter) ;
  edge_intersect(2, x2, x3, x5, &ne, inter) ;
  edge_intersect(3, x3, x1, x6, &ne, inter) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    if ( inter[i*INTER_WIDTH + INTER_TAG] != INTER_ON_EDGE )
      inter[i*INTER_WIDTH + INTER_PHI] = atan2(inter[i*INTER_WIDTH + INTER_Y],
					       inter[i*INTER_WIDTH + INTER_X]) ;
  }

  clean_intersections(inter, &ne) ;

  qsort(inter, ne, INTER_WIDTH*sizeof(REAL), compare_real) ;

  memcpy(&(inter[ne*INTER_WIDTH + INTER_PHI]), inter, 
	 INTER_WIDTH*sizeof(REAL)) ;
  inter[ne*INTER_WIDTH + INTER_PHI] =
    inter[0*INTER_WIDTH + INTER_PHI] + 2.0*M_PI ;
  ne ++ ;

  if ( hypersingular ) {
    rule = GQR_GAUSS_LEGENDRE | GQR_GAUSS_HYPERSINGULAR ;
    gqr_parameter_clear(&pgr) ;
    gqr_parameter_set_int(&pgr, 32) ;
    gqr_parameter_set_double(&pgr, -1.0) ;
    gqr_rule_select(gr, rule, 128, &pgr) ;    
  } else 
    rule = GQR_GAUSS_LEGENDRE ;

  *ngp = 0 ; wsum = 0.0 ;
  for ( i = 0 ; i < ne - 1 ; i ++ ) {
    t0 = inter[i*INTER_WIDTH + INTER_PHI] ;
    t1 = inter[(i+1)*INTER_WIDTH + INTER_PHI] ;

    ngt = rule_select_length(t0, t1, dt, ntmin, ntmax) ;
    gqr_rule_select(gt, GQR_GAUSS_LEGENDRE, ngt, NULL) ;
    gqr_rule_scale(gt, t0, t1, &thbar, &dth) ;

    if ( dth > 1e-12 ) {
      for ( j = 0 ; j < gqr_rule_length(gt) ; j ++ ) {
	th = gqr_rule_abscissa(gt, j)*dth + thbar ;
	wth = gqr_rule_weight(gt, j)*dth ;
	S = REAL_SIN(th) ; C = REAL_COS(th) ;

	nr = tri_intersections(x1, x2, x3, x4, x5, x6, th, intr) + d - 1 ;

	for ( ; nr > 0 ; nr -= 2 ) {
	  ngr = rule_select_length(r[4*(nr-1)], r[4*nr], dr, nrmin, nrmax) ;

	  if ( !hypersingular || r[4*(nr-1)] != 0)
	    gqr_rule_select(gr, GQR_GAUSS_LEGENDRE, ngr, &pgr) ;	    
	  gqr_rule_scale(gr, r[4*(nr-1)], r[4*nr], &sbar, &ds) ;

	  if ( ds > 1e-6 ) {
	    xi1 = &(xi[xistr*(*ngp)]) ; eta1 = &(eta[etastr*(*ngp)]) ;
	    *xi1 = r[4*(nr-1)+2] ; *eta1 = r[4*(nr-1)+3] ;
	    for ( k = 0 ; k < gqr_rule_length(gr) ; k ++ ) {
	      xi1 = &(xi[xistr*(*ngp)]) ; eta1 = &(eta[etastr*(*ngp)]) ;
	      s = gqr_rule_abscissa(gr, k)*ds + sbar ;
	      shape_coordinates(x1, x2, x3, x4, x5, x6, s*C, s*S, xi1, eta1) ;

	      area_shape_diffs(*xi1, *eta1, L, dLx, dLe) ;
	      area_interpolate(x1, x2, x3, x4, x5, x6, L, y) ;
	      area_interpolate(x1, x2, x3, x4, x5, x6, dLx, dyx) ;
	      area_interpolate(x1, x2, x3, x4, x5, x6, dLe, dye) ;
	  
	      K = dyx[0]*dye[1] - dyx[1]*dye[0] ;
	  	      
	      wsum += (wt[wtstr*(*ngp)] = s*gqr_rule_weight(gr,k)*ds*wth/K) ;
	      (*ngp) ++ ;
	      xi[xistr*(*ngp)] = *xi1 ; eta[etastr*(*ngp)] = *eta1 ;
	    }
	  }
	}
      }
    }
    g_debug("%s: sum of %d quadrature weights == %1.16e", 
	    __FUNCTION__, *ngp, wsum) ;  
  }

  if ( fabs(wsum - 0.5) > 0.05 ) {
    fprintf(stderr, "%s: sum of weights=%e\n", __FUNCTION__, wsum) ;
    return -1 ;
  }

  return 0 ;
}
