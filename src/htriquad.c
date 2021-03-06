/* #define _USE_SINCOS_ */

#ifdef  _USE_SINCOS_
#define _GNU_SOURCE
#endif /*_USE_SINCOS_*/

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <glib.h>

#include "htriquad.h"
#include "binomials.h"

#include "trace.h"

static gdouble polyval2(gdouble *c, gint cstr,
			gdouble *d, gint dstr,
			gint nc, gdouble x)
/*
  evaluate polynomial \sum c_{i}d_{i}x^{i}

  c_i = c[i*cstr], d_i = d[i*dstr]

  strides can be negative but note no check is performed on array
  bounds
*/

{
  gdouble p = c[nc*cstr]*d[nc*dstr] ;
  gint i ;

  for ( i = nc-1 ; i >= 0 ; i -- ) p = p*x + c[i*cstr]*d[i*dstr] ;

  return p ;
}

gint quad_delta_C_n(gdouble al, gdouble t0, gdouble t1, gint N,
		    gdouble *Q)

/*
  integrals \int_{\theta_0}^{\theta_1} (\Delta/\cos\theta)^n d \theta,
  n = -3, ..., N
*/

{
  gdouble al2, ad, ad2, u0, u1, d0, d1, D0, D1, S0, S1, I0, I1, A0, A1 ;
  /* gdouble S[2], u[2] ; */
  gdouble d1n, d0n ;
  gint n, off ;

  off = 3 ;
  al2 = al*al ; ad2 = 1.0 - al2 ; ad = sqrt(ad2) ;

#ifdef _USE_SINCOS_
  sincos(t0, &S0, &u0) ; u0 = S0/u0 ;
  sincos(t1, &S1, &u1) ; u1 = S1/u1 ;
  /* S0 = sin(t0) ; S1 = sin(t1) ; */
  /* u0 = tan(t0) ; u1 = tan(t1) ;  */
#else
  S0 = sin(t0) ; S1 = sin(t1) ;
  u0 = tan(t0) ; u1 = tan(t1) ; 
#endif /*_USE_SINCOS_*/

  /* S[0] = S0 ; S[1] = S1 ; */
  /* u[0] = u0 ; u[1] = u1 ; */

  d0 = sqrt(1.0 + ad2*u0*u0) ; d1 = sqrt(1.0 + ad2*u1*u1) ;

  I0 = (atan(ad*u1) - atan(ad*u0))/ad ;
  I1 = log((d1 + ad*u1)/(d0 + ad*u0))/ad ;

  if ( al != 0.0 ) {
    A0 = asin(al*S0) ; A1 = asin(al*S1) ;
    D0 = sqrt(1.0 - al2*S0*S0) ; D1 = sqrt(1.0 - al2*S1*S1) ; 

    Q[off-3] = -ad2*(S1/D1 - S0/D0)/al2 + (A1 - A0)/al2/al ;
    Q[off-2] = (t1 - t0 - ad2*I0)/al2 ;
    Q[off-1] = (A1 - A0)/al ;
  } else {
    Q[off-3] = S1*(1-S1*S1/3.0) - S0*(1-S0*S0/3.0) ;
    Q[off-2] = 0.5*(S1*cos(t1) - S0*cos(t0)) + 0.5*(t1 - t0) ;
    Q[off-1] = S1 - S0 ;
  }

  d1n = d0n = 1.0 ;
  for ( n = 0 ; n <= N ; n += 2 ) {
    Q[off+n] = al2*Q[off+n-2] + ad2*I0 ;
    I0 = (u1*d1n - u0*d0n + I0*n)/(n+1) ;
    Q[off+n+1] = al2*Q[off+n+1-2] + ad2*I1 ;
    d1n *= d1 ; d0n *= d0 ;
    I1 = (u1*d1n - u0*d0n + I1*(n+1))/(n+1+1) ;
    d1n *= d1 ; d0n *= d0 ;
  }

  return 0 ;
}

gint quad_delta_C_n_T(gdouble al, gdouble t0, gdouble t1, gint N,
		      gdouble *Q)

/*
  integrals 
  \int_{\theta_0}^{\theta_1} (\Delta/\cos\theta)^n \tan\theta d \theta,
  n = -3, ..., N
*/

{
  gdouble al2, D0, D1 ;
  gdouble C0, C1, S0, S1 ;
  gdouble d1n, d0n ;
  gint n, off ;

  off = 3 ;
  al2 = al*al ;

#ifdef _USE_SINCOS_
  sincos(t0, &S0, &C0) ;
  sincos(t1, &S1, &C1) ;
#else
  S0 = sin(t0) ; S1 = sin(t1) ; C0 = cos(t0) ; C1 = cos(t1) ;
#endif /*_USE_SINCOS_*/

  D0 = sqrt(1.0 - al2*S0*S0) ; D1 = sqrt(1.0 - al2*S1*S1) ;

  if ( al != 0.0 ) {
    Q[off-3] = (C1/D1 - C0/D0)/al2 - 
      (log(al*C1 + D1) - log(al*C0 + D0))/al2/al ;
    Q[off-2] = -(log(D1) - log(D0))/al2 ;
    Q[off-1] = -(log(al*C1 + D1) - log(al*C0 + D0))/al ;
  } else {
    Q[off-3] = -(C1*C1*C1 - C0*C0*C0)/3.0 ;
    Q[off-2] = (S1*S1 - S0*S0)/2.0 ;
    Q[off-1] = -(C1 - C0) ;
  }
  Q[off+0] = -(log(C1) - log(C0)) ;

  d1n = D1/C1 ; d0n = D0/C0 ;
  for ( n = 1 ; n <= N ; n ++ ) {
    Q[off+n] = al2*Q[off+n-2] + (d1n - d0n)/n ;
    d1n *= D1/C1 ; d0n *= D0/C0 ;
  }

  return 0 ;
}

gint quad_delta_C_n_alpha(gdouble al, gdouble *Q, gint q,
			  gdouble *I0, gdouble *I1, gdouble *I2, gdouble *I3)

/*
  \int (\Delta/\cos\theta-\alpha)^q (\Delta/\cos\theta)^{-s}, s=0,1,2,3
  
  Q: generated using quad_delta_C_n

*/

{
  gint qstr, bstr ;
  gdouble *bn, *cft ;

  bn = _binomial_list(q) ; bstr = 1 ;
  qstr = -1 ;

  cft = &(Q[q+3]) ;
  *I0 = polyval2(bn, bstr, cft, qstr, q, -al) ;

  cft = &(Q[q+2]) ;
  *I1 = polyval2(bn, bstr, cft, qstr, q, -al) ;

  cft = &(Q[q+1]) ;
  *I2 = polyval2(bn, bstr, cft, qstr, q, -al) ;

  cft = &(Q[q+0]) ;
  *I3 = polyval2(bn, bstr, cft, qstr, q, -al) ;

  return 0 ;
}

gint htri_tri_z_parameters(gdouble r1, gdouble r2, gdouble th, 
			   gdouble phi, gdouble s, gdouble z,
			   gdouble *S, gdouble *Rmax, gdouble *al)

{
  *S = s*s + z*z ;

  *al = sqrt(z*z/(*S)) ;

  *Rmax = MAX(r1,r2) ;
  *Rmax = sqrt((*Rmax)*(*Rmax) + z*z) ;

  *S = sqrt(*S) ;

  return 0 ;
}

gint htri_tri_parameters(gdouble r1, gdouble r2, gdouble th, 
			 gdouble *phi, gdouble *s)

{
  *phi = atan2(r1 - r2*cos(th), r2*sin(th)) ;
  *s = r1*cos(*phi) ;

  return 0 ;
}

gint htri_quad_triangle_0(gdouble r1, gdouble r2, gdouble th, gdouble z,
			  gdouble k, gdouble tol,
			  gdouble *d0G, gdouble *d1G, gdouble *d2G)

/*
 * integrals for wavenumber k with zero order source terms over
 * reference triangle with sides r1 and r2 at the origin spanning
 * angle theta and observer at (0,0,z), evaluated to tolerance tol
 *
 * Output: d0G, d1G, d2G, \partial^{i}I/\partial z^{i}, i=0,1,2,
 *
 * I = \int_{A} \exp[j k R]/4\pi R dA
 */

{
  gdouble i0[2] = {0}, i1[2] = {0}, i2[2] = {0}, i3[2] = {0}, Iq[128] ;
  gdouble al, phi, s, S, Rmax, skz, ckz, za, sgnz, kSn ;
  gdouble t0, t1, t2, t3 ;
  gdouble *cfft_cos, *cfft_sin ;
  gint ncs, ncc, q, Q ;
  /* g_assert(Q < 31) ; */

  htri_tri_parameters(r1, r2, th, &phi, &s) ;
  t0 = -phi ; t1 = th - phi ;

  za = fabs(z) ; sgnz = ( z < 0.0 ? -1.0 : 1.0) ;

  ckz = cos(k*za) ; skz = sin(k*za) ;
    
  htri_tri_z_parameters(r1, r2, th, phi, s, za, &S, &Rmax, &al) ;
  htri_select_expansion(k*(Rmax-za), tol, &cfft_cos, &ncc, &cfft_sin, &ncs) ;

  Q = MAX(ncs, ncc) ;

  /*base integrals*/
  quad_delta_C_n(al, t0, t1, Q, Iq) ;
  kSn = 1.0 ;
  for ( q = 0 ; q <= Q ; q ++ ) {
    quad_delta_C_n_alpha(al, Iq, q, &t0, &t1, &t2, &t3) ;
    i0[0] += cfft_cos[q]*kSn*t0 ; i0[1] += cfft_sin[q]*kSn*t0 ;
    i1[0] += cfft_cos[q]*kSn*t1 ; i1[1] += cfft_sin[q]*kSn*t1 ;
    i2[0] += cfft_cos[q]*kSn*t2 ; i2[1] += cfft_sin[q]*kSn*t2 ;
    i3[0] += cfft_cos[q]*kSn*t3 ; i3[1] += cfft_sin[q]*kSn*t3 ;
    kSn *= k*S ;
  }
  
  t0 = i0[1]/k ; t1 = (-i0[0] + th)/k ;
  d0G[0] = ckz*t0 - skz*t1 ; d0G[1] = ckz*t1 + skz*t0 ;

  t0 = sgnz*(al*i1[0] - th) ; t1 = sgnz*al*i1[1] ;
  d1G[0] = ckz*t0 - skz*t1 ; d1G[1] = ckz*t1 + skz*t0 ;

  t0 = (i1[0] - k*S*al*al*i2[1] - al*al*i3[0])/S ;
  t1 = (i1[1] + k*S*al*al*i2[0] - al*al*i3[1] - k*S*th)/S ;
  d2G[0] = ckz*t0 - skz*t1 ; d2G[1] = ckz*t1 + skz*t0 ;

  return 0 ;
}

gint htri_quad_triangle_1_0(gdouble r1, gdouble r2, gdouble th,
			    gdouble k, gdouble tol,
			    gdouble *d0G, gdouble *d1G, gdouble *d2G)
/*
 * integrals for wavenumber k with first order source terms over
 * reference triangle with sides r1 and r2 at the origin spanning
 * angle theta and observer at (0,0,z=0), evaluated to tolerance tol
 *
 * Output: d0G, d1G, d2G, \partial^{i}I/\partial z^{i}, i=0,1,2,
 *
 * I = \int_{A} (1,x,y) \exp[j k R]/R dA
 */

{
  gdouble i0, i1, i2, j0, j1, j2 ;
  gdouble d0g[6], d1g[6], d2g[6] ;
  gdouble al, phi, s, S, Rmax, skz, ckz, kSq ;
  gdouble tr, ti ;
  gdouble *cfft_cos, *cfft_sin ;
  gdouble Lcos, Lsin ;
  gdouble C0, C1, S0, S1, Cq1, Cq0 ;
  gint ncs, ncc, q, Q ;

  htri_tri_parameters(r1, r2, th, &phi, &s) ;

  ckz = 1.0 ; skz = 0.0 ;

  htri_tri_z_parameters(r1, r2, th, phi, s, 0.0, &S, &Rmax, &al) ;
  htri_select_expansion(k*Rmax, tol, &cfft_cos, &ncc, &cfft_sin, &ncs) ;
  
  Q = MAX(ncs, ncc) ;

  memset(d0g, 0, 6*sizeof(gdouble)) ;
  memset(d1g, 0, 6*sizeof(gdouble)) ;
  memset(d1G, 0, 6*sizeof(gdouble)) ;
  memset(d2g, 0, 6*sizeof(gdouble)) ;
  memset(d2G, 0, 6*sizeof(gdouble)) ;

  S0 = sin(-phi) ; S1 = sin(th-phi) ;
  C0 = cos(-phi) ; C1 = cos(th-phi) ;

  i0 = 0.5*(log((1.0 + S1)/(1.0 - S1)) - log((1.0 + S0)/(1.0 - S0))) ;
  i1 = th ;
  i2 = S1 - S0 ;

  j0 = 1.0/C1 - 1.0/C0 ;
  j1 = -(log(C1) - log(C0)) ;

  Lcos = -(log(s/C1)*S1 - log(s/C0)*S0) + -(S1 - S0) + i0 ;
  Lsin = (log(s/C1)*C1 - log(s/C0)*C0) + (C1 - C0) ;

  Cq1 = C1 ; Cq0 = C0 ; kSq = 1.0 ;
  for ( q = 0 ; q <= Q ; q ++ ) {
    d0g[0] += cfft_cos[q]*kSq*s/(q+1)*i0 ;
    d0g[1] += cfft_sin[q]*kSq*s/(q+1)*i0 ;
    d0g[2] += cfft_cos[q]*kSq*s*s/(q+2)*i0 ;
    d0g[3] += cfft_sin[q]*kSq*s*s/(q+2)*i0 ;
    d0g[4] += cfft_cos[q]*kSq*s*s/(q+2)*j0 ;
    d0g[5] += cfft_sin[q]*kSq*s*s/(q+2)*j0 ;

    if ( q != 1 ) {
      d2g[0] += cfft_cos[q]*kSq/s*i2 ;
      d2g[1] += cfft_sin[q]*kSq/s*i2 ;
    }
    if ( q == 0 ) {
      d2g[2] += cfft_cos[q]*Lcos ;
      d2g[3] += cfft_sin[q]*Lcos ;
      d2g[4] += cfft_cos[q]*Lsin ;
      d2g[5] += cfft_sin[q]*Lsin ;
    }

    i2 = i1 ; i1 = i0 ;
    i0 = 1.0/(q+1)*(S1/Cq1 - S0/Cq0) + (gdouble)(q)/(q+1)*i2 ;
    j2 = j1 ; j1 = j0 ;
    j0 = 1.0/(q+2)*(S1*S1/Cq1/C1 - S0*S0/Cq0/C0) + (gdouble)(q)/(q+2)*j2 ;
    Cq1 *= C1 ; Cq0 *= C0 ;
    kSq *= k*s ;
  }

  tr = d0g[0] ; ti = d0g[1] ;
  d0G[0] = tr*ckz - ti*skz ; d0G[1] = tr*skz + ti*ckz ;

  tr = C0*d0g[2] + S0*d0g[4] ; ti = C0*d0g[3] + S0*d0g[5] ;
  d0G[2] = tr*ckz - ti*skz ; d0G[3] = tr*skz + ti*ckz ;

  tr = -S0*d0g[2] + C0*d0g[4] ; ti = -S0*d0g[3] + C0*d0g[5] ;
  d0G[4] = tr*ckz - ti*skz ; d0G[5] = tr*skz + ti*ckz ;

  tr = d2g[0] ; ti = d2g[1] ;
  d2G[0] = tr*ckz - ti*skz ; d2G[1] = tr*skz + ti*ckz ;

  tr = C0*d2g[2] + S0*d2g[4] ; ti = C0*d2g[3] + S0*d2g[5] ;
  d2G[2] = tr*ckz - ti*skz ; d2G[3] = tr*skz + ti*ckz ;

  tr = -S0*d2g[2] + C0*d2g[4] ; ti = -S0*d2g[3] + C0*d2g[5] ;
  d2G[4] = tr*ckz - ti*skz ; d2G[5] = tr*skz + ti*ckz ;

  return 0 ;
}

gint htri_quad_triangle_1(gdouble r1, gdouble r2, gdouble th, gdouble z,
			  gdouble k, gdouble tol,
			  gdouble *d0G, gdouble *d1G, gdouble *d2G)

/*
 * integrals for wavenumber k with first order source terms over
 * reference triangle with sides r1 and r2 at the origin spanning
 * angle theta and observer at (0,0,z), evaluated to tolerance tol
 *
 * Output: d0G, d1G, d2G, \partial^{i}I/\partial z^{i}, i=0,1,2,
 *
 * I = \int_{A} (1,x,y) \exp[j k R]/R dA
 */

{
  gdouble i0, i1, i2, i3, j0, j1, j2, j3, i1m1, j1m1 ;
  gdouble d0g[10] = {0}, d1g[10] = {0}, d2g[10] = {0} ;
  gdouble Iq[64], Iqt[64] ;
  gdouble al, al2, phi, s, S, Rmax, skz, ckz, za, sgnz, kSq ;
  gdouble t0, t1, tr, ti ;
  gdouble *cfft_cos, *cfft_sin ;
  gdouble Jcos, Jsin, dJcos, dJsin, d2Jcos, d2Jsin, Jcosm1, Jsinm1, ad ;
  gdouble Lcos, Lsin ;
  gdouble A0, A1, C0, C1, S0, S1, D0, D1, L0, L1 ;
  gint ncs, ncc, q, Q ;

  if ( (th < 1e-3) || (fabs(M_PI-th) < 1e-3) ) {  
    memset(d0G, 0, 10*sizeof(gdouble)) ;
    memset(d1G, 0, 10*sizeof(gdouble)) ;
    memset(d2G, 0, 10*sizeof(gdouble)) ;
    
    return 0 ;
  }
  
  if ( z == 0 ) {
    htri_quad_triangle_1_0(r1, r2, th, k, tol, d0G, d1G, d2G) ;

    return 0 ;
  }
  
  htri_tri_parameters(r1, r2, th, &phi, &s) ;
  t0 = -phi ; t1 = th - phi ;

  za = fabs(z) ; sgnz = ( z < 0.0 ? -1.0 : 1.0) ;

  ckz = cos(k*za) ; skz = sin(k*za) ;
    
  htri_tri_z_parameters(r1, r2, th, phi, s, za, &S, &Rmax, &al) ;
  htri_select_expansion(k*(Rmax-za), tol, &cfft_cos, &ncc, &cfft_sin, &ncs) ;

  /* fprintf(stderr, "(%lg, %d, %d) ", k*(Rmax-za), ncc, ncs) ; */
  
  al2 = al*al ; ad = sqrt(1.0 - al2) ;
  Q = MAX(ncs, ncc) ;

  /*base integrals*/
  quad_delta_C_n(al, t0, t1, Q+1, Iq) ;
  quad_delta_C_n_T(al, t0, t1, Q+1, Iqt) ;

#ifdef _USE_SINCOS_
  sincos(-phi, &S0, &C0) ;
  sincos(th-phi, &S1, &C1) ;
#else
  S0 = sin(-phi) ; S1 = sin(th-phi) ;
  C0 = cos(-phi) ; C1 = cos(th-phi) ;
#endif /*_USE_SINCOS_*/

  D0 = sqrt(1.0-al*al*S0*S0) ; D1 = sqrt(1.0-al*al*S1*S1) ;
  L1 = log((D1-ad)/(D1+ad)) ; L0 = log((D0-ad)/(D0+ad)) ;
  A1 = asin(al*S1) ; A0 = asin(al*S0) ;

  /*\int (\cos\theta,\sin\theta)\log (\Delta-\alpha')/(\Delta+\alpha') 
    d\theta*/

  Lcos = S1*L1 - S0*L0 + 
    log((D1+ad*S1)/(D1-ad*S1)) - log((D0+ad*S0)/(D0-ad*S0)) 
    - 2.0*ad/al*(A1 - A0) ;
  Lsin = -C1*L1 + C0*L0 + 2.0*ad/al*(log(al*C1+D1) - log(al*C0+D0)) ;

  q = -1 ; kSq = 1.0 ;
  Jcosm1 = -0.25/k*Lcos ; Jsinm1 = -0.25/k*Lsin ;

  quad_delta_C_n_alpha(al, Iq,  q+1, &i0, &i1, &i2, &i3) ;
  quad_delta_C_n_alpha(al, Iqt, q+1, &j0, &j1, &j2, &j3) ;
  Jcos = 0.5*s*kSq/(q+2)*i0 - k*za*(2*q+3)/(q+2)*Jcosm1 ;
  Jsin = 0.5*s*kSq/(q+2)*j0 - k*za*(2*q+3)/(q+2)*Jsinm1 ;

  dJcos = sgnz*0.25*Lcos + sgnz*0.5*s/S*(A1 - A0)/al ;
  dJsin = sgnz*0.25*Lsin + 
    sgnz*0.25*s/S*(log((D1 - al*C1)/(D1 + al*C1)) -
		   log((D0 - al*C0)/(D0 + al*C0)))/al ;

  d2Jcos = 0.5*(s/S)*(s/S)*(s/S)*(S1/D1 - S0/D0)/za ;
  d2Jsin = 0.5*(s/S)*(s/S)*(s/S)*(-C1/D1 + C0/D0)/ad/ad/za ;

  for ( q = 0 ; q <= Q ; q ++ ) {
    i1m1 = i1 ;
    j1m1 = j1 ;
    quad_delta_C_n_alpha(al, Iq,  q+1, &i0, &i1, &i2, &i3) ;
    quad_delta_C_n_alpha(al, Iqt, q+1, &j0, &j1, &j2, &j3) ;

    d0g[0] += cfft_cos[q]*S*kSq/(q+1)*i0 ;
    d0g[1] += cfft_sin[q]*S*kSq/(q+1)*i0 ;

    tr = (s*S*kSq*i0 + 2.0*za*Jcos)/(q+2) ;
    d0g[2] += cfft_cos[q]*tr ;
    d0g[3] += cfft_sin[q]*tr ;

    tr = (s*S*kSq*j0 + 2.0*za*Jsin)/(q+2) ;
    d0g[4] += cfft_cos[q]*tr ;
    d0g[5] += cfft_sin[q]*tr ;

    tr = -kSq*s/S*i1m1 + 2.0*k*Jcosm1 ;
    d0g[6] += cfft_cos[q]*tr ;
    d0g[7] += cfft_sin[q]*tr ;

    tr = -kSq*s/S*j1m1 + 2.0*k*Jsinm1 ;
    d0g[8] += cfft_cos[q]*tr ;
    d0g[9] += cfft_sin[q]*tr ;


    d1g[0] += -sgnz*cfft_cos[q]*kSq*i1 ;
    d1g[1] += -sgnz*cfft_sin[q]*kSq*i1 ;

    tr = (-(gdouble)(q+1)*kSq*s*i1 + 2.0*Jcos + 2.0*z*dJcos)/(q+2) ;
    d1g[2] += sgnz*cfft_cos[q]*tr ;
    d1g[3] += sgnz*cfft_sin[q]*tr ;

    tr = (-(gdouble)(q+1)*kSq*s*j1 + 2.0*Jsin + 2.0*z*dJsin)/(q+2) ;
    d1g[4] += sgnz*cfft_cos[q]*tr ;
    d1g[5] += sgnz*cfft_sin[q]*tr ;

    d2g[0] += cfft_cos[q]*kSq*(al*i3 + (q+1)*i2)/S ;
    d2g[1] += cfft_sin[q]*kSq*(al*i3 + (q+1)*i2)/S ;

    tr = (gdouble)(q+1)/(q+2)*kSq*s/S*(al*i3 + (q+1)*i2)
      + sgnz*4.0/(q+2)*dJcos + 2.0*za/(q+2)*d2Jcos ;
    d2g[2] += cfft_cos[q]*tr ; d2g[3] += cfft_sin[q]*tr ;

    tr = (gdouble)(q+1)/(q+2)*kSq*s/S*(al*j3 + (q+1)*j2)
      + sgnz*4.0/(q+2)*dJsin + 2.0*za/(q+2)*d2Jsin ;
    d2g[4] += cfft_cos[q]*tr ; d2g[5] += cfft_sin[q]*tr ;

    kSq *= k*S ;

    d2Jcos = 0.5*(q+1)/(q+2)*kSq*s/S/S*(al*i3  + (q+1)*i2)
      - sgnz*2.0*k*(2*q+3)/(q+2)*dJcos - za*k*(2*q+3)/(q+2)*d2Jcos ;
    d2Jsin = 0.5*(q+1)/(q+2)*kSq*s/S/S*(al*j3 + (q+1)*j2)
      - sgnz*2.0*k*(2*q+3)/(q+2)*dJsin - za*k*(2*q+3)/(q+2)*d2Jsin ;

    dJcos = -sgnz*0.5*s/S*kSq*(q+1)/(q+2)*i1
      - sgnz*k*(2*q+3)/(q+2)*Jcos - k*za*(2*q+3)/(q+2)*dJcos ;
    dJsin = -sgnz*0.5*s/S*kSq*(q+1)/(q+2)*j1
      - sgnz*k*(2*q+3)/(q+2)*Jsin - k*za*(2*q+3)/(q+2)*dJsin ;

    Jcosm1 = Jcos ; Jsinm1 = Jsin ;
    Jcos = 0.5*s*kSq/(q+2)*i0 - k*za*(2*q+3)/(q+2)*Jcos ;
    Jsin = 0.5*s*kSq/(q+2)*j0 - k*za*(2*q+3)/(q+2)*Jsin ;
  }

  tr = d0g[0] ; ti = d0g[1] ;
  d0G[0] = tr*ckz - ti*skz ; d0G[1] = tr*skz + ti*ckz ;

  tr = C0*d0g[2] + S0*d0g[4] ; ti = C0*d0g[3] + S0*d0g[5] ;
  d0G[2] = tr*ckz - ti*skz ; d0G[3] = tr*skz + ti*ckz ;

  tr = -S0*d0g[2] + C0*d0g[4] ; ti = -S0*d0g[3] + C0*d0g[5] ;
  d0G[4] = tr*ckz - ti*skz ; d0G[5] = tr*skz + ti*ckz ;

  tr = C0*d0g[6] + S0*d0g[8] ; ti = C0*d0g[7] + S0*d0g[9] ;
  d0G[6] = tr*ckz - ti*skz ; d0G[7] = tr*skz + ti*ckz ;

  tr = -S0*d0g[6] + C0*d0g[8] ; ti = -S0*d0g[7] + C0*d0g[9] ;
  d0G[8] = tr*ckz - ti*skz ; d0G[9] = tr*skz + ti*ckz ;


  tr = d1g[0] - sgnz*k*d0g[1] ; ti = d1g[1] + sgnz*k*d0g[0] ;
  d1G[0] = tr*ckz - ti*skz ; d1G[1] = tr*skz + ti*ckz ;

  tr =  C0*d1g[2] + S0*d1g[4] - sgnz*k*(C0*d0g[3] + S0*d0g[5]) ;
  ti =  C0*d1g[3] + S0*d1g[5] + sgnz*k*(C0*d0g[2] + S0*d0g[4]) ;
  d1G[2] = tr*ckz - ti*skz ; d1G[3] = tr*skz + ti*ckz ;

  tr =  -S0*d1g[2] + C0*d1g[4] - sgnz*k*(-S0*d0g[3] + C0*d0g[5]) ;
  ti =  -S0*d1g[3] + C0*d1g[5] + sgnz*k*(-S0*d0g[2] + C0*d0g[4]) ;
  d1G[4] = tr*ckz - ti*skz ; d1G[5] = tr*skz + ti*ckz ;

  tr = d2g[0] - sgnz*2*k*d1g[1] - k*k*d0g[0] ;
  ti = d2g[1] + sgnz*2*k*d1g[0] - k*k*d0g[1] ;
  d2G[0] = tr*ckz - ti*skz ; d2G[1] = tr*skz + ti*ckz ;

  tr = -k*k*(C0*d0g[2] + S0*d0g[4])
    - sgnz*2*k*(C0*d1g[3] + S0*d1g[5]) + 
    C0*d2g[2] + S0*d2g[4] ;
  ti = -k*k*(C0*d0g[3] + S0*d0g[5])
    + sgnz*2*k*(C0*d1g[2] + S0*d1g[4]) + 
    C0*d2g[3] + S0*d2g[5] ;
  d2G[2] = tr*ckz - ti*skz ; d2G[3] = tr*skz + ti*ckz ;

  tr = -k*k*(-S0*d0g[2] + C0*d0g[4])
    - sgnz*2*k*(-S0*d1g[3] + C0*d1g[5]) + 
    -S0*d2g[2] + C0*d2g[4] ;
  ti = -k*k*(-S0*d0g[3] + C0*d0g[5])
    + sgnz*2*k*(-S0*d1g[2] + C0*d1g[4]) + 
    -S0*d2g[3] + C0*d2g[5] ;
  d2G[4] = tr*ckz - ti*skz ; d2G[5] = tr*skz + ti*ckz ;

  return 0 ;
}

/*
  Integrate Green's functions weighted on linear shape functions over
  general triangle with field point in general position
*/

gint htri_quad_shape_1(gdouble *xf, 
		       gdouble *x1, gdouble *x2, gdouble *x3, 
		       gdouble k, gdouble tol, gint qmax,
		       gdouble *g, gdouble *dg, gdouble *d2g)

{
  gdouble s[3], t[3], n[3], og[3], *y1, *y2, *y3, yf[3], smax, smin ;
  gdouble psi[3]={0}, r1[3]={0}, r2[3]={0}, th[3]={0} ;
  gdouble ia[10], iz[10], izz[10] ;
  gdouble A[10], Ai[9], I1, Ix, Iy, C, S ;
  gint ntri, j, sgn[3], q, in ;

  htri_triangle_axes(x1, x2, x3, xf, og, s, t, n) ;
  y1 = &(A[1]) ; y2 = &(A[4]) ; y3 = &(A[7]) ;

  htri_triangle_project(og, s, t, n, x1, y1) ;
  htri_triangle_project(og, s, t, n, x2, y2) ;
  htri_triangle_project(og, s, t, n, x3, y3) ;
  htri_triangle_project(og, s, t, n, xf, yf) ;

  if ( fabs(yf[2]) < 1e-9 ) yf[2] = 0.0 ;

  in = htri_triangle_decompose(y1, y2, y3, yf, psi, th, r1, r2, sgn, &ntri) ;

  smax = MAX(r1[0], r2[0]) ; smin = MIN(r1[0], r2[0]) ; 
  for ( j = 1 ; j < ntri ; j ++ ) {
    g_assert(!isnan(th[j])) ;
    smax = MAX(smax, MAX(r1[j], r2[j])) ;
    smin = MIN(smin, MIN(r1[j], r2[j])) ;
  }

  if ( in >= 0 ) smin = 0.0 ;
  q = htri_quad_order_1R(smin, smax, yf[2], qmax, tol) ;
  if ( q != -1 ) return q ;

  /* g_assert(k*smax < 1.57) ; */
  if ( k*smax > 0.5 ) return -1 ;
  if ( k*k*(smax*smax + yf[2]*yf[2]) > 8 ) return -1 ;
  
  A[0] = A[3] = A[6] = 1.0 ;
  _htri_invert3x3(Ai, A) ;

  memset(g, 0, 6*sizeof(gdouble)) ;
  memset(dg, 0, 6*sizeof(gdouble)) ;
  memset(d2g, 0, 6*sizeof(gdouble)) ;
  
  for ( j = 0 ; j < ntri ; j ++ ) {
    htri_quad_triangle_1(r1[j], r2[j], th[j], yf[2], k, tol, ia, iz, izz) ;
    /* fprintf(stderr, "%d %lg %lg %lg %lg %lg %lg %lg %lg\n", */
    /* 	    j, ia[0], ia[1], iz[0], iz[1], izz[0], izz[1], yf[2], M_PI-th[j] ) ; */
    /* fprintf(stderr, "\n") ; */
    C = cos(psi[j])*sgn[j] ; S = sin(psi[j])*sgn[j] ;
    
    I1 = ia[0]*sgn[j] ; Ix = C*ia[2] - S*ia[4] ; Iy = C*ia[4] + S*ia[2] ;
    g[0] += Ai[0]*I1 + Ai[3]*Ix + Ai[6]*Iy ;
    g[2] += Ai[1]*I1 + Ai[4]*Ix + Ai[7]*Iy ;
    g[4] += Ai[2]*I1 + Ai[5]*Ix + Ai[8]*Iy ;
    
    I1 = ia[1]*sgn[j] ; Ix = C*ia[3] - S*ia[5] ; Iy = C*ia[5] + S*ia[3] ;
    g[1] += Ai[0]*I1 + Ai[3]*Ix + Ai[6]*Iy ;
    g[3] += Ai[1]*I1 + Ai[4]*Ix + Ai[7]*Iy ;
    g[5] += Ai[2]*I1 + Ai[5]*Ix + Ai[8]*Iy ;
    
    I1 = iz[0]*sgn[j] ; Ix = C*iz[2] - S*iz[4] ; Iy = C*iz[4] + S*iz[2] ;
    dg[0] += Ai[0]*I1 + Ai[3]*Ix + Ai[6]*Iy ;
    dg[2] += Ai[1]*I1 + Ai[4]*Ix + Ai[7]*Iy ;
    dg[4] += Ai[2]*I1 + Ai[5]*Ix + Ai[8]*Iy ;
    
    I1 = iz[1]*sgn[j] ; Ix = C*iz[3] - S*iz[5] ; Iy = C*iz[5] + S*iz[3] ;
    dg[1] += Ai[0]*I1 + Ai[3]*Ix + Ai[6]*Iy ;
    dg[3] += Ai[1]*I1 + Ai[4]*Ix + Ai[7]*Iy ;
    dg[5] += Ai[2]*I1 + Ai[5]*Ix + Ai[8]*Iy ;
    
    I1 = izz[0]*sgn[j] ; Ix = C*izz[2] - S*izz[4] ; Iy = C*izz[4] + S*izz[2] ;
    d2g[0] += Ai[0]*I1 + Ai[3]*Ix + Ai[6]*Iy ;
    d2g[2] += Ai[1]*I1 + Ai[4]*Ix + Ai[7]*Iy ;
    d2g[4] += Ai[2]*I1 + Ai[5]*Ix + Ai[8]*Iy ;
    
    I1 = izz[1]*sgn[j] ; Ix = C*izz[3] - S*izz[5] ; Iy = C*izz[5] + S*izz[3] ;
    d2g[1] += Ai[0]*I1 + Ai[3]*Ix + Ai[6]*Iy ;
    d2g[3] += Ai[1]*I1 + Ai[4]*Ix + Ai[7]*Iy ;
    d2g[5] += Ai[2]*I1 + Ai[5]*Ix + Ai[8]*Iy ;

    if ( isnan(g[0]) || isnan(dg[0]) || isnan(d2g[0]) )
      g_error("%s: NaN for k = %lg\n"
	      "  g=%lg+j%lg\n"
	      "  dg=%lg+j%lg\n"
	      "  d2g=%lg+j%lg\n"
	      "  xf=(%lg,%lg,%lg), (%lg,%lg,%lg)\n"
	      "  x1=(%lg,%lg,%lg), (%lg,%lg,%lg)\n"
	      "  x2=(%lg,%lg,%lg), (%lg,%lg,%lg)\n"
	      "  x3=(%lg,%lg,%lg), (%lg,%lg,%lg)\n"
	      "  s=(%lg,%lg,%lg), t=(%lg,%lg,%lg),\n"
	      "  n=(%lg,%lg,%lg)\n"
	      "  r1=(%lg,%lg,%lg), r2=(%lg,%lg,%lg),\n"
	      "  th=(%lg,%lg,%lg), psi=(%lg,%lg,%lg),\n"
	      "  th-pi=(%lg,%lg,%lg),\n"
	      "  ntri=%d, j=%d",
	      __FUNCTION__, k,
	      g[0], g[1], dg[0], dg[1], d2g[0], d2g[1],
	      xf[0], xf[1], xf[2], yf[0], yf[1], yf[2],
	      x1[0], x1[1], x1[2], y1[0], y1[1], y1[2],
	      x2[0], x2[1], x2[2], y2[0], y2[1], y2[2],
	      x3[0], x3[1], x3[2], y3[0], y3[1], y3[2],
	      s[0], s[1], s[2],
	      t[0], t[1], t[2],
	      n[0], n[1], n[2],
	      r1[0], r1[1], r1[2],
	      r2[0], r2[1], r2[2],
	      th[0], th[1], th[2],
	      psi[0], psi[1], psi[2],
	      th[0]-M_PI, th[1]-M_PI, th[2]-M_PI,
	      ntri, j
	      ) ;

  }
  
  return 0 ;
}

/*
  Integrate Green's functions weighted on constant shape functions
  over general triangle with field point in general position
*/

gint htri_quad_shape_0(gdouble *xf, 
		       gdouble *x1, gdouble *x2, gdouble *x3, 
		       gdouble k, gdouble tol, gint qmax,
		       gdouble *g, gdouble *dg, gdouble *d2g)

{
  gdouble s[3], t[3], n[3], og[3], y1[3], y2[3], y3[3], yf[3], smax, smin ;
  gdouble psi[3]={0}, r1[3]={0}, r2[3]={0}, th[3]={0} ;
  gdouble ia[2], iz[2], izz[2] ;
  gint ntri, j, sgn[3], q, in ;

  htri_triangle_axes(x1, x2, x3, xf, og, s, t, n) ;

  htri_triangle_project(og, s, t, n, x1, y1) ;
  htri_triangle_project(og, s, t, n, x2, y2) ;
  htri_triangle_project(og, s, t, n, x3, y3) ;
  htri_triangle_project(og, s, t, n, xf, yf) ;

  /* if ( _bem3d_trace_set(0) && _bem3d_trace_set(1) ) */
  /*   fprintf(stderr, "Hello\n") ; */
  
  if ( fabs(yf[2]) < 1e-9 ) yf[2] = 0.0 ;
  
  in = htri_triangle_decompose(y1, y2, y3, yf, psi, th, r1, r2, sgn, &ntri) ;

  g_assert(ntri >= 1) ;
  
  smax = MAX(r1[0], r2[0]) ; smin = MIN(r1[0], r2[0]) ; 
  for ( j = 1 ; j < ntri ; j ++ ) {
    smax = MAX(smax, MAX(r1[j], r2[j])) ;
    smin = MIN(smin, MIN(r1[j], r2[j])) ;
  }

  if ( in >= 0 ) smin = 0.0 ;
  q = htri_quad_order_1R(smin, smax, yf[2], qmax, tol) ;
  if ( q != -1 ) return q ;

  /* g_assert(k*smax < 1.57) ; */
  if ( k*smax > 1.57 ) return -1 ;
  if ( k*k*(smax*smax + yf[2]*yf[2]) > 15 ) return -1 ;
  
  g[0] = g[1] = dg[0] = dg[1] = d2g[0] = d2g[1] = 0.0 ;
  
  for ( j = 0 ; j < ntri ; j ++ ) {
    htri_quad_triangle_0(r1[j], r2[j], th[j], yf[2], k, tol, ia, iz, izz) ;

    /* fprintf(stderr, "%d %lg %lg %lg %lg %lg %lg\n", */
    /* 	    j, ia[0], ia[1], iz[0], iz[1], izz[0], izz[1]) ; */
    
    g[0] += ia[0]*sgn[j] ;
    g[1] += ia[1]*sgn[j] ;
    
    dg[0] += iz[0]*sgn[j] ;
    dg[1] += iz[1]*sgn[j] ;
    
    d2g[0] += izz[0]*sgn[j] ;
    d2g[1] += izz[1]*sgn[j] ;
  }

  /* if ( yf[2] == 0.0 ) dg[0] = dg[1] = 0.0 ; */
  
  return 0 ;
}

gint htri_quad_order_1R(gdouble rmin, gdouble rmax, gdouble z,
			gint qmax, gdouble tol)

/*
  Estimate the order of expansion needed to approximate 1/R to
  tolerance tol with R=(r^{2}+z^{2})^{1/2}, rmin<=r<=rmax

  on exit, return order of expansion q which gives error less than
  tol, or -1 if this cannot be achieved for order less than qmax
*/

{
  gdouble P0, P1, tt, rb, dr, rho, rq ;
  gint q ;

  rb = 0.5*(rmax + rmin) ;
  rho = sqrt(z*z + rb*rb) ;
  dr = 0.5*(rmax - rmin)/rho ;
  rb /= rho ;

  rq = dr/rho ; P1 = rb ; P0 = 1.0 ;
  for ( q = 0 ; q < qmax ; q ++ ) {
    if ( (fabs(rq*P1) < tol) && (q > 3) ) return q ;
    tt = P1 ; 
    P1 = ((2.0*q+3)*rb*P1 - (gdouble)(q+1)*P0)/(q+2) ;
    P0 = tt ;
    rq *= dr ;
  }
  
  return -1 ;
}
