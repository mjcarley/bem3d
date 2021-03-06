/* series.c
 * 
 * Copyright (C) 2012, 2013, 2018 by Michael Carley
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

#define SERIES_MAX_TERMS 4096

gint quad_sin_cos_nmr_hw(gint m, gint n, gint r, gdouble k,
			 gdouble t, gdouble S, gdouble C, gdouble *I) ;
gdouble quad_sin_cos(gint m, gint n, gdouble t, gdouble S, gdouble C) ;
gdouble quad_sin_cos_nmr(gint m, gint n, gint r, gdouble k,
			 gdouble t, gdouble S, gdouble C) ;

static gdouble quad_cos_sin_n(gdouble t, gdouble S, gdouble C, gint n) 

/*\int \cos^{n}t/\sin t dt*/

{
  gdouble I = 0.0, dI, sgn ;
  gint k, n0 ;

  if ( t == 0.0 ) return 0.0 ;
  sgn = ( t < 0.0 ? -1.0 : 1.0 ) ;
  if ( 2*((n0 = n/2)) == n ) {
    I = log(tan(0.5*fabs(t))) ;
    dI = C ;
    for ( k = 1 ; k <= n0 ; ( dI *= C*C), (k ++) ) I += dI/(2*k-1) ;
    return sgn*I ;
  }

  I = log(fabs(S)) ; dI = 1.0 ;
  for ( k = 1 ; k <= n0 ; k ++ ) I += ((dI *= C*C)/(2*k)) ;

  return sgn*I ;
}

static gdouble quad_sin_cos_n(gdouble t, gdouble S, gdouble C, gint n) 

/*\int \sin^{n}t/\cos t dt*/

{
  gdouble I = 0.0, dI ;
  gint k, n0 ;

  if ( n >= 0 ) {
    if ( 2*((n0 = n/2)) == n ) {
      I = log(tan(M_PI*0.25 + 0.5*t)) ;
      dI = -S ;
      for ( k = 1 ; k <= n0 ; (dI *= S*S), (k ++) ) 
	I += dI/(2*k-1) ;
      return I ;
    }
    
    I = -log(C) ; dI = -1.0 ;
    for ( k = 1 ; k <= n0 ; k ++ ) I += ((dI *= S*S)/(2*k)) ;

    return I ;
  }

  if ( 2*((n0 = -n/2)) == -n ) {
    I = log(tan(M_PI*0.25 + 0.5*t)) ;
    dI = pow(1.0/S,2*n0+1) ;
    for ( k = 1 ; k <= n0 ; k ++ ) I -= (dI *= S*S)/(2*n0-2*k+1) ;
    return I ;
  }

  I = log(S/C) ;
  dI = pow(1.0/S,2*n0+2) ;
  for ( k = 1 ; k <= n0 ; k ++ ) I -= (dI *= S*S)/(2*n0-2*k+2) ;

  return I ;
}

static gdouble quad_cos_n(gdouble t, gdouble S, gdouble C, gint n) 

/*\int \cos^{n}t dt*/

{
  gdouble I = 0.0, tt, dI ;
  gint n0, k ;

  switch ( n ) {
  default: break ;
  case 0: return t ; break ;
  case 1: return S ; break ;
  case 2: return 0.5*(S*C + t) ; break ;
  case 3: return S*(1.0 - S*S/3.0) ; break ;
  case 4: return (S*C*(3 + 2.0*C*C) +3.0*t)/8.0 ; break ;
  case 5: return S*(12 - S*S*4 + 3*C*C*C*C)/15.0 ; break ;
  case 6: return 0.31250*(t + S*C) + S*C*C*C*(5.0 + 4*C*C)/24.0 ; break ;
  case 7: return S*(24.0 - 8*S*S + C*C*C*C*(6 + 5*C*C))/35.0 ; break ;
  }
  
  if ( 2*(n0 = n/2) == n ) {
    /*n even*/
    I = dI = S/n*pow(C,n-1) ; tt = 0.5/n0*t ;
    for ( k = 1 ; k <= n0-1 ; k ++ ) {
      I += (dI *= 0.5*(n-2*k+1)/(n0-k)/C/C) ;
      tt *= 0.5*(2*n0+1-2*k)/(n0-k) ;
    }
    I += tt ;

    return I ;
  }

  I = dI = S/n*pow(C,n-1) ;
  for ( k = 0 ; k <= n0-1 ; k ++ ) I += (dI *= 2.0*(n0-k)/(2*n0-2*k-1)/C/C) ;

  return I ;
}

static gdouble quad_cos_minus_n(gdouble t, gdouble S, gdouble C, gint n) 

/*\int \cos^{-n}t dt*/

{
  gdouble I = 0.0, T = S/C, L, dI ;
  gint n0, k ;

  switch ( n ) {
  default: break ;
  case 0: return t ; break ;
  case 1: return log(sqrt((1+S)/(1-S))) ; break ;
  case 2: return T ; break ;
  case 3: return 0.5*T/C + 0.5*log(sqrt((1+S)/(1-S))) ; break ;
  case 4: return (1.0/C/C + 2.0)*T/3.0 ; break ;
  case 5: return (0.25/C/C + 0.375)*T/C +
      0.375*log(sqrt((1+S)/(1-S))) ;
    break ;
  case 6: return T*((0.2/C/C + 4.0/15.0*S*S)/C/C + 0.8) ; break ;
  case 7: return T/C*((1.0/C/C + 1.25)/6.0/C/C + 0.3125) +
      0.3125*log(sqrt((1+S)/(1-S))) ;
    break ;
  case 8: return T*(1.0 + T*T*(1.0 + T*T*(0.6 + T*T/7.0))) ; break ;
  }

  if ( 2*(n0=n/2) == n ) {
    /*n even n=2*n0*/
    I = dI = pow(1/C,n-1)*S/(n-1) ;
    for ( k = 1 ; k <= n0-1 ; k ++ ) I += (dI *= C*C*(n-2*k)/(n-2*k-1)) ;

    return I ;
  }

  I = dI = pow(1/C,2*n0)*S/(n-1) ; ;
  L = log(sqrt((1+S)/(1-S)))*0.5*(n-2)/n0 ;
  for ( k = 1 ; k <= n0-1 ; k ++ ) {    
    I += (dI *= C*C*(n0-k+0.5)/(n0-k)) ;
    L *= (gdouble)(n0-k-0.5)/k ;
  }
  I += L ;

  return I ;
}

static gdouble quad_sin_n(gdouble t, gdouble S, gdouble C, gint n) 

/*\int \sin^{n}t dt*/

{
  gdouble I = 0.0, dI, tt ;
  gint n0, k ;

  switch ( n ) {
  default: break ;
  case 0: return t ; break ;
  case 1: return -C ; break ;
  case 2: return -0.5*(S*C - t) ; break ;
  case 3: return (C*C/3.0 - 1.0)*C ; break ;
  case 4: return (3*t - S*C*(3 + 2*S*S))/8 ; break ;
  case 5: return C*(4*C*C - (3*S*S*S*S+12))/15.0 ; break ;
  case 6: return 0.31250*t - S*C*(0.3125 + S*S*(5.0/24 + S*S/6.0)) ; break ;
  case 7: return -C*(24 - 8*C*C + S*S*S*S*(6.0 + S*S*5))/35.0 ; break ;
  }

  if ( 2*(n0=n/2) == n ) {
    /*n even n=2*n0*/
    I = dI = -pow(S,n-1)*C/n ; tt = 0.5/n0*t ;
    for ( k = 1 ; k <= n0-1 ; k ++ ) {
      I += (dI *= 0.5*(2*n0+1-2*k)/(n0-k)/S/S ) ;
      tt *= 0.5*(2*n0+1-2*k)/(n0-k) ;
    }
    I += tt ;

    return I ;
  }

  I = dI = -pow(S,2*n0)*C/n ;
  for ( k = 0 ; k <= n0-1 ; k ++ ) I += (dI *= 2.0*(n0-k)/(2*n0-2*k-1)/S/S) ;

  return I ;
}

static gdouble quad_sin_minus_n(gdouble t, gdouble S, gdouble C, gint n) 

/*\int \sin^{n}t dt*/

{
  gdouble I = 0.0, dI, tt, T = S/C ;
  gint n0, k ;

  switch ( n ) {
  default: break ;
  case 0: return t ; break ;
  case 1: return log(tan(0.5*t)) ; break ;
  case 2: return -1.0/T ; break ;
  case 3: return -0.5*C/S/S + 0.5*log(tan(0.5*t)) ; break ;
  case 4: return -C/S*(1.0/S/S + 2.0)/3.0 ; break ;
  case 5: return (-C/S/S*(3.0 + 2.0/S/S) + 3*log(tan(0.5*t)))/8.0 ; break ;
  case 6: return -(1 + (2.0/3.0 + 0.2/T/T)/T/T)/T ; break ;
  case 7: return -C/S/S/6.0*(15.0/8.0 + (1.25 + 1.0/S/S)/S/S) + 
      0.3125*log(tan(0.5*t)) ; break ;
  case 8: return -(1.0 + (1.0 + (0.6 + 1.0/7/T/T)/T/T)/T/T)/T ; break ;
  }

  if ( 2*(n0=n/2) == n ) {
    /*n even n == 2*n0*/
    I = dI = -pow(1.0/S,n-1)*C/(n-1) ;

    for ( k = 1 ; k <= n0-1 ; k ++ ) I += (dI *= 2.0*(n0-k)/(2*n0-2*k-1)*S*S) ;

    return I ;
  }

  I = dI = -pow(1.0/S,2*n0)*C/(2*n0) ;
  tt = 0.5/n0*log(tan(0.5*t)) ;
  for ( k = 1 ; k <= n0-1 ; k ++ ) {
    I += (dI *= 0.5*(2*n0+1-2*k)/(n0-k)*S*S ) ;
    tt *= 0.5*(2*n0+1-2*k)/(n0-k) ;
  }
  I += tt ;

  return I ;
}

gdouble quad_sin_cos(gint m, gint n, gdouble t, gdouble S, gdouble C)

/*
  integral of \sin^{m}t \cos^{n}t
  t:   limit of integration
  S:   \sin t
  C:   \cos t
 */

{
  gdouble I = 0.0, dI, sgn ;
  gint p, m0, n0, k ;

  if ( m >= 0 && n >= 0 ) {
    if ( 2*((m0 = m/2)) == m ) {
      I = quad_cos_n(t, S, C, n) ;
      for ( p = 2 ; p <= m ; p += 2 ) 
	I = -pow(S, p-1)*pow(C,n+1)/(p+n) + I*(p-1)/(p+n) ;
      return I ;
    }
    I = dI = -pow(S,2*m0)*pow(C,n+1)/(2*m0+n+1) ;
    for ( k = 1 ; k <= m0 ; k ++ )
      I += ( dI *= 2.0*(m0-k+1)/(2*m0+n-2*k+1)/S/S ) ;
    return I ;
  }

  if ( m == -1 && n >=0 ) return quad_cos_sin_n(t, S, C, n) ;

  if ( m < 0 && n >= 0 ) {
    if ( 2*((m0 = -m/2)) == -m ) {
      I = dI = -pow(C,n+1)*pow(1.0/S,2*m0-1)/(2*m0-1) ; 
      sgn = -(gdouble)n/(2*m0-1) ;
      for ( k = 1 ; k <= m0-1 ; k ++ ) {
	I += ( dI *= S*S*(2*m0-n-2*k)/(2*m0-2*k-1) ) ;
	sgn *= (gdouble)(2*m0-n-2*k)/(2*m0-2*k-1) ;
      }
      I += sgn*quad_cos_n(t, S, C, n) ;
      return I ;
    }

    I = dI = -pow(C,n+1)*pow(1.0,2*m0)/(2*m0) ; 
    sgn = (gdouble)(1-n)/m0*0.5 ;
    for ( k = 1 ; k <= m0-1 ; k ++ ) {
      I += ( dI *= S*S*0.5*(2*m0-n-2*k+1)/(m0-k) ) ;
      sgn *= (gdouble)(2*m0-n-2*k+1)/(m0-k)*0.5 ;
    }
    I /= pow(S,2*m0) ;
    I += sgn*quad_cos_sin_n(t, S, C, n) ;
    return I ;
  }

  if ( m ==  0 ) return quad_cos_minus_n(t,S,C,-n) ;
  if ( n == -1 ) return quad_sin_cos_n(t,S,C,m) ;

  if ( 2*((n0 = -n/2)) == -n ) {
    I = dI = pow(1/C, 2*n0-1)*pow(S, m+1)/(2*n0-1) ;
    sgn = -(gdouble)m/(2*n0-1) ;
    for ( k = 1 ; k <= n0-1 ; k ++ ) {
      I += ( dI *= C*C*(2*n0-m-2*k)/(2*n0-2*k-1) ) ;
      sgn *= (gdouble)(2*n0-m-2*k)/(2*n0-2*k-1) ;
    }
    I += sgn*( m >= 0 ? quad_sin_n(t, S, C, m) : 
	       quad_sin_minus_n(t, S, C, -m)) ;
    return I ;
  }

  I = dI = pow(S,m+1)*pow(1.0/C,2*n0)/(2*n0) ; 
  sgn = (gdouble)(1-m)/n0*0.5 ;
  for ( k = 1 ; k <= n0-1 ; k ++ ) {
    I += ( dI *= C*C*0.5*(2*n0-m-2*k+1)/(n0-k) ) ;
    sgn *= (gdouble)(2*n0-m-2*k+1)/(n0-k)*0.5 ;
  }
  I += sgn*quad_sin_cos_n(t, S, C, m) ;

  return I ;
}

gdouble quad_sin_cos_nmr(gint m, gint n, gint r, gdouble k,
			 gdouble t, gdouble S, gdouble C)

/*
  integral of \sin^{m}t \cos^{n}t (1-k^{2}\sin^{2} t)^{r/2}
  t:   limit of integration
  S:   \sin t
  C:   \cos t

  note: r = \pm 1; 0\leq k < 1 (not yet implemented for more general case),
  otherwise no restriction on input parameters.
 */

{
  gdouble f, I, tol = 1e-16, df, SC, S2 = S*S, k2 = k*k, kc = 0.2 ;
  gint q ;

  g_assert(k >= 0.0) ;
  if ( k >= 1.0 ) return 0.0 ;
  
  g_assert(r == 1 || r == -1 || r == -3) ;

  if ( r == -3 ) {
    gdouble d = sqrt(1-k2*S2) ;
    gdouble kd2 = 1.0 - k2 ;

    /*some of these are hard-wired*/
    if ( m == 1 && n == 0 ) return -C/kd2/d ;      
    if ( m == 0 && n == 1 ) return  S/d ;
    if ( m == 0 && n == 3 ) return -kd2*S/k2/d + 1.0/k2/k*asin(k*S) ;
    if ( m == 0 && n == 5 ) return  ((-k2*S2+2*k2*k2-4*k2+3)/2./d*S + 
				     (4*k2-3.)/2./k*asin(k*S))/k2/k2 ;
    if ( m == 1 && n == 4 ) 
      return ((k2*S2 + 2*k2 - 3.0)/2.0/d*C + 3*kd2/2.0/k*log(k*C+d))/k2/k2 ;

    if ( m == 2 && n == 3 )
      return ((k2*S2 + 2*k2 - 3.0)/2/d*S - (2*k2 - 3)/2/k*asin(k*S))/k2/k2 ;

    f = pow(S,m-1)*pow(C,n-1)/k2/d ;
    f -= (gdouble)(m-1)/k2*quad_sin_cos_nmr(m-2, n  , -1, k, t, S, C) ;
    f += (gdouble)(n-1)/k2*quad_sin_cos_nmr(m  , n-2, -1, k, t, S, C) ;

    return f ;
  }

  /*check for hard-wired cases if k is not too small*/
  if ( k > kc ) {
    if ( (quad_sin_cos_nmr_hw(m, n, r, k, t, S, C, &I)) == 0 )
      return I ;
  }

  f = I = quad_sin_cos(m, n, t, S, C) ;
  SC = pow(S,m+1)*pow(C,n+1) ;

  q = 1 ;
  if ( m+2*q+n == 0) 
    I = quad_sin_cos(m+2*q, n, t, S, C) ;
  else
    I = -SC/(m+2*q+n) + I*(gdouble)(m+2*q-1)/(m+n+2*q) ;

  if ( r == 1 ) {
    df = -0.5*k2 ; f += df*I ;
    
    for ( q = 2 ; (q < SERIES_MAX_TERMS) && (fabs(df*I) > tol) ; q ++ ) {
      SC *= S2 ;
      df *= k2*(gdouble)(2*q-3)/(2*q) ;
      if ( m+2*q+n == 0) 
	I = quad_sin_cos(m+2*q, n, t, S, C) ;
      else
	I = -SC/(m+2*q+n) + I*(gdouble)(m+2*q-1)/(m+n+2*q) ;
      f += df*I ;
    }

    return f ;
  }

  df = 0.5*k2 ; f += df*I ;
    
  for ( q = 2 ; (q < SERIES_MAX_TERMS) && (fabs(df*I) > tol) ; q ++ ) {
    SC *= S2 ;
    df *= k2*(gdouble)(2*q-1)/(2*q) ;
    if ( m+2*q+n == 0)
      I = quad_sin_cos(m+2*q, n, t, S, C) ;
    else
      I = -SC/(m+2*q+n) + I*(gdouble)(m+2*q-1)/(m+n+2*q) ;
    f += df*I ;
  }

  return f ;
}
