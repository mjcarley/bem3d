/* intsincos.c
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

#define FUNCTION_NOT_IMPLEMENTED return -1 

gint quad_sin_cos_nmr_hw(gint m, gint n, gint r, gdouble k,
			 gdouble t, gdouble S, gdouble C, gdouble *I) ;

gint quad_sin_cos_nmr_hw(gint m, gint n, gint r, gdouble k,
			 gdouble t, gdouble S, gdouble C, gdouble *I)

/* hand-crafted analytical evaluation of the most common integrals,
   from Gradshteyn and Ryzhik, section 2.58 */

/*
 *I = \int \sin^{m} x\cos^{n}x \Delta^{r} dx, \Delta^{2}=1-k^{2}\sin^{2}x
 
 t:    \theta
 S:    \sin\theta
 C:    \cos\theta
 
 return value: 0 on success; -1 if the integral is not implemented.
 */

{
  gdouble k2 = k*k, kd2 = 1.0-k2, kd = sqrt(kd2), S2 = S*S, C2 = C*C ;
  gdouble d = sqrt(1.0-k2*S2) ;

  g_assert(r == 1 || r == -1) ;

  if ( r == 1 ) {
    switch ( m ) {
    default: break ;
    case -5:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -1: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  0: FUNCTION_NOT_IMPLEMENTED ;
	*I = ((k2-3)*S2+2.0)/8.0/S2/S2*C*d +
	  kd2*(k2+3.0)/16.0*log((d+C)/(d-C)) ;
	return 0 ;
	break ;
      case  1: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case -4:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -1: 
	*I = -((3-k2)*S2 + 1)/3/S2/S*d - kd/2*log((d-kd*S)/(d+kd*S)) ;
	return 0 ;
	break ;
      case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  1: 
	*I = -d*d*d/3/S2/S ;
	return 0 ;
	break ;
      case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case -3:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -2:
	*I = (3*S2-1)/2/S2/C*d - (k2-3)/4*log((d-C)/(d+C)) ;
	return 0 ;
	break ;
      case -1:
	*I = -d/2/S2 + kd/2.0*log((d+kd)/(d-kd)) + (k2-2)/4*log((1+d)/(1-d)) ;
	return 0 ;
	break ;
      case  0: 
	*I = -d*C/2/S2 - kd2/4*log((d+C)/(d-C)) ;
	return 0 ;
	break ;
      case  1: 
	*I = -d*0.5/S2 + 0.25*k2*log((1+d)/(1-d)) ;
	return 0 ;
	break ;
      case  2: 
	*I = -C/2/S2*d + (k2+1)/4*log((d+C)/(d-C)) - k*log(k*C+d) ;
	return 0 ;
	break ;
      case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case -2:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -3:
	*I = (3*S2 - 2)/2/S/C2*d - 
	  (2*k2-3)/4/kd*log((d+kd*S)/(d-kd*S)) ;
	return 0 ;
	break ;
      case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -1: FUNCTION_NOT_IMPLEMENTED ;
	*I = -d/S - (1.0+k2)/2.0/kd*log((d-kd*S)/(d+kd*S)) ;
	return 0 ;
	break ;
      case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  1: 
	*I = -d/S - k*asin(k*S) ;
	return 0 ;
	break ;
      case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  3:
	*I = -(S2+2)/2/S*d - (2*k2+1)/2/k*asin(k*S) ;
	return 0 ;
	break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case -1:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: 
	*I = ((2*k2-3)*S2 - 3*k2 + 4)/3/kd2/C2/C*d - 
	  0.5*log((d+C)/(d-C)) ;
	return 0 ;
	break ;
      case -3: 
	*I = d/2/C2 - 0.5*log((1+d)/(1-d)) + (2-k2)/4/kd*log((d+kd)/(d-kd)) ;
	return 0 ;
	break ;
      case -2:
	*I = d/C - 0.5*log((d+C)/(d-C)) ;
	return 0 ;
	break ;
      case -1: 
	*I = 0.5*log((1-d)/(1+d)) + kd/2*log((d+kd)/(d-kd)) ;
	return 0 ;
	break ;
      case  0: 
	*I = -0.5*log((d+C)/(d-C)) + k*log(k*C+d) ;
	return 0 ;
	break ;
      case  1:
	*I = d + 0.5*log((1-d)/(1+d)) ;
	return 0 ;
	break ;
      case  2:
	*I = d*C/2 + (k2+1)/2/k*log(k*C+d) - 0.5*log((d+C)/(d-C)) ;
	return 0 ;
	break ;
      case  3: 
	*I = -(k2*S2 - 3*k2 -1)/3/k2*d + 0.5*log((1-d)/(1+d)) ;
	return 0 ;
	break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ;
	*I = (-2*k2*S2 + 5*k2 + 1)/8.0/k2*C*d + 
	  0.5*log((d+C)/(d-C)) + (3*k2*k2+6*k2-1)/8/k2/k*log(k*C+d) ;
	return 0 ;
	break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case 0:
      switch ( n ) {
      case -5:
	*I = ((2.0*k2-3.0)*S2 - 4*k2+5)/(8*kd2*C2*C2)*S*d
	  - (4*k2-3)/16/kd2/kd*log((d + kd*S)/(d-kd*S)) ;
	return 0 ;
	break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ;	break ;
      case -3:
	*I = d*S/2/C2 + 1.0/4/kd*log((d + kd*S)/(d-kd*S)) ;
	return 0 ;
	break ;
      case -2: FUNCTION_NOT_IMPLEMENTED ;	break ;
      case -1:
	*I = kd/2*log((d + kd*S)/(d-kd*S)) + k*asin(k*S) ;
	return 0 ;
	break ;
      case 0: FUNCTION_NOT_IMPLEMENTED ;break ;
      case 1:
	*I = d*S/2 + 0.5/k*asin(k*S) ;
	return 0 ;
	break ;
      case 2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case 3:
	*I = (2*k2*C2 + 2*k2 + 1)/8/k2*d*S + (4*k2-1)/8/k2/k*asin(k*S) ;
	return 0 ;
	break ;
      case 4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case 5:
	*I = 
	  (8*k2*k2*S2*S2 - 2*k2*(12*k2+1)*S2 + 24*k2*k2 + 
	   12*k2 -3)/48/k2/k2*d*S +
	  (8*k2*k2-4*k2+1)/16/k2/k2/k*asin(k*S) ;
	return 0 ;
	break ;
      }
    case 1:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4:
	FUNCTION_NOT_IMPLEMENTED ;
	*I = (-(2*k2+1)*k2*S2 + 3*k2*k2 - k2 + 1)/3/kd2/C2/C*d ;
	return 0 ;
	break ;
      case -3:
	*I = d/2/C2 + k2/4/kd*log((d+kd)/(d-kd)) ;
	return 0 ;
	break ;
      case -2:
	*I = d/C - k*log(k*C + d) ;
	return 0 ;
	break ;
      case -1: 
	*I = -d + kd/2*log((d+kd)/(d-kd)) ;
	return 0 ;
	break ;
      case  0: 
	*I = -d*C/2 - kd2/2/k*log(k*C + d) ;
	return 0 ;
	break ;
      case  1: 
	*I = -d*d*d/3/k2 ; 
	return 0 ;
	break ;
      case  2: 
	*I = -(2*k2*C2 + kd2)/8/k2*d*C + kd2*kd2/8/k/k2*log(k*C+d) ;
	return 0 ;
	break ;
      case  3:
	*I = -(3*k2*k2*S2*S2 - k2*(5*k2+1)*S2 + 5*k2-2)/15/k2/k2*d ;
	return 0 ;
	break ;
      case  4: 
	*I = (-8*k2*k2*S2*S2 + 2*k2*(7*k2+1)*S2 - 
	      3*k2*k2 - 8*k2 +3)/48/k2/k2*d*C -
	  kd2*kd2*kd2/16/k2/k2/k*log(k*C+d) ;
	return 0 ;
	break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case 2:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -3:
	*I = S/2/C2*d + (2*k2-1)/4/kd*log((d+kd*S)/(d-kd*S)) - k*asin(k*S) ;
	return 0 ;
	break ;
      case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -1: 
	*I = -d*S/2 + (2*k2-1)/2/k*asin(k*S) + kd/2*log((d+kd*S)/(d-kd*S)) ;
	return 0 ;
	break ;
      case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  1:
	*I = (2*k2*S2 - 1)/8/k2*d*S + 1.0/8/k2/k*asin(k*S) ;
	return 0 ;
	break ;
      case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case 3:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -2: 
	*I = -(S2-3)/2/C*d - (3*k2-1)/2/k*log(k*C+d) ;
	return 0 ;
	break ;
      case -1: 
	*I = -(k2*S2 + 3*k2 -1)/3/k2*d + kd/2*log((d+kd)/(d-kd)) ;
	return 0 ;
	break ;
      case  0:
	*I = -(2*k2*S2 + 3*k2-1)/8/k2*d*C + 
	  (3*k2*k2 - 2*k2 -1)/8/k2/k*log(k*C+d) ;
	return 0 ;
	break ;
      case  1:
	*I = (3*k2*k2*S2*S2 - k2*S2-2)/15/k2/k2*d ;
	break ;
      case  2: 
	*I = (8*k2*k2*S2*S2 - 2*k2*(k2+1)*S2 - 
	      3*k2*k2 + 2*k2 -3)/48/k2/k2*d*C + 
	  kd2*kd2*(k2+1)/16/k2/k2/k*log(k*C+d) ;
	return 0 ;
	break ;
      case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case 4:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -1: 
	*I = -(2*k2*S2 + 4*k2-1)/8/k2*d*S + 
	  (8*k2*k2 - 4*k2 -1)/8/k2/k*asin(k*S) +
	  kd/2*log((d+kd*S)/(d-kd*S)) ;
	return 0 ;
	break ;
      case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  1: 
	*I = (8*k2*k2*S2*S2 - 2*k2*S2 -3)/48/k2/k2*d*S + 
	  1.0/16/k2/k2/k*asin(k*S) ;
	return 0 ;
	break ;
      case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    case 5:
      switch (n) {
      case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case -1: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  0: 
	*I = (-8*k2*k2*S2*S2 - 2*k2*(5*k2-1)*S2 - 
	      15*k2*k2 + 4*k2 +3)/48/k2/k2*d*C + 
	  (5*k2*k2*k2 - 3*k2*k2 - k2 - 1)/16/k2/k2/k*log(k*C+d) ;
	return 0 ;
	break ;
      case  1: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
      case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
      }
      break ;
    }
    return -1 ;
  }

  /*r == -1*/
  switch ( m ) {
  default: break ;
  case -5:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  0: 
      *I = -(3*(k2+1)*S2+2)/8/S2*d*C + 
	(3*k2*k2+2*k2+3)/16*log((d+C)/(d-C)) ;
      return 0 ;
      break ;
    case  1: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case -4:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1: 
      *I = -((3+2*k2)*S2 + 1)/3/S2/S*d - log((d-kd*S)/(d+kd*S))/2/kd ;
      return 0 ;
      break ;
    case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  1: 
      *I = -(2*k2*S2+1)/3/S2/S*d ;
      return 0 ;
      break ;
    case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case -3:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -2: 
      *I = ((3-k2)*S2-kd2)/2/kd2/S2/C*d + (k2+3)/4*log((d-C)/(d+C)) ;
      return 0 ;
      break ;
    case -1: 
      *I = -d/2/S2 + log((d+kd)/(d-kd))/2/kd - (k2+2)/4*log((1+d)/(1-d)) ;
      return 0 ;
      break ;
    case  0: 
      *I = -d/2*C/S2 - (k2+1)/4*log((d+C)/(d-C)) ;
      return 0 ;
      break ;
    case  1: FUNCTION_NOT_IMPLEMENTED ;
      *I = -d/2/S2 - k2/d*log((1+d)/(1-d)) ;
      return 0 ;
      break ;
    case  2: 
      *I = -d*C/S2/2 + kd2/4*log((d+C)/(d-C)) ;
      return 0 ;
      break ;
    case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case -2:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: 
      *I = ((3-2*k2)*S2-2*kd2)/2/kd2/S/C2*d - 
	(4*k2-3)/4/kd2/kd*log((d+kd*S)/(d-kd*S)) ;
      return 0 ;
      break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1: 
      *I = -d/S - log((d-kd*S)/(d+kd*S))/2/kd ;
      return 0 ;
      break ;
    case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  1: 
      *I = -d/S ;
      return 0 ;
      break ;
    case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  3: 
      *I = -d/S - asin(k*S)/k ;
      return 0 ;
      break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case -1:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: 
      *I = ((5*k2-3)*S2-6*k2+4)/3/kd2/kd2/C2/C*d - 
	0.5*log((d+C)/(d-C)) ;
      return 0 ;
      break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ;
      *I = -d/2/kd2/C2 - 0.5*log((1+d)/(1-d)) + 
	(2.0-3*k2)/4/kd2/kd*log((d+kd)/(d-kd)) ;
      return 0 ;
      break ;
    case -2: 
      *I = d/kd2/C + 0.5*log((d-C)/(d+C)) ;
      return 0 ;
      break ;
    case -1: 
      *I = 0.5*log((1-d)/(1+d)) + log((d+kd)/(d-kd))/2/kd ;
      return 0 ;
      break ;
    case  0: 
      *I = -0.5*log((d+C)/(d-C)) ;
      return 0 ;
      break ;
    case  1: 
      *I = 0.5*log((1-d)/(1+d)) ;
      return 0 ;
      break ;
    case  2: 
      *I = -0.5*log((d+C)/(d-C)) + log(k*C+d)/k ;
      return 0 ;
      break ;
    case  3: 
      *I = d/k2 - 0.5*log((1+d)/(1-d)) ;
      return 0 ;
      break ;
    case  4: 
      *I = d*C/2/k2 - 0.5*log((d+C)/(d-C)) + (3*k2-1)/2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case  0:
    switch (n) {
    case -5: 
      *I = (3*(2*k2-1)*S2-8*k2+5)/8/kd2/kd2/C2/C2*d*S + 
	(8*k2*k2-8*k2+3)/16/kd2/kd2/kd*log((d+kd*S)/(d-kd*S)) ;
      return 0 ;
      break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: 
      *I = d*S/2/kd2/C2 + (2*k2-1)/4/kd2/kd*log((d-kd*S)/(d+kd*S)) ;
      return 0 ;
      break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1:
      *I = -0.5/kd*log((d-kd*S)/(d+kd*S)) ;
      return 0 ;
      break ;
    case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  1: 
      *I = asin(k*S)/k ;
      return 0 ;
      break ;
    case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  3: 
      *I = S*d/2/k2 + (2*k2-1)/2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: 
      *I = (2*k2*C2 + 6*k2-3)/8/k2/k2*S*d + 
	(8*k2*k2-8*k2+3)/8/k2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    case  6: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  7:
      *I = (8*k2*k2*S2*S2 - 2*k2*(18*k2-5)*S2 + 72*k2*k2 - 
	    54*k2+15)/48/k2/k2/k2*S*d + 
	(16*k2*k2*k2 - 24*k2*k2 + 18*k2 - 5)/16/k2/k2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    }
    break ;
  case  1:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ;
      *I = -(2*k2*C2-kd2)/2/kd2/kd2/C2/C*d ;
      return 0 ;
      break ;
    case -3: 
      *I = d/2/kd2/C2 - k2/4/kd2/kd*log((d+kd)/(d-kd)) ;
      return 0 ;
      break ;
    case -2:
      *I = d/kd2/C ;
      return 0 ;
      break ;
    case -1: 
      *I = log((d+kd)/(d-kd))/2/kd ;
      return 0 ;
      break ;
    case  0: 
      *I = -log(k*C+d)/k ;
      return 0 ;
      break ;
    case  1: 
      *I = -d/k2 ;
      return 0 ;
      break ;
    case  2: 
      *I = -C*d/2/k2 + kd2/2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  3:  FUNCTION_NOT_IMPLEMENTED ;
      *I = -(k2*C2-2.0*kd2)/3.0/k2/k2*d ;
      return 0 ;
      break ;
    case  4: 
      *I = (3-5*k2+2*k2*S2)/8/k2/k2*C*d - 
	(3*k2*k2-6*k2+3)/8/k2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  5: 
      *I = (-3*k2*k2*C2*C2 + 4*k2*kd2*C2 - 8*k2*k2 + 
	    16*k2 -8)/15/k2/k2/k2*d ;
      return 0 ;
      break ;
    case  6: 
      *I = (-8*k2*k2*S2*S2 + 2*k2*(13*k2-5)*S2-
	    33*k2*k2+40*k2-15)/48/k2/k2/k2*C*d + 
	5*kd2*kd2*kd2/16/k2/k2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  7: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case  2:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: 
      *I = d*S/2/kd2/C2 - log((d+kd*S)/(d-kd*S))/4/kd2/kd ;
      return 0 ;
      break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1: 
      *I = log((d+kd*S)/(d-kd*S))/2/kd - asin(k*S)/k ;
      return 0 ;
      break ;
    case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  1: 
      *I = -S*d/2/k2 + asin(k*S)/2/k2/k ;
      return 0 ;
      break ;
    case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  3: 
      *I = -(2*k2*C2+2*k2-3)/8/k2/k2*S*d + (4*k2-3)/8/k2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: 
      *I = (-8*k2*k2*S2*S2 + 2*k2*(12*k2-5)*S2 -24*k2*k2 +
	    36*k2-15)/48/k2/k2/k2*S*d +
	(8*k2*k2-12*k2+5)/16/k2/k2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    case  6: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  7: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case  3:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -2: 
      *I = d/kd2/C + log(k*C+d)/k ;
      return 0 ;
      break ;
    case -1: 
      *I = d/k2 + log((d+kd)/(d-kd))/2/kd ;
      return 0 ;
      break ;
    case  0: 
      *I = C*d/2/k2 - (1+k2)/2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  1:
      *I = -(2+k2*S2)*d/3/k2/k2 ;
      break ;
    case  2:
      *I = (2*k2*C2-k2-3)/8/k2/k2*C*d - (k2*k2+2*k2-3)/8/k2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  3: 
      *I = (3*k2*k2*S2*S2 - (5*k2*k2-4*k2)*S2-10*k2+8)/15/k2/k2/k2*d ;
      return 0 ;
      break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ;
      *I = (8*k2*k2*S2*S2-2*k2*(6*k2-5)*S2+3*k2*k2-22*k2+15)/48/k2/k2/k2*C*d -
	(k2*k2*k2+3*k2*k2-9*k2+5)/16/k2/k2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  6: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  7: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case  4:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1: 
      *I = d*S/2/k2 + log((d+kd*S)/(d-kd*S))/2/kd - (2*k2+1)/2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  1: 
      *I = -(2*k2*S2+3)/8/k2/k2*S*d + 3.0/8/k2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  3: 
      *I = (8*k2*k2*S2*S2 - 2*k2*(6*k2-5)*S2-18*k2+15)/48/k2/k2/k2*S*d +
	(6*k2-5)/16/k2/k2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  6: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  7: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  case  5:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  0: 
      *I = (2*k2*S2+3*k2+3)/8/k2/k2*C*d - 
	(3+2*k2+3*k2*k2)/8/k2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  1: 
      *I = -(3*k2*k2*S2*S2 + 4*k2*S2+8)/15/k2/k2/k2*d ;
      return 0 ;
      break ;
    case  2: 
      *I = (-8*k2*k2*S2*S2 + 2*k2*(k2-5)*S2+3*k2*k2+
	    4*k2-15.0)/48/k2/k2/k2*C*d -
	(k2*k2*k2+k2*k2+3*k2-5)/16/k2/k2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  6: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  7: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
  case  6:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  0: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  1: 
      *I = -(8*k2*k2*S2*S2+10*k2*S2+15)/48/k2/k2/k2*S*d + 
	5.0/16/k2/k2/k2/k*asin(k*S) ;
      return 0 ;
      break ;
    case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  6: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  7: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
  case  7:
    switch (n) {
    case -5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case -1: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  0: 
      *I = (8*k2*k2*S2*S2 + 10*k2*(k2+1)*S2+15*k2*k2+
	    14*k2+15)/48/k2/k2/k2*C*d -
	(5*k2*k2-2*k2+5)*(k2+1)/16/k2/k2/k2/k*log(k*C+d) ;
      return 0 ;
      break ;
    case  1: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  2: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  3: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  4: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  5: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  6: FUNCTION_NOT_IMPLEMENTED ; break ;
    case  7: FUNCTION_NOT_IMPLEMENTED ; break ;
    }
    break ;
  }

  return -1 ;
}

