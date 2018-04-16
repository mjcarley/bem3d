/* fmmlib3d_1_2.c
 * 
 * Copyright (C) 2017 Michael Carley
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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <glib.h>
#include <gts.h>

#include "bem3d.h"
#include "bem3d-private.h"

#ifdef HAVE_FMMLIB3D_1_2

#include "fmmlib3d_1_2.h"




void hfmm3dparttarg_(gint *ier, gint *iprec, _fmmlib12_c *z,
		     gint *nsource, gdouble *source,
		     gint *ifcharge, gdouble *charge,
		     gint *ifdipole, gdouble *dipstr, gdouble *dipvec,
		     gint *ifpot, gdouble *pot, 
		     gint *iffld, gdouble *fld,
		     gint *ntarget, gdouble *target,
		     gint *ifpottarg, gdouble *pottarg,
		     gint *iffldtarg, gdouble *fldtarg) ;

void lfmm3dparttarg_(gint *ier, gint *iprec,
		     gint *nsource, gdouble *source,
		     gint *ifcharge, gdouble *charge,
		     gint *ifdipole, gdouble *dipstr, gdouble *dipvec,
		     gint *ifpot, gdouble *pot, 
		     gint *iffld, gdouble *fld,
		     gint *ntarget, gdouble *target,
		     gint *ifpottarg, gdouble *pottarg,
		     gint *iffldtarg, gdouble *fldtarg) ;

void prini_(gint *ip1, gint *iq1) ;

#else /*HAVE_FMMLIB3D_1_2*/
#endif /*HAVE_FMMLIB3D_1_2*/

static gint complex_fill(gdouble *dest, gdouble *src, gint n)

{
  gint i ;

  for ( i = 0 ; i < n ; i ++ ) {
    dest[2*i+0] = src[i] ; dest[2*i+1] = 0.0 ;
  }

  return 0 ;
}

static gint real_fill(gdouble *dest, gdouble *src, gint n)

{
  gint i ;

  for ( i = 0 ; i < n ; i ++ ) dest[i] = src[2*i+0] ;

  return 0 ;
}

gint fmmlib3d_1_2_set_prec(gdouble tol)

{
  gdouble prec[] = {0.5, 0.5e-1, 0.5e-2, 0.5e-3, 0.5e-6, 0.5e-9, 0.5e-12, 
		    0.5e-15, 0} ;
  gint iprec ;

  for ( iprec = 0 ; (tol < prec[iprec]) && ( iprec < 9 ) ; iprec ++ ) ;
  
  return iprec-2 ;
}


gint _bem3d_fmm_helmholtz_fmmlib3d_1_2(BEM3DFastMultipole solver,
				       BEM3DFastMultipoleProblem problem,
				       BEM3DParameters *param,
				       BEM3DMeshSkeleton *s, 
				       gdouble tol,
				       gdouble *q, gdouble *dq,
				       gdouble *p, gdouble *dp,
				       BEM3DFMMWorkspace *w)
  
{
  gint ier, iprec ;
  gboolean ifcharge, ifdipole, ifpot, iffld, ifpottarg, iffldtarg ;
  gdouble *charge, *dipstr, *pot, *fld, *pottarg, *fldtarg ;
  _fmmlib12_c zk, dum[3] ;

  g_assert(solver == BEM3D_FMM_FMMLIB3D_1_2) ;
  g_assert(problem == BEM3D_FMM_HELMHOLTZ) ;
  
  ifcharge = FALSE ;  charge  = dum ;
  ifdipole = FALSE ;  dipstr  = dum ;
  ifpot    = FALSE ;  pot     = dum ;
  iffld    = FALSE ;  fld     = dum ;
  ifpottarg = FALSE ; pottarg = dum ;
  iffldtarg = FALSE ; fldtarg = dum ;

  if (  q != NULL ) { ifcharge = TRUE ; charge = q ; }
  if ( dq != NULL ) { ifdipole = TRUE ; dipstr = dq ; }

  if (  p != NULL ) { ifpottarg = TRUE ; pottarg = p ; }
  if ( dp != NULL ) { iffldtarg = TRUE ; fldtarg = dp ; }

  zk.r = bem3d_parameters_wavenumber(param) ; zk.i = 0.0 ;

  iprec = fmmlib3d_1_2_set_prec(tol) ;

#ifdef HAVE_FMMLIB3D_1_2

  ier = 0 ;
  hfmm3dparttarg_(&ier, &iprec, &zk,
		  &(s->ns), s->x,
		  &ifcharge, charge,
		  &ifdipole, dipstr, s->n,
		  &ifpot, dum,
		  &iffld, dum,
		  &(s->nt), &(s->x[3*(s->ns)]),
		  &ifpottarg, pottarg,
		  &iffldtarg, fldtarg) ;
  g_assert(ier == 0) ;

#else /*HAVE_FMMLIB3D_1_2*/
  g_error("%s: no fmmlib3d-1.2 support compiled in", __FUNCTION__) ;
#endif /*HAVE_FMMLIB3D_1_2*/
  
  return 0 ;
}

gint _bem3d_fmm_laplace_fmmlib3d_1_2(BEM3DFastMultipole solver,
				     BEM3DFastMultipoleProblem problem,
				     BEM3DParameters *param,
				     BEM3DMeshSkeleton *s, 
				     gdouble tol,
				     gdouble *q, gdouble *dq,
				     gdouble *p, gdouble *dp,
				     BEM3DFMMWorkspace *w)
  
/*param can be NULL here (no wavenumber to set)*/

{
  gint ier, iprec ;
  gboolean ifcharge, ifdipole, ifpot, iffld, ifpottarg, iffldtarg ;
  gdouble *charge, *dipstr, *pot, *fld, *pottarg, *fldtarg, dum[6] ;
  gdouble *buf ;

  g_assert(solver == BEM3D_FMM_FMMLIB3D_1_2) ;
  g_assert(problem == BEM3D_FMM_LAPLACE) ;

  ifcharge = FALSE ;  charge  = dum ;
  ifdipole = FALSE ;  dipstr  = dum ;
  ifpot = FALSE ;     pot     = dum ;
  iffld = FALSE ;     fld     = dum ;
  ifpottarg = FALSE ; pottarg = dum ;
  iffldtarg = FALSE ; fldtarg = dum ;

  buf = w->d ;
  if (  q != NULL ) { 
    ifcharge = TRUE ; charge = buf ;
    complex_fill(charge, q, s->ns) ;
    buf = &(charge[2*(s->ns)]) ;
  }

  if ( dq != NULL ) { 
    ifdipole = TRUE ; dipstr = buf ; 
    complex_fill(dipstr, dq, s->ns) ;
    buf = &(dipstr[2*(s->ns)]) ;
  }

  if (  p != NULL ) { 
    ifpottarg = TRUE ; pottarg = buf ;
    buf = &(pottarg[2*(s->nt)]) ;
  }

  if ( dp != NULL ) { 
    iffldtarg = TRUE ; fldtarg = buf ; 
    buf = &(fldtarg[6*(s->ns)]) ;
  }

  iprec = fmmlib3d_1_2_set_prec(tol) ;

#ifdef HAVE_FMMLIB3D_1_2

  ier = 0 ;
  lfmm3dparttarg_(&ier, &iprec,
		  &(s->ns), s->x,
		  &ifcharge, charge,
		  &ifdipole, dipstr, s->n,
		  &ifpot, dum,
		  &iffld, dum,
		  &(s->nt), &(s->x[3*(s->ns)]),
		  &ifpottarg, pottarg,
		  &iffldtarg, fldtarg) ;
  g_assert(ier == 0) ;

#else /*HAVE_FMMLIB3D_1_2*/
  g_error("%s: no fmmlib3d-1.2 support compiled in", __FUNCTION__) ;
#endif /*HAVE_FMMLIB3D_1_2*/

  if (  p != NULL ) { real_fill(p, pottarg, s->nt) ; }

  if (  dp != NULL ) { real_fill(dp, fldtarg, s->nt) ; }
  
  return 0 ;
}
