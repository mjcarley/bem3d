#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include "fmm.h"

gint main(gint argc, gchar **argv)

{
  gint ier, iprec, nsource, ntarget ;
  gboolean ifcharge, ifdipole, ifpot, iffld, ifpottarg, iffldtarg ;
  gdouble *source, *charge, *dipstr, *dipvec, *pot, *fld, *target,
    *pottarg, *fldtarg ;
  _fmm_c zk ;
  gint i, j, k ;
  FILE *output ;

  ier = 6 ; iprec = 13 ;
  prini_(&ier, &iprec) ;

  ifcharge = TRUE ;
  ifdipole = FALSE ; 
  ifpot = TRUE ; 
  iffld = FALSE ; 
  ifpottarg = FALSE ; 
  iffldtarg = FALSE ; 

  zk.r = 1.0 ; zk.i = 0.3 ;
  
  nsource = 100000 ;
  source = (gdouble *)g_malloc(3*nsource*sizeof(gdouble)) ;
  charge = (gdouble *)g_malloc(2*nsource*sizeof(gdouble)) ;
  dipstr = (gdouble *)g_malloc(2*nsource*sizeof(gdouble)) ;
  dipvec = (gdouble *)g_malloc(3*nsource*sizeof(gdouble)) ;
  pot = (gdouble *)g_malloc(2*nsource*sizeof(gdouble)) ;
  fld = (gdouble *)g_malloc(6*nsource*sizeof(gdouble)) ;

  ntarget = 10 ;
  pottarg = (gdouble *)g_malloc(ntarget*sizeof(gdouble)) ;
  fldtarg = (gdouble *)g_malloc(3*ntarget*sizeof(gdouble)) ;

  for ( i = 0 ; i < nsource ; i ++ ) {
    source[3*i+0] = drand48() ;
    source[3*i+1] = drand48() ;
    source[3*i+2] = drand48() ;
    charge[2*i+0] = i+1 ;
    charge[2*i+1] = -i+1 ;
  }

  ier = 0 ; iprec = 3 ; ntarget = 0 ;
  hfmm3dparttarg_(&ier, &iprec, &zk,
		  &nsource, source,
		  &ifcharge, charge,
		  &ifdipole, dipstr, dipvec,
		  &ifpot, pot, 
		  &iffld, fld,
		  &ntarget, target,
		  &ifpottarg, pottarg,
		  &iffldtarg, fldtarg) ;

  output = fopen("field.dat", "w") ;

  for ( i = 0 ; i < nsource ; i ++ ) {
    fprintf(output, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",
	    source[3*i+0], source[3*i+1], source[3*i+2],
	    charge[2*i+0], charge[2*i+1],
	    pot[2*i+0], pot[2*i+1]) ;
  }

  fclose(output) ;

  return 0 ;
}
