#ifndef __FMM_H_INCLUDED__
#define __FMM_H_INCLUDED__

#include <glib.h>

typedef struct {gdouble r, i ;} _fmm_c ;

void hfmm3dparttarg_(gint *ier, gint *iprec, _fmm_c *z,
		     gint *nsource, gdouble *source,
		     gint *ifcharge, gdouble *charge,
		     gint *ifdipole, gdouble *dipstr, gdouble *dipvec,
		     gint *ifpot, gdouble *pot, 
		     gint *iffld, gdouble *fld,
		     gint *ntarget, gdouble *target,
		     gint *ifpottarg, gdouble *pottarg,
		     gint *iffldtarg, gdouble *fldtarg) ;

void prini_(gint *ip1, gint *iq1) ;

#endif /*__FMM_H_INCLUDED__*/
