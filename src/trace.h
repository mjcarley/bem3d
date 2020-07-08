#ifndef BEM3D_TRACE_H_INCLUDED
#define BEM3D_TRACE_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

extern gboolean _bem3d_trace[64] ;

#define _bem3d_set_trace(_i) (_bem3d_trace[(_i)] = TRUE)
#define _bem3d_unset_trace(_i) (_bem3d_trace[(_i)] = FALSE)
#define _bem3d_trace_set(_i)   (_bem3d_trace[(_i)])

#endif /*BEM3D_TRACE_H_INCLUDED*/
