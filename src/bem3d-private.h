#ifndef _BEM_PRIVATE_H_INCLUDED_
#define _BEM_PRIVATE_H_INCLUDED_

#include <glib.h>
#include <gts.h>
#include "bem3d.h"

#define BEM3D_STARTUP_MESSAGE \
"This is free software; see the source code for copying conditions.\n\
There is ABSOLUTELY NO WARRANTY, not even for MERCHANTABILITY or\n\
FITNESS FOR A PARTICULAR PURPOSE.\n"

#define BEM3D_BMESH_DATA_WIDTH   32
#define BEM3D_BMESH_DATA_MESH     0
#define BEM3D_BMESH_DATA_N        1 
#define BEM3D_BMESH_DATA_EDGES    2
#define BEM3D_BMESH_DATA_N_EDGES  3
#define BEM3D_BMESH_DATA_HASH     4
#define BEM3D_BMESH_DATA_ELEMENT  5
#define BEM3D_BMESH_DATA_QFUNC    6
#define BEM3D_BMESH_DATA_QDATA    7
#define BEM3D_BMESH_DATA_EFUNC    9
#define BEM3D_BMESH_DATA_EDATA   10
#define BEM3D_BMESH_DATA_GFUNC   11
#define BEM3D_BMESH_DATA_GDATA   12
#define BEM3D_BMESH_DATA_NODE    13
#define BEM3D_BMESH_DATA_ROW     14
#define BEM3D_BMESH_DATA_V       15
#define BEM3D_BMESH_DATA_BFUNC   16
#define BEM3D_BMESH_DATA_BDATA   17
#define BEM3D_BMESH_DATA_IMIN    18
#define BEM3D_BMESH_DATA_IMAX    19
#define BEM3D_BMESH_DATA_LFUNC   20
#define BEM3D_BMESH_DATA_LDATA   21
#define BEM3D_BMESH_DATA_POINT   22
#define BEM3D_BMESH_DATA_F       23
#define BEM3D_BMESH_DATA_RFUNC   24
#define BEM3D_BMESH_DATA_RDATA   25
#define BEM3D_BMESH_DATA_ANGLE   26
#define BEM3D_BMESH_DATA_TREE    27
#define BEM3D_BMESH_DATA_FIELD   28
#define BEM3D_BMESH_DATA_CONFIG  29

#ifndef g_debug
#define g_debug(format...) g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, format)
#endif
#ifndef g_message
#define g_message(format...) g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, format)
#endif
#ifndef g_warning
#define g_warning(format...) g_log(G_LOG_DOMAIN, G_LOG_LEVEL_WARNING, format)
#endif
#ifndef g_error
#define g_error(format...) g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, format)
#endif

extern BEM3DGreensFunction greens_func_laplace, greens_func_helmholtz ;

/* GtsEdge *connect_vertices(GtsVertex *v1, GtsVertex *v2) ; */
void rotation_indices(gint nc, gint rot[]) ;
gint gauss_select_quadrature(guint npts, 
			     const gdouble **xk, const gdouble **wt,
			     quadrature_rule_t quad_type) ;
gint gauss_set_scaling(gdouble x0, gdouble x1, gdouble *xbar, gdouble *dx,
		       guint npts, quadrature_rule_t quad_type) ;

gint newman_tri_shape(gdouble p[], gdouble x1[], gdouble x2[], gdouble x3[],
		      gdouble Imn[], gint hmax,
		      gdouble I[], gdouble J[]) ;
gint newman_tri_quad(gdouble p[], gdouble x1[], gdouble x2[], gdouble x3[],
		     gdouble I[], gdouble J[]) ;
gint newman_tri_moments(gdouble x1[], gdouble x2[],
			gdouble x3[], gint hmax, gdouble Imn[]) ;
gint newman_tri_shape_gradient(gdouble p[], gdouble x1[], gdouble x2[], 
			       gdouble x3[], gdouble Imn[], gint hmax,
			       gdouble I[], gdouble J[]) ;

gint triangle_quad_shape(gdouble *x1, gdouble *x2, gdouble *x3,
			 gdouble *p, 
			 gdouble o12, gdouble o23, gdouble o31,
			 gint order, gboolean hs,
			 gdouble *G1, gdouble *G3, gdouble *G5) ;
gdouble quad_sin_cos_nmr(gint m, gint n, gint r, gdouble k,
			 gdouble t, gdouble S, gdouble C) ;

FILE *file_open(gchar *fname, gchar *namedefault, gchar *mode, 
		FILE *fdefault) ;
gint file_close(FILE *f) ;

gint parse_complex(gchar *v, gdouble z[]) ;

#endif /*_BEM_PRIVATE_H_INCLUDED_*/