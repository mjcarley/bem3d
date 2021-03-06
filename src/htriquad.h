#ifndef BEM3D_HTRIQUAD_H_INCLUDED
#define BEM3D_HTRIQUAD_H_INCLUDED

#include <glib.h>

#define SIGN(_x) ((_x) < 0 ? -1 : 1)

#define _htri_vector_cross(_p, _x, _y)		\
  do {						\
  _p[0] = _x[1]*_y[2] - _x[2]*_y[1] ;		\
  _p[1] = _x[2]*_y[0] - _x[0]*_y[2] ;		\
  _p[2] = _x[0]*_y[1] - _x[1]*_y[0] ;		\
  } while (0) 

#define _htri_invert3x3(_Ai,_A)						\
  do {                                                                  \
    gdouble _det = 1.0/((_A)[0]*((_A)[8]*(_A)[4]-(_A)[7]*(_A)[5]) -     \
                        (_A)[3]*((_A)[8]*(_A)[1]-(_A)[7]*(_A)[2]) +     \
                        (_A)[6]*((_A)[5]*(_A)[1]-(_A)[4]*(_A)[2])) ;    \
    (_Ai)[0] =  _det*((_A)[8]*(_A)[4] - (_A)[7]*(_A)[5]) ;              \
    (_Ai)[1] = -_det*((_A)[8]*(_A)[1] - (_A)[7]*(_A)[2]) ;              \
    (_Ai)[2] =  _det*((_A)[5]*(_A)[1] - (_A)[4]*(_A)[2]) ;              \
                                                                        \
    (_Ai)[3] = -_det*((_A)[8]*(_A)[3] - (_A)[6]*(_A)[5]) ;              \
    (_Ai)[4] =  _det*((_A)[8]*(_A)[0] - (_A)[6]*(_A)[2]) ;              \
    (_Ai)[5] = -_det*((_A)[5]*(_A)[0] - (_A)[3]*(_A)[2]) ;              \
                                                                        \
    (_Ai)[6] =  _det*((_A)[7]*(_A)[3] - (_A)[6]*(_A)[4]) ;              \
    (_Ai)[7] = -_det*((_A)[7]*(_A)[0] - (_A)[6]*(_A)[1]) ;              \
    (_Ai)[8] =  _det*((_A)[4]*(_A)[0] - (_A)[3]*(_A)[1]) ;              \
  } while (0)



gint htri_quad_triangle_1_0(gdouble r1, gdouble r2, gdouble th,
			    gdouble k, gdouble tol,
			    gdouble *d0G, gdouble *d1G, gdouble *d2G) ;
gint htri_tri_z_parameters(gdouble r1, gdouble r2, gdouble th, 
			   gdouble phi, gdouble s, gdouble z,
			   gdouble *S, gdouble *Rmax, gdouble *al) ;
gint quad_delta_C_n(gdouble al, gdouble t0, gdouble t1, gint N,
		    gdouble *Q) ;
gint quad_delta_C_n_T(gdouble al, gdouble t0, gdouble t1, gint N,
		      gdouble *Q) ;
gint quad_delta_C_n_alpha(gdouble al, gdouble *Q, gint n,
			  gdouble *I0, gdouble *I1, gdouble *I2, 
			  gdouble *I3) ;
gint quad_Rzq_cos_sin(gdouble s, gdouble z, gdouble th, gdouble phi,
		      gdouble al, 
		      gdouble *Iq, gdouble *Iqt, 
		      gint Q, 
		      gdouble *I0, gint zstr,
		      gdouble *Ic, gint cstr,
		      gdouble *Is, gint sstr) ;

gint htri_tri_parameters(gdouble r1, gdouble r2, gdouble th, 
			 gdouble *phi, gdouble *s) ;
gint htri_select_expansion(gdouble dx, gdouble tol,
			   gdouble **coefficients_cos, gint *nc,
			   gdouble **coefficients_sin, gint *ns) ;
gint htri_quad_triangle_0(gdouble r1, gdouble r2, gdouble th, gdouble z,
			  gdouble k, gdouble tol,
			  gdouble *G, gdouble *dG, gdouble *d2G) ;
gint htri_quad_triangle_0_xy(gdouble r1, gdouble r2, gdouble th, gdouble z,
			     gdouble k, gdouble tol,
			     gdouble *d0G, gdouble *d1G, gdouble *d2G) ;
gint htri_quad_triangle_1(gdouble r1, gdouble r2, gdouble th, gdouble z,
			  gdouble k, gdouble tol,
			  gdouble *d0G, gdouble *d1G, gdouble *d2G) ;

gint htri_triangle_decompose(gdouble *x1, gdouble *x2, gdouble *x3,
			     gdouble *x,
			     gdouble psi[], gdouble th[],
			     gdouble r1[], gdouble r2[],
			     gint sgn[],
			     gint *ntri) ;
gint htri_triangle_project(gdouble *og, gdouble *s, gdouble *t, gdouble *n,
			   gdouble *x, gdouble *y) ;
gint htri_triangle_axes(gdouble *x1, gdouble *x2, gdouble *x3,
			gdouble *x,
			gdouble *og, gdouble *s, gdouble *t, gdouble *n) ;
gint htri_quad_shape_1(gdouble *xf, 
		       gdouble *x1, gdouble *x2, gdouble *x3, 
		       gdouble k, gdouble tol, gint qmax,
		       gdouble *g, gdouble *dg, gdouble *d2g) ;
gint htri_quad_shape_0(gdouble *xf, 
		       gdouble *x1, gdouble *x2, gdouble *x3, 
		       gdouble k, gdouble tol, gint qmax,
		       gdouble *g, gdouble *dg, gdouble *d2g) ;

gint htri_quad_order_1R(gdouble rmin, gdouble rmax, gdouble z,
			gint qmax, gdouble tol) ;

#endif /*BEM3D_HTRIQUAD_H_INCLUDED*/
