#ifndef BEM3D_POLAR_H_INCLUDED
#define BEM3D_POLAR_H_INCLUDED

gdouble orient2d(gdouble *, gdouble *, gdouble *) ;
gint triangle_orientations(gdouble *x1, gdouble *x2, gdouble *x3,
			   gdouble *p, gdouble *o12, gdouble *o23, 
			   gdouble *o31) ;
gint triangle_rotation_matrix(gdouble *x1, gdouble *x2, gdouble *x3,
			      gdouble *A) ;
void matrix_transform(gdouble *A, gdouble *x, gdouble *p, gdouble *y) ;
gdouble quad_sin_cos(gint m, gint n, gdouble t, gdouble S, gdouble C) ;

gdouble quad_sin_cos_nmr(gint m, gint n, gint r, gdouble k,
			 gdouble t, gdouble S, gdouble C) ;
gint quad_sin_cos_nmr_hw(gint m, gint n, gint r, gdouble k,
			 gdouble t, gdouble S, gdouble C, gdouble *I) ;

#define matvecmul3x3(_y,_A,_x)				\
  (((_y)[0] = (_A)[0]*(_x)[0] + (_A)[1]*(_x)[1] + (_A)[2]*(_x)[2]),	\
   ((_y)[1] = (_A)[3]*(_x)[0] + (_A)[4]*(_x)[1] + (_A)[5]*(_x)[2]),	\
   ((_y)[2] = (_A)[6]*(_x)[0] + (_A)[7]*(_x)[1] + (_A)[8]*(_x)[2]))

#define matvecmul3x3trans(_y,_A,_x)				\
  (((_y)[0] = (_A)[0]*(_x)[0] + (_A)[3]*(_x)[1] + (_A)[6]*(_x)[2]),	\
   ((_y)[1] = (_A)[1]*(_x)[0] + (_A)[4]*(_x)[1] + (_A)[7]*(_x)[2]),	\
   ((_y)[2] = (_A)[2]*(_x)[0] + (_A)[5]*(_x)[1] + (_A)[8]*(_x)[2]))

#define invert3x3(_Ai,_A)				\
  do {									\
    gdouble _det = 1.0/((_A)[0]*((_A)[8]*(_A)[4]-(_A)[7]*(_A)[5]) -	\
			(_A)[3]*((_A)[8]*(_A)[1]-(_A)[7]*(_A)[2]) +	\
			(_A)[6]*((_A)[5]*(_A)[1]-(_A)[4]*(_A)[2])) ;	\
    (_Ai)[0] =  _det*((_A)[8]*(_A)[4] - (_A)[7]*(_A)[5]) ;		\
    (_Ai)[1] = -_det*((_A)[8]*(_A)[1] - (_A)[7]*(_A)[2]) ;		\
    (_Ai)[2] =  _det*((_A)[5]*(_A)[1] - (_A)[4]*(_A)[2]) ;		\
									\
    (_Ai)[3] = -_det*((_A)[8]*(_A)[3] - (_A)[6]*(_A)[5]) ;		\
    (_Ai)[4] =  _det*((_A)[8]*(_A)[0] - (_A)[6]*(_A)[2]) ;		\
    (_Ai)[5] = -_det*((_A)[5]*(_A)[0] - (_A)[3]*(_A)[2]) ;		\
									\
    (_Ai)[6] =  _det*((_A)[7]*(_A)[3] - (_A)[6]*(_A)[4]) ;		\
    (_Ai)[7] = -_det*((_A)[7]*(_A)[0] - (_A)[6]*(_A)[1]) ;		\
    (_Ai)[8] =  _det*((_A)[4]*(_A)[0] - (_A)[3]*(_A)[1]) ;		\
  } while (0)

gint triangle_gradient_quad_shape(gdouble *x1, gdouble *x2, gdouble *x3,
				  gdouble *p,
				  gdouble o12, gdouble o23, gdouble o31,
				  gdouble *G, gdouble *dG) ;

#endif /*BEM3D_POLAR_H_INCLUDED*/
