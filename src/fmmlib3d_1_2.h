/* fmmlib3d_1_2.c
 * 
 * Copyright (C) 2017, 2018 Michael Carley
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

#ifndef __FMMLIB3D_1_2_H_INCLUDED__
#define __FMMLIB3D_1_2_H_INCLUDED__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <glib.h>

typedef struct {gdouble r, i ;} _fmmlib12_c ;

gint _bem3d_fmm_helmholtz_fmmlib3d_1_2(BEM3DFastMultipole solver,
				       BEM3DFastMultipoleProblem problem,
				       BEM3DParameters *param,
				       BEM3DMeshSkeleton *s, 
				       gdouble tol,
				       gdouble *q, gdouble *dq,
				       gdouble *p, gdouble *dp,
				       BEM3DFMMWorkspace *w) ;
gint _bem3d_fmm_laplace_fmmlib3d_1_2(BEM3DFastMultipole solver,
				     BEM3DFastMultipoleProblem problem,
				     BEM3DParameters *param,
				     BEM3DMeshSkeleton *s, 
				     gdouble tol,
				     gdouble *q, gdouble *dq,
				     gdouble *p, gdouble *dp,
				     BEM3DFMMWorkspace *w) ;

#endif /*__FMMLIB3D_1_2_H_INCLUDED__*/
