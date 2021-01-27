/* wbfmm-bem3d.h
 * 
 * Copyright (C) 2019 Michael Carley
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

#ifndef BEM3D_WBFMM_H_INCLUDED
#define BEM3D_WBFMM_H_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <glib.h>

gint _bem3d_fmm_helmholtz_wbfmm(BEM3DFastMultipole solver,
				BEM3DFastMultipoleProblem problem,
				BEM3DParameters *param,
				BEM3DMeshSkeleton *s, 
				gdouble tol,
				gdouble *q, gdouble *dq,
				gdouble *p, gdouble *dp,
				BEM3DFMMWorkspace *w) ;
gint _bem3d_fmm_laplace_wbfmm(BEM3DFastMultipole solver,
			      BEM3DFastMultipoleProblem problem,
			      BEM3DParameters *param,
			      BEM3DMeshSkeleton *s, 
			      gdouble tol,
			      gdouble *q, gdouble *dq,
			      gdouble *p, gdouble *dp,
			      BEM3DFMMWorkspace *w) ;

#endif /*BEM3D_WBFMM_H_INCLUDED*/
