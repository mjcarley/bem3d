/* qselect.c
 * 
 * Copyright (C) 2006, 2018 Michael Carley
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

/**
 * @defgroup qselect Quadrature selection
 *
 * BEM3D has a built-in method for selecting quadrature rules (which
 * can be overridden using a \link config configuration file \endlink.
 * The selection is performed using a parameter \f$\sigma\f$,
 * calculated from geometrical properties of the element, and of the
 * point. For an element, define
 * \f$\mathbf{c}=\overline{\mathbf{x}_{i}}\f$, the mean of the element
 * vertices. Then \f$R_{s}\f$ is defined as the maximum value of
 * \f$\sqrt{2}|\mathbf{c}-\mathbf{x}_{i}|\f$, so that a sphere
 * containing the element is defined. The point-element distance is
 * defined as \f$R_{p}=|\mathbf{c}-\mathbf{x}|\f$. For field points
 * which coincide with an element vertex \f$\sigma=0\f$. For
 * \f$R_{p}<R_{s}\f$, \f$\sigma=z_{n}/R_{s}\f$ where \f$z_{n}\f$ is
 * the approximate normal distance from the field point to the plane
 * of the element. For \f$R_{p}\geq R_{s}\f$,
 * \f$\sigma=R_{p}/R_{s}\f$.
 *
 * @{
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>

#include <gts.h>
#include <glib.h>

#include "bem3d.h"
#include "bem3d-private.h"

/** 
 * Allocate a new BEM3DQuadratureSelector
 * 
 * @return pointer to new selector
 */

BEM3DQuadratureSelector *bem3d_quadrature_selector_new(void)

{
  BEM3DQuadratureSelector *s ;
  
  g_debug("%s: ", __FUNCTION__) ;

  s = (BEM3DQuadratureSelector *)
    g_malloc(sizeof(BEM3DQuadratureSelector)) ;
  
  s->f = g_ptr_array_new() ;
  s->sigma = g_array_new(FALSE, TRUE, sizeof(gdouble)) ;
  s->NM = g_array_new(FALSE, TRUE, sizeof(gint)) ;

  return s ;
}

/** 
 * Clear a quadrature selection rule
 * 
 * @param s quadrature selector
 * 
 * @return 0 on success
 */

gint bem3d_quadrature_selector_clear(BEM3DQuadratureSelector *s)

{
  g_return_val_if_fail(s != NULL, BEM3D_NULL_ARGUMENT) ;

  g_debug("%s: ", __FUNCTION__) ;

  g_ptr_array_set_size(s->f, 0) ;
  g_array_set_size(s->sigma, 0) ;
  g_array_set_size(s->NM, 0) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Add a new quadrature to a selection rule
 * 
 * @param s ::BEM3DQuadratureSelector;
 * @param f ::BEM3DQuadratureRuleFunc to add.
 * @param p parameter, the quadrature will be selected if
 * bem3d_quadrature_parameter returns a value greater than p
 * @param N number of points in angle in polar rules
 * @param M number of points in radius in polar rules
 *
 * \a N and \a M can be used for different purposes in other rules,
 * but will always be passed as an array of int of length 2, [N M]
 * 
 * @return 0 on success
 */

gint bem3d_quadrature_selector_add(BEM3DQuadratureSelector *s,
				   BEM3DQuadratureRuleFunc f,
				   gdouble p, 
				   gint N, gint M)

{
  g_debug("%s: ", __FUNCTION__) ;

  g_return_val_if_fail(s != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p >= 0.0, BEM3D_ARGUMENT_OUT_OF_RANGE) ;

  g_ptr_array_add(s->f, f) ;
  g_array_append_val(s->sigma, p) ;
  g_array_append_val(s->NM, N) ;
  g_array_append_val(s->NM, M) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Select a quadrature for a given parameter
 * 
 * @param s BEM3DQuadratureSelector
 * @param p quadrature parameter, returned by bem3d_quadrature_parameter
 * @param f BEM3DQuadratureRuleFunc to call
 * @param data to pass to f
 * 
 * @return 
 */

gint bem3d_quadrature_select(BEM3DQuadratureSelector *s,
			     gdouble p,
			     BEM3DQuadratureRuleFunc *f,
			     gpointer *data)
{
  gint i ;
  
  g_debug("%s: s=%p; p=%lg; f=%p; data=%p", 
	  __FUNCTION__, s, p, f, data) ;

  g_return_val_if_fail(s != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(p >= 0.0, BEM3D_ARGUMENT_OUT_OF_RANGE) ;

  for ( i = 0 ; 
	(i < bem3d_quadrature_selector_length(s)-1) &&
	  bem3d_quadrature_selector_sigma(s,i) < p ;
	i++ ) ;

  *f = bem3d_quadrature_selector_rule(s,i) ;
  *data = &(bem3d_quadrature_selector_data(s,i)) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Default quadrature selection rule which is reasonably good for most
 * problems. For points close to, or on, an element, polar or Hayami's
 * PART quadrature is used, switching to a symmetric high-order rule
 * further away.
 * 
 * @return pointer to quadrature selector 
 */

BEM3DQuadratureSelector *bem3d_quadrature_selector_default()

{
  BEM3DQuadratureSelector *s = NULL ;

  g_debug("%s: ", __FUNCTION__) ;

  s = bem3d_quadrature_selector_new() ;
  bem3d_quadrature_selector_add(s, bem3d_quadrature_rule_polar, 0.0, 64, 64) ;
  bem3d_quadrature_selector_add(s, bem3d_quadrature_rule_polar, 0.125, 32, 32) ;
  bem3d_quadrature_selector_add(s, bem3d_quadrature_rule_polar, 0.25, 16, 16) ;
  bem3d_quadrature_selector_add(s, bem3d_quadrature_rule_polar, 0.5, 16, 16) ;
  bem3d_quadrature_selector_add(s, bem3d_quadrature_rule_wx, 1.0, 54, 0) ;
  bem3d_quadrature_selector_add(s, bem3d_quadrature_rule_wx, 2.0, 25, 0) ;
  /* bem3d_quadrature_selector_add(s, bem3d_quadrature_rule_wx, 4.0, 7, 0) ; */
  /* bem3d_quadrature_selector_add(s, bem3d_quadrature_rule_gauss, 8.0, 1, 0) ; */

  return s ;
}

/**
 * @}
 * 
 */
