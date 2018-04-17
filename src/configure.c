/* configure.c
 * 
 * Copyright (C) 2010, 2017 by Michael Carley
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
 * @defgroup config Configuring BEM3D
 *
 * BEM3D codes can be configured at run-time using glib key files and
 * the following functions.
 *
 * @{
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <gts.h>
#include <glib.h>

#include "bem3d.h"
#include "bem3d-private.h"

GHashTable *idmap = NULL, *varmap = NULL ;

/*pre-defined structs for Green's functions*/
BEM3DGreensFunction greens_func_laplace = 
  {bem3d_greens_func_laplace, TRUE, 1} ;
BEM3DGreensFunction greens_func_gradient_laplace = 
  {bem3d_greens_func_gradient_laplace, TRUE, 4} ;
BEM3DGreensFunction greens_func_helmholtz = 
  {bem3d_greens_func_helmholtz, FALSE, 1} ;
BEM3DGreensFunction greens_func_helmholtz_hs = 
  {bem3d_greens_func_helmholtz_hs, FALSE, 1} ;
BEM3DGreensFunction greens_func_gradient_helmholtz = 
  {bem3d_greens_func_gradient_helmholtz, FALSE, 4} ;

/** 
 * Allocate a new ::BEM3DConfiguration, filled with default settings
 * (Laplace, ::bem3d_quadrature_rule_default, direct solver).
 * 
 * @return newly allocated ::BEM3DConfiguration
 */

BEM3DConfiguration *bem3d_configuration_new()

{
  BEM3DConfiguration *c ;

  c = (BEM3DConfiguration *)g_malloc(sizeof(BEM3DConfiguration)) ;

  c->job = g_string_new("BEM3D") ;
  c->keys = g_ptr_array_new() ;

  bem3d_greens_function_func(&(c->gfunc)) = bem3d_greens_func_laplace ;
  bem3d_greens_function_is_real(&(c->gfunc)) = TRUE ;
  c->qrule = bem3d_quadrature_rule_default ;
  c->qdata = bem3d_quadrature_selector_default() ;

  c->solver = BEM3D_SOLVER_DIRECT ;
  c->fmm    = BEM3D_FMM_FMMLIB3D_1_2 ;

  c->skel_order = 1 ;
  c->fmm_radius = 1.0 ;
  c->fmm_tol    = 1e-6 ;
  c->bc_default_admittance[0] = c->bc_default_admittance[1] = 0 ;

  return c ;
}

static void parse_quadrature_key(GKeyFile *key, gchar *group, 
				 gchar *k, BEM3DQuadratureSelector *s)
  
{
  gchar *v, **t ;
  GError *error = NULL ;
  gint N = 0, M = 0 ;
  gdouble sigma ;
  BEM3DQuadratureRuleFunc f ;

  v = g_key_file_get_value(key, group, k, &error) ;
  if ( error != NULL ) g_error("%s: %s", __FUNCTION__, error->message) ;

  t = g_strsplit(g_strdelimit(v, "(), ", '|'), "|", 0) ;

  f = bem3d_configuration_pointer_from_identifier(t[0]) ;

  if ( f == NULL ) 
    g_error("%s: unrecognized identifier %s in config file", 
	    __FUNCTION__, t[0]) ;
  if ( t[1] == NULL ) 
    g_error("%s: no selection parameter specified for quadrature rule %s",
	    __FUNCTION__, t[0]) ;

  errno = 0 ; sigma = strtod(t[1], NULL) ;
  if ( errno != 0 ) 
    g_error("%s: could not convert selection parameter %s for "
	    "quadrature rule %s", __FUNCTION__, t[1], t[0]) ;

  if ( t[2] == NULL ) {
    bem3d_quadrature_selector_add(s, f, sigma, N, M) ;
    return ;
  }

  errno = 0 ; N = strtol(t[2], NULL, 0) ; 
  if ( errno != 0 ) 
    g_error("%s: could not convert selection parameter %s "
	    "for quadrature rule %s", __FUNCTION__, t[2], t[0]) ;
  if ( t[3] == NULL ) {
    bem3d_quadrature_selector_add(s, f, sigma, N, M) ;
    return ;
  }

  errno = 0 ; M = strtol(t[3], NULL, 0) ; 
  if ( errno != 0 ) 
    g_error("%s: could not convert selection parameter %s "
	    "for quadrature rule %s", __FUNCTION__, t[3], t[0]) ;
  bem3d_quadrature_selector_add(s, f, sigma, N, M) ;


  return ;
}

static void quadrature_rules_set(GKeyFile *key, gchar *group, 
				 BEM3DQuadratureSelector *s)

{
  gint i, nk ;
  gchar **w ;
  GError *error = NULL ;

  if ( !g_key_file_has_group(key, "BEM3D::Quadrature") ) return ;

  w = g_key_file_get_keys(key, "BEM3D::Quadrature", (gsize *)(&nk), &error) ;
  if ( error != NULL ) g_error("%s: %s", __FUNCTION__, error->message) ;

  if ( nk == 0 || w == NULL ) return ;

  bem3d_quadrature_selector_clear(s) ;
  for ( i = 0 ; i < nk ; i ++ )
    parse_quadrature_key(key, "BEM3D::Quadrature", w[i], s) ;

  /*check that the parameters are ordered*/
  if ( bem3d_quadrature_selector_sigma(s,0) != 0.0 ) 
    g_error("%s: first quadrature selection parameter must be 0", 
	    __FUNCTION__) ;
  for ( i = 1 ; i < bem3d_quadrature_selector_length(s); i ++ ) 
    if ( !(bem3d_quadrature_selector_sigma(s,i) >
	   bem3d_quadrature_selector_sigma(s,i-1)) )
      g_error("%s: quadrature selection parameters must be in ascending order", 
	      __FUNCTION__) ;
      

  return ;
}

static void physics_set(GKeyFile *key, BEM3DConfiguration *c)

/*set Green's function and other physical data for problem*/

{
  gchar *v ;
  GError *error = NULL ;
  BEM3DGreensFunction *p ;
  gint ret ;

  if ( !g_key_file_has_group(key, "BEM3D::Physics") ) return ;

  if ( g_key_file_has_key(key, "BEM3D::Physics", "greens_function", &error) ) {
    /*set Green's function*/
    v = g_key_file_get_value(key, "BEM3D::Physics", "greens_function", NULL) ;
    if ( (p = bem3d_configuration_pointer_from_identifier(v)) == NULL ) 
      g_error("%s: unrecognized Green's function (%s) in key file", 
	      __FUNCTION__, v) ;
    c->gfunc = *p ;
  }

  if ( g_key_file_has_key(key, "BEM3D::Physics", "surface_admittance", 
			  &error) ) {
    /*default surface admittance for boundary conditions*/
    v = g_key_file_get_value(key, "BEM3D::Physics", "surface_admittance", 
			     NULL) ;
    if ( (ret = parse_complex(v, c->bc_default_admittance)) != 0 )
      g_error("%s: cannot parse default admittance \"%s\"", __FUNCTION__, v) ;  
  }

  return ;
}

static void general_set(GKeyFile *key, BEM3DConfiguration *c)

/*set Green's function and other physical data for problem*/

{
  gchar *v ;
  GError *error = NULL ;
  gpointer p ;

  if ( !g_key_file_has_group(key, "BEM3D::Solver") ) return ;

  if ( g_key_file_has_key(key, "BEM3D::Solver", "solver", &error) ) {
    /*set solver*/
    v = g_key_file_get_value(key, "BEM3D::Solver", "solver", NULL) ;
    if ( (p = bem3d_configuration_pointer_from_identifier(v)) == NULL ) 
      g_error("%s: unrecognized solver (%s) in key file", 
	      __FUNCTION__, v) ;
    c->solver = GPOINTER_TO_INT(p) ;
  }

  error = NULL ;
  if ( g_key_file_has_key(key, "BEM3D::Solver", "fmm", &error) ) {
    /*set solver*/
    v = g_key_file_get_value(key, "BEM3D::Solver", "fmm", NULL) ;
    if ( (p = bem3d_configuration_pointer_from_identifier(v)) == NULL ) 
      g_error("%s: unrecognized fast multipole method solver (%s) in key file", 
	      __FUNCTION__, v) ;
    c->fmm = GPOINTER_TO_INT(p) ;
  }

  if ( g_key_file_has_key(key, "BEM3D::Solver", "tolerance_fmm", &error) ) {
    /*set solver*/
    v = g_key_file_get_value(key, "BEM3D::Solver", "tolerance_fmm", NULL) ;
    c->fmm_tol = atof(v) ;
  }

  if ( g_key_file_has_key(key, "BEM3D::Solver", "radius_fmm", &error) ) {
    /*set solver*/
    v = g_key_file_get_value(key, "BEM3D::Solver", "radius_fmm", NULL) ;
    c->fmm_radius = atof(v) ;
  }

  if ( g_key_file_has_key(key, "BEM3D::Solver", "skeleton_order", &error) ) {
    /*set solver*/
    v = g_key_file_get_value(key, "BEM3D::Solver", "skeleton_order", NULL) ;
    c->skel_order = atoi(v) ;
  }

  return ;
}

/** 
 * Set fields of a ::BEM3DConfiguration from a key file, overwriting
 * previous settings. Settings are read from the "BEM3D" group of the
 * file, to allow codes to keep all their settings in one
 * place. Values not in the "BEM3D" group are ignored.
 * 
 * @param c a ::BEM3DConfiguration;
 * @param k a GKeyFile containing settings for \a c.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_configuration_set(BEM3DConfiguration *c, gpointer k) 

{
#if GLIB_CHECK_VERSION(2,14,0)
  GKeyFile *key = k ;

  if ( !g_key_file_has_group(key, "BEM3D") &&
       !g_key_file_has_group(key, "BEM3D::Quadrature")
       ) return BEM3D_SUCCESS ;

  /*quadrature rules*/
  quadrature_rules_set(key, "BEM3D::Quadrature", c->qdata) ;

  /*Green's function and other physical data*/
  physics_set(key, c) ;

  /*general solver settings*/
  general_set(key, c) ;

#else
  g_message("%s: key file support not available with this version of GLIB",
	    __FUNCTION__) ;
#endif /*GLIB_CHECK_VERSION(2,14,0)*/

  return BEM3D_SUCCESS ;
}

/** 
 * Read BEM3D configuration data from a key file.
 * 
 * @param c a ::BEM3DConfiguration
 * @param file a GLIB format key file containing configuration data
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_configuration_read(BEM3DConfiguration *c, gchar *file)

{
#if GLIB_CHECK_VERSION(2,14,0)
  GKeyFile *key ;
  GError *error ;
  GKeyFileFlags flags = G_KEY_FILE_KEEP_COMMENTS |
    G_KEY_FILE_KEEP_TRANSLATIONS ;

  key = g_key_file_new() ;

  error = NULL ;
  if ( !(g_key_file_load_from_file(key, file, flags, &error)) ) {
    g_error("%s: %s, reading %s", __FUNCTION__, error->message, file) ;
  }

  bem3d_configuration_set(c, key) ;			  

  g_key_file_free(key) ;  
#else
  g_message("%s: key file support not available with this version of GLIB",
	    __FUNCTION__) ;
#endif /*GLIB_CHECK_VERSION(2,14,0)*/

  return BEM3D_SUCCESS ;
}

/** 
 * Initialize internal data for ::BEM3DConfiguration. This function
 * should be called before any configurations are set or used, i.e. at
 * the start of any program.
 * 
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_configuration_init(void)

{
  idmap = g_hash_table_new(g_str_hash, g_str_equal) ;
  varmap = g_hash_table_new(NULL, NULL) ;

  /*add the built-in functions*/
  /*quadrature*/
  bem3d_configuration_add_identifier("bem3d_quadrature_newman", 
				     bem3d_quadrature_rule_newman) ;
  bem3d_configuration_add_identifier("bem3d_quadrature_newman_gradient", 
				     bem3d_quadrature_rule_newman_gradient) ;
  bem3d_configuration_add_identifier("bem3d_quadrature_gauss",
  				     bem3d_quadrature_rule_gauss) ;
  bem3d_configuration_add_identifier("bem3d_quadrature_khayat_wilton", 
				     bem3d_quadrature_rule_kw) ;
  bem3d_configuration_add_identifier("bem3d_quadrature_polar", 
				     bem3d_quadrature_rule_polar) ;
  bem3d_configuration_add_identifier("bem3d_quadrature_polar_hypersingular", 
				     bem3d_quadrature_rule_polar_hs) ;
  bem3d_configuration_add_identifier("bem3d_quadrature_wandzura_xiao", 
				     bem3d_quadrature_rule_wx) ;
  bem3d_configuration_add_identifier("bem3d_quadrature_hayami", 
				     bem3d_quadrature_rule_hayami) ;

  /*Green's functions*/
  bem3d_configuration_add_identifier("bem3d_greens_function_laplace", 
				     &greens_func_laplace) ;
  bem3d_configuration_add_identifier("bem3d_greens_function_gradient_laplace", 
				     &greens_func_gradient_laplace) ;
  bem3d_configuration_add_identifier("bem3d_greens_function_helmholtz", 
				     &greens_func_helmholtz) ;
  bem3d_configuration_add_identifier("bem3d_greens_function_gradient_helmholtz", 
				     &greens_func_gradient_helmholtz) ;
  bem3d_configuration_add_identifier("bem3d_greens_function_helmholtz_hypersingular", 
				     &greens_func_helmholtz_hs) ;
  
  /*solver settings*/
  bem3d_configuration_add_identifier("bem3d_fmm_fmmlib3d1.2",
				     GINT_TO_POINTER(BEM3D_FMM_FMMLIB3D_1_2)) ;
  bem3d_configuration_add_identifier("bem3d_solver_direct",
				     GINT_TO_POINTER(BEM3D_SOLVER_DIRECT)) ;
  bem3d_configuration_add_identifier("bem3d_solver_fmm",
				     GINT_TO_POINTER(BEM3D_SOLVER_FMM)) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Add an identifier to the configuration data. This allows the
 * identifier to be recognized by the bem3d_configuration_
 * functions. For example,
 *
 * bem3d_configuration_add_identifier("bem3d_quadrature_polar",
 * bem3d_quadrature_rule_polar) ;
 * 
 * adds the identifier "bem3d_quadrature_polar" which refers to the
 * quadrature function of the same name, allowing a key file to
 * contain a configuration line of the form
 *
 * quadrature[1] = bem3d_quadrature_polar(0.0,16,16)
 *
 * @param id name of identifier;
 * @param v a pointer to the identified object
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_configuration_add_identifier(const gchar *id, gpointer v)

{
  gboolean rt ;

  g_assert(idmap != NULL && varmap != NULL) ;

  rt = (g_hash_table_lookup(idmap, id) != NULL) ;
  if ( rt ) g_error("%s: value %s already in idmap", __FUNCTION__, id) ; 
  g_hash_table_insert(idmap, (gpointer)id, v) ;

  rt = (g_hash_table_lookup(varmap, id) != NULL) ;
  if ( rt ) g_error("%s: value %s already in varmap", __FUNCTION__, id) ;
  g_hash_table_insert(varmap, v, (gpointer)id) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Find the identifier associated with a given object, set using
 * ::bem3d_configuration_add_identifier
 * 
 * @param v pointer to an object
 * 
 * @return the identifier associated with \a v
 */

gchar *bem3d_configuration_identifier_from_pointer(gpointer v)

{
  if ( idmap == NULL ) return NULL ;

  return g_hash_table_lookup(varmap, v) ;
}

/** 
 * Look up the object associated with a given configuration
 * identifier, set using ::bem3d_configuration_add_identifier
 * 
 * @param id identifier used in configuration settings;
 * 
 * @return a pointer to the object associated with \a id
 */

gpointer bem3d_configuration_pointer_from_identifier(const gchar *id)

{
  gpointer p ;

  if ( idmap == NULL ) return NULL ;

  p = g_hash_table_lookup(idmap, id) ;

  return p ;
}

/**
 * @}
 * 
 */