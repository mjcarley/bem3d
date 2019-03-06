/* configure.c
 * 
 * Copyright (C) 2010, 2017, 2018 by Michael Carley
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
  {bem3d_greens_func_laplace, TRUE, 0, 1, 1} ;
BEM3DGreensFunction greens_func_gradient_laplace = 
  {bem3d_greens_func_gradient_laplace, TRUE, 0, 4, 4} ;
BEM3DGreensFunction greens_func_helmholtz = 
  {bem3d_greens_func_helmholtz, FALSE, 0, 1, 2} ;
BEM3DGreensFunction greens_func_helmholtz_hs = 
  {bem3d_greens_func_helmholtz_hs, FALSE,
   BEM3D_GREENS_FUNC_PRECOMPUTE_NORMAL |
   BEM3D_GREENS_FUNC_PRECOMPUTE_LAMBDA,
   1, 2} ;
BEM3DGreensFunction greens_func_helmholtz_ch = 
  {bem3d_greens_func_helmholtz_ch, FALSE,
   BEM3D_GREENS_FUNC_PRECOMPUTE_NORMAL |
   BEM3D_GREENS_FUNC_PRECOMPUTE_LAMBDA,
   3, 6} ;
BEM3DGreensFunction greens_func_gradient_helmholtz = 
  {bem3d_greens_func_gradient_helmholtz, FALSE, 0, 4, 8} ;

#define _CONFIG_STRIDE  4
#define _CONFIG_NAME    0
#define _CONFIG_SETTING 1
#define _CONFIG_DESC    2
#define _CONFIG_POINTER 3

/*
  configuration settings and descriptions, four entries per setting:

  [variable name] [setting to which variable can be assigned] 
  [setting description] [pointer to setting data]

  if pointer to setting data is NULL, the setting is handled within a
  function rather than in the hash tables
*/

static gpointer _config_quadratures[] =
  {"bem3d_quadrature_newman",
   "quadrature[]",
   "analytical quadrature method for linear shape functions on "
   "planar elements (for Laplace equation)",
   bem3d_quadrature_rule_newman,
   "bem3d_quadrature_newman_gradient",
   "quadrature[]",
   "analytical quadrature method for linear shape functions on "
   "planar elements (for evaluation of gradient of Laplace potential)",   
   bem3d_quadrature_rule_newman_gradient,
   "bem3d_quadrature_gauss",
   "quadrature[]",
   "Gaussian quadrature",
   bem3d_quadrature_rule_gauss,
   "bem3d_quadrature_khayat_wilton",
   "quadrature[]",
   "Khayat and Wilton's method for integration of 1/R potentials",
   bem3d_quadrature_rule_kw,
   "bem3d_quadrature_polar",
   "quadrature[]",
   "polar transformation on planar elements",
   bem3d_quadrature_rule_polar,
   "bem3d_quadrature_polar_hypersingular",
   "quadrature[]",
   "polar transformation on planar elements for hypersingular integrand "
   "(do not use)",
   bem3d_quadrature_rule_polar_hs,
   "bem3d_quadrature_wandzura_xiao",
   "quadrature[]",
   "Wandzura and Xiao's high-order symmetric rules for triangles",
   bem3d_quadrature_rule_wx,
   "bem3d_quadrature_hayami",
   "quadrature[]",
   "Hayami's method for near-singularities",
   bem3d_quadrature_rule_hayami,
   "bem3d_quadrature_series",
   "quadrature[]",
   "Carley's method for semi-analytical evaluation of Helmholtz potentials",
   bem3d_quadrature_rule_series,
   "bem3d_quadrature_mzht",
   "quadrature[]",
   "experimental method for hypersingular integrands",
   bem3d_quadrature_rule_mzht,
   NULL, NULL, NULL, NULL} ;

static gpointer _config_physics[] =
  {"bem3d_greens_function_laplace",
   "greens_function", "",
   &greens_func_laplace,

   "bem3d_greens_function_gradient_laplace", 
   "greens_function", "",
   &greens_func_gradient_laplace,

   "bem3d_greens_function_helmholtz", 
   "greens_function", "",
   &greens_func_helmholtz,

   "bem3d_greens_function_gradient_helmholtz", 
   "greens_function", "",
   &greens_func_gradient_helmholtz,

   "bem3d_greens_function_helmholtz_hyper", 
   "greens_function", "",
   &greens_func_helmholtz_hs,

   "bem3d_greens_function_helmholtz_chen_harris", 
   "greens_function", "",
   &greens_func_helmholtz_ch,
   
   "[real or complex admittance]",
   "surface_admittance",
   "default surface admittance, multiplied by wavenumber in Helmholtz problems",
   NULL,

   NULL, NULL, NULL, NULL} ;

static gpointer _config_solver[] =
  {"bem3d_fmm_fmmlib3d1.2",
   "fmm", "fast multipole solver of Greengard et al.",
   GINT_TO_POINTER(BEM3D_FMM_FMMLIB3D_1_2),

   "bem3d_solver_direct",
   "solver", "direct solution using full matrices",
   GINT_TO_POINTER(BEM3D_SOLVER_DIRECT),

   "bem3d_solver_fmm",
   "solver", "solution using fast multipole method (selected using fmm = )",
   GINT_TO_POINTER(BEM3D_SOLVER_FMM),

   "[exclusion distance for point sources in FMM]",
   "radius_fmm", "",
   NULL,

   "[tolerance for FMM approximations]",
   "tolerance_fmm", "",
   NULL,

   "[number of point sources per element in FMM]",   
   "skeleton_order", "",
   NULL,

   NULL, NULL, NULL, NULL} ;


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

static gint add_configuration_options(gpointer options[])

{
  gint i ;
  
  for ( i = 0 ; options[_CONFIG_STRIDE*i] != NULL ; i ++ ) {
    if ( options[_CONFIG_STRIDE*i + _CONFIG_POINTER] != NULL ) 
      bem3d_configuration_add_identifier(options[_CONFIG_STRIDE*i +
						 _CONFIG_NAME],
					 options[_CONFIG_STRIDE*i +
						 _CONFIG_POINTER]) ;
  }

  return 0 ;
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
  add_configuration_options(_config_quadratures) ;

  /*physics*/
  add_configuration_options(_config_physics) ;
  
  /*solver settings*/
  add_configuration_options(_config_solver) ;

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
 * @param description text string describing the variable (e.g. to be
 * used as a comment in a configuration file)
 * @param v a pointer to the identified object
 * 
 * @return ::BEM3D_SUCCESS on success
 */

gint bem3d_configuration_add_identifier(const gchar *id,
					gpointer v)

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
  if ( varmap == NULL ) return NULL ;

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

gint bem3d_configuration_write(BEM3DConfiguration *config,
			       gboolean write_descriptions,
			       gint line_width, FILE *f)

{
  
  return 0 ;
}

static gint write_configuration_options(gpointer options[],
					gboolean write_descriptions,
					gint line_width,
					gchar *prefix,
					gchar *comment,
					FILE *f)
{
  gint i ;

  if ( !write_descriptions ) {
    for ( i = 0 ; options[_CONFIG_STRIDE*i] != NULL ; i ++ ) {
      fprintf(f, "%s%s = %s\n",
	      comment,
	      (gchar *)(options[_CONFIG_STRIDE*i+_CONFIG_SETTING]),
	      (gchar *)(options[_CONFIG_STRIDE*i+_CONFIG_NAME])) ;
    }

    return 0 ;
  }

    for ( i = 0 ; options[_CONFIG_STRIDE*i] != NULL ; i ++ ) {
      printf_fixed_width((gchar *)(options[_CONFIG_STRIDE*i+_CONFIG_DESC]),
			 line_width, comment, f) ;
      fprintf(f, "%s%s = %s\n",
	      prefix,
	      (gchar *)(options[_CONFIG_STRIDE*i+_CONFIG_SETTING]),
	      (gchar *)(options[_CONFIG_STRIDE*i+_CONFIG_NAME])) ;
    }

  return 0 ;
}

/** 
 * Write generic configuration data containing all possible settings,
 * which can be used as a template for a configuration file. Output is
 * written in the form:
 *
 * [Group name]
 * \a prefix [setting] = [variable name]
 * \a prefix [setting] = [variable name]
 * \a prefix [setting] = [variable name]
 * 
 * @param write_descriptions insert comment with description of
 * variable before each line
 * @param line_width total number of characters per line in output, applied
 * to lines of description only
 * @param prefix string to insert at start of each variable line, such as `#'
 * to comment out entries
 * @param comment string to insert at start of descriptions, to comment out
 * text
 * @param f output file stream
 * 
 * @return 0 on success
 */

gint bem3d_configuration_write_generic(gboolean write_descriptions,
				       gint line_width, gchar *prefix,
				       gchar *comment,
				       FILE *f)

{
  fprintf(f, "\n[BEM3D::Quadrature]\n") ;
  write_configuration_options(_config_quadratures, write_descriptions,
			      line_width, prefix, comment, f) ;

  fprintf(f, "\n[BEM3D::Physics]\n") ;
  write_configuration_options(_config_physics, write_descriptions,
			      line_width, prefix, comment, f) ;

  fprintf(f, "\n[BEM3D::Solver]\n") ;
  write_configuration_options(_config_solver, write_descriptions,
			      line_width, prefix, comment, f) ;

  return 0 ;
}

/**
 * @}
 * 
 */
