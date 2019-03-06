/* functions.c
 * 
 * Copyright (C) 2006, 2010, 2018 Michael Carley
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
 * @defgroup functions Functions of surface data
 *
 * Application of analytical functions to data on meshes. The
 * underlying function handler is GNU libmatheval and functions should
 * conform to that syntax, including the built-in functions. Reserved
 * variable names are \a x, \a y, \a z for node coordinates, \a u, \a
 * v, \a w for components of node velocity, \a nx, \a ny, \a nz for
 * node normal, and \a i for node index. Functions are applied to data
 * at each node of a mesh with input and output data in
 * ::BEM3DMeshData structs. Input and output data are specified in
 * functions as \a f and \a g, with components given as \a f[0],
 * etc. Output data are given in \a f.
 *
 * @{
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <glib.h>
#include <gts.h>

#include <wmpi.h>

#ifdef HAVE_LIBMATHEVAL
#include <matheval.h>
#endif /*HAVE_LIBMATHEVAL*/

#include "bem3d.h"
#include "bem3d-private.h"

#define BEM3D_EFDATA_WIDTH  8
#define BEM3D_EFDATA_DATA   0
#define BEM3D_EFDATA_FUNC   1
#define BEM3D_EFDATA_MDATA  2

#define BEM3D_FUNCTION_NRESERVED 10
#define BEM3D_FUNCTION_X         0
#define BEM3D_FUNCTION_Y         1
#define BEM3D_FUNCTION_Z         2
#define BEM3D_FUNCTION_U         3
#define BEM3D_FUNCTION_V         4
#define BEM3D_FUNCTION_W         5
#define BEM3D_FUNCTION_NX        6
#define BEM3D_FUNCTION_NY        7
#define BEM3D_FUNCTION_NZ        8
#define BEM3D_FUNCTION_INDEX     9

gchar *BEM3D_FUNCTION_RESERVED[] = {"x", "y", "z", "u", "v", "w",
				    "nx", "ny", "nz", "i"} ;

static void bem3d_function_class_init (BEM3DFunctionClass * klass)
{
  /* define new methods and overload inherited methods here */

}

static void bem3d_function_init (BEM3DFunction * object)
{
  /* initialize object here */
}

BEM3DFunctionClass * bem3d_function_class (void)
{
  static BEM3DFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo b_e_m3_d_function_info = {
      "BEM3DFunction",
      sizeof (BEM3DFunction),
      sizeof (BEM3DFunctionClass),
      (GtsObjectClassInitFunc) bem3d_function_class_init,
      (GtsObjectInitFunc) bem3d_function_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &b_e_m3_d_function_info);
  }

  return klass;
}

BEM3DFunction * bem3d_function_new (BEM3DFunctionClass * klass)

{
  BEM3DFunction * object;

  object = BEM3D_FUNCTION (gts_object_new (GTS_OBJECT_CLASS (klass)));

  object->defs = g_hash_table_new(g_str_hash, g_str_equal) ;
  object->idx = g_array_new(TRUE, TRUE, sizeof(gint)) ;
  object->definitions = g_ptr_array_new() ;
  object->functions = g_ptr_array_new() ;
  object->expansions = g_ptr_array_new() ;
  object->evaluators = g_ptr_array_new() ;

  return object;
}

static void expand_variable_name(GString *s, gint i)

{
  gint i0, i1, len, j ;
  gchar substr[64] ;

  g_debug("%s: %s\n", __FUNCTION__, s->str) ;
  substr[0] = '_' ; len = 1 ;

  /* fprintf(stderr, "s: %s -> ", s->str) ; */
  
  /*find the start of the variable name*/
  for ( i0 = i-1 ; i0 >= 0 && g_ascii_isalpha(s->str[i0]) ; i0 -- ) ;
  i0 ++ ;

  if ( i - i0 > 32 ) 
    g_error("%s: syntax error in `%s' near `%s'", 
	    __FUNCTION__, s->str, &(s->str[i0])) ;

  len += i-i0 ;
  strncpy(&(substr[1]), &(s->str[i0]), i-i0) ; substr[len] = '\0' ;
  
  if ( strcmp(&(substr[1]), "f") &&
       strcmp(&(substr[1]), "dfdx") &&
       strcmp(&(substr[1]), "dfdy") &&
       strcmp(&(substr[1]), "dfdz") &&
       strcmp(&(substr[1]), "g") &&
       strcmp(&(substr[1]), "dgdx") &&
       strcmp(&(substr[1]), "dgdy") &&
       strcmp(&(substr[1]), "dgdz")
       )
    g_error("%s: unrecognized token `%s' in expression `%s'",
	    __FUNCTION__, substr, s->str) ;
  
  for ( i1 = i+1 ; i1 < s->len && g_ascii_isdigit(s->str[i1]) ; i1 ++ ) ;
  if ( s->str[i1] != ']' || i1-i > 32 )
    g_error("%s: syntax error in `%s' near `%s'", 
	    __FUNCTION__, s->str, &(s->str[i1])) ;
  i1 -- ;

  strncpy(&(substr[len]), &(s->str[i+1]), i1-i) ; len += i1-i ;
  substr[len] = '\0' ;

  for ( j = 0 ; j < len ; j ++ ) 
    s->str[i0+j] = g_ascii_toupper(substr[j]) ;
  g_string_erase(s, i1+1, 1) ;

  g_debug("%s: %s\n", __FUNCTION__, s->str) ;

  /* fprintf(stderr, "%s\n", s->str) ; */

  return ;
}

static void function_variables(GString *f)

{
  gint i ;

  for ( i = 0 ; i < f->len ; i ++ ) {
    if ( f->str[i] == '[' ) expand_variable_name(f, i) ;
  }

  return ;
}

gint bem3d_function_add_function(BEM3DFunction *f, gint i, gchar *def)

{
  GString *expanded ;

  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(f), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(i >= 0, BEM3D_ARGUMENT_OUT_OF_RANGE) ;
  g_return_val_if_fail(def != NULL, BEM3D_NULL_ARGUMENT) ;

  /* fprintf(stderr, "add %s -> ", def) ; */
  
  expanded = g_string_new(def) ;
  function_variables(expanded) ;

  /* fprintf(stderr, " %s\n", expanded->str) ; */

  g_ptr_array_add(f->definitions, g_strdup(def)) ;
  g_ptr_array_add(f->functions, expanded) ;
  g_ptr_array_add(f->evaluators, NULL) ;
  g_array_append_val(f->idx, i) ;

  return BEM3D_SUCCESS ;
}

gboolean bem3d_function_token_is_reserved(gchar *token)

{
  gint i ;
  
  g_return_val_if_fail(token != NULL, FALSE) ;  

  for ( i = 0 ; i < BEM3D_FUNCTION_NRESERVED ; i ++ ) 
    if ( strcmp(token, BEM3D_FUNCTION_RESERVED[i]) == 0 )
      return TRUE ;

  return FALSE ;
}

/**
 * Add a variable definition to a ::BEM3DFunction. This will overwrite
 * any existing definition of the same variable.
 *
 * @param f a ::BEM3DFunction;
 * @param var variable name;
 * @param def variable definition.
 *
 * @return 0 on success.
 */

gint bem3d_function_variable_add(BEM3DFunction *f, gchar *var, gchar *def)

{
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(f), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(var != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(def != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( bem3d_function_token_is_reserved(var) || 
       (var[0] == 'f' && !g_ascii_isalpha(var[1])) ||
       (var[0] == 'g' && !g_ascii_isalpha(var[1]))
       )
    g_error("%s: `%s' is a reserved token and may not be redefined",
	    __FUNCTION__, var) ;

  if ( g_hash_table_lookup(f->defs, var) != NULL ) 
    g_hash_table_remove(f->defs, var) ;

  g_hash_table_insert(f->defs, var, def) ;

  return BEM3D_SUCCESS ;
}

gchar *bem3d_function_variable_lookup(BEM3DFunction *f, gchar *var)

{
  g_return_val_if_fail(f != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(f), NULL) ;
  g_return_val_if_fail(var != NULL, NULL) ;

  return (g_hash_table_lookup(f->defs, var)) ;
}


static void insert_definition(GString *s, gint off, gint len, gchar *def)

{
  g_string_erase(s, off, len) ;
  g_string_insert_c(s, off, '(') ;
  g_string_insert_c(s, off+1, ')') ;
  g_string_insert(s, off+1, def) ;

  return ;
}

static void copy_and_expand_function(BEM3DFunction *f, GString *e, GString *s,
				     GString *buf) 

{
  gchar *def ;
  gint i, nc, nex, npass = 0 ;

  g_string_assign(e, s->str) ;

  do {
    nex = 0 ;
    for ( i = 0 ; i < e->len ; i ++ ) {
      g_string_assign(buf, "") ;
      nc = strcspn(&(e->str[i]), "+-*/^() ") ;
      g_string_append_len(buf, &(e->str[i]), nc) ;
      if ( (nc != 0) && 
	   (def = bem3d_function_variable_lookup(f, buf->str)) != NULL ) {
	insert_definition(e, i, nc, def) ;
	nex ++ ;
      } else i += nc ;
    }
    npass ++ ;
  } while ( nex != 0 && npass < 64 ) ;

  if ( npass > 63 ) 
    g_error("%s: variables seem to have a circular definition",
	    __FUNCTION__) ;

  return ;
}

gint bem3d_function_expand_functions(BEM3DFunction *f)

{
  gint i ;
  GString *buf ;

  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(f), BEM3D_ARGUMENT_WRONG_TYPE) ;

  buf = g_string_new("") ;

  for ( i = 0 ; i < f->functions->len ; i ++ ) {
    if ( f->expansions->len <= i )
      g_ptr_array_add(f->expansions, g_string_new("")) ;
    copy_and_expand_function(f, 
			     g_ptr_array_index(f->expansions, i),
			     g_ptr_array_index(f->functions, i), buf) ;
			     
  }

  g_string_free(buf, TRUE) ;

  return BEM3D_SUCCESS ;
}

#ifdef HAVE_LIBMATHEVAL
static void function_apply(gpointer key, gpointer val, gpointer data[])

{
  BEM3DFunction *func = data[0] ;
  BEM3DMeshData *f = data[1] ;
  gchar **vars = data[2] ;
  BEM3DMotion *m = data[3] ;
  gdouble t = *((gdouble *)data[4]) ;
  BEM3DOperator *op = data[5] ;
  BEM3DMeshData *g = data[6] ;
  gint i, k, p, off, nvals ;
  gdouble values[128], result[128], *opi, *fi, *fp, *gi ;
  GtsVertex *v ;
  GtsVector u, n ;

  nvals = BEM3D_FUNCTION_NRESERVED + 4*bem3d_mesh_data_element_number(f)
    + ( g == NULL ? 0: 4*bem3d_mesh_data_element_number(g) ) ;

  g_assert(nvals < 128) ;
  
  i = GPOINTER_TO_INT(key)-1 ;
  fi = bem3d_mesh_data_get(f, i) ;
  memcpy(result, fi, bem3d_mesh_data_element_number(f)*sizeof(gdouble)) ;
	 
  v = bem3d_mesh_node_from_index(bem3d_motion_mesh(m), i) ;
  bem3d_node_normal(bem3d_motion_mesh(m), i, n, BEM3D_AVERAGE_MWE) ;
  bem3d_motion_node_velocity(m, i, t, u) ;

  off = 0 ; 
  values[off + BEM3D_FUNCTION_X] = GTS_POINT(v)->x ;
  values[off + BEM3D_FUNCTION_Y] = GTS_POINT(v)->y ;
  values[off + BEM3D_FUNCTION_Z] = GTS_POINT(v)->z ;
  values[off + BEM3D_FUNCTION_U] = u[0] ;
  values[off + BEM3D_FUNCTION_V] = u[1] ;
  values[off + BEM3D_FUNCTION_W] = u[2] ;
  values[off + BEM3D_FUNCTION_NX] = n[0] ;
  values[off + BEM3D_FUNCTION_NY] = n[1] ;
  values[off + BEM3D_FUNCTION_NZ] = n[2] ;
  values[off + BEM3D_FUNCTION_INDEX] = (gdouble)i ;

  bem3d_operator_gradient(bem3d_motion_mesh(m), i, op, BEM3D_AVERAGE_MWA) ;

  off = BEM3D_FUNCTION_NRESERVED ;

  for ( k = 0 ; k < 4*bem3d_mesh_data_element_number(f) ; k ++ )
    values[off+k] = 0.0 ;

  for ( k = 0 ; k < bem3d_mesh_data_element_number(f) ; k ++ ) 
    values[off + 4*k+0] = fi[k] ;

  for ( p = 0 ; p < bem3d_operator_length(op) ; p ++ ) {
    opi = bem3d_operator_weight(op,p) ;
    g_assert(!isnan(opi[0])) ;
    fp = bem3d_mesh_data_get(f, bem3d_operator_index(op,p)) ;
    for ( k = 0 ; k < bem3d_mesh_data_element_number(f) ; k ++ ) {
      values[off + 4*k+1] += fp[k]*opi[0] ;
      values[off + 4*k+2] += fp[k]*opi[1] ;
      values[off + 4*k+3] += fp[k]*opi[2] ;
    }
  }

  if ( g != NULL ) {
    off = 4*bem3d_mesh_data_element_number(f) + BEM3D_FUNCTION_NRESERVED ;  

    gi = bem3d_mesh_data_get(g, i) ;
    
    for ( k = 0 ; k < 4*bem3d_mesh_data_element_number(g) ; k ++ )
      values[off+k] = 0.0 ;

    for ( k = 0 ; k < bem3d_mesh_data_element_number(g) ; k ++ ) 
      values[off + 4*k+0] = gi[k] ;

    for ( p = 0 ; p < bem3d_operator_length(op) ; p ++ ) {
      opi = bem3d_operator_weight(op,p) ;
      g_assert(!isnan(opi[0])) ;
      fp = bem3d_mesh_data_get(g, bem3d_operator_index(op,p)) ;
      for ( k = 0 ; k < bem3d_mesh_data_element_number(g) ; k ++ ) {
	values[off + 4*k+1] += fp[k]*opi[0] ;
	values[off + 4*k+2] += fp[k]*opi[1] ;
	values[off + 4*k+3] += fp[k]*opi[2] ;
      }
    }
  }

  /* for ( i = 0 ; i < nvals ; i ++ ) { */
  /*   fprintf(stderr, "%s %lg\n", vars[i], values[i]) ; */
  /* } */

  /* exit(0) ; */
  
  off = BEM3D_FUNCTION_NRESERVED ;
  for ( k = 0 ; k < bem3d_function_function_number(func) ; k ++ ) {
    result[g_array_index(func->idx,gint,k)] = 
      evaluator_evaluate(g_ptr_array_index(func->evaluators,k),
			 nvals, vars, values) ;
  }

  memcpy(fi, result, bem3d_mesh_data_element_number(f)*sizeof(gdouble)) ;

  return ;
}
#endif /*HAVE_LIBMATHEVAL*/

/**
 * Apply a function to surface data. The function may use any of the
 * reserved words which apply to a ::BEM3DFunction and entries for the
 * supplied data blocks \a f and \a g (if not NULL) with reference to
 * `f[0]', `f[1]', etc. for elements of data and `dfdx[0]' etc. for
 * gradient terms, which are computed as required at each node. If a
 * second data block is given, it may be included in the function as
 * `g[0]' etc. On output, \a f will contain the results of the applied
 * function.
 *
 * @param func ::BEM3DFunction to apply;
 * @param m ::BEM3DMotion for surface;
 * @param t time for evaluation of surface position using \a m;
 * @param f a ::BEM3DMeshData block containing data for mesh (must 
 * not be NULL);
 * @param g a ::BEM3DMeshData block containing supplementary data 
 * (may be NULL);
 *
 * @return 0 on success.
 */

gint bem3d_function_apply(BEM3DFunction *func, 
			  BEM3DMotion *m,
			  gdouble t,
			  BEM3DMeshData *f,
			  BEM3DMeshData *g)

{
  gpointer data[8] ;
  gchar **vars ;
  gint i, off ;
  BEM3DOperator *op ;

  g_return_val_if_fail(func != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(func), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

#ifdef HAVE_LIBMATHEVAL
  if ( g == NULL ) 
    vars = (gchar **)g_malloc((4*bem3d_mesh_data_element_number(f)+
			       BEM3D_FUNCTION_NRESERVED)*
			      sizeof(gchar *)) ;
  else
    vars = (gchar **)g_malloc((4*bem3d_mesh_data_element_number(f)+
			       4*bem3d_mesh_data_element_number(g)+
			       BEM3D_FUNCTION_NRESERVED)*
			      sizeof(gchar *)) ;
   
  for ( i = 0 ; i < bem3d_function_function_number(func) ; i ++ ) {
    if ( g_array_index(func->idx,gint,i) >= 
	 bem3d_mesh_data_element_number(f) )
      g_error("%s: there is no f[%d] in %d element mesh data block",
	      __FUNCTION__, g_array_index(func->idx,gint,i),
	      bem3d_mesh_data_element_number(f)) ;
  }

  /*reserved names (x, y, z, etc)*/
  off = 0 ;
  for ( i = 0 ; i < BEM3D_FUNCTION_NRESERVED ; i ++ )     
    vars[off+i] = BEM3D_FUNCTION_RESERVED[i] ;

  /*variable names packed into vars as f[0], f[1], ... 
    dfdx[0], dfdy[0], dfdz[0], dfdx[1], ...*/
  off = BEM3D_FUNCTION_NRESERVED ;
  for ( i = 0 ; i < bem3d_mesh_data_element_number(f) ; i ++ ) {
    vars[off + 4*i+0] = g_strdup_printf("_F%d", i) ;
    vars[off + 4*i+1] = g_strdup_printf("_DFDX%d", i) ;
    vars[off + 4*i+2] = g_strdup_printf("_DFDY%d", i) ;
    vars[off + 4*i+3] = g_strdup_printf("_DFDZ%d", i) ;
  }

  /*if a second set of data is supplied, include it in the vars*/
  if ( g != NULL ) {
    off = 4*bem3d_mesh_data_element_number(f) + BEM3D_FUNCTION_NRESERVED ;  
    for ( i = 0 ; i < bem3d_mesh_data_element_number(g) ; i ++ ) {
      vars[off + 4*i+0] = g_strdup_printf("_G%d", i) ;
      vars[off + 4*i+1] = g_strdup_printf("_DGDX%d", i) ;
      vars[off + 4*i+2] = g_strdup_printf("_DGDY%d", i) ;
      vars[off + 4*i+3] = g_strdup_printf("_DGDZ%d", i) ;
    }
  }

  if ( func->expansions->len != func->evaluators->len ) 
    g_error("%s: mismatch between number of evaluators (%d) and "
	    "expansions (%d)",
	    __FUNCTION__, func->evaluators->len, func->expansions->len) ;

  for ( i = 0 ; i < bem3d_function_function_number(func) ; i ++ ) {
    /* fprintf(stderr, "%s\n", */
    /* 	    ((GString *)(g_ptr_array_index(func->expansions,i)))->str) ; */
    if ( g_ptr_array_index(func->evaluators, i) != NULL ) 
      evaluator_destroy(g_ptr_array_index(func->evaluators, i)) ;
    if ( (g_ptr_array_index(func->evaluators, i)  = 
	  evaluator_create(((GString *)
			    (g_ptr_array_index(func->expansions,i)))->str))
	 == NULL ) 
      g_error("%s: evaluator_create failed for function (%s)",
	      __FUNCTION__, 
	      ((GString *)(g_ptr_array_index(func->expansions,i)))->str) ;
  }

  data[0] = func ; data[1] = f ; data[2] = vars ; data[3] = m ;
  data[4] = &t ;
  data[5] = op = bem3d_operator_new() ;
  data[6] = g ;

  g_hash_table_foreach(f->t, (GHFunc)function_apply, data) ;
#else /*HAVE_LIBMATHEVAL*/
  g_warning("%s: function evaluation not implemented (requires libmatheval)") ;
#endif /*HAVE_LIBMATHEVAL*/
  return BEM3D_SUCCESS ;
}

static void write_function_defs(gchar *var, gchar *def, FILE *f)

{
  fprintf(f, "%s = %s\n", var, def) ;

  return ;
}

gint bem3d_function_write(BEM3DFunction *f, FILE *fid)

{
  gint i ;

  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(f), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(fid != NULL, BEM3D_NULL_ARGUMENT) ;
  
  fprintf(fid, "BEM3DFunction\n") ;
  g_hash_table_foreach(f->defs, (GHFunc)write_function_defs, fid) ;
  for ( i = 0 ; i < bem3d_function_function_number(f) ; i ++ ) {
    fprintf(fid, "f[%d] = %s\n", g_array_index(f->idx,gint,i),
	    (gchar *)(g_ptr_array_index(f->definitions, i))) ;
  }

  return BEM3D_SUCCESS ;
}

static gint parse_and_insert(BEM3DFunction *f, gchar *var, gchar *def)

{
  gint i, j ;

  g_strstrip(var) ;

  /* fprintf(stderr, "%s = %s;\n", var, def) ; */
  
  for ( i = 0 ; (var[i] != '\0') && (var[i] != '[') ; i ++ ) ;
  for ( j = 0 ; (var[j] != '\0') && (var[j] != ']') ; j ++ ) ;

  if ( (var[i] == '\0') && (var[j] == '\0') ) {
    /*no square brackets*/
    bem3d_function_variable_add(f, var, def) ;
    return 0 ;
  }

  /*check for only one bracket or brackets in wrong order or nothing 
    between them*/
  if ( (var[i] == '\0') || (var[j] == '\0') ) return 1 ;
  if ( j < i || i == j-1 ) return 1 ;

  /*invalid name*/
  if ( var[0] != 'f' || g_ascii_isalnum(var[1]) ) return 2 ;

  var[j] = '\0' ;
  i = atoi(&(var[i+1])) ;
  bem3d_function_add_function(f, i, def) ;

  return 0 ;
}

gint bem3d_function_read(BEM3DFunction *fn, GtsFile *f)

{
  gchar *var, *def ;
  GString *dstr ;
  gint i ;

  g_return_val_if_fail(fn != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(fn), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;
  
  dstr = g_string_new("") ;

  if ( f->type != GTS_STRING ) {
    gts_file_error(f, "expecting a string (BEM3DFunction)") ;
    return f->line ;
  }  

  gts_file_first_token_after(f, '\n') ;
  while ( f->type != GTS_NONE ) {
    if ( f->type != GTS_STRING ) {
      gts_file_error(f, "expecting a string (variable name or output)") ;
      return f->line ;
    }
    var = g_strdup(f->token->str) ;
    gts_file_next_token(f) ;

    if ( strcmp(f->token->str, "=") != 0 ) {
      gts_file_error(f, "expecting a string (equal sign)") ;
      return f->line ;
    }
    
    g_string_assign(dstr, "") ;
    gts_file_next_token(f) ;
    while ( strcmp(f->token->str, "\n") != 0)  {
      g_string_append(dstr, f->token->str) ;
      gts_file_next_token(f) ;
    }

    for ( i = 0 ; i < dstr->len ; i ++ ) {
      if ( dstr->str[i] == '[' ) expand_variable_name(dstr, i) ;
    }
    def = g_strdup(dstr->str) ;
    if ( (i = parse_and_insert(fn, var, def)) != 0 )
      g_error("%s: parse error at line %u:\n\t%s = %s", __FUNCTION__, f->line,
	      var, def) ;
    gts_file_first_token_after(f, '\n') ;
  }  

  return BEM3D_SUCCESS ;
}

gint bem3d_function_insert_string(BEM3DFunction *fn, gchar *str)

{
  gchar **tokens ;
  gint i ;

  g_return_val_if_fail(fn != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(fn), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(str != NULL, BEM3D_NULL_ARGUMENT) ;

  tokens = g_strsplit(str, "=", 0) ;
  if ( tokens[0] == NULL )
    g_error("%s: cannot parse expression \"%s\" (no variable specified)",
	    __FUNCTION__, str) ; 

  if ( tokens[1] == NULL )
    g_error("%s: cannot parse expression \"%s\" (no right hand side)",
	    __FUNCTION__, str) ; 

  if ( tokens[2] != NULL )
    g_error("%s: cannot parse expression \"%s\" (extra =)",
	    __FUNCTION__, str) ; 
  
  if ( (i = parse_and_insert(fn, tokens[0], tokens[1])) != 0 )
    g_error("%s: syntax error, cannot parse \"%s\"", __FUNCTION__, str) ;

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
