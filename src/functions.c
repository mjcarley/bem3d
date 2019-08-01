/* functions.c
 * 
 * Copyright (C) 2006, 2010, 2018, 2019 Michael Carley
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
 * node normal, \a i for node index, and \a mesh for mesh index (when
 * multiple meshes are passed to the function evaluation). In integral
 * evaluations, \f$(x,y,z)\f$ refers to the node coordinates in an
 * integrand, \a X, \a Y, and \a Z are reserved for node coordinates
 * at an evaluation point, \a NX, \a NY, and \a NZ for the
 * corresponding normal, and \a I for the node index.
 *
 * Functions are applied to data at each node of a mesh
 * with input and output data in ::BEM3DMeshData structs. Input and
 * output data are specified in functions as \a f and \a g, with
 * components given as \a f[0], etc. Output data are given in \a f.
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

#define BEM3D_FUNCTION_DATA_WIDTH       16
#define BEM3D_FUNCTION_DATA_FUNC         0
#define BEM3D_FUNCTION_DATA_DATA_F       1
#define BEM3D_FUNCTION_DATA_VARS         2
#define BEM3D_FUNCTION_DATA_MOTION       3
#define BEM3D_FUNCTION_DATA_TIME         4
#define BEM3D_FUNCTION_DATA_OP           5
#define BEM3D_FUNCTION_DATA_DATA_G       6
#define BEM3D_FUNCTION_DATA_MESH_INDEX   7
#define BEM3D_FUNCTION_DATA_MESH         8
#define BEM3D_FUNCTION_DATA_VERTEX       9
#define BEM3D_FUNCTION_DATA_NORMAL      10
#define BEM3D_FUNCTION_DATA_INDEX       11
#define BEM3D_FUNCTION_DATA_VALUES      12
#define BEM3D_FUNCTION_DATA_QUADRATURE  13
#define BEM3D_FUNCTION_DATA_NVALS       14

#define BEM3D_FUNCTION_NRESERVED  18
#define BEM3D_FUNCTION_X           0
#define BEM3D_FUNCTION_Y           1
#define BEM3D_FUNCTION_Z           2
#define BEM3D_FUNCTION_U           3
#define BEM3D_FUNCTION_V           4
#define BEM3D_FUNCTION_W           5
#define BEM3D_FUNCTION_NX          6
#define BEM3D_FUNCTION_NY          7
#define BEM3D_FUNCTION_NZ          8
#define BEM3D_FUNCTION_INDEX       9
#define BEM3D_FUNCTION_MESH_INDEX 10
#define BEM3D_FUNCTION_X_EVAL     11
#define BEM3D_FUNCTION_Y_EVAL     12
#define BEM3D_FUNCTION_Z_EVAL     13
#define BEM3D_FUNCTION_NX_EVAL    14
#define BEM3D_FUNCTION_NY_EVAL    15
#define BEM3D_FUNCTION_NZ_EVAL    16
#define BEM3D_FUNCTION_INDEX_EVAL 17

gchar *BEM3D_FUNCTION_RESERVED[] = {"x", "y", "z", "u", "v", "w",
				    "nx", "ny", "nz", "i", "mesh",
				    "X", "Y", "Z", "NX", "NY", "NZ",
				    "I"} ;

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

  expanded = g_string_new(def) ;
  function_variables(expanded) ;

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
  gint i, nc, nex, npass = 0, npass_max = 64 ;

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
  } while ( nex != 0 && npass < npass_max ) ;

  if ( npass >= npass_max ) 
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

/** 
 * Insert a string, such as a variable definition, into a function
 * definition, overwriting any previous definition of the same
 * variable.
 * 
 * @param fn a ::BEM3DFunction to modify
 * @param str string containing definition to insert
 * 
 * @return 0 on success.
 */

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

#ifdef HAVE_LIBMATHEVAL
static void function_apply(gpointer key, gpointer val, gpointer data[])

{
  BEM3DFunction *func = data[BEM3D_FUNCTION_DATA_FUNC] ;
  BEM3DMeshData *f = data[BEM3D_FUNCTION_DATA_DATA_F] ;
  gchar **vars = data[BEM3D_FUNCTION_DATA_VARS] ;
  BEM3DMesh *mesh = data[BEM3D_FUNCTION_DATA_MESH] ;
  BEM3DMotion *motion = data[BEM3D_FUNCTION_DATA_MOTION] ;
  BEM3DOperator *op = data[BEM3D_FUNCTION_DATA_OP] ;
  BEM3DMeshData *g = data[BEM3D_FUNCTION_DATA_DATA_G] ;
  gint imesh = *(gint *)data[BEM3D_FUNCTION_DATA_MESH_INDEX] ;
  gint i, k, p, off, nvals ;
  gdouble values[128], result[128], *opi, *fi, *fp, *gi, t ;
  GtsVertex *v ;
  GtsVector u = {0.0}, n ;

  if ( mesh != NULL && motion != NULL ) {
    g_error("%s: mesh and motion should not both be non-NULL",
	    __FUNCTION__) ;
  }
  
  nvals = BEM3D_FUNCTION_NRESERVED + 4*bem3d_mesh_data_element_number(f)
    + ( g == NULL ? 0: 4*bem3d_mesh_data_element_number(g) ) ;

  g_assert(nvals < 128) ;
  i = GPOINTER_TO_INT(key)-1 ;
  fi = bem3d_mesh_data_get(f, i) ;
  memcpy(result, fi, bem3d_mesh_data_element_number(f)*sizeof(gdouble)) ;

  if ( mesh != NULL ) {
    if ( (v = bem3d_mesh_node_from_index(mesh, i)) == NULL )
      return ; /*this node is not in the mesh currently under consideration*/
    
  } else {
    g_assert_not_reached() ; /*untested code*/
    mesh = bem3d_motion_mesh(motion) ;
    if ( (v = bem3d_mesh_node_from_index(mesh, i)) == NULL ) return ;
    t = *((gdouble *)data[BEM3D_FUNCTION_DATA_TIME]) ;
    bem3d_motion_node_velocity(motion, i, t, u) ;    
  }

  bem3d_node_normal(mesh, i, n, BEM3D_AVERAGE_MWE) ;

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
  values[off + BEM3D_FUNCTION_MESH_INDEX] = (gdouble)imesh ;

  bem3d_operator_gradient(mesh, i, op, BEM3D_AVERAGE_MWA) ;

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
  
  off = BEM3D_FUNCTION_NRESERVED ;
  for ( k = 0 ; k < bem3d_function_function_number(func) ; k ++ ) {
    result[g_array_index(func->idx,gint,k)] = 
      evaluator_evaluate(g_ptr_array_index(func->evaluators,k),
			 nvals, vars, values) ;
  }

  memcpy(fi, result, bem3d_mesh_data_element_number(f)*sizeof(gdouble)) ;

  return ;
}

static gint bem3d_function_apply_index(BEM3DFunction *func, 
				       gpointer m0,
				       gdouble t,
				       BEM3DMeshData *f,
				       BEM3DMeshData *g,
				       gint imesh)
  
{
  gpointer data[BEM3D_FUNCTION_DATA_WIDTH] ;
  gchar **vars ;
  gint i, off ;
  BEM3DOperator *op ;

  if ( BEM3D_IS_MESH(m0) ) {
    data[BEM3D_FUNCTION_DATA_MESH] = m0 ;
    data[BEM3D_FUNCTION_DATA_MOTION] = NULL ;
    data[BEM3D_FUNCTION_DATA_TIME] = NULL ;
  } else {
    data[BEM3D_FUNCTION_DATA_MESH] = NULL ;
    data[BEM3D_FUNCTION_DATA_MOTION] = m0 ;
    data[BEM3D_FUNCTION_DATA_TIME] = &t ;
  }
  
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
  
  data[BEM3D_FUNCTION_DATA_FUNC] = func ;
  data[BEM3D_FUNCTION_DATA_DATA_F] = f ;
  data[BEM3D_FUNCTION_DATA_VARS] = vars ;
  data[BEM3D_FUNCTION_DATA_OP] = op = bem3d_operator_new() ;
  data[BEM3D_FUNCTION_DATA_DATA_G] = g ;
  data[BEM3D_FUNCTION_DATA_MESH_INDEX] = &imesh ;
  
  g_hash_table_foreach(f->t, (GHFunc)function_apply, data) ;

  return 0 ;
}
#endif /*HAVE_LIBMATHEVAL*/

/**
 * Apply a function to surface data specified as ::BEM3DMesh. The
 * function may use any of the reserved words which apply to a
 * ::BEM3DFunction and entries for the supplied data blocks \a f and
 * \a g (if not NULL) with reference to `f[0]', `f[1]', etc. for
 * elements of data and `dfdx[0]' etc. for gradient terms, which are
 * computed as required at each node. If a second data block is given,
 * it may be included in the function as `g[0]' etc. On output, \a f
 * will contain the results of the applied function.
 *
 * @param func ::BEM3DFunction to apply;
 * @param m ::BEM3DMesh for surface;
 * @param f a ::BEM3DMeshData block containing data for mesh (must 
 * not be NULL);
 * @param g a ::BEM3DMeshData block containing supplementary data 
 * (may be NULL);
 *
 * @return 0 on success.
 */

gint bem3d_function_apply_mesh(BEM3DFunction *func, 
			       BEM3DMesh *m,
			       BEM3DMeshData *f,
			       BEM3DMeshData *g)

{

  g_return_val_if_fail(func != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(func), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

#ifdef HAVE_LIBMATHEVAL
  bem3d_function_apply_index(func, m, 0, f, g, 0) ;
#else /*HAVE_LIBMATHEVAL*/
  g_warning("%s: function evaluation not implemented (requires libmatheval)") ;
#endif /*HAVE_LIBMATHEVAL*/
  return BEM3D_SUCCESS ;
}

/** 
 * Apply a function to a list of surface data, in the same way as
 * ::bem3d_function_apply_mesh. The list of surfaces is in the form
 * of an array of ::BEM3DMesh pointers, which are visited in turn in
 * the same manner as in ::bem3d_function_apply_mesh, with the
 * reserved variable \a mesh set to the index of the surface in the
 * list (so the ordering of the data in the list matters).
 * 
 * @param func ::BEM3DFunction to apply;
 * @param meshes array of pointers to ::BEM3DMesh for surfaces;
 * @param f a ::BEM3DMeshData block containing data for mesh (must 
 * not be NULL);
 * @param g a ::BEM3DMeshData block containing supplementary data 
 * (may be NULL);
 *
 * @return 0 on success.
 */

gint bem3d_function_apply_mesh_list(BEM3DFunction *func, 
				    GPtrArray *meshes,
				    BEM3DMeshData *f,
				    BEM3DMeshData *g)

{
  gint i ;
  BEM3DMesh *m ;
  
  g_return_val_if_fail(func != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(func), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

#ifdef HAVE_LIBMATHEVAL
  for ( i = 0 ; i < meshes->len ; i ++ ) {
    m = g_ptr_array_index(meshes, i) ;
    bem3d_function_apply_index(func, m, 0, f, g, i) ;
  }

#else /*HAVE_LIBMATHEVAL*/
  g_warning("%s: function evaluation not implemented (requires libmatheval)") ;
#endif /*HAVE_LIBMATHEVAL*/
  return BEM3D_SUCCESS ;
}

/**
 * Apply a function to surface data specified as ::BEM3DMotion,
 * computing surface velocity for passing to the function. The
 * function may use any of the reserved words which apply to a
 * ::BEM3DFunction and entries for the supplied data blocks \a f and
 * \a g (if not NULL) with reference to `f[0]', `f[1]', etc. for
 * elements of data and `dfdx[0]' etc. for gradient terms, which are
 * computed as required at each node. If a second data block is given,
 * it may be included in the function as `g[0]' etc. On output, \a f
 * will contain the results of the applied function.
 *
 * @param func ::BEM3DFunction to apply;
 * @param m ::BEM3DMotion for surface;
 * @param t time for evaluation of surface position and velocity using \a m;
 * @param f a ::BEM3DMeshData block containing data for mesh (must 
 * not be NULL);
 * @param g a ::BEM3DMeshData block containing supplementary data 
 * (may be NULL);
 *
 * @return 0 on success.
 */

gint bem3d_function_apply_motion(BEM3DFunction *func, 
				 BEM3DMotion *m,
				 gdouble t,
				 BEM3DMeshData *f,
				 BEM3DMeshData *g)

{
  g_return_val_if_fail(func != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(func), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

#ifdef HAVE_LIBMATHEVAL
  bem3d_function_apply_index(func, m, t, f, g, 0) ;
#else /*HAVE_LIBMATHEVAL*/
  g_warning("%s: function evaluation not implemented (requires libmatheval)") ;
#endif /*HAVE_LIBMATHEVAL*/
  return BEM3D_SUCCESS ;
}

/** 
 * Apply a function to a list of surface data, in the same way as
 * ::bem3d_function_apply_motion. The list of surfaces is in the form
 * of an array of ::BEM3DMotion pointers, which are visited in turn in
 * the same manner as in ::bem3d_function_apply_motion, with the
 * reserved variable \a mesh set to the index of the surface in the
 * list (so the ordering of the data in the list matters).
 * 
 * @param func ::BEM3DFunction to apply;
 * @param motions array of pointers to ::BEM3DMotion for surfaces;
 * @param t time for evaluation of surface position and velocity using \a m;
 * @param f a ::BEM3DMeshData block containing data for mesh (must 
 * not be NULL);
 * @param g a ::BEM3DMeshData block containing supplementary data 
 * (may be NULL);
 *
 * @return 0 on success.
 */

gint bem3d_function_apply_motion_list(BEM3DFunction *func, 
				      GPtrArray *motions,
				      gdouble t,
				      BEM3DMeshData *f,
				      BEM3DMeshData *g)
  
{
  gint i ;
  BEM3DMotion *m ;

  g_return_val_if_fail(func != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(motions != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_FUNCTION(func), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  for ( i = 0 ; i < motions->len ; i ++ ) {
    m = g_ptr_array_index(motions, i) ;
    bem3d_function_apply_index(func, m, t, f, g, i) ;
  }
  
  return 0 ;
}

static void integrate_func(BEM3DElement *e, gpointer *data)

{
  BEM3DFunction *func = data[BEM3D_FUNCTION_DATA_FUNC] ;
  /* BEM3DMesh *m = data[BEM3D_FUNCTION_DATA_MESH] ; */
  /* GtsVertex *x = data[BEM3D_FUNCTION_DATA_VERTEX] ; */
  /* gdouble *N = data[BEM3D_FUNCTION_DATA_NORMAL] ; */
  /* gint I = *((gint *)data[BEM3D_FUNCTION_DATA_INDEX]) ; */
  BEM3DMeshData *f = data[BEM3D_FUNCTION_DATA_DATA_F] ;
  /* gint imesh = *((gint *)data[BEM3D_FUNCTION_DATA_MESH_INDEX]) ; */
  gdouble *values = data[BEM3D_FUNCTION_DATA_VALUES] ;
  BEM3DQuadratureRule *q = data[BEM3D_FUNCTION_DATA_QUADRATURE] ;
  gchar **vars = data[BEM3D_FUNCTION_DATA_VARS] ;
  gint nvals = *((gint *)data[BEM3D_FUNCTION_DATA_NVALS]) ;
  
  gdouble L[32], dLds[32], dLdt[32], s, t, wt, n[3], J, result[128], *fj ;
  BEM3DShapeFunc shfunc = bem3d_element_shape_func(e) ;
  BEM3DShapeFunc cpfunc = bem3d_element_node_func(e) ;
  GtsPoint y ;
  gint i, j, k, off, ngp = 7 ;

  /*note no gradient operator included (all zeroed)*/
  
  bem3d_quadrature_rule_wx(NULL, e, q, NULL, NULL, &ngp, NULL) ;

  for ( i = 0 ; i < bem3d_quadrature_vertex_number(q) ; i ++ ) {
    s = bem3d_quadrature_xi(q,i) ;
    t = bem3d_quadrature_eta(q,i) ;
    wt = bem3d_quadrature_weight(q,i) ;
    shfunc(s, t, L, dLds, dLdt, NULL) ;
    bem3d_element_position(e, L, &y) ;
    bem3d_element_normal(e, dLds, dLdt, n, &J) ;

    cpfunc(s, t, L, NULL, NULL, NULL) ;

    off = 0 ;
    values[off + BEM3D_FUNCTION_X] = GTS_POINT(&y)->x ;
    values[off + BEM3D_FUNCTION_Y] = GTS_POINT(&y)->y ;
    values[off + BEM3D_FUNCTION_Z] = GTS_POINT(&y)->z ;
    values[off + BEM3D_FUNCTION_NX] = n[0] ;
    values[off + BEM3D_FUNCTION_NY] = n[1] ;
    values[off + BEM3D_FUNCTION_NZ] = n[2] ;
    
    /* off = BEM3D_FUNCTION_NRESERVED ; */

    /* for ( k = 0 ; k < 4*bem3d_mesh_data_element_number(f) ; k ++ ) */
    /*   values[off+k] = 0.0 ; */

    /* for ( k = 0 ; k < bem3d_mesh_data_element_number(f) ; k ++ )  */
    /*   values[off + 4*k+0] = fi[k] ; */

    /* off = BEM3D_FUNCTION_NRESERVED ; */

    /*evaluation of integrand at quadrature point*/
    for ( k = 0 ; k < bem3d_function_function_number(func) ; k ++ ) {
      result[g_array_index(func->idx,gint,k)] = 
	evaluator_evaluate(g_ptr_array_index(func->evaluators,k),
			 nvals, vars, values) ;
    }

    J *= wt ;
    for ( j = 0 ; j < bem3d_element_node_number(e) ; j ++ ) {
      k = bem3d_element_global_index(e, j) ;
      fj = bem3d_mesh_data_get(f, k) ;
      for ( k = 0 ; k < bem3d_function_function_number(func) ; k ++ ) {
	fj[k] += result[k]*L[j]*J ;
      }
    }
  }
  
  return ;
}

/** 
 * Fill a ::BEM3DMeshData with the integral of a function such that
 * the entries form the weights of a quadrature,
 * i.e. \f$\int_{S}w(\mathbf{x})\phi(\mathbf{x})\approx\sum_{i}\phi_{i}f_{i}\f$
 * where \f$\phi_{i}\f$ are the values of a function at the mesh nodes
 * and \f$f_{i}\f$ are the entries in the ::BEM3DMeshData
 * 
 * @param func weighting function to integrate (may have multiple terms);
 * @param m mesh to integrate on;
 * @param imesh index of \a m, used to set reserved word \a mesh in function
 * evaluation;
 * @param x used to set coordinates of evaluation point in position-dependent
 * integrands;
 * @param n normal at evaluation point;
 * @param i index of evaluation point;
 * @param q ::BEM3DQuadratureRule used for integration on elements;
 * @param f output data containing integral weights.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_function_integral_weights(BEM3DFunction *func,
				     BEM3DMesh *m, gint imesh,
				     GtsVertex *x, GtsVector n, gint i,
				     BEM3DQuadratureRule *q,
				     BEM3DMeshData *f)

{
  gpointer data[BEM3D_FUNCTION_DATA_WIDTH] ;
  gint nvals, off ;
  gdouble values[128] = {0.0} ;
  gchar **vars ;

  nvals = BEM3D_FUNCTION_NRESERVED + 4*bem3d_mesh_data_element_number(f) ;
  g_assert(nvals < 128) ;

  vars = (gchar **)g_malloc(nvals*sizeof(gchar *)) ;
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

  data[BEM3D_FUNCTION_DATA_FUNC] = func ;
  data[BEM3D_FUNCTION_DATA_MESH] = m ;
  data[BEM3D_FUNCTION_DATA_VERTEX] = x ;
  data[BEM3D_FUNCTION_DATA_NORMAL] = n ;
  data[BEM3D_FUNCTION_DATA_INDEX] = &i ;
  data[BEM3D_FUNCTION_DATA_DATA_F] = f ;
  data[BEM3D_FUNCTION_DATA_MESH_INDEX] = &imesh ;
  data[BEM3D_FUNCTION_DATA_VALUES] = values ;
  data[BEM3D_FUNCTION_DATA_QUADRATURE] = q ;
  data[BEM3D_FUNCTION_DATA_NVALS] = &nvals ;
  data[BEM3D_FUNCTION_DATA_VARS] = vars ;


  off = 0 ;
  for ( i = 0 ; i < BEM3D_FUNCTION_NRESERVED ; i ++ )     
    vars[off+i] = BEM3D_FUNCTION_RESERVED[i] ;

  if ( x != NULL ) {
    values[off + BEM3D_FUNCTION_X_EVAL] = GTS_POINT(x)->x ;
    values[off + BEM3D_FUNCTION_Y_EVAL] = GTS_POINT(x)->y ;
    values[off + BEM3D_FUNCTION_Z_EVAL] = GTS_POINT(x)->z ;
  }
  if ( n != NULL ) {
    values[off + BEM3D_FUNCTION_NX_EVAL] = n[0] ;
    values[off + BEM3D_FUNCTION_NY_EVAL] = n[1] ;
    values[off + BEM3D_FUNCTION_NZ_EVAL] = n[2] ;
  }
  values[off + BEM3D_FUNCTION_INDEX_EVAL] = (gdouble)i ;
  values[off + BEM3D_FUNCTION_MESH_INDEX] = (gdouble)imesh ;

  for ( i = 0 ; i < bem3d_function_function_number(func) ; i ++ ) {
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
  
  bem3d_mesh_foreach_element(m, (BEM3DElementFunc)integrate_func, data) ;
  
  return 0 ;
}

/** 
 * Evaluate a ::BEM3DFunction at a ::GtsPoint. The function is
 * evaluated in the same way as for a ::BEM3DMesh except that the
 * normal and index must be supplied explicitly.
 * 
 * @param func a ::BEM3DFunction;
 * @param x a ::GtsPoint;
 * @param n the normal at \a x;
 * @param idx the mesh index to be used for \a x in \a f, if required;
 * @param result the results of the evaluation of \a f;
 * @param nres number of entries in \a result;
 * 
 * @return the number of functions evaluated, on success, -1 otherwise.
 */

gint bem3d_function_eval_point(BEM3DFunction *func,
			       GtsPoint *x, GtsVector n, gint idx,
			       gdouble *result, gint nres)

{
  gint nvals, off, i, k ;
  gdouble values[128] = {0.0} ;
  gchar **vars ;

  nvals = BEM3D_FUNCTION_NRESERVED + 4*nres ;
  g_assert(nvals < 128) ;

  vars = (gchar **)g_malloc(nvals*sizeof(gchar *)) ;
  /*reserved names (x, y, z, etc)*/
  off = 0 ;
  for ( i = 0 ; i < BEM3D_FUNCTION_NRESERVED ; i ++ )     
    vars[off+i] = BEM3D_FUNCTION_RESERVED[i] ;

  /*variable names packed into vars as f[0], f[1], ... 
    dfdx[0], dfdy[0], dfdz[0], dfdx[1], ...*/
  off = BEM3D_FUNCTION_NRESERVED ;
  for ( i = 0 ; i < nres ; i ++ ) {
    vars[off + 4*i+0] = g_strdup_printf("_F%d", i) ;
    vars[off + 4*i+1] = g_strdup_printf("_DFDX%d", i) ;
    vars[off + 4*i+2] = g_strdup_printf("_DFDY%d", i) ;
    vars[off + 4*i+3] = g_strdup_printf("_DFDZ%d", i) ;
  }

  off = 0 ;
  for ( i = 0 ; i < BEM3D_FUNCTION_NRESERVED ; i ++ )     
    vars[off+i] = BEM3D_FUNCTION_RESERVED[i] ;

  
  off = 0 ;
  if ( x != NULL ) {
    values[off + BEM3D_FUNCTION_X] = GTS_POINT(x)->x ;
    values[off + BEM3D_FUNCTION_Y] = GTS_POINT(x)->y ;
    values[off + BEM3D_FUNCTION_Z] = GTS_POINT(x)->z ;
  }
  if ( n != NULL ) {
    values[off + BEM3D_FUNCTION_NX] = n[0] ;
    values[off + BEM3D_FUNCTION_NY] = n[1] ;
    values[off + BEM3D_FUNCTION_NZ] = n[2] ;
  }
  values[off + BEM3D_FUNCTION_INDEX] = (gdouble)idx ;

  for ( i = 0 ; i < bem3d_function_function_number(func) ; i ++ ) {
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

  i = -1 ;
  for ( k = 0 ; k < bem3d_function_function_number(func) ; k ++ ) {
    result[g_array_index(func->idx,gint,k)] = 
      evaluator_evaluate(g_ptr_array_index(func->evaluators,k),
			 nvals, vars, values) ;
    i = MAX(i, g_array_index(func->idx,gint,k)) ;
  }
  
  return i ;
}

/**
 * @}
 * 
 */
