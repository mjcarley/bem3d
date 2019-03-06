/* motion.c
 * 
 * Copyright (C) 2009, 2018 Michael Carley
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>
#include <string.h>

#include <bem3d.h>

#include "bem3d-private.h"

#ifdef HAVE_LIBMATHEVAL
#include <matheval.h>
#endif /*HAVE_LIBMATHEVAL*/

#define BEM3D_MOTION_NVALUES    12
#define BEM3D_MOTION_X0          0
#define BEM3D_MOTION_Y0          1
#define BEM3D_MOTION_Z0          2
#define BEM3D_MOTION_TIME        3
#define BEM3D_MOTION_INDEX       4
#define BEM3D_MOTION_SQRT_M1     5
#define BEM3D_MOTION_NX          6
#define BEM3D_MOTION_NY          7
#define BEM3D_MOTION_NZ          8
#define BEM3D_MOTION_NX0         9
#define BEM3D_MOTION_NY0        10
#define BEM3D_MOTION_NZ0        11

gchar *BEM3D_MOTION_VARIABLES[] = {"x0", "y0", "z0", "t", "i", "j",
				   "nx", "ny", "nz", "nx0", "ny0", "nz0"} ;

/* BEM3DMotion: Object */

static void bem3d_motion_class_init (BEM3DMotionClass * klass)
{
  /* define new methods and overload inherited methods here */

}

static void bem3d_motion_init (BEM3DMotion * object)
{
  /* initialize object here */
  object->m = object->m0 = NULL ;
  object->defs = g_hash_table_new(g_str_hash, g_str_equal) ;
  g_hash_table_insert(object->defs, "x", "x0") ;
  g_hash_table_insert(object->defs, "y", "y0") ;
  g_hash_table_insert(object->defs, "z", "z0") ;
  g_hash_table_insert(object->defs, "u", "0") ;
  g_hash_table_insert(object->defs, "v", "0") ;
  g_hash_table_insert(object->defs, "w", "0") ;

  object->x = g_string_new("") ;
  object->y = g_string_new("") ;
  object->z = g_string_new("") ;
  object->u = g_string_new("") ;
  object->v = g_string_new("") ;
  object->w = g_string_new("") ;

  object->fx = object->fy = object->fz = 
    object->fdx = object->fdy = object->fdz = 
    object->fd2x = object->fd2y = object->fd2z = 
    object->fvx = object->fvy = object->fvz = NULL ;
    
  return ;
}

/**
 * @defgroup motion Motion specification
 * 
 * BEM3D contains a simple interpreter, using the libmatheval library
 * (http://www.gnu.org/software/libmatheval/) of Aleksandar Samardzic,
 * which can be used to analytically specify the motion of the nodes
 * of a ::BEM3DMesh and analytically derive the surface velocity for
 * use in the boundary condition. Motion definitions are defined using
 * equations of the form `x = ...', where the right hand side may
 * include other user-defined functions and constants. See
 * ::bem3d_motion_variable_add to see how this works.
 *
 * @{
 */


/** 
 * The basic class for a ::BEM3DMotion.
 * 
 * 
 * @return the class definition for a ::BEM3DMotion.
 */

BEM3DMotionClass * bem3d_motion_class (void)

{
  static BEM3DMotionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo bem3d_motion_info = {
      "BEM3DMotion",
      sizeof (BEM3DMotion),
      sizeof (BEM3DMotionClass),
      (GtsObjectClassInitFunc) bem3d_motion_class_init,
      (GtsObjectInitFunc) bem3d_motion_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &bem3d_motion_info);
  }

  return klass;
}

/** 
 * Generate a new ::BEM3DMotion.
 * 
 * @param klass a ::BEM3DMotionClass; 
 * @param m a ::BEM3DMesh to which the motion equations will be applied;
 * @param m0 a ::BEM3DMesh which will be used as the reference
 * geometry in computing the motion. The ::BEM3DMotion is initialized
 * with the default equations `x = x0', `y = y0', `z = z0', `u = 0',
 * `v = 0', `w = 0', `nx = 0', `ny = 0', `nz = 0'.
 * 
 * @return the new ::BEM3DMotion.
 */

BEM3DMotion * bem3d_motion_new    (BEM3DMotionClass * klass,
				   BEM3DMesh *m, BEM3DMesh *m0)
{
  BEM3DMotion * object;

  object = BEM3D_MOTION (gts_object_new (GTS_OBJECT_CLASS (klass)));

  object->m = m ; object->m0 = m0 ;

  if ( m0 != NULL ) {
    if ( bem3d_mesh_node_number(m) != bem3d_mesh_node_number(m0) )
      g_error("%s: m and m0 do not have the same number of nodes",
	      __FUNCTION__) ;    
  }
 
  return object;
}

/** 
 * Add a variable definition to a ::BEM3DMotion, in the form \a var =
 * \a def. The variable name \a var may not be one of the reserved
 * tokens but may otherwise be any valid libmatheval token. The
 * variables which must be defined (and are given default definitions
 * in ::BEM3DMotion) are \f$x\f$, \f$y\f$, \f$z\f$ which will be
 * evaluated to give the node position and \f$u\f$, \f$v\f$ and
 * \f$w\f$ which will give the surface velocity, other than that due
 * to surface motion (e.g. in scattering problems). Reserved tokens
 * are \f$x0\f$, \f$y0\f$, \f$z0\f$, the position of a node on the
 * reference ::BEM3DMesh, \f$nx\f$, \f$ny\f$, \f$nz\f$, \f$nx0\f$,
 * \f$ny0\f$, \f$nz0\f$, the components of the local surface normal on
 * the mesh and the reference mesh respectively, \f$t\f$, time,
 * \f$i\f$, the node index and \f$j\f$, the complex variable.

 * 
 * @param m a ::BEM3DMotion;
 * @param var the name of the variable to add or redefine;
 * @param def the definition of \a var.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_motion_variable_add(BEM3DMotion *m, gchar *var, gchar *def)

{
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(var != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(def != NULL, BEM3D_NULL_ARGUMENT) ;

  if ( bem3d_motion_token_is_reserved(var) )
    g_error("%s: `%s' is a reserved token and may not be redefined",
	    __FUNCTION__, var) ;  

  if ( g_hash_table_lookup(m->defs, var) != NULL ) 
    g_hash_table_remove(m->defs, var) ;

  g_hash_table_insert(m->defs, var, def) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Find the definition of a variable in a ::BEM3DMotion.
 * 
 * @param m a ::BEM3DMotion.
 * @param var the name of the variable to look up.
 * 
 * @return the definition of \a var or NULL if \a var is not defined
 * in \a m.
 */

gchar *bem3d_motion_variable_lookup(BEM3DMotion *m, gchar *var)

{
  g_return_val_if_fail(m != NULL, NULL) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), NULL) ;
  g_return_val_if_fail(var != NULL, NULL) ;

  return (g_hash_table_lookup(m->defs, var)) ;
}

static void write_defs(gchar *var, gchar *def, FILE *f) 

{
  fprintf(f, "%s = %s\n", var, def) ;

  return ;
}

/** 
 * Write the definitions in a ::BEM3DMotion to file. The file format
 * is a header line `BEM3DMotion', followed by the variable
 * definitions, one per line, in the format: `x = ...'.
 * 
 * @param m the ::BEM3DMotion to write to file;
 * @param f file to which \a m is to be written.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_motion_write(BEM3DMotion *m, FILE *f)

{
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  fprintf(f, "BEM3DMotion\n") ;
  g_hash_table_foreach(m->defs, (GHFunc)write_defs, f) ;

  return BEM3D_SUCCESS ;
}

static void insert_definition(GString *s, gint off, gint len, gchar *def)

{
  g_string_erase(s, off, len) ;
  g_string_insert_c(s, off, '(') ;
  g_string_insert_c(s, off+1, ')') ;
  g_string_insert(s, off+1, def) ;

  return ;
}

static void copy_and_expand_def(BEM3DMotion *m, gchar *var, GString *s,
				GString *buf)

{
  gchar *def ;
  gint i, nc, nex, npass = 0 ;

  def = bem3d_motion_variable_lookup(m, var) ;
  g_string_assign(s, def) ;

  do {
    nex = 0 ;
    for ( i = 0 ; i < s->len ; i ++ ) {
      g_string_assign(buf, "") ;
      nc = strcspn(&(s->str[i]), "+-*/() ") ;
      g_string_append_len(buf, &(s->str[i]), nc) ;
      if ( (nc != 0) && 
	   (def = bem3d_motion_variable_lookup(m, buf->str)) != NULL ) {
	insert_definition(s, i, nc, def) ;
      } else i += nc ;
    }
    npass ++ ;
  } while ( nex != 0 && npass < 64 ) ;

  if ( npass > 63 ) 
    g_error("%s: expression (%s) seems to have a circular definition",
	    __FUNCTION__, bem3d_motion_variable_lookup(m, var)) ;

  return ;
}

/** 
 * Expand the equations for surface variables (displacement and
 * velocity) in a ::BEM3DMotion to give expressions in terms of the
 * reserved tokens only so that the variables can be evaluated using
 * ::bem3d_motion_mesh_position.
 * 
 * @param m the ::BEM3DMotion whose expressions are to be expanded.
 * 
 * @return ::BEM3D_SUCCESS on success.
 */

gint bem3d_motion_expand_defs(BEM3DMotion *m)

{
  GString *buf ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

  buf = g_string_new("") ;

  copy_and_expand_def(m, "x", m->x, buf) ;
  copy_and_expand_def(m, "y", m->y, buf) ;
  copy_and_expand_def(m, "z", m->z, buf) ;
  copy_and_expand_def(m, "u", m->u, buf) ;
  copy_and_expand_def(m, "v", m->v, buf) ;
  copy_and_expand_def(m, "w", m->w, buf) ;

  g_string_free(buf, TRUE) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Write the expanded expressions for displacement and velocity in a
 * ::BEM3DMotion to file. There is no need to do this to store a
 * ::BEM3DMotion: ::bem3d_motion_write should be used instead.
 * 
 * @param m a ::BEM3DMotion;
 * @param f a FILE to which the expansions in \a m should be written.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_motion_write_expansions(BEM3DMotion *m, FILE *f)

{
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  fprintf(f, "x = %s\n", m->x->str) ;
  fprintf(f, "y = %s\n", m->y->str) ;
  fprintf(f, "z = %s\n", m->z->str) ;
  fprintf(f, "u = %s\n", m->u->str) ;
  fprintf(f, "v = %s\n", m->v->str) ;
  fprintf(f, "w = %s\n", m->w->str) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Create the libmatheval (http://www.gnu.org/software/libmatheval/)
 * evaluators for the expressions of a ::BEM3DMotion, prior to
 * evaluating node positions for a ::BEM3DMesh. If BEM3D has been
 * compiled without libmatheval support, this function will fail.
 *
 * @param m a ::BEM3DMotion. 
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_motion_create_evaluators(BEM3DMotion *m)

{
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

#ifndef HAVE_LIBMATHEVAL
  g_error("%s: BEM3D compiled without matheval support", __FUNCTION__) ;
#else /*HAVE_LIBMATHEVAL*/
  if ( (bem3d_motion_evaluator_x(m) = 
	evaluator_create(bem3d_motion_position_func_x(m))) == NULL )
    g_error("%s: evaluator failed for x-displacement function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_x(m)) ;  
  if ( (bem3d_motion_evaluator_y(m) = 
	evaluator_create(bem3d_motion_position_func_y(m))) == NULL )
    g_error("%s: evaluator failed for y-displacement function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_y(m)) ;  
  if ( (bem3d_motion_evaluator_z(m) = 
	evaluator_create(bem3d_motion_position_func_z(m))) == NULL )
    g_error("%s: evaluator failed for y-displacement function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_z(m)) ;  

  if ( (bem3d_motion_evaluator_u(m) = 
	evaluator_create(bem3d_motion_velocity_func_u(m))) == NULL )
    g_error("%s: evaluator failed for u velocity function (%s)",
	    __FUNCTION__, bem3d_motion_velocity_func_u(m)) ;
  if ( (bem3d_motion_evaluator_v(m) = 
	evaluator_create(bem3d_motion_velocity_func_v(m))) == NULL )
    g_error("%s: evaluator failed for v velocity function (%s)",
	    __FUNCTION__, bem3d_motion_velocity_func_v(m)) ;
  if ( (bem3d_motion_evaluator_w(m) = 
	evaluator_create(bem3d_motion_velocity_func_w(m))) == NULL )
    g_error("%s: evaluator failed for w velocity function (%s)",
	    __FUNCTION__, bem3d_motion_velocity_func_w(m)) ;

  if ( (bem3d_motion_evaluator_dx(m) = 
	evaluator_derivative(bem3d_motion_evaluator_x(m), "t")) == NULL )
    g_error("%s: evaluator failed for derivative of x displacement "
	    "function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_x(m)) ;
  if ( (bem3d_motion_evaluator_dy(m) = 
	evaluator_derivative(bem3d_motion_evaluator_y(m), "t")) == NULL )
    g_error("%s: evaluator failed for derivative of y displacement "
	    "function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_y(m)) ;
  if ( (bem3d_motion_evaluator_dz(m) = 
	evaluator_derivative(bem3d_motion_evaluator_z(m), "t")) == NULL )
    g_error("%s: evaluator failed for derivative of z displacement "
	    "function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_z(m)) ;

  if ( (bem3d_motion_evaluator_d2x(m) = 
	evaluator_derivative(bem3d_motion_evaluator_dx(m), "t")) == NULL )
    g_error("%s: evaluator failed for derivative of x velocity "
	    "function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_x(m)) ;
  if ( (bem3d_motion_evaluator_d2y(m) = 
	evaluator_derivative(bem3d_motion_evaluator_dy(m), "t")) == NULL )
    g_error("%s: evaluator failed for derivative of y velocity "
	    "function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_y(m)) ;
  if ( (bem3d_motion_evaluator_d2z(m) = 
	evaluator_derivative(bem3d_motion_evaluator_dz(m), "t")) == NULL )
    g_error("%s: evaluator failed for derivative of z velocity "
	    "function (%s)",
	    __FUNCTION__, bem3d_motion_position_func_z(m)) ;

#endif /*HAVE_LIBMATHEVAL*/  

  return BEM3D_SUCCESS ;
}

/** 
 * Free the function evaluators of a ::BEM3DMotion. Note that this
 * does not affect the stored expressions, only the libmatheval
 * evaluators. This function will fail if BEM3D has been compiled
 * without libmatheval support.
 * 
 * @param m a ::BEM3DMotion.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_motion_free_evaluators(BEM3DMotion *m)

{
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

#ifndef HAVE_LIBMATHEVAL
  g_error("%s: BEM3D compiled without matheval support", __FUNCTION__) ;
#else /*HAVE_LIBMATHEVAL*/

  if ( bem3d_motion_evaluator_x(m) == NULL ) return BEM3D_SUCCESS ;

  evaluator_destroy(bem3d_motion_evaluator_x(m)) ;
  bem3d_motion_evaluator_x(m) = NULL ;
  evaluator_destroy(bem3d_motion_evaluator_y(m)) ;
  bem3d_motion_evaluator_y(m) = NULL ;
  evaluator_destroy(bem3d_motion_evaluator_z(m)) ;
  bem3d_motion_evaluator_z(m) = NULL ;

  evaluator_destroy(bem3d_motion_evaluator_dx(m)) ;
  bem3d_motion_evaluator_dx(m) = NULL ;
  evaluator_destroy(bem3d_motion_evaluator_dy(m)) ;
  bem3d_motion_evaluator_dy(m) = NULL ;
  evaluator_destroy(bem3d_motion_evaluator_dz(m)) ;
  bem3d_motion_evaluator_dz(m) = NULL ;

  evaluator_destroy(bem3d_motion_evaluator_d2x(m)) ;
  bem3d_motion_evaluator_d2x(m) = NULL ;
  evaluator_destroy(bem3d_motion_evaluator_d2y(m)) ;
  bem3d_motion_evaluator_d2y(m) = NULL ;
  evaluator_destroy(bem3d_motion_evaluator_d2z(m)) ;
  bem3d_motion_evaluator_d2z(m) = NULL ;

  evaluator_destroy(bem3d_motion_evaluator_u(m)) ;
  bem3d_motion_evaluator_u(m) = NULL ;
  evaluator_destroy(bem3d_motion_evaluator_v(m)) ;
  bem3d_motion_evaluator_v(m) = NULL ;
  evaluator_destroy(bem3d_motion_evaluator_w(m)) ;
  bem3d_motion_evaluator_w(m) = NULL ;

#endif /*HAVE_LIBMATHEVAL*/

  return BEM3D_SUCCESS ;  
}

static gint position_node(gint i, GtsVertex *v, gpointer data[])

{
  BEM3DMotion *m = data[0] ;  
  gdouble *values = data[1] ;
  GtsVertex *v0 ;
  
  v0 = bem3d_mesh_node_from_index(bem3d_motion_base_mesh(m), i) ;
  values[BEM3D_MOTION_X0] = GTS_POINT(v0)->x ;
  values[BEM3D_MOTION_Y0] = GTS_POINT(v0)->y ;
  values[BEM3D_MOTION_Z0] = GTS_POINT(v0)->z ;
  values[BEM3D_MOTION_INDEX] = (gdouble)i ;
  values[BEM3D_MOTION_SQRT_M1] = 0.0 ;

#ifdef HAVE_LIBMATHEVAL

  GTS_POINT(v)->x = evaluator_evaluate(bem3d_motion_evaluator_x(m),
				       BEM3D_MOTION_NVALUES, 
				       BEM3D_MOTION_VARIABLES,
				       values) ;
  GTS_POINT(v)->y = evaluator_evaluate(bem3d_motion_evaluator_y(m),
				       BEM3D_MOTION_NVALUES, 
				       BEM3D_MOTION_VARIABLES,
				       values) ;
  GTS_POINT(v)->z = evaluator_evaluate(bem3d_motion_evaluator_z(m),
				       BEM3D_MOTION_NVALUES, 
				       BEM3D_MOTION_VARIABLES,
				       values) ;

#endif /*HAVE_LIBMATHEVAL*/

  return 0 ;
}

/** 
 * Position the nodes of a ::BEM3DMesh linked to a ::BEM3DMotion. 
 * 
 * @param m the ::BEM3DMotion containing a mesh;
 * @param t time variable.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_motion_mesh_position(BEM3DMotion *m, gdouble t)

{
  gdouble values[BEM3D_MOTION_NVALUES] ;
  gpointer data[2] ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

#ifndef HAVE_LIBMATHEVAL
  g_error("%s: BEM3D compiled without libmatheval support", __FUNCTION__) ;
#endif /*HAVE_LIBMATHEVAL*/


  values[BEM3D_MOTION_TIME] = t ;
  values[BEM3D_MOTION_NX] = 0.0 ;
  values[BEM3D_MOTION_NY] = 0.0 ;
  values[BEM3D_MOTION_NZ] = 0.0 ;

  data[0] = m ; data[1] = values ;

  bem3d_mesh_foreach_node(bem3d_motion_mesh(m), 
			  (BEM3DNodeFunc)position_node, data) ;

  return BEM3D_SUCCESS ;
}

/** 
 * Check if a token is a reserved word for a ::BEM3DMotion function. 
 * 
 * @param token a string containing a potential token.
 * 
 * @return TRUE if \a token is a reserved word for ::BEM3DMotion
 * functions.
 */

gboolean bem3d_motion_token_is_reserved(gchar *token)

{
  gint i ;
  
  g_return_val_if_fail(token != NULL, FALSE) ;  

  for ( i = 0 ; i < BEM3D_MOTION_NVALUES ; i ++ ) 
    if ( strcmp(token, BEM3D_MOTION_VARIABLES[i]) == 0 )
      return TRUE ;
  
  return FALSE ;
}

/** 
 * Read the definition of a ::BEM3DMotion from file. The file format
 * is that used by ::bem3d_motion_write. 
 * 
 * @param m a ::BEM3DMotion to which the definitions should be added;
 * @param f a GtsFile for the input stream.
 * 
 * @return 0 on success or the line number in the input file where an
 * error was encountered.
 */

guint bem3d_motion_read(BEM3DMotion *m, GtsFile *f)

{
  gchar *var, *def ;
  GString *dstr ;

  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;
  g_return_val_if_fail(f != NULL, BEM3D_NULL_ARGUMENT) ;

  dstr = g_string_new("") ;

  if ( f->type != GTS_STRING ) {
    gts_file_error(f, "expecting a string (BEM3DMotion)") ;
    return f->line ;
  }  

  gts_file_first_token_after(f, '\n') ;
  while ( f->type != GTS_NONE ) {
    if ( f->type != GTS_STRING ) {
      gts_file_error(f, "expecting a string (variable name)") ;
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

    def = g_strdup(dstr->str) ;
    bem3d_motion_variable_add(m, var, def) ;
    gts_file_first_token_after(f, '\n') ;
  }

  return BEM3D_SUCCESS ;  
}

/** 
 * Evaluate the velocity of a node on a mesh from the position
 * specification in a ::BEM3DMotion.
 * 
 * @param m a ::BEM3DMotion which has been initialized and had its
 * evaluators set using bem3d_motion_create_evaluators;

 * @param i index of the node whose velocity is to be computed;
 * @param t time;
 * @param u on exit contains velocity of node i.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_motion_node_velocity(BEM3DMotion *m, gint i, gdouble t,
				GtsVector u)

{
  GtsVertex *v, *v0 ;
  gdouble values[BEM3D_MOTION_NVALUES] ;
  
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

#ifdef HAVE_LIBMATHEVAL
  if ( (v = bem3d_mesh_node_from_index(bem3d_motion_mesh(m), i)) == NULL) 
    g_error("%s: node %d is not on mesh of m", __FUNCTION__, i) ;
  if ( (v0 = bem3d_mesh_node_from_index(bem3d_motion_base_mesh(m), i)) == NULL) 
    g_error("%s: node %d is not on base mesh of m", __FUNCTION__, i) ;

  if ( bem3d_motion_evaluator_x(m) == NULL )
    g_error("%s: evaluators must be created using "
	    "bem3d_motion_create_evaluators", __FUNCTION__) ;

  values[BEM3D_MOTION_TIME] = t ;
  values[BEM3D_MOTION_NX] = 0.0 ;
  values[BEM3D_MOTION_NY] = 0.0 ;
  values[BEM3D_MOTION_NZ] = 0.0 ;

  values[BEM3D_MOTION_X0] = GTS_POINT(v0)->x ;
  values[BEM3D_MOTION_Y0] = GTS_POINT(v0)->y ;
  values[BEM3D_MOTION_Z0] = GTS_POINT(v0)->z ;
  values[BEM3D_MOTION_INDEX] = (gdouble)i ;
  values[BEM3D_MOTION_SQRT_M1] = 0.0 ;

  u[0] = evaluator_evaluate(bem3d_motion_evaluator_dx(m),
			    BEM3D_MOTION_NVALUES, 
			    BEM3D_MOTION_VARIABLES,
			    values) +
    evaluator_evaluate(bem3d_motion_evaluator_u(m),
		       BEM3D_MOTION_NVALUES, 
		       BEM3D_MOTION_VARIABLES,
		       values) ;
  u[1] = evaluator_evaluate(bem3d_motion_evaluator_dy(m),
			    BEM3D_MOTION_NVALUES, 
			    BEM3D_MOTION_VARIABLES,
			    values) +
    evaluator_evaluate(bem3d_motion_evaluator_v(m),
		       BEM3D_MOTION_NVALUES, 
		       BEM3D_MOTION_VARIABLES,
		       values) ;
  u[2] = evaluator_evaluate(bem3d_motion_evaluator_dz(m),
			    BEM3D_MOTION_NVALUES, 
			    BEM3D_MOTION_VARIABLES,
			    values) +
    evaluator_evaluate(bem3d_motion_evaluator_w(m),
		       BEM3D_MOTION_NVALUES, 
		       BEM3D_MOTION_VARIABLES,
		       values) ;
    
  if ( isnan(u[0]) || isnan(u[1]) || isnan(u[2]) )
    g_error("%s: NaN in velocity (%lg, %lg, %lg); check the motion is "
	    "physically valid", __FUNCTION__, u[0], u[1], u[2]) ;
#else /*HAVE_LIBMATHEVAL*/
  g_warning("%s: cannot evaluate source motion (libmatheval not installed)") ;
#endif /*HAVE_LIBMATHEVAL*/

  return BEM3D_SUCCESS ;  
}

/** 
 * Evaluate the acceleration of a node on a mesh from the position
 * specification in a ::BEM3DMotion.
 * 
 * @param m a ::BEM3DMotion which has been initialized and had its
 * evaluators set using bem3d_motion_create_evaluators;

 * @param i index of the node whose acceleration is to be computed;
 * @param t time;
 * @param a on exit contains acceleration of node i.
 * 
 * @return BEM3D_SUCCESS on success.
 */

gint bem3d_motion_node_acceleration(BEM3DMotion *m, gint i, gdouble t,
				    GtsVector a)

{
  GtsVertex *v, *v0 ;
  gdouble values[BEM3D_MOTION_NVALUES] ;
  
  g_return_val_if_fail(m != NULL, BEM3D_NULL_ARGUMENT) ;
  g_return_val_if_fail(BEM3D_IS_MOTION(m), BEM3D_ARGUMENT_WRONG_TYPE) ;

#ifdef HAVE_LIBMATHEVAL
  if ( (v = bem3d_mesh_node_from_index(bem3d_motion_mesh(m), i)) == NULL) 
    g_error("%s: node %d is not on mesh of m", __FUNCTION__, i) ;
  if ( (v0 = bem3d_mesh_node_from_index(bem3d_motion_base_mesh(m), i)) == NULL) 
    g_error("%s: node %d is not on base mesh of m", __FUNCTION__, i) ;

  if ( bem3d_motion_evaluator_x(m) == NULL )
    g_error("%s: evaluators must be created using "
	    "bem3d_motion_create_evaluators", __FUNCTION__) ;

  values[BEM3D_MOTION_TIME] = t ;
  values[BEM3D_MOTION_NX] = 0.0 ;
  values[BEM3D_MOTION_NY] = 0.0 ;
  values[BEM3D_MOTION_NZ] = 0.0 ;

  values[BEM3D_MOTION_X0] = GTS_POINT(v0)->x ;
  values[BEM3D_MOTION_Y0] = GTS_POINT(v0)->y ;
  values[BEM3D_MOTION_Z0] = GTS_POINT(v0)->z ;
  values[BEM3D_MOTION_INDEX] = (gdouble)i ;
  values[BEM3D_MOTION_SQRT_M1] = 0.0 ;

  a[0] = evaluator_evaluate(bem3d_motion_evaluator_d2x(m),
			    BEM3D_MOTION_NVALUES, 
			    BEM3D_MOTION_VARIABLES,
			    values) ;
  a[1] = evaluator_evaluate(bem3d_motion_evaluator_d2y(m),
			    BEM3D_MOTION_NVALUES, 
			    BEM3D_MOTION_VARIABLES,
			    values) ;
  a[2] = evaluator_evaluate(bem3d_motion_evaluator_d2z(m),
			    BEM3D_MOTION_NVALUES, 
			    BEM3D_MOTION_VARIABLES,
			    values) ;
    
  if ( isnan(a[0]) || isnan(a[1]) || isnan(a[2]) )
    g_error("%s: NaN in velocity (%lg, %lg, %lg); check the motion is "
	    "physically valid", __FUNCTION__, a[0], a[1], a[2]) ;
#else /*HAVE_LIBMATHEVAL*/
  g_warning("%s: cannot evaluate source acceleration (libmatheval "
	    "not installed)") ;
#endif /*HAVE_LIBMATHEVAL*/

  return BEM3D_SUCCESS ;
}

/**
 * @}
 * 
 */
