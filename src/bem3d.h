/* bem3d.h
 * 
 * Copyright (C) 2006, 2009, 2018 Michael Carley
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

#ifndef BEM3D_H_INCLUDED
#define BEM3D_H_INCLUDED

/*@file */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <glib.h>
#include <gts.h>
#include <gqr.h>
  
#include <bem3dconfig.h>


  /*configuration options for the maintainer*/
#define _FOREACH_USE_HASH_TABLE_ 1

#ifndef BEM3D_QUADRATURE_CACHING
#define BEM3D_QUADRATURE_CACHING 0
#endif /*BEM3D_QUADRATURE_CACHING*/

  /**
   * @enum bem3d_error
   * @ingroup logging
   * 
   * Return codes for BEM3D functions.
   */

  typedef enum {
    BEM3D_FAILURE = -1,		/**< unspecified failure */
    BEM3D_SUCCESS = 0,		/**< success */
    BEM3D_EINVAL = 1,		/**< invalid parameter */
    BEM3D_ITERMAX = 2,		/**< maximum number of iterations reached */
    BEM3D_UNKNOWN_FORMAT = 3,	/**< unrecognized file format */
    BEM3D_NULL_ARGUMENT = 4,	/**< argument NULL when it is required 
				   not to be */ 
    BEM3D_ARGUMENT_WRONG_TYPE = 5, /**< an argument was of an incompatible 
				      type */
    BEM3D_ARGUMENT_OUT_OF_RANGE = 6, /**< an argument was outside the valid 
					range */
    BEM3D_BLOCK_FULL = 7,	/**< a ::BEM3DMeshData was full */
    BEM3D_FILE_ERROR = 8        /**< file input or output error */
  } bem3d_error ;

  typedef enum {
    GAUSS_LINEAR = 0,
    GAUSS_LOGARITHMIC = 1,
    GAUSS_TRIANGULAR = 2, 
    GAUSS_CHEBYSHEV_2 = 3,
    GAUSS_HYPER = 4
  } quadrature_rule_t ;

  typedef enum {
    BEM3D_GMSH_SCALAR = 0,
    BEM3D_GMSH_VECTOR = 1,
    BEM3D_GMSH_TENSOR = 2
  } bem3d_gmsh_mode_t ;


  /**
   * @addtogroup operators
   * @{
   */
  
  /**
   * @struct BEM3DOperator 
   *
   * Opaque data structure for computation of surface differential
   * operators. This should only be accessed through the relevant
   * functions and macros.
   * @hideinitializer
   */

  typedef struct _BEM3DOperator BEM3DOperator ;

  struct _BEM3DOperator {
    GArray *w, *id ;
    gint nc ;
  } ;

  /**
   * The number of elements in a ::BEM3DOperator \a op, e.g. three
   * elements for a gradient operator, one for a Laplacian.
   * @hideinitializer
   */

#define bem3d_operator_size(op) ((op->nc))

  /**
   * The number of nodes used in computing a ::BEM3DOperator \a op
   * @hideinitializer
   */

#define bem3d_operator_length(op) ((op->id->len))

  /**
   * Global index of the \a i th node forming the operator \a op
   * @hideinitializer
   */

#define bem3d_operator_index(op,i) ((g_array_index(op->id,gint,i)))

  /**
   * Weight for the x-component of \a i th node forming the operator \a op
   * @hideinitializer
   */

#define bem3d_operator_weight(op,i) (&(g_array_index(op->w,gdouble,(op->nc)*i)))

  /**
   * Modes of calculating differential operators and vertex normals
   * from quantities on connected elements.
   */

typedef enum {
    BEM3D_AVERAGE_MWE,		/**< Mean weighted equally */
    BEM3D_AVERAGE_MWA,		/**< Mean weighted by angle */
    BEM3D_AVERAGE_MWSELR,	/**< Mean weighted by sine and edge 
				   length reciprocal*/ 
    BEM3D_AVERAGE_MWAAT,        /**< Mean weighted by areas of adjacent 
				   triangles*/
    BEM3D_AVERAGE_MWELR,        /**< Mean weighted by edge length reciprocals*/
    BEM3D_AVERAGE_MWRELR        /**< Mean weighted by square root of edge 
				   length reciprocals*/
  } BEM3DAverage ;

  /**
   * @}
   * 
   */

#define BEM3D_PARAMETERS_WAVENUMBER     0
#define BEM3D_PARAMETERS_MACH_NUMBER    1
#define BEM3D_PARAMETERS_NORMAL         4
#define BEM3D_PARAMETERS_LAMBDA_REAL    7
#define BEM3D_PARAMETERS_LAMBDA_IMAG    8
#define BEM3D_PARAMETERS_CONDITIONING   9
#define BEM3D_PARAMETERS_JACOBIAN       10
#define BEM3D_PARAMETERS_JACOBIAN_INV   19
#define BEM3D_PARAMETERS_COUPLING_REAL  28
#define BEM3D_PARAMETERS_COUPLING_IMAG  29
#define BEM3D_PARAMETERS_QUADRATURE_TOL 31
  
#define BEM3D_PARAMETERS_GRADIENT        0
  
#define BEM3D_PARAMETERS_REAL_SIZE      64
#define BEM3D_PARAMETERS_INT_SIZE       16
#define BEM3D_PARAMETERS_POINTER_SIZE   8

  /**
   * @typedef BEM3DShapeFunc
   * @ingroup shapefunc
   *
   * BEM3D shape function definition
   * 
   * @param s local coordinate
   * @param t local coordinate
   * @param L array to be filled with the shape function(s) at (s,t)
   * @param dLds array to be filled with the derivatives of L at (s,t)
   * @param dLdt array to be filled with the derivatives of L at (s,t)
   * @param data data to be passed to the shape function
   *
   * @return 0 on success
   */

  typedef gint (* BEM3DShapeFunc)(gdouble s, 
				  gdouble t, 
				  gdouble *L, 
				  gdouble *dLds, 
				  gdouble *dLdt,
				  gpointer data) ;


  /**
   * @struct BEM3DMesh
   * @ingroup mesh
   * Opaque data structure for a BEM3D Mesh. This should only be accessed
   * through the relevant functions and macros.
   */

  typedef struct _BEM3DMesh         BEM3DMesh ;

  struct _BEM3DMesh {
    /*< private >*/
    GtsSurface parent;
    GHashTable *e,
      *f,
      *c ;
/*     GByteArray *node ; */
    gint i0, i1 ;
    /*< public >*/
    /* add extra data here (if public) */
  } ;

  /**
   * @struct BEM3DMeshClass
   * The basic class for a BEM3DMesh.
   * @ingroup mesh
   */

  typedef struct _BEM3DMeshClass    BEM3DMeshClass;

  struct _BEM3DMeshClass {
    /*< private >*/
    GtsSurfaceClass parent_class;

    /*< public >*/
    /* add extra methods here */
  };

  /** 
   * @ingroup mesh
   * Casts \a obj to ::BEM3DMesh
   * 
   * @param obj a GtsObject
   * 
   * @return result of casting \a obj to ::BEM3DMesh.
   * @hideinitializer
   */
#define BEM3D_MESH(obj)            GTS_OBJECT_CAST (obj,		\
						    BEM3DMesh,		\
						    bem3d_mesh_class ())
  /** 
   * @ingroup mesh
   * Casts \a klass to ::BEM3DMeshClass
   * 
   * @param klass class to cast
   * 
   * @hideinitializer
   */
#define BEM3D_MESH_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,	\
							  BEM3DMeshClass, \
							  bem3d_mesh_class())

  /** 
   * @ingroup mesh
   * TRUE if \a obj is from the ::BEM3DMeshClass class
   * 
   * @hideinitializer
   */
#define BEM3D_IS_MESH(obj)         (gts_object_is_from_class \
				    (obj, bem3d_mesh_class ()))

  BEM3DMeshClass *bem3d_mesh_class  (void);
  BEM3DMesh *bem3d_mesh_new    (BEM3DMeshClass * klass,
				GtsFaceClass *face_class,
				GtsEdgeClass *edge_class,
				GtsVertexClass *vertex_class) ;

  /**
   * @ingroup mesh
   * Number of elements in a ::BEM3DMesh
   * 
   * @hideinitializer
   */

#define bem3d_mesh_element_number(m) (g_hash_table_size(m->e))

  /**
   * @ingroup mesh
   * Minimum index of node to be considered in
   * ::bem3d_mesh_foreach_node
   *
   * @hideinitializer
   */

#define bem3d_mesh_node_index_min(m) (m->i0)
  /**
   * @ingroup mesh
   * Maximum index of node to be considered in
   * ::bem3d_mesh_foreach_node
   * 
   * @hideinitializer
   */

#define bem3d_mesh_node_index_max(m) (m->i1)

  /**
   * @struct BEM3DQuadratureRule
   * @ingroup quadrature
   *
   * Quadrature rule for integration on two-dimensional elements
   * containing the nodes and weights normalized to a unit simplex,
   * and a set of `free terms'.
   */

  typedef struct {
    gint nc, nmax, n, nfree, nfree_max, wfree ;
    gdouble *rule ;
    gdouble free_g[32], free_dg[32] ;
  } BEM3DQuadratureRule ;

  /**
   * @typedef BEM3DElement 
   * Data type for BEM3D elements, containing the data
   * for geometry and collocation points. 
   *
   * The ::BEM3DElement type should usually be accessed only through the
   * provided functions and macros. In order to implement new elements,
   * you should look up the ::BEM3DElementBuildFunc type.
   *
   * @ingroup belement
   * 
   */

  typedef struct _BEM3DElement         BEM3DElement;

  struct _BEM3DElement {
    /*< private >*/
    GtsObject parent;

    /*< public >*/
    /* add extra data here (if public) */
    gint nf,
      nv,
      nc,
      ns,
      mo,
      *i,
      *s ;
    gpointer *f,
      *v,
      *c ;
    gdouble *xs,
      *xc ;
    BEM3DShapeFunc shf,		
      cpf ;
    GArray *Imn ;
    gpointer reserved ;		
  };

  typedef struct _BEM3DElementClass    BEM3DElementClass;

  /**
   * @typedef BEM3DElementClass
   * The basic class for a BEM3DElement
   * @ingroup belement
   */

  struct _BEM3DElementClass {
    /*< private >*/
    GtsObjectClass parent_class;

    /*< public >*/
    /* add extra methods here */
  };

#define BEM3D_ELEMENT(obj)            GTS_OBJECT_CAST (obj,		\
						       BEM3DElement,	\
						       bem3d_element_class ())
#define BEM3D_ELEMENT_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,	\
							     BEM3DElementClass,	\
							     bem3d_element_class())
#define BEM3D_IS_ELEMENT(obj)         (gts_object_is_from_class (obj,	\
								 bem3d_element_class ()))

  BEM3DElementClass * bem3d_element_class  (void);

  /**
   * Numerical data associated with a BEM3DMesh
   * 
   * @hideinitializer
   */

  typedef struct {
    gint nd ;
    GHashTable *t ;
    GArray *d ;
  } BEM3DMeshData ;


  
  /**
   * @struct BEM3DMeshDataEntryFunc
   *
   * gint BEM3DMeshDataEntryFunc(gint i, gdouble *d, gint nf, gpointer data)
   *
   * Function for visiting entries in a BEM3DMeshData
   *
   * @param i index of entry in a BEM3DMeshData
   * @param d row of data for block
   * @param nf number of fields in block
   * @param data user data
   */
  
  typedef gint (*BEM3DMeshDataEntryFunc)(gint i, gdouble *d, gint nf,
					 gpointer data) ;
  
  
  /**
   * @struct BEM3DLookupFunc
   *
   * gint BEM3DLookupFunc(gint i, gint j, gpointer ldata, GArray *f,
   * GArray *g) ;
   *
   * Function for looking up data
   *
   * @param i global index of collocation point
   * @param j local index of collocation point
   * @param ldata user data to pass to lookup function
   * @param f solution at collocation point
   * @param g normal derivative of solution at collocation point
   */

  typedef gint (*BEM3DLookupFunc)(gint i, gint j, gpointer ldata, 
				  GArray *f, GArray *g) ;

  /**
   * Solver selection to choose between, for example, direct and Fast
   * Multipole Method in configuration files.
   *
   * @ingroup config
   */

  typedef enum {
    BEM3D_SOLVER_DIRECT = 1,
    BEM3D_SOLVER_FMM    = 2
  } BEM3DSolver ;

  struct _BEM3DParameters {
    gdouble f[BEM3D_PARAMETERS_REAL_SIZE] ;
    gint    n[BEM3D_PARAMETERS_INT_SIZE] ;
    gpointer p[BEM3D_PARAMETERS_POINTER_SIZE] ;
    gpointer user_data ;
  } ;

  /**
   * @struct BEM3DParameters
   * @ingroup gfunc
   *
   * Basic parameter passing for Green's functions et al. Data can be
   * set or extracted using the predefined macros such as
   * ::bem3d_parameters_wavenumber and there is a field for
   * user-defined data which can be accessed using
   * ::bem3d_parameters_user_data.
   * 
   */

  typedef struct _BEM3DParameters BEM3DParameters ;

  /**
   * @typedef BEM3DGreensFunc
   * @ingroup gfunc
   *
   * gint BEM3DGreensFunc(GtsPoint *x, GtsPoint *y, GtsVector n, 
   * BEM3DParameters *p, gdouble *G, gdouble *dGdn)
   *
   * Green's function for various problems. The function should return
   * 0 on success and fill \a G with the Green's function and \a dGdn
   * with its normal derivative, or their equivalents for a particular
   * problem.
   * 
   * @param x field point
   * @param y source point
   * @param n normal at source point
   * @param param parameters to use in Green's function
   * @param G Green's function
   * @param dGdn normal derivative of Green's function
   */

  typedef gint (*BEM3DGreensFunc)(GtsPoint *x, GtsPoint *y,
				  GtsVector n, BEM3DParameters *p,
				  gdouble *G, gdouble *dGdn) ;

  typedef struct _BEM3DGreensFunction BEM3DGreensFunction ;

  /**
   * @ingroup gfunc
   * @struct BEM3DGreensFunction
   *
   * Opaque data structure for Green's functions. 
   *
   * @hideinitializer
   */
  struct _BEM3DGreensFunction {
    BEM3DGreensFunc func ;
    gboolean real ;
    guint precompute ;
    gint nc, size ;
  } ;

  /*used in precompute field to set quantities for pre-calculation before
   calling Green's functions, etc*/
#define BEM3D_GREENS_FUNC_PRECOMPUTE_NORMAL   1 << 0
#define BEM3D_GREENS_FUNC_PRECOMPUTE_JACOBIAN 1 << 1
#define BEM3D_GREENS_FUNC_PRECOMPUTE_LAMBDA   1 << 2
#define BEM3D_GREENS_FUNC_PRECOMPUTE_GRADIENT 1 << 3
  
  /**
   * @ingroup gfunc
   * The ::BEM3DGreensFunc of a ::BEM3DGreensFunction.
   *
   * @param g a pointer to a ::BEM3DGreensFunction.
   *
   * @return the ::BEM3DGreensFunc of \a g.
   * @hideinitializer
   */
#define bem3d_greens_function_func(_g)             ((_g)->func)

  /**
   * @ingroup gfunc
   * Check if a Green's function is real or complex.
   *
   * @param g a pointer to a ::BEM3DGreensFunction.
   *
   * @return TRUE if \a g is real.
   * @hideinitializer
   */
#define bem3d_greens_function_is_real(_g)          ((_g)->real)

  /**
   * @ingroup gfunc
   * Find the number of values returned by a Green's function.
   *
   * @param g a pointer to a ::BEM3DGreensFunction.
   *
   * @return the number of values returned for one call to \a g.
   * @hideinitializer
   */
#define bem3d_greens_function_component_number(_g) ((_g)->nc)

  /**
   * @ingroup gfunc
   * Check if field point normal is to be pre-computed
   *
   * @param g a pointer to a ::BEM3DGreensFunction.
   *
   * @return TRUE if normal must be computed before calling \a g
   * @hideinitializer
   */
#define bem3d_greens_function_compute_normal(_g) \
  ((_g)->precompute & BEM3D_GREENS_FUNC_PRECOMPUTE_NORMAL)

  /**
   * @ingroup gfunc
   * Check if coupling constant for hypersingular formulations is used 
   * and should be computed
   *
   * @param g a pointer to a ::BEM3DGreensFunction.
   *
   * @return TRUE if coupling constant \f$\lambda\f$ is used
   * @hideinitializer
   */
#define bem3d_greens_function_compute_lambda(_g) \
  ((_g)->precompute & BEM3D_GREENS_FUNC_PRECOMPUTE_LAMBDA)

  /**
   * @ingroup gfunc 
   * 
   * Check if gradient operator is to be computed for a Green's
   * function
   *
   * @param g a pointer to a ::BEM3DGreensFunction.
   *
   * @return TRUE if gradient operator is to be computed
   * @hideinitializer
   */
#define bem3d_greens_function_compute_gradient(_g) \
  ((_g)->precompute & BEM3D_GREENS_FUNC_PRECOMPUTE_GRADIENT)
  
  
  /**
   * @ingroup gfunc
   * Precompute field of a Green's function, defining quantities to be
   * calculated before calling the Green's function proper
   *
   * @param g a pointer to a ::BEM3DGreensFunction.
   *
   * @hideinitializer
   */
#define bem3d_greens_function_precompute(_g) ((_g)->precompute)

  /**
   * @ingroup gfunc
   * Check if field point Jacobian matrix is to be pre-computed
   *
   * @param g a pointer to a ::BEM3DGreensFunction.
   *
   * @return TRUE if Jacobian (and inverse) must be computed before calling \a g
   * @hideinitializer
   */
#define bem3d_greens_function_compute_jacobian(_g) \
  ((_g)->precompute & BEM3D_GREENS_FUNC_PRECOMPUTE_JACOBIAN)

  /**
   * @ingroup gfunc
   *
   * Set all fields of a ::BEM3Parameters struct to default values
   * (NULL for pointers, 0 for everything else)
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @hideinitializer
   */
#define bem3d_parameters_init(_p)					\
    do { memset(((_p)->f), 0, sizeof(gdouble)*BEM3D_PARAMETERS_REAL_SIZE) ; \
      memset(((_p)->n), 0, sizeof(gint)*BEM3D_PARAMETERS_INT_SIZE) ;	\
      memset(((_p)->p), 0, sizeof(gpointer)*BEM3D_PARAMETERS_POINTER_SIZE) ; \
      (_p)->f[BEM3D_PARAMETERS_COUPLING_REAL] = 1.0 ;			\
      (_p)->f[BEM3D_PARAMETERS_COUPLING_IMAG] = 0.0 ;			\
      (_p)->f[BEM3D_PARAMETERS_QUADRATURE_TOL] = 1e-3 ;			\
    }									\
    while (0)
    
  /**
   * @ingroup gfunc
   * Wavenumber in a ::BEM3DParameters struct 
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return the value of the wavenumber contained in \a p
   * @hideinitializer
   */

#define bem3d_parameters_wavenumber(_p) ((_p)->f[BEM3D_PARAMETERS_WAVENUMBER]) 

  /**
   * @ingroup gfunc
   * Mach number in a ::BEM3DParameters struct 
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return the value of the Mach number field of \a p
   * @hideinitializer
   */

#define bem3d_parameters_mach_number(_p) ((_p)->f[BEM3D_PARAMETERS_MACH_NUMBER])


#define bem3d_parameters_normal(_p) (&((_p)->f[BEM3D_PARAMETERS_NORMAL])) 

  /**
   * @ingroup gfunc
   *
   * Real part of coupling parameter for Burton and Miller type
   * hypersingular formulations.
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return real part of coupling parameter (address can be cast to
   * *gsl_complex)
   *
   * @hideinitializer
   */

#define bem3d_parameters_lambda_real(_p) ((_p)->f[BEM3D_PARAMETERS_LAMBDA_REAL])

  /**
   * @ingroup gfunc
   *
   * Imaginary part of coupling parameter for Burton and Miller type
   * hypersingular formulations.
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return imaginary part of coupling parameter
   *
   * @hideinitializer
   */

#define bem3d_parameters_lambda_imag(_p) ((_p)->f[BEM3D_PARAMETERS_LAMBDA_IMAG])

#define bem3d_parameters_coupling_real(_p)	\
  ((_p)->f[BEM3D_PARAMETERS_COUPLING_REAL])
#define bem3d_parameters_coupling_imag(_p)	\
  ((_p)->f[BEM3D_PARAMETERS_COUPLING_IMAG])

  /**
   * @ingroup quadrature
   * Quadrature tolerance in a ::BEM3DParameters struct 
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return the tolerance to be passed to ::bem3d_quadrature_parameter
   * @hideinitializer
   */

#define bem3d_parameters_quadrature_tol(_p)		\
      ((_p)->f[BEM3D_PARAMETERS_QUADRATURE_TOL])
  
  /**
   * @ingroup gfunc
   *
   * Conditioning parameter for Burton and Miller type hypersingular
   * formulations (see Wolf and Lele)
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return value of conditioning parameter
   *
   * @hideinitializer
   */
#define bem3d_parameters_conditioning(_p) ((_p)->f[BEM3D_PARAMETERS_CONDITIONING])

  /**
   * @ingroup gfunc
   *
   * Pointer to Jacobian matrix at field point for conversion between
   * physical and local element coordinates
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return pointer to nine-element array of doubles for 3x3 Jacobian
   * \f$\partial x_{i}/\partial \xi_{j}\f$
   *
   * @hideinitializer
   */
#define bem3d_parameters_jacobian(_p) (&((_p)->f[BEM3D_PARAMETERS_JACOBIAN]))

  /**
   * @ingroup gfunc
   *
   * Pointer to inverse Jacobian matrix at field point for conversion
   * between physical and local element coordinates
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return pointer to nine-element array of doubles for 3x3 inverse
   * Jacobian \f$\partial \xi_{i}/\partial x_{j}\f$
   *
   * @hideinitializer
   */
#define bem3d_parameters_jacobian_inverse(_p)		\
  (&((_p)->f[BEM3D_PARAMETERS_JACOBIAN_INV]))
    
  /**
   * @ingroup gfunc
   * Gradient operator in a ::BEM3DParameters struct
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return a ::BEM3DOperator reserved for gradient calculations
   * @hideinitializer
   */

#define bem3d_parameters_gradient(_p) ((_p)->p[BEM3D_PARAMETERS_GRADIENT]) 

  /**
   * @ingroup gfunc
   * User data in a ::BEM3DParameters struct 
   *
   * @param p a pointer to a ::BEM3DParameters struct
   *
   * @return the gpointer reserved for user data in \a p
   * @hideinitializer
   */
  
#define bem3d_parameters_user_data(_p) ((_p)->user_data)

#define BEM3D_WORKSPACE_DOUBLE_ARRAY_NUMBER 8
#define BEM3D_WORKSPACE_GTS_POINT_NUMBER    8
#define BEM3D_WORKSPACE_GQR_RULE_NUMBER     8
  
  /**
   * @typedef BEM3DWorkspace
   * @ingroup work
   *
   * Data type for workspaces passed to functions
   *
   */

  typedef struct _BEM3DWorkspace         BEM3DWorkspace ;

  struct _BEM3DWorkspace {
    BEM3DQuadratureRule *q ; /*quadrature rule*/
    GArray *doubles[BEM3D_WORKSPACE_DOUBLE_ARRAY_NUMBER] ;
    GtsPoint *points[BEM3D_WORKSPACE_GTS_POINT_NUMBER] ;
    gqr_rule_t *gqr[BEM3D_WORKSPACE_GQR_RULE_NUMBER] ;
    gboolean
    used_gqr[BEM3D_WORKSPACE_GQR_RULE_NUMBER],
    used_points[BEM3D_WORKSPACE_GTS_POINT_NUMBER],
      used_doubles[BEM3D_WORKSPACE_DOUBLE_ARRAY_NUMBER] ;
  } ;

  /**
   * @struct BEM3DQuadratureRuleFunc
   * @ingroup quadrature 
   * gint BEM3DQuadratureRuleFunc(GtsPoint *p, BEM3DElement *e,
   * BEM3DQuadratureRule *q, BEM3DGreensFunction gfunc, 
   * BEM3DParameters *param, gpointer data, 
   * BEM3DWorkspace *work)
   *
   * Function for generating quadrature rules
   * 
   * @param p field point for integration;
   * @param e element over which to integrate;
   * @param q quadrature rule to hold points;
   * @param gfunc Green's function to integrate (ignored for some rules);
   * @param param parameters for @a gfunc;
   * @param data user data passed to quadrature generation;
   * @param work a ::BEM3DWorkspace
   * 
   * @return 0 on success
   */

  typedef gint (*BEM3DQuadratureRuleFunc)(GtsPoint *p, BEM3DElement *e,
					  BEM3DQuadratureRule *q, 
					  BEM3DGreensFunction *gfunc,
					  BEM3DParameters *param,
					  gpointer data,
					  BEM3DWorkspace *work) ;					  
  /**
   * Quadrature selection rule structure. 
   * @ingroup quadrature
   */

  typedef struct {
    GPtrArray *f ;
    GArray *sigma, *NM ;
    gdouble tol ;
  } BEM3DQuadratureSelector ;

#define bem3d_quadrature_selector_length(_s) ((_s)->f->len)
#define bem3d_quadrature_selector_rule(_s,_i) \
  (g_ptr_array_index((_s)->f,(_i)))
#define bem3d_quadrature_selector_sigma(_s,_i) \
  (g_array_index((_s)->sigma,gdouble,(_i)))
#define bem3d_quadrature_selector_data(_s,_i) \
  (g_array_index((_s)->NM,gint,(2*(_i))))
#define bem3d_quadrature_selector_tol(_s)   ((_s)->tol)

  /**
   * @struct BEM3DEdge
   * 
   * A data structure for sharp edges on ::BEM3DMesh. It contains an
   * ordered list of mesh nodes with multiple indices, ordered to
   * give the correct orientation about the edge curve.
   *
   * @hideinitializer
   */

typedef struct _BEM3DEdge         BEM3DEdge;

struct _BEM3DEdge {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  /* add extra data here (if public) */
  GArray *i ;
  GPtrArray *e, *v ;
  BEM3DMesh *m ;
};

  /**
   * @struct BEM3DEdgeClass
   * The basic class for a BEM3DEdge.
   * @ingroup edge
   */

typedef struct _BEM3DEdgeClass    BEM3DEdgeClass;

struct _BEM3DEdgeClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/
  /* add extra methods here */
};

#define BEM3D_EDGE(obj)            GTS_OBJECT_CAST (obj,\
					         BEM3DEdge,\
					         bem3d_edge_class ())
#define BEM3D_EDGE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 BEM3DEdgeClass,\
						 bem3d_edge_class())
#define BEM3D_IS_EDGE(obj)         (gts_object_is_from_class (obj,\
						 bem3d_edge_class ()))

  /**
   * @struct BEM3DElementBuildFunc
   * @ingroup elements
   * Build an element from edge and vertex data
   * 
   * @param e an array of edges
   * @param v an array of vertices
   *
   * @return a pointer to the newly built element
   *
   * The function should insert the vertices into the new element and
   * construct new faces and edges as required to build up the
   * element. Note that the ordering is important for the shape functions. 
   */

  typedef BEM3DElement *(*BEM3DElementBuildFunc)(GtsEdge **e, 
						 GtsVertex **v) ;

  /**
   * @struct BEM3DNodeFunc
   *
   * A function which visits collocation points, usually called from
   * ::bem3d_mesh_foreach_node
   * 
   * @param i the collocation point index
   * @param item pointer to the vertex
   * @param data user data to pass to the function 
   *
   * @return 0 on success
   */

  typedef gint (*BEM3DNodeFunc)(gint i, gpointer item, gpointer data) ;

  /**
   * @struct BEM3DELementFunc
   *
   * A function which visits elements, usually called from
   * bem3d_mesh_foreach_element
   * 
   * @param item pointer to the element
   * @param data user data to pass to the function 
   *
   * @return 0 on success
   */

  typedef gint (*BEM3DElementFunc)(BEM3DElement *item, gpointer data) ;

  /**
   * @struct BEM3DEquationFunc
   *
   * A function to insert data in an equation matrix
   *
   * @param i row index 
   * @param j column index 
   * @param G Green's function matrix terms
   * @param dGdn normal derivative terms
   * @param n number of elements in term (e.g. 2 for complex data)
   * @param data user data to pass to function
   *
   * @return 0 on success
   */
  typedef gint (*BEM3DEquationFunc)(gint i, gint j,
				    gdouble *G, gdouble *dGdn,
				    gint n, gpointer data) ;

  /**
   * @struct BEM3DBCFunc
   *
   * A function to set boundary conditions
   * 
   * @param v the vertex at which to set the boundary condition;
   * @param n the normal at that vertex;
   * @param i the global index of the collocation point;
   * @param data user data to pass to the function.
   *
   * @return 0 on success
   */

  typedef gint (*BEM3DBCFunc)(GtsVertex *v, GtsVector n,
			      gint i, gpointer data) ;

  /**
   * @struct BEM3DMeshDataFunc
   *
   * A function to set numerical data for a mesh and insert it in a
   * BEM3DMeshData variable
   *  
   * @param i global index of node
   * @param v node vertex
   * @param data user data to pass to the function
   * @param f array of output data
   * 
   * @return 0 on success
   */

  typedef gint (*BEM3DMeshDataFunc)(gint i, GtsVertex *v, 
				    gpointer data, gdouble *f) ;

  /**
   * @ingroup matrix
   *
   * A basic definition of a 3x3 matrix.
   * 
   */

  typedef gdouble BEM3DMatrix[9] ;

  /**
   * @struct BEM3DMotion
   * 
   * Data structure containing specification of motion of mesh nodes,
   * with facilities for analytical evaluation of quantities such as
   * velocity, for use in setting boundary conditions.
   *
   * @hideinitializer
   */

typedef struct _BEM3DMotion         BEM3DMotion;

struct _BEM3DMotion {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  /* add extra data here (if public) */
  BEM3DMesh *m, *m0 ;
  GHashTable *defs ;
  GString *x, *y, *z, *u, *v, *w ;
  gpointer fx, fy, fz, fdx, fdy, fdz, fd2x, fd2y, fd2z, fvx, fvy, fvz ;
} ;

typedef struct _BEM3DMotionClass    BEM3DMotionClass;

struct _BEM3DMotionClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/
  /* add extra methods here */
};

#define BEM3D_MOTION(obj)            GTS_OBJECT_CAST (obj,\
					         BEM3DMotion,\
					         bem3d_motion_class ())
#define BEM3D_MOTION_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 BEM3DMotionClass,\
						 bem3d_motion_class())
#define BEM3D_IS_MOTION(obj)         (gts_object_is_from_class (obj,\
						bem3d_motion_class ()))

#define bem3d_motion_mesh(mt)              ((mt)->m)
#define bem3d_motion_base_mesh(mt)         ((mt)->m0)
#define bem3d_motion_position_func_x(mt)   ((mt)->x->str)
#define bem3d_motion_position_func_y(mt)   ((mt)->y->str)
#define bem3d_motion_position_func_z(mt)   ((mt)->z->str)
#define bem3d_motion_velocity_func_u(mt)   ((mt)->u->str)
#define bem3d_motion_velocity_func_v(mt)   ((mt)->v->str)
#define bem3d_motion_velocity_func_w(mt)   ((mt)->w->str)
#define bem3d_motion_evaluator_x(mt)       ((mt)->fx)
#define bem3d_motion_evaluator_y(mt)       ((mt)->fy)
#define bem3d_motion_evaluator_z(mt)       ((mt)->fz)
#define bem3d_motion_evaluator_u(mt)       ((mt)->fvx)
#define bem3d_motion_evaluator_v(mt)       ((mt)->fvy)
#define bem3d_motion_evaluator_w(mt)       ((mt)->fvz)
#define bem3d_motion_evaluator_dx(mt)      ((mt)->fdx)
#define bem3d_motion_evaluator_dy(mt)      ((mt)->fdy)
#define bem3d_motion_evaluator_dz(mt)      ((mt)->fdz)
#define bem3d_motion_evaluator_d2x(mt)     ((mt)->fd2x)
#define bem3d_motion_evaluator_d2y(mt)     ((mt)->fd2y)
#define bem3d_motion_evaluator_d2z(mt)     ((mt)->fd2z)

  /**
   * Basic function class which allows the evaluation of analytical
   * functions of mesh data, including differentiation. 
   *
   * @ingroup functions
   */

  typedef struct _BEM3DFunction         BEM3DFunction;

struct _BEM3DFunction {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  /* add extra data here (if public) */
  GHashTable *defs ;
  GPtrArray *functions, *expansions, *evaluators, *definitions ;
  GArray *idx ;
  gchar **vars ;
};

typedef struct _BEM3DFunctionClass    BEM3DFunctionClass;

struct _BEM3DFunctionClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/
  /* add extra methods here */
};

#define BEM3D_FUNCTION(obj)            GTS_OBJECT_CAST (obj,\
					         BEM3DFunction,\
					         bem3d_function_class ())
#define BEM3D_FUNCTION_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 BEM3DFunctionClass,\
						 bem3d_function_class())
#define BEM3D_IS_FUNCTION(obj)         (gts_object_is_from_class (obj,\
						 bem3d_function_class ()))

  BEM3DFunctionClass * bem3d_function_class  (void) ;
  BEM3DFunction * bem3d_function_new    (BEM3DFunctionClass * klass) ;
  gint bem3d_function_add_function(BEM3DFunction *f, gint i, gchar *def) ;

#define bem3d_function_function_number(f)      ((f->functions->len)) 

  typedef gint (*BEM3DReductionFunc)(BEM3DMesh *m,
				     BEM3DMeshData *f,
				     gint *idata, gint *ni,
				     gdouble *ddata, gint *nd,
				     gpointer data) ;
  
  /**
   * @struct BEM3DMeshSkeleton
   * 
   * A data structure which represents a ::BEM3DMesh as an array of
   * points, mainly for use in interfacing to Fast Multipole Method
   * codes.
   *
   * @hideinitializer
   */

  typedef struct _BEM3DMeshSkeleton BEM3DMeshSkeleton ;

  struct _BEM3DMeshSkeleton {
    BEM3DMesh *m ;
    gint nnodes, /*number of mesh nodes*/
      nelem,     /*number of mesh elements*/
      npts,      /*number of points allocated*/
      order,     /*order of source point interpolation*/
      ns,        /*number of source points*/
      nt,        /*number of target (collocation) points*/
      *idx,      /*indices of nodes for source interpolation*/
      ppe ;      /*(maximum) points per element*/
    guint imin,      /*minimum index on the mesh (indexing assumed contiguous)*/
      imax ;      /*minimum index on the mesh (indexing assumed contiguous)*/

    BEM3DAverage anorm ; /*averaging for normal calculations*/
    gdouble *x, *n, *w ;  /*positions and normals, and interpolation weights*/
    GHashTable *e ; /*to tie elements to their point sources*/
  } ;

  /**
   * Selection of different fast multipole implementations linked into
   * the library
   *
   * @ingroup fmm
   */

  typedef enum {
    BEM3D_FMM_FMMLIB3D_1_2 = 1,    /**< fmmlib3d-1.2, 
				https://github.com/zgimbutas/fmmlib3d*/ 
    BEM3D_FMM_WBFMM = 2,          /**< wbfmm,
				     https://github.com/mjcarley/wbfmm*/
  } BEM3DFastMultipole ;

  /**
   * Selection of problem type for fast multipole calculation. Not all
   * FMM solvers will be able to solve all problems.
   *
   * @ingroup fmm
   */

  typedef enum {
    BEM3D_FMM_LAPLACE,          /**< Laplace equation*/
    BEM3D_FMM_HELMHOLTZ         /**< Helmholtz equation*/
  } BEM3DFastMultipoleProblem ;

#define BEM3D_FMM_POINTER_NUMBER 8
#define BEM3D_FMM_INT_NUMBER     8
  
  typedef struct _BEM3DFMMWorkspace BEM3DFMMWorkspace ;

  struct _BEM3DFMMWorkspace {
    BEM3DFastMultipole solver ;
    gpointer p[BEM3D_FMM_POINTER_NUMBER] ;
    /*generic pointers made available to given FMM solvers for
      internal data*/
    gint i[BEM3D_FMM_INT_NUMBER] ;
    gint nda ;            /*number of doubles allocated in d*/
    gdouble *d ;          /*double array for workspaces*/
  } ;


  /**
   * A (non-) matrix type for using Fast Multipole Method summations
   * in solvers.
   *
   * @ingroup fmm
   */

  typedef struct _BEM3DFMMMatrix BEM3DFMMMatrix ;
  struct _BEM3DFMMMatrix {
    BEM3DFastMultipole solver ;
    BEM3DFastMultipoleProblem problem ;
    BEM3DMeshSkeleton *skel ;
    gdouble tol, /*FMM calculation tolerance*/
      *C,        /*diagonal terms*/
      scaleA[2],  /*scaling of FMM output for `A' matrix, multiplying \phi*/
      scaleB[2] ; /*scaling of FMM output for `B' matrix, multiplying d\phi*/
    GArray *gcorr, *dgcorr,  /*correction weights for G and dGdn*/
      *icorr ;               /*correction indices*/
    gint *idxcorr ;          /*node indices into correction arrays*/
  } ;

  /**
   * @struct BEM3DConfiguration
   * 
   * A data structure containing settings for BEM3D calculations,
   * including choice of solver (direct or FMM), quadrature rules, and
   * Green's function.
   *
   * @hideinitializer
   */

typedef struct _BEM3DConfiguration BEM3DConfiguration ;

struct _BEM3DConfiguration {
  GPtrArray               *keys ;
  BEM3DGreensFunction     gfunc ;
  BEM3DQuadratureRuleFunc qrule ;
  BEM3DQuadratureSelector *qdata ;
  BEM3DSolver             solver ;
  BEM3DFastMultipole      fmm ;
  gint                    skel_order, fmm_tree_depth ;
  gboolean diagonal_shortcut ; /*use short cut method to compute
				 diagonal terms*/
  gdouble fmm_radius,
    fmm_tol,
    bc_default_admittance[2],
    quad_tol ;
  GString                 *job ;
  GString                 *gfunc_comment, *qrule_comment ;
} ;

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * 
   * @return the number of GtsFaces on the element \a e
   * @hideinitializer
   */
#define bem3d_element_face_number(e) (e->nf)

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * 
   * @return the number of geometrical vertices on the element \a e
   * @hideinitializer
   */
#define bem3d_element_vertex_number(e) (e->nv)

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * 
   * @return the number of nodes (collocation points) on the element
   * \a e @hideinitializer
   */
#define bem3d_element_node_number(e) (e->nc)

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * 
   * @return the number of corners on the element \a e
   * @hideinitializer
   */
#define bem3d_element_corner_number(e) (e->ns)

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i index of a corner 0 < i < bem3d_element_corner_number(e)
   *
   * @return the GtsVertex at the \a i th corner of \a e. 
   * @hideinitializer
   */
#define bem3d_element_corner(e,i) (e->v[e->s[(i)]])

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i index of a corner 0 < i < bem3d_element_corner_number(e)
   *
   * @return the local index of the \a i th corner of \a e.
   * @hideinitializer
   */
#define bem3d_element_corner_index(e,i) (e->s[(i)])

  /** 
   * @ingroup belement
   * 
   * @param _e a ::BEM3DElement
   * @param _i local index of a collocation point
   *
   * @return the global index of the \a i th collocation point of \a e.
   * @hideinitializer
   */
#define bem3d_element_global_index(_e,_i) (_e->i[(_i)])

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   *
   * @return the geometric ::BEM3DShapeFunc of \a e
   * @hideinitializer
   */
#define bem3d_element_shape_func(e) (e->shf)

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   *
   * @return the collocation ::BEM3DShapeFunc of \a e
   * @hideinitializer
   */
#define bem3d_element_node_func(e) (e->cpf)

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i local index of a collocation point on \a e
   *
   * @return the GtsVertex of the \a i th node (collocation point) of
   * \a e.  @hideinitializer
   */
#define bem3d_element_node(e,i) (e->c[i]) 

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i local index of a node on \a e
   *
   * @return the local coordinate \f$\xi\f$ of the \a i th node
   * of \a e.
   * @hideinitializer
   */

#define bem3d_element_node_xi(e,i) ((e)->xc[2*(i)])

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i local index of a node on \a e
   *
   * @return the local coordinate \f$\eta\f$ of the \a i th node
   * of \a e.
   * @hideinitializer
   */

#define bem3d_element_node_eta(e,i) ((e)->xc[2*(i)+1])

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i local index of a geometric point on \a e
   *
   * @return the GtsVertex of the \a i th geometric point of \a e.
   * @hideinitializer
   */
#define bem3d_element_vertex(e,i) (e->v[(i)])

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i local index of a collocation point on \a e
   *
   * @return the local coordinate \f$\xi\f$ of the \a i th vertex
   * of \a e.
   * @hideinitializer
   */

#define bem3d_element_vertex_xi(e,i) ((e)->xs[2*(i)])

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i local index of a collocation point on \a e
   *
   * @return the local coordinate \f$\eta\f$ of the \a i th vertex
   * of \a e.
   * @hideinitializer
   */

#define bem3d_element_vertex_eta(e,i) ((e)->xs[2*(i)+1])

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param i local index of a GtsFace on \a e
   *
   * @return the \a i th GtsFace of \a e.
   * @hideinitializer
   */
#define bem3d_element_face(e,i) (e->f[(i)])

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   *
   * @return the order of the moments for element \a e.
   * @hideinitializer
   */
#define bem3d_element_moment_order(e) ((e)->mo)

  /** 
   * @ingroup belement
   * 
   * @param e a ::BEM3DElement
   * @param m order of \f$x\f$ in multipole moment;
   * @param n order of \f$y\f$ in multipole moment.
   *
   * @return \f$I_{mn}\f$ for element \a e.
   * @hideinitializer
   */
#define bem3d_element_moment(e,m,n)		\
  (g_array_index((e)->Imn,gdouble,(m)))

#define bem3d_element_edge_vertex_number(e) ((e->n))
#define bem3d_element_edge_vertex(e,_i) (e->v[(_i)])
#define bem3d_element_edge_index_upper(e,_i) (e->i[(2*(_i))])
#define bem3d_element_edge_index_lower(e,_i) (e->i[(2*(_i)+1)])
#define bem3d_element_edge_upper(e) (e->e1)
#define bem3d_element_edge_lower(e) (e->e2)

  /** 
   * @ingroup belement
   * 
   * @param m a ::BEM3DMesh
   * @param el a ::BEM3DElement
   *
   * @return TRUE if \a el is an element of \a m, FALSE otherwise.
   * @hideinitializer
   */

#define bem3d_element_has_parent_mesh(m,el)		\
  ((g_hash_table_lookup(((m)->e), (el)) != NULL))

  /* BEM3DElement *bem3d_element_new    (BEM3DElementClass * klass); */
  BEM3DElement *bem3d_element_new(BEM3DElementClass * klass,
				  gint nf, gint nv, gint nc, gint ns,
				  BEM3DShapeFunc shf, BEM3DShapeFunc cpf) ;
  gint bem3d_element_add_vertex(BEM3DElement *e, gpointer v, gint i) ;
  gboolean bem3d_element_point_inside(BEM3DElement *e, gdouble xi, gdouble eta) ;
  gint bem3d_element_gradient_weights(BEM3DElement *e, gdouble xi, gdouble eta,
				      gdouble *w) ;
  gint bem3d_element_add_node(BEM3DElement *e, gpointer v, gint i) ;
  gint bem3d_element_add_face(BEM3DElement *e, gpointer f, gint i) ;
  gint bem3d_element_write(BEM3DElement *e, 
			   GHashTable *v, GHashTable *t,
			   FILE *f) ;
  gint bem3d_element_set_index(BEM3DElement *e, gint i, gint j) ;
  gint bem3d_element_set_corner(BEM3DElement *e, gint i, gint j)  ;
  gint bem3d_element_find_vertex(BEM3DElement *e, GtsVertex *v) ;
  gint bem3d_element_find_face(BEM3DElement *e, GtsFace *f) ;
  gint bem3d_element_find_node(BEM3DElement *e, GtsVertex *v) ;
  gint bem3d_element_find_index(BEM3DElement *e, gint i) ;
  gdouble bem3d_element_area(BEM3DElement *e, gint ngp) ;
  gint bem3d_element_position(BEM3DElement *e, gdouble *L,
			      GtsPoint *q) ;
  gint bem3d_element_normal(BEM3DElement *e, gdouble *dLds, gdouble *dLdt,
			    GtsVector normal, gdouble *J) ;
  gdouble bem3d_element_jacobian(BEM3DElement *e, gdouble *dLds, 
				 gdouble *dLdt) ;
  gint bem3d_element_jacobian_matrix_normal(BEM3DElement *e,
					    gdouble *dLds, gdouble *dLdt,
					    GtsVector normal, gdouble *J,
					    gdouble *Ji) ;
gint bem3d_element_assemble_equations_direct(BEM3DElement *e,
					     GtsPoint *x, gint ix,
					     BEM3DConfiguration *config,
					     BEM3DParameters *gdata,
					     gdouble *a, gdouble *b,
					     BEM3DWorkspace *work) ;
					     
  gint bem3d_element_assemble_equations(BEM3DElement *e, GtsPoint *x,
					BEM3DConfiguration *config,
					BEM3DParameters *gdata,
					GArray *G, GArray *dGdn,
					BEM3DWorkspace *work) ;
					
  gint bem3d_element_nearest_vertex(BEM3DElement *e, GtsPoint *p,
				    gint *i, gdouble *R) ;
  gint bem3d_element_boundary_nearest_point(BEM3DElement *e, GtsPoint *x,
					    gdouble *xi, gdouble *eta) ;
  gint bem3d_element_edge_nearest_point(BEM3DElement *e, gint i, gint j,
					GtsPoint *x,
					gdouble *xi, gdouble *eta) ;

  GSList *bem3d_elements_common_edges(BEM3DElement *e1, BEM3DElement *e2) ;
  GSList *bem3d_elements_from_vertices(BEM3DMesh *m, GtsVertex *v1, 
				       GtsVertex *v2) ;
  gint bem3d_element_local_vectors(BEM3DElement *e,
				   gdouble dLds[], gdouble dLdt[],
				   GtsVector u1, gdouble *h1,
				   GtsVector u2, gdouble *h2) ;
  gboolean bem3d_element_has_vertex(BEM3DElement *e, GtsVertex *v) ;
  gboolean bem3d_element_has_node(BEM3DElement *e, GtsVertex *v) ;
  gint bem3d_element_replace_vertex(BEM3DElement *e, GtsVertex *v, 
				    GtsVertex *w) ;
  gint bem3d_element_reset_index(BEM3DMesh *m, BEM3DElement *e, gint i, gint j) ;
  gint bem3d_element_slopes(BEM3DElement *e, gdouble *dLds, gdouble *dLdt,
			    GtsVector dxds, GtsVector dxdt) ;
  gint bem3d_element_nearest_point(BEM3DElement *e, GtsPoint *x,
				   gdouble *xi, gdouble *eta, 
				   gboolean constrain) ;
  GSList *bem3d_element_neighbours(BEM3DElement *el, BEM3DMesh *m) ;
  gint bem3d_element_moments_make(BEM3DElement *e, gint H) ;

  BEM3DElement *bem3d_element_from_node(BEM3DMesh *m, GtsVertex *v, gint i) ;
  gint bem3d_element_index_nodes(BEM3DElement *e, BEM3DMesh *m, gint *n) ;
  gint bem3d_element_vertex_is_corner(BEM3DElement *e, GtsVertex *v) ;

  BEM3DElement *bem3d_element_build_t0(GtsEdge **e, GtsVertex **v) ;
  BEM3DElement *bem3d_element_build_t1(GtsEdge **e, GtsVertex **v) ;
  BEM3DElement *bem3d_element_build_t2(GtsEdge **e, GtsVertex **v) ;
  BEM3DElement *bem3d_element_build_t3(GtsEdge **e, GtsVertex **v) ;
  BEM3DElement *bem3d_element_build_q1(GtsEdge **e, GtsVertex **v) ;
  BEM3DElement *bem3d_element_build_q2(GtsEdge **e, GtsVertex **v) ;

  GSList *bem3d_element_common_nodes(BEM3DElement *e1, BEM3DElement *e2) ;
  GtsBBox *bem3d_element_bounding_box(GtsBBoxClass *klass,
				      BEM3DElement *e) ;

  gint bem3d_greens_func_laplace(GtsPoint *x, GtsPoint *y,
				 GtsVector n, BEM3DParameters *p, 
				 gdouble *G, gdouble *dGdn) ;
  gint bem3d_greens_func_helmholtz(GtsPoint *x, GtsPoint *y,
				   GtsVector n, BEM3DParameters *p,
				   gdouble *G, gdouble *dGdn) ;
  gint bem3d_greens_func_helmholtz_hs(GtsPoint *x, GtsPoint *y,
				      GtsVector ny, BEM3DParameters *p,
				      gdouble *G, gdouble *dGdn) ;
  gint bem3d_greens_func_helmholtz_ch(GtsPoint *x, GtsPoint *y,
				      GtsVector ny, BEM3DParameters *p,
				      gdouble *G, gdouble *dGdn) ;
  gint bem3d_greens_func_gradient_laplace(GtsPoint *x, GtsPoint *y,
					  GtsVector n, BEM3DParameters *p,
					  gdouble *G, gdouble *dGdn) ;
  gint bem3d_greens_func_gradient_helmholtz(GtsPoint *x, GtsPoint *y,
					    GtsVector n, BEM3DParameters *p,
					    gdouble *G, gdouble *dGdn) ;
  BEM3DParameters *bem3d_parameters_new(void) ;

  gint bem3d_radiation_func_laplace(GArray *G, GArray *dG,
				    GArray *phi, GArray *dphi,
				    GArray *f, gpointer data) ;
  gint bem3d_radiation_func_helmholtz(GArray *G, GArray *dG,
				      GArray *phi, GArray *dphi,
				      GArray *f, gpointer data) ;
  gint bem3d_mesh_radiation_point(BEM3DMesh *m,
				  BEM3DConfiguration *config,
				  BEM3DParameters *gdata,
				  BEM3DLookupFunc lf, gpointer ldata,
				  GtsPoint *x, GArray *f,
				  BEM3DWorkspace *work) ;
  gint bem3d_element_radiation_point(BEM3DElement *e, 
				     BEM3DConfiguration *config,
				     BEM3DParameters *gdata,	     
				     BEM3DLookupFunc lfunc, gpointer ldata,
				     GtsPoint *x, 
				     GArray *G, GArray *dGdn,
				     GArray *phi, GArray *dphi,
				     BEM3DWorkspace *work) ;
  gint bem3d_mesh_radiation_mesh(BEM3DMesh *m,
				 BEM3DConfiguration *config,
				 BEM3DParameters *gdata,
				 BEM3DLookupFunc lf, gpointer ldata,
				 BEM3DMesh *s, BEM3DMeshData *f,
				 BEM3DWorkspace *work) ;

  GtsVertex *bem3d_mesh_node_from_index(BEM3DMesh *m, gint i) ;
  gint bem3d_mesh_index_from_node(BEM3DMesh *m, GtsVertex *v) ;
  gint bem3d_mesh_write(BEM3DMesh *m, FILE *fptr) ;
  guint bem3d_mesh_read(BEM3DMesh *m, GtsFile *f) ;
  guint bem3d_gmsh_read(BEM3DMesh *m, FILE *f) ;
  gint bem3d_mesh_add_element(BEM3DMesh *m, BEM3DElement *e, gboolean force) ;
  gint bem3d_mesh_remove_element(BEM3DMesh *m, BEM3DElement *e) ;
  GSList *bem3d_mesh_vertex_elements(BEM3DMesh *m, GtsVertex *v) ;
  GSList *bem3d_mesh_node_elements(BEM3DMesh *m, gint i) ;
  gint bem3d_mesh_index_nodes(BEM3DMesh *m, gdouble angle, gint n) ;
  gint bem3d_mesh_node_number(BEM3DMesh *m) ;
  gint bem3d_mesh_discretize(GtsSurface *s, gint nne,
			     BEM3DElementBuildFunc bfunc,
			     BEM3DMesh *m) ;
  gint bem3d_mesh_assemble_equations(BEM3DMesh *m, BEM3DMesh *n,
				     BEM3DConfiguration *config,
				     BEM3DParameters *gdata,
				     BEM3DEquationFunc efunc, gpointer edata,
				     BEM3DWorkspace *work) ;

  gint bem3d_equation_func_simple(gint i, gint j,
				  gdouble *G, gdouble *dGdn, gint n,
				  gpointer *e) ;
  gint bem3d_mesh_write_nodes(BEM3DMesh *m, FILE *f) ;
  gint bem3d_mesh_foreach_element(BEM3DMesh *m, BEM3DElementFunc f, 
				  gpointer data) ;
  gint bem3d_mesh_foreach_node(BEM3DMesh *m, BEM3DNodeFunc f, gpointer data) ;
  gint bem3d_mesh_element_clear_reserved(BEM3DMesh *m) ;
  gint bem3d_mesh_set_bc(BEM3DMesh *m, BEM3DBCFunc bcf, gpointer bdata) ;
  gint bem3d_mesh_quad_dgdn(BEM3DMesh *m,
			    BEM3DConfiguration *config,
			    BEM3DParameters *gdata,
			    BEM3DLookupFunc lfunc, gpointer ldata,
			    BEM3DEquationFunc efunc, gpointer edata,
			    BEM3DWorkspace *work) ;
  gint bem3d_mesh_disjoint_quad_dgdn(BEM3DMesh *m1,
				     BEM3DMesh *m2,
				     BEM3DConfiguration *config,
				     BEM3DParameters *gdata,
				     BEM3DLookupFunc lfunc, gpointer ldata,
				     BEM3DEquationFunc efunc, gpointer edata,
				     BEM3DWorkspace *work) ;

  gint bem3d_mesh_merge(BEM3DMesh *m, BEM3DMesh *n) ;
  gint bem3d_mesh_element_moments(BEM3DMesh *m, gint H) ;
  gint bem3d_mesh_element_node_number_max(BEM3DMesh *m) ;
  gint bem3d_mesh_index_range(BEM3DMesh *m, guint *imin, guint *imax) ;
  BEM3DElement *bem3d_mesh_element_sample(BEM3DMesh *m) ;

  BEM3DElement *bem3d_mesh_face_element(BEM3DMesh *m, GtsFace *f) ;
  gdouble bem3d_mesh_surface_area(BEM3DMesh *m, gint ngp) ;

  /**
   * The number of elements (numerical entries) per node in a ::BEM3DMeshData.
   * @hideinitializer
   */

#define bem3d_mesh_data_element_number(m) ((m)->nd)

/*   /\** */
/*    * The number of nodes in a ::BEM3DMeshData. */
/*    * @hideinitializer */
/*    *\/ */

/* #define bem3d_mesh_data_node_number(m) ((m)->nc) */

  gint bem3d_mesh_data_multiproc_sum(BEM3DMeshData *m) ;
  BEM3DMeshData *bem3d_mesh_data_new(BEM3DMesh *m, gint n) ;
  gint bem3d_mesh_data_expand(BEM3DMeshData *d, gint ne) ;
  gint bem3d_mesh_data_free(BEM3DMeshData *d) ;
  BEM3DMeshData *bem3d_mesh_data_sized_new(gint n, gint m) ;
  gint bem3d_mesh_data_add_node(BEM3DMeshData *m, gint i) ;
  gdouble *bem3d_mesh_data_get(BEM3DMeshData *m, gint i) ;
  gint bem3d_mesh_data_node_number(BEM3DMeshData *d) ;
  gint bem3d_mesh_data_add_mesh(BEM3DMeshData *d, BEM3DMesh *m) ;
  gint bem3d_mesh_data_foreach(BEM3DMeshData *f, BEM3DMeshDataEntryFunc func,
			       gpointer data) ;

  gint bem3d_mesh_data_clear(BEM3DMeshData *m) ;
  gint bem3d_mesh_data_add(BEM3DMeshData *f, gint i, GArray *g) ;

  BEM3DQuadratureSelector *bem3d_quadrature_selector_new(void) ;
  gint bem3d_quadrature_selector_add(BEM3DQuadratureSelector *s,
				     BEM3DQuadratureRuleFunc f,
				     gdouble p, gint N, gint M) ;
  BEM3DQuadratureSelector *bem3d_quadrature_selector_default(void) ;
  BEM3DQuadratureSelector *bem3d_quadrature_selector_hypersingular_default(void) ;
  gint bem3d_quadrature_selector_clear(BEM3DQuadratureSelector *s) ;
  gint bem3d_quadrature_select(BEM3DQuadratureSelector *s,
			       gdouble p,
			       BEM3DQuadratureRuleFunc *f,
			       gpointer *data) ;

  gint bem3d_shfunc_t0(gdouble s, gdouble t, 
		       gdouble *L, gdouble *dLds,
		       gdouble *dLdt, gpointer data) ;
  gint bem3d_shfunc_t1(gdouble s, gdouble t, 
		       gdouble *L, gdouble *dLds,
		       gdouble *dLdt, gpointer data) ;
  gint bem3d_shfunc_t2(gdouble s, gdouble t, 
		       gdouble *L, gdouble *dLds,
		       gdouble *dLdt, gpointer data) ;
  gint bem3d_shfunc_t3(gdouble s, gdouble t, 
		       gdouble *L, gdouble *dLds,
		       gdouble *dLdt, gpointer data) ;
  gint bem3d_shfunc_q1(gdouble s, gdouble t, 
		       gdouble *L, gdouble *dLds,
		       gdouble *dLdt, gpointer data) ;
  gint bem3d_shfunc_q2(gdouble s, gdouble t, 
		       gdouble *L, gdouble *dLds,
		       gdouble *dLdt, gpointer data) ;

  gint bem3d_shapefunc_lookup_init(void) ;
  gint bem3d_shapefunc_lookup_add(BEM3DShapeFunc func, gchar *name) ;
  BEM3DShapeFunc bem3d_shapefunc_lookup_func(const gchar *name) ;
  const gchar *bem3d_shapefunc_lookup_name(BEM3DShapeFunc func) ;

  gboolean bem3d_edge_is_sharp(GtsEdge *e, gdouble angle) ;
  GSList *bem3d_edges_local(GtsVertex *v, GtsSurface *s) ;

  gint bem3d_lookup_func_unit(gint i, gint j, 
			      gpointer data,
			      GArray *s, GArray *ds) ;

  gint bem3d_lookup_func_unit_c(gint i, gint j, 
				gpointer data,
				GArray *s, GArray *ds) ;
  gint bem3d_lookup_func_both_unit(gint i, gint j, 
				   gpointer data,
				   GArray *s, GArray *ds) ;

  gint bem3d_area_coordinates_tri(GtsPoint *t1,
				  GtsPoint *t2,
				  GtsPoint *t3,
				  GtsPoint *x,
				  GtsVector xi) ;

  gint bem3d_logging_init(FILE *f, gchar *p, 
			  GLogLevelFlags log_level,
			  gpointer exit_func) ;

#define bem3d_quadrature_clear(_q) ((_q)->n=(_q)->nfree=0) 
#define bem3d_quadrature_vertex_number(_q) ((_q)->n)
#define bem3d_quadrature_vertex_number_max(_q) ((_q)->nmax)
#define bem3d_quadrature_component_number(_q) ((_q)->nc)
#define bem3d_quadrature_xi(_q,_i) ((_q)->rule[(3*(_i)+0)])
#define bem3d_quadrature_eta(_q,_i) ((_q)->rule[(3*(_i)+1)])
#define bem3d_quadrature_weight(_q,_i) ((_q)->rule[(3*(_i)+2)])

#define bem3d_quadrature_free_number(_q) ((_q)->nfree)
#define bem3d_quadrature_free_term_g(_q,_i) ((_q)->free_g[(_q)->wfree*(_i)])
#define bem3d_quadrature_free_term_dg(_q,_i) ((_q)->free_dg[(_q)->wfree*(_i)])

  BEM3DQuadratureRule *bem3d_quadrature_rule_new(gint n, gint nc) ;
  gint bem3d_quadrature_rule_free(BEM3DQuadratureRule *q) ;
  gdouble bem3d_quadrature_rule_sum_weights(BEM3DQuadratureRule *q) ;
  gint bem3d_quadrature_rule_write(BEM3DQuadratureRule *q, FILE *f) ;
  gint bem3d_quadrature_rule_realloc(BEM3DQuadratureRule *q, gint n) ;
  gint bem3d_quadrature_add_point(BEM3DQuadratureRule *q,
				  gdouble xi, gdouble eta,
				  gdouble w) ;
  gint bem3d_quadrature_rule_gauss(GtsPoint *p, BEM3DElement *e,
				   BEM3DQuadratureRule *q, 
				   BEM3DGreensFunction *gfunc, 
				   BEM3DParameters *param,
				   gpointer n,
				   BEM3DWorkspace *work) ;
				  
  gint bem3d_quadrature_rule_default(GtsPoint *p,
				     BEM3DElement *e,
				     BEM3DQuadratureRule *q,
				     BEM3DGreensFunction *gfunc, 
				     BEM3DParameters *param,
				     gpointer data,
				     BEM3DWorkspace *work) ;
				     
  gint bem3d_quadrature_rule_kw(GtsPoint *p, BEM3DElement *e,
				BEM3DQuadratureRule *q, 
				BEM3DGreensFunction *gfunc, 
				BEM3DParameters *param,
				gpointer data,
				BEM3DWorkspace *work) ;
				
  gint bem3d_quadrature_rule_polar(GtsPoint *p, BEM3DElement *e,
				   BEM3DQuadratureRule *q, 
				   BEM3DGreensFunction *gfunc, 
				   BEM3DParameters *param,
				   gpointer data,
				   BEM3DWorkspace *work) ;
				   
  gint bem3d_quadrature_rule_polar_hs(GtsPoint *p, BEM3DElement *e,
				      BEM3DQuadratureRule *q, 
				      BEM3DGreensFunction *gfunc, 
				      BEM3DParameters *param,
				      gpointer data,
				      BEM3DWorkspace *work) ;
				      
  gint bem3d_quadrature_rule_hayami(GtsPoint *xs, BEM3DElement *e,
				    BEM3DQuadratureRule *q, 
				    BEM3DGreensFunction *gfunc, 
				    BEM3DParameters *param,
				    gpointer data,
				   BEM3DWorkspace *work) ;
				    
  gint bem3d_quadrature_rule_rnvr(GtsPoint *p, BEM3DElement *e,
				  BEM3DQuadratureRule *q, 
				  BEM3DGreensFunction *gfunc, 
				  BEM3DParameters *param,
				  gpointer data,
				   BEM3DWorkspace *work) ;
				  
  gint bem3d_quadrature_rule_wx(GtsPoint *p, BEM3DElement *e,
				BEM3DQuadratureRule *q, 
				BEM3DGreensFunction *gfunc, 
				BEM3DParameters *param,
				gpointer data,
				BEM3DWorkspace *work) ;
				
  gint bem3d_quadrature_rule_newman(GtsPoint *xs, BEM3DElement *e,
				    BEM3DQuadratureRule *q, 
				    BEM3DGreensFunction *gfunc, 
				    BEM3DParameters *param,
				    gpointer data,
				   BEM3DWorkspace *work) ;
				    
  gint bem3d_quadrature_rule_newman_gradient(GtsPoint *xs, BEM3DElement *e,
					     BEM3DQuadratureRule *q, 
					     BEM3DGreensFunction *gfunc, 
					     BEM3DParameters *param,
					     gpointer data,
				   BEM3DWorkspace *work) ;
					     
  gint bem3d_quadrature_rule_series(GtsPoint *xs, BEM3DElement *e,
				    BEM3DQuadratureRule *q, 
				    BEM3DGreensFunction *gfunc,
				    BEM3DParameters *param,
				    gpointer data,
				   BEM3DWorkspace *work) ;
				    
  gint bem3d_quadrature_rule_mzht(GtsPoint *xs, BEM3DElement *e,
				  BEM3DQuadratureRule *q, 
				  BEM3DGreensFunction *gfunc,
				  BEM3DParameters *param,
				  gpointer data,
				  BEM3DWorkspace *work) ;
				  

gint bem3d_quadrature_rule_decomp(GtsPoint *xs, BEM3DElement *e,
				  BEM3DQuadratureRule *q, 
				  BEM3DGreensFunction *gfunc,
				  BEM3DParameters *param,
				  gpointer data,
				  BEM3DWorkspace *work) ;
				  
gint bem3d_quadrature_rule_decomp_gradient(GtsPoint *xs, BEM3DElement *e,
					   BEM3DQuadratureRule *q, 
					   BEM3DGreensFunction *gfunc,
					   BEM3DParameters *param,
					   gpointer data,
					   BEM3DWorkspace *work) ;
					   
  gdouble bem3d_quadrature_parameter(GtsPoint *p, BEM3DElement *e,
				     gdouble tol) ;
  gint bem3d_quadrature_rule_remap(gdouble xi0, gdouble eta0,
				   gdouble xi1, gdouble eta1,
				   gdouble xi2, gdouble eta2,
				   gdouble s, gdouble t, gdouble w,
				   gdouble *sn, gdouble *tn,
				   gdouble *wn) ;

  gint bem3d_geometry_plane(GtsSurface *s, gint ni, gint nj) ;
  GtsVertex *bem3d_vertex_from_segments(GtsSegment *s1, GtsSegment *s2) ;
  gint bem3d_mesh_data_write(BEM3DMeshData *f, FILE *fp, gchar *header) ;
  gint bem3d_mesh_data_write_weights(BEM3DMeshData *f, FILE *fp,
				     gint *fields, gint nfields) ;
  gint bem3d_mesh_data_read(BEM3DMeshData **f, FILE *fp, gint width) ;
  BEM3DMeshData *bem3d_mesh_data_merge(BEM3DMeshData *f1,
				       BEM3DMeshData *f2,
				       gboolean strict) ;
  gint bem3d_mesh_write_msh(BEM3DMesh *m, BEM3DMeshData *f, gint k,
			    gchar *view, gdouble t, bem3d_gmsh_mode_t mode,
			    gint offset, FILE *fp) ;
  gint bem3d_mesh_write_pos(BEM3DMesh *m, BEM3DMeshData *f, gint k,
			    gchar *view, bem3d_gmsh_mode_t mode,
			    FILE *fp) ;
  gint bem3d_element_write_pos(BEM3DElement *e, BEM3DMeshData *f,
			       gint nv, 
			       gint *indices_g, gint *indices_d,
			       gchar *ename,
			       gint k, gint nf, FILE *fp) ;
  gint bem3d_edge_write_pos(BEM3DEdge *edge, gchar *view, FILE *f) ;

  GSList *bem3d_mesh_index_sharp_edges(BEM3DMesh *m, gdouble angle, 
				       gint *n) ;

  gint bem3d_matrix_vector_mul(BEM3DMatrix m, GtsVector v, GtsVector w) ;
  gdouble bem3d_matrix_det(BEM3DMatrix m) ;
  gint bem3d_matrix_inverse(BEM3DMatrix m, BEM3DMatrix im) ;

  BEM3DOperator *bem3d_operator_new(void) ;
  gint bem3d_operator_gradient(BEM3DMesh *m, gint i, BEM3DOperator *op,
			       BEM3DAverage mode) ;
  gint bem3d_node_normal(BEM3DMesh *m, gint i, GtsVector n, 
			 BEM3DAverage mode) ;


#define bem3d_edge_node_number(edge)  ((edge->i->len)/2)
#define bem3d_edge_element_number(edge)  ((edge->e->len)/2)
#define bem3d_edge_node_index_upper(edge,j)  \
  ((g_array_index(edge->i,gint,2*(j)+0)))
#define bem3d_edge_node_index_lower(edge,j)  \
  ((g_array_index(edge->i,gint,2*(j)+1)))
#define bem3d_edge_element_upper(edge,j)  \
  ((g_ptr_array_index(edge->e,2*(j)+0)))
#define bem3d_edge_element_lower(edge,j)  \
  ((g_ptr_array_index(edge->e,2*(j)+1)))
#define bem3d_edge_vertex(edge,j)  ((g_ptr_array_index(edge->v,(j))))
#define bem3d_edge_mesh(edge)      (edge)->m

BEM3DEdgeClass * bem3d_edge_class  (void) ;
BEM3DEdge * bem3d_edge_new    (BEM3DEdgeClass * klass) ;

gint bem3d_mesh_vertex_index_number(BEM3DMesh *m, GtsVertex *v) ;
GSList *bem3d_mesh_sharp_vertices(BEM3DMesh *m) ;
gint bem3d_edge_add_node(BEM3DEdge *e, GtsVertex *v, gint i, gint j) ;
gint bem3d_edge_add_elements(BEM3DEdge *e, 
			     BEM3DElement *e1, BEM3DElement *e2) ;
GSList *bem3d_mesh_extract_edges(BEM3DMesh *m, GSList **e) ;
gint bem3d_edge_clear(BEM3DEdge *e) ;
  gint bem3d_edge_copy(BEM3DEdge *e, BEM3DEdge *f) ;
gint bem3d_link_edges(BEM3DMesh *m, BEM3DEdge *e1, BEM3DEdge *e2,
		      BEM3DElementBuildFunc build, gint nv, 
		      gint nc, gint *idx) ;
  gint bem3d_edge_write(BEM3DEdge *e, FILE *f) ;
  guint bem3d_edge_read(BEM3DEdge *e, GtsFile *f) ;
  gint bem3d_edge_link_to_mesh(BEM3DEdge *e, BEM3DMesh *m) ;
  gboolean bem3d_edge_is_oriented(BEM3DEdge *e, BEM3DMesh *m) ;
  gint bem3d_invert_edge(BEM3DEdge *e) ;

BEM3DMotionClass * bem3d_motion_class  (void);
BEM3DMotion * bem3d_motion_new    (BEM3DMotionClass * klass,
				   BEM3DMesh *m, BEM3DMesh *m0);

gint bem3d_motion_variable_add(BEM3DMotion *m, gchar *var, gchar *def) ;
gint bem3d_motion_write(BEM3DMotion *m, FILE *f) ;
gint bem3d_motion_expand_defs(BEM3DMotion *m) ;
gint bem3d_motion_write_expansions(BEM3DMotion *m, FILE *f) ;
gchar *bem3d_motion_variable_lookup(BEM3DMotion *m, gchar *var) ;
gint bem3d_motion_mesh_position(BEM3DMotion *m, gdouble t) ;
gint bem3d_motion_create_evaluators(BEM3DMotion *m) ;
gint bem3d_motion_free_evaluators(BEM3DMotion *m) ;
gboolean bem3d_motion_token_is_reserved(gchar *token) ;
guint bem3d_motion_read(BEM3DMotion *m, GtsFile *f) ;
gint bem3d_motion_node_velocity(BEM3DMotion *m, gint i, gdouble t,
				GtsVector u) ;
gint bem3d_motion_node_acceleration(BEM3DMotion *m, gint i, gdouble t,
				    GtsVector a) ;


  gint bem3d_function_variable_add(BEM3DFunction *f, gchar *var, gchar *def) ;
  gint bem3d_function_expand_functions(BEM3DFunction *f) ;
  gchar *bem3d_function_variable_lookup(BEM3DFunction *f, gchar *var) ;
  gint bem3d_function_apply_motion(BEM3DFunction *f, 
				   BEM3DMotion *m,
				   gdouble t,
				   BEM3DMeshData *d,
				   BEM3DMeshData *e) ;
  gint bem3d_function_apply_motion_list(BEM3DFunction *func, 
					GPtrArray *motions,
					gdouble t,
					BEM3DMeshData *f,
					BEM3DMeshData *g) ;
  gint bem3d_function_apply_mesh(BEM3DFunction *f, 
				 BEM3DMesh *m,
				 BEM3DMeshData *d,
				 BEM3DMeshData *e) ;
  gint bem3d_function_apply_mesh_list(BEM3DFunction *func, 
				      GPtrArray *meshes,
				      BEM3DMeshData *f,
				      BEM3DMeshData *g) ;
  gint bem3d_function_integral_weights(BEM3DFunction *func,
				       BEM3DMesh *m, gint imesh,
				       GtsVertex *x, GtsVector n, gint i,
				       BEM3DQuadratureRule *q,
				       BEM3DMeshData *f) ;
  gint bem3d_function_eval_point(BEM3DFunction *func,
				 GtsPoint *x, GtsVector n, gint i,
				 gdouble *result, gint nres) ;
  
  gboolean bem3d_function_token_is_reserved(gchar *var) ;
  gint bem3d_function_write(BEM3DFunction *f, FILE *fid) ; 
  gint bem3d_function_read(BEM3DFunction *f, GtsFile *fid) ;
  gint bem3d_function_insert_string(BEM3DFunction *fn, gchar *str) ;

  gint bem3d_reduction_func_max(BEM3DMesh *m, BEM3DMeshData *f,
				gint *idata, gint *ni,
				gdouble *ddata, gint *nd,
				gpointer data) ;
  gint bem3d_reduction_func_min(BEM3DMesh *m, BEM3DMeshData *f,
				gint *idata, gint *ni,
				gdouble *ddata, gint *nd,
				gpointer data) ;
  gint bem3d_reduction_func_limits(BEM3DMesh *m, BEM3DMeshData *f,
				   gint *idata, gint *ni,
				   gdouble *ddata, gint *nd,
				   gpointer data) ;
  gint bem3d_reduction_func_sum(BEM3DMesh *m, BEM3DMeshData *f,
				gint *idata, gint *ni,
				gdouble *ddata, gint *nd,
				gpointer data) ;
  gint bem3d_reduction_func_int(BEM3DMesh *m, BEM3DMeshData *f,
				gint *idata, gint *ni,
				gdouble *ddata, gint *nd,
				gpointer data) ;
  gint bem3d_reduction_func_apply(gchar *func,
				  BEM3DMesh *m, BEM3DMeshData *f,
				  gint *idata, gint *ni,
				  gdouble *ddata, gint *nd,
				  gpointer data) ;
  gchar *bem3d_reduction_func_description(gchar *func) ;
  gint bem3d_reduction_func_list(FILE *output, gchar *format,
				 gboolean describe) ;
  
  BEM3DWorkspace *bem3d_workspace_new(void) ;
  GArray *bem3d_workspace_double_array_get(BEM3DWorkspace *w) ;
  gint bem3d_workspace_double_array_put(BEM3DWorkspace *w, GArray *g) ; 
  GtsPoint *bem3d_workspace_gts_point_get(BEM3DWorkspace *w) ;
  gint bem3d_workspace_gts_point_put(BEM3DWorkspace *w, GtsPoint *p) ;
  gqr_rule_t *bem3d_workspace_gqr_rule_get(BEM3DWorkspace *w) ;
  gint bem3d_workspace_gqr_rule_put(BEM3DWorkspace *w, gqr_rule_t *g) ;
 
  BEM3DConfiguration *bem3d_configuration_new(void) ;
  gint bem3d_configuration_init(void) ;
  gint bem3d_configuration_read(BEM3DConfiguration *c, gchar *file) ;
  gint bem3d_configuration_add_identifier(const gchar *id,
					  gpointer v) ;
  gchar *bem3d_configuration_identifier_from_pointer(gpointer v) ;
  gpointer bem3d_configuration_pointer_from_identifier(const gchar *id) ;  
  gint bem3d_configuration_set(BEM3DConfiguration *c, gpointer k)  ;
  gint bem3d_configuration_write(BEM3DConfiguration *config,
				 gboolean write_descriptions,
				 gint line_width, FILE *f) ;
  gint bem3d_configuration_write_generic(gboolean write_descriptions,
					 gint line_width, gchar *prefix,
					 gchar *comment,
					 FILE *f) ;

  BEM3DMeshSkeleton *bem3d_mesh_skeleton_new(BEM3DMesh *m, gint order_max) ;
  gint bem3d_mesh_skeleton_init(BEM3DMeshSkeleton *s, BEM3DQuadratureRule *q,
				BEM3DAverage anorm) ;
  gint bem3d_mesh_skeleton_write(BEM3DMeshSkeleton *s, FILE *f) ;
  gint bem3d_mesh_skeleton_read(BEM3DMeshSkeleton *s, FILE *f) ;

  gint bem3d_fmm_calculate(BEM3DFastMultipole solver,
			   BEM3DFastMultipoleProblem problem,
			   BEM3DParameters *param,
			   BEM3DMeshSkeleton *s, 
			   gdouble tol,
			   gdouble *q, gdouble *dq,
			   gdouble *p, gdouble *dp,
			   BEM3DFMMWorkspace *w) ;
  
  BEM3DFMMWorkspace *bem3d_fmm_workspace_alloc(BEM3DFastMultipole solver,
					       BEM3DMeshSkeleton *skel,
					       BEM3DConfiguration *config) ;
					       
  BEM3DFMMMatrix *bem3d_fmm_matrix_new(BEM3DFastMultipole solver,
				       BEM3DFastMultipoleProblem problem,
				       BEM3DMeshSkeleton *skel,
				       BEM3DConfiguration *config,
				       BEM3DParameters *param,
				       gdouble r, BEM3DWorkspace *work) ;

  gchar *bem3d_solver_name(BEM3DSolver s) ;
  BEM3DSolver bem3d_solver_type(gchar *s) ;
  gchar *bem3d_fmm_name(BEM3DFastMultipole solver) ;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*BEM3D_H_INCLUDED*/
