#ifndef included_PflotranJacobianMultilevelOperator
#define included_PflotranJacobianMultilevelOperator


#ifndef included_Array
#define included_Array
#include "tbox/Array.h"
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
extern "C"{
#include "assert.h"
}
#endif

#include <vector>

#include "MultilevelLinearOperator.h"
#include "PflotranJacobianLevelOperator.h"

extern "C" {
#include "petscmat.h"
}

namespace SAMRAI{

class PflotranJacobianMultilevelOperator: public MultilevelLinearOperator
{
public:

   PflotranJacobianMultilevelOperator(MultilevelOperatorParameters *parameters);

   ~PflotranJacobianMultilevelOperator();

   /**
    * Compute forward apply operation on a level, with the default arguments for a and b
    * the residual, r=f-Au, should be calculated, else b*f+a*A*u is calculated.
    *  \param ln
    *            level number for which apply() should be done
    * \param f_id
    *            array of integer id's for right hand side
    * \param u_id
    *            array of integer id's for solution variables
    * \param r_id
    *            array of integer id's for residual
    * \param f_idx
    *            Array of integer id's for components of right hand side.
    *            These must correspond to the integer id's supplied in f_id
    * \param u_idx
    *            Array of integer id's for components of solution variables.
    *            These must correspond to the integer id's supplied in u_id
    * \param r_idx
    *            Array of integer id's for components of residual.
    *            These must correspond to the integer id's supplied in r_id
    * \param a
    *            Double scaling factor for A*u, default of -1.0
    * \param b
    *            Double scaling factor for rhs, default of 1.0
    */
   void apply(const int ln,
              const int *f_id,
              const int *u_id, 
              const int *r_id,
              const int *f_idx=NULL,
              const int *u_idx=NULL,
              const int *r_idx=NULL,
              const double a = -1.0,
              const double b = 1.0);

   /**
    * Compute forward apply operation on a range of levels, with the default arguments for a and b
    * the residual, r=f-Au, should be calculated, else b*f+a*A*u is calculated.
    *  \param coarse_ln
    *            lower bound on levels for which apply() should be done
    *  \param fine_ln
    *            upper bound on levels for which apply() should be done
    * \param f_id
    *            array of integer id's for right hand side
    * \param u_id
    *            array of integer id's for solution variables
    * \param r_id
    *            array of integer id's for residual
    * \param f_idx
    *            Array of integer id's for components of right hand side.
    *            These must correspond to the integer id's supplied in f_id
    * \param u_idx
    *            Array of integer id's for components of solution variables.
    *            These must correspond to the integer id's supplied in u_id
    * \param r_idx
    *            Array of integer id's for components of residual.
    *            These must correspond to the integer id's supplied in r_id
    * \param a
    *            Double scaling factor for A*u, default of -1.0
    * \param b
    *            Double scaling factor for rhs, default of 1.0
    */
   void apply(const int coarse_ln,
              const int fine_ln,
              const int *f_id,
              const int *u_id, 
              const int *r_id,
              const int *f_idx=NULL,
              const int *u_idx=NULL,
              const int *r_idx=NULL,
              const double a = -1.0,
              const double b = 1.0);

   /**
   * Apply a boundary condition for a subset of variable, index pairs. Periodic boundaries
   * will be automatically detected and set.
   * \pre Initialization of the boundary condition types for the different variables
   *      should be done when the constructor is called.
   * \param var_id
   *         Variable descriptor index 
   * \param var_idx
   *         Component index
   * \param var_components
   *         We assume that the order of the descriptor indices supplied in var_id and
   *         the corresponding component indices supplied in var_idx need not correspond 
   *         with the order of the variables assumed within the discrete operator. For example
   *         the discrete equations may be ax+by=f, cx+dy=g. Within the operator, x is considered
   *         to be the first variable and y the second variable. However, var_id[0] may correspond
   *         y and var_id[1] to x. In which case var_components[0]=1 and var_components[1]=0
   *         var_components specifies the mapping from (var_id, var_idx) to the order assumed
   *         within the discrete operator.
   * \param number_of_variables
   *         For scalars the default is always true, for systems, in cases where it is desirable
   *         to call applyBoundaryCondition() for a subset of the variables, this specifies the
   *         number of variables for which we are calling applyBoundaryCondition().
   * \param reset_ghost_values
   *          boolean specifying whether interpolation from coarser levels and caching
   *          internally of the data is required before applying boundary conditions, default
   *          is false to minimize interpolations
   */
   void applyBoundaryCondition(const int ln,
                               const int *var_id,
                               const int *var_idx=NULL,
                               const int *var_components=NULL,
                               const int number_of_variables=-1,
                               const bool reset_ghost_values = false);
  
   /**
   * Returns the number of primitive variables for the discretization
   */
   const int getNumberOfVariables(void){ return d_ndof; }
   
   static PetscErrorCode wrapperMatMult(Mat mat,Vec x,Vec y);

   PetscErrorCode MatMult(Mat mat,Vec x,Vec y);

   static PetscErrorCode wrapperMatZeroEntries(Mat mat);

   PetscErrorCode MatZeroEntries(Mat mat);

   static PetscErrorCode wrapperMatSetValuesLocal(Mat mat,
                                                  PetscInt nrow,const PetscInt irow[],
                                                  PetscInt ncol,const PetscInt icol[],
                                                  const PetscScalar y[],InsertMode addv);
      
   PetscErrorCode MatSetValuesLocal(Mat mat,
                                    PetscInt nrow,const PetscInt irow[],
                                    PetscInt ncol,const PetscInt icol[],
                                    const PetscScalar y[],InsertMode addv);
      
   static PetscErrorCode wrapperMatSetValuesBlockedLocal(Mat mat,
                                                         PetscInt nrow,const PetscInt irow[],
                                                         PetscInt ncol,const PetscInt icol[],
                                                         const PetscScalar y[],InsertMode addv);
      
   PetscErrorCode MatSetValuesBlockedLocal(Mat mat,
                                           PetscInt nrow,const PetscInt irow[],
                                           PetscInt ncol,const PetscInt icol[],
                                           const PetscScalar y[],InsertMode addv);
      
   tbox::Pointer<hier::PatchHierarchy<NDIM> > getHierarchy(void){ return d_hierarchy; }
   
protected:

   void initializePetscMatInterface(void);
   /**
   * Allocates space for the internally used face centered flux
   */
   void initializeInternalVariableData(void);

   void getFromInput(tbox::Pointer<tbox::Database> db);

   /**
   * Compute fluxes.
   *  \param coarse_ln
   *            lower bound on levels for which setFlux() should be called
   *  \param fine_ln
   *            upper bound on levels for which setFlux() should be called
   * \param u_id
   *        descriptor index for cell centered variable to calculate flux from
   * \param u_idx
   *        component index for cell centered variable to calculate flux from
   */

   void setFlux(const int coarse_ln,
                const int fine_ln,
                const int *u_id,
                const int *u_idx=NULL);

   /**
   * Routine used to set up internal transfer schedules
   */
   void setupTransferSchedules(void);

   void initializeBoundaryConditionStrategy(tbox::Pointer<tbox::Database> &db);

private:

   PflotranJacobianMultilevelOperator();

   bool d_adjust_cf_coefficients;
   bool d_coarsen_diffusive_fluxes;
   bool d_schedules_initialized;
   bool d_use_cf_interpolant;
   bool d_variable_order_interpolation;
   bool d_reset_ghost_values;

   std::string d_face_coarsen_op_str;
   std::string d_face_refine_op_str;
   std::string d_cell_coarsen_op_str;
   std::string d_cell_refine_op_str;

   int d_flux_id;

   int d_bdry_types[2*NDIM];

   int d_ndof;

   int d_stencilSize;

   tbox::Pointer< pdat::FaceVariable<NDIM,double> > d_flux;

   RefinementBoundaryInterpolation::InterpolationScheme d_tangent_interp_scheme;
   RefinementBoundaryInterpolation::InterpolationScheme d_normal_interp_scheme;

   tbox::Array<PflotranJacobianLevelOperator*> d_level_operators;

   Mat *d_pMatrix;
};


}

#endif //included_PflotranJacobianMultilevelOperator
