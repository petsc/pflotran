#ifndef included_PflotranJacobianLevelOperator
#define included_PflotranJacobianLevelOperator

#ifndef included_RefineSchedule
#include "RefineSchedule.h"
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
extern "C"{
#include "assert.h"
}
#endif

#include "FaceVariable.h"
#include "CCellVariable.h"
#include "LevelOperatorParameters.h"
#include "LevelLinearOperator.h"
#include <vector>
#include "RefinementBoundaryInterpolation.h"

extern "C" {
#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
}

namespace SAMRAI{

class PflotranJacobianLevelOperator: public LevelLinearOperator
{
public:

   /**
   * constructor
   * \param parameters 
   *        A LevelOperatorParameters class used to provide arguments
   *        to initialize the PflotranJacobianLevelOperator class
   */
   PflotranJacobianLevelOperator(LevelOperatorParameters *parameters);

   ~PflotranJacobianLevelOperator();

   /**
    * Compute forward apply operation, with the default arguments for a and b
    * the residual, r=f-Au, should be calculated, else b*f+a*A*u is calculated.
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
    */
   void apply(const int *f_id,
              const int *u_id, 
              const int *r_id,
              const int *f_idx=NULL,
              const int *u_idx=NULL,
              const int *r_idx=NULL,
              const double a = -1.0,
              const double b = 1.0);
   
   /**
   * Apply a boundary condition for a specific variable and index. Periodic boundaries
   * will be automatically detected and set.
   * \pre Initialization of the boundary condition types for the different variables
   *      should be done when the constructor is called.
   * \param var_id
   *         variable descriptor index 
   * \param var_idx
   *         component index
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
   */
   void applyBoundaryCondition(const int *var_id,
                               const int *var_idx=NULL,
                               const int *var_components=NULL,
                               const int number_of_variables=-1);
   
   /**
   * Return the number of non zero stencil elements for block i,j of the stencil.
   * Block i,j of the stencil contains non zero connections between the i-th and j-th variables.
   * \param i
   *        i is the first variable
   * \param j
   *        index of second variable
   */
   const int getStencilSize(const int i=0,
                            const int j=0,
                            const int k=0);
   
   /**
   * Return whether the stencil for block (i,j,k) is a constant stencil or variable coefficient
   * Block i,j of the stencil contains non zero connections between the i-th and j-th variables.
   * Returns zero if constant coefficient else returns one.
   * \param i
   *        i is the first variable
   * \param j
   *        index of second variable
   */
   const int getStencilType(const int i=0,
                            const int j=0,
                            const int k=0);
   
   /**
   * Returns a pointer to the stencil block data for block i,j of the stencil.
   * Block i,j of the stencil contains non zero connections between the i-th and j-th variables.
   * Will return a NULL pointer if all elements of that stencil block are zero.
   * \param p
   *        patch number for which to get stencil data
   * \param i
   *        index of first variable
   * \param j
   *        index of second variable
   */
   tbox::Pointer< hier::PatchData<NDIM > > getStencilBlock(const int p, 
                                                           const int i=0, 
                                                           const int j=0,
                                                           const int k=0);


   /**
   * Returns integer array of offsets for non zero elements of the i,j stencil block.
   * Block i,j of the stencil contains non zero connections between the i-th and j-th variables.
   * \param i
   *        index of first variable
   * \param j
   *        index of second variable   
   */
   std::vector<int> getStencilOffsets(const int i=0, 
                                      const int j=0, 
                                      const int k=0);

   /**
   * Returns the number of primitive variables for the discretization
   */
   const int getNumberOfVariables(void){ return d_ndof; }

   int getStencilID(void){return d_stencil_id;}

   void MatSetValuesLocal(int patchNumber,
                          PetscInt nrow,const PetscInt irow[],
                          PetscInt ncol,const PetscInt icol[],
                          const PetscScalar y[],InsertMode addv);

   void MatSetValuesBlockedLocal(int patchNumber,
                                 PetscInt nrow,const PetscInt irow[],
                                 PetscInt ncol,const PetscInt icol[],
                                 const PetscScalar y[],InsertMode addv);

   int MatMult(Vec x, Vec y);

protected:

   void getFromInput(const tbox::Pointer<tbox::Database> &db);

   void initializeInternalVariableData();
   
private:

   PflotranJacobianLevelOperator();

   int d_ndof;
   
   bool d_adjust_cf_coefficients;
   bool d_interpolate_ghost_values;
   bool d_sibling_fill_cached;
   bool d_stencil_initialized;
   bool d_coefficients_changed;
   bool d_variable_order_interpolation;
   bool d_reset_ghost_cells;

   int d_extrapolation_order;           // extrapolation order to use for ghost cells on physical boundaries
   int d_bdry_types[2*NDIM];
   int d_stencil_size;
   int d_stencil_id;

   tbox::Pointer<xfer::RefineSchedule<NDIM> > d_sibling_fill_schedule;

   tbox::Pointer<pdat::FaceVariable<NDIM,double> > d_flux;
   tbox::Pointer<pdat::CCellVariable<NDIM,double> > d_stencil;

   int d_flux_id;

   RefinementBoundaryInterpolation::InterpolationScheme d_tangent_interp_scheme;
   RefinementBoundaryInterpolation::InterpolationScheme d_normal_interp_scheme;
   
};

}
#endif //  included_PflotranJacobianLevelOperator
