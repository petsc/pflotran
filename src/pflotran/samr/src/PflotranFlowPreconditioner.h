#ifndef included_PflotranFlowPreconditioner
#define included_PflotranFlowPreconditioner


#ifndef included_PreconditionerStrategy
#include "PreconditionerStrategy.h"
#endif

#ifndef included_RefinementBoundaryInterpolation
#include "RefinementBoundaryInterpolation.h"
#endif

#ifndef included_PflotranFlowPreconditionerParameters
#include "PflotranFlowPreconditionerParameters.h"
#endif

#ifndef included_FaceVariable
#include "FaceVariable.h"
#endif

#include "MultilevelSolverFactory.h"
#include "MultilevelSolver.h"
#include "MultilevelLinearOperatorFactory.h"
#include "PflotranJacobianMultilevelOperator.h"
#include "CoarsenSchedule.h"
#include "CoarsenAlgorithm.h"
#include "tbox/List.h"
#include <vector>

typedef SAMRAI::tbox::List< SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > crsList;

namespace SAMRAI{

class PflotranFlowPreconditioner: public PreconditionerStrategy
{
public:
   PflotranFlowPreconditioner(PflotranFlowPreconditionerParameters *parameters);

   ~PflotranFlowPreconditioner();
   
   int setupPreconditioner(  PreconditionerParameters* parameters );
   
   int applyPreconditioner( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > r,
                            tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > z );

   /**
    * Functions to read data from input and restart databases. If the
    * boolean flag is true, all data members are read from restart.
    * They can later be overwritten from values in the input file.
    * When the flag is false, all data values are set from those given
    * in input.
    *
    * If assertion checking is enabled, an unrecoverable exception
    * results if the database pointer is null.
   */
   void getFromInput( tbox::Pointer<tbox::Database> &db,
                      bool is_from_restart = false);

   void setRefinementBoundaryInterpolant(RefinementBoundaryInterpolation *cf_interpolant);

   void setOperator(PflotranJacobianMultilevelOperator *op){d_pc_operator = op;}

   static PetscErrorCode wrapperSetupPreconditioner(void *ptr);

   static PetscErrorCode wrapperApplyPreconditioner(void *ptr,Vec xin,Vec xout);

protected:

private:

   // private constuctor to prevent it being called
   PflotranFlowPreconditioner();

   void preprocessPCApply( int r_id );

   void postprocessPCApply( int z_id );

   void interpolateVariable( const int src_id, 
                             const int dest_id, 
                             RefinementBoundaryInterpolation::InterpolationScheme tangential_interp_scheme,
                             RefinementBoundaryInterpolation::InterpolationScheme normal_interp_scheme );

   void initializeSolvers(tbox::Pointer<tbox::Database> &db);

   void initializePetscInterface(void);

   void coarsenVariable( const int var_id,
                         std::string coarsen_op_str);
   /*
    * Retrieve list of potentially usable sibling fill schedules for a level.
    */
   crsList* getCoarsenSchedules(int ln) { return( &d_coarsen_fill_schedules[ln] ); }

   /*
    * Retrieve list of potentially usable sibling fill schedules for a level.
    */
    tbox::Pointer<xfer::CoarsenSchedule<NDIM> > getCoarsenSchedule(int ln, xfer::CoarsenAlgorithm<NDIM> &);

   bool  d_preconditioner_print_flag;

   bool d_pc_solver_op_registered;

   /*
    * Cached list of fill schedules, to avoid creating one for every
    * sibling fill.
    */
   std::vector< crsList > d_coarsen_fill_schedules;

   tbox::Pointer<hier::PatchHierarchy<NDIM> > d_hierarchy;

   RefinementBoundaryInterpolation *d_cf_interpolant;

   MultilevelSolverFactory *d_PCSolverFactory;

   MultilevelSolver *d_pc_solver;

   PflotranJacobianMultilevelOperator *d_pc_operator;

   PC *d_pc;
   
};
}
#endif
