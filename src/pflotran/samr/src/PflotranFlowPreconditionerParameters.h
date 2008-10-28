#ifndef included_PflotranFlowPreconditionerParameters
#define included_PflotranFlowPreconditionerParameters

#ifndef included_PreconditionerParameters
#include "PreconditionerParameters.h"
#endif

#ifndef included_Database
#include "tbox/Database.h"
#endif

#ifndef included_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef included_PatchHierarchy
#include "PatchHierarchy.h"
#endif

#ifndef included_ParameterBase
#include "ParameterBase.h"
#endif

#ifndef included_RefinementBoundaryInterpolation
#include "RefinementBoundaryInterpolation.h"
#endif

extern "C"{
#include "petsc.h" 
#include "petscpc.h" 
}

namespace SAMRAI{

class PflotranFlowPreconditionerParameters: public PreconditionerParameters
{
public:
   PflotranFlowPreconditionerParameters();

   PflotranFlowPreconditionerParameters( const tbox::Pointer<tbox::Database>& database);

   ~PflotranFlowPreconditionerParameters();

   /**
   *  Database object which needs to be initialized specific to the preconditioner.
   *  Documentation for parameters required by each preconditioner can be found in the
   *  documentation for the preconditioner.
   */
   tbox::Pointer<tbox::Database> d_db;

   /**
   * Pointer to the patch hierarchy object used by the preconditioner.
   */
   tbox::Pointer<hier::PatchHierarchy<NDIM> > d_hierarchy;

   /**
   * Pointer to the RefinementBoundaryInterpolation object used by the preconditioner to initialize and set coarse-fine boundary values.
   */
   RefinementBoundaryInterpolation *d_cf_interpolant;      // object storing refinement boundary geometry and ghost value info.
   
   PC *d_pc;  // pointer to PETSc preconditioner object 

protected:

private:

};
}
#endif
