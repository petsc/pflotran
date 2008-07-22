#ifndef included_PflotranApplicationStrategy
#define included_PflotranApplicationStrategy


#include <string>
#include <vector>

#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "CellVariable.h"
#include "FaceVariable.h"
#include "PatchHierarchy.h"
#include "VariableContext.h"
#include "FaceData.h"

#include "ComponentSelector.h"
#include "ApplicationStrategy.h"
#include "PflotranApplicationParameters.h"

namespace SAMRAI{

class PflotranApplicationStrategy: public ApplicationStrategy
{
public:
   PflotranApplicationStrategy(PflotranApplicationParameters *params);
   ~PflotranApplicationStrategy();

   /**
    * Return ComponentSelector of data to allocate on a new level.
    */
   hier::ComponentSelector getDataToAllocate();

   /**
    * Return ComponentSelector of data to time stamp on a new level.
    */
   hier::ComponentSelector getDataToTimeStamp();

   /** 
    * Set initial conditions on new level.
    */
   void setInitialConditions( const double initial_time, tbox::Pointer< hier::PatchLevel<NDIM> > level );

   /**
    * Set initial conditions in specified vector.
    */
   void setInitialConditions( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > ic );

   /** 
    * Set data values on new level.
    */
   void setValuesOnNewLevel( tbox::Pointer< hier::PatchLevel<NDIM> > level );

   /**
    * Return a list of variables used by this application.
    */
   tbox::Array< tbox::Pointer< hier::Variable<NDIM> > > getVariables();

   /**
    * Create the RefineSchedules needed to transfer data to a new level.
    */
   tbox::Array< tbox::Pointer< xfer::RefineSchedule<NDIM> > > setupRegridRefineSchedules( 
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int level_number,  
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level );

   /**
    * Mark locations where additional refinement is desired.
    */
   void tagCells( const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
                          const int level_number,
                          const double error_data_time,
                          const bool initial_time,
                          const int tag_index,
                          const bool uses_richardson_extrapolation_too );

   /**
    * Update data structures that change when the grid hierarchy changes.
    */
   void resetHierarchyConfiguration( const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
                                             const int coarsest_level,
                                             const int finest_level );

   /**
    * Evaluate IVP forcing term.
    */
   int evaluateFunction( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  x,
                                 tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  f );

   /**
    * Print identifying std::string.
    */
   void printObjectName( std::ostream& os );

   /**
    * Allocate data for time integrator
    */ 
   void allocateVectorData(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > x,
                                   double time, bool flag=true ); 


   // Hierarchy
   tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;


protected:

private:

   PflotranApplicationStrategy();

   void initialize(PflotranApplicationParameters *params);
   
   // Name of application
   std::string d_object_name;

   int d_number_solution_components;
   tbox::Array< tbox::Pointer< hier::Variable<NDIM> > > d_variable_list;

   tbox::Array< tbox::Pointer< xfer::RefineSchedule<NDIM> > > d_regrid_refine_scheds;

   tbox::Pointer<hier::VariableContext> d_application_ctx;

   hier::ComponentSelector d_AllocSelector;
   hier::ComponentSelector d_TimestampSelector;
};

}

#endif // included_PflotranApplicationStrategy
