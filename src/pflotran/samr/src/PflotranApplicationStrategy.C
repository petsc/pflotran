#include "PflotranApplicationStrategy.h"

namespace SAMRAI{
  
PflotranApplicationStrategy::PflotranApplicationStrategy()
{
}

PflotranApplicationStrategy::PflotranApplicationStrategy(PflotranApplicationParameters *params)
{
   initialize(params);
}

PflotranApplicationStrategy::~PflotranApplicationStrategy()
{
}

void
PflotranApplicationStrategy::initialize(PflotranApplicationParameters *params)
{
   d_hierarchy = params->d_hierarchy;
}

/**
* Return ComponentSelector of data to allocate on a new level.
*/
hier::ComponentSelector 
PflotranApplicationStrategy::getDataToAllocate()
{
   return d_AllocSelector;
}

/**
* Return ComponentSelector of data to time stamp on a new level.
*/
hier::ComponentSelector 
PflotranApplicationStrategy::getDataToTimeStamp()
{
   return d_TimestampSelector;
}

/** 
* Set initial conditions on new level.
*/
void 
PflotranApplicationStrategy::setInitialConditions( const double initial_time, 
                                                   tbox::Pointer< hier::PatchLevel<NDIM> > level )
{
   tbox::pout << "ERROR::PflotranApplicationStrategy::setInitialConditions not implemented as yet!!" << std::endl;
}

/**
* Set initial conditions in specified vector.
*/
void 
PflotranApplicationStrategy::setInitialConditions( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > ic )
{
   tbox::pout << "ERROR::PflotranApplicationStrategy::setInitialConditions not implemented as yet!!" << std::endl;
}

/** 
* Set data values on new level.
*/
void 
PflotranApplicationStrategy::setValuesOnNewLevel( tbox::Pointer< hier::PatchLevel<NDIM> > level )
{
   tbox::pout << "ERROR::PflotranApplicationStrategy::setValuesOnNewLevel not implemented as yet!!" << std::endl;
}

/**
* Return a list of variables used by this application.
*/
tbox::Array< tbox::Pointer< hier::Variable<NDIM> > > 
PflotranApplicationStrategy::getVariables()
{
   return(d_variable_list);
}

/**
* Create the RefineSchedules needed to transfer data to a new level.
*/
tbox::Array< tbox::Pointer< xfer::RefineSchedule<NDIM> > > 
PflotranApplicationStrategy::setupRegridRefineSchedules( const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
                                                         const int level_number,  
                                                         const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level )
{

   return(d_regrid_refine_scheds);
}

/**
* Mark locations where additional refinement is desired.
*/
void 
PflotranApplicationStrategy::tagCells( const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
                                       const int level_number,
                                       const double error_data_time,
                                       const bool initial_time,
                                       const int tag_index,
                                       const bool uses_richardson_extrapolation_too )
{
   tbox::pout << "ERROR::PflotranApplicationStrategy::tagCells not implemented as yet!!" << std::endl;
}

/**
* Update data structures that change when the grid hierarchy changes.
*/
void 
PflotranApplicationStrategy::resetHierarchyConfiguration( const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
                                                          const int coarsest_level,
                                                          const int finest_level )
{
   tbox::pout << "ERROR::PflotranApplicationStrategy::resetHierarchyConfiguration not implemented as yet!!" << std::endl;
}

/**
* Evaluate IVP forcing term.
*/
int 
PflotranApplicationStrategy::evaluateFunction( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  x,
                                               tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  f )
{
   tbox::pout << "ERROR:: PflotranApplicationStrategy::evaluateFunction not implemented as yet!!" << std::endl;
   return(0);
}

/**
* Print identifying std::string.
*/
void 
PflotranApplicationStrategy::printObjectName( std::ostream& os )
{
   os << d_object_name;
}

/**
* Allocate data for time integrator
*/ 
void 
PflotranApplicationStrategy::allocateVectorData(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > x,
                                                double time, bool flag)
{
   tbox::pout << "ERROR:: PflotranApplicationStrategy::allocateVectorData not implemented as yet!!" << std::endl;

   x->allocateVectorData();
} 

}
