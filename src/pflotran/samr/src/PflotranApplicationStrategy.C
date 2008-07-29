#include "PflotranApplicationStrategy.h"
#include "tbox/TimerManager.h"

namespace SAMRAI{
  
PflotranApplicationStrategy::PflotranApplicationStrategy()
{
}

PflotranApplicationStrategy::PflotranApplicationStrategy(PflotranApplicationParameters *params)
{
   initialize(params);

   d_read_regrid_boxes                       = false;
   d_error_checkpoint                        = false;
   d_is_after_regrid                         = false;
   d_use_variable_order_interpolation        = false;
   d_coarsen_fluxes                          = true;

   d_soln_refine_op =  d_grid_geometry->lookupRefineOperator(d_solution,
                                                             "CONSTANT_REFINE");

   d_GlobalToLocalRefineSchedule.resizeArray(d_hierarchy->getNumberOfLevels());

}

PflotranApplicationStrategy::~PflotranApplicationStrategy()
{
}

void
PflotranApplicationStrategy::initialize(PflotranApplicationParameters *params)
{
   d_hierarchy = params->d_hierarchy;
   d_grid_geometry = d_hierarchy->getGridGeometry();
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

void
PflotranApplicationStrategy::interpolateLocalToLocalVector(int ndof,
                                                           tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  srcVec,
                                                           tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  destVec,
                                                           int ierr)
{

}

void
PflotranApplicationStrategy::interpolateGlobalToLocalVector(int ndof,
                                                            tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  globalVec,
                                                            tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  localVec,
                                                            int ierr)
{
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy =  d_hierarchy;

    static tbox::Pointer<tbox::Timer> t_interpolate_variable = tbox::TimerManager::getManager()->getTimer("PFlotran::PflotranApplicationStrategy::interpolateVector");

    t_interpolate_variable->start();

    int src_id = globalVec->getComponentDescriptorIndex(0);
    int dest_id = localVec->getComponentDescriptorIndex(0);   

    xfer::RefineAlgorithm<NDIM> ghost_cell_fill;

    ghost_cell_fill.registerRefine(dest_id,
				   src_id,
				   dest_id,
				   d_soln_refine_op);

    // now refine and set physical boundaries also
    for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ln++ ) 
    {

      tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

      d_refine_patch_strategy->setDataID(dest_id);

      if(ghost_cell_fill.checkConsistency(d_GlobalToLocalRefineSchedule[ln]))
      {
         ghost_cell_fill.resetSchedule(d_GlobalToLocalRefineSchedule[ln]);
      }
      else
      {
         d_GlobalToLocalRefineSchedule[ln] = ghost_cell_fill.createSchedule(
            level, 
            ln-1,
            hierarchy,
            d_refine_patch_strategy);
      }
      
      d_GlobalToLocalRefineSchedule[ln]->fillData(d_current_time);            

      if(ln>0)
      {	
          d_cf_interpolant->setVariableOrderInterpolation(d_use_variable_order_interpolation);
          d_cf_interpolant->setPhysicalCornerRefGhostValues(ln,
                                                            dest_id,
                                                            0);
          d_cf_interpolant->setGhostCellData(ln, dest_id);

	  for ( int idx=0; idx<d_number_solution_components; idx++)
	  {
             d_cf_interpolant->interpolateGhostValues(ln,
                                                      d_nl_tangential_interp_scheme,
                                                      d_nl_normal_interp_scheme,
                                                      dest_id,
                                                      idx);     
          }
      }
    }

    // should add code to coarsen variables

    t_interpolate_variable->stop();
}

void
PflotranApplicationStrategy::setRefinementBoundaryInterpolant(RefinementBoundaryInterpolation *cf_interpolant)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(cf_interpolant!=NULL);
#endif
   d_cf_interpolant = cf_interpolant;

   if(d_Jacobian.get()!=NULL)
   {
      d_Jacobian->setRefinementBoundaryInterpolant(cf_interpolant);
   }

}

}
