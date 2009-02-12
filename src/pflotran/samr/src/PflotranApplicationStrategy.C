#include "PflotranApplicationStrategy.h"
#include "tbox/TimerManager.h"
#include "CCellDataFactory.h"
#include "CSideDataFactory.h"
#include "tbox/RestartManager.h"
#include "CellVariable.h"
#include "CCellVariable.h"
#include "SideVariable.h"
#include "CSideVariable.h"
#include "SAMRAIVectorReal.h"
#include "PETSc_SAMRAIVectorReal.h"
#include "HierarchyCCellDataOpsReal.h"
#include "HierarchyCSideDataOpsReal.h"


namespace SAMRAI{
  
int PflotranApplicationStrategy::d_vec_instance_id=0;

PflotranApplicationStrategy::PflotranApplicationStrategy()
{
}

PflotranApplicationStrategy::PflotranApplicationStrategy(PflotranApplicationParameters *params)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(params!=NULL);
#endif

   d_read_regrid_boxes                       = false;
   d_is_after_regrid                         = false;
   d_use_variable_order_interpolation        = false;
   d_coarsen_fluxes                          = true;
   d_error_checkpoint                        = false;

   d_current_time                            = 0.0;

   d_object_name                             = "PflotranApplicationStrategy";
   
   d_number_solution_components              = 0;

   d_cf_interpolant                          = NULL;

   d_nl_tangential_interp_scheme             = RefinementBoundaryInterpolation::linear;

   d_nl_normal_interp_scheme                 = RefinementBoundaryInterpolation::linear;

   d_refine_patch_strategy                   = new BoundaryConditionStrategy(-1);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_refine_patch_strategy!=NULL);
#endif

   d_variable_list.resizeArray(0);

   d_soln_refine_op.setNull();

   d_soln_coarsen_op.setNull();

   d_flux_coarsen_op.setNull();

   d_regrid_refine_scheds.resizeArray(0);

   initialize(params);

   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();
   
   if(!variable_db->checkVariableExists("pflotranSolution"))
   {
      d_solution                                = new pdat::CCellVariable<NDIM,double>("pflotranSolution", d_number_solution_components);
   }
   else
   {
      d_solution = variable_db->getVariable("pflotranSolution");
   }

   if(!variable_db->checkVariableExists("pflotranWeight"))
   {
      d_pflotran_weight = new pdat::CCellVariable<NDIM,double>("pflotranWeight", d_number_solution_components);
   }
   else
   {
      d_pflotran_weight = variable_db->getVariable("pflotranWeight");
   }

   const SAMRAI::tbox::Pointer< SAMRAI::hier::VariableContext > pflotran_cxt = variable_db->getContext("PFLOTRAN");

   d_pflotran_weight_id = variable_db->registerVariableAndContext(d_pflotran_weight,
                                                                  pflotran_cxt,
                                                                  hier::IntVector<NDIM>(0));

   for(int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++)
   {
      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      level->allocatePatchData(d_pflotran_weight_id);
   }
   
   AMRUtilities::setVectorWeights(d_hierarchy, d_pflotran_weight_id);

   d_application_ctx = variable_db->getContext(d_object_name);

   d_soln_refine_op =  d_grid_geometry->lookupRefineOperator(d_solution,
                                                             "CCELL_CONSTANT_REFINE");

   d_soln_coarsen_op = d_grid_geometry->lookupCoarsenOperator(d_solution,
                                                              "CONSERVATIVE_COARSEN");
   
   d_GlobalToLocalRefineSchedule.resizeArray(d_hierarchy->getNumberOfLevels());
   d_LocalToLocalRefineSchedule.resizeArray(d_hierarchy->getNumberOfLevels());
   d_CoarsenSchedule.resizeArray(d_hierarchy->getNumberOfLevels());
   d_FluxCoarsenSchedule.resizeArray(d_hierarchy->getNumberOfLevels());
   
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_number_solution_components>=1);
#endif

   for(int ln=0;ln<d_hierarchy->getNumberOfLevels(); ln++)
   {
      d_GlobalToLocalRefineSchedule[ln].resizeArray(d_number_solution_components);
      d_LocalToLocalRefineSchedule[ln].resizeArray(d_number_solution_components);
      d_CoarsenSchedule[ln].resizeArray(d_number_solution_components);
      d_FluxCoarsenSchedule[ln].resizeArray(d_number_solution_components);

      for(int i=0;i<d_number_solution_components; i++)
      {
         d_GlobalToLocalRefineSchedule[ln][i].setNull();
         d_LocalToLocalRefineSchedule[ln][i].setNull();
         d_CoarsenSchedule[ln][i].setNull();
         d_FluxCoarsenSchedule[ln][i].setNull();
      }
   }

   d_ccell_math_op = new math::HierarchyCCellDataOpsReal< NDIM, double >(d_hierarchy,
                                                                         0, d_hierarchy->getFinestLevelNumber());

   d_cside_math_op = new math::HierarchyCSideDataOpsReal< NDIM, double >(d_hierarchy,
                                                                         0, d_hierarchy->getFinestLevelNumber());

   d_visit_writer = new appu::VisItDataWriter<NDIM>("rmhd visit writer", d_viz_directory);
   
}

PflotranApplicationStrategy::~PflotranApplicationStrategy()
{
}

void
PflotranApplicationStrategy::initialize(PflotranApplicationParameters *params)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(params!=NULL);
#endif

   d_hierarchy       = params->d_hierarchy;
   d_grid_geometry   = d_hierarchy->getGridGeometry();
   d_application_db  = params->d_database;   

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();

   PflotranApplicationStrategy::getFromInput(d_application_db, is_from_restart);

}

void
PflotranApplicationStrategy::getFromInput(tbox::Pointer<tbox::Database> db,
                                          bool is_from_restart)
{
   
   if (db->keyExists("nl_tangential_coarse_fine_scheme")) {
      d_nl_tangential_interp_scheme=RefinementBoundaryInterpolation::lookupInterpolationScheme(db->getString("nl_tangential_coarse_fine_scheme"));
   }
   
   if (db->keyExists("nl_normal_coarse_fine_scheme")) {
      d_nl_normal_interp_scheme = RefinementBoundaryInterpolation::lookupInterpolationScheme(db->getString("nl_normal_coarse_fine_scheme"));
   }
   
   if (db->keyExists("number_solution_components")) {
      d_number_solution_components = db->getInteger("number_solution_components");
   }
   
   d_viz_directory = db->getStringWithDefault("viz_directory", "viz");
   
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

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!x.isNull());
#endif

   x->allocateVectorData();
} 

void
PflotranApplicationStrategy::interpolateLocalToLocalVector(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  srcVec,
                                                           tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  dstVec,
                                                           int ierr)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!srcVec.isNull());
   assert(!dstVec.isNull());
#endif

    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy =  d_hierarchy;

    static tbox::Pointer<tbox::Timer> t_interpolate_variable = tbox::TimerManager::getManager()->getTimer("PFlotran::PflotranApplicationStrategy::interpolateVector");

    t_interpolate_variable->start();

    int src_id = srcVec->getComponentDescriptorIndex(0);
    int dest_id = dstVec->getComponentDescriptorIndex(0);   

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(src_id>=0);
    assert(dest_id>=0);
#endif

    tbox::Pointer< hier::Variable< NDIM > > srcVar = srcVec->getComponentVariable(0);

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!srcVar.isNull());
#endif

    tbox::Pointer< pdat::CCellDataFactory< NDIM, double > > srcFactory = srcVar->getPatchDataFactory(); 

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!srcFactory.isNull());
#endif

    int srcDOF = srcFactory->getDefaultDepth();

    tbox::Pointer< hier::Variable< NDIM > > dstVar = dstVec->getComponentVariable(0);

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!dstVar.isNull());
#endif

    tbox::Pointer< pdat::CCellDataFactory< NDIM, double > > dstFactory = dstVar->getPatchDataFactory(); 

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!dstFactory.isNull());
#endif

    int dstDOF = dstFactory->getDefaultDepth();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(srcDOF>=1);
    assert(srcDOF=dstDOF);
#endif

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

      if((!d_LocalToLocalRefineSchedule[ln][srcDOF-1].isNull()) && ghost_cell_fill.checkConsistency(d_LocalToLocalRefineSchedule[ln][srcDOF-1]))
      {
         ghost_cell_fill.resetSchedule(d_LocalToLocalRefineSchedule[ln][srcDOF-1]);
      }
      else
      {
         d_LocalToLocalRefineSchedule[ln][srcDOF-1] = ghost_cell_fill.createSchedule(
            level, 
            ln-1,
            hierarchy,
            d_refine_patch_strategy);

#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!d_LocalToLocalRefineSchedule[ln][srcDOF-1].isNull());
#endif

      }
      
      d_LocalToLocalRefineSchedule[ln][srcDOF-1]->fillData(d_current_time);            

      if(ln>0)
      {	
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(d_cf_interpolant!=NULL);
#endif
         
         d_cf_interpolant->setVariableOrderInterpolation(d_use_variable_order_interpolation);
         d_cf_interpolant->setPhysicalCornerRefGhostValues(ln,
                                                           dest_id,
                                                           0);
         d_cf_interpolant->setGhostCellData(ln, dest_id);
         
         for ( int idx=0; idx<srcDOF; idx++)
         {
            d_cf_interpolant->interpolateGhostValues(ln,
                                                     d_nl_tangential_interp_scheme,
                                                     d_nl_normal_interp_scheme,
                                                     dest_id,
                                                     idx);     
         }
      }
    }
#if 1
    // should add code to coarsen variables
    xfer::CoarsenAlgorithm<NDIM> cell_coarsen;
    cell_coarsen.registerCoarsen(dest_id, dest_id, d_soln_coarsen_op);

    for (int ln = hierarchy->getNumberOfLevels()-2; ln>=0; ln-- ) 
    {
      tbox::Pointer<hier::PatchLevel<NDIM> > clevel = hierarchy->getPatchLevel(ln);
      tbox::Pointer<hier::PatchLevel<NDIM> > flevel = hierarchy->getPatchLevel(ln+1);

      for ( int i=0; i<srcDOF; i++)
      {
         if((!d_CoarsenSchedule[ln][i].isNull()) && cell_coarsen.checkConsistency(d_CoarsenSchedule[ln][i]))
         {
            cell_coarsen.resetSchedule(d_CoarsenSchedule[ln][i]);
         }
         else
         {
            d_CoarsenSchedule[ln][i] = cell_coarsen.createSchedule(clevel, 
                                                                   flevel);
         }
         
         d_CoarsenSchedule[ln][i]->coarsenData();            
      }
    }
#endif

    t_interpolate_variable->stop();
}

void
PflotranApplicationStrategy::interpolateGlobalToLocalVector(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  globalVec,
                                                            tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  localVec,
                                                            int ierr)
{
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy =  d_hierarchy;

    static tbox::Pointer<tbox::Timer> t_interpolate_variable = tbox::TimerManager::getManager()->getTimer("PFlotran::PflotranApplicationStrategy::interpolateVector");

    t_interpolate_variable->start();

    int src_id = globalVec->getComponentDescriptorIndex(0);
    int dest_id = localVec->getComponentDescriptorIndex(0);   

    tbox::Pointer< hier::Variable< NDIM > > localVar = localVec->getComponentVariable(0);
    tbox::Pointer< pdat::CCellDataFactory< NDIM, double > > localFactory = localVar->getPatchDataFactory(); 
    int localDOF = localFactory->getDefaultDepth();

    tbox::Pointer< hier::Variable< NDIM > > globalVar = globalVec->getComponentVariable(0);
    tbox::Pointer< pdat::CCellDataFactory< NDIM, double > > globalFactory = globalVar->getPatchDataFactory(); 
    int globalDOF = globalFactory->getDefaultDepth();

    assert(localDOF=globalDOF);

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

      if((!d_GlobalToLocalRefineSchedule[ln][globalDOF-1].isNull()) && ghost_cell_fill.checkConsistency(d_GlobalToLocalRefineSchedule[ln][globalDOF-1]))
      {
         ghost_cell_fill.resetSchedule(d_GlobalToLocalRefineSchedule[ln][globalDOF-1]);
      }
      else
      {
         d_GlobalToLocalRefineSchedule[ln][globalDOF-1] = ghost_cell_fill.createSchedule(
            level, 
            ln-1,
            hierarchy,
            d_refine_patch_strategy);
      }
      
      d_GlobalToLocalRefineSchedule[ln][globalDOF-1]->fillData(d_current_time);            

      if(ln>0)
      {	
          d_cf_interpolant->setVariableOrderInterpolation(d_use_variable_order_interpolation);
          d_cf_interpolant->setPhysicalCornerRefGhostValues(ln,
                                                            dest_id,
                                                            0);
          d_cf_interpolant->setGhostCellData(ln, dest_id);

	  for ( int idx=0; idx<globalDOF; idx++)
	  {
             d_cf_interpolant->interpolateGhostValues(ln,
                                                      d_nl_tangential_interp_scheme,
                                                      d_nl_normal_interp_scheme,
                                                      dest_id,
                                                      idx);     
          }
      }

#if 0
      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) 
      {
         tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
         
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!patch.isNull());
#endif


         tbox::Pointer< pdat::CCellData<NDIM, double > > src = patch->getPatchData(src_id);
         tbox::Pointer< pdat::CCellData<NDIM, double > > dst = patch->getPatchData(dest_id);
         tbox::pout << "debug code" << std::endl;
         
      }
#endif
  
    }

    // should add code to coarsen variables
    xfer::CoarsenAlgorithm<NDIM> cell_coarsen;
    cell_coarsen.registerCoarsen(dest_id, dest_id, d_soln_coarsen_op);

    for (int ln = hierarchy->getNumberOfLevels()-2; ln>=0; ln-- ) 
    {
      tbox::Pointer<hier::PatchLevel<NDIM> > clevel = hierarchy->getPatchLevel(ln);
      tbox::Pointer<hier::PatchLevel<NDIM> > flevel = hierarchy->getPatchLevel(ln+1);

      for ( int i=0; i<globalDOF; i++)
      {
         if((!d_CoarsenSchedule[ln][i].isNull()) && cell_coarsen.checkConsistency(d_CoarsenSchedule[ln][i]))
         {
            cell_coarsen.resetSchedule(d_CoarsenSchedule[ln][i]);
         }
         else
         {
            d_CoarsenSchedule[ln][i] = cell_coarsen.createSchedule(clevel, 
                                                                   flevel);
         }
         
         d_CoarsenSchedule[ln][i]->coarsenData();            
      }
    }

    t_interpolate_variable->stop();
}

void
PflotranApplicationStrategy::coarsenFaceFluxes(tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > fluxVec, 
                                               int ierr)
{

    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy =  d_hierarchy;

    static tbox::Pointer<tbox::Timer> t_coarsen_flux_variable = tbox::TimerManager::getManager()->getTimer("PFlotran::PflotranApplicationStrategy::coarsenFaceFluxes");

    t_coarsen_flux_variable->start();

    int flux_id = fluxVec->getComponentDescriptorIndex(0);

    tbox::Pointer< hier::Variable< NDIM > > fluxVar = fluxVec->getComponentVariable(0);
    tbox::Pointer< pdat::CSideDataFactory< NDIM, double > > sideFactory = fluxVar->getPatchDataFactory(); 
    int ndof = sideFactory->getDefaultDepth();

    if(d_flux_coarsen_op.isNull())
    {
       d_flux_coarsen_op = d_grid_geometry->lookupRefineOperator(fluxVar,
                                                                 "CONSERVATIVE_COARSEN");
    }

    // should add code to coarsen variables
    xfer::CoarsenAlgorithm<NDIM> flux_coarsen;
    flux_coarsen.registerCoarsen(flux_id, flux_id, d_flux_coarsen_op);

    for (int ln = hierarchy->getNumberOfLevels()-2; ln>=0; ln-- ) 
    {
      tbox::Pointer<hier::PatchLevel<NDIM> > clevel = hierarchy->getPatchLevel(ln);
      tbox::Pointer<hier::PatchLevel<NDIM> > flevel = hierarchy->getPatchLevel(ln+1);

      for ( int i=0; i<ndof; i++)
      {
         if((!d_FluxCoarsenSchedule[ln][i].isNull()) && flux_coarsen.checkConsistency(d_FluxCoarsenSchedule[ln][i]))
         {
            flux_coarsen.resetSchedule(d_FluxCoarsenSchedule[ln][i]);
         }
         else
         {
            d_FluxCoarsenSchedule[ln][i] = flux_coarsen.createSchedule(clevel, 
                                                                       flevel);
         }
         
         d_FluxCoarsenSchedule[ln][i]->coarsenData();            
      }
    }

    t_coarsen_flux_variable->stop();
}

void
PflotranApplicationStrategy::setRefinementBoundaryInterpolant(RefinementBoundaryInterpolation *cf_interpolant)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(cf_interpolant!=NULL);
#endif
   d_cf_interpolant = cf_interpolant;

   if(d_FlowJacobian.get()!=NULL)
   {
      d_FlowJacobian->setRefinementBoundaryInterpolant(cf_interpolant);
   }

   if(d_TransportJacobian.get()!=NULL)
   {
      d_TransportJacobian->setRefinementBoundaryInterpolant(cf_interpolant);
   }

   if(d_FlowPreconditioner.get()!=NULL)
   {
      d_FlowPreconditioner->setRefinementBoundaryInterpolant(cf_interpolant);
   }

   if(d_TransportPreconditioner.get()!=NULL)
   {
      d_TransportPreconditioner->setRefinementBoundaryInterpolant(cf_interpolant);
   }

}

void
PflotranApplicationStrategy::createVector(int &dof, 
                                          int &centering,
                                          bool &use_ghost, 
                                          bool &use_components, 
                                          Vec *vec)
{
   std::ostringstream ibuffer;
   ibuffer<<(long)PflotranApplicationStrategy::d_vec_instance_id;
   std::string object_str=ibuffer.str();

   std::string dataName("PFLOTRAN_variable_");
   dataName+=object_str;

   const int nlevels = d_hierarchy->getNumberOfLevels();
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > samrai_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(dataName,
                                                                                                                                     d_hierarchy,
                                                                                                                                     0, nlevels-1);
   
   PflotranApplicationStrategy::d_vec_instance_id++;

   int pflotran_var_id;
   
   SAMRAI::hier::VariableDatabase<NDIM>* variable_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
   
   const SAMRAI::tbox::Pointer< SAMRAI::hier::VariableContext > pflotran_cxt = variable_db->getContext("PFLOTRAN");
   
   SAMRAI::hier::IntVector<NDIM> nghosts;

   if(use_ghost)
   {
      nghosts = SAMRAI::hier::IntVector<NDIM>(1);
   }
   else
   {
      nghosts = SAMRAI::hier::IntVector<NDIM>(0);
   }

   SAMRAI::tbox::Pointer< SAMRAI::hier::Variable<NDIM> > pflotran_var;
   
   if(use_components)
   {
      // in this case regular cell variables are used

      for(int i=0;i<dof; i++)
      {
         std::ostringstream cbuffer;
         std::string vName;

         cbuffer<<(long)i;
         object_str=cbuffer.str();
         vName = dataName+"_component_"+object_str;
        
         PflotranApplicationStrategy::createVariable(vName,centering, 0, 1, pflotran_var);
                           
         pflotran_var_id = variable_db->registerVariableAndContext(pflotran_var,
                                                                   pflotran_cxt,
                                                                   nghosts);
                  
         for(int ln=0; ln<nlevels; ln++)
         {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(pflotran_var_id);
         }
         
         samrai_vec->addComponent(pflotran_var,
                                  pflotran_var_id);
      }
   }
   else
   {
      PflotranApplicationStrategy::createVariable(dataName,centering, 1, dof, pflotran_var);            
      
      pflotran_var_id = variable_db->registerVariableAndContext(pflotran_var,
                                                                pflotran_cxt,
                                                                nghosts);
      
      
      for(int ln=0; ln<nlevels; ln++)
      {
         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
         level->allocatePatchData(pflotran_var_id);
      }

      tbox::Pointer< math::HierarchyDataOpsReal< NDIM, double > > math_op = (centering==0)? d_ccell_math_op:d_cside_math_op; 

      samrai_vec->addComponent(pflotran_var,
                               pflotran_var_id,
                               d_pflotran_weight_id,
                               math_op);
   }

   *vec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM,double>::createPETScVector(samrai_vec, PETSC_COMM_WORLD);
}

void
PflotranApplicationStrategy::writePlotData(int time_step,
                                           double sim_time)
{
   d_visit_writer->writePlotData(d_hierarchy, time_step, sim_time);
}

void
PflotranApplicationStrategy::initializePreconditioner(int *which_pc, PC *pc)
{
   if(*which_pc==0)
   {
      if(d_FlowPreconditioner.get()==NULL)
      {
         tbox::Pointer<tbox::Database> application_db = this->getDatabase();
         tbox::Pointer<tbox::Database> pc_db = application_db->getDatabase("PflotranFlowPreconditioner");
         PflotranFlowPreconditionerParameters *parameters = new PflotranFlowPreconditionerParameters(pc_db);
         parameters->d_hierarchy = d_hierarchy;
         parameters->d_pc = pc;
         parameters->d_cf_interpolant = this->getRefinementBoundaryInterpolant();

         SAMRAI::PflotranFlowPreconditioner *pFlowPC = new SAMRAI::PflotranFlowPreconditioner(parameters);
         d_FlowPreconditioner.reset(pFlowPC);
         d_FlowPreconditioner->setOperator(d_FlowJacobian.get());
      }

   }
   else
   {
      if(d_TransportPreconditioner.get()==NULL)
      {
         tbox::Pointer<tbox::Database> application_db = this->getDatabase();
         tbox::Pointer<tbox::Database> pc_db = application_db->getDatabase("PflotranTransportPreconditioner");
         PflotranTransportPreconditionerParameters *parameters = new PflotranTransportPreconditionerParameters(pc_db);
         parameters->d_hierarchy = d_hierarchy;
         parameters->d_pc = pc;
         parameters->d_cf_interpolant = this->getRefinementBoundaryInterpolant();

         SAMRAI::PflotranTransportPreconditioner *pTransportPC = new SAMRAI::PflotranTransportPreconditioner(parameters);
         d_TransportPreconditioner.reset(pTransportPC);
         d_TransportPreconditioner->setOperator(d_TransportJacobian.get());
      }

   }
}

void 
PflotranApplicationStrategy::createVariable(std::string &vname,
                                            int centering,
                                            int type,
                                            int dof,
                                            SAMRAI::tbox::Pointer< SAMRAI::hier::Variable<NDIM> > &var)
{
   if(centering==0)
   {
      if(type==0)
      {
         var =  new SAMRAI::pdat::CellVariable<NDIM,double>(vname,dof);
      }
      else
      {
         var =  new SAMRAI::pdat::CCellVariable<NDIM,double>(vname,dof);
      }
   }
   else
   {
      if(type==0)
      {
         var =  new SAMRAI::pdat::SideVariable<NDIM,double>(vname,dof);
      }
      else
      {
         var =  new SAMRAI::pdat::CSideVariable<NDIM,double>(vname,dof);
      }
   }
}
}
