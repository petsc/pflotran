#ifndef included_CartesianGridGeometry
#include "CartesianGridGeometry.h"
#endif
#ifndef included_CartesianPatchGeometry
#include "CartesianPatchGeometry.h"
#endif

#ifndef included_FaceData
#include "FaceData.h"
#endif

#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "PETSc_SAMRAIVectorReal.h"

#include "MultilevelOperatorParameters.h"

#include "PflotranTransportPreconditioner.h"
 
PflotranTransportPreconditioner::PflotranTransportPreconditioner(PflotranTransportPreconditionerParameters *parameters)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(parameters!=NULL);
#endif

   d_preconditioner_print_flag   = false;

   d_pc_solver_op_registered     = false;

   d_hierarchy                   = parameters->d_hierarchy;

   d_cf_interpolant              = parameters->d_cf_interpolant;

   d_pc                          = parameters->d_pc;

   d_pc_solver_op_registered     = false;

   d_PCSolverFactory             = NULL;

   d_pc_solver                   = NULL;

   d_pc_operator                 = NULL;

   getFromInput(parameters->d_db);   

   initializeSolvers(parameters->d_db);
}

PflotranTransportPreconditioner::~PflotranTransportPreconditioner()
{

   if(d_pc_solver)
   {
      delete d_pc_solver;
      d_pc_solver = NULL;
   }

   if(d_pc_operator)
   {
      delete d_pc_operator;
      d_pc_operator = NULL;
   }

   if (d_PCSolverFactory) 
   {
      delete d_PCSolverFactory;
      d_PCSolverFactory = NULL;
   }
}

void
PflotranTransportPreconditioner::initializePetscInterface(void)
{
   int ierr;

   if(d_pc!=NULL)
   {
      std::string pc_type = "shell";

      ierr = PCSetType(*d_pc,(PCType) pc_type.c_str());
      PETSC_SAMRAI_ERROR(ierr);

      PCShellSetContext(*d_pc, this);

      ierr = PCShellSetSetUp(*d_pc, PflotranTransportPreconditioner::wrapperSetupPreconditioner);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetApply(*d_pc, PflotranTransportPreconditioner::wrapperApplyPreconditioner); 
      PETSC_SAMRAI_ERROR(ierr);       
   }
   else
   {
      abort();
   }
}

void
PflotranTransportPreconditioner::getFromInput( tbox::Pointer<tbox::Database> &db,
                                  bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif
   if (!is_from_restart) 
   {
      if (db->keyExists("pc_solver"))
      {
         d_PCSolverFactory = new MultilevelSolverFactory();         
      } 
      else 
      {
         TBOX_ERROR("PflotranTransportPreconditioner" << " -- Key data `pc_solver'"
                    << " missing in input.");
      }        

      if (db->keyExists("preconditioner_print_flag")) 
      {
         d_preconditioner_print_flag = db->getBool("preconditioner_print_flag");
      } 
      else 
      {
         TBOX_ERROR("PflotranTransportPreconditioner" 
                    << " -- Required key `preconditioner_print_flag'"
                    << " missing in input.");
      }
   }
   else
   {
      abort();
   }
}

void
PflotranTransportPreconditioner::initializeSolvers(tbox::Pointer<tbox::Database> &db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
#endif

   tbox::Pointer<tbox::Database> solver_db = db->getDatabase("pc_solver");
   MultilevelSolverParameters *params      = new  MultilevelSolverParameters(solver_db);
   params->d_hierarchy                     = d_hierarchy;

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_PCSolverFactory!=NULL);
#endif

   std::auto_ptr<MultilevelSolver> pcs          = d_PCSolverFactory->createMultilevelSolver(params);
   d_pc_solver                             = pcs.get();
   pcs.release();
   
   delete params;

}

void 
PflotranTransportPreconditioner::setRefinementBoundaryInterpolant(RefinementBoundaryInterpolation *cf_interpolant)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(cf_interpolant!=NULL);
   assert(d_pc_operator!=NULL);
#endif

   d_cf_interpolant = cf_interpolant;
}

PetscErrorCode
PflotranTransportPreconditioner::wrapperApplyPreconditioner(void *ptr, Vec r, Vec z)
{
   
   PflotranTransportPreconditioner *pc = (PflotranTransportPreconditioner *)ptr;

   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > rVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(r);
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > zVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(z);

   pc->applyPreconditioner(rVec,zVec);

   return(0);
}

/*
*************************************************************************
*                                                                       *
* Apply the FAC Preconditioner to the right-hand side r, placing the    *
* result in the vector z.                                               *
*                                                                       *
*************************************************************************
*/
int 
PflotranTransportPreconditioner::applyPreconditioner(
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > r,
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > z)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!r.isNull());
   assert(!z.isNull());
#endif
   int r_id, z_id;

   r_id = r->getComponentDescriptorIndex(0);
   z_id = z->getComponentDescriptorIndex(0);

   d_pc_solver->solve(&r_id, &z_id);
   return (0);

}   

PetscErrorCode
PflotranTransportPreconditioner::wrapperSetupPreconditioner( void *ptr)
{
   PflotranTransportPreconditioner *pc = (PflotranTransportPreconditioner *)ptr;
   PreconditionerParameters* params = NULL;
   pc->setupPreconditioner(params);

   return(0);
}

int
PflotranTransportPreconditioner::setupPreconditioner( PreconditionerParameters* parameters )
{
   static tbox::Pointer<tbox::Timer> t_setup_pc = tbox::TimerManager::getManager()->getTimer("PflotranTransportPreconditioner::setupPreconditioner");
   t_setup_pc->start();

   if(d_pc_solver_op_registered)
   {
      d_pc_solver->reset();
   }
   else
   {
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_pc_operator!=NULL);
#endif
      d_pc_solver->registerOperator(d_pc_operator);
      d_pc_solver_op_registered = true;
   }

   t_setup_pc->stop();

   return(0);
}

void
PflotranTransportPreconditioner::coarsenVariable( const int var_id,
                                     std::string coarsen_op_str)
{
   xfer::CoarsenAlgorithm<NDIM> coarsen_var;
   tbox::Pointer<xfer::CoarsenOperator<NDIM> > soln_coarsen_op;
   tbox::Pointer< hier::Variable<NDIM> > cvar;
   tbox::Pointer<xfer::Geometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
   
   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();

   variable_db->mapIndexToVariable(var_id, cvar);
   soln_coarsen_op = grid_geometry->lookupCoarsenOperator(cvar, coarsen_op_str);
   
   coarsen_var.registerCoarsen(var_id, var_id, soln_coarsen_op);
   
   // coarsen data to be consistent across all levels
   for ( int ln = d_hierarchy->getFinestLevelNumber();
         ln >= 0;
         ln-- ) 
   {
      
      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      
      if(ln<d_hierarchy->getFinestLevelNumber())
      {
         tbox::Pointer< xfer::CoarsenSchedule<NDIM> > coarsen_schedule = getCoarsenSchedule(ln, coarsen_var);
         coarsen_schedule->coarsenData();
      }
   }   
}

tbox::Pointer<xfer::CoarsenSchedule<NDIM> >
PflotranTransportPreconditioner::getCoarsenSchedule(int ln, xfer::CoarsenAlgorithm<NDIM> &crs_fill_alg)
{
   tbox::Pointer<xfer::CoarsenSchedule<NDIM> > crs_fill_schedule;

   tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
   tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(ln+1);

   crsList* coarsen_fill_schedules = getCoarsenSchedules(ln);

   bool found_compatible_schedule = false;

   for (crsList::Iterator crs(*coarsen_fill_schedules); crs; crs++) 
   {
      if (crs_fill_alg.checkConsistency(crs())) 
      {
         crs_fill_alg.resetSchedule(crs());
         crs_fill_schedule = crs();
         found_compatible_schedule = true;
      }

      if (found_compatible_schedule) break;
   }

   if (!found_compatible_schedule) 
   {
      crs_fill_schedule = crs_fill_alg.createSchedule(level, flevel);
      coarsen_fill_schedules->appendItem(crs_fill_schedule);
   }
   
   return crs_fill_schedule;

}
