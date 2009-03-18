#include "PflotranJacobianMultilevelOperator.h"
#include "CartesianGridGeometry.h"
#include "PflotranJacobianMultilevelOperatorParameters.h"
#include "CCellData.h"
#include "PETSc_SAMRAIVectorReal.h"
#include "tbox/TimerManager.h"
#include "CCellDataFactory.h"
#include "RefineAlgorithm.h"
#include "CoarsenAlgorithm.h"
#include "HierarchyCCellDataOpsReal.h"
#include "CartesianGridGeometry.h"

namespace SAMRAI{
#ifndef iC
#define iC(fun) {CHKERRQ(fun);}
#endif

PflotranJacobianMultilevelOperator::PflotranJacobianMultilevelOperator()
{
}

PflotranJacobianMultilevelOperator::~PflotranJacobianMultilevelOperator()
{
}

PflotranJacobianMultilevelOperator::PflotranJacobianMultilevelOperator(MultilevelOperatorParameters *parameters)
   :MultilevelLinearOperator(parameters)
{
   d_flux_id                      = -1;
   d_coarsen_diffusive_fluxes     = true;
   d_schedules_initialized        = false;
   d_adjust_cf_coefficients       = false; 
   d_variable_order_interpolation = false;
   d_reset_ghost_values           = true;
   d_face_coarsen_op_str          = "SUM_COARSEN";
   d_cell_soln_coarsen_op_str     = "CONSERVATIVE_COARSEN";
   d_cell_src_coarsen_op_str      = "SUM_COARSEN";
   d_cell_refine_op_str           = "CONSTANT_REFINE";
   d_face_refine_op_str           = "CONSTANT_REFINE";
   d_flux.setNull();
   
   d_patch                        = NULL;
   d_refine_patch_strategy        = new BoundaryConditionStrategy(-1);

   const int hierarchy_size       = d_hierarchy->getNumberOfLevels();

   d_level_operators.resizeArray(hierarchy_size);

   PflotranJacobianMultilevelOperatorParameters *pm = dynamic_cast<PflotranJacobianMultilevelOperatorParameters *>(parameters);

   d_pMatrix = pm->d_pMatrix;

   MatShellSetContext(*d_pMatrix,this);

   PflotranJacobianMultilevelOperator::getFromInput(parameters->d_db);

   initializePetscMatInterface();

   // we initialize the internal variable data before
   // initializing the level operators so that the variable id
   // for the src/sink's is set and can be passed to the level operators
   initializeInternalVariableData();

   for(int ln=0; ln<hierarchy_size; ln++)
   {
      LevelOperatorParameters *params = new LevelOperatorParameters(parameters->d_db);
      // The next call is important
      // It lets the level operators know what object_id to use as a suffix
      // when creating internal data to minimize the number of variables created
      parameters->d_db->putInteger("object_id", d_object_id);
      parameters->d_db->putInteger("srcsink_id", d_srcsink_id);
      
      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

      params->d_level               = level;
      params->d_cf_interpolant      = parameters->d_cf_interpolant;
      params->d_set_boundary_ghosts = d_set_boundary_ghosts;
      d_level_operators[ln]         = new PflotranJacobianLevelOperator(params);
      delete params;
   }
   
   d_math_op = new math::HierarchyCCellDataOpsReal< NDIM, double >(d_hierarchy,
                                                                   0, d_hierarchy->getFinestLevelNumber());
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();

   d_soln_refine_op =  grid_geometry->lookupRefineOperator(d_scratch_variable, d_cell_refine_op_str);

   d_soln_coarsen_op = grid_geometry->lookupCoarsenOperator(d_scratch_variable, d_cell_soln_coarsen_op_str);

   d_src_coarsen_op = grid_geometry->lookupCoarsenOperator(d_scratch_variable, d_cell_src_coarsen_op_str);

   d_flux_coarsen_schedule.resizeArray(d_hierarchy->getNumberOfLevels());
   d_GlobalToLocalRefineSchedule.resizeArray(d_hierarchy->getNumberOfLevels());
   d_src_coarsen_schedule.resizeArray(d_hierarchy->getNumberOfLevels());
   d_soln_coarsen_schedule.resizeArray(d_hierarchy->getNumberOfLevels());
   for(int ln=0; ln<hierarchy_size; ln++)
   {
      d_GlobalToLocalRefineSchedule[ln].setNull();
      d_src_coarsen_schedule[ln].setNull();
      d_soln_coarsen_schedule[ln].setNull();
      d_flux_coarsen_schedule[ln].setNull();
   }

}

void
PflotranJacobianMultilevelOperator::initializePetscMatInterface(void)
{
   MatShellSetOperation(*d_pMatrix, MATOP_MULT, (void(*)(void))(&SAMRAI::PflotranJacobianMultilevelOperator::wrapperMatMult));

   MatShellSetOperation(*d_pMatrix, MATOP_ZERO_ENTRIES, (void(*)(void))(&SAMRAI::PflotranJacobianMultilevelOperator::wrapperMatZeroEntries));

   MatShellSetOperation(*d_pMatrix, MATOP_SET_VALUES_LOCAL, (void(*)(void))(&SAMRAI::PflotranJacobianMultilevelOperator::wrapperMatSetValuesLocal));

   MatShellSetOperation(*d_pMatrix, MATOP_SET_VALUES_BLOCKED, (void(*)(void))(&SAMRAI::PflotranJacobianMultilevelOperator::wrapperMatSetValuesBlockedLocal));
}

void
PflotranJacobianMultilevelOperator::apply(const int ln,
                                          const int *f_id,
                                          const int *u_id, 
                                          const int *r_id,
                                          const int *f_idx,
                                          const int *u_idx,
                                          const int *r_idx,
                                          const double a,
                                          const double b)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(f_id!=NULL);
   assert(u_id!=NULL);
   assert(r_id!=NULL);
   assert(d_level_operators[ln]!=NULL);
#endif

   d_level_operators[ln]->apply(f_id , u_id , r_id ,
                                f_idx, u_idx, r_idx,
                                a, b);
}

void
PflotranJacobianMultilevelOperator::apply(const int coarse_ln,
                                          const int fine_ln,
                                          const int *f_id,
                                          const int *u_id, 
                                          const int *r_id,
                                          const int *f_idx,
                                          const int *u_idx,
                                          const int *r_idx,
                                          const double a,
                                          const double b)
{
   if(d_coarsen_diffusive_fluxes)
   {
      setupTransferSchedules();
      
      for(int ln=fine_ln; ln>=coarse_ln; ln--)
      {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(d_level_operators[ln]!=NULL);
#endif
         d_level_operators[ln]->setFlux(d_flux_id, u_id, u_idx);
      }

      for(int ln=fine_ln-1; ln>=coarse_ln; ln--)
      {
         bool coarsen_rhs = (b!=0.0)?true:false;   
         bool coarsen_soln = false;
         coarsenSolutionAndSourceTerm(ln, u_id[0], f_id[0], coarsen_soln, coarsen_rhs);      
         
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(d_schedules_initialized==true);
         assert(!d_flux_coarsen_schedule[ln].isNull());
#endif
         d_flux_coarsen_schedule[ln]->coarsenData();
      }

      for(int ln=fine_ln; ln>=coarse_ln; ln--)
      {
         d_level_operators[ln]->apply(d_flux_id,
                                      f_id , u_id , r_id, 
                                      f_idx, u_idx, r_idx,
                                      a, b);
      }

//       for(int ln=fine_ln-1; ln>=coarse_ln; ln--)
//       {
//          bool coarsen_rhs = (b!=0.0)?true:false;         
//          coarsenSolutionAndSourceTerm(ln, u_id[0], r_id[0], coarsen_rhs);      
//       }
   }
   else
   {
      // this version does not compute correctly at c-f interfaces
      for(int ln=fine_ln; ln>=coarse_ln; ln--)
      {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(d_level_operators[ln]!=NULL);
#endif
         d_level_operators[ln]->apply(f_id, u_id, r_id, f_idx, u_idx, r_idx, a, b);
         
      }
   }
}

void
PflotranJacobianMultilevelOperator::applyBoundaryCondition(const int ln,
                                                           const int *var_id,
                                                           const int *var_idx,
                                                           const int *var_components,
                                                           const int number_of_variables,
                                                           const bool reset_ghost_values)
{
}

void
PflotranJacobianMultilevelOperator::initializeInternalVariableData(void)
{
   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();

   const hier::IntVector<NDIM> zero_ghosts(0);

   const tbox::Pointer< hier::VariableContext > scratch_cxt = variable_db->getContext("SCRATCH");

   std::ostringstream ibuffer;
   ibuffer<<(long)d_object_id;
   std::string object_str=ibuffer.str();

   std::string vecname("PflotranJacobianMultilevelOperator_scratchVector");
   vecname+=object_str;

   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > scratch_vector = new solv::SAMRAIVectorReal<NDIM,double>(vecname,
                                                                                                                 d_hierarchy,
                                                                                                                 0, d_hierarchy->getFinestLevelNumber());

   d_scratch_variable = new pdat::CCellVariable<NDIM,double>(vecname,d_ndof);
   int scratch_var_id = variable_db->registerVariableAndContext(d_scratch_variable,
                                                                scratch_cxt,
                                                                hier::IntVector<NDIM>(1));
   
   scratch_vector->addComponent(d_scratch_variable,
                                scratch_var_id, -1, d_math_op);
   

   std::string cellFlux("PflotranJacobianMultilevelOperator_InternalFlux");
   cellFlux+=object_str;

   d_flux = variable_db->getVariable(cellFlux);

   if (!d_flux) 
   {
      d_flux = new pdat::CSideVariable<NDIM,double>(cellFlux,1);
   }

   d_flux_id = variable_db->registerVariableAndContext(d_flux,
                                                       scratch_cxt,
                                                       zero_ghosts);

   
   std::string SrcSink("PflotranJacobianOperator_SrcSink");
   SrcSink+=object_str;

#ifdef DEBUG_CHECK_ASSERTIONS
   if(d_debug_print_info_level>4)
   {
      tbox::pout << "PflotranJacobianMultilevelOperator::SrcSink variable name " << SrcSink << std::endl;
   }
#endif

   d_srcsink = variable_db->getVariable(SrcSink);

   if (!d_srcsink) 
   {
      d_srcsink = new pdat::CCellVariable<NDIM,double>(SrcSink,d_ndof);
   }

   if(d_srcsink)
   {
      d_srcsink_id = variable_db->registerVariableAndContext(d_srcsink,
                                                             scratch_cxt,
                                                             zero_ghosts);
   }
   else
   {
      tbox::pout << "PflotranJacobianLevelOperator::ERROR could not allocate memory for internal srcsink data" << std::endl; 
      abort();
   }

   for(int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++)
   {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      // allocate storage space
      if(!level->checkAllocated(d_flux_id))
      {
         level->allocatePatchData(d_flux_id);
      }

      if(!level->checkAllocated(scratch_var_id))
      {
         level->allocatePatchData(scratch_var_id);
      }

      // allocate storage space
      if(!level->checkAllocated(d_srcsink_id))
      {
         level->allocatePatchData(d_srcsink_id);
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(level->checkAllocated(d_flux_id));
#endif

      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) 
      {
         tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
         tbox::Pointer< pdat::CCellData<NDIM, double > > data = patch->getPatchData(d_srcsink_id);
         data->fillAll(0.0);
         data = patch->getPatchData(scratch_var_id);
         data->fillAll(0.0);
      }
   }

   d_scratch_vector = solv::PETSc_SAMRAIVectorReal<NDIM,double>::createPETScVector(scratch_vector);
}

void
PflotranJacobianMultilevelOperator::getFromInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   if (db->keyExists("tangent_interp_scheme")) 
   {
      d_tangential_interp_scheme = RefinementBoundaryInterpolation::lookupInterpolationScheme(db->getString("tangent_interp_scheme"));
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianMultilevelLevelOperator" 
                 << " -- Required key `tangent_interp_scheme'"
                 << " missing in input.");
   }

   if (db->keyExists("normal_interp_scheme")) 
   {
      d_normal_interp_scheme = RefinementBoundaryInterpolation::lookupInterpolationScheme(db->getString("normal_interp_scheme"));
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianMultilevelLevelOperator" 
                 << " -- Required key `normal_interp_scheme'"
                 << " missing in input.");
   }

   if(db->keyExists("cell_src_coarsen_op"))
   {
      d_cell_src_coarsen_op_str = db->getString("cell_src_coarsen_op");
   }
   else
   {
      TBOX_ERROR("PflotranJacobianMultilevelOperator" 
                 << " -- Required key `cell_src_coarsen_op'"
                 << " missing in input.");
   }

   if(db->keyExists("cell_soln_coarsen_op"))
   {
      d_cell_soln_coarsen_op_str = db->getString("cell_soln_coarsen_op");
   }
   else
   {
      TBOX_ERROR("PflotranJacobianMultilevelOperator" 
                 << " -- Required key `cell_soln_coarsen_op'"
                 << " missing in input.");
   }

   if(db->keyExists("use_cf_interpolant"))
   {
      d_use_cf_interpolant = db->getBool("use_cf_interpolant");
      
      if(d_use_cf_interpolant)
      {
         d_cell_refine_op_str = "CONSTANT_REFINE";

         if (db->keyExists("variable_order_interpolation")) 
         {
            d_variable_order_interpolation = db->getBool("variable_order_interpolation");
         } 
         else 
         {
            TBOX_ERROR("PflotranJacobianMultilevelOperator"
                       << " -- Required key `variable_order_interpolation'"
                       << " missing in input.");
         }
      }
      else
      {
         if(db->keyExists("cell_refine_op"))
         {
            d_cell_refine_op_str = db->getString("cell_refine_op");
         }
         else
         {
            TBOX_ERROR( "PflotranJacobianMultilevelOperator"
                       << " -- Required key `cell_refine_op'"
                       << " missing in input.");
         }
      }
   }
   else
   {
      TBOX_ERROR("PflotranJacobianMultilevelOperator"
                 << " -- Required key `use_cf_interpolant'"
                 << " missing in input.");
   }

   if (db->keyExists("coarsen_diffusive_fluxes")) 
   {
      d_coarsen_diffusive_fluxes = db->getBool("coarsen_diffusive_fluxes");
   } 
   else 
   {
      TBOX_ERROR("PflotranJacobianMultilevelOperator" 
                 << " -- Required key `coarsen_diffusive_fluxes'"
                 << " missing in input.");
   }
   
   if(db->keyExists("face_coarsen_op"))
   {
      d_face_coarsen_op_str = db->getString("face_coarsen_op");
   }
   else
   {
      TBOX_ERROR("PflotranJacobianMultilevelOperator" 
                 << " -- Required key `face_coarsen_op'"
                 << " missing in input.");
   }

   if(db->keyExists("face_refine_op"))
   {
      d_face_refine_op_str = db->getString("face_refine_op");
   }
   else
   {
      TBOX_ERROR("PflotranJacobianMultilevelOperator" 
                 << " -- Required key `face_refine_op'"
                 << " missing in input.");
   }

   if (db->keyExists("adjust_cf_coefficients")) 
   {
      d_adjust_cf_coefficients = db->getBool("adjust_cf_coefficients");      
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianLevelOperator" 
                 << " -- Required key `adjust_cf_coefficients'"
                 << " missing in input.");
   }

   if (db->keyExists("boundary_conditions")) 
   {
      // get the database object for boundary conditions
      db->getIntegerArray("boundary_conditions", (int *)&d_bdry_types, 2*NDIM);
      
      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = 
         d_hierarchy->getGridGeometry();
      hier::IntVector<NDIM> shift = grid_geometry->getPeriodicShift();
      
      for (int d=0; d<NDIM; d++) 
      {
         if (shift[d]!=0)
         {
            d_bdry_types[2*d+0] = PERIODIC;
            d_bdry_types[2*d+1] = PERIODIC;
         }
      }

      initializeBoundaryConditionStrategy(db);

   } 
   else 
   {
      TBOX_ERROR("PflotranJacobianMultilevelOperator" 
                 << " -- Required key `boundary_conditions'"
                 << " missing in input.");
   } 

   if (db->keyExists("stencilsize")) 
   {
      d_stencilSize = db->getInteger("stencilsize");      
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianMultilevelOperator" 
                 << " -- Required key `stencilsize'"
                 << " missing in input.");
   }

   if (db->keyExists("ndof")) 
   {
      d_ndof = db->getInteger("ndof");      
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianMultilevelOperator" 
                 << " -- Required key `ndof'"
                 << " missing in input.");
   }
}

void
PflotranJacobianMultilevelOperator::setFlux(const int coarse_ln,
                                            const int fine_ln,
                                            const int *u_id,
                                            const int *u_idx)
{

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_flux_id>=0);
      assert(u_id!=NULL);
#endif   

   for(int ln=fine_ln; ln>=coarse_ln; ln--)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_level_operators[ln]!=NULL);
#endif
      
      d_level_operators[ln]->setFlux(d_flux_id, u_id, u_idx);

      if(ln<fine_ln && d_coarsen_diffusive_fluxes)
      {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(d_schedules_initialized==true);
         assert(!d_flux_coarsen_schedule[ln].isNull());
#endif
         d_flux_coarsen_schedule[ln]->coarsenData();
      } 
   }
}

void
PflotranJacobianMultilevelOperator::setupTransferSchedules(void)
{
   int ln;
   
   if(d_coarsen_diffusive_fluxes && (!d_schedules_initialized))
   {
      tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(0);
      tbox::Pointer<hier::PatchLevel<NDIM> > clevel;

      tbox::Pointer<xfer::Geometry<NDIM> > geometry = flevel->getGridGeometry();

      xfer::CoarsenAlgorithm<NDIM> flux_coarsen_alg;
      flux_coarsen_alg.registerCoarsen(d_flux_id, d_flux_id,
                                       geometry->lookupCoarsenOperator(d_flux,d_face_coarsen_op_str));
      
      for (ln = 1; ln < d_hierarchy->getNumberOfLevels(); ln++) 
      {
         flevel = d_hierarchy->getPatchLevel(ln);
         clevel = d_hierarchy->getPatchLevel(ln-1);
         d_flux_coarsen_schedule[ln-1] = flux_coarsen_alg.createSchedule(clevel, flevel); 
      }
      
      d_schedules_initialized=true;

   }
}

void 
PflotranJacobianMultilevelOperator::initializeBoundaryConditionStrategy(tbox::Pointer<tbox::Database> &db)
{
}

PetscErrorCode
PflotranJacobianMultilevelOperator::wrapperMatMult(Mat mat,Vec x,Vec y)
{

   PflotranJacobianMultilevelOperator *pMatrix = NULL;
   
   MatShellGetContext(mat,(void**)&pMatrix);

   return pMatrix->MatMult(mat,x,y);

}

PetscErrorCode
PflotranJacobianMultilevelOperator::MatMult(Mat mat,Vec x,Vec y)
{

   PflotranJacobianMultilevelOperator *pMatrix = NULL;
   
   MatShellGetContext(mat,(void**)&pMatrix);

   pMatrix->initializeScratchVector(x);

   Vec scratch = pMatrix->getScratchVector();

   for(int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++)
   {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      
      PflotranJacobianLevelOperator *levelMatrix = dynamic_cast<PflotranJacobianLevelOperator *>(pMatrix->getLevelOperator(ln));

      levelMatrix->MatMult(scratch,y);
   }

   return (0);

}


PetscErrorCode  
PflotranJacobianMultilevelOperator::wrapperMatZeroEntries(Mat mat)
{

  PflotranJacobianMultilevelOperator *pMatrix = NULL;
  MatShellGetContext(mat,(void**)&pMatrix);
  return pMatrix->MatZeroEntries(mat);

}

PetscErrorCode  
PflotranJacobianMultilevelOperator::MatZeroEntries(Mat mat)
{
  PflotranJacobianMultilevelOperator *pMatrix = NULL;
   
   MatShellGetContext(mat,(void**)&pMatrix);

   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy = pMatrix->getHierarchy();

   for(int ln=0; ln<hierarchy->getNumberOfLevels(); ln++)
   {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      
      PflotranJacobianLevelOperator *levelMatrix = dynamic_cast<PflotranJacobianLevelOperator *>(pMatrix->getLevelOperator(ln));
      
      int stencil_id = levelMatrix->getStencilID();

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(level->checkAllocated(stencil_id));
#endif

      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) 
      {
         tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());

#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!patch.isNull());
#endif
         tbox::Pointer< pdat::CCellData<NDIM, double > > stencil = patch->getPatchData(stencil_id);
         stencil->fillAll(0.0);
      }
   }   

   return (0);

}

PetscErrorCode  
PflotranJacobianMultilevelOperator::wrapperMatSetValuesLocal(Mat mat,
                                                             PetscInt nrow,const PetscInt irow[],
                                                             PetscInt ncol,const PetscInt icol[],
                                                             const PetscScalar y[],InsertMode addv)
{
  PflotranJacobianMultilevelOperator *pMatrix = NULL;
   
   MatShellGetContext(mat,(void**)&pMatrix);

   return pMatrix->MatSetValuesLocal(mat, nrow, irow, ncol, icol, y, addv);
}

PetscErrorCode  
PflotranJacobianMultilevelOperator::MatSetValuesLocal(Mat mat,
                                                      PetscInt nrow,const PetscInt irow[],
                                                      PetscInt ncol,const PetscInt icol[],
                                                      const PetscScalar y[],InsertMode addv)
{
  PflotranJacobianMultilevelOperator *pMatrix = NULL;
   
   MatShellGetContext(mat,(void**)&pMatrix);
   
   if(pMatrix->d_patch !=NULL)
   {
      int ln = d_patch->getPatchLevelNumber();
      int patchNumber = d_patch->getPatchNumber();

      PflotranJacobianLevelOperator *levelMatrix = dynamic_cast<PflotranJacobianLevelOperator *>(pMatrix->getLevelOperator(ln));

      levelMatrix->MatSetValuesLocal(patchNumber,
                                     nrow, irow,
                                     ncol, icol,
                                     y, addv);
      
   }

   return (0);
}

PetscErrorCode  
PflotranJacobianMultilevelOperator::wrapperMatSetValuesBlockedLocal(Mat mat,
                                                                    PetscInt nrow,const PetscInt irow[],
                                                                    PetscInt ncol,const PetscInt icol[],
                                                                    const PetscScalar y[],InsertMode addv)
{
  PflotranJacobianMultilevelOperator *pMatrix = NULL;
   
   MatShellGetContext(mat,(void**)&pMatrix);
   return pMatrix->MatSetValuesBlockedLocal(mat, nrow, irow, ncol, icol, y, addv);
}

PetscErrorCode  
PflotranJacobianMultilevelOperator::MatSetValuesBlockedLocal(Mat mat,
                                                             PetscInt nrow,const PetscInt irow[],
                                                             PetscInt ncol,const PetscInt icol[],
                                                             const PetscScalar y[],InsertMode addv)
{
  PflotranJacobianMultilevelOperator *pMatrix = NULL;
   
   MatShellGetContext(mat,(void**)&pMatrix);

   if(pMatrix->d_patch !=NULL)
   {
      int ln = d_patch->getPatchLevelNumber();
      int patchNumber = d_patch->getPatchNumber();

      PflotranJacobianLevelOperator *levelMatrix = dynamic_cast<PflotranJacobianLevelOperator *>(pMatrix->getLevelOperator(ln));
      
      levelMatrix->MatSetValuesBlockedLocal(patchNumber,
                                            nrow, irow,
                                            ncol, icol,
                                            y, addv);
      
   }
   return (0);
}

void
PflotranJacobianMultilevelOperator::initializeScratchVector( Vec x )
{
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy =  d_hierarchy;

    static tbox::Pointer<tbox::Timer> t_interpolate_variable = tbox::TimerManager::getManager()->getTimer("PFlotran::PflotranJacobianMultilevelOperator::initializeScratchVector");

    t_interpolate_variable->start();

    SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > globalVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(x);

    SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > localVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(d_scratch_vector);

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

      if((!d_GlobalToLocalRefineSchedule[ln].isNull()) && ghost_cell_fill.checkConsistency(d_GlobalToLocalRefineSchedule[ln]))
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
      
      d_GlobalToLocalRefineSchedule[ln]->fillData(0.0);            

      if(ln>0)
      {	
          d_cf_interpolant->setVariableOrderInterpolation(d_variable_order_interpolation);
          d_cf_interpolant->setPhysicalCornerRefGhostValues(ln,
                                                            dest_id,
                                                            0);
          d_cf_interpolant->setGhostCellData(ln, dest_id);

	  for ( int idx=0; idx<globalDOF; idx++)
	  {
             d_cf_interpolant->interpolateGhostValues(ln,
                                                      d_tangential_interp_scheme,
                                                      d_normal_interp_scheme,
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

      if((!d_soln_coarsen_schedule[ln].isNull()) && cell_coarsen.checkConsistency(d_soln_coarsen_schedule[ln]))
      {
         cell_coarsen.resetSchedule(d_soln_coarsen_schedule[ln]);
      }
      else
      {
         d_soln_coarsen_schedule[ln] = cell_coarsen.createSchedule(clevel, 
                                                             flevel);
      }
      
      d_soln_coarsen_schedule[ln]->coarsenData();            
      
    }
#endif
    t_interpolate_variable->stop();
}

int
PflotranJacobianMultilevelOperator::getVariableIndex(std::string &name, 
                                                     tbox::Pointer<hier::VariableContext> &context,
                                                     tbox::Pointer<hier::Variable<NDIM> > &var,
                                                     hier::IntVector<NDIM> nghosts,
                                                     int depth,
                                                     bool bOverride,
                                                     std::string centering)
{
   int var_id = -1;
   if(!bOverride)
   {
      hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();
      
      var = variable_db->getVariable(name);
      
      if(!var)
      {
         var = new pdat::CCellVariable<NDIM, double>(name, depth);         
      }

      var_id = variable_db->registerVariableAndContext(var,
                                                       context,
                                                       nghosts);
   }
   else
   {
      abort();
   }

   return var_id;

}

void
PflotranJacobianMultilevelOperator::setSourceValueOnPatch(SAMRAI::hier::Patch<NDIM> **patch, 
                                                          int *index, 
                                                          double *val)
{
   int ln = (*patch)->getPatchLevelNumber();

   d_level_operators[ln]->setSourceValueOnPatch(patch, index, val);

}

void
PflotranJacobianMultilevelOperator::setSrcCoefficientsOnPatch(SAMRAI::hier::Patch<NDIM> **patch)
{
   int ln = (*patch)->getPatchLevelNumber();

   d_level_operators[ln]->setSrcCoefficientsOnPatch(patch);

}

void
PflotranJacobianMultilevelOperator::coarsenSolutionAndSourceTerm(const int ln, 
                                                                 const int u_id,
                                                                 const int f_id, 
                                                                 const bool coarsen_soln,
                                                                 const bool coarsen_rhs)
{

   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();
   tbox::Pointer< hier::Variable< NDIM > > f;
   variable_db->mapIndexToVariable(f_id, f);
   tbox::Pointer< hier::Variable< NDIM > > u;
   variable_db->mapIndexToVariable(u_id, u);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!f.isNull());      
   assert(!u.isNull());      
   assert(ln<d_hierarchy->getFinestLevelNumber());
#endif

   tbox::Pointer< xfer::Geometry<NDIM> > geometry = d_hierarchy->getGridGeometry();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!geometry.isNull());      
#endif

   xfer::CoarsenAlgorithm<NDIM> coarsen_src_alg;
   xfer::CoarsenAlgorithm<NDIM> coarsen_soln_alg;

   if(coarsen_rhs)
   {
      coarsen_src_alg.registerCoarsen(f_id, f_id,
                                      geometry->lookupCoarsenOperator(f, d_cell_src_coarsen_op_str));
   

      // if a schedule does not currently exist create it
      if(d_src_coarsen_schedule[ln].isNull())
      {
         tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(ln+1);
         tbox::Pointer<hier::PatchLevel<NDIM> > clevel = d_hierarchy->getPatchLevel(ln);
         d_src_coarsen_schedule[ln]=coarsen_src_alg.createSchedule(clevel, flevel); 
      }
      else
      {
         // recreating a schedule when it is not consistent is not the
         // most efficient way of doing this. An improvement would be
         // to cache multiple schedules 
         if(coarsen_src_alg.checkConsistency(d_src_coarsen_schedule[ln]))
         {
            coarsen_src_alg.resetSchedule(d_src_coarsen_schedule[ln]);
         }
         else
         {
            tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(ln+1);
            tbox::Pointer<hier::PatchLevel<NDIM> > clevel = d_hierarchy->getPatchLevel(ln);
            d_src_coarsen_schedule[ln]=coarsen_src_alg.createSchedule(clevel, flevel); 
            tbox::pout << "CellDiffusionMultilevelOperator::coarsenSolutionAndSourceTerm()::Forced to recreate src coarsen schedule " << std::endl;
         }
      }
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!d_src_coarsen_schedule[ln].isNull());      
#endif
      
      d_src_coarsen_schedule[ln]->coarsenData();   
   }

   if(coarsen_soln)
   {
      coarsen_soln_alg.registerCoarsen(u_id, u_id,
                                      geometry->lookupCoarsenOperator(u, d_cell_soln_coarsen_op_str));
   

      // if a schedule does not currently exist create it
      if(d_soln_coarsen_schedule[ln].isNull())
      {
         tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(ln+1);
         tbox::Pointer<hier::PatchLevel<NDIM> > clevel = d_hierarchy->getPatchLevel(ln);
         d_soln_coarsen_schedule[ln]=coarsen_soln_alg.createSchedule(clevel, flevel); 
      }
      else
      {
         // recreating a schedule when it is not consistent is not the
         // most efficient way of doing this. An improvement would be
         // to cache multiple schedules 
         if(coarsen_soln_alg.checkConsistency(d_soln_coarsen_schedule[ln]))
         {
            coarsen_soln_alg.resetSchedule(d_soln_coarsen_schedule[ln]);
         }
         else
         {
            tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(ln+1);
            tbox::Pointer<hier::PatchLevel<NDIM> > clevel = d_hierarchy->getPatchLevel(ln);
            d_soln_coarsen_schedule[ln]=coarsen_soln_alg.createSchedule(clevel, flevel); 
            tbox::pout << "CellDiffusionMultilevelOperator::coarsenSolutionAndSourceTerm()::Forced to recreate soln coarsen schedule " << std::endl;
         }
      }
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!d_soln_coarsen_schedule[ln].isNull());      
#endif
      
      d_soln_coarsen_schedule[ln]->coarsenData();   
   }
}


}

