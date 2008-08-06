#include "PflotranJacobianMultilevelOperator.h"
#include "CartesianGridGeometry.h"
#include "PflotranJacobianMultilevelOperatorParameters.h"
#include "CCellData.h"
#include "PETSc_SAMRAIVectorReal.h"
#include "tbox/TimerManager.h"
#include "CellDataFactory.h"
#include "RefineAlgorithm.h"

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
   d_face_coarsen_op_str          = "CONSERVATIVE_COARSEN";
   d_cell_coarsen_op_str          = "CONSERVATIVE_COARSEN";
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

   for(int ln=0; ln<hierarchy_size; ln++)
   {
      LevelOperatorParameters *params = new LevelOperatorParameters(parameters->d_db);
      // The next call is important
      // It lets the level operators know what object_id to use as a suffix
      // when creating internal data to minimize the number of variables created
      parameters->d_db->putInteger("object_id", d_object_id);
      
      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

      params->d_level               = level;
      params->d_cf_interpolant      = parameters->d_cf_interpolant;
      params->d_set_boundary_ghosts = d_set_boundary_ghosts;
      d_level_operators[ln]         = new PflotranJacobianLevelOperator(params);
      delete params;
   }
   
   initializeInternalVariableData();

   d_GlobalToLocalRefineSchedule.resizeArray(d_hierarchy->getNumberOfLevels());
   for(int ln=0; ln<hierarchy_size; ln++)
   {
      d_GlobalToLocalRefineSchedule[ln].setNull();
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

   tbox::Pointer< pdat::CellVariable<NDIM,double> > scratchVar = new pdat::CellVariable<NDIM,double>(vecname,d_ndof);
   int scratch_var_id = variable_db->registerVariableAndContext(scratchVar,
                                                                scratch_cxt,
                                                                hier::IntVector<NDIM>(1));
   
   scratch_vector->addComponent(scratchVar,
                                scratch_var_id);
   

   std::string cellFlux("PflotranJacobianMultilevelOperator_InternalFlux");
   cellFlux+=object_str;

   d_flux = variable_db->getVariable(cellFlux);

   if (!d_flux) 
   {
      d_flux = new pdat::FaceVariable<NDIM,double>(cellFlux,1);
   }

   d_flux_id = variable_db->registerVariableAndContext(d_flux,
                                                       scratch_cxt,
                                                       zero_ghosts);

   

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

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(level->checkAllocated(d_flux_id));
#endif

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

   if(db->keyExists("cell_coarsen_op"))
   {
      d_cell_coarsen_op_str = db->getString("cell_coarsen_op");
   }
   else
   {
      TBOX_ERROR("PflotranJacobianMultilevelOperator" 
                 << " -- Required key `cell_coarsen_op'"
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

#if 0
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
#endif
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
}

void
PflotranJacobianMultilevelOperator::setupTransferSchedules(void)
{
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
      
      levelMatrix->MatSetValuesLocal(patchNumber,
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
    tbox::Pointer< pdat::CellDataFactory< NDIM, double > > localFactory = localVar->getPatchDataFactory(); 
    int localDOF = localFactory->getDefaultDepth();

    tbox::Pointer< hier::Variable< NDIM > > globalVar = globalVec->getComponentVariable(0);
    tbox::Pointer< pdat::CellDataFactory< NDIM, double > > globalFactory = globalVar->getPatchDataFactory(); 
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

    // should add code to coarsen variables

    t_interpolate_variable->stop();
}

}

