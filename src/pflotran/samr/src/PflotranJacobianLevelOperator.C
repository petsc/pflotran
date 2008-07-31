#include "PflotranJacobianLevelOperator.h"
#include "CartesianGridGeometry.h"

namespace SAMRAI{

PflotranJacobianLevelOperator::PflotranJacobianLevelOperator(LevelOperatorParameters *parameters):LevelLinearOperator(parameters)
{
   PflotranJacobianLevelOperator::getFromInput(parameters->d_db);

   initializeInternalVariableData();
}


void
PflotranJacobianLevelOperator::initializeInternalVariableData()
{
   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();

   const hier::IntVector<NDIM> zero_ghosts(0);

   const tbox::Pointer< hier::VariableContext > scratch_cxt = variable_db->getContext("SCRATCH");

   std::ostringstream ibuffer;
   ibuffer<<d_object_id;
   std::string object_str=ibuffer.str();
   std::string cellFlux("PflotranJacobianOperator_InternalFlux");
   cellFlux+=object_str;

#ifdef DEBUG_CHECK_ASSERTIONS
   if(d_debug_print_info_level>4)
   {
      const int ln = (d_level->inHierarchy())?d_level->getLevelNumber():-1;
      tbox::pout << "PflotranJacobianLevelOperator::level number " << ln << std::endl;
      tbox::pout << "PflotranJacobianLevelOperator::Cell Flux variable name " << cellFlux << std::endl;
      
   }
#endif

   d_flux = variable_db->getVariable(cellFlux);

   if (!d_flux) 
   {
      d_flux = new pdat::FaceVariable<NDIM,double>(cellFlux,1);
   }

   d_flux_id = variable_db->registerVariableAndContext(d_flux,
                                                       scratch_cxt,
                                                       zero_ghosts);

   // allocate storage space
   if(!d_level->checkAllocated(d_flux_id))
   {
      d_level->allocatePatchData(d_flux_id);
   }

   std::string Stencil("PflotranJacobianOperator_Stencil");
   Stencil+=object_str;

#ifdef DEBUG_CHECK_ASSERTIONS
   if(d_debug_print_info_level>4)
   {
      tbox::pout << "PflotranJacobianLevelOperator::stencil variable name " << Stencil << std::endl;
   }
#endif

   d_stencil = variable_db->getVariable(Stencil);

   if (!d_stencil) 
   {
      d_stencil = new pdat::CCellVariable<NDIM,double>(Stencil,d_stencil_size*d_ndof*d_ndof);
   }

   if(d_stencil)
   {
      d_stencil_id = variable_db->registerVariableAndContext(d_stencil,
                                                             scratch_cxt,
                                                             zero_ghosts);
   }
   else
   {
      tbox::pout << "PflotranJacobianLevelOperator::ERROR could not allocate memory for internal stencil data" << std::endl; 
      abort();
   }
   // allocate storage space
   if(!d_level->checkAllocated(d_stencil_id))
   {
      d_level->allocatePatchData(d_stencil_id);
   }
}

PflotranJacobianLevelOperator::PflotranJacobianLevelOperator()
{
}

PflotranJacobianLevelOperator::~PflotranJacobianLevelOperator()
{
}

void 
PflotranJacobianLevelOperator::apply(const int *f_id,
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
PflotranJacobianLevelOperator::applyBoundaryCondition(const int *var_id,
                                                      const int *var_idx,
                                                      const int *var_components,
                                                      const int number_of_variables)
{
}

const int 
PflotranJacobianLevelOperator::getStencilType(const int i,
                                              const int j,
                                              const int k)
{
   return(0);
}

tbox::Pointer< hier::PatchData<NDIM > >
PflotranJacobianLevelOperator::getStencilBlock(const int p, 
                                               const int i, 
                                               const int j,
                                               const int k)
{
}

const int
PflotranJacobianLevelOperator::getStencilSize(const int i,
                                              const int j,
                                              const int k)
{
   return (0);
}

std::vector<int>
PflotranJacobianLevelOperator::getStencilOffsets(const int i, 
                                                 const int j, 
                                                 const int k)
{
   std::vector<int> offsets;
  
   return offsets;
}

void 
PflotranJacobianLevelOperator::
getFromInput(const tbox::Pointer<tbox::Database> &db)
{
   LevelLinearOperator::getFromInput(db);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif
#if 0
   if (db->keyExists("tangent_interp_scheme")) 
   {
      d_tangent_interp_scheme = RefinementBoundaryInterpolation::lookupInterpolationScheme(db->getString("tangent_interp_scheme"));
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianLevelOperator" 
                 << " -- Required key `tangent_interp_scheme'"
                 << " missing in input.");
   }

   if (db->keyExists("normal_interp_scheme")) 
   {
      d_normal_interp_scheme = RefinementBoundaryInterpolation::lookupInterpolationScheme(db->getString("normal_interp_scheme"));
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianLevelOperator" 
                 << " -- Required key `normal_interp_scheme'"
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

   if (db->keyExists("variable_order_interpolation")) 
   {
      d_variable_order_interpolation = db->getBool("variable_order_interpolation");
   } 
   else 
   {
      TBOX_ERROR("PflotranJacobianLevelOperator" 
                 << " -- Required key `variable_order_interpolation'"
                 << " missing in input.");
   }
   // this key needs to be explicitly specified though on first appearance
   // it may appear that it is set to true only when adjust_cf_coefficients is false.
   // However, it is needed when ghost cell data
   // is already aligned from the coarser level and we do not
   // want to disturb those values. The situation where this arises
   // is on restricted patches underlying fine patch levels during
   // the AFAC or AFACx solution process
   if (db->keyExists("interpolate_ghost_values")) 
   {
      d_interpolate_ghost_values = db->getBool("interpolate_ghost_values");
      
   } 
   else 
   {
      TBOX_ERROR("PflotranJacobianLevelOperator" 
                 << " -- Required key `interpolate_ghost_values'"
                 << " missing in input.");
   }

   if (db->keyExists("extrapolation_order")) {
      d_extrapolation_order = db->getInteger("extrapolation_order");
   } else {
      TBOX_ERROR("PflotranJacobianLevelOperator"
                 << " -- Required key `extrapolation_order'"
                 << " missing in input.");
   }


   if (db->keyExists("boundary_conditions")) 
   {
      // get the database object for boundary conditions
      db->getIntegerArray("boundary_conditions", (int *)&d_bdry_types, 2*NDIM);

      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = 
         d_level->getGridGeometry();
      hier::IntVector<NDIM> shift = grid_geometry->getPeriodicShift();
      
      for (int d=0; d<NDIM; d++) 
      {
         if (shift[d]!=0)
         {
            d_bdry_types[2*d+0] = PERIODIC;
            d_bdry_types[2*d+1] = PERIODIC;
         }
      }
   } 
   else 
   {
      TBOX_ERROR("PflotranJacobianLevelOperator" 
                 << " -- Required key `boundary_conditions'"
                 << " missing in input.");
   }   
#endif

   if (db->keyExists("stencilsize")) 
   {
      d_stencil_size = db->getInteger("stencilsize");      
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianLevelOperator" 
                 << " -- Required key `stencilsize'"
                 << " missing in input.");
   }

   if (db->keyExists("ndof")) 
   {
      d_ndof = db->getInteger("ndof");      
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianLevelOperator" 
                 << " -- Required key `ndof'"
                 << " missing in input.");
   }

}

void 
PflotranJacobianLevelOperator::MatSetValuesLocal(int patchNumber,
                                                 PetscInt nrow,const PetscInt irow[],
                                                 PetscInt ncol,const PetscInt icol[],
                                                 const PetscScalar y[],InsertMode addv)
{
   int ngx;
   int ngxy;

   for(int i=0;i<nrow; i++)
   {
      int currentRow = irow[i];
      k= int((currentRow)/ngxy) + 1;
      j= int(mod(currentRow, ngxy)/ngx) + 1;
      i= mod(mod(currentRow, ngxy),ngx) + 1;  

      for(int j=0;j<ncol; j++)
      {
         int currentCol = icol[i];
      }
   }
}

void 
PflotranJacobianLevelOperator::MatSetValuesBlockedLocal(int patchNumber,
                                                        PetscInt nrow,const PetscInt irow[],
                                                        PetscInt ncol,const PetscInt icol[],
                                                        const PetscScalar y[],InsertMode addv)
{
   for(int i=0;i<nrow; i++)
   {
      int currentRow = irow[i];
      for(int j=0;j<ncol; j++)
      {
         int currentCol = icol[i];
      }
   }

}

}
