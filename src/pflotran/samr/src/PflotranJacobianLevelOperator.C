#include "PflotranJacobianLevelOperator.h"
#include "CartesianGridGeometry.h"
#include "CCellData.h"
#include "PETSc_SAMRAIVectorReal.h"

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
   tbox::Pointer< hier::Patch<NDIM> > patch = d_level->getPatch(patchNumber);
   // the pointer to the patch can be NULL if the patch is off processor
   if(!patch.isNull())
   {
      tbox::Pointer< pdat::CCellData<NDIM,double> > stencil = patch->getPatchData(d_stencil_id);
      double *stencilArray = stencil->getPointer();

      const hier::Box<NDIM> &box  = patch->getBox();
      hier::Box<NDIM> gbox = box;
      gbox.grow(1);

      const int ngx = gbox.numberCells(0);
      const int ngxy = ngx*(gbox.numberCells(1));

      const int nx = box.numberCells(0);
      const int nxy = nx*(box.numberCells(1));

      for(int i=0;i<nrow; i++)
      {
         int currentRow = irow[i];
         int kr= int((currentRow)/ngxy) - 1;
         int jr= int((currentRow%ngxy)/ngx) - 1;
         int ir= ((currentRow%ngxy)%ngx) - 1;  
         
         int offSet = kr*nxy+jr*nx+ir;

         for(int j=0;j<ncol; j++)
         {
            int currentCol = icol[j];
            int kl= int((currentCol)/ngxy) - 1;
            int jl= int((currentCol%ngxy)/ngx) - 1;
            int il= ((currentCol%ngxy)%ngx) - 1;  

            il = il-ir;
            jl = jl-jr;
            kl = kl-kr;
            
            if(d_stencil_size==7)
            {
               int pos;

               if((il==0)&(jl==0)&(kl==0))
               {
                  pos = 0;
               }
               else if((il==-1)&(jl==0)&(kl==0))
               {
                  pos = 1;
               }
               else if((il==1)&(jl==0)&(kl==0))
               {
                  pos = 2;
               }
               else if((il==0)&(jl==-1)&(kl==0))
               {
                  pos = 3;
               }
               else if((il==0)&(jl==1)&(kl==0))
               {
                  pos = 4;
               }
               else if((il==0)&(jl==0)&(kl==-1))
               {
                  pos = 5;
               }
               else if((il==0)&(jl==0)&(kl==1))
               {
                  pos = 6;
               }
               else
               {
                  pos=-1;
                  abort();
               }

               stencilArray[offSet+pos] = (addv==INSERT_VALUES)?y[j]:stencilArray[offSet+pos]+y[j];

            }
            else
            {

               tbox::pout << "ERROR:: PflotranJacobianLevelOperator:: invalid stencil size" << std::endl;
               abort();

            }
         }
      }
   }
}

void 
PflotranJacobianLevelOperator::MatSetValuesBlockedLocal(int patchNumber,
                                                        PetscInt nrow,const PetscInt irow[],
                                                        PetscInt ncol,const PetscInt icol[],
                                                        const PetscScalar y[],InsertMode addv)
{
   tbox::Pointer< hier::Patch<NDIM> > patch = d_level->getPatch(patchNumber);
   // the pointer to the patch can be NULL if the patch is off processor
   if(!patch.isNull())
   {
      tbox::Pointer< pdat::CCellData<NDIM,double> > stencil = patch->getPatchData(d_stencil_id);
      double *stencilArray = stencil->getPointer();

      const hier::Box<NDIM> &box  = patch->getBox();
      hier::Box<NDIM> gbox = box;
      gbox.grow(1);

      const int ngx = gbox.numberCells(0);
      const int ngxy = ngx*(gbox.numberCells(1));

      const int nx = box.numberCells(0);
      const int nxy = nx*(box.numberCells(1));

      for(int i=0;i<nrow; i++)
      {
         int currentRow = irow[i];
         int kr= int((currentRow)/ngxy) - 1;
         int jr= int((currentRow%ngxy)/ngx) - 1;
         int ir= ((currentRow%ngxy)%ngx) - 1;  
         
         int offSet = kr*nxy+jr*nx+ir;
         offSet = offSet*d_ndof;

         for(int j=0;j<ncol; j++)
         {
            int currentCol = icol[j];
            int kl= int((currentCol)/ngxy) - 1;
            int jl= int((currentCol%ngxy)/ngx) - 1;
            int il= ((currentCol%ngxy)%ngx) - 1;  

            il = il-ir;
            jl = jl-jr;
            kl = kl-kr;
            
            if(d_stencil_size==7)
            {
               int pos;

               if((il==0)&(jl==0)&(kl==0))
               {
                  pos = 0;
               }
               else if((il==-1)&(jl==0)&(kl==0))
               {
                  pos = 1;
               }
               else if((il==1)&(jl==0)&(kl==0))
               {
                  pos = 2;
               }
               else if((il==0)&(jl==-1)&(kl==0))
               {
                  pos = 3;
               }
               else if((il==0)&(jl==1)&(kl==0))
               {
                  pos = 4;
               }
               else if((il==0)&(jl==0)&(kl==-1))
               {
                  pos = 5;
               }
               else if((il==0)&(jl==0)&(kl==1))
               {
                  pos = 6;
               }
               else
               {
                  pos=-1;
                  abort();
               }

               pos = pos*d_ndof;

               if(addv==INSERT_VALUES)
               {
                  for(int m=0;m<d_ndof; m++)
                  {
                     for(int l=0;l<d_ndof; l++)
                     {
                        stencilArray[offSet+pos+m*d_ndof+l] = y[j*d_ndof+m*d_ndof+l];
                     }      
                  }
               }
               else
               {
                  for(int m=0;m<d_ndof; m++)
                  {
                     for(int l=0;l<d_ndof; l++)
                     {
                        stencilArray[offSet+pos+m*d_ndof+l] += y[j*d_ndof+m*d_ndof+l];
                     }      
                  }

               }
            }
            else
            {

               tbox::pout << "ERROR:: PflotranJacobianLevelOperator:: invalid stencil size" << std::endl;
               abort();

            }
         }
      }
   }
}

int
PflotranJacobianLevelOperator::MatMult(Vec x, Vec y )
{
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > xVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(x);
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > yVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(y);

   int src_id = xVec->getComponentDescriptorIndex(0);
   int dest_id = yVec->getComponentDescriptorIndex(0);   

   for (hier::PatchLevel<NDIM>::Iterator p(d_level); p; p++) 
   {
      tbox::Pointer<hier::Patch<NDIM> > patch = d_level->getPatch(p());
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!patch.isNull());
#endif
      tbox::Pointer< pdat::CCellData<NDIM, double > > stencil = patch->getPatchData(d_stencil_id);
   }
   
   return(0);
}

}
