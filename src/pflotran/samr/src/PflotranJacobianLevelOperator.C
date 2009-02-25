#include "PflotranJacobianLevelOperator.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CCellData.h"
#include "CSideData.h"
#include "PETSc_SAMRAIVectorReal.h"
#include "fortran/3d/prototypes.h"

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
      d_flux = new pdat::CSideVariable<NDIM,double>(cellFlux,1);
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
   for (hier::PatchLevel<NDIM>::Iterator p(d_level); p; p++) 
   {
      tbox::Pointer<hier::Patch<NDIM> > patch = d_level->getPatch(p());
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!patch.isNull());
#endif

      tbox::Pointer< pdat::CCellData<NDIM, double > > stencil = patch->getPatchData(d_stencil_id);
      if(d_ndof==1)
      {
         tbox::Pointer< pdat::CCellData<NDIM, double > > u_data = patch->getPatchData(u_id[0]);
         assert(f_id[0]>=0);
         tbox::Pointer< pdat::CCellData<NDIM, double > > f_data = patch->getPatchData(f_id[0]);
         assert(r_id[0]>=0);
         tbox::Pointer< pdat::CCellData<NDIM, double > > r_data = patch->getPatchData(r_id[0]);

         hier::Box<NDIM> box = patch->getBox();

         const hier::Index<NDIM> ifirst = box.lower();
         const hier::Index<NDIM> ilast = box.upper();

         box = u_data->getGhostBox(); 
         const hier::Index<NDIM> ufirst = box.lower();
         const hier::Index<NDIM> ulast = box.upper();

         box = f_data->getGhostBox(); 
         const hier::Index<NDIM> ffirst = box.lower();
         const hier::Index<NDIM> flast = box.upper();

         box = r_data->getGhostBox(); 
         const hier::Index<NDIM> rfirst = box.lower();
         const hier::Index<NDIM> rlast = box.upper();

         samrapply7ptstencil3d_(ifirst(0),ifirst(1),ifirst(2),
                                ilast(0),ilast(1),ilast(2),
                                stencil->getPointer(),
                                ufirst(0),ufirst(1),ufirst(2),
                                ulast(0) ,ulast(1) ,ulast(2),
                                u_data->getPointer(),
                                ffirst(0),ffirst(1),ffirst(2),
                                flast(0) ,flast(1) ,flast(2),
                                f_data->getPointer(),
                                rfirst(0),rfirst(1),rfirst(2),
                                rlast(0) ,rlast(1) ,rlast(2),
                                r_data->getPointer());

      }
      else
      {
         abort();
      }
      
   }
}

void
PflotranJacobianLevelOperator::apply(const int flux_id,
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
#endif

   const int fidx=(f_idx==NULL)?0:f_idx[0];
   const int uidx=(u_idx==NULL)?0:u_idx[0];
   const int ridx=(r_idx==NULL)?0:r_idx[0];

   for (hier::PatchLevel<NDIM>::Iterator p(d_level); p; p++) 
   {
      tbox::Pointer<hier::Patch<NDIM> > patch = d_level->getPatch(p());
      const hier::Box<NDIM> interior(patch->getBox());
      const hier::Index<NDIM> ifirst  = interior.lower();
      const hier::Index<NDIM> ilast   = interior.upper();
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(patch->checkAllocated(f_id[0]));
      assert(patch->checkAllocated(flux_id));
      assert(patch->checkAllocated(d_srcsink_id));
      assert(patch->checkAllocated(r_id[0]));
#endif

      tbox::Pointer< pdat::CCellData<NDIM,double> > f_data = patch->getPatchData(f_id[0]);
      tbox::Pointer< pdat::CCellData<NDIM,double> > src_data = patch->getPatchData(d_srcsink_id);
      tbox::Pointer< pdat::CSideData<NDIM,double> > flux_data = patch->getPatchData(flux_id);
      tbox::Pointer< pdat::CCellData<NDIM,double> > r_data = patch->getPatchData(r_id[0]);
      tbox::Pointer< pdat::CCellData<NDIM,double> > u_data = patch->getPatchData(u_id[0]);
      tbox::Pointer< pdat::CCellData<NDIM,double> > stencil = patch->getPatchData(d_stencil_id);

      const hier::Index<NDIM> ufirst = u_data->getGhostBox().lower();
      const hier::Index<NDIM> ulast  = u_data->getGhostBox().upper();
      const hier::Index<NDIM> ffirst = f_data->getGhostBox().lower();
      const hier::Index<NDIM> flast  = f_data->getGhostBox().upper();
      const hier::Index<NDIM> sfirst = src_data->getGhostBox().lower();
      const hier::Index<NDIM> slast  = src_data->getGhostBox().upper();
      const hier::Index<NDIM> rfirst = r_data->getGhostBox().lower();
      const hier::Index<NDIM> rlast  = r_data->getGhostBox().upper();

      pflotranpcapply3d_(
         ifirst(0),ifirst(1),ifirst(2),ilast(0),ilast(1),ilast(2),
         ufirst(0),ufirst(1),ufirst(2),ulast(0),ulast(1),ulast(2),
         ffirst(0),ffirst(1),ffirst(2),flast(0),flast(1),flast(2),
         sfirst(0),sfirst(1),sfirst(2),slast(0),slast(1),slast(2),
         rfirst(0),rfirst(1),rfirst(2),rlast(0),rlast(1),rlast(2),
         d_stencil_size,
         d_ndof,
         stencil->getPointer(),
         u_data->getPointer(),
         flux_data->getPointer(0), flux_data->getPointer(1), flux_data->getPointer(2),
         f_data->getPointer(),
         src_data->getPointer(),
         r_data->getPointer());

   }
}

void
PflotranJacobianLevelOperator::setFlux(const int flux_id,
                                       const int *u_id,
                                       const int *u_idx)
{
   for (hier::PatchLevel<NDIM>::Iterator p(d_level); p; p++) 
   {
      tbox::Pointer<hier::Patch<NDIM> > patch = d_level->getPatch(p());
      
      const hier::Box<NDIM> box = patch->getBox();
      const hier::Index<NDIM> ifirst = box.lower();
      const hier::Index<NDIM> ilast = box.upper();
     
      tbox::Pointer< pdat::CCellData<NDIM,double> > stencil = patch->getPatchData(d_stencil_id);
      tbox::Pointer< pdat::CSideData<NDIM,double> > flux_data = patch->getPatchData(flux_id);

      tbox::Pointer< pdat::CCellData<NDIM,double> > u_data = patch->getPatchData(u_id[0]);
      const hier::Index<NDIM> gfirst = u_data->getGhostBox().lower();
      const hier::Index<NDIM> glast = u_data->getGhostBox().upper();
      
      assert(d_stencil_size==7);
      assert(d_ndof==1);

      pflotranpcflux3d_(ifirst(0),ifirst(1),ifirst(2),
                        ilast(0),ilast(1),ilast(2),
                        gfirst(0),gfirst(1),gfirst(2),
                        glast(0),glast(1),glast(2),
                        d_stencil_size,
                        d_ndof,
                        stencil->getPointer(),
                        u_data->getPointer(),
                        flux_data->getPointer(0),
                        flux_data->getPointer(1),
                        flux_data->getPointer(2));
   }
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
   tbox::Pointer< hier::PatchData<NDIM > > stencilBlock;
   stencilBlock.setNull();

   if(i==0&&j==0&&k==0)
   {
      tbox::Pointer<hier::Patch<NDIM> > patch = d_level->getPatch(p);
      if(!patch.isNull())
      {
         stencilBlock=patch->getPatchData(d_stencil_id);
      }
   }

   return stencilBlock;
}

const int
PflotranJacobianLevelOperator::getStencilSize(const int i,
                                              const int j,
                                              const int k)
{
   return d_stencil_size;
}

std::vector<int>
PflotranJacobianLevelOperator::getStencilOffsets(const int i, 
                                                 const int j, 
                                                 const int k)  
{

   std::vector<int> offsets;
   if(d_ndof==1)
   {
   if((NDIM==1)&&i==0)
   {
      offsets.resize(d_stencil_size*NDIM);
      offsets[0]=0;
      offsets[1]=-1;
      offsets[2]=1;
   }
   else if((NDIM==2)&&i==0&&j==0)
   {
      offsets.resize(d_stencil_size*NDIM);
      offsets[0]=0;
      offsets[1]=0;
      offsets[2]=-1;
      offsets[3]=0;
      offsets[4]=1;
      offsets[5]=0;
      offsets[6]=0;
      offsets[7]=-1;
      offsets[8]=0;
      offsets[9]=1;

   }
   else if((NDIM==3)&&i==0&&j==0&&k==0)
   {
      offsets.resize(d_stencil_size*NDIM);
      offsets[0]=0;
      offsets[1]=0;
      offsets[2]=0;
      offsets[3]=-1;
      offsets[4]=0;
      offsets[5]=0;
      offsets[6]=1;
      offsets[7]=0;
      offsets[8]=0;
      offsets[9]=0;
      offsets[10]=-1;
      offsets[11]=0;
      offsets[12]=0;
      offsets[13]=1;
      offsets[14]=0;
      offsets[15]=0;
      offsets[16]=0;
      offsets[17]=-1;
      offsets[18]=0;
      offsets[19]=0;
      offsets[20]=1;
   }
   }
   else
   {
      abort();
   }

   return offsets;
}

void 
PflotranJacobianLevelOperator::getFromInput(const tbox::Pointer<tbox::Database> &db)
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

   if (db->keyExists("object_id")) 
   {
      d_object_id = db->getInteger("object_id");      
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianLevelOperator" 
                 << " -- Required key `object_id'"
                 << " missing in input.");
   }

   if (db->keyExists("srcsink_id")) 
   {
      d_srcsink_id = db->getInteger("srcsink_id");      
   } 
   else 
   {
      TBOX_ERROR( "PflotranJacobianLevelOperator" 
                 << " -- Required key `srcsink_id'"
                 << " missing in input.");
   }

}

void 
PflotranJacobianLevelOperator::MatSetValuesLocal(int patchNumber,
                                                 PetscInt nrow,const PetscInt irow[],
                                                 PetscInt ncol,const PetscInt icol[],
                                                 const PetscScalar y[],InsertMode addv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(nrow==1);
   assert(ncol==1);
#endif
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
         
         int offSet = (kr*nxy+jr*nx+ir)*d_stencil_size;

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
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(nrow==1);
   assert(ncol==1);
#endif
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
         offSet = offSet*d_stencil_size*d_ndof*d_ndof;

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

               pos = pos*d_ndof*d_ndof;

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
   int dst_id = yVec->getComponentDescriptorIndex(0);   

   for (hier::PatchLevel<NDIM>::Iterator p(d_level); p; p++) 
   {
      tbox::Pointer<hier::Patch<NDIM> > patch = d_level->getPatch(p());
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!patch.isNull());
#endif
      const hier::Box<NDIM> &box = patch->getBox();

      tbox::Pointer< pdat::CCellData<NDIM, double > > stencil = patch->getPatchData(d_stencil_id);
      // for now we will use cell data instead of CCell data
      tbox::Pointer< pdat::CCellData<NDIM, double > > src = patch->getPatchData(src_id);
      tbox::Pointer< pdat::CCellData<NDIM, double > > dst = patch->getPatchData(dst_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!stencil.isNull());
      assert(!src.isNull());
      assert(!dst.isNull());
#endif

      const hier::Index<NDIM> ifirst = box.lower();
      const hier::Index<NDIM> ilast = box.upper();
      
      hier::IntVector<NDIM> ghost_cell_width = src->getGhostCellWidth();
      const int sgcw=ghost_cell_width(0);
      ghost_cell_width = dst->getGhostCellWidth();
      const int dgcw=ghost_cell_width(0);

#ifdef DEBUG_CHECK_ASSERTIONS
        assert(sgcw>=0);
        assert(dgcw>=0);
        assert(d_stencil_size>0);
        assert(d_ndof>=1);
#endif
        samrccellmatmult3d_( ifirst(0),ifirst(1),ifirst(2),
                            ilast(0),ilast(1),ilast(2),
                            d_stencil_size,
                            d_ndof,
                            stencil->getPointer(),
                            sgcw,
                            src->getPointer(),
                            dgcw,
                            dst->getPointer());
   }
   
   return(0);
}

void
PflotranJacobianLevelOperator::setSourceValueOnPatch(SAMRAI::hier::Patch<NDIM> **patch, int *index, double *val)
{
   tbox::Pointer< pdat::CCellData<NDIM,double> > src_data = (*patch)->getPatchData(d_srcsink_id);
   double *srcArray = src_data->getPointer();

   int currentRow = (*index);
   const hier::Box<NDIM> &box  = (*patch)->getBox();
   hier::Box<NDIM> gbox = box;
   gbox.grow(1);
   
   const int ngx = gbox.numberCells(0);
   const int ngxy = ngx*(gbox.numberCells(1));
   
   const int nx = box.numberCells(0);
   const int nxy = nx*(box.numberCells(1));
   int kr= int((currentRow)/ngxy) - 1;
   int jr= int((currentRow%ngxy)/ngx) - 1;
   int ir= ((currentRow%ngxy)%ngx) - 1;  

   int offSet = (kr*nxy+jr*nx+ir)*d_ndof;
   srcArray[offSet] += (*val);
}

int
PflotranJacobianLevelOperator::getVariableIndex(std::string &name, 
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

}
