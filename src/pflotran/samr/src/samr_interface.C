#include "tbox/Pointer.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PETSc_SAMRAIVectorReal.h"
#include "CCellData.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "PflotranApplicationStrategy.h"
#include "PflotranJacobianMultilevelOperator.h"
#include "PflotranJacobianMultilevelOperatorParameters.h"

extern "C" {
#include "petscvec.h"
#include "petscmat.h"
void  cf90bridge_(void *, int*, void *);

}

#ifdef __cplusplus
extern "C" {
#endif

int hierarchy_number_levels_(SAMRAI::PflotranApplicationStrategy **application_strategy)
{
   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = (*application_strategy)->getHierarchy();
   return hierarchy->getNumberOfLevels();
}

int level_number_patches_(SAMRAI::PflotranApplicationStrategy **application_strategy, int *ln)
{
   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = (*application_strategy)->getHierarchy();
   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(*ln);
   return level->getNumberOfPatches();
}

bool is_local_patch_(SAMRAI::PflotranApplicationStrategy **application_strategy, int *ln, int *pn)
{
   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = (*application_strategy)->getHierarchy();
   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(*ln);
   return(level->getProcessorMapping().isMappingLocal(*pn));
}

void *hierarchy_get_patch_(SAMRAI::PflotranApplicationStrategy **application_strategy, int *ln, int *pn)
{
   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = (*application_strategy)->getHierarchy();
   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(*ln);
   SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(*pn);
   return (void *)(patch.getPointer());
}

void samr_physical_dimensions_(SAMRAI::PflotranApplicationStrategy **application_strategy, 
                               int *nx, int *ny, int *nz)
{
   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = (*application_strategy)->getHierarchy();
   SAMRAI::hier::BoxArray<NDIM> physicalDomain = hierarchy->getGridGeometry()->getPhysicalDomain();

   // assume one box representing the domain for now
   SAMRAI::hier::Box< NDIM > physicalBox = physicalDomain.getBox(0);
   (*nx) = physicalBox.numberCells(0);
   (*ny) = physicalBox.numberCells(1);
   (*nz) = physicalBox.numberCells(2);
}

void samr_get_origin_(SAMRAI::PflotranApplicationStrategy **application_strategy,
                     double *x0,
                     double *y0,
                     double *z0)
{
   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = (*application_strategy)->getHierarchy();
   SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry = hierarchy->getGridGeometry();
   const double* xlo = grid_geometry->getXLower();
   (*x0) = xlo[0];
   (*y0) = xlo[1];
   (*z0) = xlo[2];
}


void samr_get_upper_corner_(SAMRAI::PflotranApplicationStrategy **application_strategy,
                     double *x0,
                     double *y0,
                     double *z0)
{
   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = (*application_strategy)->getHierarchy();
   SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry = hierarchy->getGridGeometry();
   const double* xup = grid_geometry->getXUpper();
   (*x0) = xup[0];
   (*y0) = xup[1];
   (*z0) = xup[2];
}

void samr_patch_get_origin_(SAMRAI::hier::Patch<NDIM> **patch,
                            double *x0,
                            double *y0,
                            double *z0)
{
   SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geometry = (*patch)->getPatchGeometry();
   const double* xlo = patch_geometry->getXLower();
   (*x0) = xlo[0];
   (*y0) = xlo[1];
   (*z0) = xlo[2];
}

void samr_patch_get_corners_(SAMRAI::hier::Patch<NDIM> **patch, 
                       int *nxs, int *nys, int *nzs,
                       int *nlx, int *nly, int *nlz)
{
   SAMRAI::hier::Index<NDIM> lower = (*patch)->getBox().lower();
   (*nxs) = lower(0); 
   (*nys) = lower(1); 
   (*nzs) = lower(2);

   SAMRAI::hier::IntVector<NDIM> length = (*patch)->getBox().numberCells();
   (*nlx) = length(0);
   (*nly) = length(1);
   (*nlz) = length(2);
   
}

// for now we hard code a ghost cell width of 1
void samr_patch_get_ghostcorners_(SAMRAI::hier::Patch<NDIM> **patch, 
                       int *nxs, int *nys, int *nzs,
                       int *nlx, int *nly, int *nlz)
{
   SAMRAI::hier::Index<NDIM> lower = (*patch)->getBox().lower();
   (*nxs) = lower(0)-1; 
   (*nys) = lower(1)-1; 
   (*nzs) = lower(2)-1;

   SAMRAI::hier::IntVector<NDIM> length = (*patch)->getBox().numberCells();
   (*nlx) = length(0)+2;
   (*nly) = length(1)+2;
   (*nlz) = length(2)+2;
   
}

void samr_vecgetarrayf90_(SAMRAI::hier::Patch<NDIM> **patch, 
                          Vec *petscVec,
                          void **f90wrap)

{
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > sVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(*petscVec);
   SAMRAI::tbox::Pointer< SAMRAI::pdat::CCellData<NDIM, double> > pData = sVec->getComponentPatchData(0, *(*patch));
   int depth = pData->getDepth();

   int len = pData->getGhostBox().size();
   len = len*depth;

   void *p_data_ptr = pData->getPointer(0);

   cf90bridge_(p_data_ptr, &len, *f90wrap);
   
}

int samr_patch_at_bc_(SAMRAI::hier::Patch<NDIM> **patch, 
                      int *axis, int *side)
{
   int istouching = (int)(*patch)->getPatchGeometry()->getTouchesRegularBoundary(*axis, *side);
   return istouching;
}


void samr_patch_get_spacing_(SAMRAI::hier::Patch<NDIM> **patch, 
                             double *dx, double *dy, double *dz)
{
   SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geometry = (*patch)->getPatchGeometry();
   const double* ds = patch_geometry->getDx();
   *dx = ds[0];
   *dy = ds[1];
   *dz = ds[2];
}

void samrvecrestorearray_(SAMRAI::hier::Patch<NDIM> **patch, 
                         Vec *rvec,
                         void **f90wrap,
                         int *ierr)
{
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > sVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(*rvec);

}

void samrcreatematrix_(SAMRAI::PflotranApplicationStrategy **application_strategy, 
                       int *ndof,
                       int *stencilSize,
                       int *flowortransport,
                       Mat *pMatrix)    
{
   tbox::Pointer<tbox::Database> application_db = (*application_strategy)->getDatabase();
   
   tbox::Pointer<tbox::Database> operator_db = application_db->getDatabase("PflotranMultilevelOperator");

   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = (*application_strategy)->getHierarchy();

   operator_db->putInteger("ndof", *ndof);
   operator_db->putInteger("stencilsize", *stencilSize);

   SAMRAI::PflotranJacobianMultilevelOperatorParameters *parameters = new  SAMRAI::PflotranJacobianMultilevelOperatorParameters(operator_db);
   parameters->d_hierarchy = hierarchy;
   parameters->d_pMatrix = pMatrix;
   parameters->d_cf_interpolant = (*application_strategy)->getRefinementBoundaryInterpolant();
   parameters->d_set_boundary_ghosts.setNull();

   SAMRAI::PflotranJacobianMultilevelOperator *pJacobian = new SAMRAI::PflotranJacobianMultilevelOperator((SAMRAI::MultilevelOperatorParameters *)parameters);
   if(*flowortransport==0)
   {
      (*application_strategy)->setFlowJacobianMatrix(pJacobian);
   }
   else
   {
      (*application_strategy)->setTransportJacobianMatrix(pJacobian);
   }

}

void samrglobaltolocal_(SAMRAI::PflotranApplicationStrategy **application_strategy, 
                        Vec *gvec, 
                        Vec *lvec, 
                        int *ierr)
{
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > globalVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(*gvec);

   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > localVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(*lvec);

   (*application_strategy)->interpolateGlobalToLocalVector(globalVec, localVec, *ierr);
   
}

void samrlocaltolocal_(SAMRAI::PflotranApplicationStrategy **application_strategy, 
                       Vec *svec, 
                       Vec *dvec, 
                       int *ierr)
{
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > srcVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(*svec);

   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM, double > > dstVec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM, double>::getSAMRAIVector(*dvec);

   (*application_strategy)->interpolateLocalToLocalVector(srcVec, dstVec, *ierr);
   
}

void 
samrsetcurrentjacobianpatch_( Mat *mat, SAMRAI::hier::Patch<NDIM> **patch)
{
   SAMRAI::PflotranJacobianMultilevelOperator *pJacobian = NULL;

   MatShellGetContext(*mat, (void **)&pJacobian);
   
   pJacobian->setCurrentPatch(*patch);
}


void create_samrai_vec_(SAMRAI::PflotranApplicationStrategy **application_strategy,
                        int &dof, 
                        bool &use_ghost,
                        Vec *vec)
{

   (*application_strategy)->createVector(dof, use_ghost, vec);
}

void
samrpetscobjectstateincrease_(Vec *vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(vec == (Vec)NULL));
#endif
   int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*vec)); PETSC_SAMRAI_ERROR(ierr);
}

void
samr_mpi_min_(double *x, double *y, double *z)
{
   double vals[3];

   vals[0] = *x;
   vals[1] = *y;
   vals[2] = *z;

   SAMRAI::tbox::SAMRAI_MPI::minReduction(vals);

   *x = vals[0];
   *y = vals[1];
   *z = vals[2];
}

void
samr_mpi_max_(double *x, double *y, double *z)
{
   double vals[3];

   vals[0] = *x;
   vals[1] = *y;
   vals[2] = *z;

   SAMRAI::tbox::SAMRAI_MPI::maxReduction(vals);

   *x = vals[0];
   *y = vals[1];
   *z = vals[2];
}

void samrbarrier_(void)
{
   SAMRAI::tbox::SAMRAI_MPI::barrier();
}

#ifdef __cplusplus
}
#endif
