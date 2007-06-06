#include "tbox/Pointer.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"

#ifdef __cplusplus
extern "C" {
int hierarchy_number_levels_(SAMRAI::hier::PatchHierarchy<NDIM> **hierarchy)
{
   return (*hierarchy)->getNumberLevels();
}

int level_number_patches_(SAMRAI::hier::PatchHierarchy<NDIM> **hierarchy, int *ln)
{
   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = (*hierarchy)->getPatchLevel(*ln);
   return level->getNumberOfPatches();
}

bool is_local_patch_(SAMRAI::hier::PatchHierarchy<NDIM> **hierarchy, int *ln, int *pn)
{
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = (*hierarchy)->getPatchLevel(*ln);
    return(level->getProcessorMapping().isMappingLocal(*pn));
}

void *hierarchy_get_patch_(SAMRAI::hier::PatchHierarchy<NDIM> **hierarchy, int *ln, int *pn)
{
   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = (*hierarchy)->getPatchLevel(*ln);
   SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(*pn);
   return (void *)(patch.getPointer());
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


}
#endif
