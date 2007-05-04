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

}
#endif
