#include "PatchHierarchy.h"
#ifdef __cplusplus
extern "C" 
void get_refinement_levels_(SAMRAI::hier::PatchHierarchy<NDIM> **hierarchy)
{
   int nlevels = (*hierarchy)->getNumberLevels();
}
#endif
