#include "PatchHierarchy.h"
#include "PETSc_SAMRAIVectorReal.h"
#include "SAMRAIVectorReal.h"
#include <string>

#ifdef __cplusplus
extern "C" 
void create_samrai_vec_(SAMRAI::hier::PatchHierarchy<NDIM> **hierarchy, 
                        int &dof, 
                        bool &use_ghost,
                        Vec *vec)
{
   string name="";
   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > p_hierarchy = (*hierarchy);
   const int nlevels = p_hierarchy->getNumberLevels();
   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > samrai_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(name,
                                                                                                             p_hierarchy,
                                                                                                             0, nlevels-1);
   *vec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM,double>::createPETScVector(samrai_vec);

}
#endif
