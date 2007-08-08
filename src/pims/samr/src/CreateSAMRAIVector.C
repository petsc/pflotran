#include "PatchHierarchy.h"

#ifndef included_PatchLevel
#include "PatchLevel.h"
#endif

#ifndef included_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef included_PatchData
#include "PatchData.h"
#endif

#ifndef included_FaceVariable
#include "CellVariable.h"
#endif

#include "VariableDatabase.h"
#include "VariableContext.h"

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

   static int vec_i = 0;
   ostringstream ibuffer;
   ibuffer<<(long)vec_i;
   string object_str=ibuffer.str();

   string dataName("PIMS_variable_");
   dataName+=object_str;

   vec_i++;

   SAMRAI::tbox::Pointer< SAMRAI::pdat::CellVariable<NDIM,double> > pims_var;

   pims_var = new SAMRAI::pdat::CellVariable<NDIM,double>(dataName,1);

   int pims_var_id;

   SAMRAI::hier::VariableDatabase<NDIM>* variable_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

   const SAMRAI::tbox::Pointer< SAMRAI::hier::VariableContext > pims_cxt = variable_db->getContext("PIMS");

   SAMRAI::hier::IntVector<NDIM> nghosts;
   if(use_ghost)
   {
      nghosts = SAMRAI::hier::IntVector<NDIM>(1);
   }
   else
   {
      nghosts = SAMRAI::hier::IntVector<NDIM>(0);
   }

   pims_var_id = variable_db->registerVariableAndContext(pims_var,
                                                         pims_cxt,
                                                         nghosts);


   SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > p_hierarchy = (*hierarchy);
   const int nlevels = p_hierarchy->getNumberLevels();

   for(int ln=0; ln<p_hierarchy->getNumberLevels(); ln++)
   {
      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = p_hierarchy->getPatchLevel(ln);
      level->allocatePatchData(pims_var_id);
   }

   SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > samrai_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(dataName,
                                                                                                             p_hierarchy,
                                                                                                             0, nlevels-1);

   samrai_vec->addComponent(pims_var,
                            pims_var_id);

   *vec = SAMRAI::solv::PETSc_SAMRAIVectorReal<NDIM,double>::createPETScVector(samrai_vec);

}
#endif
