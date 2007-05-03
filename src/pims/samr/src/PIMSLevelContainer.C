#include "PIMSLevelContainer.h"
#include "fc_interface.h"

PIMSLevelContainer::PIMSLevelContainer()
{
   
}

PIMSLevelContainer::PIMSLevelContainer(const int n)
{
   patch_data.resize(n);

   for(int i=0; i<n;i++)
   {
      f_create_local_patch_data_(&(patch_data[i]));
   }
}

PIMSLevelContainer::~PIMSLevelContainer()
{

}
