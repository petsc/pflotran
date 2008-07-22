#ifndef included_PflotranApplicationParameters
#define included_PflotranApplicationParameters

#include <string>

#include "tbox/Database.h"
#include "PatchHierarchy.h"
#include "ApplicationParameters.h"

namespace SAMRAI{

class PflotranApplicationParameters: public ApplicationParameters
{
public:
   PflotranApplicationParameters(const tbox::Pointer<tbox::Database> &database );
   ~PflotranApplicationParameters();

   // Computational grid where problem is solved.
   tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;

protected:
private:
   PflotranApplicationParameters();
};

}
#endif //  included_PflotranApplicationParameters
