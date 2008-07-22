#ifndef included_PflotranApplicationParameters
#define included_PflotranApplicationParameters

#include "ApplicationParameters.h"

namespace SAMRAI{

class PflotranApplicationParameters: public ApplicationParameters
{
public:
   PflotranApplicationParameters(const tbox::Pointer<tbox::Database> &database );
   ~PflotranApplicationParameters();
protected:
private:
   PflotranApplicationParameters();
};

}
#endif //  included_PflotranApplicationParameters
