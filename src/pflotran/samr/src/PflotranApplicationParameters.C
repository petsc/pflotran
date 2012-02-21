#include "PflotranApplicationParameters.h"

namespace SAMRAI{

PflotranApplicationParameters::PflotranApplicationParameters()
{
   d_hierarchy.setNull();
}

PflotranApplicationParameters::~PflotranApplicationParameters()
{
}

PflotranApplicationParameters::PflotranApplicationParameters(const tbox::Pointer<tbox::Database> &database)
   :ApplicationParameters(database)
{
   d_hierarchy.setNull();
}

}
