#include "PflotranApplicationParameters.h"

namespace SAMRAI{

PflotranApplicationParameters::PflotranApplicationParameters()
{
}

PflotranApplicationParameters::~PflotranApplicationParameters()
{
}

PflotranApplicationParameters::PflotranApplicationParameters(const tbox::Pointer<tbox::Database> &database)
   :ApplicationParameters(database)
{

}

}
