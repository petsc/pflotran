#include "PflotranJacobianMultilevelOperatorParameters.h"

namespace SAMRAI{

PflotranJacobianMultilevelOperatorParameters::PflotranJacobianMultilevelOperatorParameters()
{
}

PflotranJacobianMultilevelOperatorParameters::~PflotranJacobianMultilevelOperatorParameters()
{
}

PflotranJacobianMultilevelOperatorParameters::PflotranJacobianMultilevelOperatorParameters(const tbox::Pointer<tbox::Database> &database)
   :MultilevelOperatorParameters(database)
{
}

}
