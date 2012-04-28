#include "PflotranJacobianMultilevelOperatorParameters.h"

namespace SAMRAI{
  namespace SAMRSolvers{
    
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

}
