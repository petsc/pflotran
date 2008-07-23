#include "PflotranJacobianLevelOperator.h"

namespace SAMRAI{

PflotranJacobianLevelOperator::PflotranJacobianLevelOperator(LevelOperatorParameters *parameters)
{
}

PflotranJacobianLevelOperator::PflotranJacobianLevelOperator()
{
}

PflotranJacobianLevelOperator::~PflotranJacobianLevelOperator()
{
}

void 
PflotranJacobianLevelOperator::apply(const int *f_id,
                                     const int *u_id, 
                                     const int *r_id,
                                     const int *f_idx,
                                     const int *u_idx,
                                     const int *r_idx,
                                     const double a,
                                     const double b)
{
}

void
PflotranJacobianLevelOperator::applyBoundaryCondition(const int *var_id,
                               const int *var_idx,
                               const int *var_components,
                               const int number_of_variables)
{
}

const int 
PflotranJacobianLevelOperator::getStencilType(const int i,
                                              const int j,
                                              const int k)
{
   return(0);
}

tbox::Pointer< hier::PatchData<NDIM > >
PflotranJacobianLevelOperator::getStencilBlock(const int p, 
                                               const int i, 
                                               const int j,
                                               const int k)
{
}

const int
PflotranJacobianLevelOperator::getStencilSize(const int i,
                                              const int j,
                                              const int k)
{
   return (0);
}

}
