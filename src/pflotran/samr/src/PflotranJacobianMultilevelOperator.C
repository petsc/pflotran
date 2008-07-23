#include "PflotranJacobianMultilevelOperator.h"

namespace SAMRAI{

PflotranJacobianMultilevelOperator::PflotranJacobianMultilevelOperator()
{
}

PflotranJacobianMultilevelOperator::~PflotranJacobianMultilevelOperator()
{
}

PflotranJacobianMultilevelOperator::PflotranJacobianMultilevelOperator(MultilevelOperatorParameters *parameters)
{
}

void
PflotranJacobianMultilevelOperator::apply(const int ln,
                                          const int *f_id,
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
PflotranJacobianMultilevelOperator::apply(const int coarse_ln,
                                          const int fine_ln,
                                          const int *f_id,
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
PflotranJacobianMultilevelOperator::applyBoundaryCondition(const int ln,
                                                           const int *var_id,
                                                           const int *var_idx,
                                                           const int *var_components,
                                                           const int number_of_variables,
                                                           const bool reset_ghost_values)
{
}

void
PflotranJacobianMultilevelOperator::initializeInternalVariableData(void)
{
}

void
PflotranJacobianMultilevelOperator::getFromInput(tbox::Pointer<tbox::Database> db)
{
}

void
PflotranJacobianMultilevelOperator::setFlux(const int coarse_ln,
                                            const int fine_ln,
                                            const int *u_id,
                                            const int *u_idx)
{
}

void
PflotranJacobianMultilevelOperator::setupTransferSchedules(void)
{
}

}

