#include "BoundaryConditionStrategy.h"

BoundaryConditionStrategy::BoundaryConditionStrategy(const int id)
{
  d_data_id = id;
}
 
BoundaryConditionStrategy::~BoundaryConditionStrategy()
{
}

void BoundaryConditionStrategy::setPhysicalBoundaryConditions(
                                  hier::Patch<NDIM>& patch,
                                  const double time,
                                  const hier::IntVector<NDIM>& ghost_width_to_fill)
{
}
