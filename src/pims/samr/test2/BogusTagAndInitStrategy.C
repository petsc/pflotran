
#include "BogusTagAndInitStrategy.h"

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

BogusTagAndInitStrategy::BogusTagAndInitStrategy()
{
}

BogusTagAndInitStrategy::~BogusTagAndInitStrategy()
{
}

void BogusTagAndInitStrategy::initializeLevelData(
   const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy, 
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level,
   const bool allocate_data )
{
  (void) hierarchy;
  (void) old_level;
}

void BogusTagAndInitStrategy::resetHierarchyConfiguration(
   const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy, 
   const int coarsest_level, 
   const int finest_level)
{
  (void) hierarchy;
}

void BogusTagAndInitStrategy::applyGradientDetector(
   const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy, 
   const int level_number, 
   const double time, 
   const int tag_index, 
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
  (void) hierarchy;
}
