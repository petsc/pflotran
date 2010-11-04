#include <vector>
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"
#endif
#ifndef included_solv_SAMRAIVectorReal
#include "SAMRAIVectorReal.h"
#include "SAMRAIVectorReal.C"
#endif
#include "tbox/Array.h"
#include "tbox/Array.C"
#include "tbox/List.h"
#include "tbox/List.C"
#include "CoarsenSchedule.h"
#include "CoarsenSchedule.C"
#include "BoundaryConditionStrategy.h"
#include "PflotranJacobianLevelOperator.h"

template class SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >;
template class  SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >;
template class SAMRAI::tbox::Array<SAMRAI::SAMRSolvers::PflotranJacobianLevelOperator*>;
template class SAMRAI::tbox::Pointer<BoundaryConditionStrategy>;
template class SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::xfer::CoarsenSchedule<NDIM> > >;
template class SAMRAI::tbox::Array< SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::xfer::CoarsenSchedule<NDIM> > > >;

#ifndef __APPLE__ 
template class std::_Rb_tree<int, std::pair<int const, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<3> > >, std::_Select1st<std::pair<int const, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<3> > > >, std::less<int>, std::allocator<std::pair<int const, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<3> > > > >;
#endif
