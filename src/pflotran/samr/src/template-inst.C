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
#include "CoarsenSchedule.h"
#include "CoarsenSchedule.C"
#include "BoundaryConditionStrategy.h"
#include "PflotranJacobianLevelOperator.h"

template class SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >;
template class  SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >;
template class SAMRAI::tbox::Array<SAMRAI::PflotranJacobianLevelOperator*>;
template class SAMRAI::tbox::Pointer<BoundaryConditionStrategy>;
template class SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::xfer::CoarsenSchedule<NDIM> > >;
template class SAMRAI::tbox::Array< SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::xfer::CoarsenSchedule<NDIM> > > >;
