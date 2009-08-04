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
#include "CoarsenSchedule.h"
#include "CoarsenSchedule.C"
#include "BoundaryConditionStrategy.h"
#include "PflotranJacobianLevelOperator.h"
#ifndef __PGI
#include "VisInterlacedDataStrategy.h"
#include "VisInterlacedDataStrategy.C"
#endif

template class SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >;
template class  SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >;
template class SAMRAI::tbox::Array<SAMRAI::PflotranJacobianLevelOperator*>;
template class SAMRAI::tbox::Pointer<BoundaryConditionStrategy>;
template class SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::xfer::CoarsenSchedule<NDIM> > >;
template class SAMRAI::tbox::Array< SAMRAI::tbox::Array< SAMRAI::tbox::Pointer< SAMRAI::xfer::CoarsenSchedule<NDIM> > > >;
#ifndef __PGI
template void std::vector<SAMRAI::appu::VisInterlacedDataStrategy<NDIM>*, std::allocator<SAMRAI::appu::VisInterlacedDataStrategy<NDIM>*> >::_M_insert_aux(__gnu_cxx::__normal_iterator<SAMRAI::appu::VisInterlacedDataStrategy<NDIM>**, std::vector<SAMRAI::appu::VisInterlacedDataStrategy<NDIM>*, std::allocator<SAMRAI::appu::VisInterlacedDataStrategy<NDIM>*> > >, SAMRAI::appu::VisInterlacedDataStrategy<NDIM>* const&);
#endif
