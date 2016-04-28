
#ifndef included_BogusTagAndInitStrategy
#define included_BogusTagAndInitStrategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "StandardTagAndInitStrategy.h"

using namespace SAMRAI;

/** \class BogusTagAndInitStrategy
 *
 * BogusTagAndInitStrategy is a concrete implementation of
 * mesh::StandardTagAndInitStrategy<NDIM> that does nothing (hence, it is bogus).
 */

class BogusTagAndInitStrategy :
   public mesh::StandardTagAndInitStrategy<NDIM>
{
public:

   /**
    * Empty constructor.
    */
   BogusTagAndInitStrategy();

   /**
    * Empty destructor.
    */
  ~BogusTagAndInitStrategy();

   /**
    * Empty method to initialize data on a level.
    */
   void initializeLevelData(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
                            const int level_number,
                            const double time,
                            const bool can_be_refined,
                            const bool initial_time,
                            const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level = NULL,
                            const bool allocate_data = true );

   /**
    * Empty method to reset data structures that must be updated when
    * the grid is changed.
    */
   void resetHierarchyConfiguration(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
                                    const int coarsest_level,
                                    const int finest_level);

   /**
    * Empty method to identify locations where the grid should be
    * refined.
    */
   void applyGradientDetector(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
                              const int level_number,
                              const double time,
                              const int tag_index,
                              const bool initial_time,
                              const bool uses_richardson_extrapolation_too);
};

#endif
