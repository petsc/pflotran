
#ifndef included_boundary_condition_strategy
#define included_boundary_condition_strategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "RefinePatchStrategy.h"

using namespace SAMRAI;

/** \class BoundaryConditionStrategy
 *
 * BoundaryConditionStrategy is a concrete implementation of
 * xfer::RefinePatchStrategy<NDIM> that only sets physical boundary conditions.
 */
class BoundaryConditionStrategy :
 public xfer::RefinePatchStrategy<NDIM>
{
public:
   /**
    * The constructor caches the descriptor index of the variable
    * whose boundary conditions are set.
    */
   BoundaryConditionStrategy(const int id);

   /**
    * Virtual destructor.
    */
   virtual ~BoundaryConditionStrategy();

   /**
    * Set solution ghost cell values along physical boundaries.
    * \param patch
    *        Reference to a patch that touches the physical boundary.
    * \param time
    *        Simulation time.
    * \param ghost_width_to_fill
    *        Number of ghost cells in each coordinate direction that 
    *        must be filled.
    *
    * Function is overloaded from xfer::RefinePatchStrategy<NDIM>.
    */
   void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                      const double time,
                                      const hier::IntVector<NDIM>& ghost_width_to_fill);

   /**
    * Function for applying user-defined data processing prior
    * to refinement.  
    *
    * Function is overloaded from xfer::RefinePatchStrategy<NDIM>.  An empty
    * implementation is provided here.
    */
   void preprocessRefine(hier::Patch<NDIM>& fine,
                         const hier::Patch<NDIM>& coarse,
                         const hier::Box<NDIM>& fine_box,
                         const hier::IntVector<NDIM>& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   } 

   /**
    * Function for applying user-defined data processing after
    * refinement.  
    * 
    * Function is overloaded from xfer::RefinePatchStrategy<NDIM>.  An empty
    * implementation is provided here.
    */
   void postprocessRefine(hier::Patch<NDIM>& fine,
                          const hier::Patch<NDIM>& coarse,
                          const hier::Box<NDIM>& fine_box,
                          const hier::IntVector<NDIM>& ratio) 
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   } 

   /**
    * Return maximum stencil width needed for user-defined 
    * data interpolation operations.  Default is to return 
    * zero, assuming no user-defined operations provided.
    */
   hier::IntVector<NDIM> getRefineOpStencilWidth() const 
     {
       return(hier::IntVector<NDIM>(0));
     }

private:

   int d_data_id;
};

#endif
