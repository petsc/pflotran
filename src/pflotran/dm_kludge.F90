module DM_Kludge_module
  
  use Unstructured_Grid_Aux_module, only : ugdm_type
  
  implicit none

  private
 
#include "definitions.h"

#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmshell.h90"

  type, public :: dm_ptr_type
    DM :: dm  ! PETSc DM
    type(ugdm_type), pointer :: ugdm
      ! Unstructured grid "private" dm.  This gets wrapped in a PETSc DM via 
      ! DMShell routines.
  end type dm_ptr_type

end module DM_Kludge_module
