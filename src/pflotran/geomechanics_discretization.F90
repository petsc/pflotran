#ifdef GEOMECH
module Geomechanics_Discretization_module

  use Geomech_Grid_module
  use Geomech_Grid_Aux_module
  
  implicit none

  private
 
#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmshell.h90"

  type, public :: gmdm_ptr_type
    DM :: dm  ! PETSc DM
    type(gmdm_type), pointer :: gmdm
  end type gmdm_ptr_type

  type, public :: geomech_discretization_type
    PetscInt :: itype                          ! type of discretization (e.g. structured, unstructured, etc.)
    character(len=MAXWORDLENGTH) :: ctype      ! name of discretization
    PetscReal :: origin(3)                     ! origin of global domain
    type(geomech_grid_type), pointer :: grid   ! pointer to a grid object
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: dm_index_to_ndof(3)            ! mapping between a dm_ptr to the number of degrees of freedom
    type(gmdm_ptr_type), pointer :: dm_1dof
    type(gmdm_ptr_type), pointer :: dm_ngeodof    
  end type geomech_discretization_type

  public :: GeomechDiscretizationCreate, &
            GeomechDiscretizationDestroy, &
            GeomechDiscretizationCreateVector, &
            GeomechDiscretizationDuplicateVector, &         
!            GeomechDiscretizationCreateJacobian, &
            GeomechDiscretizationGlobalToLocal, &
            GeomechDiscretizationLocalToGlobal, &
            GeomechDiscretizationLocalToLocal, &
            GeomechDiscretizationGlobalToNatural, &
            GeomechDiscretizationNaturalToGlobal, &
            GeomechDiscretizationGlobalToLocalBegin, &
            GeomechDiscretizationGlobalToLocalEnd, &
            GeomechDiscretizationLocalToLocalBegin, &
            GeomechDiscretizationLocalToLocalEnd, &
            GeomechDiscretizGlobalToNaturalBegin, &
            GeomechDiscretizGlobalToNaturalEnd, &
            GeomechDiscretizNaturalToGlobalBegin, &
            GeomechDiscretizNaturalToGlobalEnd, &
            GeomechDiscretizationCreateDMs,&
            GeomechDiscretizationGetDMPtrFromIndex, &
            GeomechDiscretAOApplicationToPetsc
            
contains

! ************************************************************************** !
!
! GeomechDiscretizationCreate: Creates a geomechanics discretization
! author: Satish Karra, LANL
! date: 05/23/2013
!
! ************************************************************************** !
function GeomechDiscretizationCreate()

  implicit none
  
  type(geomech_discretization_type), pointer :: GeomechDiscretizationCreate
  type(geomech_discretization_type), pointer :: discretization
  
  allocate(discretization)
  discretization%ctype = ''
  discretization%itype = 0
  discretization%origin = 0.d0
  discretization%filename = ''

  ! nullify DM pointers
  allocate(discretization%dm_1dof)
  allocate(discretization%dm_ngeodof)
  discretization%dm_1dof%dm = 0
  discretization%dm_ngeodof%dm = 0
  nullify(discretization%dm_1dof%gmdm)
  nullify(discretization%dm_ngeodof%gmdm)  
  nullify(discretization%grid)
  
  GeomechDiscretizationCreate => discretization

end function GeomechDiscretizationCreate

! ************************************************************************** !
!
! GeomechDiscretizationCreateDMs: creates distributed, parallel meshes/grids
! If there are multiple degrees of freedom per grid cell, this will call 
! GeomechDiscretizationCreateDM() multiple times to create the DMs corresponding 
! to one degree of freedom grid cell and those corresponding to multiple 
! degrees of freedom per cell for geomechanics.
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationCreateDMs(discretization,option)
      
  use Option_module    
      
  implicit none
  
  type(geomech_discretization_type)            :: discretization
  type(option_type)                            :: option
      
  PetscInt                                     :: ndof
  PetscErrorCode                               :: ierr
  type(geomech_grid_type), pointer             :: geomech_grid

  !-----------------------------------------------------------------------
  ! Generate the DM objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call GeomechDiscretizationCreateDM(discretization,discretization%dm_1dof, &
                                     ndof,option)
  
  if (option%ngeomechdof > 0) then
    ndof = option%ngeomechdof
    call GeomechDiscretizationCreateDM(discretization,discretization%dm_ngeodof, &
                                       ndof,option)
  endif


end subroutine GeomechDiscretizationCreateDMs

! ************************************************************************** !
!
! GeomechDiscretizationCreateDM: creates a distributed, parallel mesh/grid
! for geomechanics
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationCreateDM(discretization,dm_ptr,ndof,option)

  use Option_module
  
  implicit none
  
  type(geomech_discretization_type)                :: discretization
  type(gmdm_ptr_type), pointer                     :: dm_ptr
  type(option_type)                                :: option
  PetscInt                                         :: ndof
  PetscErrorCode                                   :: ierr

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      option%io_buffer = &
        'Geomechanics currently works only with unstructured grid.'
      call printErrMsg(option)
    case(UNSTRUCTURED_GRID)
      call GMCreateGMDM(discretization%grid, &
                        dm_ptr%gmdm,ndof,option)
      call DMShellCreate(option%mycomm,dm_ptr%dm,ierr)
      call DMShellSetGlobalToLocalVecScatter(dm_ptr%dm,dm_ptr%gmdm%scatter_gtol,ierr)
  end select

end subroutine GeomechDiscretizationCreateDM

! ************************************************************************** !
!
! GeomechDiscretizationCreateVector: Creates a PETSc vector for the nodes
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationCreateVector(discretization,dm_index,vector, &
                                             vector_type,option)
  use Option_module                                      

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(option_type)                             :: option
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  PetscInt                                      :: dm_index
  Vec                                           :: vector
  PetscInt                                      :: vector_type
  PetscInt                                      :: ndof
  PetscErrorCode                                :: ierr
  
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)

  call GMGridDMCreateVector(discretization%grid,dm_ptr%gmdm,vector, &
                            vector_type,option)
                            
  call VecSet(vector,0.d0,ierr)
  
end subroutine GeomechDiscretizationCreateVector

! ************************************************************************** !
!
! GeomechDiscretizationDuplicateVector: Duplicates a Petsc vector
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationDuplicateVector(discretization,vector1,vector2)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  Vec                                           :: vector1
  Vec                                           :: vector2
  PetscErrorCode                                :: ierr
  
  call VecDuplicate(vector1,vector2,ierr)
  call VecCopy(vector1,vector2,ierr)
  
end subroutine GeomechDiscretizationDuplicateVector

! ************************************************************************** !
!
! GeomechDiscretizationGetDMPtrFromIndex: Returns the integer pointer for 
! the Geomech DM referenced
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
function GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: GeomechDiscretizationGetDMPtrFromIndex
  PetscInt                                      :: dm_index
  
  select case (dm_index)
    case(ONEDOF)
      GeomechDiscretizationGetDMPtrFromIndex => discretization%dm_1dof
    case(NGEODOF)
      GeomechDiscretizationGetDMPtrFromIndex => discretization%dm_ngeodof
  end select  
  
end function GeomechDiscretizationGetDMPtrFromIndex

! ************************************************************************** !
!
! GeomechDiscretizationGlobalToLocal: Performs global to local communication
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationGlobalToLocal(discretization,global_vec, &
                                              local_vec,dm_index)

  implicit none

  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: global_vec
  Vec                                           :: local_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  call DMGlobalToLocalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec,ierr)
  call DMGlobalToLocalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec,ierr)
  
end subroutine GeomechDiscretizationGlobalToLocal

! ************************************************************************** !
!
! GeomechDiscretizationLocalToGlobal: Performs local to global communication
! with DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationLocalToGlobal(discretization,local_vec, &
                                              global_vec,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: local_vec
  Vec                                           :: global_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_ltog,local_vec,global_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_ltog,local_vec,global_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  
end subroutine GeomechDiscretizationLocalToGlobal

! ************************************************************************** !
!
! GeomechDiscretizationLocalToLocal: Performs local to local communication 
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationLocalToLocal(discretization,local_vec1, &
                                             local_vec2,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: local_vec1
  Vec                                           :: local_vec2
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_ltol,local_vec1,local_vec2, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_ltol,local_vec1,local_vec2, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)    
  
end subroutine GeomechDiscretizationLocalToLocal

! ************************************************************************** !
!
! GeomechDiscretizationGlobalToNatural: Performs global to natural
! communication with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationGlobalToNatural(discretization,global_vec, &
                                                natural_vec,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: global_vec
  Vec                                           :: natural_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)

  call VecScatterBegin(dm_ptr%gmdm%scatter_gton,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_gton,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)       
  
end subroutine GeomechDiscretizationGlobalToNatural

! ************************************************************************** !
!
! GeomechDiscretizationNaturalToGlobal: Performs natural to global 
! communication with DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationNaturalToGlobal(discretization,natural_vec, &
                                                global_vec,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: global_vec
  Vec                                           :: natural_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_ntog,natural_vec,global_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_ntog,natural_vec,global_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  
end subroutine GeomechDiscretizationNaturalToGlobal

! ************************************************************************** !
!
! GeomechDiscretizationGlobalToLocalBegin: Begins global to local 
! communication with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationGlobalToLocalBegin(discretization,global_vec, &
                                                   local_vec,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: local_vec
  Vec                                           :: global_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_gtol,global_vec,local_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  
end subroutine GeomechDiscretizationGlobalToLocalBegin
  
! ************************************************************************** !
!
! GeomechDiscretizationGlobalToLocalEnd: Ends global to local communication
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationGlobalToLocalEnd(discretization,global_vec, &
                                                 local_vec,dm_index)

 implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: local_vec
  Vec                                           :: global_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterEnd(dm_ptr%gmdm%scatter_gtol,global_vec,local_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  
end subroutine GeomechDiscretizationGlobalToLocalEnd

! ************************************************************************** !
!
! GeomechDiscretizationLocalToLocalBegin: Begins local to local communication 
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationLocalToLocalBegin(discretization,local_vec1, &
                                                  local_vec2,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: local_vec1
  Vec                                           :: local_vec2
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_ltol,local_vec1,local_vec2, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)

end subroutine GeomechDiscretizationLocalToLocalBegin
  
! ************************************************************************** !
!
! GeomechDiscretizationLocalToLocalEnd: Ends local to local communication 
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizationLocalToLocalEnd(discretization,local_vec1, &
                                                local_vec2,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: local_vec1
  Vec                                           :: local_vec2
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterEnd(dm_ptr%gmdm%scatter_ltol,local_vec1,local_vec2, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)    

end subroutine GeomechDiscretizationLocalToLocalEnd
  
! ************************************************************************** !
!
! GeomechDiscretizGlobalToNaturalBegin: Begins global to natural communication
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizGlobalToNaturalBegin(discretization,global_vec, &
                                                natural_vec,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: global_vec
  Vec                                           :: natural_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_gton,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  
end subroutine GeomechDiscretizGlobalToNaturalBegin

! ************************************************************************** !
!
! GeomechDiscretizGlobalToNaturalEnd: Ends global to natural communication 
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizGlobalToNaturalEnd(discretization,global_vec, &
                                              natural_vec,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: global_vec
  Vec                                           :: natural_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call VecScatterEnd(dm_ptr%gmdm%scatter_gton,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)       
  
end subroutine GeomechDiscretizGlobalToNaturalEnd

! ************************************************************************** !
!
! GeomechDiscretizNaturalToGlobalBegin: Begins natural to global communication
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizNaturalToGlobalBegin(discretization,natural_vec, & 
                                                global_vec,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: global_vec
  Vec                                           :: natural_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
end subroutine GeomechDiscretizNaturalToGlobalBegin

! ************************************************************************** !
!
! GeomechDiscretizNaturalToGlobalEnd: Ends natural to global communication
! with geomech DM
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretizNaturalToGlobalEnd(discretization,natural_vec, &
                                              global_vec,dm_index)

  implicit none
  
  type(geomech_discretization_type)             :: discretization
  type(gmdm_ptr_type), pointer                  :: dm_ptr
  Vec                                           :: global_vec
  Vec                                           :: natural_vec
  PetscInt                                      :: dm_index
  PetscErrorCode                                :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
end subroutine GeomechDiscretizNaturalToGlobalEnd

! ************************************************************************** !
!
! GeomechDiscretAOApplicationToPetsc: Maps application ordering to petsc
! author: Satish Karra, LANL
! date: 06/02/13
!
! ************************************************************************** !
subroutine GeomechDiscretAOApplicationToPetsc(discretization,int_array)

  implicit none
  
#include "finclude/petscao.h"  
  
  type(geomech_discretization_type)             :: discretization
  PetscInt                                      :: int_array(:)
  PetscErrorCode                                :: ierr
  AO                                            :: ao
  
  ao = discretization%grid%ao_natural_to_petsc_nodes
  
  call AOApplicationToPetsc(ao,size(int_array),int_array,ierr)
  
end subroutine GeomechDiscretAOApplicationToPetsc

! ************************************************************************** !
!
! GeomechDiscretizationDestroy: Deallocates a geomechanics discretization
! author: Satish Karra, LANL
! date: 05/23/2013
!
! ************************************************************************** !
subroutine GeomechDiscretizationDestroy(discretization)

  implicit none
  
  type(geomech_discretization_type), pointer :: discretization
  
  PetscErrorCode :: ierr
  PetscInt :: i
    
  if (.not.associated(discretization)) return
      
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      if (discretization%dm_1dof%dm /= 0) &
        call DMDestroy(discretization%dm_1dof%dm,ierr)
      discretization%dm_1dof%dm = 0
      if (discretization%dm_ngeodof%dm /= 0) &
        call DMDestroy(discretization%dm_ngeodof%dm,ierr)
      discretization%dm_ngeodof%dm = 0
    case(UNSTRUCTURED_GRID)
      if (associated(discretization%dm_1dof%gmdm)) &
        call GMDMDestroy(discretization%dm_1dof%gmdm)
      if (associated(discretization%dm_ngeodof%gmdm)) &
        call GMDMDestroy(discretization%dm_ngeodof%gmdm)
  end select
  if (associated(discretization%dm_1dof)) &
    deallocate(discretization%dm_1dof)
  nullify(discretization%dm_1dof)
  if (associated(discretization%dm_ngeodof)) &
    deallocate(discretization%dm_ngeodof)
  nullify(discretization%dm_ngeodof)

end subroutine GeomechDiscretizationDestroy

end module Geomechanics_Discretization_module
#endif