! added by S. Karra 07/11/12

module Secondary_Continuum_Aux_module

  use Reactive_Transport_Aux_module

  implicit none

  private

#include "definitions.h"

  type, public :: slab_type
    PetscReal :: length                       ! input - length of slab
    PetscReal :: area                         ! input - surface area
  end type slab_type
  
  type, public :: nested_cube_type
    PetscReal :: matrix_block_size            ! input - side of cube
    PetscReal :: fracture_spacing             ! input - fracture spacing
  end type nested_cube_type
  
  type, public :: nested_sphere_type
    PetscReal :: radius                       ! input - radius of sphere
  end type nested_sphere_type
  
  type, public :: sec_continuum_type
    PetscInt :: itype                         ! input - type of sec. continuum (slab, nested_cube, nested_sphere,....) 
    type(slab_type) :: slab
    type(nested_cube_type) :: nested_cube
    type(nested_sphere_type) :: nested_sphere 
  end type sec_continuum_type

  type, public :: sec_heat_type  
    PetscBool :: sec_temp_update               ! flag to check if the temp is updated
    PetscInt :: ncells                         ! number of secondary grid cells
    PetscReal :: aperture                      ! fracture aperture
    PetscReal :: epsilon                       ! vol. frac. of primary continuum
    type(sec_continuum_type) :: sec_continuum
    PetscReal, pointer :: sec_temp(:)          ! array of temp. at secondary grid cells
    PetscReal, pointer :: area(:)              ! surface area
    PetscReal, pointer :: vol(:)               ! volume     face      node       face
    PetscReal, pointer :: dm_plus(:)           ! see fig.    |----------o----------|
    PetscReal, pointer :: dm_minus(:)          ! see fig.      <dm_minus> <dm_plus>
    PetscReal :: interfacial_area              ! interfacial area between prim. and sec. per unit volume of prim.+sec.
    PetscBool :: log_spacing                   ! flag to check if log spacing is set
    PetscReal :: outer_spacing                 ! value of the outer most grid cell spacing
  end type sec_heat_type  
 
  type, public :: sec_transport_type  
    PetscInt :: ncells                         ! number of secondary grid cells
    PetscReal :: aperture                      ! fracture aperture
    PetscReal :: epsilon                       ! vol. frac. of primary continuum
    type(sec_continuum_type) :: sec_continuum
    type(reactive_transport_auxvar_type), pointer :: sec_rt_auxvar(:)  ! for each secondary grid cell
    PetscReal, pointer :: area(:)              ! surface area
    PetscReal, pointer :: vol(:)               ! volume     face      node       face
    PetscReal, pointer :: dm_plus(:)           ! see fig.    |----------o----------|
    PetscReal, pointer :: dm_minus(:)          ! see fig.      <dm_minus> <dm_plus>
    PetscReal :: interfacial_area              ! interfacial area between prim. and sec. per unit volume of prim.+sec.
    PetscBool :: log_spacing                   ! flag to check if log spacing is set
    PetscReal :: outer_spacing                 ! value of the outer most grid cell spacing
    PetscReal, pointer :: updated_conc(:,:)    ! This stores the secondary concentration update values from secondary NR iteration  (naqcomp x ncells)
    PetscReal, pointer :: sec_jac(:,:)         ! stores the secondary continuum jacobian value (naqcomp x naqcomp)
    PetscBool :: sec_jac_update                ! flag to check if secondary jacobian is updated
    PetscReal, pointer :: cxm(:,:,:)           ! stores the coeff of left diag in block triag system (ncomp x ncomp x ncells-1)
    PetscReal, pointer :: cxp(:,:,:)           ! stores the coeff of right diag in block triag system (ncomp x ncomp x ncells-1)
    PetscReal, pointer :: cdl(:,:,:)           ! stores the coeff of central diag in block triag system (ncomp x ncomp x ncells)
    PetscReal, pointer :: rhs(:)                 ! stores the solution of the forward solve
  end type sec_transport_type  
        

  type, public :: sc_heat_type
    type(sec_heat_type), pointer :: sec_heat_vars(:)
  end type sc_heat_type

  type, public :: sc_rt_type
    type(sec_transport_type), pointer :: sec_transport_vars(:)
  end type sc_rt_type

  public :: SecondaryAuxHeatCreate, SecondaryAuxHeatDestroy, &
            SecondaryAuxRTCreate, SecondaryAuxRTDestroy
            
contains
  
  
! ************************************************************************** !
!
! SecondaryAuxHeatCreate: Allocate and initialize secondary continuum heat
! auxiliary object
! author: Satish Karra
! date: 01/10/13
!
! ************************************************************************** !
function SecondaryAuxHeatCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(sc_heat_type), pointer :: SecondaryAuxHeatCreate
  
  type(sc_heat_type), pointer :: aux

  allocate(aux) 
  nullify(aux%sec_heat_vars)
  
  SecondaryAuxHeatCreate => aux
  
end function SecondaryAuxHeatCreate  
  
! ************************************************************************** !
!
! SecondaryAuxHeatDestroy: Deallocates a secondary continuum heat
! auxiliary object
! author: Satish Karra
! date: 01/10/13
!
! ************************************************************************** !
subroutine SecondaryAuxHeatDestroy(aux)

  implicit none

  type(sc_heat_type), pointer :: aux
   
  if (.not.associated(aux)) return
  
  deallocate(aux)
  nullify(aux)  

end subroutine SecondaryAuxHeatDestroy


! ************************************************************************** !
!
! SecondaryAuxRTCreate: Allocate and initialize secondary continuum
! reactive transport auxiliary object
! author: Satish Karra
! date: 01/10/13
!
! ************************************************************************** !
function SecondaryAuxRTCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(sc_rt_type), pointer :: SecondaryAuxRTCreate
  
  type(sc_rt_type), pointer :: aux

  allocate(aux) 
  nullify(aux%sec_transport_vars)
  
  SecondaryAuxRTCreate => aux
  
end function SecondaryAuxRTCreate  
  
! ************************************************************************** !
!
! SecondaryAuxRTDestroy: Deallocates a secondary continuum reactive 
! transport auxiliary object
! author: Satish Karra
! date: 01/10/13
!
! ************************************************************************** !
subroutine SecondaryAuxRTDestroy(aux)

  implicit none

  type(sc_rt_type), pointer :: aux
   
  if (.not.associated(aux)) return
  
  deallocate(aux)
  nullify(aux)  

end subroutine SecondaryAuxRTDestroy


end module Secondary_Continuum_Aux_module
            