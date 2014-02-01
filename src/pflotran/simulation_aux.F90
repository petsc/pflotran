module Simulation_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type,public :: simulation_aux_type

    ! Note: These are GLOBAL vectors (i.e. they do not contain ghost control
    !       volumes)

    ! Size of entire subsurface domain
    Vec :: subsurf_pres
    Vec :: subsurf_temp
    Vec :: subsurf_sat
    Vec :: subsurf_den

    ! Size of surface cells of subsurface domain
    Vec :: subsurf_pres_top_bc
    Vec :: subsurf_temp_top_bc
    Vec :: subsurf_mflux_exchange_with_surf
    Vec :: subsurf_hflux_exchange_with_surf

    ! Size of entire surface domain
    Vec :: surf_head
    Vec :: surf_temp
    Vec :: surf_mflux_exchange_with_subsurf
    Vec :: surf_hflux_exchange_with_subsurf

    VecScatter :: surf_to_subsurf
    VecScatter :: subsurf_to_surf
    VecScatter :: subsurf_to_hydrogeophyics

  end type simulation_aux_type

  public :: SimAuxCreate, &
            SimAuxCopyVecScatter, &
            SimAuxCopySubsurfVec, &
            SimAuxCopySubsurfTopBCVec, &
            SimAuxCopySurfVec, &
            SimAuxDestroy

contains

! ************************************************************************** !

function SimAuxCreate()
  ! 
  ! This routine allocates auxillary object.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
  ! 

  use Option_module

  implicit none

  type (simulation_aux_type),pointer :: SimAuxCreate

  type (simulation_aux_type),pointer :: aux

  allocate(aux)
  aux%subsurf_pres = 0
  aux%subsurf_temp = 0
  aux%subsurf_sat = 0
  aux%subsurf_den = 0
  aux%subsurf_pres_top_bc = 0
  aux%subsurf_temp_top_bc = 0
  aux%subsurf_mflux_exchange_with_surf = 0
  aux%subsurf_hflux_exchange_with_surf = 0

  aux%surf_head = 0
  aux%surf_temp = 0
  aux%surf_mflux_exchange_with_subsurf = 0
  aux%surf_hflux_exchange_with_subsurf = 0

  aux%surf_to_subsurf = 0
  aux%subsurf_to_surf = 0
  aux%subsurf_to_hydrogeophyics = 0

  SimAuxCreate => aux

end function SimAuxCreate

! ************************************************************************** !

subroutine SimAuxCopyVecScatter(aux, vscat, vscat_index)
  ! 
  ! This routine copies VectorScatter to an appropriate context.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 

  implicit none

  type (simulation_aux_type),pointer :: aux
  VecScatter :: vscat
  PetscInt :: vscat_index

  PetscErrorCode :: ierr

  select case (vscat_index)
    case(SURF_TO_SUBSURF)
      call VecScatterCopy(vscat, aux%surf_to_subsurf, ierr)
    case(SUBSURF_TO_SURF)
      call VecScatterCopy(vscat, aux%subsurf_to_surf, ierr)
    case(SUBSURF_TO_HYDROGEOPHY)
      call VecScatterCopy(vscat, aux%subsurf_to_hydrogeophyics, ierr)
  end select  

end subroutine SimAuxCopyVecScatter

! ************************************************************************** !

subroutine SimAuxCopySubsurfVec(aux, subsurf_vec)
  ! 
  ! This routine creates 3D vectors related with subsurface-flow.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 

  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: subsurf_vec

  PetscErrorCode :: ierr

  call VecDuplicate(subsurf_vec,aux%subsurf_pres,ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_temp,ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_sat,ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_den,ierr)

end subroutine SimAuxCopySubsurfVec

! ************************************************************************** !

subroutine SimAuxCopySubsurfTopBCVec(aux, subsurf_top_bc_vec)
  ! 
  ! This routine creates vectors associated with surface of subsurface domain
  ! related with subsurface-flow.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 

  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: subsurf_top_bc_vec

  PetscErrorCode :: ierr

  call VecDuplicate(subsurf_top_bc_vec,aux%subsurf_pres_top_bc,ierr)
  call VecDuplicate(subsurf_top_bc_vec,aux%subsurf_temp_top_bc,ierr)
  call VecDuplicate(subsurf_top_bc_vec,aux%subsurf_mflux_exchange_with_surf,ierr)
  call VecDuplicate(subsurf_top_bc_vec,aux%subsurf_hflux_exchange_with_surf,ierr)

end subroutine SimAuxCopySubsurfTopBCVec

! ************************************************************************** !

subroutine SimAuxCopySurfVec(aux, surf_head_vec)
  ! 
  ! This routine creates vectors associated with surface-flow.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 

  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: surf_head_vec

  PetscErrorCode :: ierr

  call VecDuplicate(surf_head_vec,aux%surf_head,ierr)
  call VecDuplicate(surf_head_vec,aux%surf_temp,ierr)
  call VecDuplicate(surf_head_vec,aux%surf_mflux_exchange_with_subsurf,ierr)
  call VecDuplicate(surf_head_vec,aux%surf_hflux_exchange_with_subsurf,ierr)

end subroutine SimAuxCopySurfVec

! ************************************************************************** !

subroutine SimAuxDestroy(aux)
  ! 
  ! This routine deallocates auxillary object.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
  ! 

  implicit none

  type(simulation_aux_type), pointer :: aux

  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  if (aux%subsurf_pres /= 0) call VecDestroy(aux%subsurf_pres,ierr)
  if (aux%subsurf_temp /= 0) call VecDestroy(aux%subsurf_temp,ierr)
  if (aux%subsurf_sat /= 0) call VecDestroy(aux%subsurf_sat,ierr)
  if (aux%subsurf_den /= 0) call VecDestroy(aux%subsurf_den,ierr)
  if (aux%subsurf_pres_top_bc /= 0) call VecDestroy(aux%subsurf_pres_top_bc,ierr)
  if (aux%subsurf_temp_top_bc /= 0) call VecDestroy(aux%subsurf_temp_top_bc,ierr)
  if (aux%subsurf_mflux_exchange_with_surf /= 0) &
    call VecDestroy(aux%subsurf_mflux_exchange_with_surf,ierr)
  if (aux%subsurf_hflux_exchange_with_surf /= 0) &
    call VecDestroy(aux%subsurf_hflux_exchange_with_surf,ierr)

  if (aux%surf_head /= 0) call VecDestroy(aux%surf_head,ierr)
  if (aux%surf_temp /= 0) call VecDestroy(aux%surf_temp,ierr)
  if (aux%surf_mflux_exchange_with_subsurf /= 0) &
    call VecDestroy(aux%surf_mflux_exchange_with_subsurf,ierr)
  if (aux%surf_hflux_exchange_with_subsurf /= 0) &
    call VecDestroy(aux%surf_hflux_exchange_with_subsurf,ierr)

  if (aux%surf_to_subsurf /= 0) call VecScatterDestroy(aux%surf_to_subsurf,ierr)
  if (aux%subsurf_to_surf /= 0) call VecScatterDestroy(aux%subsurf_to_surf,ierr)
  if (aux%subsurf_to_hydrogeophyics /= 0) &
    call VecScatterDestroy(aux%subsurf_to_hydrogeophyics,ierr)

  deallocate(aux)
  nullify(aux)

end subroutine SimAuxDestroy

end module Simulation_Aux_module
