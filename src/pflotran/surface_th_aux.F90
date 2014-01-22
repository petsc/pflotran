#ifdef SURFACE_FLOW

module Surface_TH_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public :: Surface_TH_auxvar_type
    PetscReal :: h        ! enthalpy -- not currently used
    PetscReal :: u        ! internal energy -- not currently used
    PetscReal :: pc       ! pressure change -- not currently used
    PetscReal :: Cw       ! Specific heat capacity of surface water
    PetscReal :: Ci       ! Specific heat capacity of surface ice
    PetscReal :: Cwi      ! Weighted average of Cw and Ci
      ! RTM: Note that I believe we can simple make Cw and Ci Fortran 
      ! parameters, but I'm keeping things as is for now.
    PetscReal :: k_therm  ! Thermal conductivity of surface water
    PetscReal :: unfrozen_fraction ! Proportion of unfrozen surface water.
    PetscReal :: den_water_kg  ! Density [kg/m^3] of liquid water ONLY.
      ! Note that we currently use the den_kg(1) field of the global aux var type 
      ! to store the density of the liquid/ice mixture.  We also need to track
      ! the liquid density because it is required in in the advective term of 
      ! the temperature equation.  
  end type Surface_TH_auxvar_type

  type, public :: Surface_TH_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: aux_vars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(Surface_TH_auxvar_type), pointer :: aux_vars(:)
    type(Surface_TH_auxvar_type), pointer :: aux_vars_bc(:)
    type(Surface_TH_auxvar_type), pointer :: aux_vars_ss(:)
  end type Surface_TH_type

  public :: SurfaceTHAuxCreate, &
            SurfaceTHAuxDestroy, &
            SurfaceTHAuxVarCompute, &
            SurfaceTHAuxVarComputeUnfrozen, &
            SurfaceTHAuxVarInit, &
            SurfaceTHAuxVarCopy

contains

! ************************************************************************** !
!> This routine creates an empty object
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
function SurfaceTHAuxCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(Surface_TH_type), pointer :: SurfaceTHAuxCreate
  
  type(Surface_TH_type), pointer :: aux

  allocate(aux)
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  SurfaceTHAuxCreate => aux
  
end function SurfaceTHAuxCreate

! ************************************************************************** !
!> This routine initilizes an object.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(Surface_TH_auxvar_type) :: aux_var
  type(option_type) :: option

  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%pc = 0.d0
  aux_var%Cw = 4.188d3     ! [J/kg/K]
  aux_var%Ci = 2.050d3     ! [J/kg/K]
  aux_var%Cwi = 4.188d3     ! [J/kg/K]
  aux_var%k_therm = 0.57d0 ! [J/s/m/K]
  aux_var%unfrozen_fraction = 1.d0
  aux_var%den_water_kg = 1.d3  ! [kg/m^3]

end subroutine SurfaceTHAuxVarInit

! ************************************************************************** !
!> This routine makes a copy of an object.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(Surface_TH_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%h = aux_var%h
  aux_var2%u = aux_var%u
  aux_var2%pc = aux_var%pc
  aux_var2%Cw = aux_var%Cw
  aux_var2%Ci = aux_var%Ci
  aux_var2%Cwi = aux_var%Cwi
  aux_var2%k_therm = aux_var%k_therm
  aux_var2%unfrozen_fraction = aux_var%unfrozen_fraction
  aux_var2%den_water_kg = aux_var%den_water_kg

end subroutine SurfaceTHAuxVarCopy

! ************************************************************************** !
!> This routine computes values for auxiliary variables.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHAuxVarCompute(xx,aux_var,global_aux_var, &
                                  option)

  use Option_module
  use Surface_Global_Aux_module
  use Water_EOS_module
  use Saturation_Function_module  
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: xx(option%nflowdof)
  type(Surface_TH_auxvar_type) :: aux_var
  type(surface_global_auxvar_type) :: global_aux_var
  PetscReal :: por, perm

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: di_kg, dwi_kg
    ! Densities of ice and water-ice mixture, respectively.
  PetscReal :: unfrozen_fraction
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: k_therm_w, k_therm_i
  
  global_aux_var%den_kg(1) = 0.d0

  aux_var%h = 0.d0
  aux_var%u = 0.d0
  kr = 0.d0
 
  global_aux_var%head(1) = xx(1)
  !global_aux_var%temp(1) = xx(2)
    ! RTM: Why is the above commented out?  Is one of these internal 
    ! energy instead of temperature?
 

!***************  Liquid phase properties **************************

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0

  call wateos_noderiv(global_aux_var%temp(1),pw,dw_kg,dw_mol,hw,option%scale,ierr)
  global_aux_var%den_kg(1) = dw_kg
  di_kg = 917.d0 ![kg/m^3]
    ! RTM: WARNING!  We are hard-coding the density of ice at atmospheric 
    ! pressure here.  We should actually compute this according to the 
    ! reference pressure in case someone wants to use PFLOTRAN for planetary 
    ! science.
  k_therm_w = 0.57d0 ! [J/s/m/K]
  k_therm_i = 2.18d0 ! [J/s/m/K]
    ! RTM: Same warning for thermal conductivities--these should be computed.
  
  ! RTM: These are being set but we are not using them now.  We should get rid 
  ! of them if we settle on not using an enthalpy formulation for the energy 
  ! equation.
  aux_var%h = hw
  aux_var%u = aux_var%h - pw / dw_mol * option%scale

  ! Compute unfrozen fraction, and then compute the weighted averages of 
  ! density, specific heat capacity, thermal conductivity
  unfrozen_fraction = SurfaceTHAuxVarComputeUnfrozen(global_aux_var%temp(1))
  aux_var%unfrozen_fraction = unfrozen_fraction
  global_aux_var%den_kg(1) = unfrozen_fraction * dw_kg + (1.d0 - unfrozen_fraction) * di_kg
  aux_var%Cwi = unfrozen_fraction * aux_var%Cw + (1.d0 - unfrozen_fraction) * aux_var%Ci
  aux_var%k_therm = unfrozen_fraction * k_therm_w + (1.d0 - unfrozen_fraction) * k_therm_i
  
  
end subroutine SurfaceTHAuxVarCompute

! ************************************************************************** !
!> This returns the unfrozen fraction, given a temperature.
!!
!! RTM: This is currently a public function, but is only called within 
!! Surface_TH_Aux_module, so we may want to make it private.
!! We may also want to move this into its own module if we add some 
!! more sophisticated ways of computing this quantity.
! ************************************************************************** !
function SurfaceTHAuxVarComputeUnfrozen(temp)

  implicit none

  PetscReal :: SurfaceTHAuxVarComputeUnfrozen

  PetscReal :: temp

  ! For now, we simply say that everything is unfrozen if the temperature is 
  ! above the nominal freezing point, and everything is frozen otherwise.
  ! We may want to turn this step function into a steep sigmoidal one for 
  ! numerical reasons.
  ! At some point, we may also want to consider a more detailed model of the 
  ! freezing and thawing process that might partition the frozen and unfrozen 
  ! phases according to a more physical mechanism (maybe using state transition
  ! theory?).
  if (temp > 0.d0) then
    SurfaceTHAuxVarComputeUnfrozen = 1.d0
  else
    SurfaceTHAuxVarComputeUnfrozen = 0.d0
  endif

end function SurfaceTHAuxVarComputeUnfrozen

! ************************************************************************** !
!> This routine deallocates an object.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHAuxDestroy(aux)

  implicit none

  type(Surface_TH_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%aux_vars)) deallocate(aux%aux_vars)
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) deallocate(aux%aux_vars_bc)
  nullify(aux%aux_vars_bc)
  if (associated(aux%aux_vars_ss)) deallocate(aux%aux_vars_ss)
  nullify(aux%aux_vars_ss)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)

  deallocate(aux)
  nullify(aux)  

  end subroutine SurfaceTHAuxDestroy

end module Surface_TH_Aux_module

#endif
