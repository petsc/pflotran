module Richards_Aux_module

#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

  type, public :: richards_auxvar_type
  
    PetscReal :: pc
#ifdef USE_ANISOTROPIC_MOBILITY
    PetscReal :: kvr_x
    PetscReal :: kvr_y
    PetscReal :: kvr_z
    PetscReal :: dkvr_x_dp
    PetscReal :: dkvr_y_dp
    PetscReal :: dkvr_z_dp
#else
    PetscReal :: kvr
    PetscReal :: dkvr_dp
#endif
    PetscReal :: dsat_dp
    PetscReal :: dden_dp

  end type richards_auxvar_type
  
  type, public :: richards_parameter_type
    PetscReal, pointer :: sir(:,:)
  end type richards_parameter_type
  
  type, public :: richards_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    PetscBool :: aux_vars_up_to_date
    PetscBool :: aux_vars_cell_pressures_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(richards_parameter_type), pointer :: richards_parameter
    type(richards_auxvar_type), pointer :: aux_vars(:)
    type(richards_auxvar_type), pointer :: aux_vars_bc(:)
    type(richards_auxvar_type), pointer :: aux_vars_ss(:)
#ifdef BUFFER_MATRIX
    type(matrix_buffer_type), pointer :: matrix_buffer
#endif
  end type richards_type

  public :: RichardsAuxCreate, RichardsAuxDestroy, &
            RichardsAuxVarCompute, RichardsAuxVarInit, &
            RichardsAuxVarCopy

contains


! ************************************************************************** !
!
! RichardsAuxCreate: Allocate and initialize auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
function RichardsAuxCreate()

  use Option_module

  implicit none
  
  type(richards_type), pointer :: RichardsAuxCreate
  
  type(richards_type), pointer :: aux

  allocate(aux) 
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%aux_vars_cell_pressures_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)
  aux%n_zero_rows = 0
  allocate(aux%richards_parameter)
  ! don't allocate richards_parameter%sir quite yet, since we don't know the
  ! number of saturation functions
  nullify(aux%richards_parameter%sir)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
#ifdef BUFFER_MATRIX
  nullify(aux%matrix_buffer)
#endif

  RichardsAuxCreate => aux
  
end function RichardsAuxCreate

! ************************************************************************** !
!
! RichardsAuxVarInit: Initialize auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RichardsAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: aux_var
  type(option_type) :: option
  
  aux_var%pc = 0.d0

#ifdef USE_ANISOTROPIC_MOBILITY
  aux_var%kvr_x = 0.d0
  aux_var%kvr_y = 0.d0
  aux_var%kvr_z = 0.d0
  aux_var%dkvr_x_dp = 0.d0
  aux_var%dkvr_y_dp = 0.d0
  aux_var%dkvr_z_dp = 0.d0
#else
  aux_var%kvr = 0.d0
  aux_var%dkvr_dp = 0.d0
#endif

  aux_var%dsat_dp = 0.d0
  aux_var%dden_dp = 0.d0

end subroutine RichardsAuxVarInit

! ************************************************************************** !
!
! RichardsAuxVarCopy: Copies an auxiliary variable
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine RichardsAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%pc = aux_var%pc

#ifdef USE_ANISOTROPIC_MOBILITY
  aux_var2%kvr_x = aux_var%kvr_x 
  aux_var2%kvr_y = aux_var%kvr_y 
  aux_var2%kvr_z = aux_var%kvr_z 
  aux_var2%dkvr_x_dp = aux_var%dkvr_x_dp 
  aux_var2%dkvr_y_dp = aux_var%dkvr_y_dp 
  aux_var2%dkvr_z_dp = aux_var%dkvr_z_dp 
#else
  aux_var2%kvr = aux_var%kvr
  aux_var2%dkvr_dp = aux_var%dkvr_dp
#endif

  aux_var2%dsat_dp = aux_var%dsat_dp
  aux_var2%dden_dp = aux_var%dden_dp
 
end subroutine RichardsAuxVarCopy
  
! ************************************************************************** !
!
! RichardsAuxVarCompute: Computes auxiliary variables for each grid cell
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RichardsAuxVarCompute(x,aux_var,global_aux_var,&
                                 saturation_function,por,perm,option)

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Saturation_Function_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(richards_auxvar_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscReal :: por, perm

  PetscInt :: i
  PetscBool :: saturated
  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: pert, pw_pert, dw_kg_pert
  PetscReal :: fs, ani_A, ani_B, ani_C, ani_n, ani_coef
  PetscReal, parameter :: tol = 1.d-3
  
  global_aux_var%sat = 0.d0
  global_aux_var%den = 0.d0
  global_aux_var%den_kg = 0.d0
  
#ifdef USE_ANISOTROPIC_MOBILITY
  aux_var%kvr_x = 0.d0
  aux_var%kvr_y = 0.d0
  aux_var%kvr_z = 0.d0
#else
  aux_var%kvr = 0.d0
#endif

  kr = 0.d0
 
  global_aux_var%pres = x(1)
  global_aux_var%temp = option%reference_temperature
 
  aux_var%pc = option%reference_pressure - global_aux_var%pres(1)
  
!***************  Liquid phase properties **************************
  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
  if (aux_var%pc > 0.d0) then
    saturated = PETSC_FALSE
    call SaturationFunctionCompute(global_aux_var%pres(1), &
                                   global_aux_var%sat(1),kr, &
                                   ds_dp,dkr_dp, &
                                   saturation_function, &
                                   por,perm, &
                                   saturated, &
                                   option)
  else
    saturated = PETSC_TRUE
  endif  
  
  ! the purpose for splitting this condition from the 'else' statement
  ! above is due to SaturationFunctionCompute switching a cell to
  ! saturated to prevent unstable (potentially infinite) derivatives when 
  ! capillary pressure is very small
  if (saturated) then
    aux_var%pc = 0.d0
    global_aux_var%sat = 1.d0  
    kr = 1.d0    
    pw = global_aux_var%pres(1)
  endif

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
#ifndef DONT_USE_WATEOS
!geh  call EOSWaterDensityEnthalpy(global_aux_var%temp(1),pw,dw_kg,dw_mol,hw, &
!                               dw_dp,dw_dt,hw_dp,hw_dt,option%scale,ierr)
  call EOSWaterDensity(global_aux_var%temp(1),pw,dw_kg,dw_mol, &
                       dw_dp,dw_dt,option%scale,ierr)
#else
  call EOSWaterdensity(global_aux_var%temp(1),pw,dw_kg)
  pert = tol*pw
  pw_pert = pw+pert
  call EOSWaterdensity(global_aux_var%temp(1),pw_pert,dw_kg_pert)
  dw_dp = (dw_kg_pert-dw_kg)/pert
  ! dw_kg = kg/m^3
  ! dw_mol = kmol/m^3
  ! FMWH2O = kg/kmol h2o
  dw_mol = dw_kg/FMWH2O
  dw_dp = dw_dp/FMWH2O
#endif
! may need to compute dpsat_dt to pass to VISW
  call EOSWaterSaturationPressure(global_aux_var%temp(1),sat_pressure,ierr)
  !geh: 0.d0 passed in for derivative of pressure w/respect to temp
  call EOSWaterViscosity(global_aux_var%temp(1),pw,sat_pressure,0.d0, &
                         visl,dvis_dt,dvis_dp,dvis_dpsat,ierr) 
!geh  dvis_dpsat = -dvis_dp   ! already handled in EOSWaterViscosity
  if (.not.saturated) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dvis_dpsat = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif
 
  global_aux_var%den = dw_mol
  global_aux_var%den_kg = dw_kg
  aux_var%dsat_dp = ds_dp
  aux_var%dden_dp = dw_dp

#ifdef USE_ANISOTROPIC_MOBILITY  
  aux_var%kvr_x = kr/visl       ! For anisotropic relative perm
  aux_var%kvr_y = kr/visl       ! For anisotropic relative perm         
  aux_var%kvr_z = kr/visl       ! For anisotropic relative perm
  if (option%ani_relative_permeability) then  
    ani_coef = 1
!     do i=1, 100
!     global_aux_var%sat(1) = 0.01*i
!     ani_A = 3
!     ani_B = 44.8
!     ani_C = -7.26
     ani_A = saturation_function%ani_A  
     ani_B = saturation_function%ani_B
     ani_C = saturation_function%ani_C
     fs = ani_A + ani_B*exp(ani_C*global_aux_var%sat(1))
     ani_n = 25 
     ani_coef  =  fs/((global_aux_var%sat(1)**ani_n) * (fs -1) + 1)
     aux_var%kvr_z = aux_var%kvr_z * ani_coef
!    write(*,*) global_aux_var%sat(1), ani_coef
!    end do
!    stop
  end if  
  aux_var%dkvr_x_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp ! For anisotropic relative perm 
  aux_var%dkvr_y_dp = (dkr_dp/visl - kr/(visl*visl)*dvis_dp) ! For anisotropic relative perm
  aux_var%dkvr_z_dp = (dkr_dp/visl - kr/(visl*visl)*dvis_dp) ! For anisotropic relative perm
  aux_var%dkvr_z_dp = aux_var%dkvr_z_dp * ani_coef
#else
  aux_var%kvr = kr/visl
  aux_var%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
#endif

end subroutine RichardsAuxVarCompute

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a richards auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(richards_auxvar_type) :: aux_var
  
end subroutine AuxVarDestroy

! ************************************************************************** !
!
! RichardsAuxDestroy: Deallocates a richards auxiliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RichardsAuxDestroy(aux)

  implicit none

  type(richards_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call AuxVarDestroy(aux%aux_vars(iaux))
    enddo  
    deallocate(aux%aux_vars)
  endif
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call AuxVarDestroy(aux%aux_vars_bc(iaux))
    enddo  
    deallocate(aux%aux_vars_bc)
  endif
  nullify(aux%aux_vars_bc)
  if (associated(aux%aux_vars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call AuxVarDestroy(aux%aux_vars_ss(iaux))
    enddo  
    deallocate(aux%aux_vars_ss)
  endif
  nullify(aux%aux_vars_ss)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%richards_parameter)) then
    if (associated(aux%richards_parameter%sir)) deallocate(aux%richards_parameter%sir)
    nullify(aux%richards_parameter%sir)
    deallocate(aux%richards_parameter)
  endif
  nullify(aux%richards_parameter)

#ifdef BUFFER_MATRIX
  if (associated(aux%matrix_buffer)) then
    call MatrixBufferDestroy(aux%matrix_buffer)
  endif
  nullify(aux%matrix_buffer)
#endif
  
  deallocate(aux)
  nullify(aux)
    
end subroutine RichardsAuxDestroy

end module Richards_Aux_module
