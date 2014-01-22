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

    PetscBool :: auxvars_up_to_date
    PetscBool :: auxvars_cell_pressures_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(richards_parameter_type), pointer :: richards_parameter
    type(richards_auxvar_type), pointer :: auxvars(:)
    type(richards_auxvar_type), pointer :: auxvars_bc(:)
    type(richards_auxvar_type), pointer :: auxvars_ss(:)
#ifdef BUFFER_MATRIX
    type(matrix_buffer_type), pointer :: matrix_buffer
#endif
  end type richards_type

  public :: RichardsAuxCreate, RichardsAuxDestroy, &
            RichardsAuxVarCompute, RichardsAuxVarInit, &
            RichardsAuxVarCopy

contains

! ************************************************************************** !

function RichardsAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(richards_type), pointer :: RichardsAuxCreate
  
  type(richards_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%auxvars_cell_pressures_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
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

subroutine RichardsAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%pc = 0.d0

#ifdef USE_ANISOTROPIC_MOBILITY
  auxvar%kvr_x = 0.d0
  auxvar%kvr_y = 0.d0
  auxvar%kvr_z = 0.d0
  auxvar%dkvr_x_dp = 0.d0
  auxvar%dkvr_y_dp = 0.d0
  auxvar%dkvr_z_dp = 0.d0
#else
  auxvar%kvr = 0.d0
  auxvar%dkvr_dp = 0.d0
#endif

  auxvar%dsat_dp = 0.d0
  auxvar%dden_dp = 0.d0

end subroutine RichardsAuxVarInit

! ************************************************************************** !

subroutine RichardsAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pc = auxvar%pc

#ifdef USE_ANISOTROPIC_MOBILITY
  auxvar2%kvr_x = auxvar%kvr_x 
  auxvar2%kvr_y = auxvar%kvr_y 
  auxvar2%kvr_z = auxvar%kvr_z 
  auxvar2%dkvr_x_dp = auxvar%dkvr_x_dp 
  auxvar2%dkvr_y_dp = auxvar%dkvr_y_dp 
  auxvar2%dkvr_z_dp = auxvar%dkvr_z_dp 
#else
  auxvar2%kvr = auxvar%kvr
  auxvar2%dkvr_dp = auxvar%dkvr_dp
#endif

  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dden_dp = auxvar%dden_dp
 
end subroutine RichardsAuxVarCopy

! ************************************************************************** !

subroutine RichardsAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                                 saturation_function,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Saturation_Function_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(richards_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
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
  
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0
  
#ifdef USE_ANISOTROPIC_MOBILITY
  auxvar%kvr_x = 0.d0
  auxvar%kvr_y = 0.d0
  auxvar%kvr_z = 0.d0
#else
  auxvar%kvr = 0.d0
#endif

  kr = 0.d0
 
  global_auxvar%pres = x(1)
  global_auxvar%temp = option%reference_temperature
 
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)
  
!***************  Liquid phase properties **************************
  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
  if (auxvar%pc > 0.d0) then
    saturated = PETSC_FALSE
    call SaturationFunctionCompute(global_auxvar%pres(1), &
                                global_auxvar%sat(1),kr, &
                                ds_dp,dkr_dp, &
                                saturation_function, &
                                material_auxvar%porosity, &
                                material_auxvar%permeability(perm_xx_index), &
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
    auxvar%pc = 0.d0
    global_auxvar%sat = 1.d0  
    kr = 1.d0    
    pw = global_auxvar%pres(1)
  endif

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
#ifndef DONT_USE_WATEOS
!geh  call EOSWaterDensityEnthalpy(global_auxvar%temp(1),pw,dw_kg,dw_mol,hw, &
!                               dw_dp,dw_dt,hw_dp,hw_dt,option%scale,ierr)
  call EOSWaterDensity(global_auxvar%temp(1),pw,dw_kg,dw_mol, &
                       dw_dp,dw_dt,option%scale,ierr)
#else
  call EOSWaterdensity(global_auxvar%temp(1),pw,dw_kg)
  pert = tol*pw
  pw_pert = pw+pert
  call EOSWaterdensity(global_auxvar%temp(1),pw_pert,dw_kg_pert)
  dw_dp = (dw_kg_pert-dw_kg)/pert
  ! dw_kg = kg/m^3
  ! dw_mol = kmol/m^3
  ! FMWH2O = kg/kmol h2o
  dw_mol = dw_kg/FMWH2O
  dw_dp = dw_dp/FMWH2O
#endif
! may need to compute dpsat_dt to pass to VISW
  call EOSWaterSaturationPressure(global_auxvar%temp(1),sat_pressure,ierr)
  !geh: 0.d0 passed in for derivative of pressure w/respect to temp
  call EOSWaterViscosity(global_auxvar%temp(1),pw,sat_pressure,0.d0, &
                         visl,dvis_dt,dvis_dp,dvis_dpsat,ierr) 
!geh  dvis_dpsat = -dvis_dp   ! already handled in EOSWaterViscosity
  if (.not.saturated) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dvis_dpsat = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif
 
  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dp = dw_dp

#ifdef USE_ANISOTROPIC_MOBILITY  
  auxvar%kvr_x = kr/visl       ! For anisotropic relative perm
  auxvar%kvr_y = kr/visl       ! For anisotropic relative perm         
  auxvar%kvr_z = kr/visl       ! For anisotropic relative perm
  if (option%ani_relative_permeability) then  
    ani_coef = 1
!     do i=1, 100
!     global_auxvar%sat(1) = 0.01*i
!     ani_A = 3
!     ani_B = 44.8
!     ani_C = -7.26
     ani_A = saturation_function%ani_A  
     ani_B = saturation_function%ani_B
     ani_C = saturation_function%ani_C
     fs = ani_A + ani_B*exp(ani_C*global_auxvar%sat(1))
     ani_n = 25 
     ani_coef  =  fs/((global_auxvar%sat(1)**ani_n) * (fs -1) + 1)
     auxvar%kvr_z = auxvar%kvr_z * ani_coef
!    write(*,*) global_auxvar%sat(1), ani_coef
!    end do
!    stop
  end if  
  auxvar%dkvr_x_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp ! For anisotropic relative perm 
  auxvar%dkvr_y_dp = (dkr_dp/visl - kr/(visl*visl)*dvis_dp) ! For anisotropic relative perm
  auxvar%dkvr_z_dp = (dkr_dp/visl - kr/(visl*visl)*dvis_dp) ! For anisotropic relative perm
  auxvar%dkvr_z_dp = auxvar%dkvr_z_dp * ani_coef
#else
  auxvar%kvr = kr/visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
#endif

end subroutine RichardsAuxVarCompute

! ************************************************************************** !

subroutine AuxVarDestroy(auxvar)
  ! 
  ! Deallocates a richards auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(richards_auxvar_type) :: auxvar
  
end subroutine AuxVarDestroy

! ************************************************************************** !

subroutine RichardsAuxDestroy(aux)
  ! 
  ! Deallocates a richards auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(richards_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call AuxVarDestroy(aux%auxvars(iaux))
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call AuxVarDestroy(aux%auxvars_bc(iaux))
    enddo  
    deallocate(aux%auxvars_bc)
  endif
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call AuxVarDestroy(aux%auxvars_ss(iaux))
    enddo  
    deallocate(aux%auxvars_ss)
  endif
  nullify(aux%auxvars_ss)
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
