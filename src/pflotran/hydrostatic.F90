module Hydrostatic_module
 
  implicit none

  private

#include "include/finclude/petsc.h"

#include "definitions.h"

  real*8, parameter ::  fmwnacl = 58.44277D0,  fmwh2o = 18.0153D0

  public :: HydrostaticUpdateCoupler
 
contains

subroutine HydrostaticUpdateCoupler(coupler,option,grid)

  use water_eos_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Region_module

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  
  integer :: local_id, ghosted_id, iconn
  integer :: num_iteration
  real*8 :: dist_x, dist_y, dist_z
  real*8 :: rho, rho1, pressure0, pressure, temperature
  real*8 :: xm_nacl, dw_kg
  
  real*8, parameter :: patm = 101325.d0
  
  type(condition_type), pointer :: condition
  
  type(connection_type), pointer :: cur_connection_set
  
  condition => coupler%condition
  
  xm_nacl = option%m_nacl * fmwnacl
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)

  do iconn=1,coupler%connection%num_connections
    local_id = coupler%connection%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)

    dist_x = grid%x(ghosted_id)-condition%datum(1)
    dist_y = grid%y(ghosted_id)-condition%datum(2)
    dist_z = dabs(grid%z(ghosted_id)-condition%datum(3))
    
    pressure0 = condition%cur_value(1) + &
                condition%gradient(1,X_DIRECTION)*dist_x + & ! gradient in Pa/m
                condition%gradient(1,Y_DIRECTION)*dist_y + &
                condition%gradient(1,Z_DIRECTION)*dist_z 

    temperature = condition%cur_value(2) + &
                  condition%gradient(2,X_DIRECTION)*dist_x + & ! gradient in K/m
                  condition%gradient(2,Y_DIRECTION)*dist_y + &
                  condition%gradient(2,Z_DIRECTION)*dist_z 

    call nacl_den(temperature,pressure*1D-6,xm_nacl,dw_kg) 
    rho = dw_kg * 1.d3
    pressure = rho * option%gravity * dist_z + pressure0

    num_iteration = 0
    do 
      call nacl_den(temperature,pressure*1D-6,xm_nacl,dw_kg) 
      rho = dw_kg * 1.d3
      if (abs(rho-rho1) < 1.d-10) exit
      pressure = rho * option%gravity * dist_z + pressure0
      rho = rho1
      num_iteration = num_iteration + 1
      if (num_iteration > 100) then
        print *,'HydrostaticInitCondition failed to converge',num_iteration,rho1,rho
        stop
      endif
    enddo

    coupler%aux_real_var(1,iconn) = pressure
    coupler%aux_real_var(2,iconn) = temperature
    coupler%aux_real_var(3,iconn) = condition%cur_value(3)

    if (pressure > patm) then
      coupler%aux_real_var(1,iconn) = condition%iphase
    else
      coupler%aux_real_var(1,iconn) = 3
    endif

  enddo
    
end subroutine HydrostaticUpdateCoupler

end module Hydrostatic_module

