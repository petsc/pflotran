module Hydrostatic_module
 
  implicit none

  private

#include "include/finclude/petsc.h"

#include "definitions.h"

  real*8, parameter ::  fmwnacl = 58.44277D0,  fmwh2o = 18.0153D0

  public :: HydrostaticUpdateCoupler, HydrostaticUpdateCouplerBetter
 
contains

! ************************************************************************** !
!
! HydrostaticUpdateCoupler: Computes the hydrostatic initial/boundary condition
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine HydrostaticUpdateCoupler(coupler,option,grid)

  use water_eos_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Region_module
  use Structured_Grid_module

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  
  integer :: local_id, ghosted_id, iconn
  integer :: num_iteration, iz
  real*8 :: dist_x, dist_y, dist_z, delta_z
  real*8 :: rho, rho1, rho0, pressure0, pressure, temperature
  real*8 :: xm_nacl, dw_kg
  
  real*8 :: structured_pressure(grid%structured_grid%nz)
  
  real*8, parameter :: patm = 101325.d0
  
  type(condition_type), pointer :: condition
  
  type(connection_type), pointer :: cur_connection_set
  
  condition => coupler%condition
  
  xm_nacl = option%m_nacl * fmwnacl
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)

  if (coupler%itype == INITIAL_COUPLER_TYPE) then

    do iconn=1,coupler%connection%num_connections
      local_id = coupler%connection%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      dist_x = grid%x(ghosted_id)-condition%datum(1)
      dist_y = grid%y(ghosted_id)-condition%datum(2)
      dist_z = grid%z(ghosted_id)-condition%datum(3)
      
      pressure0 = condition%cur_value(1) + &
                  condition%gradient(1,X_DIRECTION)*dist_x + & ! gradient in Pa/m
                  condition%gradient(1,Y_DIRECTION)*dist_y + &
                  condition%gradient(1,Z_DIRECTION)*dist_z 

      temperature = condition%cur_value(2) + &
                    condition%gradient(2,X_DIRECTION)*dist_x + & ! gradient in K/m
                    condition%gradient(2,Y_DIRECTION)*dist_y + &
                    condition%gradient(2,Z_DIRECTION)*dist_z 

      pressure = pressure0
      call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
      rho = dw_kg * 1.d3
      rho1 = rho

      num_iteration = 0
      do 
        ! for the standard +x,+y,+z grid, positive dist_z results in a pressure reduction
        pressure = pressure0 + rho * option%gravity(3) * dist_z
        call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
        rho1 = dw_kg * 1.d3
        if (abs(rho-rho1) < 1.d-10) exit
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
        coupler%aux_int_var(1,iconn) = condition%iphase
      else
        coupler%aux_int_var(1,iconn) = 3
      endif

    enddo
    
  else
  
    pressure0 = condition%cur_value(1)
    temperature = condition%cur_value(2)
    pressure = pressure0
    
    call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
    rho = dw_kg * 1.d3
    
    do iz=1,grid%structured_grid%nz
      if (iz == 1) then
        delta_z = 0.5d0*grid%structured_grid%dz0(iz)
!        dist_z = delta_z
      else
        delta_z = 0.5d0*(grid%structured_grid%dz0(iz)+grid%structured_grid%dz0(iz-1))
!        dist_z = dist_z + delta_z
      endif

      call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
      rho = dw_kg * 1.d3
      rho0 = rho

      num_iteration = 0
      do 
        if (iz == 1) then
          ! for the standard +x,+y,+z grid, positive dist_z results in a pressure reduction
          pressure = pressure0 + rho * option%gravity(3) * delta_z
          call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
          rho = dw_kg * 1.d3
          exit
        else
          pressure = pressure0 + 0.5d0*(rho*grid%structured_grid%dz0(iz)+ &
                                        rho0*grid%structured_grid%dz0(iz-1))* &
                                        option%gravity(3)
        endif
        call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
        rho1 = dw_kg * 1.d3
        if (abs(rho-rho1) < 1.d-10) exit
        rho = rho1
        num_iteration = num_iteration + 1
        if (num_iteration > 100) then
          print *,'HydrostaticInitCondition failed to converge',num_iteration,rho1,rho
          stop
        endif
      enddo
      structured_pressure(iz) = pressure
      pressure0 = pressure
      rho0 = rho  
    enddo
    
    do iconn=1,coupler%connection%num_connections
      local_id = coupler%connection%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      dist_z = 0.d0
      do iz=1,grid%structured_grid%nz
        if (grid%z(ghosted_id) > dist_z .and. &
            grid%z(ghosted_id) < dist_z + grid%structured_grid%dz0(iz)) exit
        dist_z = dist_z + grid%structured_grid%dz0(iz)
      enddo
      if (iz > grid%structured_grid%nz) then
        print *, 'Error setting up hydrostatic bc.  iz > nz'
        stop
      endif
      coupler%aux_real_var(1,iconn) = structured_pressure(iz)
      coupler%aux_real_var(2,iconn) = temperature
      coupler%aux_real_var(3,iconn) = condition%cur_value(3)

      if (structured_pressure(iz) > patm) then
        coupler%aux_int_var(1,iconn) = condition%iphase
      else
        coupler%aux_int_var(1,iconn) = 3
      endif
      
    enddo
  
  endif
    
end subroutine HydrostaticUpdateCoupler

! ************************************************************************** !
!
! HydrostaticUpdateCouplerBetter: Computes the hydrostatic initial/boundary 
!                                 condition (more accurately than before0
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
subroutine HydrostaticUpdateCouplerBetter(coupler,option,grid)

  use water_eos_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Region_module
  use Structured_Grid_module

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  
  integer :: local_id, ghosted_id, iconn
  integer :: num_iteration, iz
  real*8 :: dist_x, dist_y, dist_z, delta_z
  real*8 :: rho, rho1, rho0, pressure0, pressure, temperature
  real*8 :: xm_nacl, dw_kg
  real*8 :: sign = -1.d0
  

end subroutine HydrostaticUpdateCouplerBetter

end module Hydrostatic_module

