module Hydrostatic_module
 
  implicit none

  private

#include "include/finclude/petsc.h"

#include "definitions.h"

  real*8, parameter ::  fmwnacl = 58.44277D0,  fmwh2o = 18.0153D0

  public :: HydrostaticUpdateCoupler, HydrostaticUpdateCouplerBetter, &
            HydrostaticTest, HydrostaticUpdateCouplerBest
 
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
  integer :: num_iteration, ipressure, idatum, num_pressures
  real*8 :: dist_x, dist_y, dist_z, delta_z
  real*8 :: rho, rho1, rho0, pressure0, pressure, temperature
  real*8 :: xm_nacl, dw_kg
  real*8, pointer :: pressure_array(:), density_array(:), z(:)
  
  type(condition_type), pointer :: condition
  
  type(connection_type), pointer :: cur_connection_set
  
  condition => coupler%condition
    
  xm_nacl = option%m_nacl * fmwnacl
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)
  
  nullify(pressure_array)
  
  delta_z = 1.d0

  if (dabs(condition%gradient(1,Z_DIRECTION)) < 1.d-40) then
    ! compute the vertical gradient based on a 1 meter vertical spacing and
    ! interpolate the values from that array
    num_pressures = int((max(grid%z_max,condition%datum(Z_DIRECTION)) - &
                        min(grid%z_min,condition%datum(Z_DIRECTION))) / &
                        delta_z) + 1
    allocate(pressure_array(num_pressures))
    allocate(density_array(num_pressures))
    allocate(z(num_pressures))
    pressure_array = 0.d0
    density_array = 0.d0
    z = 0.d0
    ! place this pressure in the array
    idatum = int(condition%datum(Z_DIRECTION)/delta_z - &
                 min(grid%z_min,condition%datum(Z_DIRECTION)) / &
                 (max(grid%z_max,condition%datum(Z_DIRECTION)) - &
                  min(grid%z_min,condition%datum(Z_DIRECTION))) * &
                  dble(num_pressures))+1
    pressure0 = condition%cur_value(1)
    temperature = condition%cur_value(2)
    pressure_array(idatum) = pressure0
    call nacl_den(temperature,pressure0*1.d-6,xm_nacl,dw_kg) 
    rho = dw_kg * 1.d3
    density_array(idatum) = rho
    z(idatum) = condition%datum(Z_DIRECTION)
    ! compute pressures above datum, if any
    dist_z = 0.d0
    rho0 = rho
    do ipressure=idatum+1,num_pressures
      dist_z = dist_z + delta_z
      temperature = temperature + condition%gradient(2,Z_DIRECTION)*delta_z
      call nacl_den(temperature,pressure0*1.d-6,xm_nacl,dw_kg) 
      rho = dw_kg * 1.d3

      num_iteration = 0
      do 
        pressure = pressure0 + 0.5d0*(rho+rho0) * &
                   option%gravity(Z_DIRECTION) * delta_z
        call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
        rho1 = dw_kg * 1.d3
        if (abs(rho-rho1) < 1.d-10) exit
        rho = rho1
        num_iteration = num_iteration + 1
        if (num_iteration > 100) then
          print *,'Hydrostatic iteration failed to converge',num_iteration,rho1,rho
          print *, condition%name, idatum
          print *, pressure_array
          stop
        endif
      enddo
      rho0 = rho
      pressure_array(ipressure) = pressure
      density_array(ipressure) = rho
      z(ipressure) = z(idatum)+dist_z
      pressure0 = pressure
    enddo

    ! compute pressures above datum, if any
    pressure0 = pressure_array(idatum)
    temperature = condition%cur_value(2)
    dist_z = 0.d0
    rho0 = density_array(idatum)
    do ipressure=idatum-1,1,-1
      dist_z = dist_z + delta_z
      temperature = temperature - condition%gradient(2,Z_DIRECTION)*delta_z
      call nacl_den(temperature,pressure0*1.d-6,xm_nacl,dw_kg) 
      rho = dw_kg * 1.d3

      num_iteration = 0
      do                   ! notice the negative sign (-) here
        pressure = pressure0 - 0.5d0*(rho+rho0) * &
                   option%gravity(Z_DIRECTION) * delta_z
        call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
        rho1 = dw_kg * 1.d3
        if (abs(rho-rho1) < 1.d-10) exit
        rho = rho1
        num_iteration = num_iteration + 1
        if (num_iteration > 100) then
          print *,'Hydrostatic iteration failed to converge',num_iteration,rho1,rho
          print *, condition%name, idatum
          print *, pressure_array
          stop
        endif
      enddo
      rho0 = rho
      pressure_array(ipressure) = pressure
      density_array(ipressure) = rho
      z(ipressure) = z(idatum)-dist_z
      pressure0 = pressure
    enddo
  endif

  do iconn=1,coupler%connection%num_connections
    local_id = coupler%connection%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)

    dist_x = grid%x(ghosted_id)-condition%datum(X_DIRECTION)
    dist_y = grid%y(ghosted_id)-condition%datum(Y_DIRECTION)
    dist_z = grid%z(ghosted_id)-condition%datum(Z_DIRECTION)
    
    if (associated(pressure_array)) then
      ipressure = idatum+int(dist_z/delta_z)
      pressure = pressure_array(ipressure) + &
                 density_array(ipressure)*option%gravity(Z_DIRECTION) * &
                 (grid%z(ghosted_id)-z(ipressure)) + &
                 condition%gradient(1,X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 condition%gradient(1,Y_DIRECTION)*dist_y
    else
      pressure = condition%cur_value(1) + &
                 condition%gradient(1,X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 condition%gradient(1,Y_DIRECTION)*dist_y + &
                 condition%gradient(1,Z_DIRECTION)*dist_z 
    endif
    temperature = condition%cur_value(2) + &
                  condition%gradient(2,X_DIRECTION)*dist_x + & ! gradient in K/m
                  condition%gradient(2,Y_DIRECTION)*dist_y + &
                  condition%gradient(2,Z_DIRECTION)*dist_z 

    coupler%aux_real_var(1,iconn) = pressure
    coupler%aux_real_var(2,iconn) = temperature
    coupler%aux_real_var(3,iconn) = condition%cur_value(3)

    coupler%aux_int_var(1,iconn) = 1
!    if (structured_pressure(iz) > patm) then
!      coupler%aux_int_var(1,iconn) = condition%iphase
!    else
!      coupler%aux_int_var(1,iconn) = 3
!    endif
      
  enddo
  
  if (associated(pressure_array)) deallocate(pressure_array)
  nullify(pressure_array)
  
end subroutine HydrostaticUpdateCouplerBetter

! ************************************************************************** !
!
! HydrostaticUpdateCouplerBest: Computes the hydrostatic initial/boundary condition
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine HydrostaticUpdateCouplerBest(coupler,option,grid)

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
    
end subroutine HydrostaticUpdateCouplerBest

! ************************************************************************** !
!
! HydrostaticUpdateCouplerBetter: Computes the hydrostatic initial/boundary 
!                                 condition (more accurately than before0
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
subroutine HydrostaticTest()

  use water_eos_module
  
  implicit none
  
  integer :: iz, i, i_increment, num_increment
  integer :: max_num_pressures, i_up, i_dn, num_iteration
  real*8 :: rho, rho1, rho0, pressure0, pressure, temperature
  real*8 :: increment(4)
  real*8 :: xm_nacl, dw_kg, dist_z, dist

  real*8, pointer :: density_array(:,:), pressure_array(:,:)
  
  increment(1) = 1.d-1
  increment(2) = 1.d-0
  increment(3) = 1.d+1
  increment(4) = 1.d+2
  num_increment = size(increment)
  
  temperature = 25.d0

  xm_nacl = 0.d0
  
  max_num_pressures = int(1000.d0/increment(1)+0.5d0)+1
  
  allocate(density_array(max_num_pressures,num_increment))
  allocate(pressure_array(max_num_pressures,num_increment)) 
  density_array = 0.d0
  pressure_array = 0.d0
  
  do i_increment = 1, num_increment
    pressure = 101325.d0
    call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
    rho = dw_kg * 1.d3
    dist_z = 0.d0
    pressure_array(1,i_increment) = pressure
    density_array(1,i_increment) = rho
    do iz=1,int(1000.d0/increment(i_increment)+0.5d0)
      dist_z = dist_z + increment(i_increment)
      pressure0 = pressure
      num_iteration = 0
      do
        pressure = pressure0 + rho * 9.8068d0 * increment(i_increment)
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
      i = int(dist_z/increment(1)+0.5d0)+1
      pressure_array(i,i_increment) = pressure
      density_array(i,i_increment) = rho
      pressure0 = pressure
      rho0 = rho  
    enddo
  enddo

  do i_increment=2,num_increment
    dist_z = 0.d0
    i = 1
    i_up = 1
    do iz = 1,int(1000.d0/increment(i_increment)+0.5d0)
      i_dn = i_up + int(increment(i_increment)/increment(1)+1.d-6)
      dist = increment(1)
      do 
        i = i + 1
        if (dist >= 0.9999d0*increment(i_increment)) exit
        pressure_array(i,i_increment) = pressure_array(i_up,i_increment)*(1.d0-dist/increment(i_increment))+ &
                                        pressure_array(i_dn,i_increment)*dist/increment(i_increment)
        density_array(i,i_increment) = density_array(i_up,i_increment)*(1.d0-dist/increment(i_increment))+ &
                                        density_array(i_dn,i_increment)*dist/increment(i_increment)
        dist = dist + increment(1)
      enddo
      dist_z = dist_z + increment(i_increment)
      i_up = i_dn
    enddo
  enddo

  
  open(unit=86,file='pressures.dat')
  dist_z = 0.d0
  do iz = 1,max_num_pressures
    write(86,'(100(es16.10,x))') dist_z,(density_array(iz,i_increment),i_increment=1,num_increment), &
                       (pressure_array(iz,i_increment),i_increment=1,num_increment)
    dist_z = dist_z + increment(1)
  enddo
  close(86)
  
  deallocate(pressure_array)
  deallocate(density_array)
    
end subroutine HydrostaticTest

end module Hydrostatic_module

