module Hydrostatic_module
 
  implicit none

  private

#include "definitions.h"

  PetscReal, parameter ::  fmwnacl = 58.44277D0,  fmwh2o = 18.0153D0

  public :: HydrostaticUpdateCoupler, &
            HydrostaticTest
 
contains

! ************************************************************************** !
!
! HydrostaticUpdateCoupler: Computes the hydrostatic initial/boundary 
!                           condition 
! author: Glenn Hammond
! date: 11/28/07
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
  
  PetscInt :: local_id, ghosted_id, iconn
  PetscInt :: num_iteration, ipressure, idatum, num_pressures
  PetscReal :: dist_x, dist_y, dist_z, delta_z
  PetscReal :: dx_conn, dy_conn, dz_conn
  PetscReal :: rho, rho1, rho0, pressure, pressure0, pressure_at_datum 
  PetscReal :: temperature_at_datum, temperature
  PetscReal :: concentration_at_datum
  PetscReal :: xm_nacl, dw_kg
  PetscReal :: max_z, min_z
  PetscReal, pointer :: pressure_array(:), density_array(:), z(:)
  PetscReal :: pressure_gradient(3), piezometric_head_gradient(3), datum(3)
  PetscReal :: temperature_gradient(3), concentration_gradient(3)
  
  type(flow_condition_type), pointer :: condition
  
  type(connection_set_type), pointer :: cur_connection_set
  
  condition => coupler%flow_condition
    
  xm_nacl = option%m_nacl * fmwnacl
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)
  
  nullify(pressure_array)
  
  delta_z = min((grid%z_max_global-grid%z_min_global/500),1.d0)
  temperature_at_datum = 25.d0
  concentration_at_datum = 0.d0
  
  
  ! for now, just set it; in future need to account for a different temperature datum
  if (associated(condition%temperature)) then
    temperature_at_datum = condition%temperature%dataset%cur_value(1)
    temperature_gradient(1:3) = condition%temperature%gradient%cur_value(1:3)
  endif
  if (associated(condition%temperature)) then
    concentration_at_datum = condition%concentration%dataset%cur_value(1)
    concentration_gradient(1:3) = condition%concentration%gradient%cur_value(1:3)
  endif

  datum(1:3) = condition%pressure%datum%cur_value(1:3)
  pressure_at_datum = condition%pressure%dataset%cur_value(1)
  ! gradient is in m/m; needs conversion to Pa/m
  piezometric_head_gradient(1:3) = condition%pressure%gradient%cur_value(1:3)
  call nacl_den(temperature_at_datum,pressure_at_datum*1.d-6,xm_nacl,dw_kg) 
  rho = dw_kg * 1.d3
  pressure_gradient(1:3) = piezometric_head_gradient(1:3)* &
                           rho*dabs(option%gravity(Z_DIRECTION)) ! gravity is negative, but we just want magnitude

  if (dabs(pressure_gradient(Z_DIRECTION)) < 1.d-40) then
    ! compute the vertical gradient based on a 1 meter vertical spacing and
    ! interpolate the values from that array
    max_z = max(grid%z_max_global,datum(Z_DIRECTION))+1.d0 ! add 1m buffer
    min_z = min(grid%z_min_global,datum(Z_DIRECTION))-1.d0
    
    num_pressures = int((max_z-min_z)/delta_z) + 1
    allocate(pressure_array(num_pressures))
    allocate(density_array(num_pressures))
    allocate(z(num_pressures))
    pressure_array = 0.d0
    density_array = 0.d0
    z = 0.d0
    ! place this pressure in the array
    idatum = int((datum(Z_DIRECTION)-min_z)/(max_z-min_z) * &
                 dble(num_pressures))+1
    pressure_array(idatum) = pressure_at_datum
    call nacl_den(temperature_at_datum,pressure_at_datum*1.d-6,xm_nacl,dw_kg) 
    temperature = temperature_at_datum
    pressure0 = pressure_at_datum
    rho = dw_kg * 1.d3
    density_array(idatum) = rho
    z(idatum) = datum(Z_DIRECTION)
    ! compute pressures above datum, if any
    dist_z = 0.d0
    rho0 = rho
    do ipressure=idatum+1,num_pressures
      dist_z = dist_z + delta_z
      if (option%iflowmode /= RICHARDS_MODE) &
        temperature = temperature + temperature_gradient(Z_DIRECTION)*delta_z
      call nacl_den(temperature,pressure0*1.d-6,xm_nacl,dw_kg) 
      rho = dw_kg * 1.d3

      num_iteration = 0
      do 
        pressure = pressure0 + 0.5d0*(rho+rho0) * &
                   option%gravity(Z_DIRECTION) * delta_z
        call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
        rho1 = dw_kg * 1.d3
        if (dabs(rho-rho1) < 1.d-10) exit
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
    if (option%iflowmode /= RICHARDS_MODE) temperature = temperature_at_datum
    dist_z = 0.d0
    rho0 = density_array(idatum)
    do ipressure=idatum-1,1,-1
      dist_z = dist_z + delta_z
      if (option%iflowmode /= RICHARDS_MODE) &
        temperature = temperature - temperature_gradient(Z_DIRECTION)*delta_z
      call nacl_den(temperature,pressure0*1.d-6,xm_nacl,dw_kg) 
      rho = dw_kg * 1.d3

      num_iteration = 0
      do                   ! notice the negative sign (-) here
        pressure = pressure0 - 0.5d0*(rho+rho0) * &
                   option%gravity(Z_DIRECTION) * delta_z
        call nacl_den(temperature,pressure*1.d-6,xm_nacl,dw_kg) 
        rho1 = dw_kg * 1.d3
        if (dabs(rho-rho1) < 1.d-10) exit
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

  dx_conn = 0.d0
  dy_conn = 0.d0
  dz_conn = 0.d0

  do iconn=1,coupler%connection_set%num_connections
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)

    if (associated(coupler%connection_set%dist)) then
      dx_conn = coupler%connection_set%dist(0,iconn)*coupler%connection_set%dist(1,iconn)
      dy_conn = coupler%connection_set%dist(0,iconn)*coupler%connection_set%dist(2,iconn)
      dz_conn = coupler%connection_set%dist(0,iconn)*coupler%connection_set%dist(3,iconn)
    endif
    ! note the negative (-) d?_conn is required due to the offset of the boundary face
    dist_x = grid%x(ghosted_id)-dx_conn-datum(X_DIRECTION)
    dist_y = grid%y(ghosted_id)-dy_conn-datum(Y_DIRECTION)
    dist_z = grid%z(ghosted_id)-dz_conn-datum(Z_DIRECTION)

    if (associated(pressure_array)) then
      ipressure = idatum+int(dist_z/delta_z)
      dist_z = grid%z(ghosted_id)-dz_conn-z(ipressure)
      pressure = pressure_array(ipressure) + &
                 density_array(ipressure)*option%gravity(Z_DIRECTION) * &
                 dist_z + &
!                 (grid%z(ghosted_id)-z(ipressure)) + &
                 pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 pressure_gradient(Y_DIRECTION)*dist_y
    else
      pressure = pressure_at_datum + &
                 pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 pressure_gradient(Y_DIRECTION)*dist_y + &
                 pressure_gradient(Z_DIRECTION)*dist_z 
    endif

    if (condition%pressure%itype == SEEPAGE_BC) then
      coupler%flow_aux_real_var(1,iconn) = max(pressure,option%reference_pressure)
    else
      coupler%flow_aux_real_var(1,iconn) = pressure
    endif

    if (option%iflowmode /= RICHARDS_MODE) then
      temperature = temperature_at_datum + &
                    temperature_gradient(X_DIRECTION)*dist_x + & ! gradient in K/m
                    temperature_gradient(Y_DIRECTION)*dist_y + &
                    temperature_gradient(Z_DIRECTION)*dist_z 
      coupler%flow_aux_real_var(2,iconn) = temperature
      coupler%flow_aux_real_var(3,iconn) = concentration_at_datum
    endif

    coupler%flow_aux_int_var(1,iconn) = 1
!    if (structured_pressure(iz) > patm) then
!      coupler%flow_aux_int_var(1,iconn) = condition%iphase
!    else
!      coupler%flow_aux_int_var(1,iconn) = 3
!    endif
      
  enddo
  
  if (associated(pressure_array)) deallocate(pressure_array)
  nullify(pressure_array)
  
end subroutine HydrostaticUpdateCoupler

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
  
  PetscInt :: iz, i, i_increment, num_increment
  PetscInt :: max_num_pressures, i_up, i_dn, num_iteration
  PetscReal :: rho, rho1, rho0, pressure0, pressure, temperature
  PetscReal :: increment(4)
  PetscReal :: xm_nacl, dw_kg, dist_z, dist

  PetscReal, pointer :: density_array(:,:), pressure_array(:,:)
  
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
        if (dabs(rho-rho1) < 1.d-10) exit
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

