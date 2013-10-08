module Hydrostatic_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

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

  use Water_EOS_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Region_module
  use Structured_Grid_module
  use Utility_module, only : DotProduct
  use Dataset_Gridded_class
  use Dataset_Ascii_class
  
  use General_Aux_module
  
  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  
  PetscInt :: local_id, ghosted_id, iconn
  PetscInt :: num_iteration, ipressure, idatum, num_pressures
  PetscReal :: dist_x, dist_y, dist_z, delta_z, dist_z_for_pressure
  PetscReal :: dx_conn, dy_conn, dz_conn
  PetscReal :: rho, rho1, rho0, pressure, pressure0, pressure_at_datum 
  PetscReal :: temperature_at_datum, temperature
  PetscReal :: concentration_at_datum
  PetscReal :: xm_nacl, dw_kg
  PetscReal :: max_z, min_z, temp_real
  PetscInt  :: num_faces, face_id_ghosted, conn_id, num_regions
  type(connection_set_type), pointer :: conn_set_ptr
  PetscReal, pointer :: pressure_array(:)
  PetscReal, allocatable :: density_array(:), z(:)
  PetscReal :: pressure_gradient(3), piezometric_head_gradient(3), datum(3)
  PetscReal :: temperature_gradient(3), concentration_gradient(3)
  PetscReal :: gravity_magnitude
  PetscReal :: z_offset
  
  class(dataset_gridded_type), pointer :: datum_dataset
  PetscReal :: datum_dataset_rmax
  PetscReal :: datum_dataset_rmin
  
  type(flow_condition_type), pointer :: condition
  
  type(connection_set_type), pointer :: cur_connection_set
  
  condition => coupler%flow_condition
  
  xm_nacl = option%m_nacl * FMWNACL
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)
  
  nullify(pressure_array)
  nullify(datum_dataset)
  
  delta_z = min((grid%z_max_global-grid%z_min_global)/500,1.d0)
  temperature_at_datum = option%reference_temperature
  concentration_at_datum = 0.d0
  datum = 0.d0
  
  select case(option%iflowmode)
    case(G_MODE)
      temperature_at_datum = &
        condition%general%temperature%dataset%rarray(1)
      temperature_gradient(1:3) = &
        condition%general%temperature%gradient%rarray(1:3)
      concentration_at_datum = &
        condition%general%mole_fraction%dataset%rarray(1)
      concentration_gradient(1:3) = &
        condition%general%mole_fraction%gradient%rarray(1:3)
      datum(1:3) = condition%datum%rarray(1:3)
      pressure_at_datum = &
        condition%general%liquid_pressure%dataset%rarray(1)    
      ! gradient is in m/m; needs conversion to Pa/m
      piezometric_head_gradient(1:3) = &
        condition%general%liquid_pressure%gradient%rarray(1:3)
    case default
      ! for now, just set it; in future need to account for a different temperature datum
      if (associated(condition%temperature)) then
        if (condition%temperature%itype == DIRICHLET_BC) then
          temperature_at_datum = &
            condition%temperature%dataset%rarray(1)
          temperature_gradient(1:3) = &
            condition%temperature%gradient%rarray(1:3)
        endif
      endif
      if (associated(condition%concentration)) then
        if (condition%concentration%itype == DIRICHLET_BC .and. &
            associated(condition%concentration%dataset)) then
            concentration_at_datum = &
              condition%concentration%dataset%rarray(1)
            concentration_gradient(1:3) = &
              condition%concentration%gradient%rarray(1:3)
        else
          concentration_at_datum = -999.d0
          concentration_gradient = 0.d0
        endif
      endif

      if (associated(condition%datum)) then
        nullify(datum_dataset)
        select type(dataset=>condition%datum)
          class is (dataset_ascii_type)
            datum = dataset%rarray(1:3)
          class is (dataset_gridded_type)
            datum_dataset => dataset
            !TODO(geh): move this to FlowSubConditionUpdateDataset()
            !call DatasetLoad(datum_dataset,option)
            ! set datum here equal to estimated mid value of dataset
            datum(1:3) = -999.d0
            datum_dataset_rmax = maxval(datum_dataset%rarray)
            datum_dataset_rmin = minval(datum_dataset%rarray)
            datum(3) = 0.5d0*(datum_dataset_rmax+datum_dataset_rmin)
            ! round the number to the nearest whole number
            datum(3) = int(datum(3))
          class default
            option%io_buffer = &
              'Incorrect dataset type in HydrostaticUpdateCoupler'
            call printErrMsg(option)
        end select
      endif
      pressure_at_datum = &
        condition%pressure%dataset%rarray(1)
      ! gradient is in m/m; needs conversion to Pa/m
      piezometric_head_gradient(1:3) = &
        condition%pressure%gradient%rarray(1:3)
  end select      
      
  call nacl_den(temperature_at_datum,pressure_at_datum*1.d-6,xm_nacl,dw_kg) 
  rho = dw_kg * 1.d3
  
  gravity_magnitude = sqrt(DotProduct(option%gravity,option%gravity))
  
  if (dabs(gravity_magnitude-9.8068d0) > 0.1d0) then
    option%io_buffer = 'Magnitude of gravity vector is not near 9.81.'
    call printErrMsg(option)
  endif
  
  ! if a pressure gradient is prescribed in Z (at this point it will be a
  ! piezometric head gradient), the units of the pressure gradient are
  ! Pa/m and the pressure gradient does not need conversion
  if (dabs(piezometric_head_gradient(Z_DIRECTION)) < 1.d-40) then
    pressure_gradient(1:3) = piezometric_head_gradient(1:3)* &
                             rho*gravity_magnitude
  else
    pressure_gradient(1:3) = piezometric_head_gradient(1:3)
  endif

  if (dabs(pressure_gradient(Z_DIRECTION)) < 1.d-40) then
    ! compute the vertical gradient based on a 1 meter vertical spacing and
    ! interpolate the values from that array
    if (associated(datum_dataset)) then
      temp_real = grid%z_max_global - grid%z_min_global
      max_z = max(grid%z_max_global,datum_dataset_rmax+temp_real)+1.d0
      min_z = min(grid%z_min_global,datum_dataset_rmin-temp_real)-1.d0
    else
      max_z = max(grid%z_max_global,datum(Z_DIRECTION))+1.d0 ! add 1m buffer
      min_z = min(grid%z_min_global,datum(Z_DIRECTION))-1.d0
    endif
    
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
      select case(option%iflowmode)
        case(TH_MODE,THC_MODE,MPH_MODE,IMS_MODE,FLASH2_MODE,G_MODE, MIS_MODE)
          temperature = temperature + temperature_gradient(Z_DIRECTION)*delta_z
      end select
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
    select case(option%iflowmode)
      case(TH_MODE,THC_MODE,MPH_MODE,IMS_MODE,FLASH2_MODE,MIS_MODE,G_MODE)
        temperature = temperature_at_datum
    end select
    dist_z = 0.d0
    rho0 = density_array(idatum)
    do ipressure=idatum-1,1,-1
      dist_z = dist_z + delta_z
      select case(option%iflowmode)
        case(TH_MODE,THC_MODE,MPH_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          temperature = temperature - temperature_gradient(Z_DIRECTION)*delta_z
      end select
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

  if (grid%itype==STRUCTURED_GRID_MIMETIC .or. &
      grid%discretization_itype==UNSTRUCTURED_GRID_MIMETIC) then
      num_faces = coupler%numfaces_set
  else 
      num_faces = coupler%connection_set%num_connections
  end if


  do iconn=1, num_faces !geh: this should really be num_faces!
    if (grid%itype==STRUCTURED_GRID_MIMETIC .or. &
        grid%discretization_itype==UNSTRUCTURED_GRID_MIMETIC) then
#ifdef DASVYAT
      face_id_ghosted = coupler%faces_set(iconn)
      conn_set_ptr => grid%faces(face_id_ghosted)%conn_set_ptr
      conn_id = grid%faces(face_id_ghosted)%id

      dist_x = conn_set_ptr%cntr(1,conn_id) - datum(X_DIRECTION)
      dist_y = conn_set_ptr%cntr(2,conn_id) - datum(Y_DIRECTION)
      dist_z = conn_set_ptr%cntr(3,conn_id) - datum(Z_DIRECTION)
      z_offset = 0.d0
#endif
    else 

      local_id = coupler%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
  
      ! geh: note that this is a boundary connection, thus the entire distance is between
      ! the face and cell center
      if (associated(coupler%connection_set%dist)) then
        dx_conn = coupler%connection_set%dist(0,iconn)*coupler%connection_set%dist(1,iconn)
        dy_conn = coupler%connection_set%dist(0,iconn)*coupler%connection_set%dist(2,iconn)
        dz_conn = coupler%connection_set%dist(0,iconn)*coupler%connection_set%dist(3,iconn)
      endif
      if (associated(datum_dataset)) then
        ! correct datum based on dataset value
        ! if we interpolate in x and y, then we can use grid%x/y - dx/y_conn for x and y
        ! then we set dist_x and dist_y = 0.
        dist_x = 0.d0
        dist_y = 0.d0
        !TODO(geh): check that sign is correct for dx/y_conn
        call DatasetGriddedInterpolateReal(datum_dataset, &
                                           grid%x(ghosted_id)-dx_conn, &
                                           grid%y(ghosted_id)-dy_conn, &
                                           0.d0, &
                                           0.d0,temp_real,option)
        ! temp_real is now the real datum
        dist_z = grid%z(ghosted_id)-dz_conn-temp_real
        z_offset = temp_real-datum(Z_DIRECTION)
      else
        ! note the negative (-) d?_conn is required due to the offset of the boundary face
        dist_x = grid%x(ghosted_id)-dx_conn-datum(X_DIRECTION)
        dist_y = grid%y(ghosted_id)-dy_conn-datum(Y_DIRECTION)
        dist_z = grid%z(ghosted_id)-dz_conn-datum(Z_DIRECTION)
        z_offset = 0.d0
      endif
    end if


    if (associated(pressure_array)) then
      ipressure = idatum+int(dist_z/delta_z)
      if (grid%itype==STRUCTURED_GRID_MIMETIC.or. &
          grid%discretization_itype==UNSTRUCTURED_GRID_MIMETIC) then
        dist_z_for_pressure = conn_set_ptr%cntr(3,conn_id)-(z(ipressure) + z_offset)
      else 
        dist_z_for_pressure = grid%z(ghosted_id)-dz_conn-(z(ipressure) + z_offset)
      end if
      pressure = pressure_array(ipressure) + &
                 density_array(ipressure)*option%gravity(Z_DIRECTION) * &
                 dist_z_for_pressure + &
!                 (grid%z(ghosted_id)-z(ipressure)) + &
                 pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 pressure_gradient(Y_DIRECTION)*dist_y

!      if (grid%itype==STRUCTURED_GRID_MIMETIC) then
!         pressure = 3*conn_set_ptr%cntr(3,conn_id)    !DASVYAT WORKINGCHECK
!      end if
 
    else
      pressure = pressure_at_datum + &
                 pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 pressure_gradient(Y_DIRECTION)*dist_y + &
                 pressure_gradient(Z_DIRECTION)*dist_z 
    endif
   
 
    if (pressure < option%minimum_hydrostatic_pressure) &
      pressure = option%minimum_hydrostatic_pressure

    ! assign pressure
    select case(option%iflowmode)
      case(G_MODE)
        coupler%flow_aux_real_var(GENERAL_LIQUID_PRESSURE_DOF,iconn) = pressure
      case (MPH_MODE)
        coupler%flow_aux_real_var(1,iconn) = pressure
      case default
        if (condition%pressure%itype == SEEPAGE_BC) then
          coupler%flow_aux_real_var(1,iconn) = max(pressure,option%reference_pressure)
        else if (condition%pressure%itype == CONDUCTANCE_BC) then
          coupler%flow_aux_real_var(1,iconn) = max(pressure,option%reference_pressure)
          coupler%flow_aux_real_var(2,iconn) = condition%pressure%aux_real(1)
        else
          coupler%flow_aux_real_var(1,iconn) = pressure
        endif
    end select

    ! assign other dofs
    select case(option%iflowmode)
      case(THC_MODE,MPH_MODE,IMS_MODE,FLASH2_MODE)
        temperature = temperature_at_datum + &
                    temperature_gradient(X_DIRECTION)*dist_x + & ! gradient in K/m
                    temperature_gradient(Y_DIRECTION)*dist_y + &
                    temperature_gradient(Z_DIRECTION)*dist_z 
        coupler%flow_aux_real_var(2,iconn) = temperature
        coupler%flow_aux_real_var(3,iconn) = concentration_at_datum

        coupler%flow_aux_int_var(1,iconn) = condition%iphase

      case(TH_MODE)
        temperature = temperature_at_datum + &
                    temperature_gradient(X_DIRECTION)*dist_x + & ! gradient in K/m
                    temperature_gradient(Y_DIRECTION)*dist_y + &
                    temperature_gradient(Z_DIRECTION)*dist_z
        coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = temperature
        coupler%flow_aux_int_var(TH_PRESSURE_DOF,iconn) = condition%iphase

      case(MIS_MODE)
        temperature = temperature_at_datum + &
                    temperature_gradient(X_DIRECTION)*dist_x + & ! gradient in K/m
                    temperature_gradient(Y_DIRECTION)*dist_y + &
                    temperature_gradient(Z_DIRECTION)*dist_z 
!       coupler%flow_aux_real_var(2,iconn) = temperature
        coupler%flow_aux_real_var(2,iconn) = concentration_at_datum

        coupler%flow_aux_int_var(1,iconn) = condition%iphase

      case(G_MODE)
        temperature = temperature_at_datum + &
                    temperature_gradient(X_DIRECTION)*dist_x + & ! gradient in K/m
                    temperature_gradient(Y_DIRECTION)*dist_y + &
                    temperature_gradient(Z_DIRECTION)*dist_z 
        coupler%flow_aux_real_var(GENERAL_LIQUID_STATE_TEMPERATURE_DOF,iconn) = &
          temperature
        coupler%flow_aux_real_var(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF,iconn) = &
          concentration_at_datum

        coupler%flow_aux_int_var(GENERAL_LIQUID_PRESSURE_DOF,iconn) = condition%iphase
      case default
        coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,iconn) = 1
    end select

  enddo

  if ((grid%itype==STRUCTURED_GRID_MIMETIC.or. &
       grid%discretization_itype==UNSTRUCTURED_GRID_MIMETIC).and. &
      (coupler%itype == INITIAL_COUPLER_TYPE)) then
     num_regions = coupler%region%num_cells
     do iconn = 1, num_regions

       local_id = coupler%region%cell_ids(iconn)
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
         dist_z_for_pressure = grid%z(ghosted_id)-dz_conn-z(ipressure)
         pressure = pressure_array(ipressure) + &
                 density_array(ipressure)*option%gravity(Z_DIRECTION) * &
                 dist_z_for_pressure + &
!                 (grid%z(ghosted_id)-z(ipressure)) + &
                 pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 pressure_gradient(Y_DIRECTION)*dist_y

 
       else
         pressure = pressure_at_datum + &
                 pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 pressure_gradient(Y_DIRECTION)*dist_y + &
                 pressure_gradient(Z_DIRECTION)*dist_z 
       endif
   
 
      if (pressure < option%minimum_hydrostatic_pressure) &
          pressure = option%minimum_hydrostatic_pressure

      if (condition%pressure%itype == SEEPAGE_BC) then
        coupler%flow_aux_real_var(1,num_faces + iconn) = max(pressure,option%reference_pressure)
      else if (condition%pressure%itype == CONDUCTANCE_BC) then
        coupler%flow_aux_real_var(1,num_faces + iconn) = max(pressure,option%reference_pressure)
        coupler%flow_aux_real_var(2,num_faces + iconn) = condition%pressure%aux_real(1)
      else
        coupler%flow_aux_real_var(1,num_faces + iconn) = pressure
      endif

      select case(option%iflowmode)
        case(TH_MODE,THC_MODE,MPH_MODE,IMS_MODE,FLASH2_MODE)
           temperature = temperature_at_datum + &
                      temperature_gradient(X_DIRECTION)*dist_x + & ! gradient in K/m
                      temperature_gradient(Y_DIRECTION)*dist_y + &
                     temperature_gradient(Z_DIRECTION)*dist_z 
           coupler%flow_aux_real_var(2,num_faces + iconn) = temperature
           coupler%flow_aux_real_var(3,num_faces + iconn) = concentration_at_datum
        case(G_MODE)
      end select

       coupler%flow_aux_int_var(1,num_faces + iconn) = 1
     end do 
  end if

!   write(*,*) "End of HydrostaticUpdateCoupler" 
!   read(*,*)
  if (associated(pressure_array)) deallocate(pressure_array)
  nullify(pressure_array)
  if (allocated(z)) deallocate(z)
  if (allocated(density_array)) deallocate(density_array)
  !geh: Do not deallocate datum_dataset as it is soleley a pointer to an
  !     external dataset.
  nullify(datum_dataset)

end subroutine HydrostaticUpdateCoupler

! ************************************************************************** !
!
! HydrostaticTest: Computes the hydrostatic initial/boundary 
!                                 condition (more accurately than before0
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
subroutine HydrostaticTest()

  use Water_EOS_module
  
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

#if 0
! ************************************************************************** !
!
! ProjectAOntoUnitB: Projects vector a onto b, assuming b is a unit vector
! author: Glenn Hammond
! date: 02/20/09
!
! ************************************************************************** !
function ProjectAOntoUnitB(A,B)

  implicit none
  
  PetscReal :: A(3)
  PetscReal :: B(3)
  
  PetscReal :: ProjectAOntoUnitB(3)
  
  ProjectAOntoUnitB = dot_product(A,B)*A
  
end function ProjectAOntoUnitB
#endif

end module Hydrostatic_module

