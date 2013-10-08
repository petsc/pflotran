module Saturation_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: SaturationUpdateCoupler
 
contains

! ************************************************************************** !
!
! SaturationUpdateCoupler: Computes the pressures for a saturation 
!                          initial/boundary condition 
! author: Glenn Hammond
! date: 06/07/11
!
! ************************************************************************** !
subroutine SaturationUpdateCoupler(coupler,option,grid,saturation_functions, &
                                   sat_func_id)

  use Water_EOS_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Region_module
  use Saturation_Function_module

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  type(saturation_function_ptr_type) :: saturation_functions(:)
  PetscInt :: sat_func_id(:)

  PetscInt :: local_id, ghosted_id, iconn
  PetscReal :: saturation
  PetscReal :: capillary_pressure
  PetscReal :: liquid_pressure
  
  type(flow_condition_type), pointer :: condition
  
  type(connection_set_type), pointer :: cur_connection_set
  
  condition => coupler%flow_condition

  ! in this case, the saturation is stored within concentration dataset
  saturation = condition%saturation%flow_dataset%time_series%cur_value(1)

  do iconn = 1, coupler%connection_set%num_connections
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    call SatFuncGetCapillaryPressure(capillary_pressure,saturation, &
                     saturation_functions(sat_func_id(ghosted_id))%ptr,option)
    liquid_pressure = option%reference_pressure - capillary_pressure
    coupler%flow_aux_real_var(1,iconn) = liquid_pressure
  enddo

#if 0  
  if (grid%itype==STRUCTURED_GRID_MIMETIC) then 
      num_faces = coupler%numfaces_set
  else 
      num_faces = coupler%connection_set%num_connections
  end if

  do iconn=1, num_faces
    if (grid%itype==STRUCTURED_GRID_MIMETIC) then
#ifdef DASVYAT
      face_id_ghosted = coupler%faces_set(iconn)
      conn_set_ptr => grid%faces(face_id_ghosted)%conn_set_ptr
      conn_id = grid%faces(face_id_ghosted)%id

      dist_x = conn_set_ptr%cntr(1,conn_id) - datum(X_DIRECTION)
      dist_y = conn_set_ptr%cntr(2,conn_id) - datum(Y_DIRECTION)
      dist_z = conn_set_ptr%cntr(3,conn_id) - datum(Z_DIRECTION)
#endif
    else 

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
    end if


    if (associated(pressure_array)) then
      ipressure = idatum+int(dist_z/delta_z)
      if (grid%itype==STRUCTURED_GRID_MIMETIC) then
        dist_z_for_pressure = conn_set_ptr%cntr(3,conn_id) - z(ipressure)
      else 
        dist_z_for_pressure = grid%z(ghosted_id)-dz_conn-z(ipressure)
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

    if (condition%pressure%itype == SEEPAGE_BC) then
      coupler%flow_aux_real_var(1,iconn) = max(pressure,option%reference_pressure)
    else if (condition%pressure%itype == CONDUCTANCE_BC) then
      coupler%flow_aux_real_var(1,iconn) = max(pressure,option%reference_pressure)
      coupler%flow_aux_real_var(2,iconn) = condition%pressure%aux_real(1)
    else
      coupler%flow_aux_real_var(1,iconn) = pressure
    endif

    select case(option%iflowmode)
      case(TH_MODE,THC_MODE,MPH_MODE,IMS_MODE,FLASH2_MODE)
        temperature = temperature_at_datum + &
                    temperature_gradient(X_DIRECTION)*dist_x + & ! gradient in K/m
                    temperature_gradient(Y_DIRECTION)*dist_y + &
                    temperature_gradient(Z_DIRECTION)*dist_z 
        coupler%flow_aux_real_var(2,iconn) = temperature
        coupler%flow_aux_real_var(3,iconn) = concentration_at_datum

        coupler%flow_aux_int_var(1,iconn) = condition%iphase

      case(G_MODE)
      case default
        coupler%flow_aux_int_var(1,iconn) = 1
    end select

  enddo

  if ((grid%itype==STRUCTURED_GRID_MIMETIC).and.(coupler%itype == INITIAL_COUPLER_TYPE)) then
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
#endif
end subroutine SaturationUpdateCoupler

end module Saturation_module

