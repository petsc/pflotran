module Logging_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"

! stages
PetscInt, parameter, public :: INIT_STAGE = 1
PetscInt, parameter, public :: TS_STAGE = 2
PetscInt, parameter, public :: OUTPUT_STAGE = 3

  type, public :: logging_type 
  
    PetscInt :: stage(10)
    
    PetscCookie :: class_pflotran
    
    PetscLogEvent :: event_init
    PetscLogEvent :: event_setup

    PetscLogEvent :: event_restart
    PetscLogEvent :: event_checkpoint

    PetscLogEvent :: event_condition_read
    PetscLogEvent :: event_condition_read_values

    PetscLogEvent :: event_h5dread_f
    PetscLogEvent :: event_h5dwrite_f
    PetscLogEvent :: event_read_indices_hdf5
    PetscLogEvent :: event_map_indices_hdf5
    PetscLogEvent :: event_hash_create
    PetscLogEvent :: event_hash_map
    PetscLogEvent :: event_read_real_array_hdf5
    PetscLogEvent :: event_read_int_array_hdf5
    PetscLogEvent :: event_write_real_array_hdf5
    PetscLogEvent :: event_write_int_array_hdf5
    PetscLogEvent :: event_read_array_hdf5    
    PetscLogEvent :: event_write_struct_dataset_hdf5
    PetscLogEvent :: event_region_read_hdf5
    PetscLogEvent :: event_region_read_ascii
    PetscLogEvent :: event_material_read_hdf5

    PetscLogEvent :: event_output_tecplot
    PetscLogEvent :: event_output_hdf5
    PetscLogEvent :: event_output_str_grid_tecplot
    PetscLogEvent :: event_output_write_tecplot
    PetscLogEvent :: event_output_write_flux_tecplot
    PetscLogEvent :: event_output_get_var_from_array
    PetscLogEvent :: event_output_get_cell_vel
    PetscLogEvent :: event_output_vec_tecplot
    PetscLogEvent :: event_output_breakthrough
    PetscLogEvent :: event_output_coordinates_hdf5
    
  end type logging_type
  
  type(logging_type), pointer, public :: logging
  
  public :: LoggingCreate, &
            LoggingDestroy

contains

! ************************************************************************** !
!
! LoggingCreate: Allocates and initializes a new logging object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine LoggingCreate()

  implicit none
  
  PetscErrorCode :: ierr
  
  allocate(logging)

  call PetscLogStageRegister('Init Stage',  & 
                             logging%stage(INIT_STAGE),ierr)
  call PetscLogStageRegister('Time Step Stage', &
                             logging%stage(TS_STAGE),ierr)
  call PetscLogStageRegister('Output Stage', &
                             logging%stage(OUTPUT_STAGE),ierr)
                             
  call PetscCookieRegister('PFLTORAN',logging%class_pflotran,ierr)

  call PetscLogEventRegister('Init', &
                             logging%event_init, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('Init,Setup', &
                             logging%event_setup, &
                             logging%class_pflotran,ierr)

  call PetscLogEventRegister('Restart', &
                             logging%event_restart, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('Checkpoint', &
                             logging%event_checkpoint, &
                             logging%class_pflotran,ierr)

  call PetscLogEventRegister('ConditionRead', &
                             logging%event_condition_read, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('ConditionRdVals', &
                             logging%event_condition_read_values, &
                             logging%class_pflotran,ierr)

  call PetscLogEventRegister('H5DRead_F', &
                             logging%event_h5dread_f, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5DWrite_F', &
                             logging%event_h5dwrite_f, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('HDF5ReadIndices', &
                             logging%event_read_indices_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5MapLoc2NatIndx', &
                             logging%event_map_indices_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('GrdCrNat2GhstHsh', &
                             logging%event_hash_create, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('GrdLocGhstIdHsh', &
                             logging%event_hash_map, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5ReadRealArray', &
                             logging%event_read_real_array_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5ReadIntArray', &
                             logging%event_read_int_array_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5ReadArray', &
                             logging%event_read_array_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5WriteRealArray', &
                             logging%event_write_real_array_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5WriteIntArray', &
                             logging%event_write_int_array_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5WriteStrData', &
                             logging%event_write_struct_dataset_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5ReadRegFrmFile', &
                             logging%event_region_read_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('RegReadFrmFileId', &
                             logging%event_region_read_ascii, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('H5ReadMatFrmFile', &
                             logging%event_material_read_hdf5, &
                             logging%class_pflotran,ierr)

  call PetscLogEventRegister('OutputTecplot', &
                             logging%event_output_tecplot, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('OutputHDF5', &
                             logging%event_output_hdf5, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('WriteTecStrGrid', &
                             logging%event_output_str_grid_tecplot, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('WriteTecDataSet', &
                             logging%event_output_write_tecplot, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('OutputFluxVelTec', &
                             logging%event_output_write_flux_tecplot, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('OutputGtVrFrmArr', &
                             logging%event_output_get_var_from_array, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('GetCellCentVel', &
                             logging%event_output_get_cell_vel, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('OutputFluxVelTec', &
                             logging%event_output_vec_tecplot, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('OutputBrkthuTec', &
                             logging%event_output_breakthrough, &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister('WriteHDF5Coord', &
                             logging%event_output_coordinates_hdf5, &
                             logging%class_pflotran,ierr)
  
end subroutine LoggingCreate

! ************************************************************************** !
!
! OptionDestroy: Deallocates a logging object
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine LoggingDestroy()

  implicit none
  
  ! all kinds of stuff needs to be added here.
  
  deallocate(logging)
  nullify(logging)
  
end subroutine LoggingDestroy

end module Logging_module
