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
    
    PetscInt :: class_pflotran
    
    PetscEvent :: event_init
    PetscEvent :: event_setup

    PetscEvent :: event_restart
    PetscEvent :: event_checkpoint

    PetscEvent :: event_condition_read
    PetscEvent :: event_condition_read_values

    PetscEvent :: event_h5dread_f
    PetscEvent :: event_h5dwrite_f
    PetscEvent :: event_read_indices_hdf5
    PetscEvent :: event_map_indices_hdf5
    PetscEvent :: event_hash_create
    PetscEvent :: event_hash_map
    PetscEvent :: event_read_real_array_hdf5
    PetscEvent :: event_read_int_array_hdf5
    PetscEvent :: event_write_real_array_hdf5
    PetscEvent :: event_write_int_array_hdf5
    PetscEvent :: event_read_array_hdf5    
    PetscEvent :: event_write_struct_dataset_hdf5
    PetscEvent :: event_region_read_hdf5
    PetscEvent :: event_region_read_ascii
    PetscEvent :: event_material_read_hdf5

    PetscEvent :: event_output_tecplot
    PetscEvent :: event_output_hdf5
    PetscEvent :: event_output_str_grid_tecplot
    PetscEvent :: event_output_write_tecplot
    PetscEvent :: event_output_write_flux_tecplot
    PetscEvent :: event_output_get_var_from_array
    PetscEvent :: event_output_get_cell_vel
    PetscEvent :: event_output_vec_tecplot
    PetscEvent :: event_output_breakthrough
    PetscEvent :: event_output_coordinates_hdf5
    
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

  call PetscLogStageRegister(logging%stage(INIT_STAGE), &
                             "Init Stage",ierr)
  call PetscLogStageRegister(logging%stage(TS_STAGE), &
                             "Time Step Stage",ierr)
  call PetscLogStageRegister(logging%stage(OUTPUT_STAGE), &
                             "Output Stage",ierr)
                             
  call PetscLogClassRegister(logging%class_pflotran,'PFLTORAN',ierr)

  call PetscLogEventRegister(logging%event_init, &
                             'Init', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_setup, &
                             'Init,Setup', &
                             logging%class_pflotran,ierr)

  call PetscLogEventRegister(logging%event_restart, &
                             'Restart', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_checkpoint, &
                             'Checkpoint', &
                             logging%class_pflotran,ierr)

  call PetscLogEventRegister(logging%event_condition_read, &
                             'ConditionRead', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_condition_read_values, &
                             'ConditionRdVals', &
                             logging%class_pflotran,ierr)

  call PetscLogEventRegister(logging%event_h5dread_f, &
                             'H5DRead_F', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_h5dwrite_f, &
                             'H5DWrite_F', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_read_indices_hdf5, &
                             'HDF5ReadIndices', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_map_indices_hdf5, &
                             'H5MapLoc2NatIndx', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_hash_create, &
                             'GrdCrNat2GhstHsh', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_hash_map, &
                             'GrdLocGhstIdHsh', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_read_real_array_hdf5, &
                             'H5ReadRealArray', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_read_int_array_hdf5, &
                             'H5ReadIntArray', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_read_array_hdf5, &
                             'H5ReadArray', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_write_real_array_hdf5, &
                             'H5WriteRealArray', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_write_int_array_hdf5, &
                             'H5WriteIntArray', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_write_struct_dataset_hdf5, &
                             'H5WriteStrData', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_region_read_hdf5, &
                             'H5ReadRegFrmFile', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_region_read_ascii, &
                             'RegReadFrmFileId', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_material_read_hdf5, &
                             'H5ReadMatFrmFile', &
                             logging%class_pflotran,ierr)

  call PetscLogEventRegister(logging%event_output_tecplot, &
                             'OutputTecplot', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_hdf5, &
                             'OutputHDF5', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_str_grid_tecplot, &
                             'WriteTecStrGrid', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_write_tecplot, &
                             'WriteTecDataSet', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_write_flux_tecplot, &
                             'OutputFluxVelTec', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_get_var_from_array, &
                             'OutputGtVrFrmArr', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_get_cell_vel, &
                             'GetCellCentVel', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_vec_tecplot, &
                             'OutputFluxVelTec', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_breakthrough, &
                             'OutputBrkthuTec', &
                             logging%class_pflotran,ierr)
  call PetscLogEventRegister(logging%event_output_coordinates_hdf5, &
                             'WriteHDF5Coord', &
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
