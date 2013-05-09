module Logging_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"

! stages
PetscInt, parameter, public :: INIT_STAGE = 1
PetscInt, parameter, public :: TS_STAGE = 2
PetscInt, parameter, public :: FLOW_STAGE = 3
PetscInt, parameter, public :: TRAN_STAGE = 4
PetscInt, parameter, public :: OUTPUT_STAGE = 5

  type, public :: logging_type 
  
    PetscLogStage :: stage(10)
    
    PetscClassId :: class_pflotran
    
    PetscLogEvent :: event_init
    PetscLogEvent :: event_setup
   
    PetscLogEvent :: event_create_iogroups

    PetscLogEvent :: event_restart
    PetscLogEvent :: event_checkpoint

    PetscLogEvent :: event_flow_condition_read
    PetscLogEvent :: event_tran_condition_read
    PetscLogEvent :: event_tran_constraint_read
    PetscLogEvent :: event_flow_condition_read_values

    PetscLogEvent :: event_h5dread_f
    PetscLogEvent :: event_h5dwrite_f
    PetscLogEvent :: event_read_indices_hdf5
    PetscLogEvent :: event_map_indices_hdf5
    PetscLogEvent :: event_hash_create
    PetscLogEvent :: event_hash_map
    PetscLogEvent :: event_read_real_array_hdf5
    PetscLogEvent :: event_read_ndim_real_array_hdf5
    PetscLogEvent :: event_read_int_array_hdf5
    PetscLogEvent :: event_write_real_array_hdf5
    PetscLogEvent :: event_write_int_array_hdf5
    PetscLogEvent :: event_read_array_hdf5    
    PetscLogEvent :: event_read_xyz_dataset_hdf5    
    PetscLogEvent :: event_write_struct_dataset_hdf5
    PetscLogEvent :: event_region_read_hdf5
    PetscLogEvent :: event_region_read_ascii
    PetscLogEvent :: event_cell_indx_int_read_hdf5
    PetscLogEvent :: event_cell_indx_real_read_hdf5

    PetscLogEvent :: event_output_tecplot
    PetscLogEvent :: event_output_hdf5
    PetscLogEvent :: event_output_vtk
    PetscLogEvent :: event_output_grid_vtk
    PetscLogEvent :: event_output_write_vtk
    PetscLogEvent :: event_output_mad
    PetscLogEvent :: event_output_str_grid_tecplot
    PetscLogEvent :: event_output_write_tecplot
    PetscLogEvent :: event_output_write_flux_tecplot
    PetscLogEvent :: event_output_get_var_from_array
    PetscLogEvent :: event_output_get_cell_vel
    PetscLogEvent :: event_output_vec_tecplot
    PetscLogEvent :: event_output_observation
    PetscLogEvent :: event_output_coordinates_hdf5
    PetscLogEvent :: event_output_hydrograph
    PetscLogEvent :: event_output_secondary_tecplot

    PetscLogEvent :: event_r_residual
    PetscLogEvent :: event_r_jacobian
    PetscLogEvent :: event_r_auxvars
    PetscLogEvent :: event_r_auxvars_bc
    
    PetscLogEvent :: event_rt_residual
    PetscLogEvent :: event_rt_jacobian

    PetscLogEvent :: event_rt_jacobian_flux
    PetscLogEvent :: event_rt_jacobian_fluxbc
    PetscLogEvent :: event_rt_jacobian_accum
    PetscLogEvent :: event_rt_jacobian_zero_calc
    PetscLogEvent :: event_rt_jacobian_zero
    PetscLogEvent :: event_rt_jacobian_ss
    PetscLogEvent :: event_rt_jacobian1
    PetscLogEvent :: event_rt_jacobian2

    PetscLogEvent :: event_rt_res_reaction
    PetscLogEvent :: event_rt_jac_reaction
    PetscLogEvent :: event_rt_react
    PetscLogEvent :: event_rt_auxvars
    PetscLogEvent :: event_rt_auxvars_bc
    
    PetscLogEvent :: event_mass_balance

#ifdef DASVYAT

    PetscLogEvent :: event_flow_residual_mfd1
    PetscLogEvent :: event_flow_residual_mfd2
    PetscLogEvent :: event_flow_flux_mfd
    PetscLogEvent :: event_flow_rhs_mfd

#endif

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
  call PetscLogStageRegister('Flow Stage', &
                             logging%stage(FLOW_STAGE),ierr)
  call PetscLogStageRegister('Transport Stage', &
                             logging%stage(TRAN_STAGE),ierr)
  call PetscLogStageRegister('Output Stage', &
                             logging%stage(OUTPUT_STAGE),ierr)
                             
!!  call PetscCookieRegister('PFLOTRAN',logging%class_pflotran,ierr)
  call PetscClassIdRegister('PFLOTRAN',logging%class_pflotran,ierr)

  call PetscLogEventRegister('Init', &
                             logging%class_pflotran, &
                             logging%event_init,ierr)
  call PetscLogEventRegister('Init,Setup', &
                             logging%class_pflotran, &
                             logging%event_setup,ierr)

  call PetscLogEventRegister('Create_iogroups', &
                             logging%class_pflotran, &
                             logging%event_create_iogroups,ierr)

  call PetscLogEventRegister('Restart', &
                             logging%class_pflotran, &
                             logging%event_restart,ierr)
  call PetscLogEventRegister('Checkpoint', &
                             logging%class_pflotran, &
                             logging%event_checkpoint,ierr)

  call PetscLogEventRegister('FlowCondRead', &
                             logging%class_pflotran, &
                             logging%event_flow_condition_read,ierr)
  call PetscLogEventRegister('TranCondRead', &
                             logging%class_pflotran, &
                             logging%event_tran_condition_read,ierr)
  call PetscLogEventRegister('TranConstraintRd', &
                             logging%class_pflotran, &
                             logging%event_tran_constraint_read,ierr)
  call PetscLogEventRegister('FlowCondReadVals', &
                             logging%class_pflotran, &
                             logging%event_flow_condition_read_values,ierr)

  call PetscLogEventRegister('H5DRead_F', &
                             logging%class_pflotran, &
                             logging%event_h5dread_f,ierr)
  call PetscLogEventRegister('H5DWrite_F', &
                             logging%class_pflotran, &
                             logging%event_h5dwrite_f,ierr)
  call PetscLogEventRegister('HDF5ReadIndices', &
                             logging%class_pflotran, &
                             logging%event_read_indices_hdf5,ierr)
  call PetscLogEventRegister('H5MapLoc2NatIndx', &
                             logging%class_pflotran, &
                             logging%event_map_indices_hdf5,ierr)
  call PetscLogEventRegister('GrdCrNat2GhstHsh', &
                             logging%class_pflotran, &
                             logging%event_hash_create,ierr)
  call PetscLogEventRegister('GrdLocGhstIdHsh', &
                             logging%class_pflotran, &
                             logging%event_hash_map,ierr)
  call PetscLogEventRegister('H5ReadRealArray', &
                             logging%class_pflotran, &
                             logging%event_read_real_array_hdf5,ierr)
  call PetscLogEventRegister('H5ReadNDimRealArray', &
                             logging%class_pflotran, &
                             logging%event_read_ndim_real_array_hdf5,ierr)
  call PetscLogEventRegister('H5ReadIntArray', &
                             logging%class_pflotran, &
                             logging%event_read_int_array_hdf5,ierr)
  call PetscLogEventRegister('H5ReadArray', &
                             logging%class_pflotran, &
                             logging%event_read_array_hdf5,ierr)
  call PetscLogEventRegister('H5WriteRealArray', &
                             logging%class_pflotran, &
                             logging%event_write_real_array_hdf5,ierr)
  call PetscLogEventRegister('H5WriteIntArray', &
                             logging%class_pflotran, &
                             logging%event_write_int_array_hdf5,ierr)
  call PetscLogEventRegister('H5WriteStrData', &
                             logging%class_pflotran, &
                             logging%event_write_struct_dataset_hdf5,ierr)
  call PetscLogEventRegister('H5ReadRegFrmFile', &
                             logging%class_pflotran, &
                             logging%event_region_read_hdf5,ierr)
  call PetscLogEventRegister('RegReadFrmFileId', &
                             logging%class_pflotran, &
                             logging%event_region_read_ascii,ierr)
  call PetscLogEventRegister('H5RdCellIndxInt', &
                             logging%class_pflotran, &
                             logging%event_cell_indx_int_read_hdf5,ierr)
  call PetscLogEventRegister('H5RdCellIndxReal', &
                             logging%class_pflotran, &
                             logging%event_cell_indx_real_read_hdf5,ierr)

  call PetscLogEventRegister('OutputTecplot', &
                             logging%class_pflotran, &
                             logging%event_output_tecplot,ierr)
  call PetscLogEventRegister('OutputVTK', &
                             logging%class_pflotran, &
                             logging%event_output_vtk,ierr)
  call PetscLogEventRegister('OutputMAD', &
                             logging%class_pflotran, &
                             logging%event_output_mad,ierr)
  call PetscLogEventRegister('OutputHDF5', &
                             logging%class_pflotran, &
                             logging%event_output_hdf5,ierr)
  call PetscLogEventRegister('OutputSecondaryTecplot', &
                             logging%class_pflotran, &
                             logging%event_output_secondary_tecplot,ierr)
  call PetscLogEventRegister('WriteTecStrGrid', &
                             logging%class_pflotran, &
                             logging%event_output_str_grid_tecplot,ierr)
  call PetscLogEventRegister('WriteVTKGrid', &
                             logging%class_pflotran, &
                             logging%event_output_grid_vtk,ierr)
  call PetscLogEventRegister('WriteTecDataSet', &
                             logging%class_pflotran, &
                             logging%event_output_write_tecplot,ierr)
  call PetscLogEventRegister('WriteVTKDataSet', &
                             logging%class_pflotran, &
                             logging%event_output_write_vtk,ierr)
  call PetscLogEventRegister('OutputFluxVelTec', &
                             logging%class_pflotran, &
                             logging%event_output_write_flux_tecplot,ierr)
  call PetscLogEventRegister('OutputGtVrFrmArr', &
                             logging%class_pflotran, &
                             logging%event_output_get_var_from_array,ierr)
  call PetscLogEventRegister('GetCellCentVel', &
                             logging%class_pflotran, &
                             logging%event_output_get_cell_vel,ierr)
  call PetscLogEventRegister('OutputFluxVelTec', &
                             logging%class_pflotran, &
                             logging%event_output_vec_tecplot,ierr)
  call PetscLogEventRegister('OutputBrkthuTec', &
                             logging%class_pflotran, &
                             logging%event_output_observation,ierr)
  call PetscLogEventRegister('WriteHDF5Coord', &
                             logging%class_pflotran, &
                             logging%event_output_coordinates_hdf5,ierr)
  call PetscLogEventRegister('OutputHydrograph', &
                             logging%class_pflotran, &
                             logging%event_output_hydrograph,ierr)

  call PetscLogEventRegister('RResidual', &
                             logging%class_pflotran, &
                             logging%event_r_residual,ierr)
  call PetscLogEventRegister('RJacobian', &
                             logging%class_pflotran, &
                             logging%event_r_jacobian,ierr)
  call PetscLogEventRegister('RAuxVars', &
                             logging%class_pflotran, &
                             logging%event_r_auxvars,ierr)
  call PetscLogEventRegister('RAuxVarsBC', &
                             logging%class_pflotran, &
                             logging%event_r_auxvars_bc,ierr)

  call PetscLogEventRegister('RTResidual', &
                             logging%class_pflotran, &
                             logging%event_rt_residual,ierr)
  call PetscLogEventRegister('RTJacobian', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian,ierr)

  call PetscLogEventRegister('RTJacobianFlux', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_flux,ierr)
  call PetscLogEventRegister('RTJacobianFluxBC', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_fluxbc,ierr)
  call PetscLogEventRegister('RTJacobianAccum', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_accum,ierr)
  call PetscLogEventRegister('RTJacobianZeroCalc', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_zero_calc,ierr)
  call PetscLogEventRegister('RTJacobianZero', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_zero,ierr)
  call PetscLogEventRegister('RTJacobianSS', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_ss,ierr)

  call PetscLogEventRegister('RTJacobian1', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian1,ierr)
  call PetscLogEventRegister('RTJacobian2', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian2,ierr)

  call PetscLogEventRegister('RTResReaction', &
                             logging%class_pflotran, &
                             logging%event_rt_res_reaction,ierr)
  call PetscLogEventRegister('RTJacReaction', &
                             logging%class_pflotran, &
                             logging%event_rt_jac_reaction,ierr)
  call PetscLogEventRegister('RTReact', &
                             logging%class_pflotran, &
                             logging%event_rt_react,ierr)
  call PetscLogEventRegister('RTAuxVars', &
                             logging%class_pflotran, &
                             logging%event_rt_auxvars,ierr)
  call PetscLogEventRegister('RTAuxVarsBC', &
                             logging%class_pflotran, &
                             logging%event_rt_auxvars_bc,ierr)
                             
  call PetscLogEventRegister('MassBalance', &
                             logging%class_pflotran, &
                             logging%event_mass_balance,ierr)

#ifdef DASVYAT
  
  call PetscLogEventRegister('MFD_Residual1', &
                             logging%class_pflotran, &
                             logging%event_flow_residual_mfd1, ierr)
  
  
  call PetscLogEventRegister('MFD_Residual2', &
                             logging%class_pflotran, &
                             logging%event_flow_residual_mfd2, ierr)

  call PetscLogEventRegister('MFD_Flux', &
                             logging%class_pflotran, &
                             logging%event_flow_flux_mfd, ierr)
  call PetscLogEventRegister('MFD_Rhs', &
                             logging%class_pflotran, &
                             logging%event_flow_rhs_mfd, ierr)
#endif
  
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
