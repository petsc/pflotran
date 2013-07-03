#ifdef GEOMECH

module Geomechanics_Logging_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"

! stages
PetscInt, parameter, public :: GEOMECH_INIT_STAGE = 1
PetscInt, parameter, public :: GEOMECH_TS_STAGE = 2
PetscInt, parameter, public :: GEOMECH_SOLVE_STAGE = 3
PetscInt, parameter, public :: GEOMECH_OUTPUT_STAGE = 4

  type, public :: geomech_logging_type 
  
    PetscLogStage :: stage(10)
    
    PetscClassId :: class_pflotran
 
    PetscLogEvent :: event_geomech_condition_read
    PetscLogEvent :: event_geomech_condition_read_values
    PetscLogEvent :: event_geomech_residual
    PetscLogEvent :: event_geomech_jacobian
    PetscLogEvent :: event_output_tecplot

  end type geomech_logging_type
  
  type(geomech_logging_type), pointer, public :: geomech_logging
  
  public :: GeomechLoggingCreate, &
            GeomechLoggingDestroy

contains

! ************************************************************************** !
!
! GeomechLoggingCreate: Allocates and initializes a new geomech_logging object
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechLoggingCreate()

  implicit none
  
  PetscErrorCode :: ierr
  
  allocate(geomech_logging)

  call PetscLogStageRegister('Geomech Init Stage',  & 
                             geomech_logging%stage(GEOMECH_INIT_STAGE),ierr)
  call PetscLogStageRegister('Geomech Time Step Stage', &
                             geomech_logging%stage(GEOMECH_TS_STAGE),ierr)
  call PetscLogStageRegister('Geomech Flow Stage', &
                             geomech_logging%stage(GEOMECH_SOLVE_STAGE),ierr)
  call PetscLogStageRegister('Geomech Output Stage', &
                             geomech_logging%stage(GEOMECH_OUTPUT_STAGE),ierr)
                             
  call PetscClassIdRegister('Geomech PFLOTRAN',geomech_logging%class_pflotran,ierr)

  call PetscLogEventRegister('GeomechCondRead', &
                             geomech_logging%class_pflotran, &
                             geomech_logging%event_geomech_condition_read,ierr)
 
  call PetscLogEventRegister('GeomechCondReadVals', &
                             geomech_logging%class_pflotran, &
                             geomech_logging%event_geomech_condition_read_values,ierr)

  call PetscLogEventRegister('GeomechResidual', &
                             geomech_logging%class_pflotran, &
                             geomech_logging%event_geomech_residual,ierr)
                            
  call PetscLogEventRegister('GeomechJacobian', &
                             geomech_logging%class_pflotran, &
                             geomech_logging%event_geomech_jacobian,ierr)
                             
  call PetscLogEventRegister('GeomechOutputTecplot', &
                             geomech_logging%class_pflotran, &
                             geomech_logging%event_output_tecplot,ierr)
  
end subroutine GeomechLoggingCreate

! ************************************************************************** !
!
! OptionDestroy: Deallocates a geomech_logging object
! author: Satish Karra, LANL
! date: 06/12/13
!
! ************************************************************************** !
subroutine GeomechLoggingDestroy()

  implicit none
  
  ! all kinds of stuff needs to be added here.
  
  deallocate(geomech_logging)
  nullify(geomech_logging)
  
end subroutine GeomechLoggingDestroy

end module Geomechanics_Logging_module

#endif
