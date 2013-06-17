#ifdef GEOMECH

module Geomechanics_Logging_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"

! stages
PetscInt, parameter, public :: INIT_STAGE = 1
PetscInt, parameter, public :: TS_STAGE = 2
PetscInt, parameter, public :: GEOMECH_STAGE = 3
PetscInt, parameter, public :: OUTPUT_STAGE = 4

  type, public :: geomech_logging_type 
  
    PetscLogStage :: stage(10)
    
    PetscClassId :: class_pflotran
 
    PetscLogEvent :: event_geomech_condition_read
    PetscLogEvent :: event_geomech_condition_read_values

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

  call PetscLogStageRegister('Init Stage',  & 
                             geomech_logging%stage(INIT_STAGE),ierr)
  call PetscLogStageRegister('Time Step Stage', &
                             geomech_logging%stage(TS_STAGE),ierr)
  call PetscLogStageRegister('Flow Stage', &
                             geomech_logging%stage(GEOMECH_STAGE),ierr)
  call PetscLogStageRegister('Output Stage', &
                             geomech_logging%stage(OUTPUT_STAGE),ierr)
                             
  call PetscClassIdRegister('PFLOTRAN',geomech_logging%class_pflotran,ierr)

  call PetscLogEventRegister('GeomechCondRead', &
                             geomech_logging%class_pflotran, &
                             geomech_logging%event_geomech_condition_read,ierr)
 
  call PetscLogEventRegister('FlowCondReadVals', &
                             geomech_logging%class_pflotran, &
                             geomech_logging%event_geomech_condition_read_values,ierr)
  
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
