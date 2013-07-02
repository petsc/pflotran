#ifdef SURFACE_FLOW

module Timestepper_Surface_class

  use Timestepper_Base_class

  implicit none

#include "definitions.h"

  private

  type, public, extends(stepper_base_type) :: timestepper_surface_type
  contains
    procedure, public :: SetTargetTime => TimeStepperSurfaceSetTargetTime
  end type timestepper_surface_type

contains

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/02/13
! ************************************************************************** !
subroutine TimeStepperSurfaceSetTargetTime(timestepper,sync_time,option, &
                                        stop_flag,plot_flag, &
                                        transient_plot_flag)

  use Option_module

  implicit none

  class(timestepper_surface_type) :: timestepper
  PetscReal :: sync_time
  type(option_type) :: option
  PetscInt :: stop_flag
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag


  call printErrMsg(option,'debugging in TimeStepperSurfaceSetTargetTime')

end subroutine TimeStepperSurfaceSetTargetTime

end module Timestepper_Surface_class

#endif