module Checkpoint_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private

  public :: OpenCheckpointFile, &
            CloseCheckpointFile, &
            CheckpointFlowProcessModel, &
            RestartFlowProcessModel

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdef.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petsclog.h"
#include "finclude/petscviewer.h"
#include "finclude/petscbag.h"

contains

! ************************************************************************** !

subroutine OpenCheckpointFile(viewer,id,option,id_stamp)
  ! 
  ! Opens checkpoint file; sets format
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Option_module
  use String_module, only : StringNull
  
  implicit none

#include "finclude/petscviewer.h"

  PetscViewer :: viewer
  PetscInt :: id
  type(option_type) :: option
  character(len=MAXWORDLENGTH), optional, intent(in) :: id_stamp
  PetscErrorCode :: ierr
  
  character(len=MAXWORDLENGTH) :: id_string
  character(len=MAXSTRINGLENGTH) :: filename

  write(id_string,'(i8)') id
  if (present(id_stamp)) then
     if (.not. StringNull(id_stamp)) then
        id_string = id_stamp
     end if
  else if (id < 0) then
     id_string = 'restart'
  end if
  !else if (id >= 0) then --> use default id

  filename = trim(option%global_prefix) // &
       trim(option%group_prefix) // &
       '-' // trim(adjustl(id_string)) // '.chk'

  !geh: To skip .info file, need to split PetscViewerBinaryOpen() 
  !     into the routines it calls so that PetscViewerBinarySkipInfo()
  !     can be called after PetscViewerSetType(), but before
  !     PetscViewerFileSetName().  See note in PETSc docs.
  !call PetscViewerBinaryOpen(option%mycomm, filename, FILE_MODE_WRITE, &
  !                           viewer, ierr)
  call PetscViewerCreate(option%mycomm,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerSetType(viewer,PETSCVIEWERBINARY,ierr);CHKERRQ(ierr)
  call PetscViewerFileSetMode(viewer,FILE_MODE_WRITE,ierr);CHKERRQ(ierr)
  call PetscViewerBinarySkipInfo(viewer,ierr);CHKERRQ(ierr)
  call PetscViewerFileSetName(viewer,filename,ierr);CHKERRQ(ierr)
  
  write(option%io_buffer,'(" --> Dump checkpoint file: ", a64)') &
    trim(adjustl(filename))
  call printMsg(option)

end subroutine OpenCheckpointFile

! ************************************************************************** !

subroutine CloseCheckpointFile(viewer)
  ! 
  ! Closes checkpoint file
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Option_module
  
  implicit none

#include "finclude/petscviewer.h"

  PetscViewer :: viewer
  PetscErrorCode :: ierr

  call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)

end subroutine CloseCheckpointFile

! ************************************************************************** !

subroutine CheckpointFlowProcessModel(viewer,realization)
  ! 
  ! Checkpoints flow process model vectors
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Option_module
  use Realization_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Material_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_YZ, PERMEABILITY_XZ
  
  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  PetscViewer :: viewer
  type(realization_type) :: realization
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec
  
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid
  
  global_vec = 0
  
  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    ! grid%flow_xx is the vector into which all of the primary variables are 
    ! packed for the SNESSolve().
    call VecView(field%flow_xx, viewer, ierr);CHKERRQ(ierr)


    ! If we are running with multiple phases, we need to dump the vector 
    ! that indicates what phases are present, as well as the 'var' vector 
    ! that holds variables derived from the primary ones via the translator.
    select case(option%iflowmode)
      case(MPH_MODE,TH_MODE,RICHARDS_MODE,IMS_MODE,MIS_MODE, &
           FLASH2_MODE,G_MODE)
        call DiscretizationLocalToGlobal(realization%discretization, &
                                         field%iphas_loc,global_vec,ONEDOF)
        call VecView(global_vec, viewer, ierr);CHKERRQ(ierr)
       case default
    end select 

    ! Porosity and permeability.
    ! (We only write diagonal terms of the permeability tensor for now, 
    ! since we have yet to add the full-tensor formulation.)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_X,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Y,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Z,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)

    if (grid%itype == STRUCTURED_GRID_MIMETIC) then 
      if (option%iflowmode /= RICHARDS_MODE) then
        option%io_buffer = 'Checkpointing of mimetic not set up for outside Richards.'
        call printErrMsg(option)
      endif

      call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                   field%work_loc,PERMEABILITY_XZ,ZERO_INTEGER)
      call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                       global_vec,ONEDOF)
      call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)

      call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                   field%work_loc,PERMEABILITY_XY,ZERO_INTEGER)
      call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                       global_vec,ONEDOF)
      call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)

      call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                   field%work_loc,PERMEABILITY_YZ,ZERO_INTEGER)
      call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                       global_vec,ONEDOF)
      call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)

      call VecView(field%flow_xx_faces, viewer, ierr);CHKERRQ(ierr)
    end if

  endif
  
  if (global_vec /= 0) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif  
  
end subroutine CheckpointFlowProcessModel

! ************************************************************************** !

subroutine RestartFlowProcessModel(viewer,realization)
  ! 
  ! Restarts flow process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

  use Option_module
  use Realization_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Global_module
  use Material_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_YZ, PERMEABILITY_XZ, STATE
  
  implicit none

#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  PetscViewer :: viewer
  type(realization_type) :: realization
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec
  
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid
  
  global_vec = 0
  
  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
  ! Load the PETSc vectors.
    call VecLoad(field%flow_xx,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                     field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx,field%flow_yy,ierr);CHKERRQ(ierr)

    select case(option%iflowmode)
      case(MPH_MODE,TH_MODE,RICHARDS_MODE,IMS_MODE,MIS_MODE, &
           FLASH2_MODE,G_MODE)
        call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         field%iphas_loc,ONEDOF)
        call VecCopy(field%iphas_loc,field%iphas_old_loc,ierr);CHKERRQ(ierr)
        call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                        field%iphas_old_loc,ONEDOF)
        if (option%iflowmode == G_MODE) then
          ! need to copy iphase into global_auxvar%istate
          call GlobalSetAuxVarVecLoc(realization,field%iphas_loc,STATE, &
                                     ZERO_INTEGER)
        endif
        if (option%iflowmode == MPH_MODE) then
        ! set vardof vec in mphase
        endif
        if (option%iflowmode == IMS_MODE) then
        ! set vardof vec in mphase
        endif
        if (option%iflowmode == FLASH2_MODE) then
        ! set vardof vec in mphase
        endif
 
      case default
    end select
    
    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                      field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,POROSITY,ZERO_INTEGER)
    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                      field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_X,ZERO_INTEGER)
    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                      field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Y,ZERO_INTEGER)
    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                      field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Z,ZERO_INTEGER)
    
    if (grid%itype == STRUCTURED_GRID_MIMETIC) then
      if (option%iflowmode /= RICHARDS_MODE) then
        option%io_buffer = 'Restart of mimetic not set up for outside Richards.'
        call printErrMsg(option)
      endif

      call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
      call DiscretizationGlobalToLocal(discretization,global_vec, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                   field%work_loc,PERMEABILITY_XZ,ZERO_INTEGER)
      call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
      call DiscretizationGlobalToLocal(discretization,global_vec, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                   field%work_loc,PERMEABILITY_XY,ZERO_INTEGER)
      call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
      call DiscretizationGlobalToLocal(discretization,global_vec, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                   field%work_loc,PERMEABILITY_YZ,ZERO_INTEGER)

      call VecLoad(field%flow_xx_faces, viewer,ierr);CHKERRQ(ierr)
      call DiscretizationGlobalToLocalLP(discretization, field%flow_xx_faces, &
                                         field%flow_xx_loc_faces, NFLOWDOF)
      call VecCopy(field%flow_xx_faces,field%flow_yy_faces,ierr);CHKERRQ(ierr)
    end if
    
  endif
  
  if (global_vec /= 0) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif  
  
end subroutine RestartFlowProcessModel

end module Checkpoint_module
