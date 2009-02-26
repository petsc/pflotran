module Discretization_module

  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module
  use AMR_Grid_Module

  implicit none

  private
 
#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"

  type, public :: discretization_type
    PetscInt :: itype  ! type of discretization (e.g. structured, unstructured, etc.)
    character(len=MAXWORDLENGTH) :: ctype
    PetscReal :: origin(3) ! origin of global domain
    type(grid_type), pointer :: grid  ! pointer to a grid object
    type(amrgrid_type), pointer :: amrgrid  ! pointer to an amr grid object
    DM :: dm_1_dof, dm_nflowdof, dm_ntrandof
    DM, pointer :: dmc_nflowdof(:), dmc_ntrandof(:)
      ! Arrays containing hierarchy of coarsened DMs, for use with Galerkin 
      ! multigrid.  Element i of each array is a *finer* DM than element i-1.

    PetscInt :: dm_index_to_ndof(5) ! mapping between a dm_ptr to the number of degrees of freedom

  end type discretization_type

  public :: DiscretizationCreate, &
            DiscretizationDestroy, &
            DiscretizationRead, &
            DiscretizationCreateVector, &
            DiscretizationDuplicateVector, &         
            DiscretizationCreateJacobian, &
            DiscretizationCreateInterpolation, &
            DiscretizationCreateColoring, &
            DiscretizationGlobalToLocal, &
            DiscretizationLocalToGlobal, &
            DiscretizationLocalToLocal, &
            DiscretizationGlobalToNatural, &
            DiscretizationNaturalToGlobal, &
            DiscretizationGlobalToLocalBegin, &
            DiscretizationGlobalToLocalEnd, &
            DiscretizationLocalToLocalBegin, &
            DiscretizationLocalToLocalEnd, &
            DiscretizGlobalToNaturalBegin, &
            DiscretizGlobalToNaturalEnd, &
            DiscretizNaturalToGlobalBegin, &
            DiscretizNaturalToGlobalEnd, &
            DiscretizationCreateDMs
  
contains

! ************************************************************************** !
!
! DiscretizationCreate: Creates a structured or unstructured discretization
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function DiscretizationCreate()

  implicit none
  
  type(discretization_type), pointer :: DiscretizationCreate
  
  type(discretization_type), pointer :: discretization
  
  allocate(discretization)
  discretization%ctype = ''
  discretization%itype = 0
  discretization%origin = 0.d0

  ! nullify DM pointers
  discretization%dm_1_dof = 0
  discretization%dm_nflowdof = 0
  discretization%dm_ntrandof = 0
  nullify(discretization%dmc_nflowdof)
  nullify(discretization%dmc_ntrandof)
  
  nullify(discretization%grid)
  nullify(discretization%amrgrid)
  
  DiscretizationCreate => discretization

end function DiscretizationCreate

! ************************************************************************** !
!
! DiscretizationRead: Reads a discretization from the input file
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine DiscretizationRead(discretization,input,first_time,option)

  use Option_module
  use Input_module
  use String_module
  use AMR_Grid_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type),pointer :: discretization
  PetscTruth :: first_time
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid
  type(structured_grid_type), pointer :: str_grid
  type(amrgrid_type), pointer :: amrgrid
  character(len=MAXWORDLENGTH) :: structured_grid_ctype
  PetscInt :: structured_grid_itype
  PetscInt :: nx, ny, nz
  PetscInt :: i

  nx = 0
  ny = 0
  nz = 0

! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  do
  
    call InputReadFlotranString(input,option)
    if (input%ierr /= 0) exit

    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GRID')
    call StringToUpper(word)
      
    if (first_time) then ! first time read
      select case(trim(word))
        case('TYPE')
          call InputReadWord(input,option,discretization%ctype,PETSC_TRUE)
          call InputErrorMsg(input,option,'type','GRID')   
          call StringToLower(discretization%ctype)
          select case(trim(discretization%ctype))
            case('structured')
              discretization%itype = STRUCTURED_GRID
              call InputReadWord(input,option,structured_grid_ctype,PETSC_TRUE)
              call InputDefaultMsg(input,option,'structured_grid_type')   
              call StringToLower(structured_grid_ctype)
              select case(trim(structured_grid_ctype))
                case('cartesian')
                  structured_grid_itype = CARTESIAN_GRID
                case('cylindrical')
                  structured_grid_itype = CYLINDRICAL_GRID
                case('spherical')
                  structured_grid_itype = SPHERICAL_GRID
                case default
                  structured_grid_itype = CARTESIAN_GRID
                  structured_grid_ctype = 'cartesian'
              end select
            case('unstructured')
              discretization%itype = UNSTRUCTURED_GRID
            case('amr')
              discretization%itype = AMR_GRID
            case default
              option%io_buffer = 'Discretization type: ' // &
                                 trim(discretization%ctype) // &
                                 ' not recognized.'
              call printErrMsg(option)
          end select    
        case('NXYZ')
          call InputReadInt(input,option,nx)
          call InputErrorMsg(input,option,'nx','GRID')
          call InputReadInt(input,option,ny)
          call InputErrorMsg(input,option,'ny','GRID')
          call InputReadInt(input,option,nz)
          call InputErrorMsg(input,option,'nz','GRID')
          if (structured_grid_itype /= CARTESIAN_GRID) then
            ny = 1 ! cylindrical and spherical have 1 cell in Y
            if (structured_grid_itype /= CYLINDRICAL_GRID) nz = 1 ! spherical has 1 cell in Z
          endif
        case('ORIG','ORIGIN')
          call InputReadDouble(input,option,discretization%origin(X_DIRECTION))
          call InputErrorMsg(input,option,'X direction','Origin')
          call InputReadDouble(input,option,discretization%origin(Y_DIRECTION))
          call InputErrorMsg(input,option,'Y direction','Origin')
          call InputReadDouble(input,option,discretization%origin(Z_DIRECTION))
          call InputErrorMsg(input,option,'Z direction','Origin')        
        case('FILE')
        case ('GRAVITY')
        case ('INVERT_Z')
        case('DXYZ')
          call InputSkipToEND(input,option,word) 
        case('BOUNDS')
          call InputSkipToEND(input,option,word) 
        case default
          option%io_buffer = 'Keyword: ' // trim(word) // &
                   ' not recognized in DISCRETIZATION, first read.'
          call printErrMsg(option)          
      end select 
    else ! should be the second time it is read
      select case(trim(word))
        case('TYPE')
        case('NXYZ')
        case('ORIG','ORIGIN')
        case('FILE')
        case('DXYZ')
          select case(discretization%itype)
            case(STRUCTURED_GRID)
              call StructuredGridReadDXYZ(discretization%grid%structured_grid,input,option)
            case(AMR_GRID)
              call AMRGridReadDXYZ(discretization%amrgrid,input,option)
            case default
              call printErrMsg(option,'Keyword "DXYZ" not supported for unstructured grid')
          end select
          call InputReadFlotranString(input,option) ! read END card
          call InputReadStringErrorMsg(input,option,'DISCRETIZATION,DXYZ,END')
          if (.not.(InputCheckExit(input,option))) then
            option%io_buffer = 'Card DXYZ should include either 3 entires ' // &
                     '(one for each grid direction or NX+NY+NZ entries)'
            call printErrMsg(option)
          endif
        case('BOUNDS')
          select case(discretization%itype)
            case(STRUCTURED_GRID)
              grid => discretization%grid
              if (grid%structured_grid%itype == CARTESIAN_GRID .or. &
                  grid%structured_grid%itype == CYLINDRICAL_GRID .or. &
                  grid%structured_grid%itype == SPHERICAL_GRID) then
                call InputReadFlotranString(input,option) ! x-direction
                call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,X or R')
                call InputReadDouble(input,option,grid%structured_grid%bounds(X_DIRECTION,LOWER))
                call InputErrorMsg(input,option,'Lower X or R','BOUNDS')
                call InputReadDouble(input,option,grid%structured_grid%bounds(X_DIRECTION,UPPER))
                call InputErrorMsg(input,option,'Upper X or R','BOUNDS')
              endif
              if (grid%structured_grid%itype == CARTESIAN_GRID) then
                call InputReadFlotranString(input,option) ! y-direction
                call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,Y')
                call InputReadDouble(input,option,grid%structured_grid%bounds(Y_DIRECTION,LOWER))
                call InputErrorMsg(input,option,'Lower Y','BOUNDS')
                call InputReadDouble(input,option,grid%structured_grid%bounds(Y_DIRECTION,UPPER))
                call InputErrorMsg(input,option,'Upper Y','BOUNDS')
              endif
              if (grid%structured_grid%itype == CARTESIAN_GRID .or. &
                  grid%structured_grid%itype == CYLINDRICAL_GRID) then
                call InputReadFlotranString(input,option) ! z-direction
                call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,Z')
                call InputReadDouble(input,option,grid%structured_grid%bounds(Z_DIRECTION,LOWER))
                call InputErrorMsg(input,option,'Lower Z','BOUNDS')
                call InputReadDouble(input,option,grid%structured_grid%bounds(Z_DIRECTION,UPPER))
                call InputErrorMsg(input,option,'Upper Z','BOUNDS')
              endif
              call InputReadFlotranString(input,option) ! z-direction
              call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,Z')
              if (.not.(InputCheckExit(input,option))) then
                if (OptionPrintToScreen(option)) then
                  if (grid%structured_grid%itype == CARTESIAN_GRID) then
                    print *, 'BOUNDS card for a cartesian structured grid must include ' // &
                             '5 lines.  I.e.'
                    print *, 'BOUNDS'
                    print *, 'x_min, x_max'
                    print *, 'y_min, y_max'
                    print *, 'z_min, z_max'
                    print *, 'END'
                  else if (grid%structured_grid%itype == CYLINDRICAL_GRID) then
                    print *, 'BOUNDS card for a cylindrical structured grid must include ' // &
                             '4 lines.  I.e.'
                    print *, 'BOUNDS'
                    print *, 'r_min, r_max'
                    print *, 'z_min, z_max'
                    print *, 'END'
                  else if (grid%structured_grid%itype == SPHERICAL_GRID) then
                    print *, 'BOUNDS card for a spherical structured grid must include ' // &
                             '3 lines.  I.e.'
                    print *, 'BOUNDS'
                    print *, 'r_min, r_max'
                    print *, 'END'
                  endif
                endif
                stop
              endif            
              discretization%origin(X_DIRECTION) = grid%structured_grid%bounds(X_DIRECTION,LOWER)
              discretization%origin(Y_DIRECTION) = grid%structured_grid%bounds(Y_DIRECTION,LOWER)
              discretization%origin(Z_DIRECTION) = grid%structured_grid%bounds(Z_DIRECTION,LOWER)
           case(AMR_GRID)
              call InputSkipToEND(input,option,word) 
          end select
        case ('GRAVITY')
          select case(discretization%itype)
            case(STRUCTURED_GRID, UNSTRUCTURED_GRID)
               call InputReadDouble(input,option,discretization%grid%gravity(X_DIRECTION))
               call InputErrorMsg(input,option,'x-direction','GRAVITY')
               call InputReadDouble(input,option,discretization%grid%gravity(Y_DIRECTION))
               call InputErrorMsg(input,option,'y-direction','GRAVITY')
               call InputReadDouble(input,option,discretization%grid%gravity(Z_DIRECTION))
               call InputErrorMsg(input,option,'z-direction','GRAVITY')
               if (option%myrank == option%io_rank .and. &
                    option%print_to_screen) &
                    write(option%fid_out,'(/," *GRAV",/, &
                    & "  gravity    = "," [m/s^2]",3x,3pe12.4 &
                    & )') discretization%grid%gravity(1:3)
            case(AMR_GRID)
           end select
        case ('INVERT_Z')
          if (associated(grid%structured_grid)) then
            grid%structured_grid%invert_z_axis = PETSC_TRUE
          endif          
        case default
          option%io_buffer = 'Keyword: ' // trim(word) // &
                   ' not recognized in DISCRETIZATION, second read.'
          call printErrMsg(option)          
      end select 
    endif
  
  enddo  

  if (first_time) then ! first time read

    select case(discretization%itype)
      case(UNSTRUCTURED_GRID,STRUCTURED_GRID)
        grid => GridCreate()
        select case(discretization%itype)
          case(STRUCTURED_GRID)      
            if (nx*ny*nz == 0) &
              call printErrMsg(option,'NXYZ not set correctly for structured grid.')
            str_grid => StructuredGridCreate()
            str_grid%nx = nx
            str_grid%ny = ny
            str_grid%nz = nz
            str_grid%nxy = str_grid%nx*str_grid%ny
            str_grid%nmax = str_grid%nxy*str_grid%nz
            grid%structured_grid => str_grid
            grid%nmax = str_grid%nmax
            grid%structured_grid%itype = structured_grid_itype
            grid%structured_grid%ctype = structured_grid_ctype
        end select
        discretization%grid => grid
        grid%itype = discretization%itype
        grid%ctype = discretization%ctype
      case(AMR_GRID)
         amrgrid=>discretization%amrgrid
         call AMRGridInitialize(amrgrid)
    end select

  else
    select case(discretization%itype)
      case(STRUCTURED_GRID)
        if (discretization%grid%structured_grid%invert_z_axis) then
          discretization%grid%gravity(Z_DIRECTION) =  &
            -discretization%grid%gravity(Z_DIRECTION)
        endif
    end select
  endif

end subroutine DiscretizationRead

! ************************************************************************** !
!
! DiscretizationCreateDMs: creates distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
subroutine DiscretizationCreateDMs(discretization,option)
      
  use Option_module    
      
  implicit none
  
  type(discretization_type) :: discretization
  type(option_type) :: option
      
  PetscInt :: ndof
  PetscInt, parameter :: stencil_width = 1
  PetscErrorCode :: ierr
  PetscInt :: i

  discretization%dm_index_to_ndof(ONEDOF) = 1
  discretization%dm_index_to_ndof(NPHASEDOF) = option%nphase
  discretization%dm_index_to_ndof(NFLOWDOF) = option%nflowdof
  discretization%dm_index_to_ndof(NTRANDOF) = option%ntrandof

  !-----------------------------------------------------------------------
  ! Generate the DA objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call DiscretizationCreateDM(discretization,discretization%dm_1_dof,ndof, &
                              stencil_width,option)
  
  if (option%nflowdof > 0) then
    ndof = option%nflowdof
    call DiscretizationCreateDM(discretization,discretization%dm_nflowdof, &
                                ndof,stencil_width,option)
  endif
  
  if (option%ntrandof > 0) then
    ndof = option%ntrandof
    call DiscretizationCreateDM(discretization,discretization%dm_ntrandof, &
                                ndof,stencil_width,option)
  endif

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      ! this function must be called to set up str_grid%nxs, etc.
      call StructGridComputeLocalBounds(discretization%grid%structured_grid, &
                                        discretization%dm_1_dof)    
      discretization%grid%nlmax = discretization%grid%structured_grid%nlmax
      discretization%grid%ngmax = discretization%grid%structured_grid%ngmax
    case(UNSTRUCTURED_GRID)
    case(AMR_GRID)
      call AMRGridComputeLocalBounds(discretization%amrgrid)
  end select
  
end subroutine DiscretizationCreateDMs

! ************************************************************************** !
!
! DiscretizationCreateDM: creates a distributed, parallel mesh/grid
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
subroutine DiscretizationCreateDM(discretization,dm,ndof,stencil_width, &
                                  option)

  use Option_module
  
  implicit none
  
  type(discretization_type) :: discretization
  DM :: dm
  PetscInt :: ndof
  PetscInt :: stencil_width
  type(option_type) :: option
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call StructuredGridCreateDA(discretization%grid%structured_grid, &
                                  dm,ndof,stencil_width,option)
    case(UNSTRUCTURED_GRID)
  end select

end subroutine DiscretizationCreateDM

! ************************************************************************** !
!
! DiscretizationCreateVector: Creates a global PETSc vector
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DiscretizationCreateVector(discretization,dm_index,vector, &
                                      vector_type,option)
  use Option_module                                      

  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  Vec :: vector
  PetscInt :: vector_type
  type(option_type) :: option
  PetscInt :: ndof
  PetscErrorCode :: ierr
  
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      select case (vector_type)
        case(GLOBAL)
          call DACreateGlobalVector(dm_ptr,vector,ierr)
        case(LOCAL)
          call DACreateLocalVector(dm_ptr,vector,ierr)
        case(NATURAL)
          call DACreateNaturalVector(dm_ptr,vector,ierr)
      end select
    case(UNSTRUCTURED_GRID)
    case(AMR_GRID)
      select case(dm_index)
        case(ONEDOF)
           ndof = 1
        case(NFLOWDOF)
           ndof = option%nflowdof
        case(NTRANDOF)
           ndof = option%ntrandof
      end select
      call AMRGridCreateVector(discretization%amrgrid, ndof, vector, &
                               vector_type, PETSC_FALSE, option)
  end select
  call VecSet(vector,0.d0,ierr)
  
end subroutine DiscretizationCreateVector

! ************************************************************************** !
!
! DiscretizationDuplicateVector: Creates a global PETSc vector
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DiscretizationDuplicateVector(discretization,vector1,vector2)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: vector1
  Vec :: vector2
  
  PetscErrorCode :: ierr
  call VecDuplicate(vector1,vector2,ierr)
  call VecCopy(vector1,vector2,ierr)
  
end subroutine DiscretizationDuplicateVector

! ************************************************************************** !
!
! DiscretizationGetDMPtrFromIndex: Returns the integer pointer for the DM referenced
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
function DiscretizationGetDMPtrFromIndex(discretization,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  
  DM :: DiscretizationGetDMPtrFromIndex
  
  select case (dm_index)
    case(ONEDOF)
      DiscretizationGetDMPtrFromIndex = discretization%dm_1_dof
    case(NFLOWDOF)
      DiscretizationGetDMPtrFromIndex = discretization%dm_nflowdof
    case(NTRANDOF)
      DiscretizationGetDMPtrFromIndex = discretization%dm_ntrandof

  end select  
  
end function DiscretizationGetDMPtrFromIndex

function DiscretizationGetDMCPtrFromIndex(discretization,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  
  DM, pointer :: DiscretizationGetDMCPtrFromIndex(:)
  
  select case (dm_index)
    case(NFLOWDOF)
      DiscretizationGetDMCPtrFromIndex = discretization%dmc_nflowdof
    case(NTRANDOF)
      DiscretizationGetDMCPtrFromIndex = discretization%dmc_ntrandof
  end select  
  
end function DiscretizationGetDMCPtrFromIndex

! ************************************************************************** !
!
! DiscretizationCreateJacobian: Creates Jacobian matrix associated with discretization
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DiscretizationCreateJacobian(discretization,dm_index,mat_type,Jacobian,option)

  use Option_module
  
  implicit none
  
  interface

     subroutine SAMRCreateMatrix(p_application, ndof, stencilsize, flowortransport, p_matrix)
#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
       PetscFortranAddr :: p_application
       PetscInt :: ndof
       PetscInt :: stencilsize
       PetscInt :: flowortransport
       Mat :: p_matrix
     end subroutine SAMRCreateMatrix
     
     integer function hierarchy_number_levels(p_hierarchy)
       PetscFortranAddr, intent(inout) :: p_hierarchy
     end function hierarchy_number_levels
   
     integer function level_number_patches(p_hierarchy, ln)
       PetscFortranAddr, intent(inout) :: p_hierarchy
       integer, intent(in) :: ln
     end function level_number_patches

     logical function is_local_patch(p_hierarchy, ln, pn)
       PetscFortranAddr, intent(inout) :: p_hierarchy
       integer, intent(in) :: ln
       integer, intent(in) :: pn
     end function is_local_patch
     
  end interface

#include "finclude/petscis.h"
#include "finclude/petscis.h90"

  type(discretization_type) :: discretization
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  MatType :: mat_type
  Mat :: Jacobian
  type(option_type) :: option
  PetscInt :: ndof, stencilsize
  PetscInt, pointer :: indices(:)
  PetscInt :: ngmax
  PetscInt :: imax, nlevels, ln, npatches, pn, i
  PetscInt :: flowortransport
  DM :: dm_ptr
  ISLocalToGlobalMapping :: ptmap
  PetscTruth :: islocal

  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DAGetMatrix(dm_ptr,mat_type,Jacobian,ierr)
      call MatSetOption(Jacobian,MAT_KEEP_ZEROED_ROWS,PETSC_FALSE,ierr)
      call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
    case(UNSTRUCTURED_GRID)
    case(AMR_GRID)
       select case(dm_index)
       case(ONEDOF)
          ndof = 1
       case(NFLOWDOF)
          ndof = option%nflowdof
          flowortransport = 0
       case(NTRANDOF)
          ndof = option%ntrandof
          flowortransport = 1
       end select

       stencilsize=7;
       call MatCreateShell(option%mycomm, 0,0, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, Jacobian, ierr)
       call SAMRCreateMatrix(discretization%amrgrid%p_application, ndof, stencilsize, flowortransport, Jacobian)

       if(flowortransport==1) then
! we create a dummy mapping to satisfy the PETSc mapping requirements   
          imax=0
          nlevels =  hierarchy_number_levels(discretization%amrgrid%p_application)
          do ln=0,nlevels-1
             npatches = level_number_patches(discretization%amrgrid%p_application, ln )
             do pn=0,npatches-1
                islocal = is_local_patch(discretization%amrgrid%p_application, ln, pn);
                if(islocal) then
                   ngmax =  discretization%amrgrid%gridlevel(ln+1)%grids(pn+1)%grid_ptr%ngmax
                   imax = max(ngmax,imax) 
                endif
             enddo
          enddo
       
          allocate(indices(imax))
          do i=1,imax
             indices(i:i)=i-1
          enddo
          call ISLocalToGlobalMappingCreate(option%mycomm, imax, indices, ptmap, ierr)          
!          call ISSetIdentity(ptmap, ierr)
!          call ISLocalToGlobalMappingBlock(ptmap, ndof, bmap, ierr)
          call MatSetLocalToGlobalMappingBlock(Jacobian, ptmap, ierr)
       endif

  end select

end subroutine DiscretizationCreateJacobian

! ************************************************************************** !
!
! DiscretizationCreateInterpolation: Creates interpolation matrix associated 
! with discretization for geometric multigrid.
! author: Richard Mills
! date: 4/25/08.
!
! ************************************************************************** !
subroutine DiscretizationCreateInterpolation(discretization,dm_index, &
                                             interpolation,mg_levels_x, &
                                             mg_levels_y, mg_levels_z, &
                                             option)

  use Option_module
  
  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  Mat, pointer :: interpolation(:)
  PetscInt :: mg_levels_x, mg_levels_y, mg_levels_z
  type(option_type) :: option

  PetscInt :: mg_levels
  PetscInt :: refine_x, refine_y, refine_z
  DM, target :: dm_ptr
  DM, pointer :: dmc_ptr(:)
  PetscInt :: i
  DM, pointer :: dm_fine_ptr
    ! Used to point to finer-grid DM in the loop that constructst the 
    ! interpolation hierarchy.
  
  mg_levels = max(mg_levels_x, mg_levels_y, mg_levels_z)
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
!  dmc_ptr = DiscretizationGetDMCPtrFromIndex(discretization,dm_index)
  select case (dm_index)
    case(NFLOWDOF)
      allocate(discretization%dmc_nflowdof(mg_levels))
      dmc_ptr => discretization%dmc_nflowdof
    case(NTRANDOF)
      allocate(discretization%dmc_ntrandof(mg_levels))
      dmc_ptr => discretization%dmc_ntrandof
  end select  
   
  allocate(interpolation(mg_levels))

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      dm_fine_ptr => dm_ptr
      refine_x = 2; refine_y = 2; refine_z = 2
      do i=mg_levels-1,1,-1
        ! If number of coarsenings performed so far exceeds mg_levels_x-1, 
        ! set refine_x = 1; likewise for y and z.
        if (i <= mg_levels - mg_levels_x ) refine_x = 1
        if (i <= mg_levels - mg_levels_y ) refine_y = 1
        if (i <= mg_levels - mg_levels_z ) refine_z = 1
        call DASetRefinementFactor(dm_fine_ptr, refine_x, refine_y, refine_z, &
                                   ierr)
        call DASetInterpolationType(dm_fine_ptr, DA_Q0, ierr)
        call DACoarsen(dm_fine_ptr, option%mycomm, dmc_ptr(i), ierr)
        call DAGetInterpolation(dmc_ptr(i), dm_fine_ptr, interpolation(i), &
                                PETSC_NULL_OBJECT, ierr)
        dm_fine_ptr => dmc_ptr(i)
      enddo
    case(UNSTRUCTURED_GRID)
  end select

end subroutine DiscretizationCreateInterpolation

! ************************************************************************** !
!
! DiscretizationCreateColoring: Creates ISColoring for discretization
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DiscretizationCreateColoring(discretization,dm_index,option,coloring)

  use Option_module
  
  implicit none

#include "finclude/petscis.h"
#include "finclude/petscis.h90"
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(option_type) :: option
  ISColoring :: coloring

  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DAGetColoring(dm_ptr,IS_COLORING_GLOBAL,coloring,ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizationCreateColoring

! ************************************************************************** !
!
! DiscretizationGlobalToLocal: Performs global to local communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine DiscretizationGlobalToLocal(discretization,global_vec,local_vec,dm_index)

  implicit none
  
  interface
     subroutine SAMRGlobalToLocal(p_application, gvec, lvec, ierr)
       implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       PetscFortranAddr :: p_application
       Vec :: lvec
       Vec :: gvec
       PetscInt :: ndof
       PetscInt :: ierr
       
     end subroutine SAMRGlobalToLocal

  end interface

  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      call DMGlobalToLocalBegin(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
      call DMGlobalToLocalEnd(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
      case(AMR_GRID)
         call SAMRGlobalToLocal(discretization%amrgrid%p_application, global_vec, local_vec, ierr);
  end select
  
end subroutine DiscretizationGlobalToLocal
  
! ************************************************************************** !
!
! DiscretizationLocalToGlobal: Performs local to global communication with DM
! author: Glenn Hammond
! date: 1/02/08
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine DiscretizationLocalToGlobal(discretization,local_vec,global_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      call DMLocalToGlobal(dm_ptr,local_vec,INSERT_VALUES,global_vec,ierr)
      case(AMR_GRID)
         call VecCopy(local_vec, global_vec, ierr);
  end select
  
end subroutine DiscretizationLocalToGlobal
  
! ************************************************************************** !
!
! DiscretizationLocalToLocal: Performs local to local communication with DM
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine DiscretizationLocalToLocal(discretization,local_vec1,local_vec2,dm_index)

  implicit none
  
  interface
     subroutine SAMRLocalToLocal(p_application, gvec, lvec, ierr)
       implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       PetscFortranAddr :: p_application
       Vec :: lvec
       Vec :: gvec
       PetscInt :: ndof
       PetscInt :: ierr
       
     end subroutine SAMRLocalToLocal

  end interface
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DALocalToLocalBegin(dm_ptr,local_vec1,INSERT_VALUES,local_vec2,ierr)
      call DALocalToLocalEnd(dm_ptr,local_vec1,INSERT_VALUES,local_vec2,ierr)
    case(UNSTRUCTURED_GRID)
    case(AMR_GRID)
       call SAMRLocalToLocal(discretization%amrgrid%p_application, local_vec1, local_vec2, ierr);
  end select
  
end subroutine DiscretizationLocalToLocal
  
! ************************************************************************** !
!
! DiscretizationGlobalToNatural: Performs global to natural communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DAGlobalToNaturalBegin(dm_ptr,global_vec,INSERT_VALUES,natural_vec,ierr)
      call DAGlobalToNaturalEnd(dm_ptr,global_vec,INSERT_VALUES,natural_vec,ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizationGlobalToNatural

! ************************************************************************** !
!
! DiscretizationNaturalToGlobal: Performs natural to global communication with DM
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine DiscretizationNaturalToGlobal(discretization,natural_vec,global_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DANaturalToGlobalBegin(dm_ptr,natural_vec,INSERT_VALUES,global_vec,ierr)
      call DANaturalToGlobalEnd(dm_ptr,natural_vec,INSERT_VALUES,global_vec,ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizationNaturalToGlobal

! ************************************************************************** !
!
! DiscretizationGlobalToLocalBegin: Begins global to local communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine DiscretizationGlobalToLocalBegin(discretization,global_vec,local_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      call DMGlobalToLocalBegin(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
  end select
  
end subroutine DiscretizationGlobalToLocalBegin
  
! ************************************************************************** !
!
! DiscretizationGlobalToLocalEnd: Ends global to local communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine DiscretizationGlobalToLocalEnd(discretization,global_vec,local_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      call DMGlobalToLocalEnd(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
  end select
  
end subroutine DiscretizationGlobalToLocalEnd
  
! ************************************************************************** !
!
! DiscretizationLocalToLocalBegin: Begins local to local communication with DM
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine DiscretizationLocalToLocalBegin(discretization,local_vec1,local_vec2,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DALocalToLocalBegin(dm_ptr,local_vec1,INSERT_VALUES,local_vec2,ierr)
    case(UNSTRUCTURED_GRID)
  end select

end subroutine DiscretizationLocalToLocalBegin
  
! ************************************************************************** !
!
! DiscretizationLocalToLocalEnd: Ends local to local communication with DM
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine DiscretizationLocalToLocalEnd(discretization,local_vec1,local_vec2,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DALocalToLocalEnd(dm_ptr,local_vec1,INSERT_VALUES,local_vec2,ierr)
    case(UNSTRUCTURED_GRID)
  end select

end subroutine DiscretizationLocalToLocalEnd
  
! ************************************************************************** !
!
! DiscretizGlobalToNaturalBegin: Begins global to natural communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DiscretizGlobalToNaturalBegin(discretization,global_vec,natural_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DAGlobalToNaturalBegin(dm_ptr,global_vec,INSERT_VALUES,natural_vec,ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizGlobalToNaturalBegin

! ************************************************************************** !
!
! DiscretizGlobalToNaturalEnd: Ends global to natural communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DiscretizGlobalToNaturalEnd(discretization,global_vec,natural_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DAGlobalToNaturalEnd(dm_ptr,global_vec,INSERT_VALUES,natural_vec,ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizGlobalToNaturalEnd

! ************************************************************************** !
!
! DiscretizNaturalToGlobalBegin: Begins natural to global communication with DM
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine DiscretizNaturalToGlobalBegin(discretization,natural_vec,global_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DANaturalToGlobalBegin(dm_ptr,natural_vec,INSERT_VALUES,global_vec,ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizNaturalToGlobalBegin

! ************************************************************************** !
!
! DiscretizNaturalToGlobalEnd: Ends natural to global communication with DM
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine DiscretizNaturalToGlobalEnd(discretization,natural_vec,global_vec,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DANaturalToGlobalEnd(dm_ptr,natural_vec,INSERT_VALUES,global_vec,ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizNaturalToGlobalEnd

! ************************************************************************** !
!
! DiscretizationDestroy: Deallocates a discretization
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine DiscretizationDestroy(discretization)

  implicit none
  
  type(discretization_type), pointer :: discretization
  
  PetscErrorCode :: ierr
  integer :: i
    
  if (.not.associated(discretization)) return
      
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      if (discretization%dm_1_dof /= 0) &
        call DMDestroy(discretization%dm_1_dof,ierr)
      discretization%dm_1_dof = 0
      if (discretization%dm_nflowdof /= 0) &
        call DMDestroy(discretization%dm_nflowdof,ierr)
      discretization%dm_nflowdof = 0
      if (discretization%dm_ntrandof /= 0) &
        call DMDestroy(discretization%dm_ntrandof,ierr)
      discretization%dm_ntrandof = 0
      if (associated(discretization%dmc_nflowdof)) then
        do i=1,size(discretization%dmc_nflowdof)
          call DMDestroy(discretization%dmc_nflowdof(i),ierr)
        enddo
        deallocate(discretization%dmc_nflowdof)
      endif
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizationDestroy
 
end module Discretization_module
