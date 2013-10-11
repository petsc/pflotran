module Discretization_module

  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Explicit_module
  use MFD_Aux_module
  use MFD_module
  use DM_Kludge_module

  use PFLOTRAN_Constants_module

  implicit none

  private
 
#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmshell.h90"

  type, public :: discretization_type
    PetscInt :: itype  ! type of discretization (e.g. structured, unstructured, etc.)
    !geh: note that differentiating between implicit and explicit unstructured 
    !     grids is handled within the grid%itype variable, not discritization%itype
    character(len=MAXWORDLENGTH) :: ctype
    PetscReal :: origin(3) ! origin of global domain
    type(grid_type), pointer :: grid  ! pointer to a grid object
    character(len=MAXSTRINGLENGTH) :: filename

    type(dm_ptr_type), pointer :: dmc_nflowdof(:), dmc_ntrandof(:)
      ! Arrays containing hierarchy of coarsened DMs, for use with Galerkin 
      ! multigrid.  Element i of each array is a *finer* DM than element i-1.
    PetscInt :: dm_index_to_ndof(5) ! mapping between a dm_ptr to the number of degrees of freedom
    type(dm_ptr_type), pointer :: dm_1dof
    type(dm_ptr_type), pointer :: dm_nflowdof
    type(dm_ptr_type), pointer :: dm_ntrandof
    type(mfd_type), pointer :: MFD
    VecScatter :: tvd_ghost_scatter
    
    PetscInt :: stencil_width
    PetscInt :: stencil_type
    
    PetscInt :: hydr_flux_method
    PetscInt :: therm_flux_method
    PetscInt :: mech_flux_method
    PetscBool:: lsm_flux_method
    
  end type discretization_type

  public :: DiscretizationCreate, &
            DiscretizationDestroy, &
            DiscretizationReadRequiredCards, &
            DiscretizationRead, &
            DiscretizationCreateVector, &
            DiscretizationDuplicateVector, &         
            DiscretizationCreateJacobian, &
            DiscretizationCreateInterpolation, &
            DiscretizationCreateColoring, &
            DiscretizationGlobalToLocal, &
            DiscretizationGlobalToLocalFaces, &
            DiscretizationLocalToGlobal, &
            DiscretizationLocalToLocal, &
            DiscretizationLocalToLocalFaces, &
            DiscretizationGlobalToNatural, &
            DiscretizationNaturalToGlobal, &
            DiscretizationGlobalToLocalBegin, &
            DiscretizationGlobalToLocalEnd, &
            DiscretizGlobalToLocalFacesBegin, &
            DiscretizGlobalToLocalFacesEnd, &
            DiscretizLocalToLocalFacesBegin, &
            DiscretizLocalToLocalFacesEnd, &
            DiscretizationGlobalToLocalLP, &
            DiscretizationLocalToLocalLP, &
            DiscretizGlobalToLocalLPBegin, &
            DiscretizGlobalToLocalLPEnd, &
            DiscretizLocalToLocalLPBegin, &
            DiscretizLocalToLocalLPEnd, &
            DiscretizationLocalToLocalBegin, &
            DiscretizationLocalToLocalEnd, &
            DiscretizGlobalToNaturalBegin, &
            DiscretizGlobalToNaturalEnd, &
            DiscretizNaturalToGlobalBegin, &
            DiscretizNaturalToGlobalEnd, &
            DiscretizationCreateDMs,&
            DiscretizationGetDMPtrFromIndex, &
            DiscretizationUpdateTVDGhosts, &
            DiscretAOApplicationToPetsc
  
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
  discretization%filename = ''

  ! nullify DM pointers
  nullify(discretization%dmc_nflowdof)
  nullify(discretization%dmc_ntrandof)
  allocate(discretization%dm_1dof)
  allocate(discretization%dm_nflowdof)
  allocate(discretization%dm_ntrandof)
  discretization%dm_1dof%dm = 0
  discretization%dm_nflowdof%dm = 0
  discretization%dm_ntrandof%dm = 0
  nullify(discretization%dm_1dof%ugdm)
  nullify(discretization%dm_nflowdof%ugdm)
  nullify(discretization%dm_ntrandof%ugdm)
  
  nullify(discretization%grid)
  nullify(discretization%MFD)
  
  discretization%stencil_width = 1
  discretization%stencil_type = DMDA_STENCIL_STAR
  discretization%hydr_flux_method = TWO_POINT_FLUX
  discretization%therm_flux_method = TWO_POINT_FLUX
  discretization%mech_flux_method = TWO_POINT_FLUX
  discretization%lsm_flux_method = PETSC_FALSE

  discretization%tvd_ghost_scatter = 0
  
  DiscretizationCreate => discretization

end function DiscretizationCreate

! ************************************************************************** !
!
! DiscretizationReadRequiredCards: Reads a discretization from the input file
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine DiscretizationReadRequiredCards(discretization,input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type),pointer :: discretization
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid, grid2
  type(structured_grid_type), pointer :: str_grid
  type(unstructured_grid_type), pointer :: un_str_grid
  character(len=MAXWORDLENGTH) :: structured_grid_ctype
  character(len=MAXWORDLENGTH) :: unstructured_grid_ctype

  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: structured_grid_itype
  PetscInt :: unstructured_grid_itype
  PetscInt :: nx, ny, nz
  PetscInt :: i
  PetscReal :: tempreal

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
          case('unstructured','unstructured_explicit')
            discretization%itype = UNSTRUCTURED_GRID
            word = discretization%ctype
            discretization%ctype = 'unstructured'
            select case(word)
              case('unstructured')
                unstructured_grid_itype = IMPLICIT_UNSTRUCTURED_GRID
                unstructured_grid_ctype = 'implicit unstructured'
              case('unstructured_explicit')
                unstructured_grid_itype = EXPLICIT_UNSTRUCTURED_GRID
                unstructured_grid_ctype = 'explicit unstructured'
            end select
            call InputReadNChars(input,option,discretization%filename,MAXSTRINGLENGTH, &
                                 PETSC_TRUE)
            call InputErrorMsg(input,option,'unstructured filename','GRID')
          case('unstructured_mimetic')
            discretization%itype = UNSTRUCTURED_GRID_MIMETIC
            option%mimetic = PETSC_TRUE
            word = discretization%ctype
            discretization%ctype = 'unstructured_mimetic'
            unstructured_grid_itype = IMPLICIT_UNSTRUCTURED_GRID
            unstructured_grid_ctype = 'implicit unstructured'
            call InputReadNChars(input,option,discretization%filename,MAXSTRINGLENGTH, &
                                 PETSC_TRUE)
            call InputErrorMsg(input,option,'unstructured filename','GRID')
          case('structured_mimetic')
            discretization%itype = STRUCTURED_GRID_MIMETIC
            option%mimetic = PETSC_TRUE
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
      case('ORIGIN')
        call InputReadDouble(input,option,discretization%origin(X_DIRECTION))
        call InputErrorMsg(input,option,'X direction','Origin')
        call InputReadDouble(input,option,discretization%origin(Y_DIRECTION))
        call InputErrorMsg(input,option,'Y direction','Origin')
        call InputReadDouble(input,option,discretization%origin(Z_DIRECTION))
        call InputErrorMsg(input,option,'Z direction','Origin')        
      case('FILE','GRAVITY','INVERT_Z','MAX_CELLS_SHARING_A_VERTEX',&
           'STENCIL_WIDTH','STENCIL_TYPE','FLUX_METHOD')
      case('DXYZ','BOUNDS')
        call InputSkipToEND(input,option,word) 
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                 ' not recognized in DISCRETIZATION, first read.'
        call printErrMsg(option)          
    end select 
  enddo  

  if (discretization%itype == NULL_GRID) then
    option%io_buffer = 'Discretization type not defined under ' // &
                       'keyword GRID.' 
    call printErrMsg(option)
  endif
  
  grid => GridCreate()
  select case(discretization%itype)
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      un_str_grid => UGridCreate()
      select case(unstructured_grid_itype)
        case(IMPLICIT_UNSTRUCTURED_GRID)
          if (index(discretization%filename,'.h5') > 0) then
#if !defined(PETSC_HAVE_HDF5)
            option%io_buffer = 'PFLOTRAN must be built with HDF5 ' // &
              'support to read unstructured grid .h5 files'
            call printErrMsg(option)
#else

#ifdef SCORPIO
            call UGridReadHDF5PIOLib(un_str_grid,discretization%filename,option)
#else
            call UGridReadHDF5(un_str_grid,discretization%filename,option)
#endif
! #ifdef SCORPIO

#endif
!#if !defined(PETSC_HAVE_HDF5)

          else
            call UGridRead(un_str_grid,discretization%filename,option)
          endif
          grid%unstructured_grid => un_str_grid
        case(EXPLICIT_UNSTRUCTURED_GRID)
          un_str_grid%explicit_grid => UGridExplicitCreate()
          call ExplicitUGridRead(un_str_grid, &
                                 discretization%filename,option)
          grid%unstructured_grid => un_str_grid
      end select
      grid%itype = unstructured_grid_itype
      grid%ctype = unstructured_grid_ctype
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)      
      if (nx*ny*nz <= 0) &
        call printErrMsg(option,'NXYZ not set correctly for structured grid.')
      str_grid => StructGridCreate()
      str_grid%nx = nx
      str_grid%ny = ny
      str_grid%nz = nz
      str_grid%nxy = str_grid%nx*str_grid%ny
      str_grid%nmax = str_grid%nxy*str_grid%nz
      grid%structured_grid => str_grid
      grid%nmax = str_grid%nmax
      grid%structured_grid%itype = structured_grid_itype
      grid%structured_grid%ctype = structured_grid_ctype
      grid%itype = discretization%itype
      grid%ctype = discretization%ctype
  end select
  grid%discretization_itype=discretization%itype
  discretization%grid => grid
  nullify(grid)

end subroutine DiscretizationReadRequiredCards

! ************************************************************************** !
!
! DiscretizationRead: Reads a discretization from the input file
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine DiscretizationRead(discretization,input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type),pointer :: discretization
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid, grid2
  type(structured_grid_type), pointer :: str_grid
  type(unstructured_grid_type), pointer :: un_str_grid
  character(len=MAXWORDLENGTH) :: structured_grid_ctype
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: structured_grid_itype
  PetscInt :: nx, ny, nz
  PetscInt :: i
  PetscReal :: tempreal

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
      
    select case(trim(word))
      case('TYPE','NXYZ','ORIGIN','FILE')
      case('DXYZ')
        select case(discretization%itype)
          case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
            call StructGridReadDXYZ(discretization%grid%structured_grid,input,option)
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
          case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
            grid => discretization%grid

            ! read first line and we will split off the legacy approach vs. new
            call InputReadFlotranString(input,option)
            call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,X or Min Coordinate')
            string = input%buf

            do i = 1, 3
              call InputReadDouble(input,option,tempreal)
              if (input%ierr /= 0) exit
            enddo

            input%ierr = 0
            input%buf = string

            if (i == 3) then ! only 2 successfully read
              if (grid%structured_grid%itype == CARTESIAN_GRID .or. &
                  grid%structured_grid%itype == CYLINDRICAL_GRID .or. &
                  grid%structured_grid%itype == SPHERICAL_GRID) then
!geh                  call InputReadFlotranString(input,option) ! x-direction
!geh                  call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,X or R')
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
              else
                grid%structured_grid%bounds(Y_DIRECTION,LOWER) = 0.d0
                grid%structured_grid%bounds(Y_DIRECTION,UPPER) = 1.d0
              endif
              if (grid%structured_grid%itype == CARTESIAN_GRID .or. &
                  grid%structured_grid%itype == CYLINDRICAL_GRID) then
                call InputReadFlotranString(input,option) ! z-direction
                call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,Z')
                call InputReadDouble(input,option,grid%structured_grid%bounds(Z_DIRECTION,LOWER))
                call InputErrorMsg(input,option,'Lower Z','BOUNDS')
                call InputReadDouble(input,option,grid%structured_grid%bounds(Z_DIRECTION,UPPER))
                call InputErrorMsg(input,option,'Upper Z','BOUNDS')
              else
                grid%structured_grid%bounds(Z_DIRECTION,LOWER) = 0.d0
                grid%structured_grid%bounds(Z_DIRECTION,UPPER) = 1.d0
              endif
            else ! new min max coordinate approach
              select case(grid%structured_grid%itype)
                case(CARTESIAN_GRID)
                  i = 3
                case(CYLINDRICAL_GRID)
                  i = 2
                case(SPHERICAL_GRID)
                  i = 1
              end select
              call InputReadNDoubles(input,option, &
                                     grid%structured_grid%bounds(:,LOWER), &
                                     i)
              call InputErrorMsg(input,option,'Minimum Coordinate','BOUNDS')
              call InputReadFlotranString(input,option)
              call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,MAX COORDINATE')
              call InputReadNDoubles(input,option, &
                                     grid%structured_grid%bounds(:,UPPER), &
                                     i)
              call InputErrorMsg(input,option,'Maximum Coordinate','BOUNDS')
              if (grid%structured_grid%itype == CYLINDRICAL_GRID) then
                grid%structured_grid%bounds(Y_DIRECTION,LOWER) = 0.d0
                grid%structured_grid%bounds(Y_DIRECTION,UPPER) = 1.d0
              endif
              if (grid%structured_grid%itype == SPHERICAL_GRID) then
                grid%structured_grid%bounds(Z_DIRECTION,LOWER) = 0.d0
                grid%structured_grid%bounds(Z_DIRECTION,UPPER) = 1.d0
              endif
            endif
            call InputReadFlotranString(input,option)
            call InputReadStringErrorMsg(input,option,'DISCRETIZATION,BOUNDS,END')
            if (.not.(InputCheckExit(input,option))) then
              if (OptionPrintToScreen(option)) then
                if (grid%structured_grid%itype == CARTESIAN_GRID) then
                  print *, 'BOUNDS card for a cartesian structured grid must include ' // &
                           '4 lines.  I.e.'
                  print *, 'BOUNDS'
                  print *, '  x_min  y_min  z_min'
                  print *, '  x_max  y_max  z_max'
                  print *, 'END'
                else if (grid%structured_grid%itype == CYLINDRICAL_GRID) then
                  print *, 'BOUNDS card for a cylindrical structured grid must include ' // &
                           '4 lines.  I.e.'
                  print *, 'BOUNDS'
                  print *, '  r_min  z_min'
                  print *, '  r_max  z_max'
                  print *, 'END'
                else if (grid%structured_grid%itype == SPHERICAL_GRID) then
                  print *, 'BOUNDS card for a spherical structured grid must include ' // &
                           '4 lines.  I.e.'
                  print *, 'BOUNDS'
                  print *, '  r_min'
                  print *, '  r_max'
                  print *, 'END'
                endif
              endif
              stop
            endif            
            discretization%origin(X_DIRECTION) = grid%structured_grid%bounds(X_DIRECTION,LOWER)
            discretization%origin(Y_DIRECTION) = grid%structured_grid%bounds(Y_DIRECTION,LOWER)
            discretization%origin(Z_DIRECTION) = grid%structured_grid%bounds(Z_DIRECTION,LOWER)
        end select
      case ('GRAVITY')
        call InputReadDouble(input,option,option%gravity(X_DIRECTION))
        call InputErrorMsg(input,option,'x-direction','GRAVITY')
        call InputReadDouble(input,option,option%gravity(Y_DIRECTION))
        call InputErrorMsg(input,option,'y-direction','GRAVITY')
        call InputReadDouble(input,option,option%gravity(Z_DIRECTION))
        call InputErrorMsg(input,option,'z-direction','GRAVITY')
        if (option%myrank == option%io_rank .and. &
            option%print_to_screen) &
          write(option%fid_out,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1p3e12.4 &
            & )') option%gravity(1:3)
      case ('MAX_CELLS_SHARING_A_VERTEX')
        if (associated(discretization%grid%unstructured_grid)) then
          call InputReadInt(input,option,discretization%grid% &
                            unstructured_grid%max_cells_sharing_a_vertex)
          call InputErrorMsg(input,option,'max_cells_sharing_a_vertex', &
                             'GRID')
        endif          
      case ('INVERT_Z')
      case ('STENCIL_WIDTH')
        call InputReadInt(input,option,discretization%stencil_width)
        call InputErrorMsg(input,option,'stencil_width', &
                           'GRID')
      case ('STENCIL_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','GRID')
        call StringToUpper(word)
        select case(trim(word))
          case ('BOX')
            discretization%stencil_type = DMDA_STENCIL_BOX
          case ('STAR')
            discretization%stencil_type = DMDA_STENCIL_STAR
          case default
            option%io_buffer = 'Keyword: ' // trim(word) // &
                 ' not recognized in DISCRETIZATION, second read.'
            call printErrMsg(option)
        end select
      case ('FLUX_METHOD')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','GRID')
        call StringToUpper(word)
        select case(trim(word))
          case ('HYDRAULICS')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','GRID')
            call StringToUpper(word)
            select case(trim(word))
              case ('TWO_POINT_FLUX')
                discretization%hydr_flux_method = TWO_POINT_FLUX
              case ('LSM_FLUX')
                discretization%hydr_flux_method = LSM_FLUX
                discretization%lsm_flux_method = PETSC_TRUE
              case default
                option%io_buffer = 'Keyword: ' // trim(word) // &
                  ' not recognized in DISCRETIZATION, second read.'
                call printErrMsg(option)
            end select
          case ('THERMAL')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','GRID')
            call StringToUpper(word)
            select case(trim(word))
              case ('TWO_POINT_FLUX')
                discretization%therm_flux_method = TWO_POINT_FLUX
              case ('LSM_FLUX')
                discretization%therm_flux_method = LSM_FLUX
                discretization%lsm_flux_method = PETSC_TRUE
              case default
                option%io_buffer = 'Keyword: ' // trim(word) // &
                  ' not recognized in DISCRETIZATION, second read.'
                call printErrMsg(option)
            end select
          case default
            option%io_buffer = 'Keyword: ' // trim(word) // &
            ' not recognized in DISCRETIZATION, second read.'
            call printErrMsg(option)
        end select
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                 ' not recognized in DISCRETIZATION, second read.'
        call printErrMsg(option)          
    end select 
  enddo  

  select case(discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      if (discretization%grid%structured_grid%invert_z_axis) then
        option%gravity(Z_DIRECTION) = -1.d0*option%gravity(Z_DIRECTION)
      endif
  end select
  
  if(discretization%lsm_flux_method) then
    if(discretization%stencil_type/=DMDA_STENCIL_BOX) then
      option%io_buffer='For LSM-Flux option, the stencil type needs to be DMDA_STENCIL_BOX'
      call printMsg(option)
      discretization%stencil_type = DMDA_STENCIL_BOX
      discretization%STENCIL_WIDTH = 2
    endif
  endif

end subroutine DiscretizationRead

! ************************************************************************** !
!
! DiscretizationCreateDMs: creates distributed, parallel meshes/grids
! If there are multiple degrees of freedom per grid cell, this will call 
! DiscretizationCreateDM() multiple times to create the DMs corresponding 
! to one degree of freedom grid cell and those corresponding to multiple 
! degrees of freedom per cell.
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
  !PetscInt, parameter :: stencil_width = 1
  PetscErrorCode :: ierr
  PetscInt :: i
  type(unstructured_grid_type), pointer :: ugrid

  select case(discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      discretization%dm_index_to_ndof(ONEDOF) = 1
      discretization%dm_index_to_ndof(NPHASEDOF) = option%nphase
      discretization%dm_index_to_ndof(NFLOWDOF) = option%nflowdof
      discretization%dm_index_to_ndof(NTRANDOF) = option%ntrandof
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)

      ! petsc will call parmetis to calculate the graph/dual
#if !defined(PETSC_HAVE_PARMETIS)
      option%io_buffer = &
        'Must compile with Parmetis in order to use unstructured grids.'
      call printErrMsg(option)
#endif
    
      select case(discretization%grid%itype)
        case(IMPLICIT_UNSTRUCTURED_GRID)
          call UGridDecompose(discretization%grid%unstructured_grid, &
                              option)
          if(discretization%lsm_flux_method) then
            call UGridGrowStencilSupport(discretization%grid%unstructured_grid, &
                                         discretization%stencil_width, &
                                         discretization%grid%ghosted_level, &
                                         option)
          endif
        case(EXPLICIT_UNSTRUCTURED_GRID)
          ugrid => discretization%grid%unstructured_grid
          call ExplicitUGridDecompose(ugrid,option)
      end select
  end select


  !-----------------------------------------------------------------------
  ! Generate the DM objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call DiscretizationCreateDM(discretization,discretization%dm_1dof, &
                              ndof,discretization%stencil_width, &
                              discretization%stencil_type,option)
  
  if (option%nflowdof > 0) then
    ndof = option%nflowdof
    call DiscretizationCreateDM(discretization,discretization%dm_nflowdof, &
                                ndof,discretization%stencil_width, &
                                discretization%stencil_type,option)
  endif
  
  if (option%ntrandof > 0) then
    ndof = option%ntrandof
    call DiscretizationCreateDM(discretization,discretization%dm_ntrandof, &
                                ndof,discretization%stencil_width, &
                                discretization%stencil_type,option)
  endif


  select case(discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      ! this function must be called to set up str_grid%lxs, etc.
      call StructGridComputeLocalBounds(discretization%grid%structured_grid, &
                                        discretization%dm_1dof%dm,option)    
      discretization%grid%nlmax = discretization%grid%structured_grid%nlmax
      discretization%grid%ngmax = discretization%grid%structured_grid%ngmax
      discretization%grid%global_offset = &
        discretization%grid%structured_grid%global_offset
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      discretization%grid%nmax = discretization%grid%unstructured_grid%nmax
      discretization%grid%nlmax = discretization%grid%unstructured_grid%nlmax
      discretization%grid%ngmax = discretization%grid%unstructured_grid%ngmax
      discretization%grid%global_offset = &
        discretization%grid%unstructured_grid%global_offset
  end select

end subroutine DiscretizationCreateDMs

! ************************************************************************** !
!
! DiscretizationCreateDM: creates a distributed, parallel mesh/grid
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
subroutine DiscretizationCreateDM(discretization,dm_ptr,ndof,stencil_width, &
                                  stencil_type,option)

  use Option_module
  
  implicit none
  
  type(discretization_type) :: discretization
  type(dm_ptr_type), pointer :: dm_ptr
  PetscInt :: ndof
  PetscInt :: stencil_width,stencil_type
  type(option_type) :: option
  PetscErrorCode :: ierr

  select case(discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      call StructGridCreateDM(discretization%grid%structured_grid, &
                                  dm_ptr%dm,ndof,stencil_width,stencil_type,option)
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      call UGridCreateUGDM(discretization%grid%unstructured_grid, &
                           dm_ptr%ugdm,ndof,option)
      call DMShellCreate(option%mycomm,dm_ptr%dm,ierr)
      call DMShellSetGlobalToLocalVecScatter(dm_ptr%dm,dm_ptr%ugdm%scatter_gtol,ierr)
      call DMShellSetLocalToGlobalVecScatter(dm_ptr%dm,dm_ptr%ugdm%scatter_ltog,ierr)
      call DMShellSetLocalToLocalVecScatter(dm_ptr%dm,dm_ptr%ugdm%scatter_ltol,ierr)
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
  
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)

  select case(discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      select case (vector_type)
        case(GLOBAL)
          call DMCreateGlobalVector(dm_ptr%dm,vector,ierr)
        case(LOCAL)
          call DMCreateLocalVector(dm_ptr%dm,vector,ierr)
        case(NATURAL)
          call DMDACreateNaturalVector(dm_ptr%dm,vector,ierr)
      end select
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      call UGridDMCreateVector(discretization%grid%unstructured_grid, &
                               dm_ptr%ugdm,vector, &
                               vector_type,option)
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
  
  type(dm_ptr_type), pointer :: DiscretizationGetDMPtrFromIndex
  
  select case (dm_index)
    case(ONEDOF)
      DiscretizationGetDMPtrFromIndex => discretization%dm_1dof
    case(NFLOWDOF)
      DiscretizationGetDMPtrFromIndex => discretization%dm_nflowdof
    case(NTRANDOF)
      DiscretizationGetDMPtrFromIndex => discretization%dm_ntrandof
  end select  
  
end function DiscretizationGetDMPtrFromIndex

! ************************************************************************** !
! ************************************************************************** !
function DiscretizationGetDMCPtrFromIndex(discretization,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  
  type(dm_ptr_type), pointer :: DiscretizationGetDMCPtrFromIndex(:)
  
  select case (dm_index)
    case(NFLOWDOF)
      DiscretizationGetDMCPtrFromIndex => discretization%dmc_nflowdof
    case(NTRANDOF)
      DiscretizationGetDMCPtrFromIndex => discretization%dmc_ntrandof
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
  type(dm_ptr_type), pointer :: dm_ptr
  ISLocalToGlobalMapping :: ptmap
  PetscInt :: islocal

  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)


  select case(discretization%itype)
    case(STRUCTURED_GRID)
#ifndef DMGET
      call DMSetMatType(dm_ptr%dm,mat_type,ierr)
      call DMCreateMatrix(dm_ptr%dm,Jacobian,ierr)
#else
      call DMGetMatrix(dm_ptr%dm,mat_type,Jacobian,ierr)
#endif
      call MatSetOption(Jacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,ierr)
      call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
    case(UNSTRUCTURED_GRID)
      call UGridDMCreateJacobian(discretization%grid%unstructured_grid, &
                                 dm_ptr%ugdm,mat_type,Jacobian,option)
      call MatSetOption(Jacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,ierr)
      call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
    case(UNSTRUCTURED_GRID_MIMETIC)
      select case(dm_index)
        case(NFLOWDOF)
          call MFDCreateJacobianLP(discretization%grid, discretization%MFD, mat_type, Jacobian, option)
          call MatSetOption(Jacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,ierr)
          call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
        case(NTRANDOF)
          call UGridDMCreateJacobian(discretization%grid%unstructured_grid, &
                                     dm_ptr%ugdm,mat_type,Jacobian,option)
          call MatSetOption(Jacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,ierr)
          call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
      end select
    case(STRUCTURED_GRID_MIMETIC)
#ifdef DASVYAT
      select case(dm_index)
        case(NFLOWDOF)
!          call MFDCreateJacobian(discretization%grid, discretization%MFD, mat_type, Jacobian, option)
          call MFDCreateJacobianLP(discretization%grid, discretization%MFD, mat_type, Jacobian, option)
          call MatSetOption(Jacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,ierr)
          call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
        case(NTRANDOF)
#ifndef DMGET
          call DMSetMatType(dm_ptr%dm,mat_type,ierr)
          call DMCreateMatrix(dm_ptr%dm,Jacobian,ierr)
#else
          call DMGetMatrix(dm_ptr%dm,mat_type,Jacobian,ierr)
#endif
          call MatSetOption(Jacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,ierr)
          call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
      end select
#endif
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
  type(dm_ptr_type), pointer :: dm_ptr
  type(dm_ptr_type), pointer :: dmc_ptr(:)
  PetscInt :: i
  type(dm_ptr_type), pointer :: dm_fine_ptr
    ! Used to point to finer-grid DM in the loop that constructst the 
    ! interpolation hierarchy.
  
  mg_levels = max(mg_levels_x, mg_levels_y, mg_levels_z)
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
!  dmc_ptr = DiscretizationGetDMCPtrFromIndex(discretization,dm_index)
  select case (dm_index)
    case(NFLOWDOF)
      allocate(discretization%dmc_nflowdof(mg_levels))
      do i=1, mg_levels
        discretization%dmc_nflowdof(i)%dm = 0
        nullify(discretization%dmc_nflowdof(i)%ugdm)
      enddo
      dmc_ptr => discretization%dmc_nflowdof
    case(NTRANDOF)
      allocate(discretization%dmc_ntrandof(mg_levels))
      do i=1, mg_levels
        discretization%dmc_ntrandof(i)%dm = 0
        nullify(discretization%dmc_ntrandof(i)%ugdm)
      enddo
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
        call DMDASetRefinementFactor(dm_fine_ptr%dm, refine_x, refine_y, refine_z, &
                                   ierr)
        call DMDASetInterpolationType(dm_fine_ptr%dm, DMDA_Q0, ierr)
        call DMCoarsen(dm_fine_ptr%dm, option%mycomm, dmc_ptr(i)%dm, ierr)
#ifndef DMGET
        call DMCreateInterpolation(dmc_ptr(i)%dm, dm_fine_ptr%dm, &
                                   interpolation(i), PETSC_NULL_OBJECT, ierr)
#else
        call DMGetInterpolation(dmc_ptr(i)%dm, dm_fine_ptr%dm, &
                                interpolation(i), PETSC_NULL_OBJECT, ierr)
#endif
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

  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  select case(discretization%itype)
    case(STRUCTURED_GRID)
#ifndef DMGET
      call DMCreateColoring(dm_ptr%dm,IS_COLORING_GLOBAL,MATBAIJ,coloring,&
                            ierr)
#else
      call DMGetColoring(dm_ptr%dm,IS_COLORING_GLOBAL,MATBAIJ,coloring,ierr)
#endif
      ! I have set the above to use matrix type MATBAIJ, as that is what we 
      ! usually want (note: for DAs with 1 degree of freedom per grid cell, 
      ! the MATAIJ and MATBAIJ colorings should be equivalent).  What we should 
      ! eventually do here is query the type of the Jacobian matrix, but I'm 
      ! not sure of the best way to do this, as this is currently stashed in 
      ! the 'solver' object. --RTM
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

  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  call DMGlobalToLocalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec,ierr)
  call DMGlobalToLocalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec,ierr)
  
end subroutine DiscretizationGlobalToLocal

! ************************************************************************** !
!
! DiscretizationGlobalToLocalFaces: Performs global to local communication for MFD
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizationGlobalToLocalFaces(discretization,global_vec,local_vec,dm_index)

  implicit none
  

  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
    
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC,UNSTRUCTURED_GRID_MIMETIC)
      call DiscretizGlobalToLocalFacesBegin(discretization,global_vec,local_vec,dm_index)
      call DiscretizGlobalToLocalFacesEnd(discretization,global_vec,local_vec,dm_index)
  end select
  
end subroutine DiscretizationGlobalToLocalFaces


! ************************************************************************** !
!
! DiscretizationGlobalToLocalLP: Performs global to local communication for MFD
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizationGlobalToLocalLP(discretization,global_vec,local_vec,dm_index)

  implicit none
  

  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
    
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC,UNSTRUCTURED_GRID_MIMETIC)
      call DiscretizGlobalToLocalLPBegin(discretization,global_vec,local_vec,dm_index)
      call DiscretizGlobalToLocalLPEnd(discretization,global_vec,local_vec,dm_index)
  end select
  
end subroutine DiscretizationGlobalToLocalLP
  
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMLocalToGlobalBegin(dm_ptr%dm,local_vec,INSERT_VALUES,global_vec,ierr)
  call DMLocalToGlobalEnd(dm_ptr%dm,local_vec,INSERT_VALUES,global_vec,ierr)
 
end subroutine DiscretizationLocalToGlobal
  
! ************************************************************************** !
!
! DiscretizationLocalToLocal: Performs local to local communication with DM
! author: Glenn Hammond
! date: 11/14/07
!
! Some clarification:
! A "local to local" operation, in PETSc parlance, refers to communicating 
! values from a local ghosted vector (in which the ghost points are 
! irrelevant) and putting those values directly into another ghosted local 
! vector (in which those ghost points are set correctly).
! This uses the same communication pattern as a "global to local" operation, 
! but a in a "global to local", the originating vector is a PETSc global 
! vector, not a ghosted local vector.
!
! ************************************************************************** !
subroutine DiscretizationLocalToLocal(discretization,local_vec1,local_vec2,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMLocalToLocalBegin(dm_ptr%dm,local_vec1,INSERT_VALUES,local_vec2,ierr)
  call DMLocalToLocalEnd(dm_ptr%dm,local_vec1,INSERT_VALUES,local_vec2,ierr)
  
end subroutine DiscretizationLocalToLocal
  
! ************************************************************************** !
!
! DiscretizationLocalToLocalFaces: Performs local to local communication for face unknowns
! author: Daniil Svyatskiy
! date: 11/14/07
!
! ************************************************************************** !
subroutine DiscretizationLocalToLocalFaces(discretization,local_vec1,local_vec2,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC,UNSTRUCTURED_GRID_MIMETIC)
      call DiscretizLocalToLocalFacesBegin(discretization,local_vec1,local_vec2,dm_index)
      call DiscretizLocalToLocalFacesEnd(discretization,local_vec1,local_vec2,dm_index)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizationLocalToLocalFaces

! ************************************************************************** !
!
! DiscretizationLocalToLocalLP: Performs local to local communication for face unknowns
! author: Daniil Svyatskiy
! date: 11/14/07
!
! ************************************************************************** !
subroutine DiscretizationLocalToLocalLP(discretization,local_vec1,local_vec2,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC,UNSTRUCTURED_GRID_MIMETIC)
      call DiscretizLocalToLocalLPBegin(discretization,local_vec1,local_vec2,dm_index)
      call DiscretizLocalToLocalLPEnd(discretization,local_vec1,local_vec2,dm_index)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizationLocalToLocalLP
  
 
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call DMDAGlobalToNaturalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,natural_vec,ierr)
      call DMDAGlobalToNaturalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,natural_vec,ierr)
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      call VecScatterBegin(dm_ptr%ugdm%scatter_gton,global_vec,natural_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(dm_ptr%ugdm%scatter_gton,global_vec,natural_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)       
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call DMDANaturalToGlobalBegin(dm_ptr%dm,natural_vec,INSERT_VALUES,global_vec,ierr)
      call DMDANaturalToGlobalEnd(dm_ptr%dm,natural_vec,INSERT_VALUES,global_vec,ierr)
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      call VecScatterBegin(dm_ptr%ugdm%scatter_ntog,natural_vec,global_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(dm_ptr%ugdm%scatter_ntog,natural_vec,global_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMGlobalToLocalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec,ierr)
  
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMGlobalToLocalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec,ierr)
 
end subroutine DiscretizationGlobalToLocalEnd
  
! ************************************************************************** !
!
! DiscretizLocalToLocalFacesBegin: Begins Local to local communication with MFD
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizLocalToLocalFacesBegin(discretization,local_vec1,local_vec2,dm_index)


  use MFD_Aux_module

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
      call VecScatterBegin( discretization%MFD%scatter_ltol_faces,local_vec1,local_vec2 , &
                                INSERT_VALUES,SCATTER_FORWARD, ierr)
  end select
  
end subroutine DiscretizLocalToLocalFacesBegin
  
! ************************************************************************** !
!
! DiscretizLocalToLocalFacesEnd: Ends local  to local communication with DM
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizLocalToLocalFacesEnd(discretization,local_vec1,local_vec2,dm_index)

 
 use MFD_Aux_module

 implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
      call VecScatterEnd( discretization%MFD%scatter_ltol_faces,  local_vec1, local_vec2, &
                                INSERT_VALUES,SCATTER_FORWARD, ierr)
  end select
  
end subroutine DiscretizLocalToLocalFacesEnd
  
! DiscretizGlobalToLocalFacesBegin: Begins global to local communication with MFD
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizGlobalToLocalFacesBegin(discretization,global_vec,local_vec,dm_index)


  use MFD_Aux_module

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
      call VecScatterBegin( discretization%MFD%scatter_gtol_faces, global_vec, local_vec, &
                                INSERT_VALUES,SCATTER_FORWARD, ierr)
  end select
  
end subroutine DiscretizGlobalToLocalFacesBegin
  
! ************************************************************************** !
!
! DiscretizGlobalToLocalFacesEnd: Ends global to local communication with DM
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizGlobalToLocalFacesEnd(discretization,global_vec,local_vec,dm_index)

 
 use MFD_Aux_module

 implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
      call VecScatterEnd( discretization%MFD%scatter_gtol_faces, global_vec, local_vec, &
                                INSERT_VALUES,SCATTER_FORWARD, ierr)
  end select
  
end subroutine DiscretizGlobalToLocalFacesEnd

! ************************************************************************** !
!
! DiscretizLocalToLocalLPBegin: Begins Local to local communication with MFD
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizLocalToLocalLPBegin(discretization,local_vec1,local_vec2,dm_index)


  use MFD_Aux_module

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
      call VecScatterBegin( discretization%MFD%scatter_ltol_LP,local_vec1,local_vec2 , &
                                INSERT_VALUES,SCATTER_FORWARD, ierr)
  end select
  
end subroutine DiscretizLocalToLocalLPBegin
  
! ************************************************************************** !
!
! DiscretizLocalToLocalLPEnd: Ends local  to local communication with DM
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizLocalToLocalLPEnd(discretization,local_vec1,local_vec2,dm_index)

 
 use MFD_Aux_module

 implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC)
      call VecScatterEnd( discretization%MFD%scatter_ltol_LP,  local_vec1, local_vec2, &
                                INSERT_VALUES,SCATTER_FORWARD, ierr)
  end select
  
end subroutine DiscretizLocalToLocalLPEnd
  
! DiscretizGlobalToLocalLPBegin: Begins global to local communication with MFD
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizGlobalToLocalLPBegin(discretization,global_vec,local_vec,dm_index)


  use MFD_Aux_module

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC,UNSTRUCTURED_GRID_MIMETIC)
      call VecScatterBegin( discretization%MFD%scatter_gtol_LP, global_vec, local_vec, &
                                INSERT_VALUES,SCATTER_FORWARD, ierr)
  end select
  
end subroutine DiscretizGlobalToLocalLPBegin
  
! ************************************************************************** !
!
! DiscretizGlobalToLocalLPEnd: Ends global to local communication with DM
! author: Daniil Svyatskiy
! date: 07/20/10
!
!
! ************************************************************************** !
subroutine DiscretizGlobalToLocalLPEnd(discretization,global_vec,local_vec,dm_index)

 
 use MFD_Aux_module

 implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  
  select case(discretization%itype)
    case(STRUCTURED_GRID_MIMETIC,UNSTRUCTURED_GRID_MIMETIC)
      call VecScatterEnd( discretization%MFD%scatter_gtol_LP, global_vec, local_vec, &
                                INSERT_VALUES,SCATTER_FORWARD, ierr)
  end select
  
end subroutine DiscretizGlobalToLocalLPEnd


  
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMLocalToLocalBegin(dm_ptr%dm,local_vec1,INSERT_VALUES,local_vec2,ierr)

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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMLocalToLocalEnd(dm_ptr%dm,local_vec1,INSERT_VALUES,local_vec2,ierr)

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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call DMDAGlobalToNaturalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,natural_vec,ierr)
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      call VecScatterBegin(dm_ptr%ugdm%scatter_gton,global_vec,natural_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr)
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call DMDAGlobalToNaturalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,natural_vec,ierr)
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      call VecScatterEnd(dm_ptr%ugdm%scatter_gton,global_vec,natural_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)       
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call DMDANaturalToGlobalBegin(dm_ptr%dm,natural_vec,INSERT_VALUES,global_vec,ierr)
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
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call DMDANaturalToGlobalEnd(dm_ptr%dm,natural_vec,INSERT_VALUES,global_vec,ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizNaturalToGlobalEnd

! ************************************************************************** !
!
! DiscretizationUpdateTVDGhosts: Updates tvd extended ghost cell values
! author: Glenn Hammond
! date: 02/04/12
!
! ************************************************************************** !
subroutine DiscretizationUpdateTVDGhosts(discretization,global_vec, &
                                         tvd_ghost_vec)

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: tvd_ghost_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  call VecScatterBegin(discretization%tvd_ghost_scatter,global_vec, &
                       tvd_ghost_vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(discretization%tvd_ghost_scatter,global_vec, &
                       tvd_ghost_vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  
end subroutine DiscretizationUpdateTVDGhosts

! ************************************************************************** !
!
! DiscretAOApplicationToPetsc: Maps application ordering to petsc
! author: Glenn Hammond
! date: 10/12/12
!
! ************************************************************************** !
subroutine DiscretAOApplicationToPetsc(discretization,int_array)

  implicit none
  
#include "finclude/petscao.h"  
  
  type(discretization_type) :: discretization
  PetscInt :: int_array(:)
  PetscErrorCode :: ierr
  
  AO :: ao
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      call DMDAGetAO(discretization%dm_1dof,ao,ierr)
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      ao = discretization%grid%unstructured_grid%ao_natural_to_petsc
  end select
  call AOApplicationToPetsc(ao,size(int_array),int_array,ierr)
  
end subroutine DiscretAOApplicationToPetsc

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
  PetscInt :: i
    
  if (.not.associated(discretization)) return
      
  select case(discretization%itype)
    case(STRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      if (discretization%dm_1dof%dm /= 0) &
        call DMDestroy(discretization%dm_1dof%dm,ierr)
      discretization%dm_1dof%dm = 0
      if (discretization%dm_nflowdof%dm /= 0) &
        call DMDestroy(discretization%dm_nflowdof%dm,ierr)
      discretization%dm_nflowdof%dm = 0
      if (discretization%dm_ntrandof%dm /= 0) &
        call DMDestroy(discretization%dm_ntrandof%dm,ierr)
      discretization%dm_ntrandof%dm = 0
      if (associated(discretization%dmc_nflowdof)) then
        do i=1,size(discretization%dmc_nflowdof)
          call DMDestroy(discretization%dmc_nflowdof(i)%dm,ierr)
        enddo
        deallocate(discretization%dmc_nflowdof)
        nullify(discretization%dmc_nflowdof)
      endif
      if (associated(discretization%dmc_ntrandof)) then
        do i=1,size(discretization%dmc_ntrandof)
          call DMDestroy(discretization%dmc_ntrandof(i)%dm,ierr)
        enddo
        deallocate(discretization%dmc_ntrandof)
        nullify(discretization%dmc_ntrandof)
      endif
    case(UNSTRUCTURED_GRID,UNSTRUCTURED_GRID_MIMETIC)
      if (associated(discretization%dm_1dof%ugdm)) &
        call UGridDMDestroy(discretization%dm_1dof%ugdm)
      if (associated(discretization%dm_nflowdof%ugdm)) &
        call UGridDMDestroy(discretization%dm_nflowdof%ugdm)
      if (associated(discretization%dm_ntrandof%ugdm)) &
        call UGridDMDestroy(discretization%dm_ntrandof%ugdm)
  end select
  if (associated(discretization%dm_1dof)) &
    deallocate(discretization%dm_1dof)
  nullify(discretization%dm_1dof)
  if (associated(discretization%dm_nflowdof)) &
    deallocate(discretization%dm_nflowdof)
  nullify(discretization%dm_nflowdof)
  if (associated(discretization%dm_ntrandof)) &
    deallocate(discretization%dm_ntrandof)
  nullify(discretization%dm_ntrandof)


  if (associated(discretization%MFD)) &
        call MFDAuxDestroy(discretization%MFD) 
  nullify(discretization%MFD)

  if (discretization%tvd_ghost_scatter /= 0) &
    call VecScatterDestroy(discretization%tvd_ghost_scatter)
  
  ! solely nullify grid since destroyed in patch
  call GridDestroy(discretization%grid)
  
  deallocate(discretization)
  nullify(discretization)
  
end subroutine DiscretizationDestroy
 
end module Discretization_module
