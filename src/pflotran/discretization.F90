module Discretization_module

  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module
  use AMR_Grid_Module

  implicit none

  private
 
#include "definitions.h"

#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"

  type, public :: discretization_type
    PetscInt :: itype  ! type of discretization (e.g. structured, unstructured, etc.)
    character(len=MAXWORDLENGTH) :: ctype
    type(grid_type), pointer :: grid  ! pointer to a grid object
    type(amrgrid_type), pointer :: amrgrid  ! pointer to an amr grid object
    DM :: dm_1_dof, dm_nflowdof, dm_ntrandof
    DM, pointer :: dmc_nflowdof(:), dmc_ntrandof(:)
      ! Arrays containing hierarchy of coarsened DMs, for use with Galerkin 
      ! multigrid.  Element i of each array is a *finer* DM than element i-1.

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

  ! nullify DM pointers
  discretization%dm_1_dof = 0
  discretization%dm_nflowdof = 0
  discretization%dm_ntrandof = 0
  
  DiscretizationCreate => discretization

end function DiscretizationCreate

! ************************************************************************** !
!
! DiscretizationRead: Reads a discretization from the input file
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine DiscretizationRead(discretization,fid,option)

  use Fileio_module
  use Option_module
  use AMR_Grid_module

  implicit none
  
  type(option_type) :: option
  type(discretization_type) :: discretization
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid
  type(structured_grid_type), pointer :: str_grid
  type(amrgrid_type), pointer :: amrgrid
  PetscInt :: length
  PetscInt :: nx, ny, nz
  PetscErrorCode :: ierr

  nx = 0
  ny = 0
  nz = 0

  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','GRID', ierr)   
    length = len_trim(word)
    call fiCharsToUpper(word,length)
      
    select case(trim(word))
      case('TYPE')
        call fiReadWord(string,discretization%ctype,.true.,ierr)
        call fiErrorMsg(option%myrank,'type','GRID', ierr)   
        length = len_trim(discretization%ctype)
        call fiCharsToLower(discretization%ctype,length)
        select case(trim(discretization%ctype))
          case('structured')
            discretization%itype = STRUCTURED_GRID
          case('unstructured')
            discretization%itype = UNSTRUCTURED_GRID
          case('amr')
            discretization%itype = AMR_GRID
          case default
            call printErrMsg(option,'Discretization type: '//trim(discretization%ctype)//' not recognized.')
        end select    
      case('NXYZ')
        call fiReadInt(string,nx,ierr)
        call fiErrorMsg(option%myrank,'nx','GRID', ierr)
        call fiReadInt(string,ny,ierr)
        call fiErrorMsg(option%myrank,'ny','GRID', ierr)
        call fiReadInt(string,nz,ierr)
        call fiErrorMsg(option%myrank,'nz','GRID', ierr)
      case('FILE')
      case('END')
        exit
    end select 
  
  enddo  

  select case(discretization%itype)
    case(UNSTRUCTURED_GRID,STRUCTURED_GRID)
      grid => GridCreate()
      select case(discretization%itype)
        case(STRUCTURED_GRID)      
          if (nx*ny*nz == 0) call printErrMsg(option,'NXYZ not set correctly for structured grid.')
          str_grid => StructuredGridCreate()
          str_grid%nx = nx
          str_grid%ny = ny
          str_grid%nz = nz
          str_grid%nxy = str_grid%nx*str_grid%ny
          str_grid%nmax = str_grid%nxy*str_grid%nz
          grid%structured_grid => str_grid
          grid%nmax = str_grid%nmax
      end select
      discretization%grid => grid
      grid%itype = discretization%itype
      grid%ctype = discretization%ctype
    case(AMR_GRID)
       amrgrid=>discretization%amrgrid
       call AMRGridInitialize(amrgrid)
  end select

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

  !-----------------------------------------------------------------------
  ! Generate the DA objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call DiscretizationCreateDM(discretization,discretization%dm_1_dof,ndof, &
                              stencil_width)
  
  if (option%nflowdof > 0) then
    ndof = option%nflowdof
    call DiscretizationCreateDM(discretization,discretization%dm_nflowdof, &
                                ndof,stencil_width)
  endif
  
  if (option%ntrandof > 0) then
    ndof = option%ntrandof
    call DiscretizationCreateDM(discretization,discretization%dm_ntrandof, &
                                ndof,stencil_width)
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
subroutine DiscretizationCreateDM(discretization,dm,ndof,stencil_width)

  implicit none
  
  type(discretization_type) :: discretization
  DM :: dm
  PetscInt :: ndof
  PetscInt :: stencil_width
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call StructuredGridCreateDA(discretization%grid%structured_grid,dm, &
                                  ndof,stencil_width)
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
      call AMRGridCreateVector(discretization%amrgrid, ndof, vector, vector_type)
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
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  MatType :: mat_type
  Mat :: Jacobian
  type(option_type) :: option
  PetscFortranAddr :: p_samr_matrix
  PetscInt :: ndof, stencilsize

  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DAGetMatrix(dm_ptr,mat_type,Jacobian,ierr)
#if (PETSC_VERSION_RELEASE == 1)
      call MatSetOption(Jacobian,MAT_KEEP_ZEROED_ROWS,ierr)
      call MatSetOption(Jacobian,MAT_COLUMN_ORIENTED,ierr)
#else
      call MatSetOption(Jacobian,MAT_KEEP_ZEROED_ROWS,PETSC_FALSE,ierr)
      call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
#endif
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
       call MatCreateShell(PETSC_COMM_WORLD, 0,0, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_NULL, Jacobian, ierr)
       call SAMRCreateMatrix(discretization%amrgrid%p_application, ndof, stencilsize, Jacobian)
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
                                             mg_levels_y, mg_levels_z)

  use Option_module
  
  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  Mat, pointer :: interpolation(:)
  PetscInt :: mg_levels_x, mg_levels_y, mg_levels_z

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
      dmc_ptr = discretization%dmc_nflowdof
    case(NTRANDOF)
      dmc_ptr = discretization%dmc_ntrandof
  end select  
   
  allocate(dmc_ptr(mg_levels))
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
        call DACoarsen(dm_fine_ptr, PETSC_COMM_WORLD, dmc_ptr(i), ierr)
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

#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
  
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
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  DM :: dm_ptr
  
  dm_ptr = DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
#if (PETSC_VERSION_RELEASE == 1)    
      call DAGlobalToLocalBegin(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
      call DAGlobalToLocalEnd(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
#else
      call DMGlobalToLocalBegin(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
      call DMGlobalToLocalEnd(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
#endif
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
#if (PETSC_VERSION_RELEASE == 1)
      call DALocalToGlobal(dm_ptr,local_vec,INSERT_VALUES,global_vec,ierr)
#else    
      call DMLocalToGlobal(dm_ptr,local_vec,INSERT_VALUES,global_vec,ierr)
#endif
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
#if (PETSC_VERSION_RELEASE == 1)    
      call DAGlobalToLocalBegin(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
#else
      call DMGlobalToLocalBegin(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
#endif      
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
#if (PETSC_VERSION_RELEASE == 1)        
      call DAGlobalToLocalEnd(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
#else
      call DMGlobalToLocalEnd(dm_ptr,global_vec,INSERT_VALUES,local_vec,ierr)
#endif      
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
#if (PETSC_VERSION_RELEASE == 0)    
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
#else
      if (discretization%dm_1_dof /= 0) &
        call DADestroy(discretization%dm_1_dof,ierr)
      discretization%dm_1_dof = 0
      if (discretization%dm_nflowdof /= 0) &
        call DADestroy(discretization%dm_nflowdof,ierr)
      discretization%dm_nflowdof = 0
      if (discretization%dm_ntrandof /= 0) &
        call DADestroy(discretization%dm_ntrandof,ierr)
      discretization%dm_ntrandof = 0
#endif
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizationDestroy
 
end module Discretization_module
