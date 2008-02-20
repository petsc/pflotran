module Grid_module

  use Structured_Grid_module
  use Unstructured_Grid_module
  use Connection_module
 
  implicit none

  private
 
#include "definitions.h"

! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"

  type, public :: grid_type
  
    PetscInt :: igrid  ! type of grid (e.g. structured, unstructured, etc.)
    
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
   
    !nL2G :  not collective, local processor: local  =>  ghosted local  
    !nG2L :  not collective, local processor:  ghosted local => local  
    !nG2N :  collective,  ghosted local => global index , used for   
    !                     matsetvaluesblocked ( not matsetvaluesblockedlocal)  
    !nL2A :   collective, local => natural index, used for initialization   
    !                              and source/sink setup  
    PetscInt, pointer :: nL2G(:), nG2L(:), nL2A(:), nG2N(:)
    PetscInt, pointer :: nG2A(:)
    
    PetscReal, pointer :: x(:), y(:), z(:)
    
    PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
    PetscReal :: origin(3)

    Vec :: volume 
    
    DM :: dm_1_dof, dm_nphase_dof, dm_3np_dof, dm_nflowdof, dm_nphancomp_dof, &
          dm_nphanspec_dof, dm_nphanspecncomp_dof, dm_var_dof , &
          dm_ntrandof
    
    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash_bins
    
    PetscInt :: igeom
    type(structured_grid_type), pointer :: structured_grid
    type(unstructured_grid_type), pointer :: unstructured_grid
    
    type(connection_list_type), pointer :: internal_connection_list

  end type grid_type


  public :: GridCreate, &
            GridDestroy, &
            GridComputeInternalConnect, &
            GridCreateVector, &
            GridDuplicateVector, &         
            GridCreateJacobian, &
            GridCreateColoring, &
            GridGlobalToLocal, &
            GridLocalToGlobal, &
            GridLocalToLocal, &
            GridGlobalToNatural, &
            GridNaturalToGlobal, &
            GridGlobalToLocalBegin, &
            GridGlobalToLocalEnd, &
            GridLocalToLocalBegin, &
            GridLocalToLocalEnd, &
            GridGlobalToNaturalBegin, &
            GridGlobalToNaturalEnd, &
            GridNaturalToGlobalBegin, &
            GridNaturalToGlobalEnd, &
            GridMapIndices, &
            GridCreateDMs, &
            GridComputeSpacing, &
            GridComputeCoordinates, &
            GridComputeVolumes, &
            GridLocalizeRegions, &
            GridPopulateConnection, &
            GridCopyIntegerArrayToPetscVec, &
            GridCopyRealArrayToPetscVec, &
            GridCopyPetscVecToIntegerArray, &
            GridCopyPetscVecToRealArray, &
            GridCreateNaturalToGhostedHash, &
            GridDestroyHashTable, &
            GridGetLocalGhostedIdFromHash
  
contains

! ************************************************************************** !
!
! GridCreate: Creates a structured or unstructured grid
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function GridCreate(igeom_)

  implicit none
  
  PetscInt :: igeom_
  
  type(grid_type), pointer :: GridCreate
  
  type(grid_type), pointer :: grid
  
  allocate(grid)
  call initGrid(grid)
  grid%igeom = igeom_
  
  select case(igeom_)
    case(1,2,3)
      grid%igrid = STRUCTURED
      allocate(grid%structured_grid)
      call StructuredGridInit(grid%structured_grid)
      select case(igeom_)
        case(1)
          grid%structured_grid%igeom = STRUCTURED_CARTESIAN
        case(2)
          grid%structured_grid%igeom = STRUCTURED_CYLINDRICAL
        case(3)
          grid%structured_grid%igeom = STRUCTURED_SPHERICAL
      end select
    case(0)
      grid%igrid = UNSTRUCTURED
      allocate(grid%unstructured_grid)
      call UnstructuredGridInit(grid%unstructured_grid)
  end select
  
  GridCreate => grid

end function GridCreate

! ************************************************************************** !
!
! initGrid: Initializes the abstract grid object
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine initGrid(grid)

  implicit none
  
  type(grid_type) :: grid
  
  grid%igeom = 0
  grid%igrid = 0

  nullify(grid%structured_grid)
  nullify(grid%unstructured_grid)

  nullify(grid%internal_connection_list)

  nullify(grid%nL2G)
  nullify(grid%nG2L)
  nullify(grid%nL2A)
  nullify(grid%nG2N)
  nullify(grid%nG2A)

  nullify(grid%x)
  nullify(grid%y)
  nullify(grid%z)

  grid%x_min = 1.d20
  grid%x_max = -1.d20
  grid%y_min = 1.d20
  grid%y_max = -1.d20
  grid%z_min = 1.d20
  grid%z_max = -1.d20

  grid%origin = 0.d0

  grid%nmax = 0
  grid%nlmax = 0 
  grid%ngmax = 0
  
  ! nullify DM pointers
  grid%dm_1_dof = 0
  grid%dm_nphase_dof = 0
  grid%dm_3np_dof = 0
  grid%dm_nflowdof = 0
  grid%dm_ntrandof = 0
  grid%dm_nphancomp_dof = 0
  grid%dm_nphanspec_dof = 0
  grid%dm_nphanspecncomp_dof = 0
  grid%dm_var_dof = 0
  
  grid%volume = 0
  
  nullify(grid%hash)
  grid%num_hash_bins = 1000

end subroutine initGrid

! ************************************************************************** !
!
! GridCreateDMs: creates distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
subroutine GridCreateDMs(grid,option)
      
  use Option_module    
      
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
      
  PetscInt :: ndof
  PetscInt, parameter :: stencil_width = 1
  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------
  ! Generate the DA objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call GridCreateDM(grid,grid%dm_1_dof,ndof,stencil_width)
  
  if (option%nflowdof > 0) then
    ndof = option%nphase
    call GridCreateDM(grid,grid%dm_nphase_dof,ndof,stencil_width)

    ndof = option%nflowdof
    call GridCreateDM(grid,grid%dm_nflowdof,ndof,stencil_width)

    select case(option%iflowmode) 
      case(MPH_MODE)
        ndof = (option%nflowdof+1)*(2+7*option%nphase + 2*option%nspec*option%nphase)
        call GridCreateDM(grid,grid%dm_var_dof,ndof,stencil_width)
    end select
  endif
  
  if (option%ntrandof > 0) then
    ndof = option%ntrandof
    call GridCreateDM(grid,grid%dm_ntrandof,ndof,stencil_width)
  endif

  select case(grid%igrid)
    case(STRUCTURED)
      grid%nlmax = grid%structured_grid%nlmax
      grid%ngmax = grid%structured_grid%ngmax
    case(UNSTRUCTURED)
  end select

  ! allocate coordinate arrays  
  allocate(grid%x(grid%ngmax))
  allocate(grid%y(grid%ngmax))
  allocate(grid%z(grid%ngmax))
  
end subroutine GridCreateDMs

! ************************************************************************** !
!
! GridCreateDM: creates a distributed, parallel mesh/grid
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
subroutine GridCreateDM(grid,dm,ndof,stencil_width)

  implicit none
  
  type(grid_type) :: grid
  DM :: dm
  PetscInt :: ndof
  PetscInt :: stencil_width
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructuredGridCreateDA(grid%structured_grid,dm,ndof, &
                                  stencil_width)
    case(UNSTRUCTURED)
  end select

end subroutine GridCreateDM

! ************************************************************************** !
!
! GridCreateVector: Creates a global PETSc vector
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridCreateVector(grid,dm_index,vector,vector_type)

  implicit none
  
  type(grid_type) :: grid
  PetscInt :: dm_index
  Vec :: vector
  PetscInt :: vector_type
  
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructuredGridCreateVecFromDA(dm_ptr,vector,vector_type)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridCreateVector

! ************************************************************************** !
!
! GridDuplicateVector: Creates a global PETSc vector
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridDuplicateVector(grid,vector1,vector2)

  implicit none
  
  type(grid_type) :: grid
  Vec :: vector1
  Vec :: vector2
  
  PetscErrorCode :: ierr
  
  select case(grid%igrid)
    case(STRUCTURED,UNSTRUCTURED)
      call VecDuplicate(vector1,vector2,ierr)
  end select
  
end subroutine GridDuplicateVector

! ************************************************************************** !
!
! GridGetDMPtrFromIndex: Returns the integer pointer for the DM referenced
! author: Glenn Hammond
! date: 02/08/08
!
! ************************************************************************** !
function GridGetDMPtrFromIndex(grid,dm_index)

  implicit none
  
  type(grid_type) :: grid
  PetscInt :: dm_index
  
  DM :: GridGetDMPtrFromIndex
  
  select case (dm_index)
    case(ONEDOF)
      GridGetDMPtrFromIndex = grid%dm_1_dof
    case(NPHASEDOF)
      GridGetDMPtrFromIndex = grid%dm_nphase_dof
    case(THREENPDOF)
      GridGetDMPtrFromIndex = grid%dm_3np_dof
    case(NFLOWDOF)
      GridGetDMPtrFromIndex = grid%dm_nflowdof
    case(NTRANDOF)
      GridGetDMPtrFromIndex = grid%dm_ntrandof
    case(NPHANCOMPDOF)
      GridGetDMPtrFromIndex = grid%dm_nphancomp_dof
    case(NPHANSPECDOF)
      GridGetDMPtrFromIndex = grid%dm_nphanspec_dof
    case(NPHANSPECNCOMPDOF)
      GridGetDMPtrFromIndex = grid%dm_nphanspecncomp_dof
    case(VARDOF)
      GridGetDMPtrFromIndex = grid%dm_var_dof

  end select  
  
end function GridGetDMPtrFromIndex

! ************************************************************************** !
!
! GridComputeInternalConnect: computes internal connectivity of a grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
subroutine GridComputeInternalConnect(grid,option)

  use Connection_module
  use Option_module    
    
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  type(connection_type), pointer :: connection
  
  select case(grid%igrid)
    case(STRUCTURED)
      connection => &
        StructGridComputeInternConnect(grid%structured_grid,option)
    case(UNSTRUCTURED) 
      connection => &
        UnstGridComputeInternConnect(grid%unstructured_grid,option)
  end select
  
  allocate(grid%internal_connection_list)
  call ConnectionInitList(grid%internal_connection_list)
  call ConnectionAddToList(connection,grid%internal_connection_list)
  
end subroutine GridComputeInternalConnect

! ************************************************************************** !
!
! GridPopulateConnection: computes connectivity coupler to a grid
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine GridPopulateConnection(grid,connection,iface,iconn,cell_id_local)

  use Connection_module
  use Structured_Grid_module
  
  implicit none
 
  type(grid_type) :: grid
  type(connection_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: cell_id_local
  
  PetscInt :: cell_id_ghosted
  
  cell_id_ghosted = grid%nL2G(cell_id_local)
  ! Use ghosted index to access dx, dy, dz because we have
  ! already done a global-to-local scatter for computing the
  ! interior node connections.
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructGridPopulateConnection(grid%structured_grid,connection, &
                                        iface,iconn,cell_id_ghosted)
    case(UNSTRUCTURED)
  end select

end subroutine GridPopulateConnection

! ************************************************************************** !
!
! GridMapIndices: maps global, local and natural indices of cells 
!                 to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridMapIndices(grid)

  implicit none
  
  type(grid_type) :: grid
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructuredGridMapIndices(grid%structured_grid,grid%nG2L,grid%nL2G, &
                                    grid%nL2A,grid%nG2A,grid%nG2N)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridMapIndices

! ************************************************************************** !
!
! GridComputeSpacing: Computes grid spacing (only for structured grid
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine GridComputeSpacing(grid)

  implicit none
  
  type(grid_type) :: grid
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructuredGridComputeSpacing(grid%structured_grid,grid%nL2A)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridComputeSpacing

! ************************************************************************** !
!
! GridComputeCoordinates: Computes x,y,z coordinates of grid cells
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridComputeCoordinates(grid,option)

  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructuredGridComputeCoord(grid%structured_grid,option, &
                                      grid%origin,grid%x,grid%y,grid%z, &
                                      grid%x_min,grid%x_max,grid%y_min, &
                                      grid%y_max,grid%z_min,grid%z_max)
    case(UNSTRUCTURED)
  end select

end subroutine GridComputeCoordinates

! ************************************************************************** !
!
! GridComputeVolumes: Computes the volumes of cells in structured grid
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GridComputeVolumes(grid,option)

  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructuredGridComputeVolumes(grid%structured_grid,option, &
                                        grid%nL2G,grid%volume)
    case(UNSTRUCTURED)
  end select

end subroutine GridComputeVolumes

! ************************************************************************** !
!
! GridCreateJacobian: Creates Jacobian matrix associated with grid
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridCreateJacobian(grid,dm_index,Jacobian,option)

  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  PetscInt :: dm_index
  Mat :: Jacobian
  type(option_type) :: option

  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
    
  select case(grid%igrid)
    case(STRUCTURED)
      call StructuredGridCreateJacobian(dm_ptr,Jacobian,option)
    case(UNSTRUCTURED)
  end select

end subroutine GridCreateJacobian

! ************************************************************************** !
!
! GridCreateColoring: Creates ISColoring for grid
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridCreateColoring(grid,dm_index,option,coloring)

  use Option_module
  
  implicit none

#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
  
  type(grid_type) :: grid
  PetscInt :: dm_index
  type(option_type) :: option
  ISColoring :: coloring

  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
    
  select case(grid%igrid)
    case(STRUCTURED)
      call StructuredGridCreateColoring(dm_ptr,option,coloring)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridCreateColoring

! ************************************************************************** !
!
! GridGlobalToLocal: Performs global to local communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine GridGlobalToLocal(grid,global_vec,local_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
    
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridGlobalToLocal(dm_ptr,global_vec,local_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridGlobalToLocal
  
! ************************************************************************** !
!
! GridLocalToGlobal: Performs local to global communication with DM
! author: Glenn Hammond
! date: 1/02/08
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine GridLocalToGlobal(grid,local_vec,global_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: local_vec
  Vec :: global_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridLocalToGlobal(dm_ptr,local_vec,global_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridLocalToGlobal
  
! ************************************************************************** !
!
! GridLocalToLocal: Performs local to local communication with DM
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine GridLocalToLocal(grid,local_vec1,local_vec2,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridLocalToLocal(dm_ptr,local_vec1,local_vec2)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridLocalToLocal
  
! ************************************************************************** !
!
! GridGlobalToNatural: Performs global to natural communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridGlobalToNatural(grid,global_vec,natural_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridGlobalToNatural(dm_ptr,global_vec,natural_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridGlobalToNatural

! ************************************************************************** !
!
! GridNaturalToGlobal: Performs natural to global communication with DM
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine GridNaturalToGlobal(grid,natural_vec,global_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridNaturalToGlobal(dm_ptr,natural_vec,global_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridNaturalToGlobal

! ************************************************************************** !
!
! GridGlobalToLocalBegin: Begins global to local communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine GridGlobalToLocalBegin(grid,global_vec,local_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridGlobalToLocalBegin(dm_ptr,global_vec,local_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridGlobalToLocalBegin
  
! ************************************************************************** !
!
! GridGlobalToLocalEnd: Ends global to local communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! Note that 'dm_index' should correspond to one of the macros defined 
! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
!
! ************************************************************************** !
subroutine GridGlobalToLocalEnd(grid,global_vec,local_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridGlobalToLocalEnd(dm_ptr,global_vec,local_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridGlobalToLocalEnd
  
! ************************************************************************** !
!
! GridLocalToLocalBegin: Begins local to local communication with DM
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine GridLocalToLocalBegin(grid,local_vec1,local_vec2,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridLocalToLocalBegin(dm_ptr,local_vec1,local_vec2)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridLocalToLocalBegin
  
! ************************************************************************** !
!
! GridLocalToLocalEnd: Ends local to local communication with DM
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine GridLocalToLocalEnd(grid,local_vec1,local_vec2,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridLocalToLocalEnd(dm_ptr,local_vec1,local_vec2)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridLocalToLocalEnd
  
! ************************************************************************** !
!
! GridGlobalToNaturalBegin: Begins global to natural communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridGlobalToNaturalBegin(grid,global_vec,natural_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridGlobalToNaturBegin(dm_ptr,global_vec,natural_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridGlobalToNaturalBegin

! ************************************************************************** !
!
! GridGlobalToNaturalEnd: Ends global to natural communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridGlobalToNaturalEnd(grid,global_vec,natural_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridGlobalToNaturEnd(dm_ptr,global_vec,natural_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridGlobalToNaturalEnd

! ************************************************************************** !
!
! GridNaturalToGlobalBegin: Begins natural to global communication with DM
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine GridNaturalToGlobalBegin(grid,natural_vec,global_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridNaturToGlobalBegin(dm_ptr,natural_vec,global_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridNaturalToGlobalBegin

! ************************************************************************** !
!
! GridNaturalToGlobalEnd: Ends natural to global communication with DM
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine GridNaturalToGlobalEnd(grid,natural_vec,global_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  DM :: dm_ptr
  
  dm_ptr = GridGetDMPtrFromIndex(grid,dm_index)
  
  select case(grid%igrid)
    case(STRUCTURED)
      call StructureGridNaturToGlobalEnd(dm_ptr,natural_vec,global_vec)
    case(UNSTRUCTURED)
  end select
  
end subroutine GridNaturalToGlobalEnd

! ************************************************************************** !
!
! GridLocalizeRegions: Resticts regions to cells local to processor
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine GridLocalizeRegions(region_list,grid,option)

  use Option_module
  use Region_module

  implicit none
  
  type(region_list_type), pointer :: region_list
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  type(region_type), pointer :: region
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, local_ghosted_id, local_id
  
  region => region_list%first
  do
  
    if (.not.associated(region)) exit
    
    if (.not.associated(region%cell_ids) .and. &
        region%i1 > 0 .and. region%i2 > 0 .and. &
        region%j1 > 0 .and. region%j2 > 0 .and. &
        region%k1 > 0 .and. region%k2 > 0) then

      ! convert indexing from global (entire domain) to local processor
      region%i1 = region%i1 - grid%structured_grid%nxs
      region%i2 = region%i2 - grid%structured_grid%nxs
      region%j1 = region%j1 - grid%structured_grid%nys
      region%j2 = region%j2 - grid%structured_grid%nys
      region%k1 = region%k1 - grid%structured_grid%nzs
      region%k2 = region%k2 - grid%structured_grid%nzs
        
      ! clip region to within local processor domain
      region%i1 = max(region%i1,1)

      region%i2 = min(region%i2,grid%structured_grid%nlx)
      region%j1 = max(region%j1,1)
      region%j2 = min(region%j2,grid%structured_grid%nly)
      region%k1 = max(region%k1,1)
      region%k2 = min(region%k2,grid%structured_grid%nlz)
       
      count = 0  
      if (region%i1 <= region%i2 .and. &
          region%j1 <= region%j2 .and. &
          region%k1 <= region%k2) then
        region%num_cells = (region%i2-region%i1+1)* &
                           (region%j2-region%j1+1)* &
                           (region%k2-region%k1+1)
        allocate(region%cell_ids(region%num_cells))
        region%cell_ids = 0
        do k=region%k1,region%k2
          do j=region%j1,region%j2
            do i=region%i1,region%i2
              count = count + 1
              region%cell_ids(count) = &
                     i + (j-1)*grid%structured_grid%nlx + &
                     (k-1)*grid%structured_grid%nlxy
            enddo
          enddo
        enddo
      else
        region%num_cells = 0
      endif
     
      if (count /= region%num_cells) &
        call printErrMsg(option,"Mismatch in number of cells in block region")

    else
#if 0
      allocate(temp_int_array(region%num_cells))
      temp_int_array = 0
      if (grid%igrid == STRUCTURED) then
        do count=1,region%num_cells
          i = mod(region%cell_ids(count),grid%structured_grid%nx) - &
                grid%structured_grid%nxs
          j = mod((region%cell_ids(count)-1)/grid%structured_grid%nx, &
                  grid%structured_grid%ny)+1 - &
                grid%structured_grid%nys
          k = ((region%cell_ids(count)-1)/grid%structured_grid%nxy)+1 - &
                grid%structured_grid%nzs
          if (i > 0 .and. i <= grid%structured_grid%nlx .and. &
              j > 0 .and. j <= grid%structured_grid%nly .and. &
              k > 0 .and. k <= grid%structured_grid%nlz) then
            temp_int_array(local_count) = &
                i + (j-1)*grid%structured_grid%nlx + &
                (k-1)*grid%structured_grid%nlxy
            local_count = local_count + 1
          endif
        enddo
      else
        call GridCreateNaturalToGhostedHash(grid,option)
        do count=1,region%num_cells
          local_ghosted_id = UnstructGridGetGhostIdFromHash( &
                                                    grid%unstructured_grid, &
                                                    region%cell_ids(count))
          if (local_ghosted_id > -1) then
            local_id = grid%nG2L(local_ghosted_id)
            if (local_id > -1) then
              temp_int_array(local_count) = local_id
              local_count = local_count + 1
            endif
          endif
        enddo
        call GridDestroyHashTable(grid)
      endif
      if (local_count /= region%num_cells) then
        deallocate(region%cell_ids)
        allocate(region%cell_ids(local_count))
        region%cell_ids(1:local_count) = temp_int_array(1:local_count)
      endif
      deallocate(temp_int_array)
#endif        
    endif
    
    if (region%num_cells == 0 .and. associated(region%cell_ids)) &
      deallocate(region%cell_ids)
    region => region%next
    
  enddo

end subroutine GridLocalizeRegions

! ************************************************************************** !
!
! GridCopyIntegerArrayToPetscVec: Copies values from an integer array into a 
!                                 PETSc Vec
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyIntegerArrayToPetscVec(array,vector,num_values)

  implicit none
  
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyIntegerArrayToPetscVec

! ************************************************************************** !
!
! GridCopyRealArrayToPetscVec: Copies values from an integer array into a 
!                              PETSc Vec
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyRealArrayToPetscVec(array,vector,num_values)

  implicit none
  
  PetscReal :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyRealArrayToPetscVec

! ************************************************************************** !
!
! GridCopyPetscVecToIntegerArray: Copies values from a PETSc Vec to an  
!                                 integer array
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyPetscVecToIntegerArray(array,vector,num_values)

  implicit none
  
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  do i=1,num_values
    array(i) = int(vec_ptr(i)+1.d-4)
  enddo
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyPetscVecToIntegerArray

! ************************************************************************** !
!
! GridCopyPetscVecToRealArray: Copies values from a PETSc Vec to an integer 
!                              array
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyPetscVecToRealArray(array,vector,num_values)

  implicit none
  
  PetscReal :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  array(1:num_values) = vec_ptr(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyPetscVecToRealArray

! ************************************************************************** !
!
! GridCreateNaturalToGhostedHash: Creates a hash table for looking up the  
!                                 local ghosted id of a natural id, if it 
!                                 exists
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine GridCreateNaturalToGhostedHash(grid,option)

  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: local_ghosted_id, natural_id
  PetscInt :: num_in_hash, num_ids_per_hash, hash_id, id, ierr
  PetscInt :: max_num_ids_per_hash = 0
  PetscInt, pointer :: hash(:,:,:), temp_hash(:,:,:)

  if (associated(grid%hash)) return

  ! initial guess of 10% of ids per hash
  ! must be at least 5 so that reallocation (*1.2) works below
  num_ids_per_hash = max(grid%nlmax/(grid%num_hash_bins/10),5)

  allocate(hash(2,0:num_ids_per_hash,grid%num_hash_bins))
  hash(:,:,:) = 0

  ! recall that natural ids are zero-based
  do local_ghosted_id = 1, grid%ngmax
    natural_id = grid%nG2A(local_ghosted_id)+1
    hash_id = mod(natural_id,grid%num_hash_bins)+1 
    num_in_hash = hash(1,0,hash_id)
    num_in_hash = num_in_hash+1
    if (num_in_hash > max_num_ids_per_hash) max_num_ids_per_hash = num_in_hash
    ! if a hash runs out of space reallocate
    if (num_in_hash > num_ids_per_hash) then 
      allocate(temp_hash(2,0:num_ids_per_hash,0:grid%num_hash_bins))
      ! copy old hash
      temp_hash(1:2,0:num_ids_per_hash,grid%num_hash_bins) = &
                             hash(1:2,0:num_ids_per_hash,grid%num_hash_bins)
      deallocate(hash)
      ! recompute hash 20% larger
      num_ids_per_hash = int(dble(num_ids_per_hash)*1.2)
      allocate(hash(1:2,0:num_ids_per_hash,grid%num_hash_bins))
      ! copy old to new
      do hash_id = 1, grid%num_hash_bins
        do id = 1, temp_hash(1,0,hash_id)
          hash(1:2,id,hash_id) = temp_hash(1:2,id,hash_id)
        enddo
        hash(1,0,hash_id) = temp_hash(1,0,hash_id)
      enddo
      deallocate(temp_hash)
    endif
    hash(1,0,hash_id) = num_in_hash
    hash(1,num_in_hash,hash_id) = natural_id
    hash(2,num_in_hash,hash_id) = local_ghosted_id
  enddo

  grid%hash => hash
  
!  call GridPrintHashTable(grid)
  call mpi_allreduce(max_num_ids_per_hash,num_in_hash,ONE_INTEGER,MPI_INTEGER, &
                     MPI_MAX,PETSC_COMM_WORLD,ierr)
  if (option%myrank == 0) print *, 'max_num_ids_per_hash:', num_in_hash

end subroutine GridCreateNaturalToGhostedHash

! ************************************************************************** !
!
! GetLocalIdFromNaturalId: Returns the local id corresponding to a natural
!                          id or 0, if the natural id is off-processor
! WARNING: Extremely inefficient for large jobs
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalIdFromNaturalId(grid,natural_id)

  implicit none

  type(grid_type) :: grid

  PetscInt :: natural_id, local_id
  
  do local_id = 1, grid%nlmax
    if (natural_id == grid%nL2A(local_id)+1) then
      GridGetLocalIdFromNaturalId = local_id
      return
    endif
  enddo
  GridGetLocalIdFromNaturalId = 0

end function GridGetLocalIdFromNaturalId

! ************************************************************************** !
!
! GridGetLocalGhostedIdFromNatId: Returns the local ghosted id corresponding 
!                                 to a natural id or 0, if the natural id 
!                                 is off-processor
! WARNING: Extremely inefficient for large jobs
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalGhostedIdFromNatId(grid,natural_id)

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  
  PetscInt :: local_ghosted_id
  
  do local_ghosted_id = 1, grid%ngmax
    if (natural_id == grid%nG2A(local_ghosted_id)+1) then
      GridGetLocalGhostedIdFromNatId = local_ghosted_id
      return 
    endif
  enddo
  GridGetLocalGhostedIdFromNatId = 0

end function GridGetLocalGhostedIdFromNatId

! ************************************************************************** !
!
! GridGetLocalGhostedIdFromHash: Returns the local ghosted id of a natural 
!                                id, if it exists.  Otherwise 0 is returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalGhostedIdFromHash(grid,natural_id)

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  
  PetscInt :: hash_id, id

  GridGetLocalGhostedIdFromHash = 0
  hash_id = mod(natural_id,grid%num_hash_bins)+1 
  do id = 1, grid%hash(1,0,hash_id)
    if (grid%hash(1,id,hash_id) == natural_id) then
      GridGetLocalGhostedIdFromHash = grid%hash(2,id,hash_id)
      return
    endif
  enddo

end function GridGetLocalGhostedIdFromHash

! ************************************************************************** !
!
! GridDestroyHashTable: Deallocates the hash table
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine GridDestroyHashTable(grid)

  implicit none

  type(grid_type), pointer :: grid
  
  if (associated(grid%hash)) deallocate(grid%hash)
  nullify(grid%hash)
  grid%num_hash_bins = 100

end subroutine GridDestroyHashTable

! ************************************************************************** !
!
! UnstructGridPrintHashTable: Prints the hashtable for viewing
! author: Glenn Hammond
! date: 03/09/07
!
! ************************************************************************** !
subroutine GridPrintHashTable(grid)

  implicit none

  type(grid_type) :: grid
  
  PetscInt :: ihash, id, fid

  fid = 87 
  open(fid,file='hashtable.dat',action='write')
  do ihash=1,grid%num_hash_bins
    write(fid,'(a4,i3,a,i5,a2,x,200(i6,x))') 'Hash',ihash,'(', &
                         grid%hash(1,0,ihash), &
                         '):', &
                         (grid%hash(1,id,ihash),id=1,grid%hash(1,0,ihash))
  enddo
  close(fid)

end subroutine GridPrintHashTable

! ************************************************************************** !
!
! GridDestroy: Deallocates a grid
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine GridDestroy(grid)

  implicit none
  
  type(grid_type), pointer :: grid
    
  if (.not.associated(grid)) return
      
  if (associated(grid%nL2G)) deallocate(grid%nL2G)
  nullify(grid%nL2G)
  if (associated(grid%nG2L)) deallocate(grid%nG2L)
  nullify(grid%nG2L)
  if (associated(grid%nL2A)) deallocate(grid%nL2A)
  nullify(grid%nL2A)
  ! Since nG2N is actually a pointer to a Petsc IS, cannot deallocate
  ! unless we change it.
!  if (associated(grid%nG2N)) deallocate(grid%nG2N)
!  nullify(grid%nG2N)
  if (associated(grid%nG2A)) deallocate(grid%nG2A)
  nullify(grid%nG2A)

  if (associated(grid%x)) deallocate(grid%x)
  nullify(grid%x)
  if (associated(grid%y)) deallocate(grid%y)
  nullify(grid%y)
  if (associated(grid%z)) deallocate(grid%z)
  nullify(grid%z)
  
  if (associated(grid%hash)) call GridDestroyHashTable(grid)

  select case(grid%igrid)
    case(STRUCTURED)
      if (grid%dm_1_dof /= 0) &
        call StructuredGridDestroyDA(grid%dm_1_dof)
      if (grid%dm_nphase_dof /= 0) &
        call StructuredGridDestroyDA(grid%dm_nphase_dof)
      if (grid%dm_3np_dof /= 0) &
        call StructuredGridDestroyDA(grid%dm_3np_dof)
      if (grid%dm_nflowdof /= 0) &
        call StructuredGridDestroyDA(grid%dm_nflowdof)
      if (grid%dm_nphancomp_dof /= 0) &
        call StructuredGridDestroyDA(grid%dm_nphancomp_dof)
      if (grid%dm_nphanspec_dof /= 0) &
        call StructuredGridDestroyDA(grid%dm_nphanspec_dof)
      if (grid%dm_nphanspecncomp_dof /= 0) &
        call StructuredGridDestroyDA(grid%dm_nphanspecncomp_dof)
      if (grid%dm_var_dof /= 0) &
        call StructuredGridDestroyDA(grid%dm_var_dof)
      if (grid%dm_ntrandof /= 0) &
        call StructuredGridDestroyDA(grid%dm_ntrandof)
    case(UNSTRUCTURED)
  end select
  
  call UnstructuredGridDestroy(grid%unstructured_grid)    
  call StructuredGridDestroy(grid%structured_grid)
                                           
  call ConnectionDestroyList(grid%internal_connection_list)

end subroutine GridDestroy
  
end module Grid_module
