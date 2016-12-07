!  contains

!#include "pflowgrid_output.F90"
!#include "pflowgrid_new.F90"
!#include "pflowgrid_destroy.F90"
!#include "pflowgrid_setup.F90"
!#include "pflowgrid_setvel.F90"
!#include "pflowgrid_compute_xyz.F90"
!#include "pflowgrid_update_dt.F90"
!#include "pflowgrid_step.F90"
!!#include "pflowgrid_ptran_init.F90"
!#include "pflowLIQUID.F90"
!#include "pflowCOND.F90"
!#include "pflowTH.F90"
!#include "pflowTHC.F90"
!#include "pflow2PH.F90"
!#include "pflowjacobian.F90"
!#include "pflowgrid_readinput.F90"
!#include "pflowgrid_parse_cmdline.F90"
!#include "pflowgrid_misc.F90"
!#include "pflowgrid_monitors.F90"

  module pflow_grid_module
  
#include "include/finclude/petscsnes.h"
  use petscsnes
  use pflow_gridtype_module
  use pflow_read_module
  
  implicit none
  
  private
  real*8, parameter :: Pi=3.1415926D0 

  public pflowGrid
  public pflowGrid_new
  public pflowGrid_newpm
  public pflowGrid_setup
  public pflowGrid_read_input
 ! public porperm_out
  public pflowGrid_destroy
  public Pflow_allocate_Vec
  public pflowgrid_Setup_SNES
  public pflow_setup_index   ! should make private later on

#include "definitions.h"

!#include "pflow_gridtype.h"     
     
  contains

!======================================================================

!#include "pflowgrid_new.F90"

! pflowGrid_new() is a structure constructor function for the  
! pflowGrid type. It allocates memory for all allocatable arrays.
! It initializes all class members associated with grid topology (but 
! not physical geometry -- dx, dy, areas, etc. must be set up 
! separately). Note that it does not set up all of the topology of the  
! cell connections, as it is cleaner to do this in pflowGrid_setup()   
! when the geometry of the connections is calculated.
    subroutine pflowGrid_new(grid, timestep,igeom, nx, ny, nz, npx, npy, npz, &
      nphase )
  
      implicit none

      type(pflowGrid) :: grid
      type(time_stepping_context), intent(inout) :: timestep
      
      interface 
         subroutine allocate_patch_info(p_samr_hierarchy, patchlevel_info)
           use pflow_gridtype_module
           PetscFortranAddr :: p_samr_hierarchy
           type(PatchLevelInfoPtr), dimension(:), pointer :: patchlevel_info
         end subroutine allocate_patch_info
      end interface

      integer, intent(in) :: igeom, nx, ny, nz, npx, npy, npz
      integer, intent(in) :: nphase
!      integer, intent(in) :: equ_option
      integer :: n, ng, na
      integer :: i, j, k
      integer :: ierr
    
	 ! integer, pointer::ghostind(:)

      call MPI_Comm_rank(PETSC_COMM_WORLD, grid%myrank, ierr)
      call MPI_Comm_size(PETSC_COMM_WORLD, grid%commsize, ierr)
       !call pflowGrid_parse_cmdline(grid)
	   grid%use_ksp=PETSC_FALSE
	   grid%use_isoth=PETSC_FALSE
       grid%samrai_drive = PETSC_FALSE

      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-snes_mf", & 
      grid%use_matrix_free, ierr)
      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_analytical", &
      grid%use_analytical, ierr)
      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-print_hhistory", &
      grid%print_hhistory, ierr)
      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-monitor_h", &
      grid%monitor_h, ierr)
       call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_ksp", &
      grid%use_ksp, ierr)
       call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-SAMRAI", &
      grid%samrai_drive, ierr)
      
      grid%igeom = igeom
      grid%nx = nx
      grid%ny = ny
      grid%nz = nz
      grid%nxy = nx*ny
      grid%nmax = nx*ny*nz
      print *,"new grid 1:",grid%nx,grid%ny,grid%nz
      grid%npx = npx
      grid%npy = npy
      grid%npz = npz

      grid%nphase = nphase
      grid%ndof = nphase

    
   
!-----------------------------------------------------------------------
      ! Initialize some counter variables.
!-----------------------------------------------------------------------
      grid%t = 0.d0
      grid%flowsteps = 0
      grid%newtcum = 0
      grid%icutcum = 0

      !-----------------------------------------------------------------------
      ! Initialize some parameters to sensible values.  These are parameters
      ! which should be set via the command line or the input file, but it
      ! seems good practice to set them to sensible values when a pflowGrid
      ! is created.
!-----------------------------------------------------------------------
      allocate(timestep%tfac(13))
      
      timestep%tfac(1)  = 2.0d0; timestep%tfac(2)  = 2.0d0
      timestep%tfac(3)  = 2.0d0; timestep%tfac(4)  = 2.0d0
      timestep%tfac(5)  = 2.0d0; timestep%tfac(6)  = 1.8d0
      timestep%tfac(7)  = 1.6d0; timestep%tfac(8)  = 1.4d0
      timestep%tfac(9)  = 1.2d0; timestep%tfac(10) = 1.0d0
      timestep%tfac(11) = 1.0d0; timestep%tfac(12) = 1.0d0
      timestep%tfac(13) = 1.0d0      
      timestep%iaccel = 1
      timestep%dt_min = 1.d0
      timestep%dt_max = 3.1536d6 ! One-tenth of a year.
      timestep%icut_max = 16
      timestep%dpmxe = 5.d4
      timestep%dsmxe = 5.d-1

      print *,"new grid 2:",grid%nx,grid%ny,grid%nz, grid%nphase
      
      grid%use_matrix_free = 1
      grid%newton_max = 16
            
      grid%dt = 1.d0
      grid%atol = PETSC_DEFAULT_DOUBLE_PRECISION
      grid%rtol = PETSC_DEFAULT_DOUBLE_PRECISION
      grid%stol = PETSC_DEFAULT_DOUBLE_PRECISION
      grid%maxit = PETSC_DEFAULT_INTEGER
      grid%maxf = PETSC_DEFAULT_INTEGER
      
      grid%ihydrostatic = 0
      !grid%conc0 = 1.d-6
      
   
      !physical constants
      grid%gravity = 9.8068d0    ! m/s^2
	  grid%tref=50D0
 !    grid%gravity = 0d0    ! m/s^2
     
      !-----------------------------------------------------------------------
      ! Generate the DA objects that will manage communication.
      !-----------------------------------------------------------------------
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
           nx,ny,nz,npx,npy,npz,1,1, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           grid%da_1_dof,ierr)
  
   
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
           nx,ny,nz,npx,npy,npz,3*grid%nphase,1, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           grid%da_3np_dof,ierr)
  
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
           nx,ny,nz,npx,npy,npz,grid%ndof,1, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           grid%da_ndof,ierr)
 
	   
  !    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
  !             nx,ny,nz,npx,npy,npz,(grid%ndof+1)*(2+4*grid%nphase),1, &
  !             PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
  !             grid%da_var_dof,ierr)

	 
       !print DA info for each processor
       call DAView(grid%da_ndof,PETSC_VIEWER_STDOUT_WORLD,ierr)

      !-----------------------------------------------------------------------
      ! Create the vectors with parallel layout corresponding to the DA's,
      ! and, for vectors that need to be ghosted, create the corresponding
      ! ghosted vectors.
      !-----------------------------------------------------------------------

      ! 1 degree of freedom
   
     call Pflow_allocate_Vec(grid) 
     print *,"new grid 3:",grid%nx,grid%ny,grid%nz       
   !-----------------------------------------------------------------------
      ! Set up information about corners of local domain.
!-----------------------------------------------------------------------
     ! allocate space to store information about the patch objects
     call allocate_patch_info(grid%p_samr_hierarchy, grid%patchlevel_info)
   
   print *,"new grid 4:",grid%nx,grid%ny,grid%nz   
      !     if (grid%using_pflowGrid == PETSC_TRUE) &
!     allocate(grid%vvl_loc(locpat%nconn*grid%nphase))

      ! I don't like having a fixed number of boundary condition regions.
      ! Memory for these arrays ought to allocated by parsing the input file
      ! to determine the number of regions.  This is the lazy way... I 
      ! should fix it eventually.
      ! The same goes for the number of BC blocks.
      allocate(grid%iregbc1(MAXBCREGIONS))
      allocate(grid%iregbc2(MAXBCREGIONS))
      allocate(grid%k1bc(MAXBCBLOCKS))
      allocate(grid%k2bc(MAXBCBLOCKS))
      allocate(grid%j1bc(MAXBCBLOCKS))
      allocate(grid%j2bc(MAXBCBLOCKS))
      allocate(grid%i1bc(MAXBCBLOCKS))
      allocate(grid%i2bc(MAXBCBLOCKS))
      
	  
	  
	  
      allocate(grid%k1src(MAXSRC))
      allocate(grid%k2src(MAXSRC))
      allocate(grid%j1src(MAXSRC))
      allocate(grid%j2src(MAXSRC))
      allocate(grid%i1src(MAXSRC))
      allocate(grid%i2src(MAXSRC))
      allocate(grid%timesrc(MAXSRCTIMES,MAXSRC))
      allocate(grid%tempsrc(MAXSRCTIMES,MAXSRC))
      allocate(grid%qsrc(MAXSRCTIMES,MAXSRC, grid%nphase))
      
      
      allocate(grid%i1reg(MAXPERMREGIONS))
      allocate(grid%i2reg(MAXPERMREGIONS))
      allocate(grid%j1reg(MAXPERMREGIONS))
      allocate(grid%j2reg(MAXPERMREGIONS))
      allocate(grid%k1reg(MAXPERMREGIONS))
      allocate(grid%k2reg(MAXPERMREGIONS))
      allocate(grid%icap_reg(MAXPERMREGIONS))
      allocate(grid%ithrm_reg(MAXPERMREGIONS))
      allocate(grid%por_reg(MAXPERMREGIONS))
      allocate(grid%tor_reg(MAXPERMREGIONS))
      allocate(grid%perm_reg(MAXPERMREGIONS,4))
      
      allocate(grid%i1ini(MAXINITREGIONS))
      allocate(grid%i2ini(MAXINITREGIONS))
      allocate(grid%j1ini(MAXINITREGIONS))
      allocate(grid%j2ini(MAXINITREGIONS))
      allocate(grid%k1ini(MAXINITREGIONS))
      allocate(grid%k2ini(MAXINITREGIONS))
	  
    allocate(grid%xx_ini(grid%ndof,MAXINITREGIONS))
	    	  
      allocate(grid%i1brk(MAXINITREGIONS))
      allocate(grid%i2brk(MAXINITREGIONS))
      allocate(grid%j1brk(MAXINITREGIONS))
      allocate(grid%j2brk(MAXINITREGIONS))
      allocate(grid%k1brk(MAXINITREGIONS))
      allocate(grid%k2brk(MAXINITREGIONS))
      allocate(grid%ibrktyp(MAXINITREGIONS))
      allocate(grid%ibrkface(MAXINITREGIONS))
      
      allocate(grid%rock_density(MAXPERMREGIONS))
      allocate(grid%cpr(MAXPERMREGIONS))
      allocate(grid%dencpr(MAXPERMREGIONS))
      allocate(grid%ckdry(MAXPERMREGIONS))
      allocate(grid%ckwet(MAXPERMREGIONS))
      allocate(grid%tau(MAXPERMREGIONS))
      allocate(grid%cdiff(MAXPERMREGIONS))
      allocate(grid%cexp(MAXPERMREGIONS))

      allocate(grid%icaptype(MAXPERMREGIONS))
	  allocate(grid%sir(1:grid%nphase,MAXPERMREGIONS))
	  allocate(grid%lambda(MAXPERMREGIONS))
      allocate(grid%alpha(MAXPERMREGIONS))
      allocate(grid%pckrm(MAXPERMREGIONS))
      allocate(grid%pcwmax(MAXPERMREGIONS))
	  allocate(grid%pcbetac(MAXPERMREGIONS))
	  allocate(grid%pwrprm(MAXPERMREGIONS))
      
  !-----------------------------------------------------------------------
  ! Set up boundary condition storage on blocks
  !-----------------------------------------------------------------------
      allocate(grid%velocitybc0(grid%nphase, MAXBCREGIONS))
	  grid%velocitybc0 = 0.d0


	  allocate(grid%xxbc0(grid%ndof,MAXBCREGIONS))
	  grid%xxbc0=0.D0
	  	  
      
      !set scale factor for heat equation, i.e. use units of MJ for energy
	  grid%scale = 1.d-6
    print *,"new grid 5:",grid%nx,grid%ny,grid%nz

    end subroutine pflowGrid_new


!#include "pflowgrid_destroy.F90"

  subroutine pflowGrid_destroy(grid)
  
  implicit none
  
  type(pflowGrid), intent(inout) :: grid
  
  integer :: ierr

  ! Deallocate all of the arrays contained within grid.
  
    end subroutine pflowGrid_destroy

!======================================================================

    subroutine pflowGrid_newpm( parameters )
  
      implicit none

      type(pflowGridParameters) :: parameters
      type(pflowGrid), pointer :: grid
      type(time_stepping_context), pointer :: timestep
      
      interface 
         subroutine allocate_patch_info(p_samr_hierarchy, patchlevel_info)
           use pflow_gridtype_module
           PetscFortranAddr :: p_samr_hierarchy
           type(PatchLevelInfoPtr), dimension(:), pointer :: patchlevel_info
         end subroutine allocate_patch_info
      end interface

      integer :: igeom, nx, ny, nz, npx, npy, npz
      integer :: nphase
!      integer, intent(in) :: equ_option
      integer :: n, ng, na
      integer :: i, j, k
      integer :: ierr
      integer :: nlevels

      grid =>  parameters%grid
      timestep =>  parameters%timestep

      grid%igeom   = parameters%igeom
      grid%nx      = parameters%nx
      grid%ny      = parameters%ny
      grid%nz      = parameters%nz
      grid%nphase  = parameters%nphase
      grid%npx     = parameters%npx
      grid%npy     = parameters%npy
      grid%npz     = parameters%npz
      
      grid%Samrai_drive     = parameters%usesamrai
      grid%p_samr_hierarchy = parameters%p_samr_hierarchy

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz
      npx = grid%npx
      npy = grid%npy
      npz = grid%npz
      igeom = grid%igeom
      nphase = grid%nphase

	 ! integer, pointer::ghostind(:)

      call MPI_Comm_rank(PETSC_COMM_WORLD, grid%myrank, ierr)
      call MPI_Comm_size(PETSC_COMM_WORLD, grid%commsize, ierr)
       !call pflowGrid_parse_cmdline(grid)
      grid%use_ksp=PETSC_FALSE
      grid%use_isoth=PETSC_FALSE

      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-snes_mf", & 
      grid%use_matrix_free, ierr)
      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_analytical", &
      grid%use_analytical, ierr)
      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-print_hhistory", &
      grid%print_hhistory, ierr)
      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-monitor_h", &
      grid%monitor_h, ierr)
       call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_ksp", &
      grid%use_ksp, ierr)
      
      grid%nxy = nx*ny
      grid%nmax = nx*ny*nz
      print *,"new grid 1:",grid%nx,grid%ny,grid%nz

      grid%ndof = nphase

    
   
!-----------------------------------------------------------------------
      ! Initialize some counter variables.
!-----------------------------------------------------------------------
      grid%t = 0.d0
      grid%flowsteps = 0
      grid%newtcum = 0
      grid%icutcum = 0

      !-----------------------------------------------------------------------
      ! Initialize some parameters to sensible values.  These are parameters
      ! which should be set via the command line or the input file, but it
      ! seems good practice to set them to sensible values when a pflowGrid
      ! is created.
!-----------------------------------------------------------------------
      print *,"new grid 2:",grid%nx,grid%ny,grid%nz, grid%nphase
      
      grid%use_matrix_free = 1
      grid%newton_max = 16
            
      grid%dt = 1.d0
      grid%atol = PETSC_DEFAULT_DOUBLE_PRECISION
      grid%rtol = PETSC_DEFAULT_DOUBLE_PRECISION
      grid%stol = PETSC_DEFAULT_DOUBLE_PRECISION
      grid%maxit = PETSC_DEFAULT_INTEGER
      grid%maxf = PETSC_DEFAULT_INTEGER
      
      grid%ihydrostatic = 0
      !grid%conc0 = 1.d-6
      
   
      !physical constants
      grid%gravity = 9.8068d0    ! m/s^2
	  grid%tref=50D0
 !    grid%gravity = 0d0    ! m/s^2
     
      !-----------------------------------------------------------------------
      ! Generate the DA objects that will manage communication.
      !-----------------------------------------------------------------------
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
           nx,ny,nz,npx,npy,npz,1,1, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           grid%da_1_dof,ierr)
  
   
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
           nx,ny,nz,npx,npy,npz,3*grid%nphase,1, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           grid%da_3np_dof,ierr)
  
      call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
           nx,ny,nz,npx,npy,npz,grid%ndof,1, &
           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
           grid%da_ndof,ierr)
 
	   
  !    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
  !             nx,ny,nz,npx,npy,npz,(grid%ndof+1)*(2+4*grid%nphase),1, &
  !             PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
  !             grid%da_var_dof,ierr)

	 
       !print DA info for each processor
       call DAView(grid%da_ndof,PETSC_VIEWER_STDOUT_WORLD,ierr)

      !-----------------------------------------------------------------------
      ! Create the vectors with parallel layout corresponding to the DA's,
      ! and, for vectors that need to be ghosted, create the corresponding
      ! ghosted vectors.
      !-----------------------------------------------------------------------

      ! 1 degree of freedom
   
     call Pflow_allocate_Vec(grid) 
     print *,"new grid 3:",grid%nx,grid%ny,grid%nz       

     

   !-----------------------------------------------------------------------
      ! Set up information about corners of local domain.
!-----------------------------------------------------------------------
   
     print *,"new grid 4:",grid%nx,grid%ny,grid%nz   

     ! allocate space to store information about the patch objects
     call allocate_patch_info(grid%p_samr_hierarchy, grid%patchlevel_info)
     
      !     if (grid%using_pflowGrid == PETSC_TRUE) &
!     allocate(grid%vvl_loc(locpat%nconn*grid%nphase))

      ! I don't like having a fixed number of boundary condition regions.
      ! Memory for these arrays ought to allocated by parsing the input file
      ! to determine the number of regions.  This is the lazy way... I 
      ! should fix it eventually.
      ! The same goes for the number of BC blocks.
      allocate(grid%iregbc1(MAXBCREGIONS))
      allocate(grid%iregbc2(MAXBCREGIONS))
      allocate(grid%k1bc(MAXBCBLOCKS))
      allocate(grid%k2bc(MAXBCBLOCKS))
      allocate(grid%j1bc(MAXBCBLOCKS))
      allocate(grid%j2bc(MAXBCBLOCKS))
      allocate(grid%i1bc(MAXBCBLOCKS))
      allocate(grid%i2bc(MAXBCBLOCKS))
      
	  
	  
	  
      allocate(grid%k1src(MAXSRC))
      allocate(grid%k2src(MAXSRC))
      allocate(grid%j1src(MAXSRC))
      allocate(grid%j2src(MAXSRC))
      allocate(grid%i1src(MAXSRC))
      allocate(grid%i2src(MAXSRC))
      allocate(grid%timesrc(MAXSRCTIMES,MAXSRC))
      allocate(grid%tempsrc(MAXSRCTIMES,MAXSRC))
      allocate(grid%qsrc(MAXSRCTIMES,MAXSRC, grid%nphase))
      
      
      allocate(grid%i1reg(MAXPERMREGIONS))
      allocate(grid%i2reg(MAXPERMREGIONS))
      allocate(grid%j1reg(MAXPERMREGIONS))
      allocate(grid%j2reg(MAXPERMREGIONS))
      allocate(grid%k1reg(MAXPERMREGIONS))
      allocate(grid%k2reg(MAXPERMREGIONS))
      allocate(grid%icap_reg(MAXPERMREGIONS))
      allocate(grid%ithrm_reg(MAXPERMREGIONS))
      allocate(grid%por_reg(MAXPERMREGIONS))
      allocate(grid%tor_reg(MAXPERMREGIONS))
      allocate(grid%perm_reg(MAXPERMREGIONS,4))
      
      allocate(grid%i1ini(MAXINITREGIONS))
      allocate(grid%i2ini(MAXINITREGIONS))
      allocate(grid%j1ini(MAXINITREGIONS))
      allocate(grid%j2ini(MAXINITREGIONS))
      allocate(grid%k1ini(MAXINITREGIONS))
      allocate(grid%k2ini(MAXINITREGIONS))
	  
    allocate(grid%xx_ini(grid%ndof,MAXINITREGIONS))
	    	  
      allocate(grid%i1brk(MAXINITREGIONS))
      allocate(grid%i2brk(MAXINITREGIONS))
      allocate(grid%j1brk(MAXINITREGIONS))
      allocate(grid%j2brk(MAXINITREGIONS))
      allocate(grid%k1brk(MAXINITREGIONS))
      allocate(grid%k2brk(MAXINITREGIONS))
      allocate(grid%ibrktyp(MAXINITREGIONS))
      allocate(grid%ibrkface(MAXINITREGIONS))
      
      allocate(grid%rock_density(MAXPERMREGIONS))
      allocate(grid%cpr(MAXPERMREGIONS))
      allocate(grid%dencpr(MAXPERMREGIONS))
      allocate(grid%ckdry(MAXPERMREGIONS))
      allocate(grid%ckwet(MAXPERMREGIONS))
      allocate(grid%tau(MAXPERMREGIONS))
      allocate(grid%cdiff(MAXPERMREGIONS))
      allocate(grid%cexp(MAXPERMREGIONS))

      allocate(grid%icaptype(MAXPERMREGIONS))
      allocate(grid%sir(1:grid%nphase,MAXPERMREGIONS))
      allocate(grid%lambda(MAXPERMREGIONS))
      allocate(grid%alpha(MAXPERMREGIONS))
      allocate(grid%pckrm(MAXPERMREGIONS))
      allocate(grid%pcwmax(MAXPERMREGIONS))
      allocate(grid%pcbetac(MAXPERMREGIONS))
      allocate(grid%pwrprm(MAXPERMREGIONS))
      
  !-----------------------------------------------------------------------
  ! Set up boundary condition storage on blocks
  !-----------------------------------------------------------------------
      allocate(grid%velocitybc0(grid%nphase, MAXBCREGIONS))
	  grid%velocitybc0 = 0.d0


	  allocate(grid%xxbc0(grid%ndof,MAXBCREGIONS))
	  grid%xxbc0=0.D0
	  	  
      
      !set scale factor for heat equation, i.e. use units of MJ for energy
	  grid%scale = 1.d-6
    print *,"new grid 5:",grid%nx,grid%ny,grid%nz

    allocate(timestep%tfac(13))
    
    timestep%tfac(1)  = 2.0d0; timestep%tfac(2)  = 2.0d0
    timestep%tfac(3)  = 2.0d0; timestep%tfac(4)  = 2.0d0
    timestep%tfac(5)  = 2.0d0; timestep%tfac(6)  = 1.8d0
    timestep%tfac(7)  = 1.6d0; timestep%tfac(8)  = 1.4d0
    timestep%tfac(9)  = 1.2d0; timestep%tfac(10) = 1.0d0
    timestep%tfac(11) = 1.0d0; timestep%tfac(12) = 1.0d0
    timestep%tfac(13) = 1.0d0      
    timestep%iaccel = 1
    timestep%dt_min = 1.d0
    timestep%dt_max = 3.1536d6 ! One-tenth of a year.
    timestep%icut_max = 16
    timestep%dpmxe = 5.d4
    timestep%dsmxe = 5.d-1

    end subroutine pflowGrid_newpm

!======================================================================

!#include "pflowgrid_setup.F90"

! pflowGrid_setup():
! After the constructor function pflowGrid_new() has been used to 
! construct a new pflowGrid object, the grid's topology has been set 
! up but its geometry hasn't.  pflowGrid_setup() sets up the physical
! geometry of the grid, the initial values of the fields, and the 
! boundary conditions.

 subroutine pflowGrid_setup(grid, timestep, locpat, inputfile)
  use utilities_module
  use IMS_module
  use readfield                          
  implicit none
  
  type(pflowGrid), intent(inout) :: grid
  type(time_stepping_context), intent(inout) :: timestep
   type(pflow_localpatch_info), intent(inout) :: locpat
  
  
  character(len=*), intent(in) :: inputfile
  
  integer :: ierr
  integer :: i, j, jn1, jn2, k, ird
  integer :: mg1, mg2
  integer :: m, n, ng
  integer :: nc  ! Tracks number of connections computed.
  integer :: ibc ! Used to index boundary condition blocks.
  integer :: ir ! Used to index boundary condition regions.
  integer :: ii1, ii2, jj1, jj2, kk1, kk2
    ! Used for indexing corners of boundary regions.
  
    ! Used to access the contents of the local portion of grid%volume.
  
    ! Used for coloring the Jacobian so that we can take advantage of its
    ! sparsity when calculating it via finite differences.
  real*8, pointer :: ran_p(:),&
                     phis_p(:), icap_p(:),ithrm_p(:), por_p(:),por0_p(:),&
                     tor_p(:),perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),&
                     perm_pow_p(:)
  real*8 :: d1,d2,val,val1,val2,val3,val4,por,dw_kg
  real*8 :: dl,hl
  integer :: iseed,nx,ny,nz,na

  !,temp5_nat_vec,temp6_nat_vec,temp7_nat_vec
  
! external SNESDefaultComputeJacobian

#include "definitions.h"
  
!#ifdef DEBUG
! PetscViewer :: view_out
!#endif


  !-----------------------------------------------------------------------
  ! Parse the input file to get dx, dy, dz, fields, etc. for each cell. 
  !-----------------------------------------------------------------------
  call pflow_setup_index(grid,locpat)
  allocate(grid%ibndtyp(MAXBCREGIONS))
  allocate(grid%iface(MAXBCREGIONS))

  call pflowGrid_read_input(grid, timestep, inputfile)
  
 
! check number of dofs and phases


  if(grid%ndof .ne. (grid%nphase)) then
     write(*,*) 'Specified number of dofs or phases not correct-stop: MPH ', &
          'ndof= ',grid%ndof,' nph= ',grid%nphase
     stop
  endif
  	
  
! Calculate the x, y, z vectors that give the 
! physical coordinates of each cell.
  

  if(grid% Samrai_drive==PETSC_FALSE) then
     allocate(grid%x(grid%nmax))
     allocate(grid%y(grid%nmax))
     allocate(grid%z(grid%nmax))
     
     
     call pflowGrid_compute_xyz(grid)
     
     if (grid%myrank == 0) then
        write(*,'(/,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')
        write(*,'(" number of processors = ",i5,", npx,y,z= ",3i5)') &
             grid%commsize,grid%npx,grid%npy,grid%npz
        write(*,'(" number of dofs = ",i3,", number of phases = ",i3,i2)') &
             grid%ndof,grid%nphase
     endif
  endif


  !-----------------------------------------------------------------------
  ! Set up cell topology and geometry of interior connections, and 
  ! calculate interior interface areas and cell volumes.
  !-----------------------------------------------------------------------
  
  call pflowGrid_Setup_Geom(grid,locpat)

! set initial conditions by region for pressure, temperature, saturation
! and concentration

  call pflow_IMS_setupini(grid, locpat)
  
  ! set hydrostatic properties for initial and boundary conditions with depth
  
  if(grid%iread_init==2) call Read_init_field(grid, locpat)
  if(grid%nphase>=2) grid%velocitybc0(2,:) = grid%velocitybc0(1,:)
	  
  
 
! call VecView(grid%conc,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
!call VecView(grid%pressure,PETSC_VIEWER_STDOUT_WORLD,ierr)
!call VecView(grid%temp,PETSC_VIEWER_STDOUT_WORLD,ierr)
!call VecView(grid%xmol,PETSC_VIEWER_STDOUT_WORLD,ierr)
!call VecView(grid%sat,PETSC_VIEWER_STDOUT_WORLD,ierr)

  deallocate(grid%k1ini)
  deallocate(grid%k2ini)
  deallocate(grid%j1ini)
  deallocate(grid%j2ini)
  deallocate(grid%i1ini)
  deallocate(grid%i2ini)
  
  deallocate(grid%xx_ini)
 print *,'deallocate ini'

!************End of initial Condition Setup ***********************
 
  call pflowGrid_setup_BC(grid, locpat)

 if (grid%myrank == 0) &
  write(*,'("  Finished setting up of Geometry ")')


 call pflowGrid_Setup_Trans(grid,locpat)

 if (grid%myrank == 0) &
  write(*,'("  Finished setting up of INIT ")')


  !-----------------------------------------------------------------------
  ! Initialize field variables
  !-----------------------------------------------------------------------
      
   
     
    call pflow_update_ims(grid, locpat)
    print *, "IMS finish variable packing"
  
   
! call VecView(grid%xmol,PETSC_VIEWER_STDOUT_WORLD,ierr)
! call VecView(grid%yy,PETSC_VIEWER_STDOUT_WORLD,ierr)
! zero initial velocity
  call VecSet(grid%vl,0.d0,ierr)
   
  if (grid%myrank == 0) &
  write(*,'("  Finished setting up of INIT2 ")')

! call VecView(grid%yy,PETSC_VIEWER_STDOUT_WORLD,ierr)

  ! set phase index for each node and initialize accumulation terms
   
   call pflow_ims_initaccum(grid, locpat)
   
    
   
  if (grid%myrank == 0) &
  write(*,'("  Finished setting up ")')

  end subroutine pflowGrid_setup
!=============================================================================

! Subroutine to allocate global vector

 subroutine Pflow_allocate_Vec(grid)
 implicit none
 
 type(pflowGrid), intent(inout) :: grid
 PetscTruth :: use_ghost

 integer ierr
 
     call pflow_create_vector(grid, grid%da_1_dof, PETSC_FALSE, grid%porosity)
     call pflow_create_vector(grid, grid%da_1_dof, PETSC_TRUE, grid%dx_loc)
     call pflow_create_vector(grid, grid%da_3np_dof, PETSC_FALSE, grid%vl)
     call pflow_create_vector(grid, grid%da_ndof, PETSC_FALSE, grid%xx)
     call pflow_create_vector(grid, grid%da_ndof, PETSC_FALSE, grid%xx_loc)

! 1 degree of freedom      
      !global
      call VecDuplicate(grid%porosity, grid%porosity0, ierr)
      call VecDuplicate(grid%porosity, grid%tor, ierr)
      call VecDuplicate(grid%porosity, grid%dx, ierr)
      call VecDuplicate(grid%porosity, grid%dy, ierr)
      call VecDuplicate(grid%porosity, grid%dz, ierr)
      call VecDuplicate(grid%porosity, grid%volume, ierr)
      call VecDuplicate(grid%porosity, grid%ithrm, ierr)
      call VecDuplicate(grid%porosity, grid%icap, ierr)
      call VecDuplicate(grid%porosity, grid%temp, ierr)
	  call VecDuplicate(grid%porosity, grid%ttemp, ierr)
    
      call VecDuplicate(grid%porosity, grid%perm_xx, ierr)
      call VecDuplicate(grid%porosity, grid%perm_yy, ierr)
      call VecDuplicate(grid%porosity, grid%perm_zz, ierr)
	  call VecDuplicate(grid%porosity, grid%perm0_xx, ierr)
      call VecDuplicate(grid%porosity, grid%perm0_yy, ierr)
      call VecDuplicate(grid%porosity, grid%perm0_zz, ierr)
	  call VecDuplicate(grid%porosity, grid%perm_pow, ierr)
           
      call VecDuplicate(grid%dx_loc, grid%dy_loc, ierr)
      call VecDuplicate(grid%dx_loc, grid%dz_loc, ierr)
      call VecDuplicate(grid%dx_loc, grid%porosity_loc, ierr)
	  call VecDuplicate(grid%dx_loc, grid%tor_loc, ierr)
      call VecDuplicate(grid%dx_loc, grid%ithrm_loc, ierr)
      call VecDuplicate(grid%dx_loc, grid%icap_loc, ierr)
      
      call VecDuplicate(grid%dx_loc, grid%perm_xx_loc, ierr)
      call VecDuplicate(grid%dx_loc, grid%perm_yy_loc, ierr)
      call VecDuplicate(grid%dx_loc, grid%perm_zz_loc, ierr)

  
      ! ndof degrees of freedom
      call VecDuplicate(grid%xx, grid%yy, ierr)
      call VecDuplicate(grid%xx, grid%dxx, ierr)
      call VecDuplicate(grid%xx, grid%r, ierr)
      call VecDuplicate(grid%xx, grid%accum, ierr)
      
      call VecSetBlocksize(grid%dxx, grid%ndof, ierr)
   print *,"nVec 1:",grid%nx,grid%ny,grid%nz
     
 end subroutine Pflow_allocate_Vec
 
subroutine pflow_create_vector(grid, da_info, use_ghost, vec)
implicit none

  type(pflowGrid), intent(inout) :: grid
  DA :: da_info
  Vec :: vec
  PetscTruth :: use_ghost
  integer :: ierr
  integer :: dof

  if(grid%Samrai_drive==PETSC_TRUE) then
! in this case create a SAMRAI Vec
     call DAGetInfo(da_info, &
                    PETSC_NULL_INTEGER, &
                    PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                    PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                    dof, & 
                    PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                    ierr)
     call create_samrai_vec(grid%p_samr_hierarchy, dof, use_ghost, vec)
     
  else
! create a PETSc Vec
     if(use_ghost == PETSC_FALSE) then
        call DACreateGlobalVector(da_info, vec, ierr)
     else
        call DACreateLocalVector(da_info, vec, ierr)
     endif
  endif
  
end subroutine pflow_create_vector

subroutine grid_get_corners(grid, locpat, use_ghost)
implicit none

  type(pflowGrid), intent(inout) :: grid
  DA :: da_info
  type(pflow_localpatch_info), intent(inout) :: locpat
  PetscTruth, intent(in) :: use_ghost
  integer :: ierr

  if(grid%Samrai_drive==PETSC_FALSE) then
     if(use_ghost==PETSC_FALSE) then
        call DAGetCorners(grid%da_1_dof, locpat%nxs, &
                          locpat%nys, locpat%nzs, locpat%nlx, &
                          locpat%nly, locpat%nlz, ierr)
     else
        call DAGetGhostCorners(grid%da_1_dof, locpat%ngxs, &
                               locpat%ngys, locpat%ngzs, locpat%ngx, &
                               locpat%ngy, locpat%ngz, ierr)

     endif
  else
     if(use_ghost==PETSC_FALSE) then
        call samr_patch_get_corners(locpat%p_samr_patch, &
                               locpat%nxs, locpat%nys, locpat%nzs, &
                               locpat%nlx, locpat%nly, locpat%nlz)
     else
        call samr_patch_get_ghostcorners(locpat%p_samr_patch, &
                                    locpat%ngxs, locpat%ngys, locpat%ngzs, &
                                    locpat%ngx, locpat%ngy, locpat%ngz)
     endif
  endif

end subroutine grid_get_corners

function pims_patch_at_boundary(grid,locpat, axis, dim)
  implicit none
  integer :: pims_patch_at_boundary
  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info), intent(inout) :: locpat
  integer :: axis
  integer :: dim
  integer :: samr_patch_at_bc

  if(grid%Samrai_drive==PETSC_TRUE) then
     pims_patch_at_boundary  =  samr_patch_at_bc(locpat%p_samr_patch, axis, dim)
  else
     pims_patch_at_boundary = 1
  endif

end function pims_patch_at_boundary

!Subroutine to setup index system

subroutine pflow_setup_index(grid, locpat)
implicit none

type(pflowGrid), intent(inout) :: grid
type(pflow_localpatch_info), intent(inout) :: locpat

   integer i,j,k,n,ng, na, ierr
   call grid_get_corners(grid, locpat, PETSC_FALSE)

!   call DAGetCorners(grid%da_1_dof, locpat%nxs, &
!        locpat%nys, locpat%nzs, locpat%nlx, &
!        locpat%nly, locpat%nlz, ierr)

   locpat%nxe = locpat%nxs + locpat%nlx
   locpat%nye = locpat%nys + locpat%nly
   locpat%nze = locpat%nzs + locpat%nlz
   locpat%nlxy = locpat%nlx * locpat%nly
   locpat%nlxz = locpat%nlx * locpat%nlz
   locpat%nlyz = locpat%nly * locpat%nlz
   locpat%nlmax = locpat%nlx * locpat%nly * locpat%nlz
   locpat%nldof = locpat%nlmax * grid%nphase
   
   call grid_get_corners(grid, locpat, PETSC_TRUE)

!   call DAGetGhostCorners(grid%da_1_dof, locpat%ngxs, &
!        locpat%ngys, locpat%ngzs, locpat%ngx, &
!        locpat%ngy, locpat%ngz, ierr)
   
   locpat%ngxe = locpat%ngxs + locpat%ngx
   locpat%ngye = locpat%ngys + locpat%ngy
   locpat%ngze = locpat%ngzs + locpat%ngz
   locpat%ngxy = locpat%ngx * locpat%ngy
   locpat%ngxz = locpat%ngx * locpat%ngz
   locpat%ngyz = locpat%ngy * locpat%ngz
   locpat%ngmax = locpat%ngx * locpat%ngy * locpat%ngz
   locpat%ngdof = locpat%ngmax * grid%nphase

!-----------------------------------------------------------------------
      ! Determine number of local connections.
!-----------------------------------------------------------------------
      locpat%nconn = (locpat%ngx-1) * locpat%nly &
        * locpat%nlz + locpat%nlx * (locpat%ngy-1) &
        * locpat%nlz + locpat%nlx * locpat%nly &
        * (locpat%ngz-1)

  if(grid%samrai_drive == PETSC_TRUE) then
     if(pims_patch_at_boundary(grid,locpat, 0, 0) ==1 ) then !x, left
        locpat%nconn = locpat%nconn - locpat%nlyz
     endif  
     if(pims_patch_at_boundary(grid,locpat, 0, 1) ==1) then 
        locpat%nconn = locpat%nconn - locpat%nlyz
     endif  
     if(pims_patch_at_boundary(grid,locpat, 1, 0) ==1) then
        locpat%nconn = locpat%nconn - locpat%nlxz
     endif  
     if(pims_patch_at_boundary(grid,locpat, 1, 1) ==1) then
        locpat%nconn = locpat%nconn - locpat%nlxz
     endif  
     if(pims_patch_at_boundary(grid,locpat, 2, 0) ==1) then
        locpat%nconn = locpat%nconn - locpat%nlxy
     endif  
     if(pims_patch_at_boundary(grid,locpat, 2, 1) ==1) then
        locpat%nconn = locpat%nconn - locpat%nlxy
     endif  
   endif

!-----------------------------------------------------------------------
      ! Allocate memory for allocatable arrays.
!-----------------------------------------------------------------------
      allocate(locpat%nL2G(locpat%nlmax))
      allocate(locpat%nG2L(locpat%ngmax))

      if(grid%Samrai_drive==PETSC_FALSE) then
         allocate(locpat%nL2A(locpat%nlmax))
         allocate(locpat%nG2N(locpat%ngmax))
      endif

 !-----------------------------------------------------------------------
      ! Compute arrays for indexing between local ghosted and non-ghosted 
      ! arrays.  I think the PETSc DA facilities may make these redundant,
      ! but since the ptran code uses them, I think it is easiest to just
      ! use these.
!-----------------------------------------------------------------------
      ! Local <-> Ghosted Transformation (Two ways)
      locpat%nG2L(:) = 0  ! Must initialize this to zero!
      !grid%nL2N(:) = 0
      n = 0
      locpat%istart = locpat%nxs-locpat%ngxs
      locpat%jstart = locpat%nys-locpat%ngys
      locpat%kstart = locpat%nzs-locpat%ngzs
      locpat%iend = locpat%istart+locpat%nlx-1
      locpat%jend = locpat%jstart+locpat%nly-1
      locpat%kend = locpat%kstart+locpat%nlz-1

      do k=locpat%kstart,locpat%kend
        do j=locpat%jstart,locpat%jend
          do i=locpat%istart,locpat%iend
            n = n + 1
            ng = i+j*locpat%ngx+k*locpat%ngxy+1
			locpat%nL2G(n) = ng
            locpat%nG2L(ng) = n
		  enddo
        enddo
      enddo

      do i=1,locpat%ngmax
        j = locpat%nG2L(i)
        if (j > 0) then
          k = locpat%nL2G(j)
          if(i /= k) then
            print *,'Error in ghost-local node numbering for ghost node =', i
            print *,'node_id_gtol(i) =', j
            print *,'node_id_ltog(node_id_gtol(i)) =', k
            stop
          endif
        endif
      enddo
	  
      if(grid%Samrai_drive==PETSC_FALSE) then
         ! Local(non ghosted)->Natural(natural order starts from 0) (one way)
         n=0
         do k=1,locpat%nlz
            do j=1,locpat%nly
               do i=1,locpat%nlx
                  n = n + 1
                  na = i-1+locpat%nxs+(j-1+locpat%nys)*grid%nx+(k-1+locpat%nzs)*grid%nxy
                  if(na>(grid%nmax-1)) print *,'Wrong Nature order....'
                  locpat%nL2A(n) = na
                  !print *,grid%myrank, k,j,i,n,na
                  !grid%nG2N(ng) = na
               enddo
            enddo
         enddo
         print *,grid%myrank, locpat%nxs,locpat%ngxs,locpat%nys,&
              locpat%ngys,locpat%nzs,locpat%ngzs
         
         ! Loacal(include ghosted) to global mapping (one way)
         call DAGetGlobalIndicesF90(grid%da_1_dof,locpat%ngmax,locpat%nG2N, ierr)
      endif

end subroutine pflow_setup_index												             
      	  


! Subroutine to setup geomtry and connections

subroutine pflowGrid_setup_Geom(grid, locpat)
use readfield
implicit none
type(pflowGrid) grid
type(pflow_localpatch_info) :: locpat

interface
   subroutine pims_vecgetarrayf90(grid, patch, vec, f90ptr, ierr)
     use pflow_gridtype_module
     implicit none
#include "include/finclude/petsc.h"

     type(pflowGrid), intent(inout) :: grid
     type(pflow_localpatch_info) :: patch
     Vec :: vec
     PetscScalar, dimension(:), pointer :: f90ptr
     integer :: ierr
   end subroutine pims_vecgetarrayf90
end interface
integer i,j,k, nc, mg1, mg2, n, ierr, ng, leng
real*8 val, d1,d2
PetscScalar, pointer ::  dx_loc_p(:), dy_loc_p(:), dz_loc_p(:), volume_p(:)
PetscScalar, pointer ::  dx_p(:), dy_p(:), dz_p(:)
! Vec :: temp0_nat_vec, temp1_nat_vec, temp2_nat_vec, temp3_nat_vec, &
!          temp4_nat_vec 

 allocate(locpat%nd1(locpat%nconn))
 allocate(locpat%nd2(locpat%nconn))
 allocate(locpat%dist1(locpat%nconn))
 allocate(locpat%dist2(locpat%nconn))
 allocate(locpat%area(locpat%nconn))
 allocate(locpat%delz(locpat%nconn))
 
 allocate(locpat%iperm1(locpat%nconn))
 allocate(locpat%iperm2(locpat%nconn))


! old version **********************
#if 0
  call DACreateNaturalVector(grid%da_1_dof,temp1_nat_vec,ierr)
  call VecDuplicate(temp1_nat_vec, temp2_nat_vec, ierr)
  call VecDuplicate(temp1_nat_vec, temp3_nat_vec, ierr)
  if (grid%myrank == 0) then
!---set dx, dy, dz
    do k = 1,grid%nz
      do j = 1,grid%ny
        do i = 1,grid%nx
          n = i+(j-1)*grid%nx+(k-1)*grid%nxy-1
          val = grid%dx0(i)
          call VecSetValue(temp1_nat_vec,n,val,INSERT_VALUES,ierr)
          val = grid%dy0(j)
          call VecSetValue(temp2_nat_vec,n,val,INSERT_VALUES,ierr)
          val = grid%dz0(k)
          call VecSetValue(temp3_nat_vec,n,val,INSERT_VALUES,ierr)
        enddo
      enddo
    enddo
  endif



  call VecAssemblyBegin(temp1_nat_vec,ierr)
  call VecAssemblyEnd(temp1_nat_vec,ierr)
  call VecAssemblyBegin(temp2_nat_vec,ierr)
  call VecAssemblyEnd(temp2_nat_vec,ierr)
  call VecAssemblyBegin(temp3_nat_vec,ierr)
  call VecAssemblyEnd(temp3_nat_vec,ierr)
  
!  print *,'Pflow setup Geom :: Na Vec end assembly'
  
  call DANaturalToGlobalBegin(grid%da_1_dof,temp1_nat_vec,INSERT_VALUES, &
                                grid%dx,ierr)
  call DANaturalToGlobalEnd(grid%da_1_dof,temp1_nat_vec,INSERT_VALUES, &
                              grid%dx,ierr)
  call DANaturalToGlobalBegin(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
                                grid%dy,ierr)
  call DANaturalToGlobalEnd(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
                              grid%dy,ierr)
  call DANaturalToGlobalBegin(grid%da_1_dof,temp3_nat_vec,INSERT_VALUES, &
                                grid%dz,ierr)
  call DANaturalToGlobalEnd(grid%da_1_dof,temp3_nat_vec,INSERT_VALUES, &
                              grid%dz,ierr)
  
  call VecDestroy(temp1_nat_vec,ierr)
  call VecDestroy(temp2_nat_vec,ierr)
  call VecDestroy(temp3_nat_vec,ierr)
#endif
!old version end

!#if 0
  call pims_vecgetarrayf90(grid, locpat, grid%dx, dx_p, ierr)
  call pims_vecgetarrayf90(grid, locpat, grid%dy, dy_p, ierr)
  call pims_vecgetarrayf90(grid, locpat, grid%dz, dz_p, ierr)
   


      n=0
      do k=1,locpat%nlz
            do j=1,locpat%nly
               do i=1,locpat%nlx
                  n = n + 1
                  dx_p(n) = grid%dx0(i + locpat%nxs)
                  dy_p(n)=  grid%dy0(j + locpat%nys)
                  dz_p(n)=  grid%dz0(k + locpat%nzs)
                enddo
            enddo      
      enddo 
      
  call VecRestoreArrayF90( grid%dx, dx_p, ierr)   
  call VecRestoreArrayF90( grid%dy, dy_p, ierr)
  call VecRestoreArrayF90( grid%dz, dz_p, ierr)
!  call VecView(grid%dx,PETSC_VIEWER_STDOUT_WORLD,ierr)

!#endif
  
  
  
  ! Extract local, ghosted portions of dx, dy, dz vectors.
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%dx, INSERT_VALUES, &
                            grid%dx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%dx, INSERT_VALUES, &
                          grid%dx_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%dy, INSERT_VALUES, &
                            grid%dy_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%dy, INSERT_VALUES, &
                          grid%dy_loc,ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%dz, INSERT_VALUES, &
                            grid%dz_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%dz, INSERT_VALUES, &
                          grid%dz_loc,ierr)

  call pims_vecgetarrayf90(grid, locpat, grid%dx_loc, dx_loc_p, ierr)
  call pims_vecgetarrayf90(grid, locpat, grid%dy_loc, dy_loc_p, ierr)
  call pims_vecgetarrayf90(grid, locpat, grid%dz_loc, dz_loc_p, ierr)

  nc = 0
  locpat%nconnx = 0
  locpat%nconny = 0

  if(grid%igeom ==2) then
    allocate(grid%rd(0:grid%nx))
    grid%rd=0.D0
	grid%rd(0)=grid%Radius_0 
    do i=1, grid%nx
      grid%rd(i)=grid%rd(i-1) + grid%dx0(i)
	enddo
  endif	  
  
! print *,'setup-RD:: ',grid%myrank,grid%rd, grid%dx0
  
  ! x-connections
  if(locpat%nlx > 1) then
   
    do k = locpat%kstart, locpat%kend
      do j = locpat%jstart, locpat%jend
         leng = locpat%ngx - 1
         if(grid%Samrai_drive==PETSC_TRUE) then
           if(pims_patch_at_boundary(grid,locpat, 0, 1) ==1) leng = leng - 1
         endif
        do i = 1, leng
          mg1 = i + j * locpat%ngx + k * locpat%ngxy
          if(grid%Samrai_drive == PETSC_TRUE) mg1 = mg1 + 1 
          mg2 = mg1 + 1
          nc = nc + 1
          locpat%nd1(nc) = mg1
          locpat%nd2(nc) = mg2
          locpat%dist1(nc) = 0.5d0 * dx_loc_p(mg1)
          locpat%dist2(nc) = 0.5d0 * dx_loc_p(mg2)
          locpat%delz(nc) = 0.d0
        !  abs_x = rd(i-1+grid%nxs)
          if (grid%igeom == 1) then
            locpat%area(nc) = dy_loc_p(mg1) * dz_loc_p(mg1)
          else if (grid%igeom == 2) then
            locpat%area(nc) = 2.D0 * Pi * grid%rd(i+locpat%nxs) * dz_loc_p(mg1)
            print *,grid%myrank,nc,i,grid%rd(i+locpat%nxs)
          else if (grid%igeom == 3) then
          endif
          locpat%iperm1(nc) = 1
          locpat%iperm2(nc) = 1
        enddo
      enddo
    enddo





    locpat%nconnx = nc
  endif

  ! y-connections
  if(locpat%ngy > 1) then
    do k = locpat%kstart, locpat%kend
      do i = locpat%istart, locpat%iend
         leng = locpat%ngy - 1
         if(grid%Samrai_drive==PETSC_TRUE) then
           if(pims_patch_at_boundary(grid,locpat, 1, 1) ==1) leng = leng - 1
         endif
        do j = 1, leng
          mg1 = i + 1 + (j-1) * locpat%ngx + k * locpat%ngxy
          if(grid%Samrai_drive == PETSC_TRUE) mg1 = mg1 +  locpat%ngx
          mg2 = mg1 + locpat%ngx
          nc = nc + 1
          locpat%nd1(nc) = mg1
          locpat%nd2(nc) = mg2
          locpat%dist1(nc) = 0.5d0 * dy_loc_p(mg1)
          locpat%dist2(nc) = 0.5d0 * dy_loc_p(mg2)
          locpat%delz(nc) = 0.d0
          locpat%area(nc) = dx_loc_p(mg1) * dz_loc_p(mg1)
          locpat%iperm1(nc) = 2
          locpat%iperm2(nc) = 2
        enddo
      enddo
    enddo
    locpat%nconny = nc
  endif
      
  ! z-connections
  if(locpat%ngz > 1) then
    do j = locpat%jstart, locpat%jend
      do i = locpat%istart, locpat%iend
         leng = locpat%ngz - 1
         if(grid%Samrai_drive == PETSC_TRUE) then
           if(pims_patch_at_boundary(grid,locpat, 2, 1) ==1) leng = leng - 1
         endif
        do k = 1, leng
          mg1 = i + 1 + j * locpat%ngx + (k-1) * locpat%ngxy
          if(grid%Samrai_drive == PETSC_TRUE) mg1 = mg1 + locpat%ngxy
          mg2 = mg1 + locpat%ngxy
          nc = nc + 1
          locpat%nd1(nc) = mg1
          locpat%nd2(nc) = mg2
          d1 = 0.5d0 * dz_loc_p(mg1)
          d2 = 0.5d0 * dz_loc_p(mg2)
          locpat%dist1(nc) = d1
          locpat%dist2(nc) = d2
          locpat%delz(nc) = d1 + d2
          if (grid%igeom == 1) then
          locpat%area(nc) = dx_loc_p(mg1) * dy_loc_p(mg1)
          else if (grid%igeom == 2) then
            locpat%area(nc) =  Pi * (grid%rd(i+locpat%nxs)+ grid%rd(i-1+locpat%nxs))* &
            (grid%rd(i+locpat%nxs) - grid%rd(i-1+locpat%nxs))	
            print *, 'area nc ',grid%myrank, nc, i,  locpat%area(nc)
          else if (grid%igeom == 3) then
          endif
          locpat%iperm1(nc) = 3
          locpat%iperm2(nc) = 3
        enddo
      enddo
    enddo
  endif
  
  ! Calculate cell volumes for local cells.
  call pims_vecgetarrayf90(grid, locpat, grid%volume, volume_p, ierr)

  do n=1, locpat%nlmax
    ng = locpat%nL2G(n)
	if (grid%igeom == 1) then
      volume_p(n) = dx_loc_p(ng) * dy_loc_p(ng) * dz_loc_p(ng)
      print *, 'setup_geom: Vol ', grid%myrank, n,volume_p(n), dx_loc_p(ng),dy_loc_p(ng),dz_loc_p(ng)
    else if (grid%igeom == 2) then
      i= mod(mod((n),locpat%nlxy),locpat%nlx)!+(grid%ngxs-grid%nxs)
      if(i==0) i=  locpat%nlx
      volume_p(n) = Pi * (grid%rd(i+locpat%nxs)+ grid%rd(i-1+locpat%nxs))*&
      (grid%rd(i+locpat%nxs) - grid%rd(i-1+locpat%nxs)) * dz_loc_p(ng)
	  print *, 'setup: Vol ', grid%myrank, n,i, grid%rd(i+locpat%nxs),volume_p(n)
    endif
  enddo

  if(grid%Samrai_drive==PETSC_FALSE) then
     call VecRestoreArrayF90(grid%volume, volume_p, ierr)
     call VecRestoreArrayF90(grid%dx_loc, dx_loc_p, ierr)
     call VecRestoreArrayF90(grid%dy_loc, dy_loc_p, ierr)
     call VecRestoreArrayF90(grid%dz_loc, dz_loc_p, ierr)
  endif
  
  
  write(*,'(" myrank= ",i3,", nlmax= ",i6,", nlx,y,z= ",3i4, &
  & ", nxs,e = ",2i4,", nys,e = ",2i4,", nzs,e = ",2i4)') &
  grid%myrank,locpat%nlmax,locpat%nlx,locpat%nly,locpat%nlz, &
  locpat%nxs,locpat%nxe,locpat%nys,locpat%nye,&
  locpat%nzs,locpat%nze

  write(*,'(" myrank= ",i3,", ngmax= ",i6,", ngx,y,z= ",3i4, &
  & ", ngxs,e= ",2i4,", ngys,e= ",2i4,", ngzs,e= ",2i4)') &
  grid%myrank,locpat%ngmax,locpat%ngx,locpat%ngy,locpat%ngz, &
  locpat%ngxs,locpat%ngxe,locpat%ngys,locpat%ngye,&
  locpat%ngzs,locpat%ngze

end subroutine pflowGrid_Setup_Geom
!=========================================================================

! 
subroutine pflowGrid_Setup_Trans(grid, locpat)
  !-----------------------------------------------------------------------
  ! Set up the transformation from physical coordinates
  ! to local domains.
  !-----------------------------------------------------------------------
  use utilities_module
  use readfield
  implicit none

 type(pflowGrid), intent(inout) :: grid
 type(pflow_localpatch_info) :: locpat
 
integer :: iseed, n, nla, nx, ny, nz, ir, ierr
real*8, pointer :: icap_p(:),por_p(:),por0_p(:), ithrm_p(:), &
                   tor_p(:),perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),&
                   perm_pow_p(:), ran_p(:)
 real*8 val, val1,val2,val3,val4,por, random_nr, frand
 Vec :: temp0_nat_vec


! old version
#if 0
  call DACreateNaturalVector(grid%da_1_dof,temp0_nat_vec,ierr)
!  call VecDuplicate(temp0_nat_vec, temp1_nat_vec, ierr)
!  call VecDuplicate(temp0_nat_vec, temp2_nat_vec, ierr)
!  call VecDuplicate(temp0_nat_vec, temp3_nat_vec, ierr)
!  call VecDuplicate(temp0_nat_vec, temp4_nat_vec, ierr)
!  call VecDuplicate(temp0_nat_vec, temp5_nat_vec, ierr)
!  call VecDuplicate(temp0_nat_vec, temp6_nat_vec, ierr)
!  call VecDuplicate(temp0_nat_vec, temp7_nat_vec, ierr)

! call DACreateNaturalVector(grid%da_3_dof,temp3_nat_vec,ierr)

! set capillary index, thermal index, porosity and permeability by region
  
  if (grid%myrank == 0) then 

    iseed = 345678912
    do n=1,grid%nmax
      val=1.D0
      if (grid%ran_fac > 0.d0) random_nr = ran1(n)
      val = random_nr
      call VecSetValue(temp0_nat_vec,n-1,val,INSERT_VALUES,ierr)
     enddo
   endif
     call VecAssemblyBegin(temp0_nat_vec,ierr)
     call VecAssemblyEnd(temp0_nat_vec,ierr)
     call DANaturalToGlobalBegin(grid%da_1_dof,temp0_nat_vec,INSERT_VALUES, &
                                grid%ttemp,ierr)
     call DANaturalToGlobalEnd(grid%da_1_dof,temp0_nat_vec,INSERT_VALUES, &
                              grid%ttemp,ierr)
#endif
! old version end

     call VecGetArrayF90(grid%ttemp,ran_p,ierr)
 ! call pims_vecgetarrayf90(grid, locpat, grid%ttemp, ran_p, ierr)
 ! Clu change, Bobby please note the  pims_vecgetarrayf90 not working for me
!#if 0
    do n = 1,locpat%nlmax
       random_nr=1.D0
       if (grid%ran_fac > 0.d0) random_nr = ran1(n)
       val = random_nr
       ran_p(n)=val
    enddo   
  
!#endif
     print *, 'end random'
     call VecGetArrayF90(grid%icap,icap_p,ierr)
     call VecGetArrayF90(grid%ithrm,ithrm_p,ierr)
     call VecGetArrayF90(grid%porosity,por_p,ierr)
     call VecGetArrayF90(grid%porosity0,por0_p,ierr)
     call VecGetArrayF90(grid%perm_xx,perm_xx_p,ierr)
     call VecGetArrayF90(grid%perm_yy,perm_yy_p,ierr)
     call VecGetArrayF90(grid%perm_zz,perm_zz_p,ierr)
     call VecGetArrayF90(grid%perm_pow,perm_pow_p,ierr)
     call VecGetArrayF90(grid%tor,tor_p,ierr)
     do n = 1,locpat%nlmax
          !na = locpat%nL2A(n)
          !nz= int(na/grid%nxy) + 1
          !ny= int(mod(na,grid%nxy)/grid%nx) + 1
          !nx= mod(mod(na,grid%nxy),grid%nx) + 1
          nla=n-1
          nz = int(nla/locpat%nlxy) + 1 + locpat%nzs
          ny = int(mod(nla,locpat%nlxy)/locpat%nlx) + 1 + locpat%nys
          nx = mod(mod(nla,locpat%nlxy),locpat%nlx) + 1 + locpat%nxs


          do ir = 1,grid%iregperm        
          if ((nz>=grid%k1reg(ir)) .and. (nz<=grid%k2reg(ir)) .and.&
           (ny>=grid%j1reg(ir)) .and. (ny<=grid%j2reg(ir)) .and.&
           (nx>= grid%i1reg(ir)) .and. (nx<=grid%i2reg(ir)))then
                                
            val = grid%icap_reg(ir)
           ! call VecSetValue(temp0_nat_vec,n,val,INSERT_VALUES,ierr)
            icap_p(n)=val
            
            val = grid%ithrm_reg(ir)
           ! call VecSetValue(temp1_nat_vec,n,val,INSERT_VALUES,ierr)
           ithrm_p(n)=val
           
            random_nr=1.D0
            if(grid%ran_fac > 0.d0) then
!             frand = rand(0)
              frand = ran1(n)
!             frand = ran_p(n)
!             random_nr = 1.d0+2.d0*grid%ran_fac*(frand-0.5D0)
              random_nr = grid%ran_fac*frand+1.d-6
              
!             print *,'pflowgrid_mod: ',n,frand,random_nr,grid%ran_fac
            endif

            por = grid%por_reg(ir)
            por0_p(n)=por
            if(grid%iran_por==1) then
              por=por*(2.D0**0.666667D0*(frand)**1.5D0)
              if(por<1D-2) por=1D-2 
			  
            endif
           ! call VecSetValue(temp2_nat_vec,n,por,INSERT_VALUES,ierr)
            por_p(n)= por
            
!           nn = 3*n
!           val1 = grid%perm_reg(ir,1)
!           call VecSetValue(temp3_nat_vec,nn,val1,INSERT_VALUES,ierr)
            
!           val2 = grid%perm_reg(ir,2)
!           call VecSetValue(temp3_nat_vec,nn+1,val2,INSERT_VALUES,ierr)
            
!           val3 = grid%perm_reg(ir,3)
!           call VecSetValue(temp3_nat_vec,nn+2,val3,INSERT_VALUES,ierr)
            val1 = grid%perm_reg(ir,1)
            val1 = val1*random_nr
            !call VecSetValue(temp3_nat_vec,n,val1,INSERT_VALUES,ierr)
            perm_xx_p(n)= val1
            
            val2 = grid%perm_reg(ir,2)
            val2 = val2*random_nr
           ! call VecSetValue(temp4_nat_vec,n,val2,INSERT_VALUES,ierr)
            perm_yy_p(n)= val2
              
            val3 = grid%perm_reg(ir,3)
            val3 = val3*random_nr
            !call VecSetValue(temp5_nat_vec,n,val3,INSERT_VALUES,ierr)
            perm_zz_p(n)= val3
            
            val4 = grid%perm_reg(ir,4) ! permpower
		  !	call VecSetValue(temp6_nat_vec,n,val4,INSERT_VALUES,ierr)
            perm_pow_p(n)=val4

	 		val3 = grid%tor_reg(ir)
	 		!call VecSetValue(temp7_nat_vec,n,val3,INSERT_VALUES,ierr)
            tor_p(n)=val3
            
   !       print *,'setup trans: ',ir, n,nx,ny,nz,locpat%nxs,perm_xx_p(n),   tor_p(n)
   !         exit
          endif
          
          enddo
      enddo
      
 
     call VecRestoreArrayF90(grid%ttemp,ran_p,ierr)
     call VecRestoreArrayF90(grid%icap,icap_p,ierr)
     call VecRestoreArrayF90(grid%ithrm,ithrm_p,ierr)
     call VecRestoreArrayF90(grid%porosity,por_p,ierr)
     call VecRestoreArrayF90(grid%porosity0,por0_p,ierr)
     call VecRestoreArrayF90(grid%perm_xx,perm_xx_p,ierr)
     call VecRestoreArrayF90(grid%perm_yy,perm_yy_p,ierr)
     call VecRestoreArrayF90(grid%perm_zz,perm_zz_p,ierr)
     call VecRestoreArrayF90(grid%perm_pow,perm_pow_p,ierr)
     call VecRestoreArrayF90(grid%tor,tor_p,ierr)

!  call VecDestroy(temp0_nat_vec,ierr)
 ! call VecDestroy(temp1_nat_vec,ierr)
 ! call VecDestroy(temp2_nat_vec,ierr)
 ! call VecDestroy(temp3_nat_vec,ierr)
 ! call VecDestroy(temp4_nat_vec,ierr)
 ! call VecDestroy(temp5_nat_vec,ierr)
 ! call VecDestroy(temp6_nat_vec,ierr)
 ! call VecDestroy(temp7_nat_vec,ierr)
#if 0
  call VecSet(grid%icap,1.d0,ierr)
  call VecSet(grid%ithrm,1.d0,ierr)
  call VecSet(grid%porosity,0.5,ierr)
! call VecSet(grid%perm,1.d-15,ierr)
#endif

 if(grid%iread_perm == 1)then
   call Read_perm_field(grid, locpat)
 endif  

 !if(grid%using_pflowGrid==0) call VecCopy(grid%Porosity, grid%Porosity0, ierr)
 call VecCopy(grid%perm_xx, grid%perm0_xx, ierr) 
 call VecCopy(grid%perm_yy, grid%perm0_yy, ierr) 
 call VecCopy(grid%perm_zz, grid%perm0_zz, ierr) 
 

end subroutine pflowGrid_Setup_Trans
!==================================================================================

!======================
subroutine pflowgrid_Setup_SNES(grid, pflowsolv)
use IMS_module
use pflow_solv_module
use pflow_step

implicit none
type(pflowGrid) :: grid
type(pflow_solver_context) :: pflowsolv
external SNESDefaultComputeJacobianColor
ISColoring :: iscoloring
integer myrank, ierr, maxstep
real*8 alapha, steptol


 myrank=grid%myrank

!-----------------------------------------------------------------------
      ! Set up PETSc nonlinear solver context.
!-----------------------------------------------------------------------
      call SNESCreate(PETSC_COMM_WORLD, pflowsolv%snes, ierr)
      CHKERRQ(ierr)


 if(grid%use_analytical == PETSC_TRUE) then
  
    grid%ideriv = 1
  
   if (myrank == 0) write(*,'(" analytical jacobian as ")');print *, grid%iblkfmt

    ! call DAGetMatrix(grid%da_ndof, MATMPIAIJ, grid%J, ierr)
    ! PETSc 2.1.6 introduces the MATAIJ matrix type.
    if (grid%iblkfmt == 0) then
      call DAGetMatrix(grid%da_ndof, MATAIJ, grid%J, ierr)
    else
      call DAGetMatrix(grid%da_ndof, MATBAIJ, grid%J, ierr)
   !   call  MatSetBlocksize(grid%J,grid%ndof,ierr)
    endif
   ! 
    call MatSetOption(grid%J,MAT_COLUMN_ORIENTED,ierr)
    
    call SNESSetJacobian(pflowsolv%snes, grid%J, grid%J, IMSJacobin, &
                         grid, ierr); CHKERRQ(ierr)
    if(grid%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid,pflowsolv)
	
  else if(grid%use_matrix_free == PETSC_TRUE) then
  
    grid%ideriv = 0
  
    if (myrank == 0) write(*,'(" Using matrix-free Newton-Krylov")')
    
    call MatCreateSNESMF(pflowsolv%snes, grid%xx, grid%J, ierr)

    
    ! It seems that I ought to call SNESSetJacobian here now, but I don't know
    ! what function I am supposed to pass to it to get it to use one of the 
    ! finite-difference routines for computing the Jacobian.  It might not
    ! actually matter if -snes_mf has been specified.
    ! Pernice thinks that perhaps the I need to provide a function which 
    ! simply calls MatAssemblyBegin/End.
    call SNESSetJacobian(pflowsolv%snes, grid%J, grid%J, &
                         ComputeMFJacobian, 0, ierr)

    ! Use "Walker-Pernice" differencing.
#if ((PETSC_VERSION_RELEASE==1)&&(PETSC_VERSION_MAJOR==2)&&(PETSC_VERSION_MINOR==3)&&(PETSC_VERSION_SUBMINOR==2)&&(PETSC_VERSION_PATCH<=10))
    call MatSNESMFSetType(grid%J, MATSNESMF_WP, ierr)
    if(grid%print_hhistory == PETSC_TRUE) then
      allocate(grid%hhistory(HHISTORY_LENGTH))
      call MatSNESMFSetHHistory(grid%J, grid%hhistory, HHISTORY_LENGTH, ierr)
   endif
#else
    call MatMFFDSetType(grid%J, MATMFFD_WP, ierr)
    if(grid%print_hhistory == PETSC_TRUE) then
      allocate(grid%hhistory(HHISTORY_LENGTH))
      call MatMFFDSetHHistory(grid%J, grid%hhistory, HHISTORY_LENGTH, ierr)
    endif
#endif

    if(grid%monitor_h == PETSC_TRUE) then
#if ((PETSC_VERSION_RELEASE)&&(PETSC_VERSION_MAJOR==2)&&(PETSC_VERSION_MINOR==3)&&(PETSC_VERSION_SUBMINOR==2)&&(PETSC_VERSION_PATCH<=10))
      call SNESSetMonitor(pflowsolv%snes, pflowgrid_MonitorH, grid, &
                          0, ierr)
#else
      call SNESMonitorSet(pflowsolv%snes, pflowgrid_MonitorH, grid, &
                          0    , ierr)
#endif
    endif
    
  else
  
    grid%ideriv = 0
  
    if (myrank == 0) write(*,'(" numerical jacobian")')
    ! We will compute the Jacobian via finite differences and store it.
    
    ! Create matrix with correct parallel layout and nonzero structure to 
    ! hold the Jacobian.
      ! MatFDColoringCreate() currently does not support MATMPIBAIJ.
    ! call DAGetMatrix(grid%da_ndof, MATMPIAIJ, grid%J, ierr)
    ! PETSc 2.1.6 introduces the MATAIJ matrix type.

!    if (grid%iblkfmt == 0) then
     call DAGetMatrix(grid%da_ndof, MATAIJ, grid%J, ierr)
 !   else
 !     call DAGetMatrix(grid%da_ndof, MATBAIJ, grid%J, ierr)
 !   endif
!   call MatSetOption(grid%J,MAT_COLUMN_ORIENTED,ierr)
        
#if ((PETSC_VERSION_RELEASE)&&(PETSC_VERSION_MAJOR==2)&&(PETSC_VERSION_MINOR==3)&&(PETSC_VERSION_SUBMINOR==2)&&(PETSC_VERSION_PATCH<=10))
    call DAGetColoring(grid%da_ndof, IS_COLORING_LOCAL, iscoloring, ierr)
#else
    call DAGetColoring(grid%da_ndof, IS_COLORING_GLOBAL, iscoloring, ierr)
#endif    
    call MatFDColoringCreate(grid%J, iscoloring, grid%matfdcoloring, ierr)
    
    call ISColoringDestroy(iscoloring, ierr)

    call MatFDColoringSetFunctionSNES(grid%matfdcoloring, imsResidual, &
                                      grid, ierr)
	
	 
    call MatFDColoringSetFromOptions(grid%matfdcoloring, ierr)
    call SNESSetJacobian(pflowsolv%snes, grid%J, grid%J, &
                         SNESDefaultComputeJacobianColor,  &
                         grid%matfdcoloring, ierr)
  endif

  if (myrank == 0) &
  write(*,'("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",/)')

  call SNESSetFunction(pflowsolv%snes, grid%r, IMSResidual, grid, ierr)


  ! Set the tolerances for the Newton solver.
  call SNESSetTolerances(pflowsolv%snes, grid%atol, grid%rtol, grid%stol, & 
                         grid%maxit, grid%maxf, ierr)

 

  call SNESSetFromOptions(pflowsolv%snes, ierr)
! call SNESLineSearchGetParams(grid%snes, alapha, maxstep, steptol, ierr) 
! call SNESLineSearchSetParams(grid%snes, alapha, maxstep, grid%stol, ierr) 
! print *, 'SNES setup ::', alapha,maxstep, steptol
 if (myrank == 0) write(*,'("  Finished setting up of SNES 1")')
  
!  call SNESLineSearchGetParams(grid%snes, alapha, maxstep, steptol, ierr) 
! if (myrank == 0) write(*,'("  Finished setting up of SNES 2")')
!  call SNESLineSearchSetParams(grid%snes, alapha, maxstep, grid%stol, ierr) 
! if (myrank == 0) write(*,'("  Finished setting up of SNES 3")')

!  call SNESGetKSP(grid%snes, grid%ksp, ierr)
!  call KSPSetTolerances(grid%ksp,grid%rtol,grid%atol,grid%dtol, &
!      10000,ierr)


 if (myrank == 0) &
  write(*,'("  Finished setting up of SNES ")')

end subroutine pflowgrid_Setup_SNES
!===================================================================================


  !-----------------------------------------------------------------------
  ! Set up boundary connections.
  !-----------------------------------------------------------------------
 
 subroutine pflowGrid_setup_BC(grid, locpat)	  
 implicit none
 type(pflowGrid) :: grid
 type(pflow_localpatch_info) :: locpat
 
interface
   subroutine pims_vecgetarrayf90(grid, patch, vec, f90ptr, ierr)
     use pflow_gridtype_module
     implicit none
#include "include/finclude/petsc.h"

     type(pflowGrid), intent(inout) :: grid
     type(pflow_localpatch_info) :: patch
     Vec :: vec
     PetscScalar, dimension(:), pointer :: f90ptr
     integer :: ierr
   end subroutine pims_vecgetarrayf90
end interface

 integer :: myrank, nc, ibc, ir, ii1,ii2,jj1, jj2, kk1, kk2, i,j,k
 integer :: m,n,ng,ird, ierr
 PetscScalar, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
								 						        
 
  myrank=grid%myrank
  
  
   locpat%nconnbc = 0
  if (grid%nx > 1 .and. grid%ny == 1 .and. grid%nz == 1) then
    if(grid%Samrai_drive == PETSC_TRUE)then
      if(pims_patch_at_boundary(grid,locpat, 0, 0) ==1)locpat%nconnbc = locpat%nconnbc + 1
      if(pims_patch_at_boundary(grid,locpat, 0, 1) ==1)locpat%nconnbc = locpat%nconnbc + 1
     else  
      if (locpat%nxs == locpat%ngxs) locpat%nconnbc = locpat%nconnbc + 1
      if (locpat%nxe == locpat%ngxe) locpat%nconnbc = locpat%nconnbc + 1
    endif
  else if (grid%nx == 1 .and. grid%ny == 1 .and. grid%nz > 1) then
    if(grid%Samrai_drive == PETSC_TRUE)then
      if(pims_patch_at_boundary(grid,locpat, 2, 0) ==1)locpat%nconnbc = locpat%nconnbc + 1
      if(pims_patch_at_boundary(grid,locpat, 2, 1) ==1)locpat%nconnbc = locpat%nconnbc + 1
     else  
      if (locpat%nzs == locpat%ngzs) locpat%nconnbc = locpat%nconnbc + 1
      if (locpat%nze == locpat%ngze) locpat%nconnbc = locpat%nconnbc + 1
     endif 
  elseif(grid%nx == 1 .and. grid%ny > 1 .and. grid%nz == 1) then
     print *, '1D problem: limited to x or z direction. STOP'; stop
  else
    if (grid%nx > 1) then
      if(grid%Samrai_drive == PETSC_TRUE)then
         if(pims_patch_at_boundary(grid,locpat, 0, 0) ==1)locpat%nconnbc = locpat%nconnbc + locpat%nlyz  
         if(pims_patch_at_boundary(grid,locpat, 0, 1) ==1)locpat%nconnbc = locpat%nconnbc + locpat%nlyz
       else 
        if (locpat%nxs == locpat%ngxs) locpat%nconnbc = locpat%nconnbc + locpat%nlyz
        if (locpat%nxe == locpat%ngxe) locpat%nconnbc = locpat%nconnbc + locpat%nlyz
       endif
    endif
    if (grid%ny > 1) then
      if(grid%Samrai_drive == PETSC_TRUE)then
         if(pims_patch_at_boundary(grid,locpat, 1, 0) ==1)locpat%nconnbc = locpat%nconnbc + locpat%nlxz  
         if(pims_patch_at_boundary(grid,locpat, 1, 1) ==1)locpat%nconnbc = locpat%nconnbc + locpat%nlxz
       else 
         if (locpat%nys == locpat%ngys) locpat%nconnbc = locpat%nconnbc + locpat%nlxz
         if (locpat%nye == locpat%ngye) locpat%nconnbc = locpat%nconnbc + locpat%nlxz
      endif 
    endif
    if (grid%nz > 1) then
       if(grid%Samrai_drive == PETSC_TRUE)then
         if(pims_patch_at_boundary(grid,locpat, 2, 0) ==1)locpat%nconnbc = locpat%nconnbc + locpat%nlxy  
         if(pims_patch_at_boundary(grid,locpat, 2, 1) ==1)locpat%nconnbc = locpat%nconnbc + locpat%nlxy
       else 
         if (locpat%nzs == locpat%ngzs) locpat%nconnbc = locpat%nconnbc + locpat%nlxy
         if (locpat%nze == locpat%ngze) locpat%nconnbc = locpat%nconnbc + locpat%nlxy
       endif
    endif
  endif
      
  write(*,'(" --> pflowconn: rank = ",i4, &
 &", boundary connections =", i6)') myrank,locpat%nconnbc



  if (locpat%nconnbc > 0) then
    allocate(locpat%mblkbc(locpat%nconnbc))
    allocate(locpat%ibconn(locpat%nconnbc))
    allocate(locpat%distbc(locpat%nconnbc))
    allocate(locpat%areabc(locpat%nconnbc))
    allocate(locpat%ipermbc(locpat%nconnbc))
    allocate(locpat%delzbc(locpat%nconnbc))
    
!   allocate(grid%velocitybc(grid%nphase,locpat%nconnbc))
    
    allocate(locpat%vlbc(locpat%nconnbc))
    allocate(locpat%vvlbc(locpat%nconnbc))
    allocate(locpat%vgbc(locpat%nconnbc))
    allocate(locpat%vvgbc(locpat%nconnbc))
	
	locpat%vlbc=0.D0
	locpat%vgbc=0.D0
  
    allocate(locpat%varbc(1:(grid%ndof+1)*(2+4*grid%nphase )))
    endif

    

 call pims_vecgetarrayf90(grid, locpat, grid%dx_loc, dx_loc_p, ierr)
  call pims_vecgetarrayf90(grid, locpat, grid%dy_loc, dy_loc_p, ierr)
  call pims_vecgetarrayf90(grid, locpat, grid%dz_loc, dz_loc_p, ierr)




  nc = 0 
  if ( locpat%nconnbc > 0) then

    ! calculate boundary conditions locally on only those processors which 
    ! contain a boundary!

    do ibc = 1, grid%nblkbc
      do ir = grid%iregbc1(ibc), grid%iregbc2(ibc)
        kk1 = grid%k1bc(ir) - locpat%nzs
        kk2 = grid%k2bc(ir) - locpat%nzs
        jj1 = grid%j1bc(ir) - locpat%nys
        jj2 = grid%j2bc(ir) - locpat%nys
        ii1 = grid%i1bc(ir) - locpat%nxs
        ii2 = grid%i2bc(ir) - locpat%nxs

        kk1 = max(1,kk1)
        kk2 = min(locpat%nlz,kk2)
        jj1 = max(1,jj1)
        jj2 = min(locpat%nly,jj2)
        ii1 = max(1,ii1)
        ii2 = min(locpat%nlx,ii2)

        if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle 

        do k = kk1,kk2
          do j = jj1,jj2
            do i = ii1,ii2
              nc = nc + 1
                m = i+(j-1)*locpat%nlx+(k-1)*locpat%nlxy
              locpat%mblkbc(nc) = m  ! m is a local index
              locpat%ibconn(nc) = ibc
              ng = locpat%nL2G(m)
              ! Use ghosted index to access dx, dy, dz because we have
              ! already done a global-to-local scatter for computing the
              ! interior node connections.
        
!       print *,'pflowgrid_mod: ',nc,ibc,ir,m,ng,ii1,ii2,kk1,kk2,locpat%nblkbc,grid%igeom
        
              select case(grid%igeom)
              case(1) ! cartesian
                if (grid%iface(ibc) == 1) then
                  locpat%distbc(nc) = 0.5d0*dx_loc_p(ng)
                  locpat%areabc(nc) = dy_loc_p(ng)*dz_loc_p(ng)
                  locpat%ipermbc(nc) = 1
                  locpat%delzbc(nc) = 0.d0
                else if (grid%iface(ibc) == 2) then
                  locpat%distbc(nc) = 0.5d0*dx_loc_p(ng)
                  locpat%areabc(nc) = dy_loc_p(ng)*dz_loc_p(ng)
                  locpat%ipermbc(nc) = 1
                  locpat%delzbc(nc) = 0.d0
                else if (grid%iface(ibc) == 3) then
                  locpat%distbc(nc) = 0.5d0*dz_loc_p(ng)
                  locpat%areabc(nc) = dx_loc_p(ng)*dy_loc_p(ng)
                  locpat%ipermbc(nc) = 2
                  locpat%delzbc(nc) = locpat%distbc(nc)
                else if (grid%iface(ibc) == 4) then
                  locpat%distbc(nc) = 0.5d0*dz_loc_p(ng)
                  locpat%areabc(nc) = dx_loc_p(ng)*dy_loc_p(ng)
                  locpat%ipermbc(nc) = 2
                  locpat%delzbc(nc) = -locpat%distbc(nc)
                else if (grid%iface(ibc) == 5) then
                  locpat%distbc(nc) = 0.5d0*dy_loc_p(ng)
                  locpat%areabc(nc) = dx_loc_p(ng)*dz_loc_p(ng)
                  locpat%ipermbc(nc) = 3
                  locpat%delzbc(nc) = 0.d0
                else if (grid%iface(ibc) == 6) then
                  locpat%distbc(nc) = 0.5d0*dy_loc_p(ng)
                  locpat%areabc(nc) = dx_loc_p(ng)*dz_loc_p(ng)
                  locpat%ipermbc(nc) = 3
                  locpat%delzbc(nc) = 0.d0
                endif
              case(2) ! cylindrical
                ird= mod(mod((m),locpat%nlxy),locpat%nlx) + locpat%nxs 
                if (grid%iface(ibc) == 1) then
                  locpat%distbc(nc) = 0.5d0*dx_loc_p(ng)
                  locpat%areabc(nc) = 2.0D0* Pi * grid%rd(ird-1)*dz_loc_p(ng)
                  locpat%ipermbc(nc) = 1
                  locpat%delzbc(nc) = 0.d0
                else if (grid%iface(ibc) == 2) then
                  locpat%distbc(nc) = 0.5d0*dx_loc_p(ng)
                  locpat%areabc(nc) =  2.0D0* Pi * grid%rd(ird)*dz_loc_p(ng)
                  locpat%ipermbc(nc) = 1
                  locpat%delzbc(nc) = 0.d0
                else if (grid%iface(ibc) == 3) then
                  locpat%distbc(nc) = 0.5d0*dz_loc_p(ng)
                  locpat%areabc(nc) =  Pi * (grid%rd(ird)+ grid%rd(ird-1))*&
                       (grid%rd(ird) - grid%rd(ird-1))	
                  locpat%ipermbc(nc) = 2
                  locpat%delzbc(nc) = locpat%distbc(nc)
                else if (grid%iface(ibc) == 4) then
                  locpat%distbc(nc) = 0.5d0*dz_loc_p(ng)
                  locpat%areabc(nc) = Pi * (grid%rd(ird)+ grid%rd(ird-1))*&
                       (grid%rd(ird) - grid%rd(ird-1))
                  locpat%ipermbc(nc) = 2
                  locpat%delzbc(nc) = -locpat%distbc(nc)
                endif	
              case(3) ! spherical
              end select
            enddo ! i
          enddo ! j
        enddo ! k
      enddo ! ir
    enddo ! ibc
  endif
  
  
  if(grid%Samrai_drive==PETSC_FALSE) then
     call VecRestoreArrayF90(grid%dx_loc, dx_loc_p, ierr)
     call VecRestoreArrayF90(grid%dy_loc, dy_loc_p, ierr)
     call VecRestoreArrayF90(grid%dz_loc, dz_loc_p, ierr)
  endif

  if(locpat%nconnbc .ne. nc) then
    write(*,*) 'Error in computing boundary connections: ', &
      'rank = ',myrank,' nconnbc = ',locpat%nconnbc,' nc = ',nc
    stop
  endif
 
  
   !Hydro here  
!   print *,'Finished Hydro'

  !-----------------------------------------------------------------------
  ! Set up boundary conditions at interfaces
  !-----------------------------------------------------------------------
  allocate(locpat%velocitybc(grid%nphase, locpat%nconnbc))
     allocate(locpat%xxbc(grid%ndof,locpat%nconnbc))
!     allocate(grid%iphasebc(locpat%nconnbc))
	 do nc = 1, locpat%nconnbc
        ibc = locpat%ibconn(nc)
        locpat%xxbc(:,nc)=grid%xxbc0(:,ibc)
        locpat%velocitybc(:,nc) = grid%velocitybc0(:,ibc)
     enddo
   ! if(grid%using_pflowGrid == PETSC_FALSE)then
	 deallocate(grid%xxbc0)
   ! endif
	!if(grid%iread_init==2) call Boundary_adjustment(grid)
	
   deallocate(grid%velocitybc0)

end subroutine 
!======================================================================


!#include "pflowgrid_compute_xyz.F90"

  subroutine pflowGrid_compute_xyz(grid)
  
  implicit none
  
  type(pflowGrid), intent(inout) :: grid

! integer :: ierr
  integer :: i, j, k, n
  integer :: prevnode
    ! prevnode is used to hold the natural numbering of the node that is the
    ! previous node in either the x, y, or z direction.
! Vec :: dx_nat, dy_nat, dz_nat
! Vec :: dx_all, dy_all, dz_all  ! Holds contents of dx_nat et al. on proc 0.
! real*8, pointer :: dx_p(:), dy_p(:), dz_p(:)


!  if(grid%myrank == 0) then
!   note changed by clu 2005/08/24/17:06, now every node have coordinate of the whole domain
!    ordered in natural     
!   call VecGetArrayF90(dx_all, dx_p, ierr)
!   call VecGetArrayF90(dy_all, dy_p, ierr)
!   call VecGetArrayF90(dz_all, dz_p, ierr)

!   num_nodes = grid%nx * grid%ny * grid%nz
!   allocate(grid%x(num_nodes))
!   allocate(grid%y(num_nodes))
!   allocate(grid%z(num_nodes))
    n = 0
    do k=1, grid%nz
      do j=1, grid%ny
        do i=1, grid%nx
          n = n + 1
    
!   print *,'compute-xyz: ',n,i,j,k,dx_p(n),dy_p(n),dz_p(n)

          if(i == 1) then
!           grid%x(n) = 0.5d0 * dx_p(n) 
            grid%x(n) = 0.5d0 * grid%dx0(i) 
          else
            prevnode = n-1
!           grid%x(n) = grid%x(prevnode) + 0.5d0*(dx_p(prevnode) + dx_p(n))
            grid%x(n) = grid%x(prevnode) + 0.5d0*(grid%dx0(i-1) + grid%dx0(i))
          endif

          if(j == 1) then
!           grid%y(n) = 0.5d0 * dy_p(n)
            grid%y(n) = 0.5d0 * grid%dy0(j)
          else
!           prevnode = i + (j-2)*grid%nx + (k-1)*grid%nx*grid%ny
            prevnode = n - grid%nx
!           grid%y(n) = grid%y(prevnode) + 0.5d0*(dy_p(prevnode) + dy_p(n))
            grid%y(n) = grid%y(prevnode) + 0.5d0*(grid%dy0(j-1) + grid%dy0(j))
          endif

          if(k == 1) then
!           grid%z(n) = 0.5d0 * dz_p(n)
            grid%z(n) = 0.5d0 * grid%dz0(k)
          else
!           prevnode = i + (j-1)*grid%nx + (k-2)*grid%nx*grid%ny
            prevnode = n - grid%nx*grid%ny
!           grid%z(n) = grid%z(prevnode) + 0.5d0*(dz_p(prevnode) + dz_p(n))
            grid%z(n) = grid%z(prevnode) + 0.5d0*(grid%dz0(k-1) + grid%dz0(k))
          endif
          
!         print *,'compute-xyz: ',n,i,j,k,grid%x(n),grid%y(n),grid%z(n), &
!         dx_p(n),dy_p(n),dz_p(n)
        enddo
      enddo
    enddo
!   call VecRestoreArrayF90(dx_all, dx_p, ierr)
!   call VecRestoreArrayF90(dy_all, dy_p, ierr)
!   call VecRestoreArrayF90(dz_all, dz_p, ierr)
!  endif
  
 
  
  
  
  
! call VecDestroy(dx_nat, ierr)
! call VecDestroy(dy_nat, ierr)
! call VecDestroy(dz_nat, ierr)
! call VecDestroy(dx_all, ierr)
! call VecDestroy(dy_all, ierr)
! call VecDestroy(dz_all, ierr)
  
  end subroutine pflowGrid_compute_xyz

!================================================================================



!#include "pflowjacobian.F90"

subroutine ComputeMFJacobian(snes, x, J, B, flag, ctx, ierr)
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: x
  Mat, intent(out) :: J, B
  MatStructure, intent(in) :: flag
  integer, intent(inout) :: ctx(*)
  integer, intent(out) :: ierr

  call MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY, ierr)
  B = J
end subroutine ComputeMFJacobian

!======================================================================

!#include "pflowgrid_monitors.F90"

subroutine pflowgrid_MonitorH(snes, its, norm, grid)
  
  implicit none

  SNES, intent(in) :: snes
  integer, intent(in) :: its
  PetscReal, intent(in) :: norm
  type(pflowGrid), intent(in) :: grid
  
  integer :: ierr
  integer :: myrank
  PetscScalar :: h

#if ((PETSC_VERSION_RELEASE)&&(PETSC_VERSION_MAJOR==2)&&(PETSC_VERSION_MINOR==3)&&(PETSC_VERSION_SUBMINOR==2)&&(PETSC_VERSION_PATCH<=10))  
  call MatSNESMFGetH(grid%J, h, ierr)
#else
  call MatMFFDGetH(grid%J, h, ierr)
#endif
  call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)

  if(myrank == 0) then
    write(*,*) "#At SNES iteration ", its, "h is ", h
  endif
end subroutine pflowgrid_MonitorH

!======================================================================





end module pflow_grid_module
