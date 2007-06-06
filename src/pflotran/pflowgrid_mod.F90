module pflow_grid_module
  
  use pflow_gridtype_module
  
  implicit none
  
  private
  real*8, parameter :: Pi=3.1415926D0 

  public pflowGrid
  public pflowGrid_new
  public pflowGrid_setup
  public pflowgrid_update_dt
  public pflowGrid_step
  public pflowGrid_update
  public pflowGrid_read_input
  public pflowGrid_get_t
  public pflowGrid_setvel
  public pflowGrid_destroy
  public pflowGrid_parse_cmdline
  
#include "definitions.h"

! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petsc.h"
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
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

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
type(pflowGrid) function pflowGrid_new(igeom, nx, ny, nz, npx, npy, npz, &
                                       nphase, nspec,npricomp,ndof, icouple, &
                                       idcdm, itable)

  implicit none

  integer*4, intent(in) :: igeom, nx, ny, nz, npx, npy, npz
  integer, intent(in) :: nphase, ndof,icouple,idcdm,itable,nspec,npricomp
!  integer, intent(in) :: equ_option
  integer*4 :: n, ng, na
  integer*4 :: i, j, k
  integer :: ierr
  type(pflowGrid) grid
  
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
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_liquid", &
                           grid%use_liquid, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_cond", &
                           grid%use_cond, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_th", &
                           grid%use_th, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_thc", &
                           grid%use_thc, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_2ph", &
                           grid%use_2ph, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_mph", &
                           grid%use_mph, ierr)
   call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_flash", &
                           grid%use_flash, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_owg", &
                           grid%use_owg, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_vadose", &
                           grid%use_vadose, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_richard", &
                           grid%use_richard, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_debug", &
                           grid%use_debug, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_ksp", &
                           grid%use_ksp, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_isoth", &
                           grid%use_isoth, ierr)

  if (icouple == 0) then
    grid%using_pflowGrid = PETSC_FALSE
  else
    grid%using_pflowGrid = PETSC_TRUE
  endif
      
  grid%igeom = igeom
  grid%nx = nx
  grid%ny = ny
  grid%nz = nz
  grid%nxy = nx*ny
  grid%nmax = nx*ny*nz
      
  grid%npx = npx
  grid%npy = npy
  grid%npz = npz

  grid%nphase = nphase
  grid%ndof = ndof
  grid%nspec = nspec
  grid%npricomp = npricomp

  grid%idcdm = idcdm
      
  grid%itable = itable
      
  !set specific phase indices
  grid%jh2o = 1
  select case(grid%nphase)
    case(2)
      if (grid%use_2ph == PETSC_TRUE) then
        grid%jgas=  2
        grid%jco2 = 3
      else
        grid%jco2 = 2
        grid%jgas =3 
      endif
    case(3)
      grid%jco2 = 2
      grid%jgas =3 
  end select 

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
  allocate(grid%tfac(13))
      
  grid%tfac(1)  = 2.0d0; grid%tfac(2)  = 2.0d0
  grid%tfac(3)  = 2.0d0; grid%tfac(4)  = 2.0d0
  grid%tfac(5)  = 2.0d0; grid%tfac(6)  = 1.8d0
  grid%tfac(7)  = 1.6d0; grid%tfac(8)  = 1.4d0
  grid%tfac(9)  = 1.2d0; grid%tfac(10) = 1.0d0
  grid%tfac(11) = 1.0d0; grid%tfac(12) = 1.0d0
  grid%tfac(13) = 1.0d0      
  grid%use_matrix_free = 1
  grid%newton_max = 16
  grid%icut_max = 16
  grid%iaccel = 1
  grid%dt = 1.d0
  grid%dt_min = 1.d0
  grid%dt_max = 3.1536d6 ! One-tenth of a year.
  grid%atol = PETSC_DEFAULT_DOUBLE_PRECISION
  grid%rtol = PETSC_DEFAULT_DOUBLE_PRECISION
  grid%stol = PETSC_DEFAULT_DOUBLE_PRECISION
  grid%maxit = PETSC_DEFAULT_INTEGER
  grid%maxf = PETSC_DEFAULT_INTEGER
  
  grid%ihydrostatic = 0
  grid%conc0 = 1.d-6
  
  grid%dpmxe = 5.d4
  grid%dtmpmxe = 2.d0
  grid%dsmxe = 5.d0
  grid%dcmxe = 5.d0

  !physical constants
  grid%difaq = 1.d-9 ! m^2/s read from input file
  grid%delhaq = 12.6d0 ! kJ/mol read from input file
  grid%gravity = 9.8068d0    ! m/s^2
  grid%tref   = 50.D0
  grid%fmwh2o = 18.01534d0 ! kg H2O/mol H2O
  grid%fmwco2 = 44.0098d0
  grid%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  !-----------------------------------------------------------------------
  ! Generate the DA objects that will manage communication.
  !-----------------------------------------------------------------------
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       nx,ny,nz,npx,npy,npz,1,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       grid%da_1_dof,ierr)

! call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
!      nx,ny,nz,npx,npy,npz,3,1, &
!      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
!      grid%da_3_dof,ierr)
 
  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       nx,ny,nz,npx,npy,npz,grid%nphase,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       grid%da_nphase_dof,ierr)

  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       nx,ny,nz,npx,npy,npz,3*grid%nphase,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       grid%da_3np_dof,ierr)

  call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
       nx,ny,nz,npx,npy,npz,grid%ndof,1, &
       PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
       grid%da_ndof,ierr)

  if (grid%use_2ph == PETSC_TRUE) then
    print *,' 2ph create DA'
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    nx,ny,nz,npx,npy,npz,grid%nphase*grid%npricomp,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_nphancomp_dof,ierr)

    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    nx,ny,nz,npx,npy,npz,grid%nphase*grid%nspec,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_nphanspec_dof,ierr)

    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    nx,ny,nz,npx,npy,npz, &
                    grid%nphase*grid%nspec*grid%npricomp,1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_nphanspecncomp_dof,ierr)
  endif
       
  if (grid%use_mph == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    nx,ny,nz,npx,npy,npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)

  endif
  
 if (grid%use_richard == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    nx,ny,nz,npx,npy,npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)

  endif
    
  if (grid%use_vadose == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    nx,ny,nz,npx,npy,npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)
  endif

  if (grid%use_flash == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    nx,ny,nz,npx,npy,npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)
  endif
     
  if (grid%use_owg == PETSC_TRUE) then
    call DACreate3D(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR, &
                    nx,ny,nz,npx,npy,npz,(grid%ndof+1)*(2+7*grid%nphase + 2* &
                    grid%nspec*grid%nphase),1, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    grid%da_var_dof,ierr)
  endif

  !print DA info for each processor
  call DAView(grid%da_ndof,PETSC_VIEWER_STDOUT_WORLD,ierr)

 !-----------------------------------------------------------------------
 ! Create the vectors with parallel layout corresponding to the DA's,
 ! and, for vectors that need to be ghosted, create the corresponding
 ! ghosted vectors.
 !-----------------------------------------------------------------------

  ! 1 degree of freedom
  call DACreateGlobalVector(grid%da_1_dof, grid%porosity, ierr)
  call VecDuplicate(grid%porosity, grid%porosity0, ierr)
  call VecDuplicate(grid%porosity, grid%tor, ierr)
  call VecDuplicate(grid%porosity, grid%conc, ierr)
  call VecDuplicate(grid%porosity, grid%dx, ierr)
  call VecDuplicate(grid%porosity, grid%dy, ierr)
  call VecDuplicate(grid%porosity, grid%dz, ierr)
  call VecDuplicate(grid%porosity, grid%volume, ierr)
  call VecDuplicate(grid%porosity, grid%ithrm, ierr)
  call VecDuplicate(grid%porosity, grid%icap, ierr)
  call VecDuplicate(grid%porosity, grid%iphas, ierr)
  call VecDuplicate(grid%porosity, grid%iphas_old, ierr)
  call VecDuplicate(grid%porosity, grid%temp, ierr)
  call VecDuplicate(grid%porosity, grid%ttemp, ierr)
  call VecDuplicate(grid%porosity, grid%phis, ierr)

  call VecDuplicate(grid%porosity, grid%perm_xx, ierr)
  call VecDuplicate(grid%porosity, grid%perm_yy, ierr)
  call VecDuplicate(grid%porosity, grid%perm_zz, ierr)
  call VecDuplicate(grid%porosity, grid%perm0_xx, ierr)
  call VecDuplicate(grid%porosity, grid%perm0_yy, ierr)
  call VecDuplicate(grid%porosity, grid%perm0_zz, ierr)
  call VecDuplicate(grid%porosity, grid%perm_pow, ierr)
      
  call DACreateLocalVector(grid%da_1_dof,grid%dx_loc,ierr)
  call VecDuplicate(grid%dx_loc, grid%dy_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%dz_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%porosity_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%tor_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%ithrm_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%icap_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%iphas_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%ttemp_loc, ierr)

  call VecDuplicate(grid%dx_loc, grid%perm_xx_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%perm_yy_loc, ierr)
  call VecDuplicate(grid%dx_loc, grid%perm_zz_loc, ierr)

  ! 3 degrees of freedom
! call DACreateGlobalVector(grid%da_3_dof, grid%perm, ierr)
! call DACreateLocalVector(grid%da_3_dof, grid%perm_loc, ierr)

  ! nphase degrees of freedom
  call DACreateGlobalVector(grid%da_nphase_dof, grid%pressure, ierr)
  call VecDuplicate(grid%pressure, grid%ppressure, ierr)
  call VecDuplicate(grid%pressure, grid%ssat, ierr)
  call VecDuplicate(grid%pressure, grid%sat, ierr)
  call VecDuplicate(grid%pressure, grid%dp, ierr)
  call VecDuplicate(grid%pressure, grid%density, ierr)
  call VecDuplicate(grid%pressure, grid%ddensity, ierr)
  call VecDuplicate(grid%pressure, grid%avgmw, ierr)
  call VecDuplicate(grid%pressure, grid%d_p, ierr)
  call VecDuplicate(grid%pressure, grid%d_t, ierr)
  call VecDuplicate(grid%pressure, grid%h, ierr)
  call VecDuplicate(grid%pressure, grid%hh, ierr)
  call VecDuplicate(grid%pressure, grid%h_p, ierr)
  call VecDuplicate(grid%pressure, grid%h_t, ierr)
  call VecDuplicate(grid%pressure, grid%viscosity, ierr)
  call VecDuplicate(grid%pressure, grid%v_p, ierr)
  call VecDuplicate(grid%pressure, grid%v_t, ierr)
 ! xmol may not be nphase DOF, need change later 
  call VecDuplicate(grid%pressure, grid%xxmol, ierr)
  call VecDuplicate(grid%pressure, grid%xmol, ierr)


  call DACreateLocalVector(grid%da_nphase_dof, grid%ppressure_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%ssat_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%xxmol_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%ddensity_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%avgmw_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%d_p_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%d_t_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%hh_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%h_p_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%h_t_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%viscosity_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%v_p_loc, ierr)
  call VecDuplicate(grid%ppressure_loc, grid%v_t_loc, ierr)
      
  ! 3 * nphase degrees of freedom (velocity vector)
  call DACreateGlobalVector(grid%da_3np_dof, grid%vl, ierr)
      
! print *,'pflowgrid_new: ',grid%using_pflowGrid
      
  if (grid%using_pflowGrid == PETSC_TRUE) &
    call VecDuplicate(grid%vl, grid%vvl, ierr)
      
  ! nvar * nphase degrees of freedom
!  call DACreateGlobalVector(grid%da_ncnp_dof, grid%xmol, ierr)


  if (grid%use_2ph == PETSC_TRUE) then
    print *,'2ph add var'
    call VecDuplicate(grid%sat, grid%d_s, ierr)
    call VecDuplicate(grid%sat, grid%h_s, ierr)
    call VecDuplicate(grid%sat, grid%u, ierr)
    call VecDuplicate(grid%sat, grid%uu, ierr)
    call VecDuplicate(grid%sat, grid%u_s, ierr)
    call VecDuplicate(grid%sat, grid%u_p, ierr)
    call VecDuplicate(grid%sat, grid%u_t, ierr)

    call VecDuplicate(grid%ssat_loc, grid%d_s_loc, ierr)
    call VecDuplicate(grid%ssat_loc, grid%h_s_loc, ierr)
    call VecDuplicate(grid%ssat_loc, grid%u_s_loc, ierr)

    call VecDuplicate(grid%sat, grid%pcw, ierr)
    call VecDuplicate(grid%sat, grid%pc_p, ierr)
    call VecDuplicate(grid%sat, grid%pc_t, ierr)
    call VecDuplicate(grid%sat, grid%pc_s, ierr)
    call VecDuplicate(grid%ssat_loc, grid%pcw_loc, ierr)
    call VecDuplicate(grid%ssat_loc, grid%pc_p_loc, ierr)
    call VecDuplicate(grid%ssat_loc, grid%pc_t_loc, ierr)
    call VecDuplicate(grid%ssat_loc, grid%pc_s_loc, ierr)


    call VecDuplicate(grid%sat, grid%kvr, ierr)
    call VecDuplicate(grid%sat, grid%kvr_p, ierr)
    call VecDuplicate(grid%sat, grid%kvr_t, ierr)
    call VecDuplicate(grid%sat, grid%kvr_s, ierr)
    call VecDuplicate(grid%ssat_loc, grid%kvr_loc, ierr)
    call VecDuplicate(grid%ssat_loc, grid%kvr_p_loc, ierr)
    call VecDuplicate(grid%ssat_loc, grid%kvr_t_loc, ierr)
    call VecDuplicate(grid%ssat_loc, grid%kvr_s_loc, ierr)

    call DACreateGlobalVector(grid%da_nphancomp_dof, grid%h_c, ierr)
    call VecDuplicate(grid%h_c, grid%u_c, ierr)
    call VecDuplicate(grid%h_c, grid%avgmw_c, ierr)
    call VecDuplicate(grid%h_c, grid%d_c, ierr)
    call VecDuplicate(grid%h_c, grid%pc_c, ierr)
    call VecDuplicate(grid%h_c, grid%kvr_c, ierr)
        
    call DACreateLocalVector(grid%da_nphancomp_dof, grid%h_c_loc, ierr)
    call VecDuplicate(grid%h_c_loc, grid%avgmw_c_loc, ierr)
    call VecDuplicate(grid%h_c_loc, grid%d_c_loc, ierr)
    call VecDuplicate(grid%h_c_loc, grid%pc_c_loc, ierr)
    call VecDuplicate(grid%h_c_loc, grid%kvr_c_loc, ierr)


    call DACreateGlobalVector(grid%da_nphanspec_dof, grid%hen, ierr)
    call DACreateLocalVector(grid%da_nphanspec_dof, grid%hen_loc, ierr)
    call VecDuplicate(grid%hen, grid%hen_p, ierr)
    call VecDuplicate(grid%hen_loc, grid%hen_p_loc, ierr)
    call VecDuplicate(grid%hen, grid%hen_t, ierr)
    call VecDuplicate(grid%hen_loc, grid%hen_t_loc, ierr)
    call VecDuplicate(grid%hen, grid%hen_s, ierr)
    call VecDuplicate(grid%hen_loc, grid%hen_s_loc, ierr)
    call VecDuplicate(grid%hen, grid%df, ierr)
    call VecDuplicate(grid%hen_loc, grid%df_loc, ierr)
    call VecDuplicate(grid%df, grid%df_p, ierr)
    call VecDuplicate(grid%df_loc, grid%df_p_loc, ierr)
    call VecDuplicate(grid%df, grid%df_t, ierr)
    call VecDuplicate(grid%df_loc, grid%df_t_loc, ierr)
    call VecDuplicate(grid%df, grid%df_s, ierr)
    call VecDuplicate(grid%df_loc, grid%df_s_loc, ierr)
  
    call DACreateGlobalVector(grid%da_nphanspecncomp_dof, grid%hen_c, ierr)
    call DACreateLocalVector(grid%da_nphanspecncomp_dof, grid%hen_c_loc, ierr)
      
    call VecDuplicate(grid%hen_c, grid%df_c, ierr)
    call VecDuplicate(grid%hen_c_loc, grid%df_c_loc, ierr)
  end if  

  if (grid%use_mph == PETSC_TRUE) then
    call DACreateGlobalVector(grid%da_var_dof, grid%var, ierr)
    call DACreateLocalVector(grid%da_var_dof, grid%var_loc, ierr)
  endif

  if (grid%use_richard == PETSC_TRUE) then
    call DACreateGlobalVector(grid%da_var_dof, grid%var, ierr)
    call DACreateLocalVector(grid%da_var_dof, grid%var_loc, ierr)
  endif
 
 if (grid%use_flash == PETSC_TRUE) then
    call DACreateGlobalVector(grid%da_var_dof, grid%var, ierr)
    call DACreateLocalVector(grid%da_var_dof, grid%var_loc, ierr)
  endif
     
  if (grid%use_owg == PETSC_TRUE) then
    call DACreateGlobalVector(grid%da_var_dof, grid%var, ierr)
    call DACreateLocalVector(grid%da_var_dof, grid%var_loc, ierr)
  endif  
  if (grid%use_vadose == PETSC_TRUE) then
    call DACreateGlobalVector(grid%da_var_dof, grid%var, ierr)
    call DACreateLocalVector(grid%da_var_dof, grid%var_loc, ierr)
  endif



      ! ndof degrees of freedom
  call DACreateGlobalVector(grid%da_ndof, grid%xx, ierr)
  call VecDuplicate(grid%xx, grid%yy, ierr)
  call VecDuplicate(grid%xx, grid%dxx, ierr)
  call VecDuplicate(grid%xx, grid%r, ierr)
  call VecDuplicate(grid%xx, grid%accum, ierr)
  
  call VecSetBlocksize(grid%dxx, grid%ndof, ierr)

  call DACreateLocalVector(grid%da_ndof, grid%xx_loc, ierr)
  
  ! Create Natural Vec for output: use VecDuplicate here?
!  call DACreateNaturalVector(grid%da_1_dof,      grid%c_nat,    ierr)
!  call DACreateNaturalVector(grid%da_1_dof,      grid%phis_nat, ierr)
!  call DACreateNaturalVector(grid%da_1_dof,      grid%t_nat,    ierr)
!  call DACreateNaturalVector(grid%da_1_dof,      grid%por_nat,  ierr)
!  call DACreateNaturalVector(grid%da_nphase_dof, grid%p_nat,    ierr)
!  call DACreateNaturalVector(grid%da_nphase_dof, grid%s_nat,    ierr)
!  call DACreateNaturalVector(grid%da_3np_dof,    grid%vl_nat,   ierr)

!   if (grid%use_2ph == PETSC_TRUE) &
!     call DACreateNaturalVector(grid%da_nphase_dof, grid%x_nat, ierr)
!-----------------------------------------------------------------------
  ! Set up PETSc nonlinear solver context.
!-----------------------------------------------------------------------
  call SNESCreate(PETSC_COMM_WORLD, grid%snes, ierr)
  CHKERRQ(ierr)
!-----------------------------------------------------------------------
  ! Set up information about corners of local domain.
!-----------------------------------------------------------------------

  call DAGetCorners(grid%da_nphase_dof, grid%nxs, &
                    grid%nys, grid%nzs, grid%nlx, &
                    grid%nly, grid%nlz, ierr)

  grid%nxe = grid%nxs + grid%nlx
  grid%nye = grid%nys + grid%nly
  grid%nze = grid%nzs + grid%nlz
  grid%nlxy = grid%nlx * grid%nly
  grid%nlxz = grid%nlx * grid%nlz
  grid%nlyz = grid%nly * grid%nlz
  grid%nlmax = grid%nlx * grid%nly * grid%nlz
  grid%nldof = grid%nlmax * grid%nphase

  call DAGetGhostCorners(grid%da_nphase_dof, grid%ngxs, &
                         grid%ngys, grid%ngzs, grid%ngx, &
                         grid%ngy, grid%ngz, ierr)

  grid%ngxe = grid%ngxs + grid%ngx
  grid%ngye = grid%ngys + grid%ngy
  grid%ngze = grid%ngzs + grid%ngz
  grid%ngxy = grid%ngx * grid%ngy
  grid%ngxz = grid%ngx * grid%ngz
  grid%ngyz = grid%ngy * grid%ngz
  grid%ngmax = grid%ngx * grid%ngy * grid%ngz
  grid%ngdof = grid%ngmax * grid%nphase

!-----------------------------------------------------------------------
  ! Determine number of local connections.
!-----------------------------------------------------------------------
  grid%nconn = (grid%ngx-1) * grid%nly * &
               grid%nlz + grid%nlx * (grid%ngy-1) * &
               grid%nlz + grid%nlx * grid%nly * &
               (grid%ngz-1)

!-----------------------------------------------------------------------
      ! Allocate memory for allocatable arrays.
!-----------------------------------------------------------------------
  allocate(grid%nL2G(grid%nlmax))
  allocate(grid%nG2L(grid%ngmax))
  allocate(grid%nL2A(grid%nlmax))
  allocate(grid%nG2N(grid%ngmax))
  allocate(grid%nd1(grid%nconn))
  allocate(grid%nd2(grid%nconn))
  allocate(grid%dist1(grid%nconn))
  allocate(grid%dist2(grid%nconn))
  allocate(grid%area(grid%nconn))
  allocate(grid%delz(grid%nconn))
  allocate(grid%grav_ang(grid%nconn))

  allocate(grid%iperm1(grid%nconn))
  allocate(grid%iperm2(grid%nconn))

  allocate(grid%vl_loc(grid%nconn))
  allocate(grid%vvl_loc(grid%nconn))
  allocate(grid%vg_loc(grid%nconn))
  allocate(grid%vvg_loc(grid%nconn))
    
    
  grid%vl_loc = 0.D0
  grid%vg_loc = 0.D0
      
  allocate(grid%xphi_co2(grid%nlmax))
  allocate(grid%xxphi_co2(grid%nlmax))
  allocate(grid%den_co2(grid%nlmax))
  allocate(grid%dden_co2(grid%nlmax))

  grid%xphi_co2 = 1.d0
  grid%xxphi_co2 = 1.d0
  grid%den_co2 = 1.d0
  grid%dden_co2 = 1.d0

! if (grid%using_pflowGrid == PETSC_TRUE) &
! allocate(grid%vvl_loc(grid%nconn*grid%nphase))

  ! I don't like having a fixed number of boundary condition regions.
  ! Memory for these arrays ought to allocated by parsing the input file
  ! to determine the number of regions.  This is the lazy way... I 
  ! should fix it eventually.
  ! The same goes for the number of BC blocks.
  allocate(grid%iregbc1(MAXBCREGIONS))
  allocate(grid%iregbc2(MAXBCREGIONS))
  allocate(grid%ibndtyp(MAXBCREGIONS))
  allocate(grid%iface(MAXBCREGIONS))
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
  allocate(grid%qsrc(MAXSRCTIMES,MAXSRC))
  allocate(grid%csrc(MAXSRCTIMES,MAXSRC))
      
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

  if (grid%use_mph == PETSC_TRUE .or. grid%use_owg == PETSC_TRUE &
      .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
      .or. grid%use_richard == PETSC_TRUE) then
    allocate(grid%xx_ini(grid%ndof,MAXINITREGIONS))
    allocate(grid%iphas_ini(MAXINITREGIONS))
  else
    allocate(grid%pres_ini(MAXINITREGIONS))
    allocate(grid%temp_ini(MAXINITREGIONS))
    allocate(grid%sat_ini(MAXINITREGIONS))
    allocate(grid%xmol_ini(MAXINITREGIONS))
    allocate(grid%conc_ini(MAXINITREGIONS))
  endif
    
  allocate(grid%i1brk(MAXINITREGIONS))
  allocate(grid%i2brk(MAXINITREGIONS))
  allocate(grid%j1brk(MAXINITREGIONS))
  allocate(grid%j2brk(MAXINITREGIONS))
  allocate(grid%k1brk(MAXINITREGIONS))
  allocate(grid%k2brk(MAXINITREGIONS))
  allocate(grid%ibrktyp(MAXINITREGIONS))
  allocate(grid%ibrkface(MAXINITREGIONS))
  
  if (idcdm == 1) then
    allocate(grid%i1dcm(MAXINITREGIONS))
    allocate(grid%i2dcm(MAXINITREGIONS))
    allocate(grid%j1dcm(MAXINITREGIONS))
    allocate(grid%j2dcm(MAXINITREGIONS))
    allocate(grid%k1dcm(MAXINITREGIONS))
    allocate(grid%k2dcm(MAXINITREGIONS))
    allocate(grid%fracture_aperture(MAXINITREGIONS))
    allocate(grid%matrix_block(MAXINITREGIONS))
  endif
      
  allocate(grid%rock_density(MAXPERMREGIONS))
  allocate(grid%cpr(MAXPERMREGIONS))
  allocate(grid%dencpr(MAXPERMREGIONS))
  allocate(grid%ckdry(MAXPERMREGIONS))
  allocate(grid%ckwet(MAXPERMREGIONS))
  allocate(grid%tau(MAXPERMREGIONS))
  allocate(grid%cdiff(MAXPERMREGIONS))
  allocate(grid%cexp(MAXPERMREGIONS))

  allocate(grid%icaptype(MAXPERMREGIONS))
  if (grid%use_mph == PETSC_TRUE .or. grid%use_owg == PETSC_TRUE &
      .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
      .or. grid%use_richard == PETSC_TRUE) then
    allocate(grid%sir(1:grid%nphase,MAXPERMREGIONS))
  else
    allocate(grid%swir(MAXPERMREGIONS))
  endif
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
  if (grid%use_mph==PETSC_TRUE .or. grid%use_owg == PETSC_TRUE &
      .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE &
      .or. grid%use_richard == PETSC_TRUE) then
    allocate(grid%xxbc0(grid%ndof,MAXBCREGIONS))
    allocate(grid%iphasebc0(MAXBCREGIONS))
!    allocate(grid%dmax(0:grid%ndof-1))
    grid%xxbc0=0.D0
    grid%iphasebc0=0
  else
    allocate(grid%pressurebc0(grid%nphase, MAXBCREGIONS))
    allocate(grid%tempbc0(MAXBCREGIONS))
    allocate(grid%sgbc0(MAXBCREGIONS))
    allocate(grid%concbc0(MAXBCREGIONS))
    
      
    ! initialize
    grid%pressurebc0 = 0.d0
    grid%tempbc0 = 0.d0
    grid%concbc0 = 0.d0
    grid%sgbc0 = 0.d0
     
  endif
  grid%velocitybc0 = 0.d0
  !set scale factor for heat equation, i.e. use units of MJ for energy
  grid%scale = 1.d-6


  allocate(grid%rtot(grid%nlmax,2))
  !allocate(grid%rrtot(grid%nlmax,grid%ncomp))

  grid%rtot=0.D0
!  grid%qu_rate=0.D0

!-----------------------------------------------------------------------
  ! Compute arrays for indexing between local ghosted and non-ghosted 
  ! arrays.  I think the PETSc DA facilities may make these redundant,
  ! but since the ptran code uses them, I think it is easiest to just
  ! use these.
!-----------------------------------------------------------------------

  grid%istart = grid%nxs-grid%ngxs
  grid%jstart = grid%nys-grid%ngys
  grid%kstart = grid%nzs-grid%ngzs
  grid%iend = grid%istart+grid%nlx-1
  grid%jend = grid%jstart+grid%nly-1
  grid%kend = grid%kstart+grid%nlz-1

  ! Local <-> Ghosted Transformation
  grid%nG2L(:) = 0  ! Must initialize this to zero!
  !grid%nL2N(:) = 0
  n = 0
  do k=grid%kstart,grid%kend
    do j=grid%jstart,grid%jend
      do i=grid%istart,grid%iend
        n = n + 1
        ng = i+j*grid%ngx+k*grid%ngxy+1
        grid%nL2G(n) = ng
        grid%nG2L(ng) = n
      enddo
    enddo
  enddo

  do i=1,grid%ngmax
    j = grid%nG2L(i)
    if (j > 0) then
      k = grid%nL2G(j)
      if (i /= k) then
        print *,'Error in ghost-local node numbering for ghost node =', i
        print *,'node_id_gtol(i) =', j
        print *,'node_id_ltog(node_id_gtol(i)) =', k
        stop
      endif
    endif
  enddo
  ! Local(non ghosted)->Natural(natural order starts from 0)
  n=0
  do k=1,grid%nlz
    do j=1,grid%nly
      do i=1,grid%nlx
        n = n + 1
        na = i-1+grid%nxs+(j-1+grid%nys)*grid%nx+(k-1+grid%nzs)*grid%nxy
        if (na>(grid%nmax-1)) print *,'Wrong Nature order....'
        grid%nL2A(n) = na
        !print *,grid%myrank, k,j,i,n,na
        !grid%nG2N(ng) = na
      enddo
    enddo
  enddo
  print *,grid%myrank, grid%nxs,grid%ngxs,grid%nys,grid%ngys,grid%nzs,grid%ngzs
   
  call DAGetGlobalIndicesF90(grid%da_1_dof,grid%ngmax,grid%nG2N, ierr)

  pflowGrid_new = grid

end function pflowGrid_new


!#include "pflowgrid_destroy.F90"

subroutine pflowGrid_destroy(grid)
  
  implicit none
  
  type(pflowGrid), intent(inout) :: grid
  
  integer :: ierr

  ! Deallocate all of the arrays contained within grid.
  
  call DADestroy(grid%da_1_dof,ierr)
! call DADestroy(grid%da_3_dof,ierr)
  call DADestroy(grid%da_nphase_dof,ierr)
  call DADestroy(grid%da_3np_dof,ierr)
  call DADestroy(grid%da_ndof,ierr)

  call SNESDestroy(grid%snes,ierr)
  
! if (grid%use_numerical == PETSC_TRUE) then
!   call MatFDColoringDestroy(grid%matfdcoloring,ierr)
! endif
  
  call MatDestroy(grid%J,ierr)

  call VecDestroy(grid%porosity,ierr)
  call VecDestroy(grid%porosity0,ierr)
  call VecDestroy(grid%tor,ierr)
  call VecDestroy(grid%porosity_loc,ierr)
  call VecDestroy(grid%conc,ierr)
  call VecDestroy(grid%phis,ierr)
  call VecDestroy(grid%dx,ierr)
  call VecDestroy(grid%dy,ierr)
  call VecDestroy(grid%dz,ierr)
  call VecDestroy(grid%dx_loc,ierr)
  call VecDestroy(grid%dy_loc,ierr)
  call VecDestroy(grid%dz_loc,ierr)
  call VecDestroy(grid%volume,ierr)

  call VecDestroy(grid%ithrm,ierr)
  call VecDestroy(grid%ithrm_loc,ierr)
  call VecDestroy(grid%icap,ierr)
  call VecDestroy(grid%icap_loc,ierr)
  call VecDestroy(grid%iphas,ierr)
  call VecDestroy(grid%iphas_old,ierr)
  call VecDestroy(grid%iphas_loc,ierr)
  call VecDestroy(grid%temp,ierr)
  call VecDestroy(grid%ttemp,ierr)
  call VecDestroy(grid%ttemp_loc,ierr)
  call VecDestroy(grid%perm_xx,ierr)
  call VecDestroy(grid%perm_yy,ierr)
  call VecDestroy(grid%perm_zz,ierr)
  call VecDestroy(grid%perm0_xx,ierr)
  call VecDestroy(grid%perm0_yy,ierr)
  call VecDestroy(grid%perm0_zz,ierr)
  call VecDestroy(grid%perm_pow,ierr)
  call VecDestroy(grid%perm_xx_loc,ierr)
  call VecDestroy(grid%perm_yy_loc,ierr)
  call VecDestroy(grid%perm_zz_loc,ierr)
  call VecDestroy(grid%ppressure,ierr)
  call VecDestroy(grid%ppressure_loc,ierr)
  call VecDestroy(grid%ssat,ierr)
  call VecDestroy(grid%ssat_loc,ierr)
  call VecDestroy(grid%sat,ierr)
  call VecDestroy(grid%xxmol,ierr)
  call VecDestroy(grid%xxmol_loc,ierr)
  call VecDestroy(grid%xmol,ierr)
  call VecDestroy(grid%dp,ierr)
  call VecDestroy(grid%ddensity,ierr)
  call VecDestroy(grid%ddensity_loc,ierr)
  call VecDestroy(grid%density,ierr)
  call VecDestroy(grid%d_p,ierr)
  call VecDestroy(grid%d_t,ierr)
  call VecDestroy(grid%d_p_loc,ierr)
  call VecDestroy(grid%d_t_loc,ierr)
  call VecDestroy(grid%h,ierr)
  call VecDestroy(grid%hh,ierr)
  call VecDestroy(grid%hh_loc,ierr)
  call VecDestroy(grid%h_p,ierr)
  call VecDestroy(grid%h_t,ierr)
  call VecDestroy(grid%h_p_loc,ierr)
  call VecDestroy(grid%h_t_loc,ierr)
  call VecDestroy(grid%viscosity,ierr)
  call VecDestroy(grid%viscosity_loc,ierr)
  call VecDestroy(grid%v_p,ierr)
  call VecDestroy(grid%v_t,ierr)
  call VecDestroy(grid%v_p_loc,ierr)
  call VecDestroy(grid%v_t_loc,ierr)
  call VecDestroy(grid%vl,ierr)
  if (grid%using_pflowGrid == PETSC_TRUE) call VecDestroy(grid%vvl,ierr)
  call VecDestroy(grid%xx,ierr)
  call VecDestroy(grid%xx_loc,ierr)
  call VecDestroy(grid%yy,ierr)
  call VecDestroy(grid%dxx,ierr)
  call VecDestroy(grid%r,ierr)
  call VecDestroy(grid%accum,ierr)

  end subroutine pflowGrid_destroy

!======================================================================

!#include "pflowgrid_setup.F90"

! pflowGrid_setup():
! After the constructor function pflowGrid_new() has been used to 
! construct a new pflowGrid object, the grid's topology has been set 
! up but its geometry hasn't.  pflowGrid_setup() sets up the physical
! geometry of the grid, the initial values of the fields, and the 
! boundary conditions.
subroutine pflowGrid_setup(grid, inputfile)
  use hydrostat_module
  use water_eos_module
  use pflow_vector_ops_module
  use COND_module

  use LIQUID_module
  use TH_module
  use THC_module
  use TTPHASE_module
  use Flash_module
  use MPHASE_module
  use OWG_module
  use Vadose_module
  use Richard_module
  
  use utilities_module
  use readfield
  use Unstructured_Grid_module
  use pflow_solv_module  
                            
  implicit none
  
  type(pflowGrid), intent(inout) :: grid
  character(len=*), intent(in) :: inputfile
  
  integer :: ierr
  integer :: myrank
  integer*4 :: i, j, jn1, jn2, k, ird
  integer*4 :: mg1, mg2
  integer*4 :: m, n, ng
  integer*4 :: nc  ! Tracks number of connections computed.
  integer*4 :: ibc ! Used to index boundary condition blocks.
  integer*4 :: ir ! Used to index boundary condition regions.
  integer*4 :: ii1, ii2, jj1, jj2, kk1, kk2
    ! Used for indexing corners of boundary regions.
  PetscScalar, pointer :: dx_loc_p(:), dy_loc_p(:), dz_loc_p(:)
  real*8 alpha, maxstep, steptol
    ! Pointers used to access the arrays in grid%dx_loc et al.
  real*8, pointer :: volume_p(:)
    ! Used to access the contents of the local portion of grid%volume.
  ISColoring :: iscoloring
    ! Used for coloring the Jacobian so that we can take advantage of its
    ! sparsity when calculating it via finite differences.
  real*8, pointer :: den_p(:), &
                     pressure_p(:), temp_p(:), h_p(:), ran_p(:),&
                     phis_p(:), icap_p(:),ithrm_p(:), por_p(:),por0_p(:),&
                     tor_p(:),perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),&
                     perm_pow_p(:)
  real*8 :: d1,d2,val,val1,val2,val3,val4,por,dw_kg,random_nr,frand
  real*8 :: dl,hl
  integer*4 :: iseed,nx,ny,nz,na

  Vec :: temp0_nat_vec, temp1_nat_vec, temp2_nat_vec, temp3_nat_vec, &
          temp4_nat_vec !,temp5_nat_vec,temp6_nat_vec,temp7_nat_vec
  
  PetscViewer :: viewer
  Vec :: temp_vec
  
  ! Need to declare this function as external or else gfortran complains.
  ! Not sure why it complains and other compilers don't.
  external SNESDefaultComputeJacobianColor

#include "definitions.h"
  
!#ifdef DEBUG
! PetscViewer :: view_out
!#endif

  call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)

  !-----------------------------------------------------------------------
  ! Parse the input file to get dx, dy, dz, fields, etc. for each cell. 
  !-----------------------------------------------------------------------
  
  call pflowGrid_read_input(grid, inputfile)
  
 
! check number of dofs and phases
  if (grid%use_cond == PETSC_TRUE) then
    if (grid%ndof .ne. 1 .or. grid%nphase .ne. 1) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: COND ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  else if (grid%use_th == PETSC_TRUE) then
    if (grid%ndof .ne. 2 .or. grid%nphase .ne. 1) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: TH ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  else if (grid%use_thc == PETSC_TRUE) then
    if (grid%ndof .ne. 3 .or. grid%nphase .ne. 1) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: THC ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  else if (grid%use_2ph == PETSC_TRUE) then
    if (grid%ndof .ne. 4 .or. grid%nphase .ne. 2) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: 2PH ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  else if (grid%use_mph == PETSC_TRUE) then
    if (grid%ndof .ne. (grid%nspec+1)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: MPH ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  else if (grid%use_richard == PETSC_TRUE) then
    if (grid%ndof .ne. (grid%nspec+1)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: Richard ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  
    else if (grid%use_flash == PETSC_TRUE) then
    if (grid%ndof .ne. (grid%nspec+1)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: FLA ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  else if (grid%use_vadose == PETSC_TRUE) then
print *, grid%ndof, grid%nspec
    if (grid%ndof .ne. (grid%nspec+1)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: VAD ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  else if (grid%use_owg == PETSC_TRUE) then
    if (grid%ndof .ne. (grid%nspec)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: OWG ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  
  
  
  else
    if (grid%ndof .ne. 1 .or. grid%nphase .ne. 1) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: &
        &LIQUID ', &
        'ndof= ',grid%ndof,' nph= ',grid%nphase
      stop
    endif
  endif
  
! Calculate the x, y, z vectors that give the 
! physical coordinates of each cell.
  
  
  allocate(grid%x(grid%nmax))
  allocate(grid%y(grid%nmax))
  allocate(grid%z(grid%nmax))
 
  
  call pflowGrid_compute_xyz(grid)

  if (myrank == 0) then
    write(*,'(/,"++++++++++++++++++++++++++++++++++++++++++++++++++++&
      &++++++++")')
    write(*,'(" number of processors = ",i5,", npx,y,z= ",3i5)') &
      grid%commsize,grid%npx,grid%npy,grid%npz
    write(*,'(" number of dofs = ",i3,", number of phases = ",i3,i2)') &
      grid%ndof,grid%nphase
    if (grid%use_cond == PETSC_TRUE) then
      write(*,'(" mode = Conduction: T")')
    else if (grid%use_th == PETSC_TRUE) then
      write(*,'(" mode = TH: p, T")')
    else if (grid%use_thc == PETSC_TRUE) then
      write(*,'(" mode = THC: p, T, C")')
    else if (grid%use_2ph == PETSC_TRUE) then
      write(*,'(" mode = 2PH: p, T, s, C")')
    else if (grid%use_mph == PETSC_TRUE) then
      write(*,'(" mode = mPH: p, T, s/C")')
    else if (grid%use_flash == PETSC_TRUE) then
      write(*,'(" mode = flash: p, T, z")')
    else if (grid%use_vadose == PETSC_TRUE) then
      write(*,'(" mode = VAD: p, T, s/C")')
    else if (grid%use_richard == PETSC_TRUE) then
      write(*,'(" mode = Richard: p, T, s/C")')
    else if (grid%use_owg == PETSC_TRUE) then
      write(*,'(" mode = O+W+G: p, T, s/C")')
    else
      write(*,'(" mode = Single Liquid Phase: p")')
    endif
  endif

  !-----------------------------------------------------------------------
  ! Set up the Jacobian matrix.  We do this here instead of in 
  ! pflowgrid_new() because we may have to parse the input file to 
  ! determine how we want to do the Jacobian (matrix vs. matrix-free, for
  ! example).
  !-----------------------------------------------------------------------
  if (grid%use_analytical == PETSC_TRUE) then
  
    grid%ideriv = 1
  
!   if (myrank == 0) write(*,'(" analytical jacobian as ")'); &
!                    print *, grid%iblkfmt

    ! call DAGetMatrix(grid%da_ndof, MATMPIAIJ, grid%J, ierr)
    ! PETSc 2.1.6 introduces the MATAIJ matrix type.
    if (grid%iblkfmt == 0) then
      call DAGetMatrix(grid%da_ndof, MATAIJ, grid%J, ierr)
    else
      call DAGetMatrix(grid%da_ndof, MATBAIJ, grid%J, ierr)
    endif
   ! call  MatSetBlocksize(grid%J,grid%ndof,ierr)
    call MatSetOption(grid%J,MAT_COLUMN_ORIENTED,ierr)
    
    if (grid%use_cond == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, CondJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    else if (grid%use_th == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, THJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    else if (grid%use_thc == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, THCJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    else if (grid%use_2ph == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, TTPHASEJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    else if (grid%use_mph == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, MPHASEJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (grid%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
    else if (grid%use_flash == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, FlashJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (grid%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
    else if (grid%use_vadose == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, VADOSEJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (grid%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
    else if (grid%use_richard == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, RichardJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (grid%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
    else if (grid%use_owg == PETSC_TRUE) then
      call SNESSetJacobian(grid%snes, grid%J, grid%J, OWGJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (grid%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
    else
      call SNESSetJacobian(grid%snes, grid%J, grid%J, LiquidJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    endif

  else if (grid%use_matrix_free == PETSC_TRUE) then
  
    grid%ideriv = 0
  
    if (myrank == 0) write(*,'(" Using matrix-free Newton-Krylov")')
    
    if (grid%use_cond == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%ttemp, grid%J, ierr)
    else if (grid%use_th == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%xx, grid%J, ierr)
    else if (grid%use_thc == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%xx, grid%J, ierr)
    else if (grid%use_2ph == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%xx, grid%J, ierr)
    else if (grid%use_mph == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%xx, grid%J, ierr)
    else if (grid%use_flash == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%xx, grid%J, ierr)
    else if (grid%use_richard == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%xx, grid%J, ierr)
    else if (grid%use_vadose == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%xx, grid%J, ierr)
    else if (grid%use_owg == PETSC_TRUE) then
      call MatCreateMFFD(grid%snes, grid%xx, grid%J, ierr)
    else
      call MatCreateMFFD(grid%snes, grid%ppressure, grid%J, ierr)
    endif
    
    ! It seems that I ought to call SNESSetJacobian here now, but I don't know
    ! what function I am supposed to pass to it to get it to use one of the 
    ! finite-difference routines for computing the Jacobian.  It might not
    ! actually matter if -snes_mf has been specified.
    ! Pernice thinks that perhaps the I need to provide a function which 
    ! simply calls MatAssemblyBegin/End.
    call SNESSetJacobian(grid%snes, grid%J, grid%J, &
                         ComputeMFJacobian, PETSC_NULL_OBJECT, ierr)

    ! Use "Walker-Pernice" differencing.
    call MatMFFDSetType(grid%J, MATMFFD_WP, ierr)

    if (grid%print_hhistory == PETSC_TRUE) then
      allocate(grid%hhistory(HHISTORY_LENGTH))
      call MatMFFDSetHHistory(grid%J, grid%hhistory, HHISTORY_LENGTH, ierr)
    endif

    if (grid%monitor_h == PETSC_TRUE) then
      call SNESMonitorSet(grid%snes, pflowgrid_MonitorH, grid, &
                          PETSC_NULL_OBJECT, ierr)
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
        
    call DAGetColoring(grid%da_ndof, IS_COLORING_GLOBAL, iscoloring, ierr)
    
    call MatFDColoringCreate(grid%J, iscoloring, grid%matfdcoloring, ierr)
    
    call ISColoringDestroy(iscoloring, ierr)

    if (grid%use_cond == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, CondResidual, &
                                      grid, ierr)
    else if (grid%use_th == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, THResidual, &
                                      grid, ierr)
    else if (grid%use_thc == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, THCResidual, &
                                      grid, ierr)
    else if (grid%use_2ph == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, TTPHASEResidual, &
                                      grid, ierr)
    else if (grid%use_mph == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, MPHASEResidual, &
                                      grid, ierr)
    else if (grid%use_flash == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, FlashResidual, &
                                      grid, ierr)
    else if (grid%use_vadose == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, VADOSEResidual, &
                                      grid, ierr)
    else if (grid%use_richard == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, RichardResidual, &
                                      grid, ierr)
    else if (grid%use_owg == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, OWGResidual, &
                                      grid, ierr)
    else
      call MatFDColoringSetFunctionSNES(grid%matfdcoloring, LiquidResidual, &
                                      grid, ierr)
    endif
    
    call MatFDColoringSetFromOptions(grid%matfdcoloring, ierr)
    
    call SNESSetJacobian(grid%snes, grid%J, grid%J, &
                         SNESDefaultComputeJacobianColor,  &
                         grid%matfdcoloring, ierr)
  endif

  if (myrank == 0) write(*,'("++++++++++++++++++++++++++++++++&
                     &++++++++++++++++++++++++++++",/)')

  if (grid%use_cond == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, CondResidual, grid, ierr)
  else if (grid%use_th == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, THResidual, grid, ierr)
  else if (grid%use_thc == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, THCResidual, grid, ierr)
  else if (grid%use_2ph == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, TTPHASEResidual, grid, ierr)
  else if (grid%use_mph == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, MPHASEResidual, grid, ierr)
  else if (grid%use_flash == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, FlashResidual, grid, ierr)
  else if (grid%use_vadose == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, VADOSEResidual, grid, ierr)
  else if (grid%use_Richard == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, RichardResidual, grid, ierr)
  else if (grid%use_owg == PETSC_TRUE) then
    call SNESSetFunction(grid%snes, grid%r, OWGResidual, grid, ierr)
  else
    call SNESSetFunction(grid%snes, grid%r, LiquidResidual, grid, ierr)
  endif

  ! Set the tolerances for the Newton solver.
  call SNESSetTolerances(grid%snes, grid%atol, grid%rtol, grid%stol, & 
                         grid%maxit, grid%maxf, ierr)

 

  call SNESSetFromOptions(grid%snes, ierr)
  if (myrank == 0) write(*,'("  Finished setting up of SNES 1")')
  
  call SNESLineSearchGetParams(grid%snes, alpha, maxstep, steptol, ierr) 
  if (myrank == 0) write(*,'("  Finished setting up of SNES 2")')
  call SNESLineSearchSetParams(grid%snes, alpha, maxstep, grid%stol, ierr) 
  if (myrank == 0) write(*,'("  Finished setting up of SNES 3")')

  call SNESGetKSP(grid%snes, grid%ksp, ierr)
  call KSPSetTolerances(grid%ksp,grid%rtol,grid%atol,grid%dtol, &
      10000,ierr)

  

 if (myrank == 0) write(*,'("  Finished setting up of SNES ")')

  !-----------------------------------------------------------------------
  ! Set up cell topology and geometry of interior connections, and 
  ! calculate interior interface areas and cell volumes.
  !-----------------------------------------------------------------------
  
  call DACreateNaturalVector(grid%da_1_dof,temp1_nat_vec,ierr)
  call VecDuplicate(temp1_nat_vec, temp2_nat_vec, ierr)
  call VecDuplicate(temp1_nat_vec, temp3_nat_vec, ierr)
  if (myrank == 0) then
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
  call VecGetArrayF90(grid%dx_loc, dx_loc_p, ierr)
  call VecGetArrayF90(grid%dy_loc, dy_loc_p, ierr)
  call VecGetArrayF90(grid%dz_loc, dz_loc_p, ierr)

  nc = 0
  grid%nconnx = 0
  grid%nconny = 0

  if (grid%igeom == 2) then
    allocate(grid%rd(0:grid%nx))
    grid%rd = 0.D0
    grid%rd(0) = grid%Radius_0 
    do i = 1, grid%nx
      grid%rd(i) = grid%rd(i-1) + grid%dx0(i)
    enddo
  endif    
  
! print *,'setup-RD:: ',grid%myrank,grid%rd, grid%dx0
  
  ! x-connections
  if (grid%ngx > 1) then
    do k = grid%kstart, grid%kend
      do j = grid%jstart, grid%jend
        do i = 1, grid%ngx - 1
          mg1 = i + j * grid%ngx + k * grid%ngxy
          mg2 = mg1 + 1
          nc = nc + 1
          grid%nd1(nc) = mg1
          grid%nd2(nc) = mg2
          grid%dist1(nc) = 0.5d0 * dx_loc_p(mg1)
          grid%dist2(nc) = 0.5d0 * dx_loc_p(mg2)
          grid%delz(nc) = 0.d0
          grid%grav_ang(nc)=0.D0
        !  abs_x = rd(i-1+grid%nxs)
          if (grid%igeom == 1) then
            grid%area(nc) = dy_loc_p(mg1) * dz_loc_p(mg1)
          else if (grid%igeom == 2) then
            grid%area(nc) = 2.D0 * Pi * grid%rd(i+grid%nxs) * dz_loc_p(mg1)
            print *,'area nc: ',grid%myrank,nc,i,grid%rd(i+grid%nxs)
          else if (grid%igeom == 3) then
            grid%area(nc) = 4.D0 * Pi * grid%rd(i+grid%nxs)**2
          endif
          grid%iperm1(nc) = 1
          grid%iperm2(nc) = 1
        enddo
      enddo
    enddo
    grid%nconnx = nc
  endif

  ! y-connections
  if (grid%ngy > 1) then
    do k = grid%kstart, grid%kend
      do i = grid%istart, grid%iend
        do j = 1, grid%ngy - 1
          mg1 = i + 1 + (j-1) * grid%ngx + k * grid%ngxy
          mg2 = mg1 + grid%ngx
          nc = nc + 1
          grid%nd1(nc) = mg1
          grid%nd2(nc) = mg2
          grid%dist1(nc) = 0.5d0 * dy_loc_p(mg1)
          grid%dist2(nc) = 0.5d0 * dy_loc_p(mg2)
          grid%delz(nc) = 0.d0
          grid%grav_ang(nc)=0.D0
          grid%area(nc) = dx_loc_p(mg1) * dz_loc_p(mg1)
          grid%iperm1(nc) = 2
          grid%iperm2(nc) = 2
        enddo
      enddo
    enddo
    grid%nconny = nc
  endif
      
  ! z-connections
  if (grid%ngz > 1) then
    do j = grid%jstart, grid%jend
      do i = grid%istart, grid%iend
        do k = 1, grid%ngz - 1
          mg1 = i + 1 + j * grid%ngx + (k-1) * grid%ngxy
          mg2 = mg1 + grid%ngxy
          nc = nc + 1
          grid%nd1(nc) = mg1
          grid%nd2(nc) = mg2
          d1 = 0.5d0 * dz_loc_p(mg1)
          d2 = 0.5d0 * dz_loc_p(mg2)
          grid%dist1(nc) = d1
          grid%dist2(nc) = d2
          grid%delz(nc) = d1 + d2
          grid%grav_ang(nc)=1.D0
          if (grid%igeom == 1) then
            grid%area(nc) = dx_loc_p(mg1) * dy_loc_p(mg1)
          else if (grid%igeom == 2) then
            grid%area(nc) =  Pi * (grid%rd(i+grid%nxs)+  &
                             grid%rd(i-1+grid%nxs))* &
                             (grid%rd(i+grid%nxs) - grid%rd(i-1+grid%nxs))  
            print *, 'area nc ',grid%myrank, nc, i,  grid%area(nc)
          endif
          grid%iperm1(nc) = 3
          grid%iperm2(nc) = 3
        enddo
      enddo
    enddo
  endif
  
  ! Calculate cell volumes for local cells.
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  do n=1, grid%nlmax
    ng = grid%nL2G(n)
    if (grid%igeom == 1) then
      volume_p(n) = dx_loc_p(ng) * dy_loc_p(ng) * dz_loc_p(ng)
    else if (grid%igeom == 2) then
      i = mod(mod((n),grid%nlxy),grid%nlx)!+(grid%ngxs-grid%nxs)
      if (i==0) i = grid%nlx
      volume_p(n) = Pi * (grid%rd(i+grid%nxs) + grid%rd(i-1+grid%nxs))*&
      (grid%rd(i+grid%nxs) - grid%rd(i-1+grid%nxs)) * dz_loc_p(ng)
      print *, 'setup: Vol ', grid%myrank, n,i, grid%rd(i+grid%nxs),volume_p(n)
    else if (grid%igeom == 3) then
    endif
  enddo
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  
  write(*,'(" myrank= ",i3,", nlmax= ",i6,", nlx,y,z= ",3i4, &
    & ", nxs,e = ",2i4,", nys,e = ",2i4,", nzs,e = ",2i4)') &
    myrank,grid%nlmax,grid%nlx,grid%nly,grid%nlz, &
    grid%nxs,grid%nxe,grid%nys,grid%nye,grid%nzs,grid%nze

  write(*,'(" myrank= ",i3,", ngmax= ",i6,", ngx,y,z= ",3i4, &
    & ", ngxs,e= ",2i4,", ngys,e= ",2i4,", ngzs,e= ",2i4)') &
    myrank,grid%ngmax,grid%ngx,grid%ngy,grid%ngz, &
    grid%ngxs,grid%ngxe,grid%ngys,grid%ngye,grid%ngzs,grid%ngze

  !-----------------------------------------------------------------------
  ! Set up boundary connections.
  !-----------------------------------------------------------------------
      
  grid%nconnbc = 0
  if (grid%nx > 1 .and. grid%ny == 1 .and. grid%nz == 1) then
    if (grid%nxs == grid%ngxs) grid%nconnbc = grid%nconnbc + 1
    if (grid%nxe == grid%ngxe) grid%nconnbc = grid%nconnbc + 1
  else if (grid%nx == 1 .and. grid%ny == 1 .and. grid%nz > 1) then
    if (grid%nzs == grid%ngzs) grid%nconnbc = grid%nconnbc + 1
    if (grid%nze == grid%ngze) grid%nconnbc = grid%nconnbc + 1
  else
    if (grid%nx > 1) then
      if (grid%nxs == grid%ngxs) grid%nconnbc = grid%nconnbc + grid%nlyz
      if (grid%nxe == grid%ngxe) grid%nconnbc = grid%nconnbc + grid%nlyz
    endif
    if (grid%ny > 1) then
      if (grid%nys == grid%ngys) grid%nconnbc = grid%nconnbc + grid%nlxz
      if (grid%nye == grid%ngye) grid%nconnbc = grid%nconnbc + grid%nlxz
    endif
    if (grid%nz > 1) then
      if (grid%nzs == grid%ngzs) grid%nconnbc = grid%nconnbc + grid%nlxy
      if (grid%nze == grid%ngze) grid%nconnbc = grid%nconnbc + grid%nlxy
    endif
  endif
      
  write(*,'(" --> pflowconn: rank = ",i4, &
       &", boundary connections =", i6)') myrank,grid%nconnbc
  
! set initial conditions by region for pressure, temperature, saturation
! and concentration

  if (grid%use_mph == PETSC_TRUE) then
    call pflow_mphase_setupini(grid)
   else if (grid%use_flash == PETSC_TRUE) then
    call pflow_flash_setupini(grid)
   else if (grid%use_richard == PETSC_TRUE) then
    call pflow_richard_setupini(grid)
   else if (grid%use_owg == PETSC_TRUE) then
    call pflow_owg_setupini(grid)
   else if (grid%use_vadose == PETSC_TRUE) then
    call pflow_vadose_setupini(grid)
   else 
    call DACreateNaturalVector(grid%da_nphase_dof,temp1_nat_vec,ierr)
    call VecDuplicate(temp1_nat_vec, temp4_nat_vec, ierr)
    call DACreateNaturalVector(grid%da_1_dof,temp2_nat_vec,ierr)
    if  (grid%ndof == 3) call VecDuplicate(temp2_nat_vec, temp3_nat_vec, ierr)
    if  (grid%ndof == 4) call VecDuplicate(temp1_nat_vec, temp3_nat_vec, ierr)
 
    if (myrank == 0) then
      do ir = 1,grid%iregini
        do k = grid%k1ini(ir),grid%k2ini(ir)
          do j = grid%j1ini(ir),grid%j2ini(ir)
            do i = grid%i1ini(ir),grid%i2ini(ir)
              n = i+(j-1)*grid%nx+(k-1)*grid%nxy-1
              jn1 = 1+n*grid%nphase-1
              jn2 = 2+n*grid%nphase-1
              
              val = grid%pres_ini(ir)
              call VecSetValue(temp1_nat_vec,jn1,val,INSERT_VALUES,ierr)
              if (grid%nphase>1) &
                call VecSetValue(temp1_nat_vec,jn2,val,INSERT_VALUES,ierr)
  
              val = grid%temp_ini(ir)
              call VecSetValue(temp2_nat_vec,n,val,INSERT_VALUES,ierr)
             
              val = grid%sat_ini(ir)
              if (grid%nphase == 1) then
                call VecSetValue(temp4_nat_vec,jn1,val,INSERT_VALUES,ierr)
              else
                ! need to decide if sg or sl is to be read in-use sl for now!
                call VecSetValue(temp4_nat_vec,jn1,val,INSERT_VALUES,ierr)
                val=1.D0-val
                call VecSetValue(temp4_nat_vec,jn2,val,INSERT_VALUES,ierr)
              endif
              
              if (grid%ndof == 3) then
                val = grid%conc_ini(ir)
 !              print *,'pflowgrid_setup: ',n+1,i,j,k,ir,val
                call VecSetValue(temp3_nat_vec,n,val,INSERT_VALUES,ierr)
              endif
 
              if (grid%ndof == 4) then
                val = grid%conc_ini(ir)
 !              print *,'pflowgrid_setup: ',n+1,i,j,k,ir,val
                call VecSetValue(temp3_nat_vec,jn2,val,INSERT_VALUES,ierr)
              endif
            enddo
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
    call VecAssemblyBegin(temp4_nat_vec,ierr)
    call VecAssemblyEnd(temp4_nat_vec,ierr)
    call DANaturalToGlobalBegin(grid%da_nphase_dof,temp1_nat_vec, &
                                INSERT_VALUES,grid%pressure,ierr)
    call DANaturalToGlobalEnd(grid%da_nphase_dof,temp1_nat_vec,INSERT_VALUES, &
                              grid%pressure,ierr)
    call DANaturalToGlobalBegin(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
                                grid%temp,ierr)
    call DANaturalToGlobalEnd(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
                              grid%temp,ierr)
    call DANaturalToGlobalBegin(grid%da_nphase_dof,temp4_nat_vec, &
                                INSERT_VALUES,grid%sat,ierr)
    call DANaturalToGlobalEnd(grid%da_nphase_dof,temp4_nat_vec,INSERT_VALUES, &
                              grid%sat,ierr)
!    call VecView(temp1_nat_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
!    call VecView(temp2_nat_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
!    call VecView(temp4_nat_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)

    call VecDestroy(temp1_nat_vec,ierr)
    call VecDestroy(temp2_nat_vec,ierr)
    call VecDestroy(temp4_nat_vec,ierr)
  
    if (grid%ndof == 3) then
      call VecAssemblyBegin(temp3_nat_vec,ierr)
      call VecAssemblyEnd(temp3_nat_vec,ierr)
      call DANaturalToGlobalBegin(grid%da_1_dof,temp3_nat_vec,INSERT_VALUES, &
                                  grid%conc,ierr)
      call DANaturalToGlobalEnd(grid%da_1_dof,temp3_nat_vec,INSERT_VALUES, &
                                grid%conc,ierr)
      call VecDestroy(temp3_nat_vec,ierr)
    
!     call VecView(grid%conc,PETSC_VIEWER_STDOUT_WORLD,ierr)
    endif
 
   if (grid%ndof == 4) then
      call VecAssemblyBegin(temp3_nat_vec,ierr)
      call VecAssemblyEnd(temp3_nat_vec,ierr)
 !     call VecView(temp3_nat_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
      call DANaturalToGlobalBegin(grid%da_nphase_dof,temp3_nat_vec, &
                                  INSERT_VALUES,grid%xmol,ierr)
      call DANaturalToGlobalEnd(grid%da_nphase_dof,temp3_nat_vec, &
                                INSERT_VALUES,grid%xmol,ierr)
      call VecDestroy(temp3_nat_vec,ierr)
    
!     call VecView(grid%conc,PETSC_VIEWER_STDOUT_WORLD,ierr)
    endif
 
  endif 
 
 
  ! set hydrostatic properties for initial and boundary conditions with depth
  if (grid%ihydrostatic == 1) then
    if (grid%use_mph == PETSC_TRUE .or. grid%use_vadose == PETSC_TRUE &
     .or. grid%use_flash == PETSC_TRUE .or. grid%use_richard == PETSC_TRUE) then
      print *,'in hydro'
      call mhydrostatic(grid)
    elseif (grid%use_owg == PETSC_TRUE) then
      print *,'in hydro'
      call owghydrostatic(grid)
    else
      call hydrostatic(grid)
    endif
  endif
  
  if (grid%iread_init==2) call Read_init_field(grid)
  
  print *,'Finished Hydro'
  if (grid%use_mph == PETSC_TRUE .or. grid%use_owg == PETSC_TRUE &
      .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
       .or. grid%use_richard == PETSC_TRUE ) then
  !  grid%pressurebc0(2,:) = grid%pressurebc0(1,:)
    grid%velocitybc0(2,:) = grid%velocitybc0(1,:)
  endif
 
  if (grid%use_2ph == PETSC_TRUE) then
    grid%pressurebc0(2,:) = grid%pressurebc0(1,:)
    grid%velocitybc0(2,:) = grid%velocitybc0(1,:)
  endif

!  call VecView(grid%conc,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
!  call VecView(grid%pressure,PETSC_VIEWER_STDOUT_WORLD,ierr)
!  call VecView(grid%temp,PETSC_VIEWER_STDOUT_WORLD,ierr)
!  call VecView(grid%xmol,PETSC_VIEWER_STDOUT_WORLD,ierr)
!  call VecView(grid%sat,PETSC_VIEWER_STDOUT_WORLD,ierr)

  deallocate(grid%k1ini)
  deallocate(grid%k2ini)
  deallocate(grid%j1ini)
  deallocate(grid%j2ini)
  deallocate(grid%i1ini)
  deallocate(grid%i2ini)
  
  if (grid%use_mph == PETSC_TRUE .or. grid%use_owg==PETSC_TRUE &
      .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
      .or. grid%use_richard == PETSC_TRUE) then
    deallocate(grid%xx_ini)
    deallocate(grid%iphas_ini)
  else 
    deallocate(grid%pres_ini)
    deallocate(grid%temp_ini)
    deallocate(grid%sat_ini)
    deallocate(grid%xmol_ini)
    deallocate(grid%conc_ini)
  endif
  print *,'deallocate ini'

!************End of initial Condition Setup ***********************


!geh
  if (grid%iread_geom > -1) then

    if (grid%nconnbc > 0) then
      allocate(grid%mblkbc(grid%nconnbc))
      allocate(grid%ibconn(grid%nconnbc))
      allocate(grid%distbc(grid%nconnbc))
      allocate(grid%areabc(grid%nconnbc))
      allocate(grid%ipermbc(grid%nconnbc))
      allocate(grid%delzbc(grid%nconnbc))
    
!    allocate(grid%velocitybc(grid%nphase,grid%nconnbc))
    
      allocate(grid%vlbc(grid%nconnbc))
      allocate(grid%vvlbc(grid%nconnbc))
      allocate(grid%vgbc(grid%nconnbc))
      allocate(grid%vvgbc(grid%nconnbc))
  
      grid%vlbc=0.D0
      grid%vgbc=0.D0
    endif
!geh
  endif

  if (grid%use_mph == PETSC_TRUE .or. grid%use_owg == PETSC_TRUE &
      .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
      .or. grid%use_richard == PETSC_TRUE) then
    allocate(grid%varbc(1:(grid%ndof+1)*(2+7*grid%nphase + 2 *  &
                                       grid%nphase*grid%nspec)))
  else  
    allocate(grid%density_bc(grid%nphase))
    allocate(grid%d_p_bc(grid%nphase))
    allocate(grid%d_t_bc(grid%nphase))
    allocate(grid%d_c_bc(grid%nphase))
    allocate(grid%d_s_bc(grid%nphase))
    allocate(grid%avgmw_bc(grid%nphase))
    allocate(grid%avgmw_c_bc(grid%nphase*grid%npricomp))
    allocate(grid%hh_bc(grid%nphase))
    allocate(grid%h_p_bc(grid%nphase))
    allocate(grid%h_t_bc(grid%nphase))
    allocate(grid%h_c_bc(grid%nphase*grid%npricomp))
    allocate(grid%h_s_bc(grid%nphase))
    allocate(grid%uu_bc(grid%nphase))
    allocate(grid%u_p_bc(grid%nphase))
    allocate(grid%u_t_bc(grid%nphase))
    allocate(grid%u_c_bc(grid%nphase*grid%npricomp))
    allocate(grid%u_s_bc(grid%nphase))
    allocate(grid%df_bc(grid%nphase*grid%nspec))
    allocate(grid%df_p_bc(grid%nphase*grid%nspec))
    allocate(grid%df_t_bc(grid%nphase*grid%nspec))
    allocate(grid%df_c_bc(grid%nphase*grid%nspec*grid%npricomp))
    allocate(grid%df_s_bc(grid%nphase*grid%nspec))
    allocate(grid%hen_bc(grid%nphase*grid%nspec))
    allocate(grid%hen_p_bc(grid%nphase*grid%nspec))
    allocate(grid%hen_t_bc(grid%nphase*grid%nspec))
    allocate(grid%hen_c_bc(grid%nphase*grid%nspec*grid%npricomp))
    allocate(grid%hen_s_bc(grid%nphase*grid%nspec))
    allocate(grid%viscosity_bc(grid%nphase))
    allocate(grid%v_p_bc(grid%nphase))
    allocate(grid%v_t_bc(grid%nphase))
    allocate(grid%pc_bc(grid%nphase))
    allocate(grid%pc_p_bc(grid%nphase))
    allocate(grid%pc_t_bc(grid%nphase))
    allocate(grid%pc_c_bc(grid%nphase*grid%npricomp))
    allocate(grid%pc_s_bc(grid%nphase))
    allocate(grid%kvr_bc(grid%nphase))
    allocate(grid%kvr_p_bc(grid%nphase))
    allocate(grid%kvr_t_bc(grid%nphase))
    allocate(grid%kvr_c_bc(grid%nphase*grid%npricomp))
    allocate(grid%kvr_s_bc(grid%nphase))
  endif

!geh
  if (grid%iread_geom > -1) then

    nc = 0 
    if (grid%nxs == grid%ngxs .or. grid%nxe == grid%ngxe &
        .or. grid%nys == grid%ngys .or. grid%nye == grid%ngye &
        .or. grid%nzs == grid%ngzs .or. grid%nze == grid%ngze) then

    ! calculate boundary conditions locally on only those processors which 
    ! contain a boundary!

      do ibc = 1, grid%nblkbc
        do ir = grid%iregbc1(ibc), grid%iregbc2(ibc)
          kk1 = grid%k1bc(ir) - grid%nzs
          kk2 = grid%k2bc(ir) - grid%nzs
          jj1 = grid%j1bc(ir) - grid%nys
          jj2 = grid%j2bc(ir) - grid%nys
          ii1 = grid%i1bc(ir) - grid%nxs
          ii2 = grid%i2bc(ir) - grid%nxs

          kk1 = max(1,kk1)
          kk2 = min(grid%nlz,kk2)
          jj1 = max(1,jj1)
          jj2 = min(grid%nly,jj2)
          ii1 = max(1,ii1)
          ii2 = min(grid%nlx,ii2)

          if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle 

          do k = kk1,kk2
            do j = jj1,jj2
              do i = ii1,ii2
                nc = nc + 1
                m = i+(j-1)*grid%nlx+(k-1)*grid%nlxy
                grid%mblkbc(nc) = m  ! m is a local index
                grid%ibconn(nc) = ibc
                ng = grid%nL2G(m)
              ! Use ghosted index to access dx, dy, dz because we have
              ! already done a global-to-local scatter for computing the
              ! interior node connections.
        
!               print *,'pflowgrid_mod: ',nc,ibc,ir,m,ng,ii1,ii2,kk1,kk2, &
!                        grid%nblkbc,grid%igeom
        
                select case(grid%igeom)
                  case(1) ! cartesian
                    if (grid%iface(ibc) == 1) then
                      grid%distbc(nc) = 0.5d0*dx_loc_p(ng)
                      grid%areabc(nc) = dy_loc_p(ng)*dz_loc_p(ng)
                      grid%ipermbc(nc) = 1
                      grid%delzbc(nc) = 0.d0
                    else if (grid%iface(ibc) == 2) then
                      grid%distbc(nc) = 0.5d0*dx_loc_p(ng)
                      grid%areabc(nc) = dy_loc_p(ng)*dz_loc_p(ng)
                      grid%ipermbc(nc) = 1
                      grid%delzbc(nc) = 0.d0
                    else if (grid%iface(ibc) == 3) then
                      grid%distbc(nc) = 0.5d0*dz_loc_p(ng)
                      grid%areabc(nc) = dx_loc_p(ng)*dy_loc_p(ng)
                      grid%ipermbc(nc) = 2
                      grid%delzbc(nc) = grid%distbc(nc)
                    else if (grid%iface(ibc) == 4) then
                      grid%distbc(nc) = 0.5d0*dz_loc_p(ng)
                      grid%areabc(nc) = dx_loc_p(ng)*dy_loc_p(ng)
                      grid%ipermbc(nc) = 2
                      grid%delzbc(nc) = -grid%distbc(nc)
                    else if (grid%iface(ibc) == 5) then
                      grid%distbc(nc) = 0.5d0*dy_loc_p(ng)
                      grid%areabc(nc) = dx_loc_p(ng)*dz_loc_p(ng)
                      grid%ipermbc(nc) = 3
                      grid%delzbc(nc) = 0.d0
                    else if (grid%iface(ibc) == 6) then
                      grid%distbc(nc) = 0.5d0*dy_loc_p(ng)
                      grid%areabc(nc) = dx_loc_p(ng)*dz_loc_p(ng)
                      grid%ipermbc(nc) = 3
                      grid%delzbc(nc) = 0.d0
                    endif
                  case(2) ! cylindrical
                    ird = mod(mod((m),grid%nlxy),grid%nlx) + grid%nxs 
                    if (grid%iface(ibc) == 1) then
                      grid%distbc(nc) = 0.5d0*dx_loc_p(ng)
                      grid%areabc(nc) = 2.0D0*Pi*grid%rd(ird-1)*dz_loc_p(ng)
                      grid%ipermbc(nc) = 1
                      grid%delzbc(nc) = 0.d0
                    else if (grid%iface(ibc) == 2) then
                      grid%distbc(nc) = 0.5d0*dx_loc_p(ng)
                      grid%areabc(nc) = 2.0D0*Pi*grid%rd(ird)*dz_loc_p(ng)
                      grid%ipermbc(nc) = 1
                      grid%delzbc(nc) = 0.d0
                    else if (grid%iface(ibc) == 3) then
                      grid%distbc(nc) = 0.5d0*dz_loc_p(ng)
                      grid%areabc(nc) = Pi*(grid%rd(ird) + grid%rd(ird-1))* &
                                        (grid%rd(ird) - grid%rd(ird-1))  
                      grid%ipermbc(nc) = 2
                      grid%delzbc(nc) = grid%distbc(nc)
                    else if (grid%iface(ibc) == 4) then
                      grid%distbc(nc) = 0.5d0*dz_loc_p(ng)
                      grid%areabc(nc) = Pi*(grid%rd(ird) + grid%rd(ird-1))* &
                                        (grid%rd(ird) - grid%rd(ird-1))
                      grid%ipermbc(nc) = 2
                      grid%delzbc(nc) = -grid%distbc(nc)
                    endif  
                  case(3) ! spherical
                end select
              enddo ! i
            enddo ! j
          enddo ! k
        enddo ! ir
      enddo ! ibc
    endif
    call VecRestoreArrayF90(grid%dx_loc, dx_loc_p, ierr)
    call VecRestoreArrayF90(grid%dy_loc, dy_loc_p, ierr)
    call VecRestoreArrayF90(grid%dz_loc, dz_loc_p, ierr)

    if (grid%nconnbc .ne. nc) then
      write(*,*) 'Error in computing boundary connections: ', &
        'rank = ',myrank,' nconnbc = ',grid%nconnbc,' nc = ',nc
      stop
    endif
   
    do nc=1,grid%nconnbc
      print *, 'BC:', nc, grid%areabc(nc),grid%distbc(nc),  grid%mblkbc(nc) 
    enddo

  !-----------------------------------------------------------------------
  ! Set up boundary conditions at interfaces
  !-----------------------------------------------------------------------
    allocate(grid%velocitybc(grid%nphase, grid%nconnbc))
    if (grid%use_mph == PETSC_TRUE .or. grid%use_owg==PETSC_TRUE &
        .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
        .or. grid%use_richard == PETSC_TRUE) then
      allocate(grid%xxbc(grid%ndof,grid%nconnbc))
      allocate(grid%iphasebc(grid%nconnbc))
      allocate(grid%xphi_co2_bc(grid%nconnbc))
      allocate(grid%xxphi_co2_bc(grid%nconnbc))

      do nc = 1, grid%nconnbc
        ibc = grid%ibconn(nc)
        grid%xxbc(:,nc)=grid%xxbc0(:,ibc)
        grid%iphasebc(nc)=grid%iphasebc0(ibc)
        grid%velocitybc(:,nc) = grid%velocitybc0(:,ibc)
      enddo
      if (grid%using_pflowGrid == PETSC_FALSE) then
        deallocate(grid%xxbc0)
        deallocate(grid%iphasebc0)
      endif
      if (grid%iread_init==2) call Boundary_adjustment(grid)
  
    else
   !   allocate(grid%velocitybc(grid%nphase, grid%nconnbc))
      allocate(grid%pressurebc(grid%nphase, grid%nconnbc))
      allocate(grid%tempbc(grid%nconnbc))
      allocate(grid%sgbc(grid%nconnbc))
      allocate(grid%concbc(grid%nconnbc))
      allocate(grid%xphi_co2_bc(grid%nconnbc))
      allocate(grid%xxphi_co2_bc(grid%nconnbc))
      
    ! initialize
      grid%pressurebc = 0.d0
      grid%tempbc = 0.d0
      grid%concbc = 0.d0
      grid%sgbc = 0.d0
      grid%velocitybc = 0.d0
      
      do nc = 1, grid%nconnbc
        ibc = grid%ibconn(nc)
        grid%pressurebc(:,nc) = grid%pressurebc0(:,ibc)
        grid%tempbc(nc)       = grid%tempbc0(ibc)
        grid%concbc(nc)       = grid%concbc0(ibc)
        grid%sgbc(nc)         = grid%sgbc0(ibc)
        grid%velocitybc(:,nc) = grid%velocitybc0(:,ibc)
        
!       print *,'pflowgrid_mod: bc- ',nc,ibc,grid%pressurebc(1,nc), &
!                grid%pressurebc(2,nc), &
!                grid%tempbc(nc),grid%sgbc(nc),grid%concbc(nc), &
!                grid%velocitybc(1,nc),grid%velocitybc(2,nc)
      enddo
     
      deallocate(grid%pressurebc0)
      deallocate(grid%tempbc0)
      deallocate(grid%concbc0)
      deallocate(grid%sgbc0)
    endif
    deallocate(grid%velocitybc0)

    if (myrank == 0) write(*,'("  Finished setting up of Geometry ")')

!geh
  endif
  
  !-----------------------------------------------------------------------
  ! Set up the transformation from physical coordinates
  ! to local domains.
  !-----------------------------------------------------------------------

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
  
  if (myrank == 0) then 

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

  call VecGetArrayF90(grid%ttemp,ran_p,ierr)
  call VecGetArrayF90(grid%icap,icap_p,ierr)
  call VecGetArrayF90(grid%ithrm,ithrm_p,ierr)
  call VecGetArrayF90(grid%porosity,por_p,ierr)
  call VecGetArrayF90(grid%porosity0,por0_p,ierr)
  call VecGetArrayF90(grid%perm_xx,perm_xx_p,ierr)
  call VecGetArrayF90(grid%perm_yy,perm_yy_p,ierr)
  call VecGetArrayF90(grid%perm_zz,perm_zz_p,ierr)
  call VecGetArrayF90(grid%perm_pow,perm_pow_p,ierr)
  call VecGetArrayF90(grid%tor,tor_p,ierr)
  do n = 1,grid%nlmax
    na = grid%nL2A(n)
    nz= int(na/grid%nxy) + 1
    ny= int(mod(na,grid%nxy)/grid%nx) + 1
    nx= mod(mod(na,grid%nxy),grid%nx) + 1

    do ir = 1,grid%iregperm        
      if ((nz>=grid%k1reg(ir)) .and. (nz<=grid%k2reg(ir)) .and.&
          (ny>=grid%j1reg(ir)) .and. (ny<=grid%j2reg(ir)) .and.&
          (nx>= grid%i1reg(ir)) .and. (nx<=grid%i2reg(ir))) then
                                
        val = grid%icap_reg(ir)
       ! call VecSetValue(temp0_nat_vec,n,val,INSERT_VALUES,ierr)
        icap_p(n)=val
            
        val = grid%ithrm_reg(ir)
       ! call VecSetValue(temp1_nat_vec,n,val,INSERT_VALUES,ierr)
        ithrm_p(n)=val
           
        random_nr=1.D0
        if (grid%ran_fac > 0.d0) then
!          frand = rand(0)
          frand = ran1(n)
!          frand = ran_p(n)
!          random_nr = 1.d0+2.d0*grid%ran_fac*(frand-0.5D0)
          random_nr = grid%ran_fac*frand+1.d-6
              
!          print *,'pflowgrid_mod: ',n,frand,random_nr,grid%ran_fac
        endif

        por = grid%por_reg(ir)
        por0_p(n)=por
        if (grid%iran_por==1) then
          por=por*(2.D0**0.666667D0*(frand)**1.5D0)
          if (por<1D-2) por=1D-2 
        endif
       ! call VecSetValue(temp2_nat_vec,n,por,INSERT_VALUES,ierr)
        por_p(n)= por
            
!       nn = 3*n
!       val1 = grid%perm_reg(ir,1)
!       call VecSetValue(temp3_nat_vec,nn,val1,INSERT_VALUES,ierr)
            
!       val2 = grid%perm_reg(ir,2)
!       call VecSetValue(temp3_nat_vec,nn+1,val2,INSERT_VALUES,ierr)
            
!       val3 = grid%perm_reg(ir,3)
!       call VecSetValue(temp3_nat_vec,nn+2,val3,INSERT_VALUES,ierr)
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
  !      call VecSetValue(temp6_nat_vec,n,val4,INSERT_VALUES,ierr)
        perm_pow_p(n)=val4

        val3 = grid%tor_reg(ir)
       ! call VecSetValue(temp7_nat_vec,n,val3,INSERT_VALUES,ierr)
        tor_p(n)=val3
            
!        print *,'setup: ',n+1,ir,i,j,k,val1,val2,val3,random_nr
   !     exit
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
 
 
 ! call VecAssemblyBegin(temp0_nat_vec,ierr)
 ! call VecAssemblyEnd(temp0_nat_vec,ierr)
 ! call VecAssemblyBegin(temp1_nat_vec,ierr)
 ! call VecAssemblyEnd(temp1_nat_vec,ierr)
 ! call VecAssemblyBegin(temp2_nat_vec,ierr)
 ! call VecAssemblyEnd(temp2_nat_vec,ierr)
 ! call VecAssemblyBegin(temp3_nat_vec,ierr)
 ! call VecAssemblyEnd(temp3_nat_vec,ierr)
 ! call VecAssemblyBegin(temp4_nat_vec,ierr)
 ! call VecAssemblyEnd(temp4_nat_vec,ierr)
 ! call VecAssemblyBegin(temp5_nat_vec,ierr)
 ! call VecAssemblyEnd(temp5_nat_vec,ierr)
 ! call VecAssemblyBegin(temp6_nat_vec,ierr)
 ! call VecAssemblyEnd(temp6_nat_vec,ierr)
 ! call VecAssemblyBegin(temp7_nat_vec,ierr)
 ! call VecAssemblyEnd(temp7_nat_vec,ierr)

  
 ! call DANaturalToGlobalBegin(grid%da_1_dof,temp0_nat_vec,INSERT_VALUES, &
 !                               grid%icap,ierr)
 ! call DANaturalToGlobalEnd(grid%da_1_dof,temp0_nat_vec,INSERT_VALUES, &
 !                             grid%icap,ierr)
 ! call DANaturalToGlobalBegin(grid%da_1_dof,temp1_nat_vec,INSERT_VALUES, &
 !                                   grid%ithrm,ierr)
 ! call DANaturalToGlobalEnd(grid%da_1_dof,temp1_nat_vec,INSERT_VALUES, &
 !                             grid%ithrm,ierr)
 ! call DANaturalToGlobalBegin(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
 !                               grid%porosity,ierr)
 ! call DANaturalToGlobalEnd(grid%da_1_dof,temp2_nat_vec,INSERT_VALUES, &
 !                             grid%porosity,ierr)
 ! call DANaturalToGlobalBegin(grid%da_1_dof,temp3_nat_vec,INSERT_VALUES, &
 !                               grid%perm_xx,ierr)
 ! call DANaturalToGlobalEnd(grid%da_1_dof,temp3_nat_vec,INSERT_VALUES, &
 !                             grid%perm_xx,ierr)
 ! call DANaturalToGlobalBegin(grid%da_1_dof,temp4_nat_vec,INSERT_VALUES, &
 !                               grid%perm_yy,ierr)
 ! call DANaturalToGlobalEnd(grid%da_1_dof,temp4_nat_vec,INSERT_VALUES, &
 !                             grid%perm_yy,ierr)
 ! call DANaturalToGlobalBegin(grid%da_1_dof,temp5_nat_vec,INSERT_VALUES, &
 !                               grid%perm_zz,ierr)
 ! call DANaturalToGlobalEnd(grid%da_1_dof,temp5_nat_vec,INSERT_VALUES, &
 !                            grid%perm_zz,ierr)
 ! call DANaturalToGlobalBegin(grid%da_1_dof,temp6_nat_vec,INSERT_VALUES, &
 !                               grid%perm_pow,ierr)
 ! call DANaturalToGlobalEnd(grid%da_1_dof,temp6_nat_vec,INSERT_VALUES, &
 !                             grid%perm_pow,ierr)
 ! call DANaturalToGlobalBegin(grid%da_1_dof,temp7_nat_vec,INSERT_VALUES, &
  !                            grid%tor,ierr)
 
 ! call DANaturalToGlobalEnd(grid%da_1_dof,temp7_nat_vec,INSERT_VALUES, &
  !                            grid%tor,ierr)
                          
! call DANaturalToGlobalBegin(grid%da_3_dof,temp3_nat_vec,INSERT_VALUES, &
!                               grid%perm,ierr)
! call DANaturalToGlobalEnd(grid%da_3_dof,temp3_nat_vec,INSERT_VALUES, &
!                             grid%perm,ierr)

  call VecDestroy(temp0_nat_vec,ierr)
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

  if (grid%iread_perm == 1) then
    call Read_perm_field(grid)
  endif

!geh
  nullify(grid%imat)

  if (grid%iread_geom == 1) then
    call Read_Geom_field(grid)
!geh
  else if (grid%iread_geom == -1) then 
    if (myrank == 0) print *, 'Reading unstructured grid' 
    allocate(grid%pressurebc(grid%nphase,grid%nconnbc)) 

    allocate(grid%imat(grid%nlmax))  ! allocate material id array
    grid%imat = 0      

    call ReadUnstructuredGrid(grid) 
    call pflow_Vadose_initadj(grid)
    call pflow_update_vadose(grid)
    
    ! dump material ids to file in natural ordering
    call DACreateGlobalVector(grid%da_1_dof,temp_vec,ierr)
    call VecGetArrayF90(temp_vec,temp_p,ierr)
    do i=1, grid%nlmax
      temp_p(i) = grid%imat(i)*1.d0      
    enddo
    call VecRestoreArrayF90(temp_vec,temp_p,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'materials.dat',viewer,ierr)
    call VecView(temp_vec,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
    call VecDestroy(temp_vec,ierr)
#if 0    
    print *, 'nconnbc:', grid%nconnbc 
    do nc = 1, grid%nconnbc 
      ibc = grid%ibndtyp(grid%ibconn(nc))
!     print *, 'ibc:', ibc, grid%ibconn(nc) 
      if (ibc == 2) then 
        print *, grid%velocitybc(:,nc) 
      else 
        print *, grid%pressurebc(:,nc), grid%z(grid%nL2A(grid%mblkbc(nc))+1) 
      endif 
    enddo 
#endif
  endif  


  if (grid%using_pflowGrid==0) call VecCopy(grid%Porosity, grid%Porosity0, ierr)
  call VecCopy(grid%perm_xx, grid%perm0_xx, ierr) 
  call VecCopy(grid%perm_yy, grid%perm0_yy, ierr) 
  call VecCopy(grid%perm_zz, grid%perm0_zz, ierr) 
 
 
! Note: VecAssemblyBegin/End needed to run on the Mac - pcl (11/21/03)!
  call VecAssemblyBegin(grid%conc,ierr)
  call VecAssemblyEnd(grid%conc,ierr)

  call VecAssemblyBegin(grid%xmol,ierr)
  call VecAssemblyEnd(grid%xmol,ierr)
! call VecView(grid%conc,PETSC_VIEWER_STDOUT_WORLD,ierr)

  if (myrank == 0) write(*,'("  Finished setting up of INIT ")')

  !-----------------------------------------------------------------------
  ! Initialize field variables
  !-----------------------------------------------------------------------
  select case(grid%ndof)
    case(1)
      call VecCopy(grid%pressure, grid%ppressure, ierr)
      call VecCopy(grid%temp, grid%ttemp, ierr)
    case(2) 
      call pflow_pack_xx2(grid%yy, grid%pressure, grid%nphase, grid%temp, 1, &
                          ierr)
      call VecCopy(grid%yy, grid%xx, ierr)     
    case(3)
      if (grid%use_mph /= PETSC_TRUE .and. grid%use_owg /= PETSC_TRUE &
          .and. grid%use_vadose /= PETSC_TRUE .and. grid%use_flash /= PETSC_TRUE&
          .and. grid%use_richard /= PETSC_TRUE) then   
        call pflow_pack_xx3(grid%yy, grid%pressure, grid%nphase, grid%temp, 1, &
                            grid%conc, 1, ierr)
     
        call VecCopy(grid%yy, grid%xx, ierr)      
      endif 
  end select
   
  if (grid%use_2ph == PETSC_TRUE) then 
    print *, "2 ph begin var arrange"
    call pflow_2phase_initadj(grid)
    print *, "2 ph finish variable initadj"
    call VecCopy(grid%iphas, grid%iphas_old,ierr)
    call pflow_pack_xx4(grid%yy, grid%pressure, grid%nphase, grid%temp, 1, &
         grid%xmol,grid%nphase , grid%sat,grid%nphase , ierr)
    print *, "2 ph finish variable packing"      
    call VecCopy(grid%yy, grid%xx, ierr)
    call pflow_update_2phase(grid)
    print *, "2 ph finish variable update"
  endif 
  
  if (grid%use_mph == PETSC_TRUE) then 
    call pflow_mphase_initadj(grid)
    call VecCopy(grid%iphas, grid%iphas_old,ierr)
    call VecCopy(grid%xx, grid%yy, ierr)
    print *, "m ph finish variable packing"
    call pflow_update_mphase(grid)
  endif 

  if (grid%use_richard == PETSC_TRUE) then 
    call pflow_richard_initadj(grid)
    call VecCopy(grid%iphas, grid%iphas_old,ierr)
    call VecCopy(grid%xx, grid%yy, ierr)
    print *, "richard finish variable packing"
    call pflow_update_richard(grid)
  endif 

  if (grid%use_flash == PETSC_TRUE) then 
    call pflow_flash_initadj(grid)
    call VecCopy(grid%iphas, grid%iphas_old,ierr)
    call VecCopy(grid%xx, grid%yy, ierr)
    print *, "flash finish variable packing"
    !call VecView(grid%xx,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call pflow_update_flash(grid)
  endif 

  if (grid%use_owg == PETSC_TRUE) then 
    call pflow_owg_initadj(grid)
    print*, 'finished owg initadj'
    call VecCopy(grid%iphas, grid%iphas_old,ierr)
    call VecCopy(grid%xx, grid%yy, ierr)
    print *, "OWG finish variable packing"
    call pflow_update_owg(grid)
  endif 

  if (grid%use_vadose == PETSC_TRUE) then 
    call pflow_vadose_initadj(grid)
    call VecCopy(grid%iphas, grid%iphas_old,ierr)
    call VecCopy(grid%xx, grid%yy, ierr)
    print *, "vadose finish variable packing"
    call pflow_update_vadose(grid)
   endif 

  
  
  
!  call VecView(grid%xmol,PETSC_VIEWER_STDOUT_WORLD,ierr)
!  call VecView(grid%yy,PETSC_VIEWER_STDOUT_WORLD,ierr)
! zero initial velocity
  call VecSet(grid%vl,0.d0,ierr)
  if (grid%using_pflowGrid == PETSC_TRUE) call VecSet(grid%vvl,0.d0,ierr)
 
  if (myrank == 0) &
    write(*,'("  Finished setting up of INIT2 ")')

! call VecView(grid%yy,PETSC_VIEWER_STDOUT_WORLD,ierr)

  ! set phase index for each node and initialize accumulation terms
   
  if ( grid%use_owg == PETSC_TRUE) then
    call pflow_owg_initaccum(grid)
  else if (grid%use_mph == PETSC_TRUE) then
    call pflow_mphase_initaccum(grid)
  else if (grid%use_richard == PETSC_TRUE) then
    call pflow_richard_initaccum(grid)
  else if (grid%use_flash == PETSC_TRUE) then
    call pflow_flash_initaccum(grid)
  else if (grid%use_vadose == PETSC_TRUE) then
    call pflow_vadose_initaccum(grid)
  else if (grid%use_2ph == PETSC_TRUE) then
    call pflow_2phase_initaccum(grid)
  else if (grid%ndof > 1) then
 !   call VecSet(grid%iphas,1.d0,ierr)
    call VecGetArrayF90(grid%pressure, pressure_p, ierr)
    call VecGetArrayF90(grid%temp, temp_p, ierr)
    call VecGetArrayF90(grid%density, den_p, ierr)
    call VecGetArrayF90(grid%h, h_p, ierr)
    do m = 1, grid%nlmax
      call wateos_noderiv(temp_p(m),pressure_p(m),dw_kg,dl,hl,grid%scale,ierr)
      den_p(m) = dl
      h_p(m) = hl
    enddo
    call VecRestoreArrayF90(grid%pressure, pressure_p, ierr)
    call VecRestoreArrayF90(grid%temp, temp_p, ierr)
    call VecRestoreArrayF90(grid%density, den_p, ierr)
    call VecRestoreArrayF90(grid%h, h_p, ierr)

  else

    call VecGetArrayF90(grid%pressure, pressure_p, ierr)
    call VecGetArrayF90(grid%temp, temp_p, ierr)
    call VecGetArrayF90(grid%density, den_p, ierr)
    do m = 1, grid%nlmax
      call wateos_noderiv(temp_p(m),pressure_p(m),dw_kg,dl,hl,grid%scale,ierr)
      den_p(m) = dl
    enddo
    call VecRestoreArrayF90(grid%pressure, pressure_p, ierr)
    call VecRestoreArrayF90(grid%temp, temp_p, ierr)
    call VecRestoreArrayF90(grid%density, den_p, ierr)
  endif
  
  !initial solid reaction  
  if (grid%rk > 0.d0) then
    allocate(grid%area_var(grid%nlmax))
    allocate(grid%rate(grid%nlmax))
    call VecGetArrayF90(grid%phis,phis_p,ierr)
    do n = 1, grid%nlmax
      phis_p(n) = grid%phis0
      grid%area_var(n) = 1.d0
    enddo
    call VecRestoreArrayF90(grid%phis,phis_p,ierr)
  endif


   
  
  if (myrank == 0) write(*,'("  Finished setting up ")')

end subroutine pflowGrid_setup



!======================================================================


!#include "pflowgrid_setvel.F90"

!#define PPRESSURE_LOC(j,n) xx_loc_p(j+(n-1)*grid%ndof)

subroutine pflowGrid_setvel (grid, vl_loc, vlbc, ibconn, ibndtyp)

  implicit none

  type(pflowGrid), intent(inout) :: grid

  real*8, pointer :: vl_loc(:), vlbc(:)
  
  integer*4 :: ibndtyp(*), ibconn(*)

  real*8, pointer :: xx_loc_p(:), ddensity_loc_p(:), viscosity_loc_p(:), &
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)

  integer :: ierr 
  integer*4 ibc, ibc_ptran, ip1, ip2, j, jm, jm1, jm2, jng, m, m1, m2, &
             n1, n2, nc, ng
  real*8 :: f1, f2, gravity, v_darcy, perm1, perm2
  real*8 :: D, D1, D2, dd, dd1, dd2, density_ave

  call DAGlobalToLocalBegin(grid%da_ndof, grid%xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, grid%xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)

#if 0
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                            grid%ddensity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                          grid%ddensity_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%viscosity, INSERT_VALUES, &
                            grid%viscosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%viscosity, INSERT_VALUES, &
                          grid%viscosity_loc, ierr)
#endif

  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecGetArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecGetArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
! call VecGetArrayF90(vl, vl_p, ierr)

  do nc = 1, grid%nconn  ! For each interior connection...
    m1 = grid%nd1(nc)
    m2 = grid%nd2(nc)

    n1 = grid%nG2L(m1)
    n2 = grid%nG2L(m2)
    
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)
    ip2 = grid%iperm2(nc)

!   perm1 = perm_loc_p(ip1+3*(m1-1))
!   perm2 = perm_loc_p(ip2+3*(m2-1))

    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(m1)
    else if (ip1 == 2) then
      perm1 = perm_yy_loc_p(m1)
    else
      perm1 = perm_zz_loc_p(m1)
    endif

    if (ip2 == 1) then
      perm2 = perm_xx_loc_p(m2)
    else if (ip2 == 2) then
      perm2 = perm_yy_loc_p(m2)
    else
      perm2 = perm_zz_loc_p(m2)
    endif

    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd
    gravity = grid%fmwh2o * grid%gravity * grid%delz(nc)
    
    do j = 1, grid%nphase
      jm1 = j + (m1-1) * grid%nphase
      jm2 = j + (m2-1) * grid%nphase
      
      ! We need to calculate the "diffusion" constant D at the interface;
      ! it is defined at the cell centers.  We use the harmonic mean of 
      ! the values from the two cells.
      
!     print *, 'setvel: ',grid%myrank,nc,n1,m1,m2,ip1,ip2,perm1,perm2

!     D1 = perm1 / viscosity_loc_p(jm1)
!     D2 = perm2 / viscosity_loc_p(jm2)
!     D = (D1 * D2) / (dd2*D1 + dd1*D2)

      D1 = perm1 * viscosity_loc_p(jm2)
      D2 = perm2 * viscosity_loc_p(jm1)
      D = (perm1 * perm2) / (dd2*D1 + dd1*D2)

      density_ave = f2 * ddensity_loc_p(jm1) + f1 * ddensity_loc_p(jm2)

!     v_darcy = -D * (PPRESSURE_LOC(j,m2) - PPRESSURE_LOC(j,m1) &
      v_darcy = -D*(xx_loc_p(1+(m2-1)*grid%ndof)-xx_loc_p(1+(m1-1)*grid%ndof) &
                - gravity * density_ave)

!     print *,'ptran_init_setvel: ',grid%myrank,nc,m1,m2,n1,n2,v_darcy,D,&
!     ip1,ip2

      ! store average velocity
!     if (n1 > 0) vl_p(ip1+3*(n1-1)) = v_darcy

      vl_loc(nc) = v_darcy
    enddo
  enddo

!---------------------------------------------------------------------------
! Flux terms for boundary nodes.
!---------------------------------------------------------------------------

  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)

    ibc = grid%ibconn(nc)
    ibc_ptran = ibconn(nc)
    ip1 = grid%ipermbc(nc)
    
!   perm1 = perm_loc_p(ip1+3*(ng-1))
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(ng)
    else if (ip1 == 2) then
      perm1 = perm_yy_loc_p(ng)
    else
      perm1 = perm_zz_loc_p(ng)
    endif
    gravity = grid%fmwh2o * grid%gravity * grid%delzbc(nc)
    
!   print *,'pflowgrid_setvel: ',nc,grid%nconnbc,ibc,ibc_ptran,ibndtyp(ibc_ptran)

    if (ibndtyp(ibc_ptran) == 1) then
    
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase

        D = perm1 / grid%viscosity_bc(j) / grid%distbc(nc)
        
        !note: darcy vel. is positive for flow INTO boundary node
!       v_darcy = -D * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,ibc) &
        v_darcy = -D * (xx_loc_p(1+(ng-1)*grid%ndof) - grid%pressurebc(j,nc) &
                - gravity * ddensity_loc_p(jng))
        vlbc(nc) = v_darcy
        
!       print *,'pflowgrid_setvel: ',nc,m,ng,ibc_ptran,ibndtyp(ibc_ptran), &
!       v_darcy,D,gravity
      enddo
      
    else if (ibndtyp(ibc_ptran) == 3) then
    
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase

        D = perm1 / viscosity_loc_p(jng) / grid%distbc(nc)
        
        !note: darcy vel. is positive for flow INTO boundary node
!       v_darcy = -D * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,ibc) &
        v_darcy = -D * (xx_loc_p(1+(ng-1)*grid%ndof) - grid%pressurebc(j,nc) &
                - gravity * ddensity_loc_p(jng))
        vlbc(nc) = v_darcy
        
!       if (grid%t/grid%tconv > 0.3) &
!       print *,'pflowgrid_setvel: ',grid%myrank,nc,m,ng,ibc_ptran, &
!       ibndtyp(ibc_ptran),ddensity_loc_p(jng),viscosity_loc_p(jng),perm1, &
!       v_darcy,D,gravity,xx_loc_p(1+(ng-1)*grid%ndof),grid%pressurebc(j,ibc), &
!       grid%distbc(nc)
      enddo
    
    else if (ibndtyp(ibc_ptran) == 2) then
      vlbc(nc) = 0.d0
    endif
    
!   if (grid%t/grid%tconv > 45.) &
!   print *,'pflotgrid_setvel: ',grid%myrank,nc,m,ng,ibc_ptran, &
!   ibndtyp(ibc_ptran),vlbc(nc)
  enddo

  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecRestoreArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
! call VecRestoreArrayF90(vl, vl_p, ierr)  

  end subroutine pflowGrid_setvel

!#undef PPRESSURE_LOC

!======================================================================

!#include "pflowgrid_compute_xyz.F90"

subroutine pflowGrid_compute_xyz(grid)
  
  implicit none
  
  type(pflowGrid), intent(inout) :: grid

! integer :: ierr
  integer*4 :: i, j, k, n
  integer :: prevnode
    ! prevnode is used to hold the natural numbering of the node that is the
    ! previous node in either the x, y, or z direction.
! Vec :: dx_nat, dy_nat, dz_nat
! Vec :: dx_all, dy_all, dz_all  ! Holds contents of dx_nat et al. on proc 0.
! real*8, pointer :: dx_p(:), dy_p(:), dz_p(:)

! call DACreateNaturalVector(grid%da_1_dof, dx_nat, ierr)
! call DACreateNaturalVector(grid%da_1_dof, dy_nat, ierr)
! call DACreateNaturalVector(grid%da_1_dof, dz_nat, ierr)

! call DAGlobalToNaturalBegin(grid%da_1_dof,grid%dx,INSERT_VALUES,dx_nat,ierr)
! call DAGlobalToNaturalEnd(grid%da_1_dof,grid%dx,INSERT_VALUES,dx_nat,ierr)
! call DAGlobalToNaturalBegin(grid%da_1_dof,grid%dy,INSERT_VALUES,dy_nat,ierr)
! call DAGlobalToNaturalEnd(grid%da_1_dof,grid%dy,INSERT_VALUES,dy_nat,ierr)
! call DAGlobalToNaturalBegin(grid%da_1_dof,grid%dz,INSERT_VALUES,dz_nat,ierr)
! call DAGlobalToNaturalEnd(grid%da_1_dof,grid%dz,INSERT_VALUES,dz_nat,ierr)

#ifdef HAVE_MPITOMPIZERO
! call VecConvertMPIToMPIZero(dx_nat, dx_all, ierr)
! call VecConvertMPIToMPIZero(dy_nat, dy_all, ierr)
! call VecConvertMPIToMPIZero(dz_nat, dz_all, ierr)
#else
! call VecConvertMPIToSeqAll(dx_nat, dx_all, ierr)
! call VecConvertMPIToSeqAll(dy_nat, dy_all, ierr)
! call VecConvertMPIToSeqAll(dz_nat, dz_all, ierr)
#endif

!  if (grid%myrank == 0) then
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

          if (i == 1) then
!           grid%x(n) = 0.5d0 * dx_p(n) 
            grid%x(n) = 0.5d0 * grid%dx0(i) 
          else
            prevnode = n-1
!           grid%x(n) = grid%x(prevnode) + 0.5d0*(dx_p(prevnode) + dx_p(n))
            grid%x(n) = grid%x(prevnode) + 0.5d0*(grid%dx0(i-1) + grid%dx0(i))
          endif

          if (j == 1) then
!           grid%y(n) = 0.5d0 * dy_p(n)
            grid%y(n) = 0.5d0 * grid%dy0(j)
          else
!           prevnode = i + (j-2)*grid%nx + (k-1)*grid%nx*grid%ny
            prevnode = n - grid%nx
!           grid%y(n) = grid%y(prevnode) + 0.5d0*(dy_p(prevnode) + dy_p(n))
            grid%y(n) = grid%y(prevnode) + 0.5d0*(grid%dy0(j-1) + grid%dy0(j))
          endif

          if (k == 1) then
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

!======================================================================

!#include "pflowgrid_update_dt.F90"

subroutine pflowgrid_update_dt(grid, its)

  implicit none

  type(pflowGrid), intent(inout) :: grid
  integer, intent(in) :: its

  real*8 :: fac,dtt,up,utmp,uc,ut,uus
  
  if (grid%iaccel == 0) return

  if (grid%use_thc == PETSC_TRUE) then
  
    fac = 0.5d0
    if (its >= grid%iaccel) fac = 0.33d0
    up = grid%dpmxe/(grid%dpmax+0.1)
    utmp = grid%dtmpmxe/(grid%dtmpmax+1.d-5)
    uc = grid%dcmxe/(grid%dcmax+1.d-6)
    ut = min(up,utmp,uc)
    dtt = fac * grid%dt * (1.d0 + ut)

  else if (grid%use_2ph == PETSC_TRUE) then
  
    fac = 0.5d0
    if (its >= grid%iaccel) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = grid%dpmxe/(grid%dpmax+0.1)
      utmp = grid%dtmpmxe/(grid%dtmpmax+1.d-5)
      uc = grid%dcmxe/(grid%dcmax+1.d-6)
      uus=(0.01D0/(grid%dsmax+1.d-6))**2
      ut = min(up,utmp,uc)
    endif
    dtt = fac * grid%dt * (1.d0 + ut)
   
   else if (grid%use_mph == PETSC_TRUE) then
  
    fac = 0.5d0
    if (its >= grid%iaccel) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = grid%dpmxe/(grid%dpmax+0.1)
      utmp = grid%dtmpmxe/(grid%dtmpmax+1.d-5)
      uc = grid%dcmxe/(grid%dcmax+1.d-6)
      uus= grid%dsmxe/(grid%dsmax+1.d-6)
      ut = min(up,utmp,uc,uus)
    endif
    dtt = fac * grid%dt * (1.d0 + ut)
    
  else if (grid%use_richard == PETSC_TRUE) then
  
    fac = 0.5d0
    if (its >= grid%iaccel) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = grid%dpmxe/(grid%dpmax+0.1)
      utmp = grid%dtmpmxe/(grid%dtmpmax+1.d-5)
      uus= grid%dsmxe/(grid%dsmax+1.d-6)
      ut = min(up,utmp,uus)
    endif
    dtt = fac * grid%dt * (1.d0 + ut)
 
!************** FLASH *************************
   else if (grid%use_flash == PETSC_TRUE) then
  
    fac = 0.5d0
    if (its >= grid%iaccel) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = grid%dpmxe/(grid%dpmax+0.1)
      utmp = grid%dtmpmxe/(grid%dtmpmax+1.d-5)
      uc = grid%dcmxe/(grid%dcmax+1.d-6)
      uus= grid%dsmxe/(grid%dsmax+1.d-6)
      ut = min(up,utmp,uc,uus)
    endif
    dtt = fac * grid%dt * (1.d0 + ut)


   else if (grid%use_owg == PETSC_TRUE) then
  
    fac = 0.5d0
    if (its >= grid%iaccel) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = grid%dpmxe/(grid%dpmax+0.1)
      utmp = grid%dtmpmxe/(grid%dtmpmax+1.d-5)
      uc = grid%dcmxe/(grid%dcmax+1.d-6)
      uus= grid%dsmxe/(grid%dsmax+1.d-6)
      ut = min(up,utmp,uc,uus)
    endif
    dtt = fac * grid%dt * (1.d0 + ut)

   else if (grid%use_vadose == PETSC_TRUE) then
  
    fac = 0.5d0
    if (its >= grid%iaccel) then
      fac = 0.33d0
      ut = 0.d0

    else
      up = grid%dpmxe/(grid%dpmax+0.1)
      utmp = grid%dtmpmxe/(grid%dtmpmax+1.d-5)
      uc = grid%dcmxe/(grid%dcmax+1.d-6)
      uus= grid%dsmxe/(grid%dsmax+1.d-6)
      ut = min(up,utmp,uc,uus)
    endif
    dtt = fac * grid%dt * (1.d0 + ut)
  

  else if (grid%use_th == PETSC_TRUE) then
  
    fac = 0.5d0
    if (its >= grid%iaccel) fac = 0.33d0
    up = grid%dpmxe/(grid%dpmax+0.1)
    utmp = grid%dtmpmxe/(grid%dtmpmax+1.d-5)
    ut = min(up,utmp)
    dtt = fac * grid%dt * (1.d0 + ut)

  else if (grid%use_cond == PETSC_TRUE) then

    fac = 0.5d0
    if (its >= grid%iaccel) fac = 0.33d0
    ut = grid%dtmpmxe/(grid%dpmax+1.e-5)
    dtt = fac * grid%dt * (1.d0 + ut)
    if (dtt > 2.d0 * grid%dt) dtt = 2.d0 * grid%dt 

  else
  
    if (its <= grid%iaccel .and. its <= size(grid%tfac)) then
      if (its == 0) then
        dtt = grid%tfac(1) * grid%dt
      else
        dtt = grid%tfac(its) * grid%dt
      endif
    endif
  endif
  
  if (dtt > 2.d0 * grid%dt) dtt = 2.d0 * grid%dt 
  if (dtt > grid%dt_max) dtt = grid%dt_max
  if (dtt>.25d0*grid%t .and. grid%t>1.d-2) dtt=.25d0*grid%t
  grid%dt = dtt

  end subroutine pflowgrid_update_dt

!======================================================================

!#include "pflowgrid_step.F90"

subroutine pflowGrid_step(grid,ntstep,kplt,iplot,iflgcut,ihalcnt,its)
  
  use translator_mph_module, only : translator_mph_step_maxchange
  use translator_owg_module, only : translator_owg_step_maxchange
  use translator_vad_module, only : translator_vad_step_maxchange
  use translator_flash_module, only : translator_flash_step_maxchange
  use translator_Richard_module, only : translator_ric_step_maxchange
  use pflow_output_module
  use TTPHASE_module
  use MPHASE_module
  use Flash_module
  use OWG_module
  use vadose_module
  use Richard_module
  use pflow_solv_module
  implicit none

  type(pflowGrid), intent(inout) :: grid

  integer :: ierr, ihalcnt, ns !i,m,n,ix,jy,kz,j
  integer :: its,kplt,iplot,ntstep !,idpmax,idtmpmax,idcmax
  integer :: icut, iflgcut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  integer update_reason
  real*8 :: tsrc
! real*8, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  its = 0
  icut = 0

  ! Perform some global-to-local scatters to update the ghosted vectors.
  ! We have to do this so that the routines for calculating the residual
  ! and the Jacobian will have the ghost points they need.
  ! Note that we don't do the global-to-local scatter for the ppressure 
  ! vector, as that needs to be done within the residual calculation routine
  ! because that routine may get called several times during one Newton step
  ! if a method such as line search is being used.
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                            grid%porosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%porosity, INSERT_VALUES, &
                          grid%porosity_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%tor, INSERT_VALUES, &
                            grid%tor_loc, ierr)

  call DAGlobalToLocalEnd(grid%da_1_dof, grid%tor, INSERT_VALUES, &
                          grid%tor_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_1_dof, grid%icap, INSERT_VALUES, &
                            grid%icap_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%icap, INSERT_VALUES, &
                          grid%icap_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_1_dof, grid%ithrm, INSERT_VALUES, &
                            grid%ithrm_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%ithrm, INSERT_VALUES, &
                          grid%ithrm_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_1_dof, grid%iphas, INSERT_VALUES, &
                            grid%iphas_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%iphas, INSERT_VALUES, &
                          grid%iphas_loc, ierr)

  grid%t = grid%t + grid%dt
  grid%flowsteps = grid%flowsteps + 1

! Adjust time step to plot time
  if (grid%t + 0.2*grid%dt >= grid%tplot(kplt)) then
    grid%t = grid%t - grid%dt
    grid%dt = grid%tplot(kplt) - grid%t
    if (grid%dt > grid%dt_max) then
      grid%dt = grid%dt_max
      grid%t = grid%t + grid%dt
    else
      grid%t = grid%tplot(kplt)
      iplot = 1
    endif
  else if (grid%flowsteps == grid%stepmax) then
    iplot = 1
  endif

! source/sink time step control
  if (grid%nblksrc > 0) then
    ns = 1
    tsrc = grid%timesrc(grid%isrc1,ns)
    if (grid%t >= tsrc ) then
      if (grid%t > tsrc +1D2) then
        grid%t = grid%t - grid%dt
        grid%dt = tsrc - grid%t
        grid%t = tsrc
      endif
      grid%isrc1 = grid%isrc1 + 1
    endif
  endif
  
  if (iflgcut == 1) then
    ihalcnt = ihalcnt + 1
    if (ihalcnt > grid%ndtcmx) then
      iflgcut = 0
      ihalcnt = 0
    endif
  endif
  
  !set maximum time step
  if (ntstep > grid%nstpmax) then
  ! do nothing - keep current dt_max
  else if (grid%t > grid%tstep(ntstep)) then
    ntstep = ntstep + 1
    if (ntstep <= grid%nstpmax) then
      grid%dt_max = grid%dtstep(ntstep)
    else
      grid%dt_max = grid%dtstep(grid%nstpmax)
    endif
  endif
  
  do
    
    grid%iphch=0
    if (grid%use_cond == PETSC_TRUE) then
      call SNESSolve(grid%snes, PETSC_NULL, grid%ttemp, ierr)
    else if (grid%use_th == PETSC_TRUE) then
      call SNESSolve(grid%snes, PETSC_NULL, grid%xx, ierr)
    else if (grid%use_thc == PETSC_TRUE) then
      call SNESSolve(grid%snes, PETSC_NULL, grid%xx, ierr)
    else if (grid%use_2ph == PETSC_TRUE) then
   ! call  TTPhase_Update(grid%xx,grid)
      call SNESSolve(grid%snes, PETSC_NULL, grid%xx, ierr)
    else if (grid%use_mph == PETSC_TRUE) then
      if (grid%use_ksp == PETSC_TRUE) then
        call pflow_solve(grid,its,snes_reason,ierr)
      else 
        call SNESSolve(grid%snes, PETSC_NULL, grid%xx, ierr)
      endif
    else if (grid%use_richard == PETSC_TRUE) then
      if (grid%use_ksp == PETSC_TRUE) then
        call pflow_solve(grid,its,snes_reason,ierr)
      else 
        call SNESSolve(grid%snes, PETSC_NULL, grid%xx, ierr)
      endif
   else if (grid%use_flash == PETSC_TRUE) then
      if (grid%use_ksp == PETSC_TRUE) then
        call pflow_solve(grid,its,snes_reason,ierr)
      else 
        call SNESSolve(grid%snes, PETSC_NULL, grid%xx, ierr)
      endif
    else if (grid%use_vadose == PETSC_TRUE) then
      if (grid%use_ksp == PETSC_TRUE) then
        call pflow_solve(grid,its,snes_reason,ierr)
      else 
        call SNESSolve(grid%snes, PETSC_NULL, grid%xx, ierr)
      endif
    else if (grid%use_owg == PETSC_TRUE) then
      if (grid%use_ksp == PETSC_TRUE) then
        call pflow_solve(grid,its,snes_reason,ierr)
      else 
        call SNESSolve(grid%snes, PETSC_NULL, grid%xx, ierr)
      endif
    else
      call SNESSolve(grid%snes, PETSC_NULL, grid%ppressure, ierr)
    endif
  ! print *,'pflow_step, finish SNESSolve'
   call MPI_Barrier(PETSC_COMM_WORLD,ierr)
   if (grid%use_ksp /= PETSC_TRUE) &
     call SNESGetIterationNumber(grid%snes, its, ierr)
!   call KSPGetIterationNumber(grid%ksp, kspits, ierr)
  
!   printout residual and jacobian to screen for testing purposes
!   call VecView(grid%r,PETSC_VIEWER_STDOUT_WORLD,ierr)
!   call MatView(grid%J,PETSC_VIEWER_STDOUT_WORLD,ierr)

!   call SNESDefaultMonitor(grid%snes, its, r2norm, ierr)

    if (grid%use_ksp /= PETSC_TRUE) &
      call SNESGetConvergedReason(grid%snes, snes_reason, ierr)

!   parameter (SNES_CONVERGED_ITERATING         =  0)
!   parameter (SNES_CONVERGED_FNORM_ABS         =  2)
!   parameter (SNES_CONVERGED_FNORM_RELATIVE    =  3)
!   parameter (SNES_CONVERGED_PNORM_RELATIVE    =  4)
!   parameter (SNES_CONVERGED_GNORM_ABS         =  5)
!   parameter (SNES_CONVERGED_TR_REDUCTION      =  6)
!   parameter (SNES_CONVERGED_TR_DELTA          =  7)

!   parameter (SNES_DIVERGED_FUNCTION_COUNT     = -2)
!   parameter (SNES_DIVERGED_FNORM_NAN          = -4)
!   parameter (SNES_DIVERGED_MAX_IT             = -5)
!   parameter (SNES_DIVERGED_LS_FAILURE         = -6)
!   parameter (SNES_DIVERGED_TR_REDUCTION       = -7)
!   parameter (SNES_DIVERGED_LOCAL_MIN          = -8)


!******************************************************************
! since calculation on saturation has special requirments, need special treatment
! update_reason
    !call PETScBarrier(PETSC_NULL_OBJECT, ierr)
    update_reason = 1
    if ((grid%use_2ph == PETSC_TRUE).and.(snes_reason >= 0)) then
      call TTPhase_Update_Reason(update_reason, grid)
      if (grid%myrank==0) print *,'update_reason: ',update_reason
    endif

    if ((grid%use_mph == PETSC_TRUE).and.(snes_reason >= 0)) then
      call MPhase_Update_Reason(update_reason, grid)
      if (grid%myrank==0) print *,'update_reason: ',update_reason
    endif

    if ((grid%use_flash == PETSC_TRUE).and.(snes_reason >= 0)) then
      call flash_Update_Reason(update_reason, grid)
      if (grid%myrank==0) print *,'update_reason: ',update_reason
    endif

    if ((grid%use_vadose == PETSC_TRUE).and.(snes_reason >= 0)) then
      call Vadose_Update_Reason(update_reason, grid)
      if (grid%myrank==0) print *,'update_reason: ',update_reason
    endif

    if ((grid%use_richard == PETSC_TRUE).and.(snes_reason >= 0)) then
      call Vadose_Update_Reason(update_reason, grid)
      if (grid%myrank==0) print *,'update_reason: ',update_reason
    endif


    if ((grid%use_owg == PETSC_TRUE).and.(snes_reason >= 0)) then
      call OWG_Update_Reason(update_reason, grid)
      if (grid%myrank==0) print *,'update_reason: ',update_reason
    endif



!******************************************************************
    
    if (snes_reason < 0 .or. update_reason <= 0) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      iflgcut = 1

      if (icut > grid%icut_max .or. grid%dt<1d-9 ) then
!       call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)
        if (grid%myrank == 0) then
!         t = pflowgrid_get_t(grid)
          print *,"icut_max exceeded. icut= ",icut, "t= ",grid%t/grid%tconv, &
                   grid%dt, ". Stopping execution."
          print *, "Final state:"
        endif
        iplot = 1
 !       call pflow_output(grid%ppressure,grid%ttemp,grid%conc,grid%phis,grid%porosity, &
 !       grid%perm_xx, grid%perm_yy, grid%perm_zz, &
 !       grid%porosity, grid%sat, grid%vl, &
 !       grid%c_nat,grid%vl_nat,grid%p_nat,grid%t_nat,grid%s_nat,grid%phis_nat,grid%por_nat, &
 !       grid%ibrkface, grid%jh2o, grid%nphase, grid%nmax, &
 !       grid%snes, &
 !       grid%t, grid%dt, grid%tconv, grid%flowsteps, grid%rk, &
 !       grid%k1brk,grid%k2brk,grid%j1brk,grid%j2brk,grid%i1brk,grid%i2brk, &
 !       grid%nx,grid%ny,grid%nz,grid%nxy,grid%dx0,grid%dy0,grid%dz0,grid%x,grid%y,grid%z, &
 !       grid%da_nphase_dof,grid%da_1_dof,grid%da_3np_dof,grid%da_ndof,grid%ndof, &
 !       kplt,iplot,grid%iprint,grid%ibrkcrv, &
 !       grid%itecplot,grid%write_init,grid%myrank)
    
  !      call pflow_output(grid,kplt,iplot)
        ! The above line won't work when scope is restricted properly!
        ! Replace this with a different function call!
        call pflowgrid_destroy(grid)
        call PetscFinalize(ierr)
        stop
      endif

      grid%t = grid%t - grid%dt
      grid%dt = 0.5d0 * grid%dt
      grid%t = grid%t + grid%dt
    
      if (grid%myrank == 0) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.4, '' dt= '', &
        &   1pe12.4,i2)')  snes_reason,icut,grid%icutcum, &
            grid%t/grid%tconv,grid%dt/grid%tconv,iflgcut

      if (grid%ndof == 1) then
        ! VecCopy(x,y): y=x
        call VecCopy(grid%pressure, grid%ppressure, ierr)
        call VecCopy(grid%temp, grid%ttemp, ierr)
      else
        if (grid%use_owg==PETSC_TRUE) then
          call pflow_owg_timecut(grid)
        elseif (grid%use_mph==PETSC_TRUE) then
          call pflow_mphase_timecut(grid)
        elseif (grid%use_flash==PETSC_TRUE) then
          call pflow_flash_timecut(grid)
        elseif (grid%use_richard==PETSC_TRUE) then
          call pflow_richard_timecut(grid)
        elseif (grid%use_vadose==PETSC_TRUE) then
          call pflow_vadose_timecut(grid)
        else
          call VecCopy(grid%h, grid%hh, ierr)
          call VecCopy(grid%yy, grid%xx, ierr)
          call VecCopy(grid%density, grid%ddensity, ierr)
        endif
        call VecCopy(grid%iphas_old, grid%iphas, ierr)
      endif

    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  grid%newtcum = grid%newtcum + its
  grid%icutcum = grid%icutcum + icut

! print screen output
  if (grid%myrank == 0) then
    if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
      write(*, '(/," FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1,"]", &
        & " snes: ",i4,/,"  newt= ",i2," [",i6,"]"," cut= ",i2," [",i4,"]")') &
        grid%flowsteps,grid%t/grid%tconv,grid%dt/grid%tconv,grid%tunit, &
        snes_reason,its,grid%newtcum,icut,grid%icutcum
    
      write(IUNIT2, '(" FLOW ",i6," Time= ",1pe12.4," Dt= ",1pe12.4," [",a1, &
        & "]"," snes: ",i4,/,"  newt= ",i2," [",i6,"]"," cut= ",i2," [",i4, &
        & "]")') grid%flowsteps,grid%t/grid%tconv,grid%dt/grid%tconv, &
        grid%tunit, snes_reason,its,grid%newtcum,icut,grid%icutcum
    endif
  endif
  
  ! calculate maxium changes in fields over a time step

  if (grid%ndof == 1 .and. grid%use_cond == PETSC_TRUE) then
  
    call VecWAXPY(grid%dp,-1.d0,grid%ttemp,grid%temp,ierr)
    call VecStrideNorm(grid%dp,0,NORM_INFINITY,grid%dpmax,ierr)
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0) then
        write(*,'("  --> max chng: dTmx= ",1pe12.4)') grid%dpmax
        write(IUNIT2,'("  --> max chng: dTmx= ",1pe12.4)') grid%dpmax
      endif
    endif
    
  else if (grid%ndof == 2 .and. grid%use_th == PETSC_TRUE) then
  
    call VecWAXPY(grid%dxx,-1.d0,grid%xx,grid%yy,ierr)
!   call VecAbs(grid%dxx,ierr)
!   call VecMax(grid%dxx,grid%idxxmax,grid%dxxmax,ierr)
!   call VecStrideMax(grid%dxx,1,PETSC_NULL,grid%dxxmax,ierr)
    call VecStrideNorm(grid%dxx,0,NORM_INFINITY,grid%dpmax,ierr)
    call VecStrideNorm(grid%dxx,1,NORM_INFINITY,grid%dtmpmax,ierr)
!   call VecStrideNorm(grid%dxx,2,NORM_INFINITY,grid%dcmax,ierr)
!   n = grid%idxxmax/grid%ndof+1
!   kz = n/grid%nxy+1
!   jy = n/grid%nx+1 - (kz-1)*grid%ny
!   ix = n - (jy-1)*grid%nx - (kz-1)*grid%nxy
!   j = grid%idxxmax - (n-1)*grid%ndof
!   if (grid%myrank==0) print *,'  --> max chng: ',grid%idxxmax,grid%dxxmax, &
!     ' @ node ',ix,jy,kz,j,grid%ndof,grid%nx,grid%nxy
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4)') grid%dpmax,grid%dtmpmax
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4)') grid%dpmax,grid%dtmpmax
      endif
    endif
    
  else if (grid%ndof == 3 .and. grid%use_thc == PETSC_TRUE) then
  
    call VecWAXPY(grid%dxx,-1.d0,grid%xx,grid%yy,ierr)
!   call VecAbs(grid%dxx,ierr)
!   call VecMax(grid%dxx,grid%idxxmax,grid%dxxmax,ierr)
!   call VecStrideMax(grid%dxx,1,PETSC_NULL,grid%dxxmax,ierr)
    call VecStrideNorm(grid%dxx,0,NORM_INFINITY,grid%dpmax,ierr)
    call VecStrideNorm(grid%dxx,1,NORM_INFINITY,grid%dtmpmax,ierr)
    call VecStrideNorm(grid%dxx,2,NORM_INFINITY,grid%dcmax,ierr)
!   n = grid%idxxmax/grid%ndof+1
!   kz = n/grid%nxy+1
!   jy = n/grid%nx+1 - (kz-1)*grid%ny
!   ix = n - (jy-1)*grid%nx - (kz-1)*grid%nxy
!   j = grid%idxxmax - (n-1)*grid%ndof
!   if (grid%myrank==0) print *,'  --> max chng: ',grid%idxxmax,grid%dxxmax, &
!     ' @ node ',ix,jy,kz,j,grid%ndof,grid%nx,grid%nxy
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax
      endif
    endif
    
  else if (grid%ndof == 4 .and. grid%use_2ph == PETSC_TRUE) then
  
    call VecWAXPY(grid%dxx,-1.d0,grid%xx,grid%yy,ierr)
    call VecStrideNorm(grid%dxx,0,NORM_INFINITY,grid%dpmax,ierr)
    call VecStrideNorm(grid%dxx,1,NORM_INFINITY,grid%dtmpmax,ierr)
    call VecStrideNorm(grid%dxx,2,NORM_INFINITY,grid%dcmax,ierr)
    call VecStrideNorm(grid%dxx,3,NORM_INFINITY,grid%dsmax,ierr)
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
      endif
    endif

  else if (grid%use_mph == PETSC_TRUE) then
     call translator_mph_step_maxchange(grid)
    ! note use mph will use variable switching, the x and s change is not meaningful 
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
      endif
    endif


  else if (grid%use_richard == PETSC_TRUE) then
     call translator_ric_step_maxchange(grid)
    ! note use mph will use variable switching, the x and s change is not meaningful 
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax, grid%dsmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dsmax
      endif
    endif


  else if (grid%use_flash == PETSC_TRUE) then
     call translator_flash_step_maxchange(grid)
    ! note use mph will use variable switching, the x and s change is not meaningful 
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
      endif
    endif


  else if (grid%use_vadose == PETSC_TRUE) then
     call translator_vad_step_maxchange(grid)
    ! note use mph will use variable switching, the x and s change is not meaningful 
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
      endif
    endif

  else if (grid%use_owg == PETSC_TRUE) then
    call translator_owg_step_maxchange(grid)
   
      
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
        
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          grid%dpmax,grid%dtmpmax,grid%dcmax,grid%dsmax
      endif
    endif
    
  else ! use_liquid
  
    call VecWAXPY(grid%dp,-1.d0,grid%ppressure,grid%pressure,ierr)
!   call VecAbs(grid%dp,ierr)
!   call VecMax(grid%dp,idpmax,dpmax,ierr)
    call VecStrideNorm(grid%dp,0,NORM_INFINITY,grid%dpmax,ierr)
    if (grid%myrank==0) then
      if (mod(grid%flowsteps,grid%imod) == 0 .or. grid%flowsteps == 1) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4)') grid%dpmax
        write(IUNIT2,'("  --> max chng: dpmx= ",1pe12.4)') grid%dpmax
      endif
    endif
  endif

  if (grid%myrank == 0 .and. mod(grid%flowsteps,grid%imod) == 0) then
    print *, ""
  endif

end subroutine pflowGrid_step

!==========================================================================
  
subroutine pflowGrid_update (grid)
  
  use pflow_vector_ops_module
  use TTPHASE_module
  use MPHASE_module
  use Flash_module
  use OWG_module
  use Vadose_module
  use Richard_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
  
  integer :: ierr
  integer*4 m, n
  real*8, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:), phis_p(:)

! update solution vector and physical properties (VecCopy(x,y): y=x)
  if (grid%ndof == 1) then
    call VecCopy(grid%ppressure, grid%pressure, ierr)
    call VecCopy(grid%ttemp, grid%temp, ierr)
    call VecCopy(grid%ddensity, grid%density, ierr)
  else
    if (grid%ndof <= 3 .and. grid%use_mph/=PETSC_TRUE .and.  &
        grid%use_owg/=PETSC_TRUE .and. grid%use_vadose/=PETSC_TRUE &
         .and. grid%use_flash/=PETSC_TRUE .and. grid%use_richard/=PETSC_TRUE) then
      call VecCopy(grid%xx, grid%yy, ierr)
      call VecCopy(grid%hh, grid%h, ierr)
      call VecCopy(grid%ddensity, grid%density, ierr)
    endif    
  endif
 

! call VecView(grid%ppressure,PETSC_VIEWER_STDOUT_WORLD,ierr)
  if (grid%use_owg == PETSC_TRUE) then
    call pflow_update_owg(grid)
  elseif (grid%use_mph == PETSC_TRUE) then
    call pflow_update_mphase(grid)
  elseif (grid%use_richard == PETSC_TRUE) then
    call pflow_update_richard(grid)
  elseif (grid%use_flash == PETSC_TRUE) then
    call pflow_update_flash(grid)
  elseif (grid%use_vadose == PETSC_TRUE) then
    call pflow_update_vadose(grid)
  elseif (grid%use_2ph == PETSC_TRUE) then
!   call pflow_update_fldvar (grid%yy, grid%pressure, grid%temp, grid%sat, &
!   grid%xmol, grid%density, grid%porosity, grid%h, grid%accum, &
!   grid%dencpr, grid%ithrm, grid%iphas, grid%scale, &
!   grid%eqkair, grid%nlmax, grid%ndof, grid%nphase, grid%jgas, grid%jh2o)
   !print *,' Into pflow_update_2phase'
    call pflow_update_2phase(grid)  
   !print *,' out pflow_update_2phase'
  else if (grid%ndof > 1) then

    call VecGetArrayF90(grid%xx, xx_p, ierr)
    call VecGetArrayF90(grid%pressure, press_p, ierr)
    call VecGetArrayF90(grid%temp, temp_p, ierr)
    if (grid%ndof == 3) call VecGetArrayF90(grid%conc, conc_p, ierr)
    do m = 1, grid%nlmax
      press_p(m) = xx_p(1+(m-1)*grid%ndof)
      temp_p(m) = xx_p(2+(m-1)*grid%ndof)
      if (grid%ndof == 3) conc_p(m) = xx_p(3+(m-1)*grid%ndof)
    enddo
    call VecRestoreArrayF90(grid%xx, xx_p, ierr)
    call VecRestoreArrayF90(grid%pressure, press_p, ierr)
    call VecRestoreArrayF90(grid%temp, temp_p, ierr)
    if (grid%ndof == 3) call VecRestoreArrayF90(grid%conc, conc_p, ierr)
  endif

  if (grid%using_pflowGrid == PETSC_TRUE) then
!   call VecCopy(grid%vvl,grid%vl,ierr)
    grid%vl_loc = grid%vvl_loc
    grid%vlbc = grid%vvlbc
    grid%vg_loc = grid%vvg_loc
    grid%vgbc = grid%vvgbc
    grid%xphi_co2 = grid%xxphi_co2
    grid%den_co2=grid%dden_co2
  endif
  
  !integrate solid volume fraction using explicit finite difference
  if (grid%rk > 0.d0) then
    call VecGetArrayF90(grid%phis,phis_p,ierr)
    do n = 1, grid%nlmax
      phis_p(n) = phis_p(n) + grid%dt * grid%vbars * grid%rate(n)
      if (phis_p(n) < 0.d0) phis_p(n) = 0.d0
      grid%area_var(n) = (phis_p(n)/grid%phis0)**grid%pwrsrf
      
!     print *,'update: ',n,phis_p(n),grid%rate(n),grid%area_var(n)
    enddo
    call VecRestoreArrayF90(grid%phis,phis_p,ierr)
  endif

end subroutine pflowGrid_update

!======================================================================


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


!#include "pflowgrid_readinput.F90"

subroutine pflowGrid_read_input(grid, inputfile)
  
  ! keywords: GRID, PROC, COUP, GRAV, OPTS, TOLR, DXYZ, DIFF, RADN, HYDR,  
  !           SOLV, THRM, PCKR, PHIK, INIT, TIME, DTST, BCON, SOUR, BRK, RCTR
  
  use fileio_module
  
  implicit none

  type(pflowGrid), intent(inout) :: grid
  character(len=*), intent(in) :: inputfile
  
  integer :: ierr
#include "definitions.h"
  character(len=MAXSTRINGLENGTH) :: string 
  character(len=MAXWORDLENGTH) :: word, strtim
  character(len=MAXCARDLENGTH) :: card
  integer :: i, i1, i2, idum, ireg, isrc, j
  integer :: ibc, ibrk, ir,np
  
  
 

  open(IUNIT1, file=inputfile, action="read", status="old")
  
  do
    call fiReadFlotranString(IUNIT1, string, ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    call fiCharsToUpper(word,len_trim(word))
    call fiReadCard(word,card,ierr)

    if (grid%myrank == 0) print *, card

    select case(card)

!....................

      case ('GRID')

!....................

      case ('PROC')

!....................

      case ('COUP')

        call fiReadStringErrorMsg('COUP',ierr)

        call fiReadInt(string,grid%isync,ierr)
        call fiDefaultMsg('isync',ierr)

        if (grid%myrank == 0) &
          write(IUNIT2,'(/," *COUP",/, &
            & "  isync      = ",3x,i2 &
            & )') grid%isync

!....................

      case ('GRAV')

        call fiReadStringErrorMsg('GRAV',ierr)

        call fiReadDouble(string,grid%gravity,ierr)
        call fiDefaultMsg('gravity',ierr)

        if (grid%myrank == 0) &
          write(IUNIT2,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1pe12.4 &
            & )') grid%gravity

!....................

      case ('OPTS')

        call fiReadStringErrorMsg('OPTS',ierr)

        call fiReadInt(string,grid%write_init,ierr)
        call fiDefaultMsg('write_init',ierr)

        call fiReadInt(string,grid%iprint,ierr)
        call fiDefaultMsg('iprint',ierr)

        call fiReadInt(string,grid%imod,ierr)
        call fiDefaultMsg('mod',ierr)

        call fiReadInt(string,grid%itecplot,ierr)
        call fiDefaultMsg('itecplot',ierr)

        call fiReadInt(string,grid%iblkfmt,ierr)
        call fiDefaultMsg('iblkfmt',ierr)

        call fiReadInt(string,grid%ndtcmx,ierr)
        call fiDefaultMsg('ndtcmx',ierr)

        call fiReadInt(string,grid%iran_por,ierr)
        call fiDefaultMsg('iran_por',ierr)
  
        call fiReadDouble(string,grid%ran_fac,ierr)
        call fiDefaultMsg('ran_fac',ierr)
    
        call fiReadInt(string,grid%iread_perm,ierr)
        call fiDefaultMsg('iread_perm',ierr)
    
        call fiReadInt(string,grid%iread_geom,ierr)
        call fiDefaultMsg('iread_geom',ierr)


        if (grid%myrank == 0) &
          write(IUNIT2,'(/," *OPTS",/, &
            & "  write_init = ",3x,i2,/ &
            & "  iprint     = ",3x,i2,/, &
            & "  imod       = ",3x,i2,/, &
            & "  itecplot   = ",3x,i2,/, &
            & "  iblkfmt    = ",3x,i2,/, &
            & "  ndtcmx     = ",3x,i2,/, &
            & "  iran_por   = ",3x,i2,/, &
            & "  ran_fac    = ",3x,1pe12.4,/, &
            & "  iread_perm = ",3x,i2,/, &
            & "  iread_geom = ",3x,i2 &
            & )') grid%write_init,grid%iprint,grid%imod,grid%itecplot, &
            grid%iblkfmt,grid%ndtcmx,grid%iran_por,grid%ran_fac, &
            grid%iread_perm,grid%iread_geom

!....................

      case ('TOLR')

        call fiReadStringErrorMsg('TOLR',ierr)

        call fiReadInt(string,grid%stepmax,ierr)
        call fiDefaultMsg('stepmax',ierr)
  
        call fiReadInt(string,grid%iaccel,ierr)
        call fiDefaultMsg('iaccel',ierr)

        call fiReadInt(string,grid%newton_max,ierr)
        call fiDefaultMsg('newton_max',ierr)

        call fiReadInt(string,grid%icut_max,ierr)
        call fiDefaultMsg('icut_max',ierr)

        call fiReadDouble(string,grid%dpmxe,ierr)
        call fiDefaultMsg('dpmxe',ierr)

        call fiReadDouble(string,grid%dtmpmxe,ierr)
        call fiDefaultMsg('dtmpmxe',ierr)
  
        call fiReadDouble(string,grid%dcmxe,ierr)
        call fiDefaultMsg('dcmxe',ierr)

        call fiReadDouble(string,grid%dsmxe,ierr)
        call fiDefaultMsg('dsmxe',ierr)

        if (grid%myrank==0) write(IUNIT2,'(/," *TOLR ",/, &
          &"  flowsteps  = ",i6,/,      &
          &"  iaccel     = ",i3,/,      &
          &"  newtmx     = ",i3,/,      &
          &"  icutmx     = ",i3,/,      &
          &"  dpmxe      = ",1pe12.4,/, &
          &"  dtmpmxe    = ",1pe12.4,/, &
          &     "  dcmxe      = ",1pe12.4,/, &
          &"  dsmxe      = ",1pe12.4)') &
! For commented-out lines to work with the Sun f95 compiler, we have to 
! terminate the string in the line above; otherwise, the compiler tries to
! include the commented-out line as part of the continued string.
          grid%stepmax,grid%iaccel,grid%newton_max,grid%icut_max, &
          grid%dpmxe,grid%dtmpmxe,grid%dcmxe, grid%dsmxe

!....................

      case ('DXYZ')

        allocate(grid%dx0(grid%nx))
        allocate(grid%dy0(grid%ny))
        allocate(grid%dz0(grid%nz))
        
        call readxyz (grid%dx0,grid%nx)
        call readxyz (grid%dy0,grid%ny)
        call readxyz (grid%dz0,grid%nz)
    
        
        if (grid%myrank==0) then
          write(IUNIT2,'(/," *DXYZ ")')
          write(IUNIT2,'("  dx  ",/,(1p10e12.4))') (grid%dx0(i),i=1,grid%nx)
          write(IUNIT2,'("  dy  ",/,(1p10e12.4))') (grid%dy0(i),i=1,grid%ny)
          write(IUNIT2,'("  dz  ",/,(1p10e12.4))') (grid%dz0(i),i=1,grid%nz)
        endif

!....................


      case('RAD0')
    
        call fiReadDouble(string,grid%Radius_0,ierr)
        call fiDefaultMsg('R_0',ierr)


      case ('DIFF')

        call fiReadStringErrorMsg('DIFF',ierr)

        call fiReadDouble(string,grid%difaq,ierr)
        call fiDefaultMsg('difaq',ierr)

        call fiReadDouble(string,grid%delhaq,ierr)
        call fiDefaultMsg('delhaq',ierr)

        if (grid%myrank==0) write(IUNIT2,'(/," *DIFF ",/, &
          &"  difaq       = ",1pe12.4,"[m^2/s]",/, &
          &"  delhaq      = ",1pe12.4,"[kJ/mol]")') &
          grid%difaq,grid%delhaq

!....................

      case ('RCTR')

        call fiReadStringErrorMsg('RCTR',ierr)

        call fiReadInt(string,grid%ityprxn,ierr)
        call fiDefaultMsg('ityprxn',ierr)

        call fiReadDouble(string,grid%rk,ierr)
        call fiDefaultMsg('rk',ierr)

        call fiReadDouble(string,grid%phis0,ierr)
        call fiDefaultMsg('phis0',ierr)

        call fiReadDouble(string,grid%areas0,ierr)
        call fiDefaultMsg('areas0',ierr)

        call fiReadDouble(string,grid%pwrsrf,ierr)
        call fiDefaultMsg('pwrsrf',ierr)

        call fiReadDouble(string,grid%vbars,ierr)
        call fiDefaultMsg('vbars',ierr)

        call fiReadDouble(string,grid%ceq,ierr)
        call fiDefaultMsg('ceq',ierr)

        call fiReadDouble(string,grid%delHs,ierr)
        call fiDefaultMsg('delHs',ierr)

        call fiReadDouble(string,grid%delEs,ierr)
        call fiDefaultMsg('delEs',ierr)

        call fiReadDouble(string,grid%wfmts,ierr)
        call fiDefaultMsg('wfmts',ierr)

        if (grid%myrank == 0) &
        write(IUNIT2,'(/," *RCTR",/, &
          & "  ityp   = ",3x,i3,/, &
          & "  rk     = ",3x,1pe12.4," [mol/cm^2/s]",/, &
          & "  phis0  = ",3x,1pe12.4," [-]",/, &
          & "  areas0 = ",3x,1pe12.4," [1/cm]",/, &
          & "  pwrsrf = ",3x,1pe12.4," [-]",/, &
          & "  vbars  = ",3x,1pe12.4," [cm^3/mol]",/, &
          & "  ceq    = ",3x,1pe12.4," [mol/L]",/, &
          & "  delHs  = ",3x,1pe12.4," [J/kg]",/, &
          & "  delEs  = ",3x,1pe12.4," [J/kg]",/, &
          & "  wfmts  = ",3x,1pe12.4," [g/mol]" &
          & )') grid%ityprxn,grid%rk,grid%phis0,grid%areas0,grid%pwrsrf, &
          grid%vbars,grid%ceq,grid%delHs,grid%delEs,grid%wfmts

 ! convert: mol/cm^2 -> mol/cm^3 -> mol/dm^3 (note area 1/cm)          
        grid%rk = grid%rk * grid%areas0 * 1.d3
        grid%vbars = grid%vbars * 1.d-3 ! convert: cm^3/mol -> L/mol
      
        grid%delHs = grid%delHs * grid%wfmts * 1.d-3 ! convert kJ/kg -> kJ/mol
!        grid%delHs = grid%delHs * grid%scale ! convert J/kmol -> MJ/kmol

!....................

      case ('RADN')

        call fiReadStringErrorMsg('RADN',ierr)

        call fiReadDouble(string,grid%ret,ierr)
        call fiDefaultMsg('ret',ierr)

        call fiReadDouble(string,grid%fc,ierr)
        call fiDefaultMsg('fc',ierr)

        if (grid%myrank==0) write(IUNIT2,'(/," *RADN ",/, &
          &"  ret     = ",1pe12.4,/, &
          &"  fc      = ",1pe12.4)') &
          grid%ret,grid%fc

!....................


      case ('PHAR')

        call fiReadStringErrorMsg('PHAR',ierr)

        call fiReadDouble(string,grid%qu_kin,ierr)
        call fiDefaultMsg('TransReaction',ierr)
        if (grid%myrank==0) write(IUNIT2,'(/," *PHAR ",1pe12.4)')grid%qu_kin
        grid%yh2o_in_co2 = 0.d0
        if (grid%qu_kin > 0.d0) grid%yh2o_in_co2 = 1.d-2 ! check this number!
     
!......................

       case('RICH')
           call fiReadStringErrorMsg('PHAR',ierr)
           call fiReadDouble(string,grid%pref,ierr)
           call fiDefaultMsg('Ref. Pressure',ierr) 
     
!......................
      case ('HYDR')

        call fiReadStringErrorMsg('HYDR',ierr)

        call fiReadDouble(string,grid%dTdz,ierr)
        call fiDefaultMsg('dTdz',ierr)

        call fiReadDouble(string,grid%beta,ierr)
        call fiDefaultMsg('beta',ierr)

        call fiReadDouble(string,grid%tref,ierr)
        call fiDefaultMsg('tref',ierr)

        call fiReadDouble(string,grid%pref,ierr)
        call fiDefaultMsg('pref',ierr)

        call fiReadDouble(string,grid%conc0,ierr)
        call fiDefaultMsg('conc0',ierr)

        grid%ihydrostatic = 1
      
        if (grid%myrank==0) write(IUNIT2,'(/," *HYDR ",/, &
          &"  ihydro      = ",i3,/, &
          &"  dT/dz       = ",1pe12.4,/, &
          &"  beta        = ",1pe12.4,/, &
          &"  tref        = ",1pe12.4,/, &
          &"  pref        = ",1pe12.4,/, &
          &"  conc        = ",1pe12.4 &
          &)') &
          grid%ihydrostatic,grid%dTdz,grid%beta,grid%tref,grid%pref, &
          grid%conc0

!....................

      case ('SOLV')
    
        call fiReadStringErrorMsg('SOLV',ierr)

!       call fiReadDouble(string,eps,ierr)
!       call fiDefaultMsg('eps',ierr)

        call fiReadDouble(string,grid%atol,ierr)
        call fiDefaultMsg('atol_petsc',ierr)

        call fiReadDouble(string,grid%rtol,ierr)
        call fiDefaultMsg('rtol_petsc',ierr)

        call fiReadDouble(string,grid%stol,ierr)
        call fiDefaultMsg('stol_petsc',ierr)
      
        grid%dtol=1.D5
!       if (grid%use_ksp == 1) then
        call fiReadDouble(string,grid%dtol,ierr)
        call fiDefaultMsg('dtol_petsc',ierr)
!       endif
   
        call fiReadInt(string,grid%maxit,ierr)
        call fiDefaultMsg('maxit',ierr)
      
        call fiReadInt(string,grid%maxf,ierr)
        call fiDefaultMsg('maxf',ierr)

        if (grid%myrank==0) write(IUNIT2,'(/," *SOLV ",/, &
          &"  atol_petsc   = ",1pe12.4,/, &
          &"  rtol_petsc   = ",1pe12.4,/, &
          &"  stol_petsc   = ",1pe12.4,/, &
          &"  dtol_petsc   = ",1pe12.4,/, &
          &"  maxit        = ",8x,i5,/, &
          &"  maxf         = ",8x,i5 &
          &    )') &
           grid%atol,grid%rtol,grid%stol,grid%dtol,grid%maxit,grid%maxf

! The line below is a commented-out portion of the format string above.
! We have to put it here because of the stupid Sun compiler.
!    &"  eps          = ",1pe12.4,/, &

!....................

      case ('THRM')

        ireg = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('THRM',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ireg = ireg + 1
      
          call fiReadInt(string,idum,ierr)
          call fiDefaultMsg('idum',ierr)

          call fiReadDouble(string,grid%rock_density(ireg),ierr)
          call fiDefaultMsg('rock_density',ierr)

          call fiReadDouble(string,grid%cpr(ireg),ierr)
          call fiDefaultMsg('cpr',ierr)
        
          call fiReadDouble(string,grid%ckdry(ireg),ierr)
          call fiDefaultMsg('ckdry',ierr)
        
          call fiReadDouble(string,grid%ckwet(ireg),ierr)
          call fiDefaultMsg('ckwet',ierr)
        
          call fiReadDouble(string,grid%tau(ireg),ierr)
          call fiDefaultMsg('tau',ierr)

          call fiReadDouble(string,grid%cdiff(ireg),ierr)
          call fiDefaultMsg('cdiff',ierr)

          call fiReadDouble(string,grid%cexp(ireg),ierr)
          call fiDefaultMsg('cexp',ierr)

        !scale thermal properties
          grid%cpr(ireg) = grid%scale * grid%cpr(ireg)
          grid%dencpr(ireg) = grid%rock_density(ireg) * grid%cpr(ireg)
          grid%ckdry(ireg) = grid%scale * grid%ckdry(ireg)
          grid%ckwet(ireg) = grid%scale * grid%ckwet(ireg)
        enddo
      
        if (grid%myrank==0) then
          write(IUNIT2,'(/," *THRM: ",i3)') ireg
          write(IUNIT2,'("  itm rock_density  cpr        ckdry", &
            &                 "     ckwet       tau       cdiff     cexp")')
          write(IUNIT2,'("        [kg/m^3]  [J/kg/K]   [J/m/K/s]", &
            &              "     [J/m/K/s]     [-]        [m^2/s]       [-]")')
          do i = 1, ireg
            write(IUNIT2,'(i4,1p7e11.4)') i,grid%rock_density(i), &
            grid%cpr(i),grid%ckdry(i),grid%ckwet(i), &
            grid%tau(i),grid%cdiff(i),grid%cexp(i)
          enddo
        endif

!....................

      case ('PCKR')

        ireg = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('PCKR',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ireg = ireg + 1
       
          call fiReadInt(string,idum,ierr)
          call fiDefaultMsg('idum',ierr)
          
          call fiReadInt(string,grid%icaptype(idum),ierr)
          call fiDefaultMsg('icaptype',ierr)
      
          if (grid%use_mph == PETSC_TRUE .or. grid%use_owg == PETSC_TRUE &
              .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
              .or. grid%use_richard == PETSC_TRUE) then
            do np=1, grid%nphase
              call fiReadDouble(string,grid%sir(np,idum),ierr)
              call fiDefaultMsg('sir',ierr)
            enddo 
          else
            call fiReadDouble(string,grid%swir(idum),ierr)
            call fiDefaultMsg('swir',ierr)
          endif
        
          call fiReadDouble(string,grid%pckrm(idum),ierr)
          call fiDefaultMsg('lambda',ierr)
          grid%lambda(idum) = grid%pckrm(idum)/(-grid%pckrm(idum) +1.D0)
! Here the lambda is assigned as the same value of m

          call fiReadDouble(string,grid%alpha(idum),ierr)
          call fiDefaultMsg('alpha',ierr)

          call fiReadDouble(string,grid%pcwmax(idum),ierr)
          call fiDefaultMsg('pcwmax',ierr)
      
          call fiReadDouble(string,grid%pcbetac(idum),ierr)
          call fiDefaultMsg('pcbetac',ierr)
      
          call fiReadDouble(string,grid%pwrprm(idum),ierr)
          call fiDefaultMsg('pwrprm',ierr)

        enddo
      
        if (grid%myrank==0) then
          write(IUNIT2,'(/," *PCKR: ",i3)') ireg
          write(IUNIT2,'("  icp swir    lambda         alpha")')
          do j = 1, ireg
            i=grid%icaptype(j)
            if (grid%use_mph==PETSC_TRUE .or. grid%use_owg==PETSC_TRUE &
                 .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
                 .or. grid%use_richard == PETSC_TRUE) then
              write(IUNIT2,'(i4,1p8e12.4)') i,(grid%sir(np,i),np=1, &
                grid%nphase),grid%lambda(i),grid%alpha(i),grid%pcwmax(i), &
                grid%pcbetac(i),grid%pwrprm(i)
            else
              write(IUNIT2,'(i4,1p7e12.4)') i,grid%swir(i), &
                grid%lambda(i),grid%alpha(i),grid%pcwmax(i),grid%pcbetac(i), &
                grid%pwrprm(i)
            endif
          enddo
        end if

!....................
      
      case ('PHIK')

        ireg = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('PHIK',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ireg = ireg + 1
        
          if (ireg > MAXPERMREGIONS) then
            print *,'Error reading PHIK keyword: too many regions-stop',ireg
            stop
          endif
      
          call fiReadInt(string,grid%i1reg(ireg),ierr) 
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,grid%i2reg(ireg),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,grid%j1reg(ireg),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,grid%j2reg(ireg),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,grid%k1reg(ireg),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,grid%k2reg(ireg),ierr)
          call fiDefaultMsg('k2',ierr)
  
          call fiReadInt(string,grid%icap_reg(ireg),ierr)
          call fiDefaultMsg('icap',ierr)
  
          call fiReadInt(string,grid%ithrm_reg(ireg),ierr)
          call fiDefaultMsg('ithrm',ierr)
  
          call fiReadDouble(string,grid%por_reg(ireg),ierr)
          call fiDefaultMsg('por',ierr)
  
          call fiReadDouble(string,grid%tor_reg(ireg),ierr)
          call fiDefaultMsg('tor',ierr)
  
          call fiReadDouble(string,grid%perm_reg(ireg,1),ierr)
          call fiDefaultMsg('permx',ierr)
  
          call fiReadDouble(string,grid%perm_reg(ireg,2),ierr)
          call fiDefaultMsg('permy',ierr)
  
          call fiReadDouble(string,grid%perm_reg(ireg,3),ierr)
          call fiDefaultMsg('permz',ierr)
  
          call fiReadDouble(string,grid%perm_reg(ireg,4),ierr)
          call fiDefaultMsg('permpwr',ierr)

!          call fiReadDouble(string,grid%Perm_reg(ireg,5),ierr)
!          call fiDefaultMsg('porokin',ierr)

    
        enddo
        grid%iregperm = ireg
            
        if (grid%myrank==0) then
          write(IUNIT2,'(/," *PHIK: ireg = ",i4)') grid%iregperm
          write(IUNIT2,'("  i1  i2  j1  j2  k1  k2 icap ithrm  por      tor  &
            &",   "     permx      permy      permz [m^2]   permpwr")')
          do ireg = 1, grid%iregperm
            write(IUNIT2,'(6i4,2i4,1p6e11.4)') &
              grid%i1reg(ireg),grid%i2reg(ireg), &
              grid%j1reg(ireg),grid%j2reg(ireg), &
              grid%k1reg(ireg),grid%k2reg(ireg), &
              grid%icap_reg(ireg),grid%ithrm_reg(ireg), &
              grid%por_reg(ireg),grid%tor_reg(ireg), &
              (grid%perm_reg(ireg,i),i=1,4)
          enddo
        endif

!....................
      
      case ('INIT')
    
        call fiReadInt(string,grid%iread_init,ierr) 
        call fiDefaultMsg('iread_init',ierr)
      
        if (grid%myrank==0) then
          write(IUNIT2,'(/," *INIT: iread = ",i2)') grid%iread_init
        endif
      
        if (grid%iread_init == 0 .or. grid%iread_init == 2) then
      
          ireg = 0
          do
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('INIT',ierr)

            if (string(1:1) == '.' .or. string(1:1) == '/') exit
            ireg = ireg + 1

            call fiReadInt(string,grid%i1ini(ireg),ierr) 
            call fiDefaultMsg('i1',ierr)
            call fiReadInt(string,grid%i2ini(ireg),ierr)
            call fiDefaultMsg('i2',ierr)
            call fiReadInt(string,grid%j1ini(ireg),ierr)
            call fiDefaultMsg('j1',ierr)
            call fiReadInt(string,grid%j2ini(ireg),ierr)
            call fiDefaultMsg('j2',ierr)
            call fiReadInt(string,grid%k1ini(ireg),ierr)
            call fiDefaultMsg('k1',ierr)
            call fiReadInt(string,grid%k2ini(ireg),ierr)
            call fiDefaultMsg('k2',ierr)
         
            if (grid%use_mph==PETSC_TRUE .or. grid%use_owg==PETSC_TRUE &
                 .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
                 .or. grid%use_richard == PETSC_TRUE) then
              call fiReadInt(string,grid%iphas_ini(ireg),ierr)
              call fiDefaultMsg('iphase',ierr)
         
              do j=1,grid%ndof
                call fiReadDouble(string,grid%xx_ini(j,ireg),ierr)
                call fiDefaultMsg('xxini',ierr)
              enddo
            else
              call fiReadDouble(string,grid%pres_ini(ireg),ierr)
              call fiDefaultMsg('pres',ierr)
  
              call fiReadDouble(string,grid%temp_ini(ireg),ierr)
              call fiDefaultMsg('temp',ierr)
  
              call fiReadDouble(string,grid%sat_ini(ireg),ierr)
              call fiDefaultMsg('sat',ierr)
!              grid%sat_ini(ireg)=1.D0 - grid%sat_ini(ireg)
  
              call fiReadDouble(string,grid%conc_ini(ireg),ierr)
              call fiDefaultMsg('conc',ierr)
            endif
          enddo
      
          grid%iregini = ireg
      
          if (grid%myrank==0) then
            write(IUNIT2,'("  ireg = ",i4)') grid%iregini
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]   &
              &   ",    "sl [-]      c [mol/L]")')
            do ireg = 1, grid%iregini
              if (grid%use_mph==PETSC_TRUE .or. grid%use_owg==PETSC_TRUE &
                  .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
                  .or. grid%use_richard == PETSC_TRUE) then
                write(IUNIT2,'(7i4,1p10e12.4)') &
                  grid%i1ini(ireg),grid%i2ini(ireg), &
                  grid%j1ini(ireg),grid%j2ini(ireg), &
                  grid%k1ini(ireg),grid%k2ini(ireg), &
                  grid%iphas_ini(ireg),(grid%xx_ini(np,ireg),np =1,grid%ndof)
              else
                write(IUNIT2,'(6i4,1p10e12.4)') &
                  grid%i1ini(ireg),grid%i2ini(ireg), &
                  grid%j1ini(ireg),grid%j2ini(ireg), &
                  grid%k1ini(ireg),grid%k2ini(ireg), &
                  grid%pres_ini(ireg),grid%temp_ini(ireg),grid%sat_ini(ireg), &
                  grid%conc_ini(ireg)
              endif
            enddo
          endif

        else if (grid%iread_init == 1) then
    
!     read in initial conditions from file: pflow_init.dat
          if (grid%myrank == 0) then
            write(*,*) '--> read in initial conditions from file: &
                        &pflow_init.dat'
  
            open(IUNIT3, file='pflow_init.dat', action="read", status="old")

            ireg = 0
            do
              call fiReadFlotranString(IUNIT3,string,ierr)
!             call fiReadStringErrorMsg('INIT',ierr)

              if (string(1:1) == '.' .or. string(1:1) == '/') exit
              ireg = ireg + 1

              call fiReadInt(string,grid%i1ini(ireg),ierr) 
              call fiDefaultMsg('i1',ierr)
              call fiReadInt(string,grid%i2ini(ireg),ierr)
              call fiDefaultMsg('i2',ierr)
              call fiReadInt(string,grid%j1ini(ireg),ierr)
              call fiDefaultMsg('j1',ierr)
              call fiReadInt(string,grid%j2ini(ireg),ierr)
              call fiDefaultMsg('j2',ierr)
              call fiReadInt(string,grid%k1ini(ireg),ierr)
              call fiDefaultMsg('k1',ierr)
              call fiReadInt(string,grid%k2ini(ireg),ierr)
              call fiDefaultMsg('k2',ierr)
  
              if (grid%use_mph==PETSC_TRUE .or. grid%use_owg==PETSC_TRUE &
                  .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
                  .or. grid%use_richard == PETSC_TRUE) then
                call fiReadInt(string,grid%iphas_ini(ireg),ierr)
                call fiDefaultMsg('iphase_ini',ierr)
            
                do j=1,grid%ndof
                  call fiReadDouble(string,grid%xx_ini(j,ireg),ierr)
                  call fiDefaultMsg('xx_ini',ierr)
                enddo
              else
  
                call fiReadDouble(string,grid%pres_ini(ireg),ierr)
                call fiDefaultMsg('pres',ierr)
  
                call fiReadDouble(string,grid%temp_ini(ireg),ierr)
                call fiDefaultMsg('temp',ierr)
  
                call fiReadDouble(string,grid%sat_ini(ireg),ierr)
                call fiDefaultMsg('sat',ierr)

                call fiReadDouble(string,grid%conc_ini(ireg),ierr)
                call fiDefaultMsg('conc',ierr)
              endif
       
            enddo
            grid%iregini = ireg
            close(IUNIT3)
          endif
        endif

!....................

      case ('TIME')

        call fiReadStringErrorMsg('TIME',ierr)
      
        call fiReadWord(string,strtim,.false.,ierr)
      
        grid%tunit = strtim

        if (grid%tunit == 's') then
          grid%tconv = 1.d0
        else if (grid%tunit == 'm') then
          grid%tconv = 60.d0
        else if (grid%tunit == 'h') then
          grid%tconv = 60.d0 * 60.d0
        else if (grid%tunit == 'd') then
          grid%tconv = 60.d0 * 60.d0 * 24.d0
        else if (grid%tunit == 'mo') then
          grid%tconv = 60.d0 * 60.d0 * 24.d0 * 30.d0
        else if (grid%tunit == 'y') then
          grid%tconv = 60.d0 * 60.d0 * 24.d0 * 365.d0
        else
          if (grid%myrank == 0) then
            write(*,'(" Time unit: ",a3,/, &
              &" Error: time units must be one of ",/, &
              &"   s -seconds",/,"   m -minutes",/,"   h -hours",/, &
              &"   d -days", /, "  mo -months",/,"   y -years")') grid%tunit
          endif
          stop
        endif

        call fiReadInt(string,grid%kplot,ierr) 
        call fiDefaultMsg('kplot',ierr)
      
        allocate(grid%tplot(grid%kplot))
      
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('TIME',ierr)
        i2 = 0
        do
          i1 = i2 + 1
          i2 = i2+10
          if (i2 > grid%kplot) i2 = grid%kplot
          do i = i1, i2
            call fiReadDouble(string,grid%tplot(i),ierr)
            call fiDefaultMsg('tplot',ierr)
          enddo
          if (i2 == grid%kplot) exit
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('TIME',ierr)
        enddo

!       call fiReadFlotranString(IUNIT1,string,ierr)
!       call fiReadStringErrorMsg('TIME',ierr)
      
!       call fiReadDouble(string,grid%dt,ierr)
!       call fiDefaultMsg('dt',ierr)

!       call fiReadDouble(string,grid%dt_max,ierr)
!       call fiDefaultMsg('dt_max',ierr)

        if (grid%myrank==0) then
          write(IUNIT2,'(/," *TIME ",a3,1x,i4,/,(1p10e12.4))') grid%tunit, &
          grid%kplot,(grid%tplot(i),i=1,grid%kplot)
!         write(IUNIT2,'("  dt= ",1pe12.4,", dtmax= ",1pe12.4,/)') &
!         grid%dt,grid%dt_max
        endif
      
        ! convert time units to seconds
        do i = 1, grid%kplot
          grid%tplot(i) = grid%tconv * grid%tplot(i)
        enddo
!       grid%dt = grid%tconv * grid%dt
!       grid%dt_max = grid%tconv * grid%dt_max

!....................

      case ('DTST')

        call fiReadStringErrorMsg('DTST',ierr)
  
        call fiReadInt(string,grid%nstpmax,ierr)
        call fiDefaultMsg('nstpmax',ierr)
  
        allocate(grid%tstep(grid%nstpmax))
        allocate(grid%dtstep(grid%nstpmax))
  
        do i = 1, grid%nstpmax
          call fiReadDouble(string,grid%tstep(i),ierr)
          call fiDefaultMsg('tstep',ierr)
        enddo

        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('DTST',ierr)
        call fiReadDouble(string,grid%dt_min,ierr)
        call fiDefaultMsg('dt_min',ierr)
        do i = 1, grid%nstpmax
          call fiReadDouble(string,grid%dtstep(i),ierr)
          call fiDefaultMsg('dtstep',ierr)
        enddo
        
        grid%dt_max = grid%dtstep(1)
        
        grid%dt = grid%dt_min
      
        if (grid%myrank==0) then
          write(IUNIT2,'(/," *DTST ",i4,/," tstep= ",(1p10e12.4))')  &
            grid%nstpmax, (grid%tstep(i),i=1,grid%nstpmax)
          write(IUNIT2,'(" dtstep= ",1p10e12.4,/)') &
            grid%dt_min,(grid%dtstep(i),i=1,grid%nstpmax)
        endif
      
        ! convert time units to seconds
        do i = 1, grid%nstpmax
          grid%tstep(i) = grid%tconv * grid%tstep(i)
          grid%dtstep(i) = grid%tconv * grid%dtstep(i)
        enddo
        grid%dt = grid%tconv * grid%dt
        grid%dt_min = grid%tconv * grid%dt_min
        grid%dt_max = grid%tconv * grid%dt_max

!....................

      case ('BCON')

!-----------------------------------------------------------------------
!-----boundary conditions:  ibnd:  
!                   1-left,    2-right
!          3-top,    4-bottom
!          5-front,  6-back
!
!  ibndtyp:  1-Dirichlet         (p, T, C)
!  ibndtyp:  2-Neumann/Dirichlet (q, grad T=0, grad C=0)
!  ibndtyp:  3-Dirichlet/Neumann (p, grad T=0, grad C=0)
!-----------------------------------------------------------------------
        ibc = 0
        ir = 0
        grid%iregbc1(1) = 1
        do ! loop over blocks
        
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('BCON',ierr)
        
          if (string(1:1) == '.' .or. string(1:1) == '/') exit

          ibc = ibc + 1  ! RTM: Number of boundary conditions
        
          call fiReadInt(string,grid%ibndtyp(ibc),ierr)
          call fiDefaultMsg('ibndtyp',ierr)

          call fiReadInt(string,grid%iface(ibc),ierr)
          call fiDefaultMsg('iface',ierr)

          do ! loop over regions
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('BCON',ierr)
        
            if (string(1:1) == '.' .or. string(1:1) == '/') exit
            ir = ir + 1

            call fiReadInt(string,grid%i1bc(ir),ierr)
            call fiDefaultMsg('i1',ierr)
            call fiReadInt(string,grid%i2bc(ir),ierr)
            call fiDefaultMsg('i2',ierr)
            call fiReadInt(string,grid%j1bc(ir),ierr)
            call fiDefaultMsg('j1',ierr)
            call fiReadInt(string,grid%j2bc(ir),ierr)
            call fiDefaultMsg('j2',ierr)
            call fiReadInt(string,grid%k1bc(ir),ierr)
            call fiDefaultMsg('k1',ierr)
            call fiReadInt(string,grid%k2bc(ir),ierr)
            call fiDefaultMsg('k2',ierr)    

            ! Now read the velocities or pressures, depending on the BC type
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('BCON',ierr)
! due to this commented out loop, the below is left indented - geh
     !      do j=1, grid%nphase 
   
              if (grid%use_mph /= PETSC_TRUE .and.  &
                  grid%use_owg /= PETSC_TRUE .and. &
                  grid%use_vadose /= PETSC_TRUE .and. &
                  grid%use_flash /= PETSC_TRUE .and. &
                  grid%use_richard /= PETSC_TRUE) then  
                j=1
                if (grid%nphase>1) j=2
                if (grid%ibndtyp(ibc) == 1) then 
                  call fiReadDouble(string, grid%pressurebc0(j,ibc), ierr)
                  call fiDefaultMsg("Error reading pressure BCs:", ierr)
                else if (grid%ibndtyp(ibc) == 2) then
                  call fiReadDouble(string, grid%velocitybc0(j,ibc), ierr)
                  call fiDefaultMsg("Error reading velocity BCs:", ierr)
                else
                  call fiReadDouble(string, grid%pressurebc0(j,ibc), ierr)
                  call fiDefaultMsg("Error reading pressure BCs:", ierr)
                endif
                !enddo
                if (grid%nphase>1) grid%pressurebc0(1,ibc) = &
                  grid%pressurebc0(2,ibc)
                ! For simple input
                call fiReadDouble(string,grid%tempbc0(ibc),ierr)
                call fiDefaultMsg('tempbc',ierr)

                call fiReadDouble(string,grid%sgbc0(ibc),ierr)
                call fiDefaultMsg('sgbc',ierr)
                grid%sgbc0(ibc) = 1.D0 - grid%sgbc0(ibc) ! read in sl

                call fiReadDouble(string,grid%concbc0(ibc),ierr)
                call fiDefaultMsg('concbc',ierr)
      
              else
     
                call fiReadInt(string,grid%iphasebc0(ibc),ierr)
                call fiDefaultMsg('iphase',ierr)
       
                if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then 
                  do j=1,grid%ndof
                    call fiReadDouble(string,grid%xxbc0(j,ibc),ierr)
                    call fiDefaultMsg('xxbc',ierr)
                  enddo
                elseif (grid%ibndtyp(ibc) == 2) then
                  do j=1, grid%nphase       
                    call fiReadDouble(string, grid%velocitybc0(j,ibc), ierr)
                    call fiDefaultMsg("Error reading velocity BCs:", ierr)
                  enddo
                  do j=2,grid%ndof
                    call fiReadDouble(string,grid%xxbc0(j,ibc),ierr)
                    call fiDefaultMsg('xxbc',ierr)
                  enddo
                endif
            endif               
          enddo ! End loop over regions.
        
          grid%iregbc2(ibc) = ir
          if (ibc+1 > MAXBCBLOCKS) then
            write(*,*) 'Too many boundary condition blocks specified--stop: ', &
              ibc+1, MAXBCBLOCKS
            stop
          else
            grid%iregbc1(ibc+1) = grid%iregbc2(ibc)+1
          endif
        enddo ! End loop over blocks.
      
        grid%nblkbc = ibc
      
        if (grid%myrank == 0) then
          write(IUNIT2,'(/," *BCON: nblkbc = ",i4)') grid%nblkbc
          do ibc = 1, grid%nblkbc
            write(IUNIT2,'("  ibndtyp = ",i3," iface = ",i2)') &
              grid%ibndtyp(ibc), grid%iface(ibc)
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]  &
              &  c",     " [mol/L]")')
            do ireg = grid%iregbc1(ibc), grid%iregbc2(ibc)
              if (grid%use_mph== PETSC_TRUE .or. grid%use_owg== PETSC_TRUE &
                  .or. grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE&
                  .or. grid%use_richard == PETSC_TRUE) then
                if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
                  write(IUNIT2,'(7i4,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    grid%iphasebc0(ireg),(grid%xxbc0(j,ireg),j=1,grid%ndof)
                else if (grid%ibndtyp(ibc) == 2) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    (grid%velocitybc0(j,ireg),j=1,grid%nphase),&
                    (grid%xxbc0(j,ireg),j=2,grid%ndof)
                endif
              else
                if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    (grid%pressurebc0(j,ireg),j=1,grid%nphase), &
                    grid%tempbc0(ireg), &
                    grid%concbc0(ireg)
                else if (grid%ibndtyp(ibc) == 2) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    (grid%velocitybc0(j,ireg),j=1,grid%nphase), &
                    grid%tempbc0(ireg), &
                    grid%concbc0(ireg)
                endif
              endif
            enddo
          enddo
        endif

!....................

      case ('SOUR')

        isrc = 0
        ir = 0
      
        do ! loop over sources
      
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('SOUR',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit

          isrc = isrc + 1  ! Number of sources

          ir = ir + 1

          call fiReadInt(string,grid%i1src(ir),ierr)
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,grid%i2src(ir),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,grid%j1src(ir),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,grid%j2src(ir),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,grid%k1src(ir),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,grid%k2src(ir),ierr)
          call fiDefaultMsg('k2',ierr)    
          print *,'Source', isrc, ir   
          ! Read time, temperature, q-source
          i = 0
          do ! loop over time intervals
        
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('SOUR',ierr)
            if (string(1:1) == '.' .or. string(1:1) == '/') exit
        
            i = i + 1
        
            if (i+1 > 10) then
              write(*,*) 'Too many times specified in SOURce--stop: ', i+1, 10
              stop
            endif
        
            call fiReadDouble(string,grid%timesrc(i,isrc),ierr)
            call fiDefaultMsg('timesrc',ierr)

            call fiReadDouble(string,grid%tempsrc(i,isrc),ierr)
            call fiDefaultMsg('tempsrc',ierr)
      
            call fiReadDouble(string,grid%qsrc(i,isrc),ierr)
            call fiDefaultMsg('qsrc',ierr)
      
            call fiReadDouble(string,grid%csrc(i,isrc),ierr)
            call fiDefaultMsg('csrc',ierr)

          enddo ! End loop over time.

          grid%ntimsrc = i

          if (grid%ntimsrc > MAXSRCTIMES) then
            write(*,*) 'Too many source times specified--stop: ', &
              grid%ntimsrc, MAXSRCTIMES
            stop
          endif

          if (isrc+1 > MAXSRC) then
            write(*,*) 'Too many source blocks specified--stop: ', &
            isrc+1, MAXSRC
            stop
          endif
        enddo ! End loop over sources.

        grid%nblksrc = isrc

        if (grid%myrank == 0) then
          write(IUNIT2,'(/," *SOURce: nblksrc = ",i4)') grid%nblksrc
          do isrc = 1, grid%nblksrc
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2")')
            write(IUNIT2,'(6i4)') &
              grid%i1src(isrc),grid%i2src(isrc), &
              grid%j1src(isrc),grid%j2src(isrc), &
              grid%k1src(isrc),grid%k2src(isrc)
            write(IUNIT2,'("    t [s]        T [C]    QH2O [kg/s]    &
              &QCO2 [kg/s]")')
            do ir = 1, grid%ntimsrc
              write(IUNIT2,'(1p10e12.4)') &
                grid%timesrc(ir,isrc),grid%tempsrc(ir,isrc), &
                grid%qsrc(ir,isrc), &
                grid%csrc(ir,isrc)
            enddo
          enddo
        endif

!....................
      
      case ('BRK')

        ibrk = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('BRK',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ibrk = ibrk + 1
      
          call fiReadInt(string,grid%i1brk(ibrk),ierr) 
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,grid%i2brk(ibrk),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,grid%j1brk(ibrk),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,grid%j2brk(ibrk),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,grid%k1brk(ibrk),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,grid%k2brk(ibrk),ierr)
          call fiDefaultMsg('k2',ierr)

          call fiReadInt(string,grid%ibrktyp(ibrk),ierr)
          call fiDefaultMsg('ibrktyp',ierr)

          call fiReadInt(string,grid%ibrkface(ibrk),ierr)
          call fiDefaultMsg('ibrkface',ierr)
        enddo
        grid%ibrkcrv = ibrk
            
        if (grid%myrank==0) then
          write(IUNIT2,'(/," *BRK: ibrk = ",i4)') grid%ibrkcrv
          write(IUNIT2,'("  i1  i2  j1  j2  k1  k2  ibrktyp  ibrkface  ")')
          do ibrk = 1, grid%ibrkcrv
            write(IUNIT2,'(6i4,4x,i2,7x,i2)') &
              grid%i1brk(ibrk),grid%i2brk(ibrk), &
              grid%j1brk(ibrk),grid%j2brk(ibrk), &
              grid%k1brk(ibrk),grid%k2brk(ibrk),grid%ibrktyp(ibrk), &
              grid%ibrkface(ibrk)
          enddo
        endif

        if (grid%ndof == 1) grid%ibrkcrv = 0

!....................

      case default
    
        if (grid%myrank == 0) then
          print *, "Error reading input file: keyword not found. Terminating."
        endif
        call PetscFinalize(ierr)
        stop

    end select

  enddo

  close(IUNIT1)
  
end subroutine pflowgrid_read_input

!======================================================================

subroutine readxyz (a,n)

  use fileio_module
  
  implicit none
  
  integer*4, intent(in) :: n
  integer*4 :: i, i1, i2, m
  integer ::  ierr, nvalue=10
  real*8, intent(inout) :: a(*)
#include "definitions.h"
  character(len=MAXSTRINGLENGTH) :: string 

  save nvalue

!  call fiReadStringErrorMsg('DXYZ',ierr)

!  call fiReadDouble(string,grid%radius_0,ierr)
!  call fiDefaultMsg('radius_0',ierr)

  i2 = 0
  do
    i1 = i2+1
    i2 = i2+nvalue
    if (i2.gt.n) i2 = n
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg('DXYZ',ierr)
    do i = i1, i2
      call fiReadDouble(string, a(i), ierr)
      if (ierr .eq. 1) a(i) = 0.d0
!     print *,i,i1,i2,nvalue,a(i),n,ierr
!     call fiDefaultMsg("Error reading grid spacing", ierr)
    enddo
    do i = i1,i2
      if (a(i).eq.0.d0) then

!---------if less than nx non-zero values are read, set all the zero
!         values to the last non zero value read. Only for cartesian 
!         system

        do m = i,n
          a(m) = a(i-1)
        enddo
        return
      endif
    enddo
    if (i2.ge.n) exit
  enddo
    
end subroutine readxyz

!======================================================================

!#include "pflowgrid_parse_cmdline.F90"

subroutine pflowGrid_parse_cmdline(grid)
  
  implicit none

  type(pflowGrid), intent(in) :: grid

  integer :: ierr

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-snes_mf", & 
    grid%use_matrix_free, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_analytical", &
    grid%use_analytical, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-print_hhistory", &
    grid%print_hhistory, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-monitor_h", &
    grid%monitor_h, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_liquid", &
    grid%use_liquid, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_cond", &
    grid%use_cond, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_th", &
    grid%use_th, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_thc", &
    grid%use_thc, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_2ph", &
    grid%use_2ph, ierr)
   call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_mph", &
    grid%use_mph, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_richard", &
    grid%use_richard, ierr)
   call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_flash", &
    grid%use_flash, ierr)

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_debug", &
    grid%use_debug, ierr)
 ! debug option

end subroutine pflowGrid_parse_cmdline

!======================================================================


!#include "pflowgrid_misc.F90"

real*8 function pflowGrid_get_t(grid)
  
  implicit none

  type(pflowGrid), intent(in) :: grid

  pflowGrid_get_t = grid%t

end function pflowGrid_get_t

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
  
  call MatMFFDGetH(grid%J, h, ierr)

  call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)

  if (myrank == 0) then
    write(*,*) "#At SNES iteration ", its, "h is ", h
  endif

end subroutine pflowgrid_MonitorH

!======================================================================

end module pflow_grid_module
