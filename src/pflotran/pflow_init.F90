module Init_module

  implicit none

  private

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

  public :: initPFLOW
 
contains

! ************************************************************************** !
!
! initPFLOW: Initializes a pflow grid object
! author: Glenn Hammond
! date: 10/23/07
!
! **************************************************************************
subroutine initPFLOW(simulation,filename)

  use Simulation_module
  use Option_module
  use Grid_module
  use Solver_module
  use Solution_module
  use Material_module
  use Solver_module
  use Timestepper_module

!  use MPHASE_module
  use Richards_module
  use pflow_convergence_module
  use pflow_solv_module 
    
  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXWORDLENGTH) :: filename

  type(solver_type), pointer :: solver
  type(solution_type), pointer :: solution
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  ISColoring :: iscoloring

  integer :: nphase, ndof, nspec, npricomp
  integer :: mcomp, mphas, idcdm, itable
  integer :: i
  
  real*8 :: alpha, maxstep, steptol
  
  PetscErrorCode :: ierr
  
  solver => simulation%stepper%solver
  solution => simulation%solution
  option => solution%option
  grid => solution%grid
  
  call readSelectCardsFromInput(solution,filename,mcomp,mphas)
  
  option%use_ksp=PETSC_FALSE
  option%use_isoth=PETSC_FALSE
   
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-snes_mf", & 
                           option%use_matrix_free, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_analytical", &
                           option%use_analytical, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-print_hhistory", &
                           option%print_hhistory, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-monitor_h", &
                           option%monitor_h, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_liquid", &
                           option%use_liquid, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_cond", &
                           option%use_cond, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_th", &
                           option%use_th, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_thc", &
                           option%use_thc, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_2ph", &
                           option%use_2ph, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_mph", &
                           option%use_mph, ierr)
   call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_flash", &
                           option%use_flash, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_owg", &
                           option%use_owg, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_vadose", &
                           option%use_vadose, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_richards", &
                           option%use_richards, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_debug", &
                           option%use_debug, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_ksp", &
                           option%use_ksp, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_isoth", &
                           option%use_isoth, ierr)

  option%nphase = nphase
  option%ndof = ndof
  option%nspec = nspec
  option%npricomp = npricomp

if( option%use_liquid ==PETSC_FALSE .and.&
    option%use_cond==PETSC_FALSE  .and.&
    option%use_th==PETSC_FALSE  .and.&
    option%use_thc==PETSC_FALSE  .and.&
    option%use_2ph==PETSC_FALSE  .and.&
    option%use_mph==PETSC_FALSE  .and.&
    option%use_flash==PETSC_FALSE  .and.&
    option%use_owg==PETSC_FALSE  .and.&
    option%use_vadose==PETSC_FALSE  .and.&
    option%use_richards==PETSC_FALSE) then

 if(mcomp >0 .and. mphas>0)then
   if( option%use_liquid == PETSC_FALSE .and. mcomp ==1)then
       option%use_liquid = PETSC_TRUE
       option%nphase = 1; option%ndof =1
   endif
   if( option%use_cond == PETSC_FALSE .and. mcomp == 32)then
       option%use_cond = PETSC_TRUE
       option%nphase = 1; option%ndof =1
   endif
   if( option%use_th == PETSC_FALSE .and. mcomp == 33 .and. mphas == 3)then
       option%use_th = PETSC_TRUE
       option%nphase = 1; option%ndof =2
   endif
   if( option%use_thc == PETSC_FALSE .and. mcomp == 37)then
       option%use_thc = PETSC_TRUE
       option%nphase = 1; option%ndof =3
   endif
   if( option%use_mph == PETSC_FALSE .and. mcomp == 35)then
       option%use_mph = PETSC_TRUE
       option%nphase = 2; option%ndof =3; option%nspec =2 
   endif
   if( option%use_vadose == PETSC_FALSE .and. mcomp == 49)then
       option%use_vadose = PETSC_TRUE
       option%nphase = 2; option%ndof =3; option%nspec =2 
   endif
   if( option%use_richards == PETSC_FALSE .and. mcomp == 33 .and. mphas == 11)then
       option%use_richards = PETSC_TRUE
       option%nphase = 1; option%ndof = 2;  option%nspec =1
       if(nspec > 1) then
           option%ndof = nspec +1 ;  option%nspec = nspec
       endif
   endif
    if( option%use_owg == PETSC_FALSE .and. mcomp == 11)then
       option%use_owg = PETSC_TRUE
       option%nphase = 3; option%ndof =3; option%nspec =3 
   endif
  endif
 if( option%use_liquid ==PETSC_FALSE .and.&
    option%use_cond==PETSC_FALSE  .and.&
    option%use_th==PETSC_FALSE  .and.&
    option%use_thc==PETSC_FALSE  .and.&
    option%use_2ph==PETSC_FALSE  .and.&
    option%use_mph==PETSC_FALSE  .and.&
    option%use_flash==PETSC_FALSE  .and.&
    option%use_owg==PETSC_FALSE  .and.&
    option%use_vadose==PETSC_FALSE  .and.&
    option%use_richards==PETSC_FALSE) then
     print *,'No method determined, stop:'
     stop
 endif       
                         
endif  

! hardwire to uncoupled for now
!  if (icouple == 0) then
    option%using_pflowGrid = PETSC_FALSE
!  else
!    option%using_pflowGrid = PETSC_TRUE
!  endif

  option%idcdm = idcdm
      
  option%itable = itable
      
  !set specific phase indices
  option%jh2o = 1; option%jgas =1
  select case(option%nphase)
    case(2)
      if (option%use_2ph == PETSC_TRUE) then
        option%jgas=  2
        option%jco2 = 3
      else
        option%jco2 = 2
        option%jgas =3 
      endif
    case(3)
      option%jco2 = 2
      option%jgas =3 
  end select 

!-----------------------------------------------------------------------
      ! Initialize some counter variables.
!-----------------------------------------------------------------------
  option%t = 0.d0
  option%flowsteps = 0
  option%newtcum = 0
  option%icutcum = 0

      !-----------------------------------------------------------------------
      ! Initialize some parameters to sensible values.  These are parameters
      ! which should be set via the command line or the input file, but it
      ! seems good practice to set them to sensible values when a pflowGrid
      ! is created.
!-----------------------------------------------------------------------
  allocate(option%tfac(13))
      
  option%tfac(1)  = 2.0d0; option%tfac(2)  = 2.0d0
  option%tfac(3)  = 2.0d0; option%tfac(4)  = 2.0d0
  option%tfac(5)  = 2.0d0; option%tfac(6)  = 1.8d0
  option%tfac(7)  = 1.6d0; option%tfac(8)  = 1.4d0
  option%tfac(9)  = 1.2d0; option%tfac(10) = 1.0d0
  option%tfac(11) = 1.0d0; option%tfac(12) = 1.0d0
  option%tfac(13) = 1.0d0      
  option%use_matrix_free = 1
  option%newton_max = 16
  option%icut_max = 16
  option%iaccel = 1
  option%dt = 1.d0
  option%dt_min = 1.d0
  option%dt_max = 3.1536d6 ! One-tenth of a year.
  solver%atol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%rtol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%stol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%maxit = PETSC_DEFAULT_INTEGER
  solver%maxf = PETSC_DEFAULT_INTEGER
  solver%idt_switch = 0
  
  option%ihydrostatic = 0
  option%conc0 = 1.d-6
  
  option%dpmxe = 5.d4
  option%dtmpmxe = 2.d0
  option%dsmxe = 5.d0
  option%dcmxe = 5.d0

  !physical constants and defult variables
  option%difaq = 1.d-9 ! m^2/s read from input file
  option%delhaq = 12.6d0 ! kJ/mol read from input file
  option%gravity = 9.8068d0    ! m/s^2
  option%tref   = 50.D0
  option%fmwh2o = 18.01534d0 ! kg H2O/mol H2O
  option%fmwco2 = 44.0098d0
  option%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  allocate(option%steady_eps(ndof))
  option%steady_eps = -1.D0

 ! initialize default values  
  option%m_nacl =0.D0  ! default brine concentration
  
 ! default output variables
  if(option%use_2ph == PETSC_TRUE .or. option%use_mph == PETSC_TRUE &
              .or. option%use_vadose == PETSC_TRUE &
              .or. option%use_flash == PETSC_TRUE )then
    option%var_plot_num = 11
    allocate(option%var_plot_nm(11))
    option%var_plot_nm(1) = 'x'
    option%var_plot_nm(2) = 'y'
    option%var_plot_nm(3) = 'z'
    option%var_plot_nm(4) = 'iphase'
    option%var_plot_nm(5) = 'pl'
    option%var_plot_nm(6) = 'pg'
    option%var_plot_nm(7) = 'temp'
    option%var_plot_nm(8) = 'Sg'
    option%var_plot_nm(9) = 'Xg_Aq'
    option%var_plot_nm(10) = 'Xg_G'
    option%var_plot_nm(11) = 'Vf'
    if(option%use_mph == PETSC_TRUE &
              .or. option%use_vadose == PETSC_TRUE &
              .or. option%use_flash == PETSC_TRUE) then 
      allocate(option%var_plot_ind(5:10))  
      option%var_plot_ind(5)= 2
      option%var_plot_ind(6)= -5
      option%var_plot_ind(7)= 1
         
    endif
  else 
    option%var_plot_num = 8
    allocate(option%var_plot_nm(8))
    option%var_plot_nm(1) = 'x'
    option%var_plot_nm(2) = 'y'
    option%var_plot_nm(3) = 'z'
    option%var_plot_nm(4) = 'p'
    option%var_plot_nm(5) = 'T'
    option%var_plot_nm(6) = 'Sl'
    option%var_plot_nm(7) = 'Conc'
    option%var_plot_nm(8) = 'Vf'
 endif 

  
  call createDMs(grid,option)
  
   !-----------------------------------------------------------------------
 ! Create the vectors with parallel layout corresponding to the DA's,
 ! and, for vectors that need to be ghosted, create the corresponding
 ! ghosted vectors.
 !-----------------------------------------------------------------------

  ! 1 degree of freedom
  call createPetscVector(grid,ONEDOF,option%porosity,GLOBAL)
  call VecDuplicate(option%porosity, option%porosity0, ierr)
  call VecDuplicate(option%porosity, option%tor, ierr)
  call VecDuplicate(option%porosity, option%conc, ierr)
  call VecDuplicate(option%porosity, grid%structured_grid%dx, ierr)
  call VecDuplicate(option%porosity, grid%structured_grid%dy, ierr)
  call VecDuplicate(option%porosity, grid%structured_grid%dz, ierr)
  call VecDuplicate(option%porosity, option%volume, ierr)
  call VecDuplicate(option%porosity, option%ithrm, ierr)
  call VecDuplicate(option%porosity, option%icap, ierr)
  call VecDuplicate(option%porosity, option%iphas, ierr)
  call VecDuplicate(option%porosity, option%iphas_old, ierr)
  call VecDuplicate(option%porosity, option%temp, ierr)
  call VecDuplicate(option%porosity, option%ttemp, ierr)
  call VecDuplicate(option%porosity, option%phis, ierr)

  call VecDuplicate(option%porosity, option%perm_xx, ierr)
  call VecDuplicate(option%porosity, option%perm_yy, ierr)
  call VecDuplicate(option%porosity, option%perm_zz, ierr)
  call VecDuplicate(option%porosity, option%perm0_xx, ierr)
  call VecDuplicate(option%porosity, option%perm0_yy, ierr)
  call VecDuplicate(option%porosity, option%perm0_zz, ierr)
  call VecDuplicate(option%porosity, option%perm_pow, ierr)
      
  call createPetscVector(grid,ONEDOF,grid%structured_grid%dx_loc,LOCAL)
  call VecDuplicate(grid%structured_grid%dx_loc, grid%structured_grid%dy_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, grid%structured_grid%dz_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, option%porosity_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, option%tor_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, option%ithrm_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, option%icap_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, option%iphas_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, option%ttemp_loc, ierr)

  call VecDuplicate(grid%structured_grid%dx_loc, option%perm_xx_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, option%perm_yy_loc, ierr)
  call VecDuplicate(grid%structured_grid%dx_loc, option%perm_zz_loc, ierr)

  ! 3 degrees of freedom
! call DACreateGlobalVector(option%da_3_dof, option%perm, ierr)
! call DACreateLocalVector(option%da_3_dof, option%perm_loc, ierr)

  ! nphase degrees of freedom
  call createPetscVector(grid,NPHASEDOF,option%pressure,GLOBAL)
  call VecDuplicate(option%pressure, option%ppressure, ierr)
  call VecDuplicate(option%pressure, option%ssat, ierr)
  call VecDuplicate(option%pressure, option%sat, ierr)
  call VecDuplicate(option%pressure, option%dp, ierr)
  call VecDuplicate(option%pressure, option%density, ierr)
  call VecDuplicate(option%pressure, option%ddensity, ierr)
  call VecDuplicate(option%pressure, option%avgmw, ierr)
  call VecDuplicate(option%pressure, option%d_p, ierr)
  call VecDuplicate(option%pressure, option%d_t, ierr)
  call VecDuplicate(option%pressure, option%h, ierr)
  call VecDuplicate(option%pressure, option%hh, ierr)
  call VecDuplicate(option%pressure, option%h_p, ierr)
  call VecDuplicate(option%pressure, option%h_t, ierr)
  call VecDuplicate(option%pressure, option%viscosity, ierr)
  call VecDuplicate(option%pressure, option%v_p, ierr)
  call VecDuplicate(option%pressure, option%v_t, ierr)
 ! xmol may not be nphase DOF, need change later 
  call VecDuplicate(option%pressure, option%xxmol, ierr)
  call VecDuplicate(option%pressure, option%xmol, ierr)


  call createPetscVector(grid,NPHASEDOF, option%ppressure_loc, LOCAL)
  call VecDuplicate(option%ppressure_loc, option%ssat_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%xxmol_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%ddensity_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%avgmw_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%d_p_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%d_t_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%hh_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%h_p_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%h_t_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%viscosity_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%v_p_loc, ierr)
  call VecDuplicate(option%ppressure_loc, option%v_t_loc, ierr)
      
  ! 3 * nphase degrees of freedom (velocity vector)
  call createPetscVector(grid,THREENPDOF, option%vl, GLOBAL)
      
! print *,'pflowgrid_new: ',option%using_pflowGrid
      
  if (option%using_pflowGrid == PETSC_TRUE) &
    call VecDuplicate(option%vl, option%vvl, ierr)
      
  ! nvar * nphase degrees of freedom
!  call DACreateGlobalVector(option%da_ncnp_dof, option%xmol, ierr)


  if (option%use_2ph == PETSC_TRUE) then
    print *,'2ph add var'
    call VecDuplicate(option%sat, option%d_s, ierr)
    call VecDuplicate(option%sat, option%h_s, ierr)
    call VecDuplicate(option%sat, option%u, ierr)
    call VecDuplicate(option%sat, option%uu, ierr)
    call VecDuplicate(option%sat, option%u_s, ierr)
    call VecDuplicate(option%sat, option%u_p, ierr)
    call VecDuplicate(option%sat, option%u_t, ierr)

    call VecDuplicate(option%ssat_loc, option%d_s_loc, ierr)
    call VecDuplicate(option%ssat_loc, option%h_s_loc, ierr)
    call VecDuplicate(option%ssat_loc, option%u_s_loc, ierr)

    call VecDuplicate(option%sat, option%pcw, ierr)
    call VecDuplicate(option%sat, option%pc_p, ierr)
    call VecDuplicate(option%sat, option%pc_t, ierr)
    call VecDuplicate(option%sat, option%pc_s, ierr)
    call VecDuplicate(option%ssat_loc, option%pcw_loc, ierr)
    call VecDuplicate(option%ssat_loc, option%pc_p_loc, ierr)
    call VecDuplicate(option%ssat_loc, option%pc_t_loc, ierr)
    call VecDuplicate(option%ssat_loc, option%pc_s_loc, ierr)


    call VecDuplicate(option%sat, option%kvr, ierr)
    call VecDuplicate(option%sat, option%kvr_p, ierr)
    call VecDuplicate(option%sat, option%kvr_t, ierr)
    call VecDuplicate(option%sat, option%kvr_s, ierr)
    call VecDuplicate(option%ssat_loc, option%kvr_loc, ierr)
    call VecDuplicate(option%ssat_loc, option%kvr_p_loc, ierr)
    call VecDuplicate(option%ssat_loc, option%kvr_t_loc, ierr)
    call VecDuplicate(option%ssat_loc, option%kvr_s_loc, ierr)

    call createPetscVector(grid,NPHANCOMPDOF, option%h_c,GLOBAL)
    call VecDuplicate(option%h_c, option%u_c, ierr)
    call VecDuplicate(option%h_c, option%avgmw_c, ierr)
    call VecDuplicate(option%h_c, option%d_c, ierr)
    call VecDuplicate(option%h_c, option%pc_c, ierr)
    call VecDuplicate(option%h_c, option%kvr_c, ierr)
        
    call createPetscVector(grid,NPHANCOMPDOF, option%h_c_loc,LOCAL)
    call VecDuplicate(option%h_c_loc, option%avgmw_c_loc, ierr)
    call VecDuplicate(option%h_c_loc, option%d_c_loc, ierr)
    call VecDuplicate(option%h_c_loc, option%pc_c_loc, ierr)
    call VecDuplicate(option%h_c_loc, option%kvr_c_loc, ierr)


    call createPetscVector(grid,NPHANSPECDOF, option%hen,GLOBAL)
    call createPetscVector(grid,NPHANSPECDOF, option%hen_loc,LOCAL)
    call VecDuplicate(option%hen, option%hen_p, ierr)
    call VecDuplicate(option%hen_loc, option%hen_p_loc, ierr)
    call VecDuplicate(option%hen, option%hen_t, ierr)
    call VecDuplicate(option%hen_loc, option%hen_t_loc, ierr)
    call VecDuplicate(option%hen, option%hen_s, ierr)
    call VecDuplicate(option%hen_loc, option%hen_s_loc, ierr)
    call VecDuplicate(option%hen, option%df, ierr)
    call VecDuplicate(option%hen_loc, option%df_loc, ierr)
    call VecDuplicate(option%df, option%df_p, ierr)
    call VecDuplicate(option%df_loc, option%df_p_loc, ierr)
    call VecDuplicate(option%df, option%df_t, ierr)
    call VecDuplicate(option%df_loc, option%df_t_loc, ierr)
    call VecDuplicate(option%df, option%df_s, ierr)
    call VecDuplicate(option%df_loc, option%df_s_loc, ierr)
  
    call createPetscVector(grid,NPHANSPECNCOMPDOF, option%hen_c, GLOBAL)
    call createPetscVector(grid,NPHANSPECNCOMPDOF, option%hen_c_loc,LOCAL)
      
    call VecDuplicate(option%hen_c, option%df_c, ierr)
    call VecDuplicate(option%hen_c_loc, option%df_c_loc, ierr)
  end if  

  if (option%use_mph == PETSC_TRUE) then 
    call createPetscVector(grid,VARDOF, option%var,GLOBAL)
    call createPetscVector(grid,VARDOF, option%var_loc,LOCAL)
  endif

  if (option%use_richards == PETSC_TRUE) then
    call createPetscVector(grid,VARDOF, option%var,GLOBAL)
    call createPetscVector(grid,VARDOF, option%var_loc, LOCAL)
  endif
 
 if (option%use_flash == PETSC_TRUE) then
    call createPetscVector(grid,VARDOF, option%var, GLOBAL)
    call createPetscVector(grid,VARDOF, option%var_loc, LOCAL)
  endif
     
  if (option%use_owg == PETSC_TRUE) then
    call createPetscVector(grid,VARDOF, option%var, GLOBAL)
    call createPetscVector(grid,VARDOF, option%var_loc, LOCAL)
  endif  
  if (option%use_vadose == PETSC_TRUE) then
    call createPetscVector(grid,VARDOF, option%var, GLOBAL)
    call createPetscVector(grid,VARDOF, option%var_loc, LOCAL)
  endif



      ! ndof degrees of freedom
  call createPetscVector(grid,NDOF, option%xx, GLOBAL)
  call VecDuplicate(option%xx, option%yy, ierr)
  call VecDuplicate(option%xx, option%dxx, ierr)
  call VecDuplicate(option%xx, option%r, ierr)
  call VecDuplicate(option%xx, option%accum, ierr)

     
  call VecSetBlocksize(option%dxx, option%ndof, ierr)

  call createPetscVector(grid,NDOF, option%xx_loc, LOCAL)
  
  ! Create Natural Vec for output: use VecDuplicate here?
!  call DACreateNaturalVector(option%da_1_dof,      option%c_nat,    ierr)
!  call DACreateNaturalVector(option%da_1_dof,      option%phis_nat, ierr)
!  call DACreateNaturalVector(option%da_1_dof,      option%t_nat,    ierr)
!  call DACreateNaturalVector(option%da_1_dof,      option%por_nat,  ierr)
!  call DACreateNaturalVector(option%da_nphase_dof, option%p_nat,    ierr)
!  call DACreateNaturalVector(option%da_nphase_dof, option%s_nat,    ierr)
!  call DACreateNaturalVector(option%da_3np_dof,    option%vl_nat,   ierr)

!   if (option%use_2ph == PETSC_TRUE) &
!     call DACreateNaturalVector(option%da_nphase_dof, option%x_nat, ierr)
!-----------------------------------------------------------------------
  ! Set up PETSc nonlinear solver context.
!-----------------------------------------------------------------------
  call SNESCreate(PETSC_COMM_WORLD, option%snes, ierr)
  CHKERRQ(ierr)
!-----------------------------------------------------------------------
  ! Set up indexing of grid ids (local to global, global to local, etc
!-----------------------------------------------------------------------
  call mapGridIndices(grid)

  option%nldof = grid%nlmax * option%nphase
  option%ngdof = grid%ngmax * option%nphase

!-----------------------------------------------------------------------
      ! Allocate memory for allocatable arrays.
!-----------------------------------------------------------------------
  i = grid%internal_connection_list%first%num_connections
  allocate(option%vl_loc(i))
  allocate(option%vvl_loc(i))
  allocate(option%vg_loc(i))
  allocate(option%vvg_loc(i))
    
    
  option%vl_loc = 0.D0
  option%vg_loc = 0.D0
      
  allocate(option%xphi_co2(grid%nlmax))
  allocate(option%xxphi_co2(grid%nlmax))
  allocate(option%den_co2(grid%nlmax))
  allocate(option%dden_co2(grid%nlmax))

  option%xphi_co2 = 1.d0
  option%xxphi_co2 = 1.d0
  option%den_co2 = 1.d0
  option%dden_co2 = 1.d0

! if (grid%using_pflowGrid == PETSC_TRUE) &
! allocate(grid%vvl_loc(grid%nconn*grid%nphase))

  ! I don't like having a fixed number of boundary condition regions.
  ! Memory for these arrays ought to allocated by parsing the input file
  ! to determine the number of regions.  This is the lazy way... I 
  ! should fix it eventually.
  ! The same goes for the number of BC blocks.
  allocate(option%iregbc1(MAXBCREGIONS))
  allocate(option%iregbc2(MAXBCREGIONS))
  allocate(option%ibndtyp(MAXBCREGIONS))
!GEH - Structured Grid Dependence - Begin
  allocate(option%iface(MAXBCREGIONS))
  allocate(option%k1bc(MAXBCBLOCKS))
  allocate(option%k2bc(MAXBCBLOCKS))
  allocate(option%j1bc(MAXBCBLOCKS))
  allocate(option%j2bc(MAXBCBLOCKS))
  allocate(option%i1bc(MAXBCBLOCKS))
  allocate(option%i2bc(MAXBCBLOCKS))
    
  allocate(option%k1src(MAXSRC))
  allocate(option%k2src(MAXSRC))
  allocate(option%j1src(MAXSRC))
  allocate(option%j2src(MAXSRC))
  allocate(option%i1src(MAXSRC))
  allocate(option%i2src(MAXSRC))
!GEH - Structured Grid Dependence - End
  allocate(option%timesrc(MAXSRCTIMES,MAXSRC))
  allocate(option%tempsrc(MAXSRCTIMES,MAXSRC))
  allocate(option%qsrc(MAXSRCTIMES,MAXSRC))
  allocate(option%csrc(MAXSRCTIMES,MAXSRC))
  allocate(option%hsrc(MAXSRCTIMES,MAXSRC))
  option%qsrc =0.D0; option%csrc =0.D0; option%hsrc =0.D0

!GEH - Structured Grid Dependence - Begin          
  allocate(option%i1reg(MAXPERMREGIONS))
  allocate(option%i2reg(MAXPERMREGIONS))
  allocate(option%j1reg(MAXPERMREGIONS))
  allocate(option%j2reg(MAXPERMREGIONS))
  allocate(option%k1reg(MAXPERMREGIONS))
  allocate(option%k2reg(MAXPERMREGIONS))
!GEH - Structured Grid Dependence - End
  allocate(option%icap_reg(MAXPERMREGIONS))
  allocate(option%ithrm_reg(MAXPERMREGIONS))
  allocate(option%por_reg(MAXPERMREGIONS))
  allocate(option%tor_reg(MAXPERMREGIONS))
  allocate(option%perm_reg(MAXPERMREGIONS,4))
  
  allocate(option%i1ini(MAXINITREGIONS))
  allocate(option%i2ini(MAXINITREGIONS))
  allocate(option%j1ini(MAXINITREGIONS))
  allocate(option%j2ini(MAXINITREGIONS))
  allocate(option%k1ini(MAXINITREGIONS))
  allocate(option%k2ini(MAXINITREGIONS))

  if (option%use_mph == PETSC_TRUE .or. option%use_owg == PETSC_TRUE .or. &
      option%use_vadose == PETSC_TRUE .or. option%use_flash == PETSC_TRUE .or. &
      option%use_richards == PETSC_TRUE) then
    allocate(option%xx_ini(option%ndof,MAXINITREGIONS))
    allocate(option%iphas_ini(MAXINITREGIONS))
  else
    allocate(option%pres_ini(MAXINITREGIONS))
    allocate(option%temp_ini(MAXINITREGIONS))
    allocate(option%sat_ini(MAXINITREGIONS))
    allocate(option%xmol_ini(MAXINITREGIONS))
    allocate(option%conc_ini(MAXINITREGIONS))
  endif

!GEH - Structured Grid Dependence - Begin    
  allocate(option%i1brk(MAXINITREGIONS))
  allocate(option%i2brk(MAXINITREGIONS))
  allocate(option%j1brk(MAXINITREGIONS))
  allocate(option%j2brk(MAXINITREGIONS))
  allocate(option%k1brk(MAXINITREGIONS))
  allocate(option%k2brk(MAXINITREGIONS))
!GEH - Structured Grid Dependence - End
  allocate(option%ibrktyp(MAXINITREGIONS))
  allocate(option%ibrkface(MAXINITREGIONS))
  
  if (idcdm == 1) then
!GEH - Structured Grid Dependence - Begin
    allocate(option%i1dcm(MAXINITREGIONS))
    allocate(option%i2dcm(MAXINITREGIONS))
    allocate(option%j1dcm(MAXINITREGIONS))
    allocate(option%j2dcm(MAXINITREGIONS))
    allocate(option%k1dcm(MAXINITREGIONS))
    allocate(option%k2dcm(MAXINITREGIONS))
!GEH - Structured Grid Dependence - End
    allocate(option%fracture_aperture(MAXINITREGIONS))
    allocate(option%matrix_block(MAXINITREGIONS))
  endif
      
  allocate(option%rock_density(MAXPERMREGIONS))
  allocate(option%cpr(MAXPERMREGIONS))
  allocate(option%dencpr(MAXPERMREGIONS))
  allocate(option%ckdry(MAXPERMREGIONS))
  allocate(option%ckwet(MAXPERMREGIONS))
  allocate(option%tau(MAXPERMREGIONS))
  allocate(option%cdiff(MAXPERMREGIONS))
  allocate(option%cexp(MAXPERMREGIONS))

  allocate(option%icaptype(MAXPERMREGIONS))
  if (option%use_mph == PETSC_TRUE .or. option%use_owg == PETSC_TRUE .or. &
      option%use_vadose == PETSC_TRUE .or. option%use_flash == PETSC_TRUE .or. &
      option%use_richards == PETSC_TRUE) then
    allocate(option%sir(1:option%nphase,MAXPERMREGIONS))
  else
    allocate(option%swir(MAXPERMREGIONS))
  endif
  allocate(option%lambda(MAXPERMREGIONS))
  allocate(option%alpha(MAXPERMREGIONS))
  allocate(option%pckrm(MAXPERMREGIONS))
  allocate(option%pcwmax(MAXPERMREGIONS))
  allocate(option%pcbetac(MAXPERMREGIONS))
  allocate(option%pwrprm(MAXPERMREGIONS))
  
  !-----------------------------------------------------------------------
  ! Set up boundary condition storage on blocks
  !-----------------------------------------------------------------------
  allocate(option%velocitybc0(option%nphase, MAXBCREGIONS))
  if (option%use_mph==PETSC_TRUE .or. option%use_owg == PETSC_TRUE &
      .or. option%use_vadose == PETSC_TRUE .or. option%use_flash == PETSC_TRUE &
      .or. option%use_richards == PETSC_TRUE) then
    allocate(option%xxbc0(option%ndof,MAXBCREGIONS))
    allocate(option%iphasebc0(MAXBCREGIONS))
!    allocate(option%dmax(0:option%ndof-1))
    option%xxbc0=0.D0
    option%iphasebc0=0
  else
    allocate(option%pressurebc0(option%nphase, MAXBCREGIONS))
    allocate(option%tempbc0(MAXBCREGIONS))
    allocate(option%sgbc0(MAXBCREGIONS))
    allocate(option%concbc0(MAXBCREGIONS))
    
      
    ! initialize
    option%pressurebc0 = 0.d0
    option%tempbc0 = 0.d0
    option%concbc0 = 0.d0
    option%sgbc0 = 0.d0
     
  endif
  option%velocitybc0 = 0.d0
  !set scale factor for heat equation, i.e. use units of MJ for energy
  option%scale = 1.d-6


  allocate(option%rtot(grid%nlmax,2))
  !allocate(grid%rrtot(grid%nlmax,grid%ncomp))

  option%rtot=0.D0
!  grid%qu_rate=0.D0


  call readInput(simulation,filename)

! check number of dofs and phases
  if (option%use_cond == PETSC_TRUE) then
    if (option%ndof .ne. 1 .or. option%nphase .ne. 1) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: COND ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else if (option%use_th == PETSC_TRUE) then
    if (option%ndof .ne. 2 .or. option%nphase .ne. 1) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: TH ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else if (option%use_thc == PETSC_TRUE) then
    if (option%ndof .ne. 3 .or. option%nphase .ne. 1) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: THC ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else if (option%use_2ph == PETSC_TRUE) then
    if (option%ndof .ne. 4 .or. option%nphase .ne. 2) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: 2PH ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else if (option%use_mph == PETSC_TRUE) then
    if (option%ndof .ne. (option%nspec+1)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: MPH ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else if (option%use_richards == PETSC_TRUE) then
    if (option%ndof .ne. (option%nspec+1)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: Richards ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else if (option%use_flash == PETSC_TRUE) then
    if (option%ndof .ne. (option%nspec+1)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: FLA ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else if (option%use_vadose == PETSC_TRUE) then
    if (option%ndof .ne. (option%nspec+1)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: VAD ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else if (option%use_owg == PETSC_TRUE) then
    if (option%ndof .ne. (option%nspec)) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: OWG ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  else
    if (option%ndof .ne. 1 .or. option%nphase .ne. 1) then
      write(*,*) 'Specified number of dofs or phases not correct-stop: &
        &LIQUID ', &
        'ndof= ',option%ndof,' nph= ',option%nphase
      stop
    endif
  endif
  
  call computeGridCoordinates(grid,option)


  if (option%myrank == 0) then
    write(*,'(/,"++++++++++++++++++++++++++++++++++++++++++++++++++++&
      &++++++++")')
    if (grid%is_structured) then
      write(*,'(" number of processors = ",i5,", npx,y,z= ",3i5)') &
        option%commsize,grid%structured_grid%npx,grid%structured_grid%npy, &
        grid%structured_grid%npz
    endif
    write(*,'(" number of dofs = ",i3,", number of phases = ",i3,i2)') &
      option%ndof,option%nphase
    if (option%use_cond == PETSC_TRUE) then
      write(*,'(" mode = Conduction: T")')
    else if (option%use_th == PETSC_TRUE) then
      write(*,'(" mode = TH: p, T")')
    else if (option%use_thc == PETSC_TRUE) then
      write(*,'(" mode = THC: p, T, C")')
    else if (option%use_2ph == PETSC_TRUE) then
      write(*,'(" mode = 2-PH: p, T, s, C")')
    else if (option%use_mph == PETSC_TRUE) then
      write(*,'(" mode = MPH: p, T, s/C")')
    else if (option%use_flash == PETSC_TRUE) then
      write(*,'(" mode = flash: p, T, z")')
    else if (option%use_vadose == PETSC_TRUE) then
      write(*,'(" mode = VAD: p, T, s/C")')
    else if (option%use_richards == PETSC_TRUE) then
      write(*,'(" mode = Richards: p, T, s/C")')
    else if (option%use_owg == PETSC_TRUE) then
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
  if (option%use_analytical == PETSC_TRUE) then
  
    option%ideriv = 1
  
    call createJacobian(grid,option)
  
!   if (myrank == 0) write(*,'(" analytical jacobian as ")'); &
!                    print *, grid%iblkfmt

#if 0    
    if (option%use_cond == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, CondJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    else if (option%use_th == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, THJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    else if (option%use_thc == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, THCJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    else if (option%use_2ph == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, TTPHASEJacobian, &
                         grid, ierr); CHKERRQ(ierr)
    else if (option%use_mph == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, MPHASEJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
    else if (option%use_flash == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, FlashJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
    else if (option%use_vadose == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, VADOSEJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
#endif
!    else if (option%use_richards == PETSC_TRUE) then
    if (option%use_richards == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, RichardsJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(solution,solver)
#if 0      
    else if (option%use_owg == PETSC_TRUE) then
      call SNESSetJacobian(option%snes, option%J, option%J, OWGJacobian, &
                         grid, ierr); CHKERRQ(ierr)
      if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
    else
      call SNESSetJacobian(option%snes, option%J, option%J, LiquidJacobian, &
                         grid, ierr); CHKERRQ(ierr)
#endif                         
    endif

  else if (option%use_matrix_free == PETSC_TRUE) then
  
    option%ideriv = 0
  
    if (option%myrank == 0) write(*,'(" Using matrix-free Newton-Krylov")')
    
    if (option%use_cond == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%ttemp, option%J, ierr)
    else if (option%use_th == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%xx, option%J, ierr)
    else if (option%use_thc == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%xx, option%J, ierr)
    else if (option%use_2ph == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%xx, option%J, ierr)
    else if (option%use_mph == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%xx, option%J, ierr)
    else if (option%use_flash == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%xx, option%J, ierr)
    else if (option%use_richards == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%xx, option%J, ierr)
    else if (option%use_vadose == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%xx, option%J, ierr)
    else if (option%use_owg == PETSC_TRUE) then
      call MatCreateMFFD(option%snes, option%xx, option%J, ierr)
    else
      call MatCreateMFFD(option%snes, option%ppressure, option%J, ierr)
    endif
    
    ! It seems that I ought to call SNESSetJacobian here now, but I don't know
    ! what function I am supposed to pass to it to get it to use one of the 
    ! finite-difference routines for computing the Jacobian.  It might not
    ! actually matter if -snes_mf has been specified.
    ! Pernice thinks that perhaps the I need to provide a function which 
    ! simply calls MatAssemblyBegin/End.
    call SNESSetJacobian(option%snes, option%J, option%J, &
                         ComputeMFJacobian, PETSC_NULL_OBJECT, ierr)

    ! Use "Walker-Pernice" differencing.
    call MatMFFDSetType(option%J, MATMFFD_WP, ierr)

    if (option%print_hhistory == PETSC_TRUE) then
      allocate(option%hhistory(HHISTORY_LENGTH))
      call MatMFFDSetHHistory(option%J, option%hhistory, HHISTORY_LENGTH, ierr)
    endif

    if (option%monitor_h == PETSC_TRUE) then
      call SNESMonitorSet(option%snes, MonitorH, grid, &
                          PETSC_NULL_OBJECT, ierr)
    endif
    
  else
  
    option%ideriv = 0
  
    if (option%myrank == 0) write(*,'(" numerical jacobian")')
    ! We will compute the Jacobian via finite differences and store it.
    
    ! Create matrix with correct parallel layout and nonzero structure to 
    ! hold the Jacobian.
      ! MatFDColoringCreate() currently does not support MATMPIBAIJ.
      
    i = option%iblkfmt
    option%iblkfmt = 0   ! to turn off MATMPIBAIJ
    call createJacobian(grid,option)
    option%iblkfmt = i

    call createColoring(grid,option,iscoloring)
        
    call MatFDColoringCreate(option%J, iscoloring, option%matfdcoloring, ierr)
    
    call ISColoringDestroy(iscoloring, ierr)

#if 0
    if (option%use_cond == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, CondResidual, &
                                        option, ierr)
    else if (option%use_th == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, THResidual, &
                                      option, ierr)
    else if (option%use_thc == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, THCResidual, &
                                      option, ierr)
    else if (option%use_2ph == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, TTPHASEResidual, &
                                      option, ierr)
    else if (option%use_mph == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, MPHASEResidual, &
                                      option, ierr)
    else if (option%use_flash == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, FlashResidual, &
                                      option, ierr)
    else if (option%use_vadose == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, VADOSEResidual, &
                                      option, ierr)
#endif                                      
!    else if (option%use_richards == PETSC_TRUE) then
    if (option%use_richards == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, RichardsResidual, &
                                      option, ierr)
#if 0                                      
    else if (option%use_owg == PETSC_TRUE) then
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, OWGResidual, &
                                      option, ierr)
    else
      call MatFDColoringSetFunctionSNES(option%matfdcoloring, LiquidResidual, &
                                      option, ierr)
#endif                                      
    endif
    
    call MatFDColoringSetFromOptions(option%matfdcoloring, ierr)
    
    call SNESSetJacobian(option%snes, option%J, option%J, &
                         SNESDefaultComputeJacobianColor,  &
                         option%matfdcoloring, ierr)
  endif

  if (option%myrank == 0) write(*,'("++++++++++++++++++++++++++++++++&
                     &++++++++++++++++++++++++++++",/)')

#if 0
  if (option%use_cond == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, CondResidual, grid, ierr)
  else if (option%use_th == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, THResidual, grid, ierr)
  else if (option%use_thc == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, THCResidual, grid, ierr)
  else if (option%use_2ph == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, TTPHASEResidual, grid, ierr)
  else if (option%use_mph == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, MPHASEResidual, grid, ierr)
  else if (option%use_flash == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, FlashResidual, grid, ierr)
  else if (option%use_vadose == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, VADOSEResidual, grid, ierr)
#endif    
  !else if (option%use_Richards == PETSC_TRUE) then
  if (option%use_Richards == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, RichardsResidual, grid, ierr)
#if 0    
  else if (option%use_owg == PETSC_TRUE) then
    call SNESSetFunction(option%snes, option%r, OWGResidual, grid, ierr)
  else
    call SNESSetFunction(option%snes, option%r, LiquidResidual, grid, ierr)
#endif
  endif

  ! Set the tolerances for the Newton solver.
  call SNESSetTolerances(option%snes, option%atol, option%rtol, option%stol, & 
                         option%maxit, option%maxf, ierr)

  call SNESSetFromOptions(option%snes, ierr)
  
! shell for custom convergence test.  The default SNES convergence test 
! is call within this function.
  call SNESSetConvergenceTest(option%snes,PFLOWConvergenceTest, &
                              option,ierr)
                              
  if (option%myrank == 0) write(*,'("  Finished setting up of SNES 1")')
  
  call SNESLineSearchGetParams(option%snes, alpha, maxstep, steptol, ierr) 
  if (option%myrank == 0) write(*,'("  Finished setting up of SNES 2")')
  call SNESLineSearchSetParams(option%snes, alpha, maxstep, option%stol, ierr) 
  if (option%myrank == 0) write(*,'("  Finished setting up of SNES 3")')

  call SNESGetKSP(option%snes, option%ksp, ierr)
  call KSPSetTolerances(option%ksp,option%rtol,option%atol,option%dtol, &
      10000,ierr)

  

 if (option%myrank == 0) write(*,'("  Finished setting up of SNES ")')


  
end subroutine initPFLOW

! ************************************************************************** !
!
! readSelectCardsFromInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! **************************************************************************
subroutine readSelectCardsFromInput(solution,filename,mcomp,mphas)

  use Option_module
  use Grid_module
  use fileio_module
  use Solution_module
  
  implicit none

  type(solution_type) :: solution
  character(len=MAXWORDLENGTH) :: filename
  integer :: mcomp, mphas

  integer :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXNAMELENGTH) :: name
  character(len=MAXCARDLENGTH) :: card
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  integer :: igeom
  
  grid => solution%grid
  option => solution%option
  
  call MPI_Comm_rank(PETSC_COMM_WORLD,option%myrank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,option%commsize,ierr)
  
  open(IUNIT1, file=filename, action="read", status="old") 
  open(IUNIT2, file='pflow.out', action="write", status="unknown")

  mphas=0
  mcomp = 0

! Read in select required cards
!.........................................................................

  ! GRID information
  string = "GRID"
  call fiFindStringInFile(IUNIT1,string,ierr)
  call fiFindStringErrorMsg(string,ierr)

  ! strip card from front of string
  call fiReadWord(string,word,.false.,ierr)
 
  ! key off igeom for structured vs unstructured 
  call fiReadInt(string,igeom,ierr)
  call fiDefaultMsg('igeom',ierr)
  
  grid => createGrid(igeom) 

  if (grid%is_structured) then ! structured
    call fiReadInt(string,grid%structured_grid%nx,ierr)
    call fiDefaultMsg('nx',ierr)
    
    call fiReadInt(string,grid%structured_grid%ny,ierr)
    call fiDefaultMsg('ny',ierr)
    
    call fiReadInt(string,grid%structured_grid%nz,ierr)
    call fiDefaultMsg('nz',ierr)
    
    grid%structured_grid%nxy = grid%structured_grid%nx*grid%structured_grid%ny
    grid%structured_grid%nmax = grid%structured_grid%nxy * &
                                grid%structured_grid%nz
    grid%nmax = grid%structured_grid%nmax
  else ! unstructured
  endif
      
  call fiReadInt(string,option%nphase,ierr)
  call fiDefaultMsg('nphase',ierr)

  call fiReadInt(string,option%nspec,ierr)
  call fiDefaultMsg('nspec',ierr)

  call fiReadInt(string,option%npricomp,ierr)
  call fiDefaultMsg('npricomp',ierr)

  call fiReadInt(string,option%ndof,ierr)
  call fiDefaultMsg('ndof',ierr)
      
  call fiReadInt(string,option%idcdm,ierr)
  call fiDefaultMsg('idcdm',ierr)

  call fiReadInt(string,option%itable,ierr)
  call fiDefaultMsg('itable',ierr)

  if (option%myrank==0) then
    if (grid%is_structured) then
      write(IUNIT2,'(/," *GRID ",/, &
        &"  igeom   = ",i4,/, &
        &"  nx      = ",i4,/, &
        &"  ny      = ",i4,/, &
        &"  nz      = ",i4,/, &
        &"  nphase  = ",i4,/, &
        &"  ndof    = ",i4,/, &
        &"  idcdm   = ",i4,/, &
        &"  itable  = ",i4    &
        &   )') grid%igeom, grid%structured_grid%nx, &
                grid%structured_grid%ny, grid%structured_grid%nz, &
                option%nphase, option%ndof, option%idcdm, option%itable
    else
    endif
  endif

!.........................................................................

  if (grid%is_structured) then  ! look for processor decomposition
    
    ! PROC information
    string = "PROC"
    call fiFindStringInFile(IUNIT1,string,ierr)

    if (ierr /= 0) then

      ! strip card from front of string
      call fiReadWord(string,word,.false.,ierr)
      call fiReadInt(string,grid%structured_grid%nx,ierr)
      call fiDefaultMsg('npx',ierr)
      call fiReadInt(string,grid%structured_grid%ny,ierr)
      call fiDefaultMsg('npy',ierr)
      call fiReadInt(string,grid%structured_grid%nz,ierr)
      call fiDefaultMsg('npz',ierr)
 
      if (option%myrank == 0) &
        write(IUNIT2,'(/," *PROC",/, &
          & "  npx   = ",3x,i4,/, &
          & "  npy   = ",3x,i4,/, &
          & "  npz   = ",3x,i4)') grid%structured_grid%nx, &
            grid%structured_grid%ny, grid%structured_grid%nz
  
        call MPI_Comm_size(PETSC_COMM_WORLD,option%commsize,ierr)
        if (option%commsize /= grid%structured_grid%nx * &
                                 grid%structured_grid%ny * &
                                 grid%structured_grid%nz) then
          if (option%myrank==0) &
            write(*,*) 'Incorrect number of processors specified: ', &
                       grid%structured_grid%nx*grid%structured_grid%ny* &
                       grid%structured_grid%nz,' commsize = ',option%commsize
        stop
      endif
    endif
  endif
  
!.........................................................................

  ! COMP information
  string = "COMP"
  call fiFindStringInFile(IUNIT1,string,ierr)

  if (ierr /= 0) then

    mcomp = 0
    do
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('COMP',ierr)
      
      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadWord(string,name,.true.,ierr)
      call fiErrorMsg('namcx','GAS',ierr)
        
      call fiWordToUpper(name) 
      select case(name(1:len_trim(name)))
        case('H2O')
            mcomp = mcomp +1
        case('TRACER')     
            mcomp = mcomp +4
        case('CO2')
            mcomp = mcomp +2   
        case('C10H22')
            mcomp = mcomp +8
        case('AIR')              
            mcomp = mcomp +16
        case('ENG')
            mcomp = mcomp +32
        case default               
           print *,' Wrong comp name::  Stop:' 
           stop
      end select 
    enddo
  endif          

!....................................................................

  ! COMP information
  string = "PHAS"
  call fiFindStringInFile(IUNIT1,string,ierr)

  if (ierr /= 0) then

    mphas = 0
    do
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('phase',ierr)
    
      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadWord(string,name,.true.,ierr)
      call fiErrorMsg('namcx','phase',ierr)
        
      call fiWordToUpper(name) 
         
      select case(name(1:len_trim(name)))
        case('ROCK')
          mphas = mphas +1
        case('H2O')     
          mphas = mphas +2
        case('CO2')
          mphas = mphas +4   
        case('GAS')
          mphas = mphas +8
        case('NAPL1')              
          mphas = mphas +16
        case('NAPL2')
          mphas = mphas +32
        case default               
          print *,' Wrong phase name::  Stop:' 
          stop
      end select 
    enddo
  endif
    
end subroutine readSelectCardsFromInput

! ************************************************************************** !
!
! readInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine readInput(simulation,filename)

  use Simulation_module
  use Option_module
  use Grid_module
  use Structured_Grid_module
  use Solver_module
  use Material_module
  use fileio_module
  use Solution_module
  use Timestepper_module

  use pckr_module 
  
  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXWORDLENGTH) :: filename

  integer :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXNAMELENGTH) :: name
  character(len=MAXCARDLENGTH) :: card
    
  real*8, parameter:: fmwnacl = 58.44277D0, fmwh2o  = 18.01534d0
  integer :: i, i1, i2, idum, ireg, isrc, j
  integer :: ibc, ibrk, ir,np  
  
! keywords: GRID, PROC, COUP, GRAV, OPTS, TOLR, DXYZ, DIFF, RADN, HYDR,  
!           SOLV, THRM, PCKR, PHIK, INIT, TIME, DTST, BCON, SOUR, BRK, RCTR

  type(solution_type), pointer :: solution
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(solver_type), pointer :: solver
  
  solution => simulation%solution
  grid => solution%grid
  option => solution%option
  solver => simulation%stepper%solver
    
  do
    call fiReadFlotranString(IUNIT1, string, ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    call fiCharsToUpper(word,len_trim(word))
    call fiReadCard(word,card,ierr)

    if (option%myrank == 0) print *, 'pflow_read:: ',card

    select case(card)

!....................

      case ('GRID')

!....................

      case ('PROC')
!.....................
      case ('COMP') 
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('COMP',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        enddo
!.....................
      case ('PHAS')
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('PHASE',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        enddo 
      
!....................

      case ('COUP')

        call fiReadStringErrorMsg('COUP',ierr)

        call fiReadInt(string,option%isync,ierr)
        call fiDefaultMsg('isync',ierr)

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *COUP",/, &
            & "  isync      = ",3x,i2 &
            & )') option%isync

!....................

      case ('GRAV')

        call fiReadStringErrorMsg('GRAV',ierr)

        call fiReadDouble(string,option%gravity,ierr)
        call fiDefaultMsg('gravity',ierr)

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1pe12.4 &
            & )') option%gravity

!....................

      case ('HDF5')
        option%print_hdf5 = .true.
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          call fiCharsToUpper(word,len_trim(word))
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              option%print_hdf5_velocities = .true.
            case('FLUX')
              option%print_hdf5_flux_velocities = .true.
            case default
          end select
            
        enddo

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *HDF5",10x,i1,/)') option%print_hdf5

!....................

      case ('TECP')
        option%print_tecplot = .true.
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          call fiCharsToUpper(word,len_trim(word))
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              option%print_tecplot_velocities = .true.
            case('FLUX')
              option%print_tecplot_flux_velocities = .true.
            case default
          end select
          
        enddo

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *TECP",10x,i1,/)') option%print_tecplot

!....................


      case ('OPTS')

        call fiReadStringErrorMsg('OPTS',ierr)

        call fiReadInt(string,option%write_init,ierr)
        call fiDefaultMsg('write_init',ierr)

        call fiReadInt(string,option%iprint,ierr)
        call fiDefaultMsg('iprint',ierr)

        call fiReadInt(string,option%imod,ierr)
        call fiDefaultMsg('mod',ierr)

        call fiReadInt(string,option%itecplot,ierr)
        call fiDefaultMsg('itecplot',ierr)

        call fiReadInt(string,option%iblkfmt,ierr)
        call fiDefaultMsg('iblkfmt',ierr)

        call fiReadInt(string,option%ndtcmx,ierr)
        call fiDefaultMsg('ndtcmx',ierr)

        call fiReadInt(string,option%iran_por,ierr)
        call fiDefaultMsg('iran_por',ierr)
  
        call fiReadDouble(string,option%ran_fac,ierr)
        call fiDefaultMsg('ran_fac',ierr)
    
        call fiReadInt(string,option%iread_perm,ierr)
        call fiDefaultMsg('iread_perm',ierr)
    
        call fiReadInt(string,option%iread_geom,ierr)
        call fiDefaultMsg('iread_geom',ierr)


        if (option%myrank == 0) &
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
            & )') option%write_init,option%iprint,option%imod,option%itecplot, &
            option%iblkfmt,option%ndtcmx,option%iran_por,option%ran_fac, &
            option%iread_perm,option%iread_geom

!....................

      case ('TOLR')

        call fiReadStringErrorMsg('TOLR',ierr)

        call fiReadInt(string,option%stepmax,ierr)
        call fiDefaultMsg('stepmax',ierr)
  
        call fiReadInt(string,option%iaccel,ierr)
        call fiDefaultMsg('iaccel',ierr)

        call fiReadInt(string,option%newton_max,ierr)
        call fiDefaultMsg('newton_max',ierr)

        call fiReadInt(string,option%icut_max,ierr)
        call fiDefaultMsg('icut_max',ierr)

        call fiReadDouble(string,option%dpmxe,ierr)
        call fiDefaultMsg('dpmxe',ierr)

        call fiReadDouble(string,option%dtmpmxe,ierr)
        call fiDefaultMsg('dtmpmxe',ierr)
  
        call fiReadDouble(string,option%dcmxe,ierr)
        call fiDefaultMsg('dcmxe',ierr)

        call fiReadDouble(string,option%dsmxe,ierr)
        call fiDefaultMsg('dsmxe',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *TOLR ",/, &
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
          option%stepmax,option%iaccel,option%newton_max,option%icut_max, &
          option%dpmxe,option%dtmpmxe,option%dcmxe, option%dsmxe

!....................

      case ('DXYZ')
      
        if (grid%is_structured) then  ! look for processor decomposition
          call readStructuredDXYZ(grid%structured_grid,option)
        else
          if (option%myrank == 0) &
            print *, 'ERROR: Keyword "DXYZ" not supported for unstructured grid'
            stop
        endif

!....................


      case('RAD0')
    
        if (grid%is_structured) then  ! look for processor decomposition
          call fiReadDouble(string,grid%structured_grid%Radius_0,ierr)
          call fiDefaultMsg('R_0',ierr)
        else
          if (option%myrank == 0) &
            print *, 'ERROR: Keyword "RAD0" not supported for unstructured grid'
            stop
        endif


      case ('DIFF')

        call fiReadStringErrorMsg('DIFF',ierr)

        call fiReadDouble(string,option%difaq,ierr)
        call fiDefaultMsg('difaq',ierr)

        call fiReadDouble(string,option%delhaq,ierr)
        call fiDefaultMsg('delhaq',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *DIFF ",/, &
          &"  difaq       = ",1pe12.4,"[m^2/s]",/, &
          &"  delhaq      = ",1pe12.4,"[kJ/mol]")') &
          option%difaq,option%delhaq

!....................

      case ('RCTR')

        call fiReadStringErrorMsg('RCTR',ierr)

        call fiReadInt(string,option%ityprxn,ierr)
        call fiDefaultMsg('ityprxn',ierr)

        call fiReadDouble(string,option%rk,ierr)
        call fiDefaultMsg('rk',ierr)

        call fiReadDouble(string,option%phis0,ierr)
        call fiDefaultMsg('phis0',ierr)

        call fiReadDouble(string,option%areas0,ierr)
        call fiDefaultMsg('areas0',ierr)

        call fiReadDouble(string,option%pwrsrf,ierr)
        call fiDefaultMsg('pwrsrf',ierr)

        call fiReadDouble(string,option%vbars,ierr)
        call fiDefaultMsg('vbars',ierr)

        call fiReadDouble(string,option%ceq,ierr)
        call fiDefaultMsg('ceq',ierr)

        call fiReadDouble(string,option%delHs,ierr)
        call fiDefaultMsg('delHs',ierr)

        call fiReadDouble(string,option%delEs,ierr)
        call fiDefaultMsg('delEs',ierr)

        call fiReadDouble(string,option%wfmts,ierr)
        call fiDefaultMsg('wfmts',ierr)

        if (option%myrank == 0) &
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
          & )') option%ityprxn,option%rk,option%phis0,option%areas0,option%pwrsrf, &
          option%vbars,option%ceq,option%delHs,option%delEs,option%wfmts

 ! convert: mol/cm^2 -> mol/cm^3 -> mol/dm^3 (note area 1/cm)          
        option%rk = option%rk * option%areas0 * 1.d3
        option%vbars = option%vbars * 1.d-3 ! convert: cm^3/mol -> L/mol
      
        option%delHs = option%delHs * option%wfmts * 1.d-3 ! convert kJ/kg -> kJ/mol
!        option%delHs = option%delHs * option%scale ! convert J/kmol -> MJ/kmol

!....................

      case ('RADN')

        call fiReadStringErrorMsg('RADN',ierr)

        call fiReadDouble(string,option%ret,ierr)
        call fiDefaultMsg('ret',ierr)

        call fiReadDouble(string,option%fc,ierr)
        call fiDefaultMsg('fc',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *RADN ",/, &
          &"  ret     = ",1pe12.4,/, &
          &"  fc      = ",1pe12.4)') &
          option%ret,option%fc

!....................


      case ('PHAR')

        call fiReadStringErrorMsg('PHAR',ierr)

        call fiReadDouble(string,option%qu_kin,ierr)
        call fiDefaultMsg('TransReaction',ierr)
        if (option%myrank==0) write(IUNIT2,'(/," *PHAR ",1pe12.4)')option%qu_kin
        option%yh2o_in_co2 = 0.d0
        if (option%qu_kin > 0.d0) option%yh2o_in_co2 = 1.d-2 ! check this number!
     
!......................

      case('RICH')
        call fiReadStringErrorMsg('RICH',ierr)
        call fiReadDouble(string,option%pref,ierr)
        call fiDefaultMsg('Ref. Pressure',ierr) 

!......................

      case('BRIN')
        call fiReadStringErrorMsg('BRIN',ierr)
        call fiReadDouble(string,option%m_nacl,ierr)
        call fiDefaultMsg('NaCl Concentration',ierr) 

        call fiReadWord(string,word,.false.,ierr)
        call fiWordToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            option%m_nacl = option%m_nacl /fmwnacl/(1.D0-option%m_nacl)
          case('MOLE')    
            option%m_nacl = option%m_nacl /fmwh2o/(1.D0-option%m_nacl)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select 
         if (option%myrank == 0) print *, option%m_nacl
!......................

      case ('HYDR')

        call fiReadStringErrorMsg('HYDR',ierr)
  
        call fiReadInt(string,option%ihydrostatic,ierr)
        call fiDefaultMsg('ihydrostatic',ierr)
       
        call fiReadDouble(string,option%dTdz,ierr)
        call fiDefaultMsg('dTdz',ierr)

        call fiReadDouble(string,option%beta,ierr)
        call fiDefaultMsg('beta',ierr)

        call fiReadDouble(string,option%tref,ierr)
        call fiDefaultMsg('tref',ierr)

        call fiReadDouble(string,option%pref,ierr)
        call fiDefaultMsg('pref',ierr)

        call fiReadDouble(string,option%conc0,ierr)
        call fiDefaultMsg('conc0',ierr)

        if (option%ihydrostatic < 1) option%ihydrostatic =1
      
        if (option%myrank==0) write(IUNIT2,'(/," *HYDR ",/, &
          &"  ihydro      = ",i3,/, &
          &"  dT/dz       = ",1pe12.4,/, &
          &"  beta        = ",1pe12.4,/, &
          &"  tref        = ",1pe12.4,/, &
          &"  pref        = ",1pe12.4,/, &
          &"  conc        = ",1pe12.4 &
          &)') &
          option%ihydrostatic,option%dTdz,option%beta,option%tref,option%pref, &
          option%conc0

!....................

      case ('SOLV')
    
        call fiReadStringErrorMsg('SOLV',ierr)

!       call fiReadDouble(string,eps,ierr)
!       call fiDefaultMsg('eps',ierr)

        call fiReadDouble(string,solver%atol,ierr)
        call fiDefaultMsg('atol_petsc',ierr)

        call fiReadDouble(string,solver%rtol,ierr)
        call fiDefaultMsg('rtol_petsc',ierr)

        call fiReadDouble(string,solver%stol,ierr)
        call fiDefaultMsg('stol_petsc',ierr)
      
        solver%dtol=1.D5
!       if (option%use_ksp == 1) then
        call fiReadDouble(string,solver%dtol,ierr)
        call fiDefaultMsg('dtol_petsc',ierr)
!       endif
   
        call fiReadInt(string,solver%maxit,ierr)
        call fiDefaultMsg('maxit',ierr)
      
        call fiReadInt(string,solver%maxf,ierr)
        call fiDefaultMsg('maxf',ierr)
       
        call fiReadInt(string,solver%idt_switch,ierr)
        call fiDefaultMsg('idt',ierr)
        
        solver%inf_tol = solver%atol
        call fiReadDouble(string,solver%inf_tol,ierr)
        call fiDefaultMsg('inf_tol_pflow',ierr)
 
        if (option%myrank==0) write(IUNIT2,'(/," *SOLV ",/, &
          &"  atol_petsc   = ",1pe12.4,/, &
          &"  rtol_petsc   = ",1pe12.4,/, &
          &"  stol_petsc   = ",1pe12.4,/, &
          &"  dtol_petsc   = ",1pe12.4,/, &
          &"  maxit        = ",8x,i5,/, &
          &"  maxf        = ",8x,i5,/, &
          &"  idt         = ",8x,i5 &
          &    )') &
           solver%atol,solver%rtol,solver%stol,solver%dtol,solver%maxit, &
           solver%maxf,solver%idt_switch

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

          call fiReadDouble(string,option%rock_density(ireg),ierr)
          call fiDefaultMsg('rock_density',ierr)

          call fiReadDouble(string,option%cpr(ireg),ierr)
          call fiDefaultMsg('cpr',ierr)
        
          call fiReadDouble(string,option%ckdry(ireg),ierr)
          call fiDefaultMsg('ckdry',ierr)
        
          call fiReadDouble(string,option%ckwet(ireg),ierr)
          call fiDefaultMsg('ckwet',ierr)
        
          call fiReadDouble(string,option%tau(ireg),ierr)
          call fiDefaultMsg('tau',ierr)

          call fiReadDouble(string,option%cdiff(ireg),ierr)
          call fiDefaultMsg('cdiff',ierr)

          call fiReadDouble(string,option%cexp(ireg),ierr)
          call fiDefaultMsg('cexp',ierr)

        !scale thermal properties
          option%cpr(ireg) = option%scale * option%cpr(ireg)
          option%dencpr(ireg) = option%rock_density(ireg) * option%cpr(ireg)
          option%ckdry(ireg) = option%scale * option%ckdry(ireg)
          option%ckwet(ireg) = option%scale * option%ckwet(ireg)
        enddo
      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *THRM: ",i3)') ireg
          write(IUNIT2,'("  itm rock_density  cpr        ckdry", &
            &                 "     ckwet       tau       cdiff     cexp")')
          write(IUNIT2,'("        [kg/m^3]  [J/kg/K]   [J/m/K/s]", &
            &              "     [J/m/K/s]     [-]        [m^2/s]       [-]")')
          do i = 1, ireg
            write(IUNIT2,'(i4,1p7e11.4)') i,option%rock_density(i), &
            option%cpr(i),option%ckdry(i),option%ckwet(i), &
            option%tau(i),option%cdiff(i),option%cexp(i)
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
          
          call fiReadInt(string,option%icaptype(idum),ierr)
          call fiDefaultMsg('icaptype',ierr)
      
          if (option%use_mph == PETSC_TRUE .or. option%use_owg == PETSC_TRUE &
              .or. option%use_vadose == PETSC_TRUE .or. option%use_flash == PETSC_TRUE&
              .or. option%use_richards == PETSC_TRUE) then
            do np=1, option%nphase
              call fiReadDouble(string,option%sir(np,idum),ierr)
              call fiDefaultMsg('sir',ierr)
            enddo 
          else
            call fiReadDouble(string,option%swir(idum),ierr)
            call fiDefaultMsg('swir',ierr)
          endif
        
          call fiReadDouble(string,option%pckrm(idum),ierr)
          call fiDefaultMsg('lambda',ierr)
          option%lambda(idum) = option%pckrm(idum)/(-option%pckrm(idum) +1.D0)
! Here the lambda is assigned as the same value of m

          call fiReadDouble(string,option%alpha(idum),ierr)
          call fiDefaultMsg('alpha',ierr)

          call fiReadDouble(string,option%pcwmax(idum),ierr)
          call fiDefaultMsg('pcwmax',ierr)
      
          call fiReadDouble(string,option%pcbetac(idum),ierr)
          call fiDefaultMsg('pcbetac',ierr)
      
          call fiReadDouble(string,option%pwrprm(idum),ierr)
          call fiDefaultMsg('pwrprm',ierr)

        enddo

          if (option%use_mph == PETSC_TRUE .or. &
              option%use_vadose == PETSC_TRUE .or. &
              option%use_flash == PETSC_TRUE .or. &
              option%use_richards == PETSC_TRUE) then
              call pckr_init(option%nphase,ireg,grid%nlmax, &
                             option%icaptype,option%sir, option%pckrm, &
                             option%lambda,option%alpha,option%pcwmax, &
                             option%pcbetac,option%pwrprm)
          endif 

      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *PCKR: ",i3)') ireg
          write(IUNIT2,'("  icp swir    lambda         alpha")')
          do j = 1, ireg
            i=option%icaptype(j)
            if (option%use_mph==PETSC_TRUE .or. &
                option%use_owg==PETSC_TRUE .or. &
                option%use_vadose == PETSC_TRUE .or. &
                option%use_flash == PETSC_TRUE .or. &
                option%use_richards == PETSC_TRUE) then
              write(IUNIT2,'(i4,1p8e12.4)') i,(option%sir(np,i),np=1, &
                option%nphase),option%lambda(i),option%alpha(i), &
                option%pcwmax(i),option%pcbetac(i),option%pwrprm(i)
            else
              write(IUNIT2,'(i4,1p7e12.4)') i,option%swir(i), &
                option%lambda(i),option%alpha(i),option%pcwmax(i), &
                option%pcbetac(i),option%pwrprm(i)
            endif
          enddo
        end if

        if (option%use_mph == PETSC_TRUE .or. &
            option%use_vadose == PETSC_TRUE .or. &
            option%use_flash == PETSC_TRUE .or. &
            option%use_richards == PETSC_TRUE) then
          deallocate(option%icaptype, option%pckrm, option%lambda, &
                     option%alpha,option%pcwmax, option%pcbetac, &
                     option%pwrprm)
        endif 
 

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
      
!GEH - Structured Grid Dependence - Begin      
          call fiReadInt(string,option%i1reg(ireg),ierr) 
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,option%i2reg(ireg),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,option%j1reg(ireg),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,option%j2reg(ireg),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,option%k1reg(ireg),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,option%k2reg(ireg),ierr)
          call fiDefaultMsg('k2',ierr)
!GEH - Structured Grid Dependence - End
  
          call fiReadInt(string,option%icap_reg(ireg),ierr)
          call fiDefaultMsg('icap',ierr)
  
          call fiReadInt(string,option%ithrm_reg(ireg),ierr)
          call fiDefaultMsg('ithrm',ierr)
  
          call fiReadDouble(string,option%por_reg(ireg),ierr)
          call fiDefaultMsg('por',ierr)
  
          call fiReadDouble(string,option%tor_reg(ireg),ierr)
          call fiDefaultMsg('tor',ierr)
  
          call fiReadDouble(string,option%perm_reg(ireg,1),ierr)
          call fiDefaultMsg('permx',ierr)
  
          call fiReadDouble(string,option%perm_reg(ireg,2),ierr)
          call fiDefaultMsg('permy',ierr)
  
          call fiReadDouble(string,option%perm_reg(ireg,3),ierr)
          call fiDefaultMsg('permz',ierr)
  
          call fiReadDouble(string,option%perm_reg(ireg,4),ierr)
          call fiDefaultMsg('permpwr',ierr)

!          call fiReadDouble(string,option%Perm_reg(ireg,5),ierr)
!          call fiDefaultMsg('porokin',ierr)

    
#if 0
          call readRegion(string)
          word = adjustl(string(1:min(31,len_trim(string))))  ! remove leading spaces
          call fiCharsToLower(string,4)
          if (fiStringCompare(string,'file',4)) then ! unstructured
            id = readUnstructuredRegion()
          else
            id = readStructuredRegion(string)
          endif
#endif          
    
        enddo
        option%iregperm = ireg
        
        if (ireg > MAXINITREGIONS) then
          print *,'error: increase MAXINITREGIONS: regions ',ireg,' > ',MAXINITREGIONS
          stop
        endif

        if (option%myrank==0) then
          write(IUNIT2,'(/," *PHIK: ireg = ",i4)') option%iregperm
          write(IUNIT2,'("  i1  i2  j1  j2  k1  k2 icap ithrm  por      tor  &
            &",   "     permx      permy      permz [m^2]   permpwr")')
          do ireg = 1, option%iregperm
!GEH - Structured Grid Dependence - Begin
            write(IUNIT2,'(6i4,2i4,1p6e11.4)') &
                  option%i1reg(ireg),option%i2reg(ireg), &
                  option%j1reg(ireg),option%j2reg(ireg), &
                  option%k1reg(ireg),option%k2reg(ireg), &
                  option%icap_reg(ireg),option%ithrm_reg(ireg), &
                  option%por_reg(ireg),option%tor_reg(ireg), &
                  (option%perm_reg(ireg,i),i=1,4)
!GEH - Structured Grid Dependence - End
          enddo
        endif

!....................
      
      case ('INIT')
    
        call fiReadInt(string,option%iread_init,ierr) 
        call fiDefaultMsg('iread_init',ierr)
      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *INIT: iread = ",i2)') option%iread_init
        endif
      
        if (option%iread_init == 0 .or. option%iread_init == 2) then
      
          ireg = 0
          do
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('INIT',ierr)

            if (string(1:1) == '.' .or. string(1:1) == '/') exit
            ireg = ireg + 1
            
!GEH - Structured Grid Dependence - Begin
            call fiReadInt(string,option%i1ini(ireg),ierr) 
            call fiDefaultMsg('i1',ierr)
            call fiReadInt(string,option%i2ini(ireg),ierr)
            call fiDefaultMsg('i2',ierr)
            call fiReadInt(string,option%j1ini(ireg),ierr)
            call fiDefaultMsg('j1',ierr)
            call fiReadInt(string,option%j2ini(ireg),ierr)
            call fiDefaultMsg('j2',ierr)
            call fiReadInt(string,option%k1ini(ireg),ierr)
            call fiDefaultMsg('k1',ierr)
            call fiReadInt(string,option%k2ini(ireg),ierr)
            call fiDefaultMsg('k2',ierr)
!GEH - Structured Grid Dependence - End

            if (option%use_mph==PETSC_TRUE .or. &
                option%use_owg==PETSC_TRUE .or. &
                option%use_vadose == PETSC_TRUE .or. &
                option%use_flash == PETSC_TRUE .or. &
                option%use_richards == PETSC_TRUE) then
              call fiReadInt(string,option%iphas_ini(ireg),ierr)
              call fiDefaultMsg('iphase',ierr)
         
              do j=1,option%ndof
                call fiReadDouble(string,option%xx_ini(j,ireg),ierr)
                call fiDefaultMsg('xxini',ierr)
              enddo
            else
              call fiReadDouble(string,option%pres_ini(ireg),ierr)
              call fiDefaultMsg('pres',ierr)
  
              call fiReadDouble(string,option%temp_ini(ireg),ierr)
              call fiDefaultMsg('temp',ierr)
  
              call fiReadDouble(string,option%sat_ini(ireg),ierr)
              call fiDefaultMsg('sat',ierr)
!              option%sat_ini(ireg)=1.D0 - option%sat_ini(ireg)
  
              call fiReadDouble(string,option%conc_ini(ireg),ierr)
              call fiDefaultMsg('conc',ierr)
            endif
          enddo
      
          option%iregini = ireg
      
          if (option%myrank==0) then
            write(IUNIT2,'("  ireg = ",i4)') option%iregini
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]   &
              &   ",    "sl [-]      c [mol/L]")')
            do ireg = 1, option%iregini
!GEH - Structured Grid Dependence - Begin
              if (option%use_mph==PETSC_TRUE .or. option%use_owg==PETSC_TRUE &
                  .or. option%use_vadose == PETSC_TRUE .or. option%use_flash == PETSC_TRUE&
                  .or. option%use_richards == PETSC_TRUE) then
                write(IUNIT2,'(7i4,1p10e12.4)') &
                  option%i1ini(ireg),option%i2ini(ireg), &
                  option%j1ini(ireg),option%j2ini(ireg), &
                  option%k1ini(ireg),option%k2ini(ireg), &
                  option%iphas_ini(ireg),(option%xx_ini(np,ireg),np =1,option%ndof)
              else
                write(IUNIT2,'(6i4,1p10e12.4)') &
                  option%i1ini(ireg),option%i2ini(ireg), &
                  option%j1ini(ireg),option%j2ini(ireg), &
                  option%k1ini(ireg),option%k2ini(ireg), &
                  option%pres_ini(ireg),option%temp_ini(ireg),option%sat_ini(ireg), &
                  option%conc_ini(ireg)
              endif
!GEH - Structured Grid Dependence - End
            enddo
          endif

        else if (option%iread_init == 1) then
    
!     read in initial conditions from file: pflow_init.dat
          if (option%myrank == 0) then
            write(*,*) '--> read in initial conditions from file: &
                        &pflow_init.dat'
  
            open(IUNIT3, file='pflow_init.dat', action="read", status="old")

            ireg = 0
            do
              call fiReadFlotranString(IUNIT3,string,ierr)
!             call fiReadStringErrorMsg('INIT',ierr)

              if (string(1:1) == '.' .or. string(1:1) == '/') exit
              ireg = ireg + 1

!GEH - Structured Grid Dependence - Begin
              call fiReadInt(string,option%i1ini(ireg),ierr) 
              call fiDefaultMsg('i1',ierr)
              call fiReadInt(string,option%i2ini(ireg),ierr)
              call fiDefaultMsg('i2',ierr)
              call fiReadInt(string,option%j1ini(ireg),ierr)
              call fiDefaultMsg('j1',ierr)
              call fiReadInt(string,option%j2ini(ireg),ierr)
              call fiDefaultMsg('j2',ierr)
              call fiReadInt(string,option%k1ini(ireg),ierr)
              call fiDefaultMsg('k1',ierr)
              call fiReadInt(string,option%k2ini(ireg),ierr)
              call fiDefaultMsg('k2',ierr)
!GEH - Structured Grid Dependence - End

              if (option%use_mph==PETSC_TRUE .or. option%use_owg==PETSC_TRUE &
                  .or. option%use_vadose == PETSC_TRUE .or. option%use_flash == PETSC_TRUE&
                  .or. option%use_richards == PETSC_TRUE) then
                call fiReadInt(string,option%iphas_ini(ireg),ierr)
                call fiDefaultMsg('iphase_ini',ierr)
            
                do j=1,option%ndof
                  call fiReadDouble(string,option%xx_ini(j,ireg),ierr)
                  call fiDefaultMsg('xx_ini',ierr)
                enddo
              else
  
                call fiReadDouble(string,option%pres_ini(ireg),ierr)
                call fiDefaultMsg('pres',ierr)
  
                call fiReadDouble(string,option%temp_ini(ireg),ierr)
                call fiDefaultMsg('temp',ierr)
  
                call fiReadDouble(string,option%sat_ini(ireg),ierr)
                call fiDefaultMsg('sat',ierr)

                call fiReadDouble(string,option%conc_ini(ireg),ierr)
                call fiDefaultMsg('conc',ierr)
              endif
       
            enddo
            option%iregini = ireg
            close(IUNIT3)
          endif
        endif

!....................

      case ('TIME')

        call fiReadStringErrorMsg('TIME',ierr)
      
        call fiReadWord(string,word,.false.,ierr)
      
        option%tunit = trim(word)

        if (option%tunit == 's') then
          option%tconv = 1.d0
        else if (option%tunit == 'm') then
          option%tconv = 60.d0
        else if (option%tunit == 'h') then
          option%tconv = 60.d0 * 60.d0
        else if (option%tunit == 'd') then
          option%tconv = 60.d0 * 60.d0 * 24.d0
        else if (option%tunit == 'mo') then
          option%tconv = 60.d0 * 60.d0 * 24.d0 * 30.d0
        else if (option%tunit == 'y') then
          option%tconv = 60.d0 * 60.d0 * 24.d0 * 365.d0
        else
          if (option%myrank == 0) then
            write(*,'(" Time unit: ",a3,/, &
              &" Error: time units must be one of ",/, &
              &"   s -seconds",/,"   m -minutes",/,"   h -hours",/, &
              &"   d -days", /, "  mo -months",/,"   y -years")') option%tunit
          endif
          stop
        endif

        call fiReadInt(string,option%kplot,ierr) 
        call fiDefaultMsg('kplot',ierr)
      
        allocate(option%tplot(option%kplot))
      
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('TIME',ierr)
        i2 = 0
        do
          i1 = i2 + 1
          i2 = i2+10
          if (i2 > option%kplot) i2 = option%kplot
          do i = i1, i2
            call fiReadDouble(string,option%tplot(i),ierr)
            call fiDefaultMsg('tplot',ierr)
          enddo
          if (i2 == option%kplot) exit
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('TIME',ierr)
        enddo

!       call fiReadFlotranString(IUNIT1,string,ierr)
!       call fiReadStringErrorMsg('TIME',ierr)
      
!       call fiReadDouble(string,option%dt,ierr)
!       call fiDefaultMsg('dt',ierr)

!       call fiReadDouble(string,option%dt_max,ierr)
!       call fiDefaultMsg('dt_max',ierr)

        if (option%myrank==0) then
          write(IUNIT2,'(/," *TIME ",a3,1x,i4,/,(1p10e12.4))') option%tunit, &
          option%kplot,(option%tplot(i),i=1,option%kplot)
!         write(IUNIT2,'("  dt= ",1pe12.4,", dtmax= ",1pe12.4,/)') &
!         option%dt,option%dt_max
        endif
      
        ! convert time units to seconds
        do i = 1, option%kplot
          option%tplot(i) = option%tconv * option%tplot(i)
        enddo
!       option%dt = option%tconv * option%dt
!       option%dt_max = option%tconv * option%dt_max

!....................

      case ('DTST')

        call fiReadStringErrorMsg('DTST',ierr)
  
        call fiReadInt(string,option%nstpmax,ierr)
        call fiDefaultMsg('nstpmax',ierr)
  
        allocate(option%tstep(option%nstpmax))
        allocate(option%dtstep(option%nstpmax))
  
        do i = 1, option%nstpmax
          call fiReadDouble(string,option%tstep(i),ierr)
          call fiDefaultMsg('tstep',ierr)
        enddo

        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('DTST',ierr)
        call fiReadDouble(string,option%dt_min,ierr)
        call fiDefaultMsg('dt_min',ierr)
        do i = 1, option%nstpmax
          call fiReadDouble(string,option%dtstep(i),ierr)
          call fiDefaultMsg('dtstep',ierr)
        enddo
        
        option%dt_max = option%dtstep(1)
        
        option%dt = option%dt_min
      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *DTST ",i4,/," tstep= ",(1p10e12.4))')  &
            option%nstpmax, (option%tstep(i),i=1,option%nstpmax)
          write(IUNIT2,'(" dtstep= ",1p10e12.4,/)') &
            option%dt_min,(option%dtstep(i),i=1,option%nstpmax)
        endif
      
        ! convert time units to seconds
        do i = 1, option%nstpmax
          option%tstep(i) = option%tconv * option%tstep(i)
          option%dtstep(i) = option%tconv * option%dtstep(i)
        enddo
        option%dt = option%tconv * option%dt
        option%dt_min = option%tconv * option%dt_min
        option%dt_max = option%tconv * option%dt_max

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
        option%iregbc1(1) = 1
        do ! loop over blocks
        
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('BCON',ierr)
        
          if (string(1:1) == '.' .or. string(1:1) == '/') exit

          ibc = ibc + 1  ! RTM: Number of boundary conditions
        
          call fiReadInt(string,option%ibndtyp(ibc),ierr)
          call fiDefaultMsg('ibndtyp',ierr)

          call fiReadInt(string,option%iface(ibc),ierr)
          call fiDefaultMsg('iface',ierr)

          do ! loop over regions
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('BCON',ierr)
        
            if (string(1:1) == '.' .or. string(1:1) == '/') exit
            ir = ir + 1

!GEH - Structured Grid Dependence - Begin
            call fiReadInt(string,option%i1bc(ir),ierr)
            call fiDefaultMsg('i1',ierr)
            call fiReadInt(string,option%i2bc(ir),ierr)
            call fiDefaultMsg('i2',ierr)
            call fiReadInt(string,option%j1bc(ir),ierr)
            call fiDefaultMsg('j1',ierr)
            call fiReadInt(string,option%j2bc(ir),ierr)
            call fiDefaultMsg('j2',ierr)
            call fiReadInt(string,option%k1bc(ir),ierr)
            call fiDefaultMsg('k1',ierr)
            call fiReadInt(string,option%k2bc(ir),ierr)
            call fiDefaultMsg('k2',ierr)    
!GEH - Structured Grid Dependence - End

            ! Now read the velocities or pressures, depending on the BC type
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('BCON',ierr)
   
            if (option%use_mph /= PETSC_TRUE .and.  &
                option%use_owg /= PETSC_TRUE .and. &
                option%use_vadose /= PETSC_TRUE .and. &
                option%use_flash /= PETSC_TRUE .and. &
                option%use_richards /= PETSC_TRUE) then  
              j=1
              if (option%nphase>1) j=2
              if (option%ibndtyp(ibc) == 1) then 
                call fiReadDouble(string, option%pressurebc0(j,ibc), ierr)
                call fiDefaultMsg("Error reading pressure BCs:", ierr)
              else if (option%ibndtyp(ibc) == 2) then
                call fiReadDouble(string, option%velocitybc0(j,ibc), ierr)
                call fiDefaultMsg("Error reading velocity BCs:", ierr)
              else if (option%ibndtyp(ibc) == 4) then
                call fiReadDouble(string, option%velocitybc0(j,ibc), ierr)
                call fiDefaultMsg("Error reading velocity BCs:", ierr)  
              else
                call fiReadDouble(string, option%pressurebc0(j,ibc), ierr)
                call fiDefaultMsg("Error reading pressure BCs:", ierr)
              endif

              if (option%nphase>1) option%pressurebc0(1,ibc) = &
                  option%pressurebc0(2,ibc)
                ! For simple input
              call fiReadDouble(string,option%tempbc0(ibc),ierr)
              call fiDefaultMsg('tempbc',ierr)

              call fiReadDouble(string,option%sgbc0(ibc),ierr)
              call fiDefaultMsg('sgbc',ierr)
              option%sgbc0(ibc) = 1.D0 - option%sgbc0(ibc) ! read in sl

              call fiReadDouble(string,option%concbc0(ibc),ierr)
              call fiDefaultMsg('concbc',ierr)
      
            else
     
              call fiReadInt(string,option%iphasebc0(ir),ierr)
              call fiDefaultMsg('iphase',ierr)
       
              if (option%ibndtyp(ir) == 1 .or. option%ibndtyp(ir) == 3 &
                  .or. option%ibndtyp(ir) == 4) then 
                do j=1,option%ndof
                  call fiReadDouble(string,option%xxbc0(j,ir),ierr)
                  call fiDefaultMsg('xxbc',ierr)
                enddo
              elseif (option%ibndtyp(ir) == 2) then
                do j=1, option%nphase       
                  call fiReadDouble(string, option%velocitybc0(j,ir), ierr)
                  call fiDefaultMsg("Error reading velocity BCs:", ierr)
                enddo
                do j=2,option%ndof
                  call fiReadDouble(string,option%xxbc0(j,ir),ierr)
                  call fiDefaultMsg('xxbc',ierr)
                enddo
              endif
            endif               
          enddo ! End loop over regions.
        
          option%iregbc2(ibc) = ir
          if (ibc+1 > MAXBCBLOCKS) then
            write(*,*) 'Too many boundary condition blocks specified--stop: ', &
              ibc+1, MAXBCBLOCKS
            stop
          else
            option%iregbc1(ibc+1) = option%iregbc2(ibc)+1
          endif
        enddo ! End loop over blocks.
      
        option%nblkbc = ibc

!GEH - Structured Grid Dependence - Begin
        if (option%myrank == 0) then
          write(IUNIT2,'(/," *BCON: nblkbc = ",i4)') option%nblkbc
          do ibc = 1, option%nblkbc
            write(IUNIT2,'("  ibndtyp = ",i3," iface = ",i2)') &
              option%ibndtyp(ibc), option%iface(ibc)
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]  &
              &  c",     " [mol/L]")')
            do ireg = option%iregbc1(ibc), option%iregbc2(ibc)
              if (option%use_mph== PETSC_TRUE .or. option%use_owg== PETSC_TRUE &
                  .or. option%use_vadose == PETSC_TRUE .or. option%use_flash == PETSC_TRUE&
                  .or. option%use_richards == PETSC_TRUE) then
                if (option%ibndtyp(ibc) == 1 .or. option%ibndtyp(ibc) == 3) then
                  write(IUNIT2,'(7i4,1p10e12.4)') &
                    option%i1bc(ireg),option%i2bc(ireg), &
                    option%j1bc(ireg),option%j2bc(ireg), &
                    option%k1bc(ireg),option%k2bc(ireg), &
                    option%iphasebc0(ireg),(option%xxbc0(j,ireg),j=1,option%ndof)
                else if (option%ibndtyp(ibc) == 2 .or. option%ibndtyp(ibc) == 4) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    option%i1bc(ireg),option%i2bc(ireg), &
                    option%j1bc(ireg),option%j2bc(ireg), &
                    option%k1bc(ireg),option%k2bc(ireg), &
                    (option%velocitybc0(j,ireg),j=1,option%nphase),&
                    (option%xxbc0(j,ireg),j=2,option%ndof)
                endif
              else
                if (option%ibndtyp(ibc) == 1 .or. option%ibndtyp(ibc) == 3) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    option%i1bc(ireg),option%i2bc(ireg), &
                    option%j1bc(ireg),option%j2bc(ireg), &
                    option%k1bc(ireg),option%k2bc(ireg), &
                    (option%pressurebc0(j,ireg),j=1,option%nphase), &
                    option%tempbc0(ireg), &
                    option%concbc0(ireg)
                else if (option%ibndtyp(ibc) == 2 .or. option%ibndtyp(ibc) == 4) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    option%i1bc(ireg),option%i2bc(ireg), &
                    option%j1bc(ireg),option%j2bc(ireg), &
                    option%k1bc(ireg),option%k2bc(ireg), &
                    (option%velocitybc0(j,ireg),j=1,option%nphase), &
                    option%tempbc0(ireg), &
                    option%concbc0(ireg)
                endif
              endif
            enddo
          enddo
        endif
!GEH - Structured Grid Dependence - End

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

!GEH - Structured Grid Dependence - Begin
          call fiReadInt(string,option%i1src(ir),ierr)
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,option%i2src(ir),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,option%j1src(ir),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,option%j2src(ir),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,option%k1src(ir),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,option%k2src(ir),ierr)
          call fiDefaultMsg('k2',ierr)    
!GEH - Structured Grid Dependence - End

!         print *,'pflowgrid_mod: Source', isrc, ir   
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
        
            call fiReadDouble(string,option%timesrc(i,isrc),ierr)
            call fiDefaultMsg('timesrc',ierr)

            call fiReadDouble(string,option%tempsrc(i,isrc),ierr)
            call fiDefaultMsg('tempsrc',ierr)
      
            call fiReadDouble(string,option%qsrc(i,isrc),ierr)
            call fiDefaultMsg('qsrc',ierr)
      
            call fiReadDouble(string,option%csrc(i,isrc),ierr)
            call fiDefaultMsg('csrc',ierr)
         
            call fiReadDouble(string,option%hsrc(i,isrc),ierr)
            call fiDefaultMsg('hsrc',ierr)


          enddo ! End loop over time.

          option%ntimsrc = i

          if (option%ntimsrc > MAXSRCTIMES) then
            write(*,*) 'Too many source times specified--stop: ', &
              option%ntimsrc, MAXSRCTIMES
            stop
          endif

          if (isrc+1 > MAXSRC) then
            write(*,*) 'Too many source blocks specified--stop: ', &
            isrc+1, MAXSRC
            stop
          endif
        enddo ! End loop over sources.

        option%nblksrc = isrc

!GEH - Structured Grid Dependence - Begin
        if (option%myrank == 0) then
          write(IUNIT2,'(/," *SOURce: nblksrc = ",i4)') option%nblksrc
          do isrc = 1, option%nblksrc
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2")')
            write(IUNIT2,'(6i4)') &
              option%i1src(isrc),option%i2src(isrc), &
              option%j1src(isrc),option%j2src(isrc), &
              option%k1src(isrc),option%k2src(isrc)
            write(IUNIT2,'("    t [s]        T [C]    QH2O [kg/s]    &
              &QCO2 [kg/s]")')
            do ir = 1, option%ntimsrc
              write(IUNIT2,'(1p10e12.4)') &
                option%timesrc(ir,isrc),option%tempsrc(ir,isrc), &
                option%qsrc(ir,isrc), &
                option%csrc(ir,isrc)
            enddo
          enddo
        endif
!GEH - Structured Grid Dependence - End

!....................
      
      case ('BRK')

        ibrk = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('BRK',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ibrk = ibrk + 1

!GEH - Structured Grid Dependence - Begin
          call fiReadInt(string,option%i1brk(ibrk),ierr) 
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,option%i2brk(ibrk),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,option%j1brk(ibrk),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,option%j2brk(ibrk),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,option%k1brk(ibrk),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,option%k2brk(ibrk),ierr)
          call fiDefaultMsg('k2',ierr)
!GEH - Structured Grid Dependence - End

          call fiReadInt(string,option%ibrktyp(ibrk),ierr)
          call fiDefaultMsg('ibrktyp',ierr)

          call fiReadInt(string,option%ibrkface(ibrk),ierr)
          call fiDefaultMsg('ibrkface',ierr)
        enddo
        option%ibrkcrv = ibrk
            
        if (option%myrank==0) then
          write(IUNIT2,'(/," *BRK: ibrk = ",i4)') option%ibrkcrv
          write(IUNIT2,'("  i1  i2  j1  j2  k1  k2  ibrktyp  ibrkface  ")')
          do ibrk = 1, option%ibrkcrv
            write(IUNIT2,'(6i4,4x,i2,7x,i2)') &
              option%i1brk(ibrk),option%i2brk(ibrk), &
              option%j1brk(ibrk),option%j2brk(ibrk), &
              option%k1brk(ibrk),option%k2brk(ibrk),option%ibrktyp(ibrk), &
              option%ibrkface(ibrk)
          enddo
        endif

        if (option%ndof == 1) option%ibrkcrv = 0

!....................
      case('SDST')
         
          
          do j=1,option%ndof
            call fiReadDouble(string,option%steady_eps(j),ierr)
            call fiDefaultMsg('steady tol',ierr)
          enddo
        if (option%myrank==0) write(IUNIT2,'(/," *SDST ",/, &
          &"  dpdt        = ",1pe12.4,/, &
          &"  dtmpdt        = ",1pe12.4,/, &
          &"  dcdt        = ",1pe12.4)') &
          option%steady_eps

!....................
      case default
    
        if (option%myrank == 0) then
          print *, "Error reading input file: keyword not found. Terminating."
        endif
        call PetscFinalize(ierr)
        stop

    end select

  enddo

  close(IUNIT1)
  
end subroutine readInput

end module Init_module
