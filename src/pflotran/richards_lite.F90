! introduced grid variables: e_total :: 1 dof
!Richards translator_module                           c_total :: option%nspec dof
!                            p_total :: 1 dof
!                            s_total :: (option%nphase-1) dof
!  stands for the accumulation term at last time step, except the /Dt part 
!  should be updated in pflowgrid_mod.F90 :: pflowgrid_step          

               
module Richards_Lite_module

  implicit none
  
  private 

#include "definitions.h"
  
#include "include/finclude/petsc.h"
!#include "include/petscf90.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
!#ifdef USE_PETSC216
!#include "include/finclude/petscsles.h"
!#endif
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"


! Cutoff parameters
  real*8, parameter :: eps       = 1.D-8
  real*8, parameter :: floweps   = 1.D-24
  real*8, parameter :: perturbation_tolerance = 1.d-5

  public RichardsLiteResidual,RichardsLiteJacobian, &
         RichardsLiteUpdateFixedAccum,RichardsLiteTimeCut,&
         RichardsLiteSetup, RichardsLiteNumericalJacTest, &
         RichardsLiteGetVarFromArray, RichardsLiteMaxChange

  public :: createRichardsLiteZeroArray
  integer, save :: n_zero_rows = 0
  logical, save :: inactive_cells_exist = .false.
  integer, pointer, save :: zero_rows_local(:)  ! 1-based indexing
  integer, pointer, save :: zero_rows_local_ghosted(:) ! 0-based indexing

  type richards_type
    real*8 :: pres
    real*8 :: temp
    real*8 :: sat
    real*8 :: den
    real*8 :: avgmw
    real*8 :: pc
!    real*8 :: vis
!    real*8 :: dvis_dp
!    real*8 :: kr
!    real*8 :: dkr_dp
    real*8 :: kvr
    real*8 :: dsat_dp
    real*8 :: dden_dp
    real*8 :: dkvr_dp
  end type richards_type
  
  type(richards_type), pointer :: aux_vars(:), aux_vars_bc(:)

contains

! ************************************************************************** !
!
! RichardsLiteTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsLiteTimeCut(realization)
 
  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  PetscScalar, pointer :: xx_p(:),yy_p(:)
  integer :: ierr
  integer :: local_id

  grid => realization%grid
  option => realization%option
  field => realization%field
 
  call VecGetArrayF90(field%xx, xx_p, ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr)

  do local_id=1, grid%nlmax
    xx_p(local_id)= yy_p(local_id)
  enddo 
  call VecRestoreArrayF90(field%xx, xx_p, ierr) 
  call VecRestoreArrayF90(field%yy, yy_p, ierr)
 
end subroutine RichardsLiteTimeCut
  
! ************************************************************************** !
!
! RichardsLiteSetup: Creates arrays for auxilliary variables
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsLiteSetup(realization)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  use Region_module
  use Structured_Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition

  integer :: ghosted_id, iconn, sum_connection
  
  grid => realization%grid
  option => realization%option
  
  allocate(aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call initAuxVar(aux_vars(ghosted_id),option)
  enddo
  
  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection%num_connections
    boundary_condition => boundary_condition%next
  enddo
  allocate(aux_vars_bc(sum_connection))
  do iconn = 1, sum_connection
    call initAuxVar(aux_vars_bc(iconn),option)
  enddo
  
end subroutine RichardsLiteSetup

! ************************************************************************** !
!
! initAuxVar: Zeros the auxilliary variable for a grid cell
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine initAuxVar(aux_var,option)

  use Option_module

  implicit none
  
  type(richards_type) :: aux_var
  type(option_type) :: option

  aux_var%pres = 0.d0
  aux_var%temp = 0.d0
  aux_var%sat = 0.d0
  aux_var%den = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%pc = 0.d0
!  aux_var%kr = 0.d0
!  aux_var%dkr_dp = 0.d0
!  aux_var%vis = 0.d0
!  aux_var%dvis_dp = 0.d0
  aux_var%kvr = 0.d0
  aux_var%dsat_dp = 0.d0
  aux_var%dden_dp = 0.d0
  aux_var%dkvr_dp = 0.d0

end subroutine initAuxVar

! ************************************************************************** !
!
! copyAuxVar: Copies an auxilliary variable
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine copyAuxVar(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(richards_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%pres = aux_var%pres
  aux_var2%temp = aux_var%temp
  aux_var2%sat = aux_var%sat
  aux_var2%den = aux_var%den
  aux_var2%avgmw = aux_var%avgmw
  aux_var2%pc = aux_var%pc
!  aux_var2%kr = aux_var%kr
!  aux_var2%dkr_dp = aux_var%dkr_dp
!  aux_var2%vis = aux_var%vis
!  aux_var2%dvis_dp = aux_var%dvis_dp
  aux_var2%kvr = aux_var%kvr
  aux_var2%dsat_dp = aux_var%dsat_dp
  aux_var2%dden_dp = aux_var%dden_dp
  aux_var2%dkvr_dp = aux_var%dkvr_dp

end subroutine copyAuxVar
  
! ************************************************************************** !
!
! RichardsLiteUpdateAuxVars: Updates the auxilliary variables associated with 
!                        the Richards problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsLiteUpdateAuxVars(realization)

  use Realization_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  
  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_type), pointer :: cur_connection_set

  integer :: ghosted_id, local_id, sum_connection, idof, iconn
  integer :: iphasebc, iphase
  PetscScalar, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  real*8 :: xxbc(realization%option%ndof)
  PetscErrorCode :: ierr
  
  option => realization%option
  grid => realization%grid
  field => realization%field
  
  call VecGetArrayF90(field%xx_loc,xx_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif
    iphase = int(iphase_loc_p(ghosted_id))
   
    call computeAuxVarLite(xx_loc_p(ghosted_id:ghosted_id),aux_vars(ghosted_id), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       option)
    iphase_loc_p(ghosted_id) = iphase
  enddo

  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      select case(boundary_condition%condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC)
          xxbc(1) = boundary_condition%aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          xxbc(idof) = xx_loc_p(ghosted_id)
      end select
      
      select case(boundary_condition%condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC)
          iphasebc = boundary_condition%aux_int_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

      call computeAuxVarLite(xxbc(1),aux_vars_bc(sum_connection), &
                         iphasebc, &
                         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                         option)
    enddo
    boundary_condition => boundary_condition%next
  enddo


  call VecRestoreArrayF90(field%xx_loc,xx_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)

end subroutine RichardsLiteUpdateAuxVars

! ************************************************************************** !
!
! RichardsLiteUpdateFixedAccum: Updates the fixed portion of the 
!                                  accumulation term
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsLiteUpdateFixedAccum(realization)

  use Realization_module
  use Option_module
  use Field_module
  use Grid_module

  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field

  integer :: ghosted_id, local_id, iphase
  PetscScalar, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscScalar, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          ithrm_loc_p(:), accum_p(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  grid => realization%grid
  field => realization%field
  
  call VecCopy(field%xx, field%yy, ierr)   
  
  call VecGetArrayF90(field%xx,xx_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayF90(field%tor_loc,tor_loc_p,ierr)
  call VecGetArrayF90(grid%volume,volume_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)

  call VecGetArrayF90(field%accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif
!    iend = local_id*option%ndof
!    istart = iend-option%ndof+1
    iphase = int(iphase_loc_p(ghosted_id))
    call computeAuxVarLite(xx_p(local_id:local_id),aux_vars(ghosted_id), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       option)
    iphase_loc_p(ghosted_id) = iphase
    call RichardsLiteAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                                  volume_p(local_id), &
                                  option%dencpr(int(ithrm_loc_p(ghosted_id))), &
                                  option,accum_p(local_id:local_id)) 
  enddo

  call VecRestoreArrayF90(field%xx,xx_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecRestoreArrayF90(field%tor_loc,tor_loc_p,ierr)
  call VecRestoreArrayF90(grid%volume,volume_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)

  call VecRestoreArrayF90(field%accum, accum_p, ierr)

#if 0
!  call RichardsLiteNumericalJacTest(field%xx,realization)
#endif

end subroutine RichardsLiteUpdateFixedAccum

! ************************************************************************** !
!
! RichardsLiteNumericalJacTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsLiteNumericalJacTest(xx,realization)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  type(realization_type) :: realization

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  real*8 :: derivative, perturbation
  
  PetscScalar, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  integer :: idof, idof2, icell

  grid => realization%grid
  option => realization%option
  field => realization%field
  
  call VecDuplicate(xx,xx_pert,ierr)
  call VecDuplicate(xx,res,ierr)
  call VecDuplicate(xx,res_pert,ierr)
  
  call MatCreate(PETSC_COMM_WORLD,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%ndof,grid%nlmax*option%ndof,ierr)
  call MatSetType(A,MATAIJ,ierr)
  call MatSetFromOptions(A,ierr)
    
  call RichardsLiteResidual(PETSC_NULL_OBJECT,xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr)
  do icell = 1,grid%nlmax
    if (associated(field%imat)) then
      if (field%imat(grid%nL2G(icell)) <= 0) cycle
    endif
     idof = icell
!    do idof = (icell-1)*option%ndof+1,icell*option%ndof 
      call veccopy(xx,xx_pert,ierr)
      call vecgetarrayf90(xx_pert,vec_p,ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call vecrestorearrayf90(xx_pert,vec_p,ierr)
      call richardsliteresidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
      call vecgetarrayf90(res_pert,vec_p,ierr)
      do idof2 = 1, grid%nlmax*option%ndof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values,ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr)
!    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'numerical_jacobian.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call MatDestroy(A,ierr)
  
  call VecDestroy(xx_pert,ierr)
  call VecDestroy(res,ierr)
  call VecDestroy(res_pert,ierr)
  
end subroutine RichardsLiteNumericalJacTest

! ************************************************************************** !
!
! RichardsLiteAccumDerivative: Computes derivatives of the accumulation 
!                                 term for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsLiteAccumDerivative(aux_var,por,vol,rock_dencpr,option, &
                                          sat_func,J)

  use Option_module
  use Material_module
  
  implicit none

  type(richards_type) :: aux_var
  type(option_type) :: option
  real*8 vol,por,rock_dencpr
  type(saturation_function_type) :: sat_func
  real*8 :: J(option%ndof,option%ndof)
     
  integer :: ispec !, iireac=1
  real*8 :: porXvol, mol(option%nspec), eng

  integer :: iphase, ideriv
  type(richards_type) :: aux_var_pert
  real*8 :: x(1), x_pert(1), pert, res(1), res_pert(1), J_pert(1,1)

  porXvol = por*vol
      
  J(1,1) = (aux_var%sat*aux_var%dden_dp+aux_var%dsat_dp*aux_var%den)*porXvol

  if (option%numerical_derivatives) then
    call copyAuxVar(aux_var,aux_var_pert,option)
    x(1) = aux_var%pres
    call RichardsLiteAccumulation(aux_var,por,vol,rock_dencpr,option,res)
    ideriv = 1
    pert = x(ideriv)*perturbation_tolerance
    x_pert = x
    x_pert(ideriv) = x_pert(ideriv) + pert
    call computeAuxVarLite(x_pert(1),aux_var_pert,iphase,sat_func,option)
#if 0      
      select case(ideriv)
        case(1)
!         print *, 'dvis_dp:', aux_var%dvis_dp, (aux_var_pert%vis-aux_var%vis)/pert(ideriv)
!         print *, 'dkr_dp:', aux_var%dkr_dp, (aux_var_pert%kr-aux_var%kr)/pert(ideriv)
          print *, 'dsat_dp:', aux_var%dsat_dp, (aux_var_pert%sat-aux_var%sat)/pert
          print *, 'dden_dp:', aux_var%dden_dp, (aux_var_pert%den-aux_var%den)/pert
          print *, 'dkvr_dp:', aux_var%dkvr_dp, (aux_var_pert%kvr-aux_var%kvr)/pert
      end select     
#endif     
    call RichardsLiteAccumulation(aux_var_pert,por,vol,rock_dencpr,option,res_pert)
    J_pert(1,1) = (res_pert(1)-res(1))/pert
    J = J_pert
  endif
   
end subroutine RichardsLiteAccumDerivative

! ************************************************************************** !
!
! RichardsLiteAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine RichardsLiteAccumulation(aux_var,por,vol,rock_dencpr,option,Res)

  use Option_module
  
  implicit none

  type(richards_type) :: aux_var
  type(option_type) :: option
  real*8 Res(1:option%ndof) 
  real*8 vol,por,rock_dencpr
       
  Res(1) = aux_var%sat * aux_var%den * por * vol

end subroutine RichardsLiteAccumulation

! ************************************************************************** !
!
! RichardsLiteFluxDerivative: Computes the derivatives of the internal flux terms
!                         for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine RichardsLiteFluxDerivative(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,sat_func_up,sat_func_dn,Jup,Jdn)
  use Option_module 
  use Material_module                             
  
  implicit none
  
  type(richards_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  real*8 :: sir_up, sir_dn
  real*8 :: por_up, por_dn
  real*8 :: tor_up, tor_dn
  real*8 :: dd_up, dd_dn
  real*8 :: perm_up, perm_dn
  real*8 :: Dk_up, Dk_dn
  real*8 :: v_darcy,area
  real*8 :: dist_gravity  ! distance along gravity vector
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  real*8 :: Jup(option%ndof,option%ndof), Jdn(option%ndof,option%ndof)
     
  real*8 :: fluxm(option%nspec),fluxe,q
  real*8 :: ukvr,DK,Dq,diffdp
  real*8 :: upweight,density_ave,cond,gravity,dphi
  
  real*8 :: dden_ave_dp_up, dden_ave_dp_dn
  real*8 :: dgravity_dden_up, dgravity_dden_dn
  real*8 :: dphi_dp_up, dphi_dp_dn
  real*8 :: dukvr_dp_up, dukvr_dp_dn
  real*8 :: dq_dp_up, dq_dp_dn

  integer :: iphase, ideriv
  type(richards_type) :: aux_var_pert_up, aux_var_pert_dn
  real*8 :: x_up(1), x_dn(1), x_pert_up(1), x_pert_dn(1), pert_up, pert_dn, &
            res(1), res_pert_up(1), res_pert_dn(1), J_pert_up(1,1), J_pert_dn(1,1)
  
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0 
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  
! Flow term
  if (aux_var_up%sat > sir_up .or. aux_var_dn%sat > sir_dn) then
    if (aux_var_up%sat <eps) then 
      upweight=0.d0
    else if (aux_var_dn%sat <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*aux_var_up%den+(1.D0-upweight)*aux_var_dn%den
    dden_ave_dp_up = upweight*aux_var_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*aux_var_dn%dden_dp

    gravity = (upweight*aux_var_up%den*aux_var_up%avgmw + &
              (1.D0-upweight)*aux_var_dn%den*aux_var_dn%avgmw) &
              * dist_gravity
    dgravity_dden_up = upweight*aux_var_up%avgmw*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*aux_var_dn%avgmw*dist_gravity

    dphi = aux_var_up%pres - aux_var_dn%pres  + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*aux_var_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*aux_var_dn%dden_dp

! note uxmol only contains one phase xmol
    if (dphi>=0.D0) then
      ukvr = aux_var_up%kvr
      dukvr_dp_up = aux_var_up%dkvr_dp
!      dukvr_dp_dn = 0.d0
    else
      ukvr = aux_var_dn%kvr
!      dukvr_dp_up = 0.d0
      dukvr_dp_dn = aux_var_dn%dkvr_dp
    endif      
   
    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
      
      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
    endif
  endif 

  Jup = Jup*option%dt
  Jdn = Jdn*option%dt
 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%numerical_derivatives) then
    call copyAuxVar(aux_var_up,aux_var_pert_up,option)
    call copyAuxVar(aux_var_dn,aux_var_pert_dn,option)
    x_up(1) = aux_var_up%pres
    x_dn(1) = aux_var_dn%pres
    call RichardsLiteFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                      aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                      area,dist_gravity,upweight, &
                      option,v_darcy,res)
    ideriv = 1
    pert_up = x_up(ideriv)*perturbation_tolerance
    pert_dn = x_dn(ideriv)*perturbation_tolerance
    x_pert_up = x_up
    x_pert_dn = x_dn
    x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    call computeAuxVarLite(x_pert_up(1),aux_var_pert_up,iphase,sat_func_up,option)
    call computeAuxVarLite(x_pert_dn(1),aux_var_pert_dn,iphase,sat_func_dn,option)
    call RichardsLiteFlux(aux_var_pert_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                      aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                      area,dist_gravity,upweight, &
                      option,v_darcy,res_pert_up)
    call RichardsLiteFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                      aux_var_pert_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                      area,dist_gravity,upweight, &
                      option,v_darcy,res_pert_dn)
    J_pert_up(1,ideriv) = (res_pert_up(1)-res(1))/pert_up
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jup = J_pert_up
    Jdn = J_pert_dn
  endif

end subroutine RichardsLiteFluxDerivative

! ************************************************************************** !
!
! RichardsLiteFlux: Computes the internal flux terms for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine RichardsLiteFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,v_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(richards_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  real*8 :: sir_up, sir_dn
  real*8 :: por_up, por_dn
  real*8 :: tor_up, tor_dn
  real*8 :: dd_up, dd_dn
  real*8 :: perm_up, perm_dn
  real*8 :: Dk_up, Dk_dn
  real*8 :: v_darcy,area
  real*8 :: Res(1:option%ndof) 
  real*8 :: dist_gravity  ! distance along gravity vector
     
  integer :: ispec
  real*8 :: fluxm(option%nspec),fluxe,q
  real*8 :: uh,uxmol(1:option%nspec),ukvr,difff,diffdp, DK,Dq
  real*8 :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0  
  
! Flow term
  if (aux_var_up%sat > sir_up .or. aux_var_dn%sat > sir_dn) then
    if (aux_var_up%sat <eps) then 
      upweight=0.d0
    else if (aux_var_dn%sat <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*aux_var_up%den+(1.D0-upweight)*aux_var_dn%den 

    gravity = (upweight*aux_var_up%den*aux_var_up%avgmw + &
              (1.D0-upweight)*aux_var_dn%den*aux_var_dn%avgmw) &
              * dist_gravity

    dphi = aux_var_up%pres - aux_var_dn%pres  + gravity

! note uxmol only contains one phase xmol
    if (dphi>=0.D0) then
      ukvr = aux_var_up%kvr
    else
      ukvr = aux_var_dn%kvr
    endif      
   

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area

      fluxm(1) = q*density_ave       
    endif
  endif 

  Res(1) = fluxm(1) * option%dt
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine RichardsLiteFlux

! ************************************************************************** !
!
! RichardsLiteBCFluxDerivative: Computes the derivatives of the boundary flux 
!                           terms for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsLiteBCFluxDerivative(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
                                    por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                                    area,dist_gravity,option, &
                                    sat_func_dn,Jdn)
  use Option_module
  use Material_module
 
  implicit none
  
  integer :: ibndtype(:)
  type(richards_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  real*8 :: dd_up, sir_dn
  real*8 :: aux_vars(:) ! from aux_real_var array in boundary condition
  real*8 :: por_dn,perm_dn,Dk_dn,tor_dn
  real*8 :: area
  type(saturation_function_type) :: sat_func_dn  
  real*8 :: Jdn(option%ndof,option%ndof)
  
  real*8 :: dist_gravity  ! distance along gravity vector
          
  real*8 :: v_darcy
  real*8 :: fluxm(option%nspec),q,density_ave
  real*8 :: ukvr,diffdp,DK,Dq
  real*8 :: upweight,cond,gravity,dphi

  real*8 :: dden_ave_dp_dn
  real*8 :: dgravity_dden_dn
  real*8 :: dphi_dp_dn
  real*8 :: dukvr_dp_dn
  real*8 :: dq_dp_dn

  integer :: iphase, ideriv
  type(richards_type) :: aux_var_pert_dn, aux_var_pert_up
  real*8 :: perturbation
  real*8 :: x_dn(1), x_up(1), x_pert_dn(1), x_pert_up(1), pert_dn, res(1), &
            res_pert_dn(1), J_pert_dn(1,1)
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_dn = 0.d0
        
  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(RICHARDS_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC)
      Dq = perm_dn / dd_up
      ! Flow term
      if (aux_var_up%sat > sir_dn .or. aux_var_dn%sat > sir_dn) then
        upweight=1.D0
        if (aux_var_up%sat < eps) then 
          upweight=0.d0
        else if (aux_var_dn%sat < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*aux_var_up%den+(1.D0-upweight)*aux_var_dn%den
        dden_ave_dp_dn = (1.D0-upweight)*aux_var_dn%dden_dp

        gravity = (upweight*aux_var_up%den*aux_var_up%avgmw + &
                  (1.D0-upweight)*aux_var_dn%den*aux_var_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*aux_var_dn%avgmw*dist_gravity
        
        dphi = aux_var_up%pres - aux_var_dn%pres + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*aux_var_dn%dden_dp
        
        if (dphi>=0.D0) then
          ukvr = aux_var_up%kvr
        else
          ukvr = aux_var_dn%kvr
          dukvr_dp_dn = aux_var_dn%dkvr_dp
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(RICHARDS_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = aux_var_up%den
        else 
          density_ave = aux_var_dn%den
          dden_ave_dp_dn = aux_var_dn%dden_dp
        endif 
        q = v_darcy * area
      endif

  end select

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)

  Jdn = Jdn * option%dt

  if (option%numerical_derivatives) then
    call copyAuxVar(aux_var_up,aux_var_pert_up,option)
    call copyAuxVar(aux_var_dn,aux_var_pert_dn,option)
    x_up(1) = aux_var_up%pres
    x_dn(1) = aux_var_dn%pres
    ideriv = 1
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_up(ideriv) = x_dn(ideriv)
    endif
    call RichardsLiteBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
                        por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                        area,dist_gravity,option,v_darcy,res)
    if (ibndtype(RICHARDS_PRESSURE_DOF) == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    ideriv = 1
    pert_dn = x_dn(ideriv)*perturbation_tolerance    
    x_pert_dn = x_dn
    x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
    x_pert_up = x_up
    if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
      x_pert_up(ideriv) = x_pert_dn(ideriv)
    endif   
    call computeAuxVarLite(x_pert_dn(1),aux_var_pert_dn,iphase,sat_func_dn,option)
    call computeAuxVarLite(x_pert_up(1),aux_var_pert_up,iphase,sat_func_dn,option)
    call RichardsLiteBCFlux(ibndtype,aux_vars,aux_var_pert_up,aux_var_pert_dn, &
                        por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                        area,dist_gravity,option,v_darcy,res_pert_dn)
    J_pert_dn(1,ideriv) = (res_pert_dn(1)-res(1))/pert_dn
    Jdn = J_pert_dn
  endif

end subroutine RichardsLiteBCFluxDerivative

! ************************************************************************** !
!
! RichardsLiteBCFlux: Computes the  boundary flux terms for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsLiteBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
                          por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                          area,dist_gravity,option,v_darcy,Res)
  use Option_module
 
  implicit none
  
  integer :: ibndtype(:)
  type(richards_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  real*8 :: dd_up, sir_dn
  real*8 :: aux_vars(:) ! from aux_real_var array
  real*8 :: por_dn,perm_dn,Dk_dn,tor_dn
  real*8 :: v_darcy, area
  real*8 :: Res(1:option%ndof) 
  
  real*8 :: dist_gravity  ! distance along gravity vector
          
  integer ispec
  real*8 :: fluxm(option%nspec),q,density_ave
  real*8 :: ukvr,diffdp,DK,Dq
  real*8 :: upweight,cond,gravity,dphi
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(RICHARDS_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC)
      Dq = perm_dn / dd_up
      ! Flow term
      if (aux_var_up%sat > sir_dn .or. aux_var_dn%sat > sir_dn) then
        upweight=1.D0
        if (aux_var_up%sat < eps) then 
          upweight=0.d0
        else if (aux_var_dn%sat < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*aux_var_up%den+(1.D0-upweight)*aux_var_dn%den
   
        gravity = (upweight*aux_var_up%den*aux_var_up%avgmw + &
                  (1.D0-upweight)*aux_var_dn%den*aux_var_dn%avgmw) &
                  * dist_gravity
       
        dphi = aux_var_up%pres - aux_var_dn%pres + gravity
   
        if (dphi>=0.D0) then
          ukvr = aux_var_up%kvr
        else
          ukvr = aux_var_dn%kvr
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(RICHARDS_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(RICHARDS_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = aux_var_up%den
        else 
          density_ave = aux_var_dn%den
        endif 
      endif

  end select

  q = v_darcy * area

  fluxm(1) = q*density_ave

  Res(1)=fluxm(1)* option%dt

end subroutine RichardsLiteBCFlux

! ************************************************************************** !
!
! RichardsLiteResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine RichardsLiteResidual(snes,xx,r,realization,ierr)

  use water_eos_module
  use Gas_Eos_Module

  use Connection_module
  use Realization_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization

  integer :: ierr
  integer :: i, jn
  integer :: ip1, ip2
  integer :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscScalar, pointer ::accum_p(:)

  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               vl_p(:) 
                          
               
  PetscScalar, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)

  integer :: iphase, index_var_begin, index_var_end,np
  integer :: icap_up, icap_dn, ithrm_up, ithrm_dn
  real*8 :: dd_up, dd_dn, &
            accum
  real*8 :: dd, f_up, f_dn, ff
  real*8 :: perm_up, perm_dn
  real*8 :: D_up, D_dn  ! "Diffusion" constants at upstream, downstream faces.
  real*8 :: dw_kg, dw_mol, dif(realization%option%nphase)
  real*8 :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  real*8 :: tmp, upweight
  real*8 :: rho
  real*8 :: xxbc(realization%option%ndof)
  integer :: iphasebc
  real*8 :: Res(realization%option%ndof), v_darcy
 PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  logical :: enthalpy_flag
  integer :: iconn, idof
  integer :: sum_connection
  real*8 :: distance, fraction_upwind
  real*8 :: distance_gravity
  
  grid => realization%grid
  option => realization%option
  field => realization%field

  ! Communication -----------------------------------------
  ! These 3 must be called before RichardsUpdateAuxVars()
  call GridGlobalToLocal(grid,xx,field%xx_loc,NDOF)
  call GridLocalToLocal(grid,field%iphas_loc,field%iphas_loc,ONEDOF)
  call GridLocalToLocal(grid,field%icap_loc,field%icap_loc,ONEDOF)

  call GridLocalToLocal(grid,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call GridLocalToLocal(grid,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call GridLocalToLocal(grid,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call GridLocalToLocal(grid,field%ithrm_loc,field%ithrm_loc,ONEDOF)

  call RichardsLiteUpdateAuxVars(realization)

! now assign access pointer to local variables
  call VecGetArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%accum, accum_p, ierr)
 
  call VecGetArrayF90(field%yy,yy_p,ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%vl, vl_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  !print *,' Finished scattering non deriv'

  if (option%rk > 0.d0) then
    call VecGetArrayF90(field%phis,phis_p,ierr)
  endif

  r_p = 0.d0
#if 1
  ! Accumulation terms ------------------------------------
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif
    call RichardsLiteAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                                  volume_p(local_id), &
                                  option%dencpr(int(ithrm_loc_p(ghosted_id))), &
                                  option,Res) 
    r_p(local_id) = r_p(local_id) + Res(1)
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => realization%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    qsrc1 = source_sink%condition%cur_value(RICHARDS_PRESSURE_DOF)
    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
      
    cur_connection_set => source_sink%connection
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif
      
      if (qsrc1 > 0.d0) then ! injection
        tsrc1 = 25.d0
        call wateos_noderiv(tsrc1,aux_vars(ghosted_id)%pres, &
                            dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        r_p(local_id) = r_p(local_id) - qsrc1 *option%dt
      endif  
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
#if 1
  ! Interior Flux Terms -----------------------------------
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(field%imat)) then
        if (field%imat(ghosted_id_up) <= 0 .or.  &
            field%imat(ghosted_id_dn) <= 0) cycle
      endif

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
        
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
   
      D_up = option%ckwet(ithrm_up)
      D_dn = option%ckwet(ithrm_dn)

      call RichardsLiteFlux(aux_vars(ghosted_id_up),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),option%sir(1,icap_up), &
                          dd_up,perm_up,D_up, &
                        aux_vars(ghosted_id_dn),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),option%sir(1,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                        cur_connection_set%area(iconn),distance_gravity, &
                        upweight,option,v_darcy,Res)

      field%internal_velocities(1,iconn) = v_darcy

      if (local_id_up > 0) then               ! If the upstream node is not a ghost node...
        vl_p(1+3*(local_id_up-1)) = &
                           v_darcy*abs(cur_connection_set%dist(1,iconn) )
        vl_p(2+3*(local_id_up-1)) = &
                           v_darcy*abs(cur_connection_set%dist(2,iconn))
        vl_p(3+3*(local_id_up-1)) = &
                             v_darcy*abs(cur_connection_set%dist(3,iconn))
      endif
     
      if (local_id_up>0) then
        r_p(local_id_up) = r_p(local_id_up) + Res(1)
      endif
   
      if (local_id_dn>0) then
        r_p(local_id_dn) = r_p(local_id_dn) - Res(1)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = option%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))

      do idof=1,option%ndof
        select case(boundary_condition%condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC)
            xxbc(idof) = boundary_condition%aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%ndof+idof)
        end select
      enddo
      
      select case(boundary_condition%condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC)
          iphasebc = boundary_condition%aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

      icap_dn = int(icap_loc_p(ghosted_id))  

      call RichardsLiteBCFlux(boundary_condition%condition%itype, &
                                boundary_condition%aux_real_var(:,iconn), &
                                aux_vars_bc(sum_connection), &
                                aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                tor_loc_p(ghosted_id), &
                                option%sir(1,icap_dn), &
                                cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
                                cur_connection_set%area(iconn), &
                                distance_gravity,option, &
                                v_darcy,Res)
      field%boundary_velocities(1,iconn) = v_darcy

      r_p(local_id)= r_p(local_id) - Res(1)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif  

  if (inactive_cells_exist) then
    do i=1,n_zero_rows
      r_p(zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%yy, yy_p, ierr)
  call VecRestoreArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%accum, accum_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%vl, vl_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
  if (option%rk > 0.d0) then
    call VecRestoreArrayF90(field%phis,phis_p,ierr)
  endif

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'residual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'xx.out',viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

end subroutine RichardsLiteResidual
                
! ************************************************************************** !
!
! RichardsLiteJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine RichardsLiteJacobian(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module
  use gas_eos_module

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_module
  use Coupler_module
  use Field_module
  use Debug_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(realization_type) :: realization
  MatStructure flag

  integer :: ierr
  integer :: nvar,neq,nr
  integer :: ithrm_up, ithrm_dn, i
  integer :: ip1, ip2 

  PetscScalar, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), phis_p(:),  tor_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscScalar, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
  integer :: icap,iphas,iphas_up,iphas_dn,icap_up,icap_dn
  integer :: ii, jj
  integer :: index_var_begin, index_var_end
  real*8 :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
  real*8 :: vv_darcy(realization%option%nphase),voldt,pvoldt
  real*8 :: ff,dif(1:realization%option%nphase)
  real*8 :: tsrc1,qsrc1,csrc1,hsrc1
  real*8 :: dd_up, dd_dn, dd, f_up, f_dn
  real*8 :: perm_up, perm_dn
  real*8 :: dw_dp,dw_dt,hw_dp,hw_dt,dresT_dp,dresT_dt
  real*8 :: D_up, D_dn  ! "Diffusion" constants upstream and downstream of a face.
  real*8 :: zero, norm
  real*8 :: ra(1:realization%option%ndof,1:2*realization%option%ndof)  
  real*8 :: tmp, upweight
  real*8 :: delxbc(1:realization%option%ndof)
  real*8 :: blkmat11(1:realization%option%ndof,1:realization%option%ndof), &
            blkmat12(1:realization%option%ndof,1:realization%option%ndof),&
            blkmat21(1:realization%option%ndof,1:realization%option%ndof),&
            blkmat22(1:realization%option%ndof,1:realization%option%ndof)
  real*8 :: ResInc(1:realization%grid%nlmax, 1:realization%option%ndof, 1:realization%option%ndof),res(1:realization%option%ndof)  
  real*8 :: max_dev  
  real*8 :: xxbc(realization%option%ndof)
  integer :: iphasebc
  integer :: local_id, ghosted_id
  integer :: local_id_up, local_id_dn
  integer :: ghosted_id_up, ghosted_id_dn
  integer ::  natural_id_up,natural_id_dn
  
  real*8 :: Jup(realization%option%ndof,realization%option%ndof), &
            Jdn(realization%option%ndof,realization%option%ndof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  logical :: enthalpy_flag
  integer :: iconn, idof
  integer :: sum_connection  
  real*8 :: distance, fraction_upwind
  real*8 :: distance_gravity 
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option 
  type(field_type), pointer :: field  
  
  PetscViewer :: viewer
  Vec :: debug_vec
!-----------------------------------------------------------------------
! R stand for residual
!  ra       1              2              3              4          5              6            7      8
! 1: p     dR/dpi         dR/dTi          dR/dci        dR/dsi   dR/dpim        dR/dTim
! 2: T
! 3: c
! 4  s         
!-----------------------------------------------------------------------

  grid => realization%grid
  option => realization%option
  field => realization%field

! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
  flag = SAME_NONZERO_PATTERN

#if 0
!  call RichardsLiteNumericalJacTest(xx,realization)
#endif

 ! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr)

  call VecGetArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)

#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(field%imat)) then
      if (field%imat(ghosted_id) <= 0) cycle
    endif
    icap = int(icap_loc_p(ghosted_id))
    call RichardsLiteAccumDerivative(aux_vars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option, &
                              realization%saturation_function_array(icap)%ptr,&
                              Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_accum.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => realization%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (size(source_sink%condition%cur_value) > RICHARDS_CONCENTRATION_DOF) then
      enthalpy_flag = .true.
    else
      enthalpy_flag = .false.
    endif

    qsrc1 = source_sink%condition%cur_value(RICHARDS_PRESSURE_DOF)
    tsrc1 = source_sink%condition%cur_value(RICHARDS_TEMPERATURE_DOF)
    csrc1 = source_sink%condition%cur_value(RICHARDS_CONCENTRATION_DOF)
    if (enthalpy_flag) hsrc1 = source_sink%condition%cur_value(RICHARDS_ENTHALPY_DOF)

    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_set => source_sink%connection
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif
      
!      if (enthalpy_flag) then
!        r_p(local_id*option%ndof) = r_p(local_id*option%ndof) - hsrc1 * option%dt   
!      endif         

      if (qsrc1 > 0.d0) then ! injection
        call wateos(tsrc1,aux_vars(ghosted_id)%pres,dw_kg,dw_mol,dw_dp,dw_dt, &
              enth_src_h2o,hw_dp,hw_dt,option%scale,ierr)        
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        ! base on r_p() = r_p() - qsrc1*enth_src_h2o*option%dt
!        dresT_dp = -qsrc1*hw_dp*option%dt
        ! dresT_dt = -qsrc1*hw_dt*option%dt ! since tsrc1 is prescribed, there is no derivative
!        call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,dresT_dp,ADD_VALUES,ierr)
        ! call MatSetValuesLocal(A,1,istart-1,1,istart-1,dresT_dt,ADD_VALUES,ierr)
      endif  
    
      if (csrc1 > 0.d0) then ! injection
        call printErrMsg(option,"concentration source not yet implemented in RichardsLite")
      endif
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Interior Flux Terms -----------------------------------  
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id_up) <= 0 .or. &
            field%imat(ghosted_id_dn) <= 0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
      natural_id_up = grid%nG2N(ghosted_id_up)
      natural_id_dn = grid%nG2N(ghosted_id_dn)
   
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
    
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))
    
      iphas_up = iphase_loc_p(ghosted_id_up)
      iphas_dn = iphase_loc_p(ghosted_id_dn)

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      D_up = option%ckwet(ithrm_up)
      D_dn = option%ckwet(ithrm_dn)
    
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
                              
      call RichardsLiteFluxDerivative(aux_vars(ghosted_id_up),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),option%sir(1,icap_up), &
                          dd_up,perm_up,D_up, &
                        aux_vars(ghosted_id_dn),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),option%sir(1,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                        cur_connection_set%area(iconn),distance_gravity, &
                        upweight,option,&
                        realization%saturation_function_array(icap_up)%ptr,&
                        realization%saturation_function_array(icap_dn)%ptr,&
                        Jup,Jdn)
      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -1.d0*Jup
        Jdn = -1.d0*Jdn
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_flux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => realization%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(field%imat)) then
        if (field%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = option%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         OptionDotProduct(option%gravity, &
                                          cur_connection_set%dist(1:3,iconn))

      do idof=1,option%ndof
        select case(boundary_condition%condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC)
            xxbc(idof) = boundary_condition%aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%ndof+idof)
        end select
      enddo
      
      select case(boundary_condition%condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC)
          iphasebc = boundary_condition%aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

      icap_dn = int(icap_loc_p(ghosted_id))  

      call RichardsLiteBCFluxDerivative(boundary_condition%condition%itype, &
                                boundary_condition%aux_real_var(:,iconn), &
                                aux_vars_bc(sum_connection), &
                                aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                tor_loc_p(ghosted_id), &
                                option%sir(1,icap_dn), &
                                cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
                                cur_connection_set%area(iconn), &
                                distance_gravity,option, &
                                realization%saturation_function_array(icap_dn)%ptr,&
                                Jdn)
      Jdn = -1.d0*Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call VecRestoreArrayF90(field%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tor_loc, tor_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  if (option%rk > 0.d0) then
    call VecRestoreArrayF90(field%phis,phis_p,ierr)
  endif
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL
  zero = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero,ierr) 
  do i=1, n_zero_rows
    ii = mod(zero_rows_local(i),option%ndof)
    ip1 = zero_rows_local_ghosted(i)
    if (ii == 0) then
      ip2 = ip1-1
    elseif (ii == option%ndof-1) then
      ip2 = ip1+1
    else
      ip2 = ip1
    endif
    call MatSetValuesLocal(A,1,ip1,1,ip2,1.d0,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
#else
  if (inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,f_up,ierr) 
  endif
#endif

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'jacobian.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    call MatNorm(A,NORM_1,norm,ierr)
    if (option%myrank == 0) print *, '1 norm:', norm
    call MatNorm(A,NORM_FROBENIUS,norm,ierr)
    if (option%myrank == 0) print *, '2 norm:', norm
    call MatNorm(A,NORM_INFINITY,norm,ierr)
    if (option%myrank == 0) print *, 'inf norm:', norm
!    call GridCreateVector(grid,ONEDOF,debug_vec,GLOBAL)
!    call MatGetRowMaxAbs(A,debug_vec,PETSC_NULL_INTEGER,ierr)
!    call VecMax(debug_vec,i,norm,ierr)
!    call VecDestroy(debug_vec,ierr)
!    if (option%myrank == 0) print *, 'max:', i, norm
  endif

end subroutine RichardsLiteJacobian

! ************************************************************************** !
!
! createRichardsLiteZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine createRichardsLiteZeroArray(realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module
  
  implicit none

  type(realization_type) :: realization
  integer :: ncount, idof
  integer :: local_id, ghosted_id

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  integer :: flag = 0, ierr
    
  grid => realization%grid
  option => realization%option
  field => realization%field
  
  n_zero_rows = 0

  if (associated(field%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (field%imat(ghosted_id) <= 0) then
        n_zero_rows = n_zero_rows + option%ndof
      endif
    enddo
  endif

  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))
  zero_rows_local = 0
  zero_rows_local_ghosted = 0
  ncount = 0

  if (associated(field%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (field%imat(ghosted_id) <= 0) then
        do idof = 1, option%ndof
          ncount = ncount + 1
          zero_rows_local(ncount) = (local_id-1)*option%ndof+idof
          zero_rows_local_ghosted(ncount) = (ghosted_id-1)*option%ndof+idof-1
        enddo
      endif
    enddo
  endif

  call MPI_Allreduce(n_zero_rows,flag,1,MPI_INTEGER,MPI_MAX, &
                     PETSC_COMM_WORLD,ierr)
  if (flag > 0) inactive_cells_exist = .true.

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine createRichardsLiteZeroArray

! ************************************************************************** !
!
! RichardsLiteMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 01/15/08
!
! ************************************************************************** !
subroutine RichardsLiteMaxChange(realization)

  use Realization_module
  use Option_module
  use Field_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  
  integer :: ierr
  
  option => realization%option
  field => realization%field

  call VecWAXPY(field%dxx,-1.d0,field%xx,field%yy,ierr)
  call VecStrideNorm(field%dxx,0,NORM_INFINITY,option%dpmax,ierr)

end subroutine RichardsLiteMaxChange

! ************************************************************************** !
!
! RichardsLiteGetVarFromArray: Extracts variables indexed by ivar and isubvar
!                          from Richards type
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine RichardsLiteGetVarFromArray(realization,vec,ivar,isubvar)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

#define TEMPERATURE 4
#define PRESSURE 5
#define LIQUID_SATURATION 6
#define GAS_SATURATION 7
#define LIQUID_ENERGY 8
#define GAS_ENERGY 9
#define LIQUID_MOLE_FRACTION 10
#define GAS_MOLE_FRACTION 11
#define VOLUME_FRACTION 12
#define PHASE 13  
#define MATERIAL_ID 14

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  integer :: ivar
  integer :: isubvar

  integer :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscScalar, pointer :: vec_ptr(:), vec2_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  grid => realization%grid
  field => realization%field

  call RichardsLiteUpdateAuxVars(realization)

  call VecGetArrayF90(vec,vec_ptr,ierr)

  select case(ivar)
    case(TEMPERATURE,GAS_SATURATION,LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION, &
         GAS_ENERGY)
      select case(ivar)
        case(TEMPERATURE)
          call printErrMsg(option,'TEMPERATURE not supported by RichardsLite')
        case(GAS_SATURATION)
          call printErrMsg(option,'GAS_SATURATION not supported by RichardsLite')
        case(LIQUID_MOLE_FRACTION)
          call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by RichardsLite')
        case(GAS_MOLE_FRACTION)
          call printErrMsg(option,'GAS_MOLE_FRACTION not supported by RichardsLite')
        case(LIQUID_ENERGY)
          call printErrMsg(option,'LIQUID_ENERGY not supported by RichardsLite')
        case(GAS_ENERGY)
          call printErrMsg(option,'GAS_ENERGY not supported by RichardsLite')
      end select
    case(PRESSURE,LIQUID_SATURATION)
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)    
        select case(ivar)
          case(PRESSURE)
            vec_ptr(local_id) = aux_vars(ghosted_id)%pres
          case(LIQUID_SATURATION)
            vec_ptr(local_id) = aux_vars(ghosted_id)%sat
        end select
      enddo
    case(VOLUME_FRACTION)
      ! need to set minimum to 0.
      call VecGetArrayF90(field%phis,vec2_ptr,ierr)
      vec_ptr(1:grid%nlmax) = vec2_ptr(1:grid%nlmax)
      call VecRestoreArrayF90(field%phis,vec2_ptr,ierr)
    case(PHASE)
      call VecGetArrayF90(field%iphas_loc,vec2_ptr,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec2_ptr(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%iphas_loc,vec2_ptr,ierr)
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = field%imat(grid%nL2G(local_id))
      enddo
  end select
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

end subroutine RichardsLiteGetVarFromArray

! ************************************************************************** !
!
! computeAuxVarLite: Computes secondary variables for each grid cell
! author: Glenn Hammond
! date: 12/15/07
!
! ************************************************************************** !
subroutine computeAuxVarLite(x,aux_var,iphase,saturation_function,option)

  use Option_module
  use water_eos_module
  use gas_eos_module  
  use pckr_module
  use Material_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  real*8 :: x(1)
  type(richards_type) :: aux_var
  integer :: iphase

  integer :: ierr
  real*8 :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  real*8 :: kr, ds_dp, dkr_dp
  real*8 :: dvis_dt, dvis_dp, dvis_dpsat
  real*8 :: dw_dp, dw_dt, hw_dp, hw_dt
  real*8 :: dpw_dp
  
  aux_var%sat = 0.d0
  aux_var%den = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%kvr = 0.d0
  kr = 0.d0
 
  aux_var%pres = x(1)
  aux_var%temp = 25.d0
 
  aux_var%pc = option%pref - aux_var%pres

!***************  Liquid phase properties **************************
  !geh aux_var%avgmw = option%fmwh2o  ! hardwire for comparison with old code
  aux_var%avgmw = 18.0153d0

  pw = option%pref
  ds_dp = 0.d0
  dkr_dp = 0.d0
!  if (aux_var%pc > 0.d0) then
  if (aux_var%pc > 1.d0) then
    iphase = 3
    call SaturationFunctionCompute(aux_var%pres,aux_var%sat,kr, &
                                   ds_dp,dkr_dp, &
                                   saturation_function, &
                                   option)
    dpw_dp = 0
!    call pflow_pckr_richards(ipckr,aux_var%sat,aux_var%pc,kr)
  else
    iphase = 1
    aux_var%pc = 0.d0
    aux_var%sat = 1.d0  
    kr = 1.d0    
    pw = aux_var%pres
    dpw_dp = 1.d0
  endif  

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
  call wateos(aux_var%temp,pw,dw_kg,dw_mol,dw_dp,dw_dt,hw,hw_dp,hw_dt, &
              option%scale,ierr)

! may need to compute dpsat_dt to pass to VISW
  call psat(aux_var%temp,sat_pressure,ierr)
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  call VISW(aux_var%temp,pw,sat_pressure,visl,dvis_dt,dvis_dp,ierr) 
  dvis_dpsat = -dvis_dp 
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif
 
  aux_var%den = dw_mol
  aux_var%kvr = kr/visl
  
!  aux_var%vis = visl
!  aux_var%dvis_dp = dvis_dp
!  aux_var%kr = kr
!  aux_var%dkr_dp = dkr_dp
  aux_var%dsat_dp = ds_dp

  aux_var%dden_dp = dw_dp
  
  aux_var%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp

end subroutine computeAuxVarLite

end module Richards_Lite_module
