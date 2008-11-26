module THC_module

  use THC_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"


! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-5

  public THCResidual,THCJacobian, &
         THCUpdateFixedAccumulation,THCTimeCut,&
         THCNumericalJacobianTest, &
         THCMaxChange, THCUpdateSolution, &
         THCGetTecplotHeader, THCInitializeTimestep, &
         THCSetup, THCResidualToMass, &
         THCUpdateAuxVars

  PetscInt, parameter :: jh2o = 1

contains

! ************************************************************************** !
!
! THCTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine THCTimeCut(realization)
 
  use Realization_module
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
 
  call VecCopy(field%flow_yy,field%flow_xx,ierr)
 
end subroutine THCTimeCut


! ************************************************************************** !
!
! THCSetup: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine THCSetup(realization)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call THCSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THCSetup
  
! ************************************************************************** !
!
! THCSetupPatch: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine THCSetupPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Condition_module
  use Connection_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(thc_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)

  PetscInt :: ghosted_id, iconn, sum_connection
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
    
  patch%aux%THC => THCAuxCreate()
    
  ! allocate aux_var data structures for all grid cells
  allocate(aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call THCAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%aux%THC%aux_vars => aux_vars
  patch%aux%THC%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! aux_var data structures for them
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo
  allocate(aux_vars_bc(sum_connection))
  do iconn = 1, sum_connection
    call THCAuxVarInit(aux_vars_bc(iconn),option)
  enddo
  patch%aux%THC%aux_vars_bc => aux_vars_bc
  patch%aux%THC%num_aux_bc = sum_connection

  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call THCCreateZeroArray(patch,option)

end subroutine THCSetupPatch

! ************************************************************************** !
!
! THCUpdateAuxVars: Updates the auxilliary variables associated with 
!                        the THC problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine THCUpdateAuxVars(realization)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call THCUpdateAuxVarsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THCUpdateAuxVars

! ************************************************************************** !
!
! THCUpdateAuxVarsPatch: Updates the auxilliary variables associated with 
!                        the THC problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine THCUpdateAuxVarsPatch(realization)

  use Realization_module
  use Patch_module
  use Field_module
  use Option_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  
  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(thc_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  aux_vars => patch%aux%THC%aux_vars
  aux_vars_bc => patch%aux%THC%aux_vars_bc
  
  call GridVecGetArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))
   
    call THCAuxVarCompute(xx_loc_p(istart:iend), &
                       aux_vars(ghosted_id), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
                       option)
    iphase_loc_p(ghosted_id) = iphase
  enddo

  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      do idof=1,option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo
      
      select case(boundary_condition%flow_condition%itype(THC_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
          iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

      call THCAuxVarCompute(xxbc,aux_vars_bc(sum_connection), &
                         iphasebc, &
                         realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                         porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &
                         option)
    enddo
    boundary_condition => boundary_condition%next
  enddo


  call GridVecRestoreArrayF90(grid,field%flow_xx_loc,xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  
  patch%aux%THC%aux_vars_up_to_date = PETSC_TRUE

end subroutine THCUpdateAuxVarsPatch

! ************************************************************************** !
!
! THCInitializeTimestep: Update data in module prior to time step
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine THCInitializeTimestep(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  call THCUpdateFixedAccumulation(realization)

end subroutine THCInitializeTimestep

! ************************************************************************** !
!
! THCUpdateSolution: Updates data in module after a successful time step
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
subroutine THCUpdateSolution(realization)

  use Realization_module
  
  implicit none
  
  type(realization_type) :: realization

  PetscErrorCode :: ierr
  
  call VecCopy(realization%field%flow_xx,realization%field%flow_yy,ierr)   

end subroutine THCUpdateSolution

! ************************************************************************** !
!
! THCUpdateFixedAccumulation: Updates the fixed portion of the 
!                                  accumulation term
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine THCUpdateFixedAccumulation(realization)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call THCUpdateFixedAccumPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THCUpdateFixedAccumulation

! ************************************************************************** !
!
! THCUpdateFixedAccumPatch: Updates the fixed portion of the 
!                                  accumulation term
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine THCUpdateFixedAccumPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module

  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(thc_auxvar_type), pointer :: aux_vars(:)

  PetscInt :: ghosted_id, local_id, istart, iend, iphase
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tor_loc_p(:), volume_p(:), &
                          ithrm_loc_p(:), accum_p(:), perm_xx_loc_p(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  aux_vars => patch%aux%THC%aux_vars
    
  call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecGetArrayF90(grid,field%ithrm_loc,ithrm_loc_p,ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)

  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))
    call THCAuxVarCompute(xx_p(istart:iend), &
                       aux_vars(ghosted_id), &
                       iphase, &
                       realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       porosity_loc_p(ghosted_id),perm_xx_loc_p(ghosted_id), &                       
                       option)
    iphase_loc_p(ghosted_id) = iphase
    call THCAccumulation(aux_vars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,accum_p(istart:iend)) 
  enddo

  call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call GridVecRestoreArrayF90(grid,field%ithrm_loc,ithrm_loc_p,ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)

  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)

#if 0
!  call THCNumericalJacobianTest(field%flow_xx,realization)
#endif

end subroutine THCUpdateFixedAccumPatch

! ************************************************************************** !
!
! THCNumericalJacobianTest: Computes the a test numerical jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine THCNumericalJacobianTest(xx,realization)

  use Realization_module
  use Patch_module
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
  
  PetscReal :: derivative, perturbation
  
  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  
  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  
  call VecDuplicate(xx,xx_pert,ierr)
  call VecDuplicate(xx,res,ierr)
  call VecDuplicate(xx,res_pert,ierr)
  
  call MatCreate(option%comm,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof,grid%nlmax*option%nflowdof,ierr)
  call MatSetType(A,MATAIJ,ierr)
  call MatSetFromOptions(A,ierr)
    
  call THCResidual(PETSC_NULL_OBJECT,xx,res,realization,ierr)
  call GridVecGetArrayF90(grid,res,vec2_p,ierr)
  do icell = 1,grid%nlmax
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    endif
    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call veccopy(xx,xx_pert,ierr)
      call vecgetarrayf90(xx_pert,vec_p,ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call vecrestorearrayf90(xx_pert,vec_p,ierr)
      call THCResidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
      call vecgetarrayf90(res_pert,vec_p,ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values,ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr)
    enddo
  enddo
  call GridVecRestoreArrayF90(grid,res,vec2_p,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call PetscViewerASCIIOpen(option%comm,'numerical_jacobian.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call MatDestroy(A,ierr)
  
  call VecDestroy(xx_pert,ierr)
  call VecDestroy(res,ierr)
  call VecDestroy(res_pert,ierr)
  
end subroutine THCNumericalJacobianTest

! ************************************************************************** !
!
! THCAccumulationDerivative: Computes derivatives of the accumulation 
!                                 term for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine THCAccumulationDerivative(aux_var,por,vol,rock_dencpr,option, &
                                          sat_func,J)

  use Option_module
  use Material_module
  
  implicit none

  type(thc_auxvar_type) :: aux_var
  type(option_type) :: option
  PetscReal vol,por,rock_dencpr
  type(saturation_function_type) :: sat_func
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec !, iireac=1
  PetscReal :: porXvol, mol(option%nflowspec), eng

  PetscInt :: iphase, ideriv
  type(thc_auxvar_type) :: aux_var_pert
  PetscReal :: x(3), x_pert(3), pert, res(3), res_pert(3), J_pert(3,3)

  porXvol = por*vol
      
  J(1,1) = (aux_var%sat*aux_var%dden_dp+aux_var%dsat_dp*aux_var%den)*porXvol*aux_var%xmol(1)
  J(1,2) = (aux_var%sat*aux_var%dden_dt)*porXvol*aux_var%xmol(1)
  J(1,3) = 0.d0
  J(2,1) = (aux_var%sat*aux_var%dden_dp+aux_var%dsat_dp*aux_var%den)*porXvol*aux_var%xmol(2)
  J(2,2) = (aux_var%sat*aux_var%dden_dt)*porXvol*aux_var%xmol(2)
  J(2,3) = (aux_var%sat*aux_var%den)*porXvol
  J(3,1) = (aux_var%dsat_dp*aux_var%den*aux_var%u+ &
            aux_var%sat*aux_var%dden_dp*aux_var%u+ &
            aux_var%sat*aux_var%den*aux_var%du_dp)* &
           porXvol
  J(3,2) = aux_var%sat* &
           (aux_var%dden_dt*aux_var%u+ &  ! pull %sat outside
            aux_var%den*aux_var%du_dt)* &
           porXvol + (1.d0 - por)* vol * rock_dencpr 
  J(3,3) = 0.d0 

  if (option%numerical_derivatives) then
    allocate(aux_var_pert%xmol(option%nflowspec),aux_var_pert%diff(option%nflowspec))
    call THCAuxVarCopy(aux_var,aux_var_pert,option)
    x(1) = aux_var%pres
    x(2) = aux_var%temp
    x(3) = aux_var%xmol(2)
    call THCAccumulation(aux_var,por,vol,rock_dencpr,option,res)
    do ideriv = 1,3
      pert = x(ideriv)*perturbation_tolerance
      x_pert = x
      x_pert(ideriv) = x_pert(ideriv) + pert
      call THCAuxVarCompute(x_pert,aux_var_pert,iphase,sat_func, &
                                 0.d0,0.d0,option)
#if 0      
      select case(ideriv)
        case(1)
!         print *, 'dvis_dp:', aux_var%dvis_dp, (aux_var_pert%vis-aux_var%vis)/pert(ideriv)
!         print *, 'dkr_dp:', aux_var%dkr_dp, (aux_var_pert%kr-aux_var%kr)/pert(ideriv)
          print *, 'dsat_dp:', aux_var%dsat_dp, (aux_var_pert%sat-aux_var%sat)/pert
          print *, 'dden_dp:', aux_var%dden_dp, (aux_var_pert%den-aux_var%den)/pert
          print *, 'dkvr_dp:', aux_var%dkvr_dp, (aux_var_pert%kvr-aux_var%kvr)/pert
          print *, 'dh_dp:', aux_var%dh_dp, (aux_var_pert%h-aux_var%h)/pert
          print *, 'du_dp:', aux_var%du_dp, (aux_var_pert%u-aux_var%u)/pert
        case(2)
          print *, 'dden_dt:', aux_var%dden_dt, (aux_var_pert%den-aux_var%den)/pert
          print *, 'dkvr_dt:', aux_var%dkvr_dt, (aux_var_pert%kvr-aux_var%kvr)/pert
          print *, 'dh_dt:', aux_var%dh_dt, (aux_var_pert%h-aux_var%h)/pert
          print *, 'du_dt:', aux_var%du_dt, (aux_var_pert%u-aux_var%u)/pert
      end select     
#endif     
      call THCAccumulation(aux_var_pert,por,vol,rock_dencpr,option,res_pert)
      J_pert(:,ideriv) = (res_pert(:)-res(:))/pert
    enddo
    deallocate(aux_var_pert%xmol,aux_var_pert%diff)
    J = J_pert
  endif
   
end subroutine THCAccumulationDerivative

! ************************************************************************** !
!
! THCAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine THCAccumulation(aux_var,por,vol,rock_dencpr,option,Res)

  use Option_module
  
  implicit none

  type(thc_auxvar_type) :: aux_var
  type(option_type) :: option
  PetscReal Res(1:option%nflowdof) 
  PetscReal vol,por,rock_dencpr
     
  PetscInt :: ispec !, iireac=1
  PetscReal :: porXvol, mol(option%nflowspec), eng
  
 ! if (present(ireac)) iireac=ireac

  porXvol = por*vol
      
  mol=0.d0
  do ispec=1, option%nflowspec  
    mol(ispec) = mol(ispec) + aux_var%sat * &
                              aux_var%den * &
                              aux_var%xmol(ispec)
  enddo
  mol = mol * porXvol
  eng = aux_var%sat * &
        aux_var%den * &
        aux_var%u * &
        porXvol + (1.d0 - por)* vol * rock_dencpr * aux_var%temp 
 
! Reaction terms here
!  if (option%run_coupled .and. iireac>0) then
!H2O
!    mol(1)= mol(1) - option%flow_dt * option%rtot(node_no,1)
!  endif
  Res(1:option%nflowdof-1)=mol(:)
  Res(option%nflowdof)=eng

end subroutine THCAccumulation

! ************************************************************************** !
!
! THCFluxDerivative: Computes the derivatives of the internal flux terms
!                         for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine THCFluxDerivative(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,sat_func_up,sat_func_dn,Jup,Jdn)
  use Option_module 
  use Material_module                             
  
  implicit none
  
  type(thc_auxvar_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: v_darcy,area
  PetscReal :: dist_gravity  ! distance along gravity vector
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec
  PetscReal :: fluxm(option%nflowspec),fluxe,q
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  
  PetscReal :: ddifff_dp_up, ddifff_dp_dn, ddifff_dt_up, ddifff_dt_dn
  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn, dden_ave_dt_up, dden_ave_dt_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn, dphi_dt_up, dphi_dt_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn, dukvr_dt_up, dukvr_dt_dn
  PetscReal :: duh_dp_up, duh_dp_dn, duh_dt_up, duh_dt_dn
  PetscReal :: dq_dp_up, dq_dp_dn, dq_dt_up, dq_dt_dn
  PetscReal :: duxmol_dxmol_up, duxmol_dxmol_dn

  PetscInt :: iphase, ideriv
  type(thc_auxvar_type) :: aux_var_pert_up, aux_var_pert_dn
  PetscReal :: x_up(3), x_dn(3), x_pert_up(3), x_pert_dn(3), pert_up, pert_dn, &
            res(3), res_pert_up(3), res_pert_dn(3), J_pert_up(3,3), J_pert_dn(3,3)
  
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0 
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dt_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dden_ave_dt_dn = 0.d0
  ddifff_dp_up = 0.d0
  ddifff_dp_dn = 0.d0
  ddifff_dt_up = 0.d0
  ddifff_dt_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dt_up = 0.d0
  dphi_dt_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dt_up = 0.d0
  dukvr_dt_dn = 0.d0
  duh_dp_up = 0.d0
  duh_dp_dn = 0.d0
  duh_dt_up = 0.d0
  duh_dt_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  dq_dt_up = 0.d0
  dq_dt_dn = 0.d0
  duxmol_dxmol_up = 0.d0
  duxmol_dxmol_dn = 0.d0
  
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
    dden_ave_dt_up = upweight*aux_var_up%dden_dt
    dden_ave_dt_dn = (1.D0-upweight)*aux_var_dn%dden_dt

    gravity = (upweight*aux_var_up%den*aux_var_up%avgmw + &
              (1.D0-upweight)*aux_var_dn%den*aux_var_dn%avgmw) &
              * dist_gravity
    dgravity_dden_up = upweight*aux_var_up%avgmw*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*aux_var_dn%avgmw*dist_gravity

    dphi = aux_var_up%pres - aux_var_dn%pres  + gravity
    dphi_dp_up = 1.d0 + dgravity_dden_up*aux_var_up%dden_dp
    dphi_dp_dn = -1.d0 + dgravity_dden_dn*aux_var_dn%dden_dp
    dphi_dt_up = dgravity_dden_up*aux_var_up%dden_dt
    dphi_dt_dn = dgravity_dden_dn*aux_var_dn%dden_dt

! note uxmol only contains one phase xmol
    if (dphi>=0.D0) then
      ukvr = aux_var_up%kvr
      dukvr_dp_up = aux_var_up%dkvr_dp
!      dukvr_dp_dn = 0.d0
      dukvr_dt_up = aux_var_up%dkvr_dt
!      dukvr_dt_dn = 0.d0
      
      uh = aux_var_up%h
      duh_dp_up = aux_var_up%dh_dp
!      duh_dp_dn = 0.d0
      duh_dt_up = aux_var_up%dh_dt
!      duh_dt_dn = 0.d0
      
      uxmol(1:option%nflowspec) = aux_var_up%xmol(1:option%nflowspec)
      duxmol_dxmol_up = 1.d0
!      duxmol_dxmol_dn = 0.d0
    else
      ukvr = aux_var_dn%kvr
!      dukvr_dp_up = 0.d0
      dukvr_dp_dn = aux_var_dn%dkvr_dp
!      dukvr_dt_up = 0.d0
      dukvr_dt_dn = aux_var_dn%dkvr_dt
      
      uh = aux_var_dn%h
!      duh_dp_up = 0.d0
      duh_dp_dn = aux_var_dn%dh_dp
!      duh_dt_up = 0.d0
      duh_dt_dn = aux_var_dn%dh_dt
      
      uxmol(1:option%nflowspec) = aux_var_dn%xmol(1:option%nflowspec)
!      duxmol_dxmol_up = 0.d0
      duxmol_dxmol_dn = 1.d0
    endif      
   
    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
      
      dq_dt_up = Dq*(dukvr_dt_up*dphi+ukvr*dphi_dt_up)*area
      dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
        
      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uxmol(1)
      Jup(1,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uxmol(1)
!      Jup(1,3:option%nflowdof) = 0d.0

      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(1)
      Jdn(1,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(1)
!      Jdn(1,3:option%nflowdof) = 0.d0
      do ispec=2,option%nflowspec 
        ! based on flux = q*density_ave*uxmol
        Jup(ispec,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uxmol(ispec)
        Jup(ispec,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uxmol(ispec)
        Jup(ispec,ispec+1) = q*density_ave*duxmol_dxmol_up

        Jdn(ispec,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(ispec)
        Jdn(ispec,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(ispec)
        Jdn(ispec,ispec+1) = q*density_ave*duxmol_dxmol_dn
      enddo
      ! based on flux = q*density_ave*uh
      Jup(option%nflowdof,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uh+q*density_ave*duh_dp_up
      Jup(option%nflowdof,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uh+q*density_ave*duh_dt_up
!      Jup(option%nflowdof,3:option%nflowdof) = 0d.0

      Jdn(option%nflowdof,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn
      Jdn(option%nflowdof,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn
!      Jdn(option%nflowdof,3:option%nflowdof) = 0.d0

    endif
  endif 
! Diffusion term   
! Note : average rule may not be correct  
  if ((aux_var_up%sat > eps) .and. (aux_var_dn%sat > eps)) then
!    difff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)* &
!                            (aux_var_up%den+aux_var_dn%den)
!    do ispec=1, option%nflowspec
!      fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
!                 (aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))* &
!                 (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
!    enddo 
    difff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)* &
                            (aux_var_up%den+aux_var_dn%den)
    ddifff_dp_up = diffdp * 0.25D0*(aux_var_up%dsat_dp*(aux_var_up%den+aux_var_dn%den)+ &
                                    (aux_var_up%sat+aux_var_dn%sat)*aux_var_up%dden_dp)
    ddifff_dt_up = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*aux_var_up%dden_dt
    ddifff_dp_dn = diffdp * 0.25D0*(aux_var_dn%dsat_dp*(aux_var_up%den+aux_var_dn%den)+ &
                                    (aux_var_up%sat+aux_var_dn%sat)*aux_var_dn%dden_dp)
    ddifff_dt_dn = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*aux_var_dn%dden_dt
                                    
    Jup(1,1) = Jup(1,1)+ddifff_dp_up*0.5d0*(aux_var_up%diff(1) + aux_var_dn%diff(1))*&
                                           (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
    Jup(1,2) = Jup(1,2)+ddifff_dt_up*0.5d0*(aux_var_up%diff(1) + aux_var_dn%diff(1))*&
                                           (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
                                           
    Jdn(1,1) = Jdn(1,1)+ddifff_dp_dn*0.5d0*(aux_var_up%diff(1) + aux_var_dn%diff(1))*&
                                           (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
    Jdn(1,2) = Jdn(1,2)+ddifff_dt_dn*0.5d0*(aux_var_up%diff(1) + aux_var_dn%diff(1))*&
                                           (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
    do ispec=2, option%nflowspec
      Jup(ispec,1) = Jup(ispec,1)+ddifff_dp_up*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*&
                                                     (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
      Jup(ispec,2) = Jup(ispec,2)+ddifff_dt_up*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*&
                                                     (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
      Jup(ispec,ispec+1) = Jup(ispec,ispec+1)+difff*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))

      Jdn(ispec,1) = Jdn(ispec,1)+ddifff_dp_dn*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*&
                                             (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
      Jdn(ispec,2) = Jdn(ispec,2)+ddifff_dt_dn*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*&
                                             (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
      Jdn(ispec,ispec+1) = Jdn(ispec,ispec+1)+difff*0.5d0*(aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))*(-1.d0)
    enddo  
  endif 

! conduction term
        
  Dk = (Dk_up * Dk_dn) / (dd_dn*Dk_up + dd_up*Dk_dn)
!  cond = Dk*area*(aux_var_up%temp-aux_var_dn%temp) 
  Jup(option%nflowdof,2) = Jup(option%nflowdof,2)+Dk*area
  Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2)+Dk*area*(-1.d0)
  Jup = Jup*option%flow_dt
  Jdn = Jdn*option%flow_dt
 ! note: Res is the flux contribution, for node up J = J + Jup
 !                                              dn J = J - Jdn  

  if (option%numerical_derivatives) then
    allocate(aux_var_pert_up%xmol(option%nflowspec),aux_var_pert_up%diff(option%nflowspec))
    allocate(aux_var_pert_dn%xmol(option%nflowspec),aux_var_pert_dn%diff(option%nflowspec))
    call THCAuxVarCopy(aux_var_up,aux_var_pert_up,option)
    call THCAuxVarCopy(aux_var_dn,aux_var_pert_dn,option)
    x_up(1) = aux_var_up%pres
    x_up(2) = aux_var_up%temp
    x_up(3) = aux_var_up%xmol(2)
    x_dn(1) = aux_var_dn%pres
    x_dn(2) = aux_var_dn%temp
    x_dn(3) = aux_var_dn%xmol(2)
    call THCFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                      aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                      area,dist_gravity,upweight, &
                      option,v_darcy,res)
    do ideriv = 1,3
      pert_up = x_up(ideriv)*perturbation_tolerance
      pert_dn = x_dn(ideriv)*perturbation_tolerance
      x_pert_up = x_up
      x_pert_dn = x_dn
      x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
      x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      call THCAuxVarCompute(x_pert_up,aux_var_pert_up,iphase,sat_func_up, &
                                 0.d0,0.d0,option)
      call THCAuxVarCompute(x_pert_dn,aux_var_pert_dn,iphase,sat_func_dn, &
                                 0.d0,0.d0,option)
      call THCFlux(aux_var_pert_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,v_darcy,res_pert_up)
      call THCFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_pert_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,v_darcy,res_pert_dn)
      J_pert_up(:,ideriv) = (res_pert_up(:)-res(:))/pert_up
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    deallocate(aux_var_pert_up%xmol,aux_var_pert_up%diff)
    deallocate(aux_var_pert_dn%xmol,aux_var_pert_dn%diff)
    Jup = J_pert_up
    Jdn = J_pert_dn
  endif

end subroutine THCFluxDerivative

! ************************************************************************** !
!
! THCFlux: Computes the internal flux terms for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** ! 
subroutine THCFlux(aux_var_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        aux_var_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,v_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(thc_auxvar_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: v_darcy,area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec
  PetscReal :: fluxm(option%nflowspec),fluxe,q
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
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
      uh = aux_var_up%h
      uxmol(1:option%nflowspec) = aux_var_up%xmol(1:option%nflowspec)
    else
      ukvr = aux_var_dn%kvr
      uh = aux_var_dn%h
      uxmol(1:option%nflowspec) = aux_var_dn%xmol(1:option%nflowspec)
    endif      
   

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
        
      do ispec=1, option%nflowspec 
        fluxm(ispec)=fluxm(ispec) + q*density_ave*uxmol(ispec)
      enddo
      fluxe = fluxe + q*density_ave*uh 
    endif
  endif 

! Diffusion term   
! Note : average rule may not be correct  
  if ((aux_var_up%sat > eps) .and. (aux_var_dn%sat > eps)) then
   
    difff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)* &
                            (aux_var_up%den+aux_var_dn%den)
    do ispec=1, option%nflowspec
      fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
                 (aux_var_up%diff(ispec) + aux_var_dn%diff(ispec))* &
                 (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
    enddo  
  endif 

! conduction term
        
  Dk = (Dk_up * Dk_dn) / (dd_dn*Dk_up + dd_up*Dk_dn)
  cond = Dk*area*(aux_var_up%temp-aux_var_dn%temp) 
  fluxe=fluxe + cond

  Res(1:option%nflowdof-1) = fluxm(:) * option%flow_dt
  Res(option%nflowdof) = fluxe * option%flow_dt
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine THCFlux

! ************************************************************************** !
!
! THCBCFluxDerivative: Computes the derivatives of the boundary flux 
!                           terms for the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine THCBCFluxDerivative(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
                                    por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                                    area,dist_gravity,option, &
                                    sat_func_dn,Jdn)
  use Option_module
  use Material_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(thc_auxvar_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array in boundary condition
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: area
  type(saturation_function_type) :: sat_func_dn  
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: v_darcy
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi

  PetscReal :: ddiff_dp_dn, ddiff_dt_dn
  PetscReal :: dden_ave_dp_dn, dden_ave_dt_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn, dphi_dt_dn
  PetscReal :: dukvr_dp_dn, dukvr_dt_dn
  PetscReal :: duh_dp_dn, duh_dt_dn
  PetscReal :: dq_dp_dn, dq_dt_dn
  PetscReal :: duxmol_dxmol_dn

  PetscInt :: iphase, ideriv
  type(thc_auxvar_type) :: aux_var_pert_dn, aux_var_pert_up
  PetscReal :: perturbation
  PetscReal :: x_dn(3), x_up(3), x_pert_dn(3), x_pert_up(3), pert_dn, res(3), &
            res_pert_dn(3), J_pert_dn(3,3)
  
  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dden_ave_dt_dn = 0.d0
  ddiff_dp_dn = 0.d0
  ddiff_dt_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dt_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dt_dn = 0.d0
  duh_dp_dn = 0.d0
  duh_dt_dn = 0.d0
  dq_dp_dn = 0.d0
  dq_dt_dn = 0.d0
  duxmol_dxmol_dn = 0.d0
        
  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(THC_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
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
        dden_ave_dt_dn = (1.D0-upweight)*aux_var_dn%dden_dt

        if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
          dden_ave_dt_dn = dden_ave_dt_dn + upweight*aux_var_up%dden_dt
        endif
        
        gravity = (upweight*aux_var_up%den*aux_var_up%avgmw + &
                  (1.D0-upweight)*aux_var_dn%den*aux_var_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*aux_var_dn%avgmw*dist_gravity

        dphi = aux_var_up%pres - aux_var_dn%pres + gravity
        dphi_dp_dn = -1.d0 + dgravity_dden_dn*aux_var_dn%dden_dp
        dphi_dt_dn = dgravity_dden_dn*aux_var_dn%dden_dt

        if (ibndtype(THC_PRESSURE_DOF) == SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. aux_var_up%pres-option%reference_pressure < eps) then
            dphi = 0.d0
            dphi_dp_dn = 0.d0
            dphi_dt_dn = 0.d0
          endif
        endif        
        
        if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
                                   !( dgravity_dden_up                   ) (dden_dt_up)
          dphi_dt_dn = dphi_dt_dn + upweight*aux_var_up%avgmw*dist_gravity*aux_var_up%dden_dt
        endif
        
        if (dphi>=0.D0) then
          ukvr = aux_var_up%kvr
          if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dukvr_dt_dn = aux_var_up%dkvr_dt
          endif
        else
          ukvr = aux_var_dn%kvr
          dukvr_dp_dn = aux_var_dn%dkvr_dp
          dukvr_dt_dn = aux_var_dn%dkvr_dt
        endif      
     
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
          dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
        endif
      endif 

    case(NEUMANN_BC)
      if (dabs(aux_vars(THC_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(THC_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = aux_var_up%den
          if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dden_ave_dt_dn = aux_var_up%dden_dt
          endif
        else 
          density_ave = aux_var_dn%den
          dden_ave_dp_dn = aux_var_dn%dden_dp
          dden_ave_dt_dn = aux_var_dn%dden_dt
        endif 
        q = v_darcy * area
      endif

  end select

  if (v_darcy >= 0.D0) then
    uh = aux_var_up%h
    uxmol(:)=aux_var_up%xmol(1:option%nflowspec)
    if (ibndtype(THC_PRESSURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dp_dn = aux_var_up%dh_dp
    endif
    if (ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dt_dn = aux_var_up%dh_dt
    endif
    if (ibndtype(THC_CONCENTRATION_DOF) == ZERO_GRADIENT_BC) then
      duxmol_dxmol_dn = 1.d0
    endif
  else
    uh = aux_var_dn%h
    duh_dp_dn = aux_var_dn%dh_dp
    duh_dt_dn = aux_var_dn%dh_dt

    uxmol(:)=aux_var_dn%xmol(1:option%nflowspec)
    duxmol_dxmol_dn = 1.d0
  endif      

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(1)
  Jdn(1,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(1)
!  Jdn(1,3:option%nflowdof) = 0.d0
  do ispec=2,option%nflowspec 
    ! based on flux = q*density_ave*uxmol
    Jdn(ispec,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uxmol(ispec)
    Jdn(ispec,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uxmol(ispec)
    Jdn(ispec,ispec+1) = q*density_ave*duxmol_dxmol_dn
  enddo
      ! based on flux = q*density_ave*uh
  Jdn(option%nflowdof,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn
  Jdn(option%nflowdof,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn
!  Jdn(option%nflowdof,3:option%nflowdof) = 0.d0

  ! Diffusion term   
  select case(ibndtype(THC_CONCENTRATION_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC) 
!      if (aux_var_up%sat > eps .and. aux_var_dn%sat > eps) then
!        diff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*(aux_var_up%den+aux_var_dn%den)
!        ddiff_dp_dn = diffdp * 0.25D0*(aux_var_dn%dsat_dp*(aux_var_up%den+aux_var_dn%den)+ &
!                                      (aux_var_up%sat+aux_var_dn%sat)*aux_var_dn%dden_dp)
!        ddiff_dt_dn = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*aux_var_dn%dden_dt
      if (aux_var_dn%sat > eps) then
        diff = diffdp * aux_var_dn%sat*aux_var_dn%den
        ddiff_dp_dn = diffdp * (aux_var_dn%dsat_dp*aux_var_dn%den+ &
                                aux_var_dn%sat*aux_var_dn%dden_dp)
        ddiff_dt_dn = diffdp * aux_var_dn%sat*aux_var_dn%dden_dt
                                    
        Jdn(1,1) = Jdn(1,1)+ddiff_dp_dn*aux_var_dn%diff(1)*&
                                        (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
        Jdn(1,2) = Jdn(1,2)+ddiff_dt_dn*aux_var_dn%diff(1)*&
                                        (aux_var_up%xmol(1) - aux_var_dn%xmol(1))
        do ispec=2, option%nflowspec
          Jdn(ispec,1) = Jdn(ispec,1)+ddiff_dp_dn*aux_var_dn%diff(ispec)*&
                                                (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
          Jdn(ispec,2) = Jdn(ispec,2)+ddiff_dt_dn*aux_var_dn%diff(ispec)*&
                                                (aux_var_up%xmol(ispec) - aux_var_dn%xmol(ispec))
          Jdn(ispec,ispec+1) = Jdn(ispec,ispec+1)+diff*aux_var_dn%diff(ispec)*(-1.d0)
        enddo  
      endif
  end select

  ! Conduction term
  select case(ibndtype(THC_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dk =  Dk_dn / dd_up
      !cond = Dk*area*(aux_var_up%temp-aux_var_dn%temp) 
      Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2)+Dk*area*(-1.d0)
  end select

  Jdn = Jdn * option%flow_dt

  if (option%numerical_derivatives) then
    allocate(aux_var_pert_dn%xmol(option%nflowspec),aux_var_pert_dn%diff(option%nflowspec))
    allocate(aux_var_pert_up%xmol(option%nflowspec),aux_var_pert_up%diff(option%nflowspec))
    call THCAuxVarCopy(aux_var_up,aux_var_pert_up,option)
    call THCAuxVarCopy(aux_var_dn,aux_var_pert_dn,option)
    x_up(1) = aux_var_up%pres
    x_up(2) = aux_var_up%temp
    x_up(3) = aux_var_up%xmol(2)
    x_dn(1) = aux_var_dn%pres
    x_dn(2) = aux_var_dn%temp
    x_dn(3) = aux_var_dn%xmol(2)
    do ideriv = 1,3
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_up(ideriv) = x_dn(ideriv)
      endif
    enddo
    call THCBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
                        por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                        area,dist_gravity,option,v_darcy,res)
    if (ibndtype(THC_PRESSURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(THC_TEMPERATURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(THC_CONCENTRATION_DOF) == ZERO_GRADIENT_BC) then
      x_pert_up = x_up
    endif
    do ideriv = 1,3
      pert_dn = x_dn(ideriv)*perturbation_tolerance    
      x_pert_dn = x_dn
      x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      x_pert_up = x_up
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_pert_up(ideriv) = x_pert_dn(ideriv)
      endif   
      call THCAuxVarCompute(x_pert_dn,aux_var_pert_dn,iphase,sat_func_dn, &
                                 0.d0,0.d0,option)
      call THCAuxVarCompute(x_pert_up,aux_var_pert_up,iphase,sat_func_dn, &
                                 0.d0,0.d0,option)
      call THCBCFlux(ibndtype,aux_vars,aux_var_pert_up,aux_var_pert_dn, &
                          por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                          area,dist_gravity,option,v_darcy,res_pert_dn)
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    deallocate(aux_var_pert_dn%xmol,aux_var_pert_dn%diff)
    Jdn = J_pert_dn
  endif

end subroutine THCBCFluxDerivative

! ************************************************************************** !
!
! THCBCFlux: Computes the  boundary flux terms for the residual
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine THCBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
                          por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
                          area,dist_gravity,option,v_darcy,Res)
  use Option_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(thc_auxvar_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: v_darcy, area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  
  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  ! Flow   
  diffdp = por_dn*tor_dn/dd_up*area
  select case(ibndtype(THC_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
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

        if (ibndtype(THC_PRESSURE_DOF) == SEEPAGE_BC) then
              ! flow in         ! boundary cell is <= pref
          if (dphi > 0.d0 .and. aux_var_up%pres-option%reference_pressure < eps) then
            dphi = 0.d0
          endif
        endif
        
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
      if (dabs(aux_vars(THC_PRESSURE_DOF)) > floweps) then
        v_darcy = aux_vars(THC_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = aux_var_up%den
        else 
          density_ave = aux_var_dn%den
        endif 
      endif

  end select

  q = v_darcy * area

  if (v_darcy >= 0.D0) then
    uh = aux_var_up%h
    uxmol(:)=aux_var_up%xmol(1:option%nflowspec)
  else
    uh = aux_var_dn%h
    uxmol(:)=aux_var_dn%xmol(1:option%nflowspec)
  endif      
    
  do ispec=1, option%nflowspec 
    fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
  enddo
  fluxe = fluxe + q*density_ave*uh

  ! Diffusion term   
  select case(ibndtype(THC_CONCENTRATION_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC) 
!      if (aux_var_up%sat > eps .and. aux_var_dn%sat > eps) then
!        diff = diffdp * 0.25D0*(aux_var_up%sat+aux_var_dn%sat)*(aux_var_up%den+aux_var_dn%den)
      if (aux_var_dn%sat > eps) then
        diff = diffdp * aux_var_dn%sat*aux_var_dn%den
        do ispec = 1, option%nflowspec
          fluxm(ispec) = fluxm(ispec) + diff * aux_var_dn%diff(ispec)* &
                           (aux_var_up%xmol(ispec)-aux_var_dn%xmol(ispec))
        enddo  
      endif
  end select

  ! Conduction term
  select case(ibndtype(THC_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dk =  Dk_dn / dd_up
      cond = Dk*area*(aux_var_up%temp-aux_var_dn%temp) 
      fluxe=fluxe + cond
  end select

  Res(1:option%nflowspec)=fluxm(:)* option%flow_dt
  Res(option%nflowdof)=fluxe * option%flow_dt

end subroutine THCBCFlux

! ************************************************************************** !
!
! THCResidual: Computes the residual equation 
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine THCResidual(snes,xx,r,realization,ierr)

  use Realization_module
  use Level_module
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module

  implicit none

  interface
     subroutine samrpetscobjectstateincrease(vec)
       implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
       Vec :: vec
     end subroutine samrpetscobjectstateincrease
     
  end interface

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  field => realization%field
  discretization => realization%discretization
  
  ! Communication -----------------------------------------
  ! These 3 must be called before THCUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call THCResidualPatch(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if(discretization%itype==AMR_GRID) then
     call samrpetscobjectstateincrease(r)
  endif

end subroutine THCResidual

! ************************************************************************** !
!
! THCResidualPatch: Computes the residual equation at patch level
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine THCResidualPatch(snes,xx,r,realization,ierr)

  use water_eos_module

  use Connection_module
  use Realization_module
  use Patch_module
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

  PetscErrorCode :: ierr
  PetscInt :: i, jn
  PetscInt :: ip1, ip2
  PetscInt :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer ::accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
                          
               
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: iphase
  PetscInt :: icap_up, icap_dn, ithrm_up, ithrm_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: dd, f_up, f_dn, ff
  PetscReal :: perm_up, perm_dn
  PetscReal :: D_up, D_dn  ! "Diffusion" constants at upstream, downstream faces.
  PetscReal :: dw_kg, dw_mol
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy
  PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(thc_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscTruth :: enthalpy_flag
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  aux_vars => patch%aux%THC%aux_vars
  aux_vars_bc => patch%aux%THC%aux_vars_bc
  
  call THCUpdateAuxVarsPatch(realization)
  ! override flags since they will soon be out of date  
  patch%aux%THC%aux_vars_up_to_date = PETSC_FALSE

! now assign access pointer to local variables
  call GridVecGetArrayF90(grid,field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecGetArrayF90(grid, r, r_p, ierr)
  call GridVecGetArrayF90(grid,field%flow_accum, accum_p, ierr)
 
  call GridVecGetArrayF90(grid,field%flow_yy,yy_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecGetArrayF90(grid,field%ithrm_loc, ithrm_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr)
  !print *,' Finished scattering non deriv'

  r_p = 0.d0
#if 1
  ! Accumulation terms ------------------------------------
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    call THCAccumulation(aux_vars(ghosted_id),porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              option%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%flow_condition%num_sub_conditions > THC_CONCENTRATION_DOF) then
      enthalpy_flag = PETSC_TRUE
    else
      enthalpy_flag = PETSC_FALSE
    endif

    qsrc1 = source_sink%flow_condition%pressure%dataset%cur_value(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%cur_value(1)

    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      if (enthalpy_flag) then
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1 * option%flow_dt   
      endif         

      if (qsrc1 > 0.d0) then ! injection
        call wateos_noderiv(tsrc1,aux_vars(ghosted_id)%pres, &
                            dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        r_p((local_id-1)*option%nflowdof + jh2o) = r_p((local_id-1)*option%nflowdof + jh2o) &
                                               - qsrc1 *option%flow_dt
        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - qsrc1*enth_src_h2o*option%flow_dt
      endif  
    
      if (csrc1 > 0.d0) then ! injection
        call printErrMsg(option,"concentration source not yet implemented in THC")
      endif
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif
#if 1
  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle
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

      call THCFlux(aux_vars(ghosted_id_up),porosity_loc_p(ghosted_id_up), &
                          tor_loc_p(ghosted_id_up),option%sir(1,icap_up), &
                          dd_up,perm_up,D_up, &
                        aux_vars(ghosted_id_dn),porosity_loc_p(ghosted_id_dn), &
                          tor_loc_p(ghosted_id_dn),option%sir(1,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                        cur_connection_set%area(iconn),distance_gravity, &
                        upweight,option,v_darcy,Res)

      patch%internal_velocities(1,sum_connection) = v_darcy
     
      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif
   
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
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

      icap_dn = int(icap_loc_p(ghosted_id))  

      call THCBCFlux(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                aux_vars_bc(sum_connection), &
                                aux_vars(ghosted_id), &
                                porosity_loc_p(ghosted_id), &
                                tor_loc_p(ghosted_id), &
                                option%sir(1,icap_dn), &
                                cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
                                cur_connection_set%area(iconn), &
                                distance_gravity,option, &
                                v_darcy,Res)
      patch%boundary_velocities(1,sum_connection) = v_darcy

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif  
  if (option%use_isoth) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      istart = 3 + (local_id-1)*option%nflowdof
      r_p(istart)=xx_loc_p(2 + (ghosted_id-1)*option%nflowdof)-yy_p(istart-1)
    enddo
  endif

  if (patch%aux%THC%inactive_cells_exist) then
    do i=1,patch%aux%THC%n_zero_rows
      r_p(patch%aux%THC%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_yy, yy_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%flow_accum, accum_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)
  call GridVecRestoreArrayF90(grid,field%ithrm_loc, ithrm_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr)

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(option%comm,'THCresidual.out',viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(option%comm,'THCxx.out',viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

end subroutine THCResidualPatch

! ************************************************************************** !
!
! THCJacobian: Computes the Jacobian
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine THCJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_module
  use Patch_module
  use Level_module
  use Grid_module
  use Option_module

  implicit none

  interface
     subroutine SAMRSetCurrentJacobianPatch(mat,patch) 
#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
       
       Mat :: mat
       PetscFortranAddr :: patch
     end subroutine SAMRSetCurrentJacobianPatch
  end interface

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  PetscReal :: norm
  
  flag = SAME_NONZERO_PATTERN
  call MatGetType(A,mat_type,ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr)

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      grid => cur_patch%grid
      ! need to set the current patch in the Jacobian operator
      ! so that entries will be set correctly
      if(associated(grid%structured_grid) .and. &
        (.not.(grid%structured_grid%p_samr_patch.eq.0))) then
         call SAMRSetCurrentJacobianPatch(J, grid%structured_grid%p_samr_patch)
      endif
      call THCJacobianPatch(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(realization%option%comm,'THCjacobian.out', &
                              viewer,ierr)
    call MatView(J,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    call MatNorm(J,NORM_1,norm,ierr)
    if (realization%option%myrank == 0) print *, '1 norm:', norm
    call MatNorm(J,NORM_FROBENIUS,norm,ierr)
    if (realization%option%myrank == 0) print *, '2 norm:', norm
    call MatNorm(J,NORM_INFINITY,norm,ierr)
    if (realization%option%myrank == 0) print *, 'inf norm:', norm
  endif
  
end subroutine THCJacobian

! ************************************************************************** !
!
! THCJacobianPatch: Computes the Jacobian
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine THCJacobianPatch(snes,xx,A,B,flag,realization,ierr)
       
  use water_eos_module

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_module
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_type) :: realization
  MatStructure flag

  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: ithrm_up, ithrm_dn, i
  PetscInt :: ip1, ip2 

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), tor_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
  PetscInt :: icap,iphas,iphas_up,iphas_dn,icap_up,icap_dn
  PetscInt :: ii, jj
  PetscReal :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
  PetscReal :: tsrc1,qsrc1,csrc1,hsrc1
  PetscReal :: dd_up, dd_dn, dd, f_up, f_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: dw_dp,dw_dt,hw_dp,hw_dt,dresT_dp,dresT_dt
  PetscReal :: D_up, D_dn  ! "Diffusion" constants upstream and downstream of a face.
  PetscReal :: zero, norm
  PetscReal :: upweight
  PetscReal :: max_dev  
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
            Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscTruth :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(thc_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  
  PetscViewer :: viewer
  Vec :: debug_vec

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  aux_vars => patch%aux%THC%aux_vars
  aux_vars_bc => patch%aux%THC%aux_vars_bc
  
#if 0
!  call THCNumericalJacobianTest(xx,realization)
#endif

  call GridVecGetArrayF90(grid,field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%volume, volume_p, ierr)

  call GridVecGetArrayF90(grid,field%ithrm_loc, ithrm_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
  call GridVecGetArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr)

#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    icap = int(icap_loc_p(ghosted_id))
    call THCAccumulationDerivative(aux_vars(ghosted_id), &
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
    call PetscViewerASCIIOpen(option%comm,'jacobian_accum.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
    if (source_sink%flow_condition%num_sub_conditions > THC_CONCENTRATION_DOF) then
      enthalpy_flag = PETSC_TRUE
    else
      enthalpy_flag = PETSC_FALSE
    endif

    qsrc1 = source_sink%flow_condition%pressure%dataset%cur_value(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%cur_value(1)

    qsrc1 = qsrc1 / option%fmwh2o ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2
      
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
!      if (enthalpy_flag) then
!        r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1 * option%flow_dt   
!      endif         

      if (qsrc1 > 0.d0) then ! injection
        call wateos(tsrc1,aux_vars(ghosted_id)%pres,dw_kg,dw_mol,dw_dp,dw_dt, &
              enth_src_h2o,hw_dp,hw_dt,option%scale,ierr)        
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        ! base on r_p() = r_p() - qsrc1*enth_src_h2o*option%flow_dt
        dresT_dp = -qsrc1*hw_dp*option%flow_dt
        ! dresT_dt = -qsrc1*hw_dt*option%flow_dt ! since tsrc1 is prescribed, there is no derivative
        istart = ghosted_id*option%nflowdof
        call MatSetValuesLocal(A,1,istart-1,1,istart-option%nflowdof,dresT_dp,ADD_VALUES,ierr)
        ! call MatSetValuesLocal(A,1,istart-1,1,istart-1,dresT_dt,ADD_VALUES,ierr)
      endif  
    
      if (csrc1 > 0.d0) then ! injection
        call printErrMsg(option,"concentration source not yet implemented in THC")
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
    call PetscViewerASCIIOpen(option%comm,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or. &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
   
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
                              
      call THCFluxDerivative(aux_vars(ghosted_id_up),porosity_loc_p(ghosted_id_up), &
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
        Jup = -Jup
        Jdn = -Jdn
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
    call PetscViewerASCIIOpen(option%comm,'jacobian_flux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
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
      icap_dn = int(icap_loc_p(ghosted_id))  

      call THCBCFluxDerivative(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
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
      Jdn = -Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES,ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%comm,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call GridVecRestoreArrayF90(grid,field%flow_xx_loc, xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%tor_loc, tor_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%volume, volume_p, ierr)

   
  call GridVecRestoreArrayF90(grid,field%ithrm_loc, ithrm_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%icap_loc, icap_loc_p, ierr)
  call GridVecRestoreArrayF90(grid,field%iphas_loc, iphase_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL
  zero = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero,ierr) 
  do i=1, n_zero_rows
    ii = mod(zero_rows_local(i),option%nflowdof)
    ip1 = zero_rows_local_ghosted(i)
    if (ii == 0) then
      ip2 = ip1-1
    elseif (ii == option%nflowdof-1) then
      ip2 = ip1+1
    else
      ip2 = ip1
    endif
    call MatSetValuesLocal(A,1,ip1,1,ip2,1.d0,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
#else
  if (patch%aux%THC%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%THC%n_zero_rows, &
                          patch%aux%THC%zero_rows_local_ghosted,f_up,ierr) 
  endif
#endif

end subroutine THCJacobianPatch

! ************************************************************************** !
!
! THCCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine THCCreateZeroArray(patch,option)

  use Patch_module
  use Grid_module
  use Option_module
  
  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  
  PetscInt :: ncount, idof
  PetscInt :: local_id, ghosted_id

  type(grid_type), pointer :: grid
  PetscInt :: flag = 0
  PetscInt :: n_zero_rows
  PetscInt, pointer :: zero_rows_local(:)
  PetscInt, pointer :: zero_rows_local_ghosted(:)
  PetscErrorCode :: ierr
    
  grid => patch%grid
  
  n_zero_rows = 0

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
        n_zero_rows = n_zero_rows + option%nflowdof
      else
#ifdef ISOTHERMAL
        n_zero_rows = n_zero_rows + 1
#endif
      endif
    enddo
  else
#ifdef ISOTHERMAL
    n_zero_rows = n_zero_rows + grid%nlmax
#endif
  endif

  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))

  zero_rows_local = 0
  zero_rows_local_ghosted = 0
  ncount = 0

  if (associated(patch%imat)) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) then
        do idof = 1, option%nflowdof
          ncount = ncount + 1
          zero_rows_local(ncount) = (local_id-1)*option%nflowdof+idof
          zero_rows_local_ghosted(ncount) = (ghosted_id-1)*option%nflowdof+idof-1
        enddo
      else
#ifdef ISOTHERMAL
        ncount = ncount + 1
        zero_rows_local(ncount) = local_id*option%nflowdof
        zero_rows_local_ghosted(ncount) = ghosted_id*option%nflowdof-1
#endif
      endif
    enddo
  else
#ifdef ISOTHERMAL
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ncount = ncount + 1
      zero_rows_local(ncount) = local_id*option%nflowdof
      zero_rows_local_ghosted(ncount) = ghosted_id*option%nflowdof-1
    enddo
#endif
  endif

  patch%aux%THC%zero_rows_local => zero_rows_local
  patch%aux%THC%zero_rows_local_ghosted => zero_rows_local_ghosted
  patch%aux%THC%n_zero_rows = n_zero_rows

  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                     option%comm,ierr)
  if (flag > 0) patch%aux%THC%inactive_cells_exist = PETSC_TRUE

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif

end subroutine THCCreateZeroArray

! ************************************************************************** !
!
! THCMaxChange: Computes the maximum change in the solution vector
! author: Glenn Hammond
! date: 01/15/08
!
! ************************************************************************** !
subroutine THCMaxChange(realization)

  use Realization_module
  use Option_module
  use Field_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field

  option%dcmax=0.D0
  
  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,option%dtmpmax,ierr)
  if (option%nflowdof > 2) &
    call VecStrideNorm(field%flow_dxx,TWO_INTEGER,NORM_INFINITY,option%dcmax,ierr)
    
end subroutine THCMaxChange

! ************************************************************************** !
!
! THCResidualToMass: Computes mass balance from residual equation
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine THCResidualToMass(realization)

  use Realization_module
  use Level_module
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Grid_module

  implicit none

  Vec :: ts_mass_balance
  type(realization_type) :: realization
  
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  PetscReal, pointer :: mass_balance_p(:)
  type(thc_auxvar_type), pointer :: aux_vars(:) 
  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart
  
  option => realization%option
  field => realization%field

  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid
      aux_vars => cur_patch%aux%THC%aux_vars

      call GridVecGetArrayF90(grid,field%flow_ts_mass_balance,mass_balance_p, ierr)
  
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (associated(cur_patch%imat)) then
          if (cur_patch%imat(ghosted_id) <= 0) cycle
        endif
        
        istart = (ghosted_id-1)*option%nflowdof+1
        mass_balance_p(istart) = mass_balance_p(istart)/ &
                                 aux_vars(ghosted_id)%den* &
                                 aux_vars(ghosted_id)%den_kg
      enddo

      call GridVecRestoreArrayF90(grid,field%flow_ts_mass_balance,mass_balance_p, ierr)

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine THCResidualToMass

! ************************************************************************** !
!
! THCLiteGetTecplotHeader: Returns THC contribution to 
!                               Tecplot file header
! author: Glenn Hammond
! date: 02/13/08
!
! ************************************************************************** !
function THCGetTecplotHeader(realization)

  use Realization_module
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: THCGetTecplotHeader
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i
  
  option => realization%option
  field => realization%field
  
  string = ',' // &
           '"T [C]",' // &
           '"P [Pa]",' // &
           '"sl",' // &
           '"Ul"' 
  do i=1,option%nflowspec
    write(string2,'('',"Xl('',i2,'')"'')') i
    string = trim(string) // trim(string2)
  enddo
  
  THCGetTecplotHeader = string

end function THCGetTecplotHeader

! ************************************************************************** !
!
! THCDestroy: Deallocates variables associated with Richard
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine THCDestroy(patch)

  use Patch_module

  implicit none
  
  type(patch_type) :: patch
  
  ! need to free array in aux vars
  call THCAuxDestroy(patch%aux%THC)

end subroutine THCDestroy

end module THC_module
