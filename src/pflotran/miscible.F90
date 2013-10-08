module Miscible_module
  
  use Miscible_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"
  
!#include "include/petscf90.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
!#ifdef USE_PETSC216
!#include "finclude/petscsles.h"
!#endif
#include "finclude/petscsnes.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsysdef.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petsclog.h"
#include "finclude/petscerror.h"

! Cutoff parameters
  PetscReal, parameter :: formeps = 1.D-4
  PetscReal, parameter :: eps = 1.D-5 
  PetscReal, parameter :: dfac = 1D-8
  PetscReal, parameter :: floweps = 1.D-24
  PetscReal, parameter :: zerocut =0.D0 

  PetscInt, parameter :: jh2o=1, jglyc=2
  
  public MiscibleResidual,MiscibleJacobian, &
         MiscibleUpdateFixedAccumulation,MiscibleTimeCut, &
         MiscibleSetup, &
         MiscibleMaxChange, MiscibleUpdateSolution, &
         MiscibleGetTecplotHeader, MiscibleInitializeTimestep, &
         MiscibleUpdateAuxVars,MiscibleComputeMassBalance

contains

! ************************************************************************** !
!
! MiscibleTimeCut: Resets arrays for time step cut
! author: Chuan Lu
! date: 9/13/08
!
! ************************************************************************** !
subroutine MiscibleTimeCut(realization)
 
  use Realization_class
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  PetscReal, pointer :: xx_p(:),yy_p(:)
  PetscErrorCode :: ierr
  PetscInt :: local_id

  option => realization%option
  field => realization%field

  call VecCopy(field%flow_yy,field%flow_xx,ierr)

end subroutine MiscibleTimeCut

! ************************************************************************** !
!
! MiscibleSetup: 
! author: Chuan Lu
! date: 9/13/08
!
! ************************************************************************** !
subroutine MiscibleSetup(realization)

  use Realization_class
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
      call MiscibleSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  call MiscibleSetPlotVariables(realization)

end subroutine MiscibleSetup

! ************************************************************************** !
!
! MiscibleSetupPatch: Creates arrays for auxiliary variables
! author: Chuan Lu
! date: 10/1/08
!
! ************************************************************************** !
subroutine MiscibleSetupPatch(realization)

  use Realization_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: ghosted_id, iconn, sum_connection, ipara
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  patch%aux%Miscible => MiscibleAuxCreate()
  
!  option%io_buffer = 'Before Miscible can be run, the Miscible_parameter object ' // &
!                     'must be initialized with the proper variables ' // &
!                     'MiscibleAuxCreate() is called anyhwere.'
!  call printErrMsg(option)
     
! Miscible_parameters create *********************************************
  allocate(patch%aux%Miscible%Miscible_parameter%sir(option%nphase, &
                                  size(realization%saturation_function_array)))
  do ipara = 1, size(realization%saturation_function_array)
    patch%aux%Miscible%Miscible_parameter%sir(:,realization%saturation_function_array(ipara)%ptr%id) = &
      realization%saturation_function_array(ipara)%ptr%Sr(:)
  enddo
  
! dencpr  
  allocate(patch%aux%Miscible%Miscible_parameter%dencpr(size(realization%material_property_array)))
  do ipara = 1, size(realization%material_property_array)
    patch%aux%Miscible%Miscible_parameter%dencpr(realization%material_property_array(ipara)%ptr%id) = &
      realization%material_property_array(ipara)%ptr%rock_density*option%scale*&
      realization%material_property_array(ipara)%ptr%specific_heat
  enddo
  
! ckwet
  allocate(patch%aux%Miscible%Miscible_parameter%ckwet(size(realization%material_property_array)))
  do ipara = 1, size(realization%material_property_array)
    patch%aux%Miscible%Miscible_parameter%ckwet(realization%material_property_array(ipara)%ptr%id) = &
      realization%material_property_array(ipara)%ptr%thermal_conductivity_wet*option%scale
  enddo
! Miscible_parameters create_end *****************************************

! allocate aux_var data structures for all grid cells  
  allocate(aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call MiscibleAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%aux%Miscible%aux_vars => aux_vars
  patch%aux%Miscible%num_aux = grid%ngmax
  
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
    call MiscibleAuxVarInit(aux_vars_bc(iconn),option)
  enddo
  patch%aux%Miscible%aux_vars_bc => aux_vars_bc
  patch%aux%Miscible%num_aux_bc = sum_connection
  option%numerical_derivatives_flow = PETSC_TRUE
  
  allocate(patch%aux%Miscible%delx(option%nflowdof, grid%ngmax))
  allocate(patch%aux%Miscible%Resold_AR(grid%nlmax,option%nflowdof))
  allocate(patch%aux%Miscible%Resold_BC(grid%nlmax,option%nflowdof))
  ! should be allocated by the number of BC connections, just for debug now
  allocate(patch%aux%Miscible%Resold_FL(ConnectionGetNumberInList(patch%grid%&
           internal_connection_set_list),option%nflowdof))
  
  ! create zero array for zeroing residual and Jacobian (1 on diagonal)
  ! for inactive cells (and isothermal)
  call MiscibleCreateZeroArray(patch,option)

end subroutine MiscibleSetupPatch

! ************************************************************************** !
!
! MiscibleComputeMassBalance: 
!                        
! author: Jitendra Kumar 
! date: 07/21/2010
! Adapted from RichardsComputeMassBalance: need to be checked
! ************************************************************************** !
subroutine MiscibleComputeMassBalance(realization,mass_balance)

  use Realization_class
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec,1)
   
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  mass_balance = 0.d0

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MiscibleComputeMassBalancePatch(realization,mass_balance)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MiscibleComputeMassBalance    

! ************************************************************************** !
!
! MiscibleComputeMassBalancePatch: 
!                        
! author: Jitendra Kumar 
! date: 07/21/2010
! Adapted from RichardsComputeMassBalancePatch: need to be checked
! ************************************************************************** !
subroutine MiscibleComputeMassBalancePatch(realization,mass_balance)
 
  use Realization_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec,1)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(miscible_auxvar_type), pointer :: miscible_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  PetscReal, pointer :: volume_p(:), porosity_loc_p(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ispec
  PetscInt :: jh2o=1,jglyc=2
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  global_aux_vars => patch%aux%Global%aux_vars
  miscible_aux_vars => patch%aux%Miscible%aux_vars

  call VecGetArrayF90(field%volume,volume_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    ! mass = saturation * density * mole fraction * volume
    do ispec = 1,option%nflowspec
      mass_balance(ispec,1) = mass_balance(ispec,1) + &
        miscible_aux_vars(ghosted_id)%aux_var_elem(0)%xmol(ispec)* &
        global_aux_vars(ghosted_id)%den(1)* &
        porosity_loc_p(ghosted_id)*volume_p(local_id)
    enddo
  enddo
  mass_balance(jh2o,1) = mass_balance(jh2o,1)*FMWH2O
  mass_balance(jglyc,1) = mass_balance(jglyc,1)*FMWGLYC

  call VecRestoreArrayF90(field%volume,volume_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  
end subroutine MiscibleComputeMassBalancePatch

! ************************************************************************** !
!
! MiscibleZeroMassBalDeltaPatch: Zeros mass balance delta array
! author: Satish Karra, LANL
! date: 12/13/11
!
! ************************************************************************** !
subroutine MiscibleZeroMassBalDeltaPatch(realization)
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Miscible%num_aux
    patch%aux%Global%aux_vars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Miscible%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Miscible%num_aux_bc
      global_aux_vars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  if (patch%aux%Miscible%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Miscible%num_aux_ss
      global_aux_vars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
 
end subroutine MiscibleZeroMassBalDeltaPatch

! ************************************************************************** !
!
! MiscibleUpdateMassBalancePatch: Updates mass balance
! author: Glenn Hammond
! date: 12/13/11
!
! ************************************************************************** !
subroutine MiscibleUpdateMassBalancePatch(realization)
 
  use Realization_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Miscible%num_aux
    patch%aux%Global%aux_vars(iconn)%mass_balance = &
      patch%aux%Global%aux_vars(iconn)%mass_balance + &
      patch%aux%Global%aux_vars(iconn)%mass_balance_delta* &
      option%flow_dt
  enddo
#endif

  if (patch%aux%Miscible%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Miscible%num_aux_bc
      global_aux_vars_bc(iconn)%mass_balance = &
        global_aux_vars_bc(iconn)%mass_balance + &
        global_aux_vars_bc(iconn)%mass_balance_delta*option%flow_dt
    enddo
  endif

  if (patch%aux%Miscible%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Miscible%num_aux_ss
      global_aux_vars_bc(iconn)%mass_balance = &
        global_aux_vars_bc(iconn)%mass_balance + &
        global_aux_vars_bc(iconn)%mass_balance_delta*option%flow_dt
    enddo
  endif


end subroutine MiscibleUpdateMassBalancePatch

! ************************************************************************** !
! MiscibleInitGuessCheckPatch: 
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
  function MiscibleInitGuessCheck(realization)
 
  use Realization_class
  use Level_module
  use Patch_module
  use Option_module
  
  PetscInt ::  MiscibleInitGuessCheck
  type(realization_type) :: realization
  type(option_type), pointer:: option
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscInt :: ipass, ipass0
  PetscErrorCode :: ierr

  option => realization%option
  cur_level => realization%level_list%first
  ipass = 1
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      ipass= MiscibleInitGuessCheckPatch(realization)
      if(ipass<=0)then
        nullify(cur_level)
        exit 
      endif
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

   call MPI_Barrier(option%mycomm,ierr)
   if(option%mycommsize >1)then
      call MPI_Allreduce(ipass,ipass0,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                         option%mycomm,ierr)
      if(ipass0 < option%mycommsize) ipass=-1
   endif
   MiscibleInitGuessCheck =ipass
 end function MiscibleInitGuessCheck



! ************************************************************************** !
! MiscibleInitGuessCheckPatch: 
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
  function MiscibleInitGuessCheckPatch(realization)
   
    use co2_span_wagner_module
     
    use Realization_class
    use Patch_module
    use Field_module
    use Grid_module
    use Option_module
    implicit none
    
    PetscInt :: MiscibleInitGuessCheckPatch 
    type(realization_type) :: realization
    type(grid_type), pointer :: grid
    type(patch_type), pointer :: patch
    type(option_type), pointer :: option
    type(field_type), pointer :: field
      
    PetscInt :: local_id, ghosted_id, ipass
    PetscErrorCode :: ierr
    PetscReal, pointer :: xx_p(:)

    patch => realization%patch
    grid => patch%grid
    option => realization%option
    field => realization%field
    
    call VecGetArrayF90(field%flow_xx,xx_p, ierr)
    
    ipass=1
    MiscibleInitGuessCheckPatch = ipass
  end function MiscibleInitGuessCheckPatch

! ***************************************************************************
!
! MiscibleUpdateAuxVars: Updates the auxiliary variables associated with 
!                        the Miscible problem
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine MiscibleUpdateAuxVars(realization)

  use Realization_class
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
      call MiscibleUpdateAuxVarsPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MiscibleUpdateAuxVars

! ************************************************************************** !
!
! MiscibleUpdateAuxVarsPatch: Updates the auxiliary variables associated with 
!                        the Miscible problem
! author: Chuan Lu
! date: 12/10/07
!
! ************************************************************************** !
subroutine MiscibleUpdateAuxVarsPatch(realization)

  use Realization_class
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
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscReal :: mnacl, ynacl, xphi
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  aux_vars => patch%aux%Miscible%aux_vars
  aux_vars_bc => patch%aux%Miscible%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    if(.not. associated(realization%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr))then
       print*, 'error!!! saturation function not allocated', ghosted_id,icap_loc_p(ghosted_id)
    endif
    
    call MiscibleAuxVarCompute_NINC(xx_loc_p(istart:iend), &
                       aux_vars(ghosted_id)%aux_var_elem(0), &
                       global_aux_vars(ghosted_id), &
                       realization%fluid_properties,option)
                      
!   update global variables
    if( associated(global_aux_vars))then
      global_aux_vars(ghosted_id)%pres = aux_vars(ghosted_id)%aux_var_elem(0)%pres
      global_aux_vars(ghosted_id)%sat(:) = 1D0
      global_aux_vars(ghosted_id)%den(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:)
      global_aux_vars(ghosted_id)%den_kg(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:) &
                                          *aux_vars(ghosted_id)%aux_var_elem(0)%avgmw(:)
    else
      print *,'Not associated global for Miscible'
    endif

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
          case(DIRICHLET_BC,HYDROSTATIC_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(:) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
        end select
      enddo
 
      call MiscibleAuxVarCompute_NINC(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0), &
                         global_aux_vars_bc(sum_connection), &
                         realization%fluid_properties, option)

      if (associated(global_aux_vars_bc)) then
        global_aux_vars_bc(sum_connection)%pres(:)= aux_vars_bc(sum_connection)%aux_var_elem(0)%pres
        global_aux_vars_bc(sum_connection)%sat(:)=1.D0
        global_aux_vars_bc(sum_connection)%den(:)=aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:)
        global_aux_vars_bc(sum_connection)%den_kg = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:) &
            * aux_vars_bc(sum_connection)%aux_var_elem(0)%avgmw(:)
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  
  patch%aux%Miscible%aux_vars_up_to_date = PETSC_TRUE

end subroutine MiscibleUpdateAuxVarsPatch

! ************************************************************************** !
!
! MiscibleInitializeTimestep: Update data in module prior to time step
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine MiscibleInitializeTimestep(realization)

  use Realization_class
  
  implicit none
  
  type(realization_type) :: realization

  call MiscibleUpdateFixedAccumulation(realization)

end subroutine MiscibleInitializeTimestep

! ************************************************************************** !
!
! MiscibleUpdateSolution: Updates data in module after a successful time step
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
subroutine MiscibleUpdateSolution(realization)

  use Realization_class
  use Field_module
  use Level_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  PetscErrorCode :: ierr
  
  field => realization%field
  
  call VecCopy(realization%field%flow_xx,realization%field%flow_yy,ierr)   
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do 
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MiscibleUpdateSolutionPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MiscibleUpdateSolution


! ************************************************************************** !
!
! MiscibleUpdateSolutionPatch: Updates data in module after a successful time 
!                             step 
! author: Satish Karra, LANL
! written based on RichardsUpdateSolutionPatch
! date: 08/23/11
!
! ************************************************************************** !
subroutine MiscibleUpdateSolutionPatch(realization)

  use Realization_class
    
  implicit none
  
  type(realization_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call MiscibleUpdateMassBalancePatch(realization)
  endif

end subroutine MiscibleUpdateSolutionPatch

! ************************************************************************** !
!
! MiscibleUpdateFixedAccumulation: Updates the fixed portion of the 
!                                  accumulation term
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine MiscibleUpdateFixedAccumulation(realization)

  use Realization_class
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
      call MiscibleUpdateFixedAccumPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine MiscibleUpdateFixedAccumulation

! ************************************************************************** !
!
! MiscibleUpdateFixedAccumPatch: Updates the fixed portion of the 
!                                  accumulation term
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine MiscibleUpdateFixedAccumPatch(realization)

  use Realization_class
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
  type(Miscible_parameter_type), pointer :: Miscible_parameter
  type(Miscible_auxvar_type), pointer :: aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)

  PetscInt :: ghosted_id, local_id, istart, iend, iphase
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tortuosity_loc_p(:), volume_p(:), &
                        ithrm_loc_p(:), accum_p(:)
                          
  PetscErrorCode :: ierr
  
  call MiscibleUpdateAuxVarsPatch(realization) 
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  global_aux_vars => patch%aux%Global%aux_vars
  Miscible_parameter => patch%aux%Miscible%Miscible_parameter
  aux_vars => patch%aux%Miscible%aux_vars
    
  call VecGetArrayF90(field%flow_xx,xx_p, ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecGetArrayF90(field%tortuosity_loc,tortuosity_loc_p,ierr)
  call VecGetArrayF90(field%volume,volume_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

    call MiscibleAccumulation(aux_vars(ghosted_id)%aux_var_elem(0), &
                              global_aux_vars(ghosted_id), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              Miscible_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
  call VecRestoreArrayF90(field%tortuosity_loc,tortuosity_loc_p,ierr)
  call VecRestoreArrayF90(field%volume,volume_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)

end subroutine MiscibleUpdateFixedAccumPatch


! ************************************************************************** !
!
! MiscibleAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !  
subroutine MiscibleAccumulation(aux_var,global_aux_var,por,vol,rock_dencpr,option,Res)

  use Option_module
  
  implicit none

  type(Miscible_auxvar_elem_type) :: aux_var
  type(option_type) :: option
  type(global_auxvar_type) :: global_aux_var
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: vol,por,rock_dencpr
     
  PetscInt :: ispec, np
  PetscReal :: porXvol, mol(option%nflowspec)

  porXvol = por*vol
  mol=0.d0
  do np = 1, option%nphase
    do ispec = 1, option%nflowspec  
      mol(ispec) = mol(ispec) + &
        aux_var%den(np) * aux_var%xmol(ispec + (np-1)*option%nflowspec)
    enddo
  enddo
  mol = mol*porXvol
  Res(1:option%nflowspec) = mol(:)

end subroutine MiscibleAccumulation

! ************************************************************************** !
!
! MiscibleAccumulation: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !  
subroutine MiscibleAccumulation_Xp(aux_var,global_aux_var,por,vol,rock_dencpr,option,Res)

use Option_module

implicit none

type(Miscible_auxvar_elem_type) :: aux_var
type(option_type) :: option
type(global_auxvar_type) :: global_aux_var
PetscReal :: Res(1:option%nflowdof) 
PetscReal :: vol,por,rock_dencpr

PetscInt :: ispec, np
PetscReal :: porXvol, mol(option%nflowspec)

porXvol = por*vol
mol = 0.d0
np = 1

! pressure equation
ispec = 1
mol(ispec) = mol(ispec) + aux_var%den(np)

!PPG
ispec = 2
mol(ispec) = mol(ispec) + &
aux_var%den(np) * aux_var%xmol(ispec + (np-1)*option%nflowspec)

mol = mol*porXvol
Res(1:option%nflowspec) = mol(:)

end subroutine MiscibleAccumulation_Xp

! ************************************************************************** !
!
! MiscibleSourceSink: Computes the non-fixed portion of the accumulation
!                       term for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !  
subroutine MiscibleSourceSink(mmsrc,nsrcpara,psrc,tsrc,hsrc,csrc,aux_var,isrctype,Res,&
                            qsrc_phase,energy_flag, option)

  use Option_module
  
  use Water_EOS_module
  use co2eos_module
  use co2_span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use co2_span_wagner_module 
  
  implicit none

  type(Miscible_auxvar_elem_type) :: aux_var
  type(option_type) :: option
  PetscReal Res(1:option%nflowdof) 
  PetscReal, pointer :: mmsrc(:)
  PetscReal psrc(option%nphase),tsrc,hsrc, csrc 
  PetscInt isrctype, nsrcpara
  PetscBool :: energy_flag
  PetscReal :: qsrc_phase(:) 
     
  PetscReal, pointer :: msrc(:)
  PetscReal dw_kg, dw_mol,dddt,dddp
  PetscReal :: enth_src_h2o, enth_src_co2 
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: ukvr, v_darcy, dq, dphi
  PetscReal :: well_status, well_diameter
  PetscReal :: pressure_bh, well_factor, pressure_max, pressure_min
  PetscReal :: well_inj_water, well_inj_co2
  PetscInt  :: np
  PetscInt :: iflag
  PetscErrorCode :: ierr
  
  Res=0D0
  allocate(msrc(nsrcpara))
  msrc = mmsrc(1:nsrcpara)
  qsrc_phase = 0.d0
 ! if (present(ireac)) iireac=ireac
!  if (energy_flag) then
!    Res(option%nflowdof) = Res(option%nflowdof) + hsrc * option%flow_dt   
!  endif         
 
  select case(isrctype)
    case(MASS_RATE_SS)
      msrc(1) = msrc(1) / FMWH2O
      msrc(2) = msrc(2) / FMWGLYC
      if (msrc(1) /= 0.d0 .or. msrc(2) /= 0.d0) then ! H2O injection
!        call wateos_noderiv(tsrc,aux_var%pres,dw_kg,dw_mol,enth_src_h2o,option%scale,ierr)
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
        Res(jh2o) = Res(jh2o) + msrc(1)*option%flow_dt
        Res(jglyc) = Res(jglyc) + msrc(2)*option%flow_dt
      endif  
#if 0   
 !  End of mass rate inplementation
    case(WELL_SS) ! production well
     !if node pessure is lower than the given extraction pressure, shut it down
    ! Flow term
!  well parameter explaination
!   1. well status. 1 injection; -1 production; 0 shut in
!                   2 rate controled injection (same as rate_ss, with max pressure control, not completed yet) 
!                  -2 rate controled production(not implemented for now) 
!
!   2. well factor [m^3],  the effective permeability [m^2/s]
!   3. bottomhole pressure:  [Pa]
!   4. max pressure: [Pa]
!   5. min pressure: [Pa]   
!   6. preferred mass flux of water [kg/s]
!   7. preferred mass flux of Co2 [kg/s]
!   8. well diameter, not used now
!   9. skin factor, not used now

      well_status = msrc(1)
      well_factor = msrc(2)
      pressure_bh = msrc(3)
      pressure_max = msrc(4)
      pressure_min = msrc(5)
      well_inj_water = msrc(6)
      well_inj_co2 = msrc(7)
    
!     if(pressure_min < 0D0) pressure_min = 0D0 !not limited by pressure lower bound   

    ! production well (well status = -1)
      if (dabs(well_status + 1D0) < 1D-1) then 
        if (aux_var%pres > pressure_min) then
          Dq = well_factor 
          do np = 1, option%nphase
            dphi = aux_var%pres - aux_var%pc(np) - pressure_bh
            if (dphi>=0.D0) then ! outflow only
              ukvr = aux_var%kvr(np)
              if(ukvr<1e-20) ukvr=0D0
              v_darcy=0D0
              if (ukvr*Dq>floweps) then
                v_darcy = Dq * ukvr * dphi
                ! store volumetric rate for ss_fluid_fluxes()
                qsrc_phase(1) = -1.d0*v_darcy
                Res(1) = Res(1) - v_darcy* aux_var%den(np)* &
                  aux_var%xmol((np-1)*option%nflowspec+1)*option%flow_dt
                Res(2) = Res(2) - v_darcy* aux_var%den(np)* &
                  aux_var%xmol((np-1)*option%nflowspec+2)*option%flow_dt
                if(energy_flag) Res(3) = Res(3) - v_darcy * aux_var%den(np)* &
                  aux_var%h(np)*option%flow_dt
              ! print *,'produce: ',np,v_darcy
              endif
            endif
          enddo
        endif
      endif 
     !print *,'well-prod: ',  aux_var%pres,psrc(1), res
    ! injection well (well status = 2)
      if (dabs(well_status - 2D0) < 1D-1) then 

        call wateos_noderiv(tsrc,aux_var%pres,dw_kg,dw_mol,enth_src_h2o, &
          option%scale,ierr)

        Dq = msrc(2) ! well parameter, read in input file
                      ! Take the place of 2nd parameter 
        ! Flow term
        if (aux_var%pres < pressure_max) then  
          do np = 1, option%nphase
            dphi = pressure_bh - aux_var%pres + aux_var%pc(np)
            if (dphi>=0.D0) then ! outflow only
              ukvr = aux_var%kvr(np)
              v_darcy=0.D0
              if (ukvr*Dq>floweps) then
                v_darcy = Dq * ukvr * dphi
                ! store volumetric rate for ss_fluid_fluxes()
                qsrc_phase(1) = v_darcy
                Res(1) = Res(1) + v_darcy* aux_var%den(np)* &
!                 aux_var%xmol((np-1)*option%nflowspec+1) * option%flow_dt
                  (1.d0-csrc) * option%flow_dt
                Res(2) = Res(2) + v_darcy* aux_var%den(np)* &
!                 aux_var%xmol((np-1)*option%nflowspec+2) * option%flow_dt
                  csrc * option%flow_dt
!               if(energy_flag) Res(3) = Res(3) + v_darcy*aux_var%den(np)*aux_var%h(np)*option%flow_dt
                if(energy_flag) Res(3) = Res(3) + v_darcy*aux_var%den(np)* &
                  enth_src_h2o*option%flow_dt
                
!               print *,'inject: ',np,v_darcy
              endif
            endif
          enddo
        endif
      endif
#endif          
    case default
    print *,'Unrecognized Source/Sink condition: ', isrctype 
  end select      
!  deallocate(msrc)
  
end subroutine MiscibleSourceSink


! ************************************************************************** !
!
! MiscibleFlux: Computes the internal flux terms for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** ! 
subroutine MiscibleFlux(aux_var_up,por_up,tor_up,dd_up,perm_up, &
                        aux_var_dn,por_dn,tor_dn,dd_dn,perm_dn, &
                        area,dist_gravity,upweight, &
                        option,vv_darcy,Res)
  use Option_module                              
  
  implicit none
  
  type(Miscible_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: vv_darcy(:),area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec, np, ind
  PetscReal :: fluxm(option%nflowspec),q,v_darcy
  PetscReal :: uxmol(1:option%nflowspec),ukvr,Dq
  PetscReal :: upweight,density_ave,gravity,dphi
  PetscReal :: portor_up,portor_dn,dendif_up,dendif_dn,dif_harmonic

  Dq = (perm_up*perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
  
! harmonic average porosity and tortuosity
  portor_up = por_up*tor_up
  portor_dn = por_dn*tor_dn
  
  fluxm = 0.D0
  vv_darcy = 0.D0 
  
  do np = 1, option%nphase
      
!   Flow term
      
    density_ave = upweight*aux_var_up%den(np) + (1.D0-upweight)*aux_var_dn%den(np) 
        
    gravity = (upweight*aux_var_up%den(np)*aux_var_up%avgmw(np) + &
             (1.D0-upweight)*aux_var_dn%den(np)*aux_var_dn%avgmw(np)) &
             *dist_gravity

    dphi = aux_var_up%pres - aux_var_dn%pres + gravity

    v_darcy = 0.D0
    ukvr = 0.D0
    uxmol = 0.D0

!   note uxmol only contains one phase xmol
    if (dphi >= 0.D0) then
      ukvr = aux_var_up%kvr(np)
      uxmol(:) = aux_var_up%xmol((np-1)*option%nflowspec+1:np*option%nflowspec)
    else
      ukvr = aux_var_dn%kvr(np)
      uxmol(:) = aux_var_dn%xmol((np-1)*option%nflowspec+1:np*option%nflowspec)
    endif

    v_darcy = Dq * ukvr * dphi
    vv_darcy(np) = v_darcy
    q = v_darcy * area
    do ispec = 1, option%nflowspec
      fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
    enddo  

!   Diffusion term   
!   Note : use harmonic average rule for diffusion
    do ispec=1, option%nflowspec
      ind = ispec + (np-1)*option%nflowspec
      dendif_up = aux_var_up%den(1)*aux_var_up%diff(ispec)
      dendif_dn = aux_var_dn%den(1)*aux_var_dn%diff(ispec)
      dif_harmonic = (portor_up*portor_dn*dendif_up*dendif_dn) &
        /(dd_dn*portor_up*dendif_up + dd_up*portor_dn*dendif_dn)*area
      fluxm(ispec) = fluxm(ispec) + dif_harmonic* &
        (aux_var_up%xmol(ind) - aux_var_dn%xmol(ind))
    enddo
  enddo

 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  
  Res(1:option%nflowspec) = fluxm(:)*option%flow_dt

end subroutine MiscibleFlux

! ************************************************************************** !
!
! MiscibleBCFlux: Computes the  boundary flux terms for the residual
! author: Chuan Lu
! date: 10/12/08
!
! ************************************************************************** !
subroutine MiscibleBCFlux(ibndtype,aux_vars,aux_var_up,aux_var_dn, &
     por_dn,tor_dn,dd_up,perm_dn, &
     area,dist_gravity,option,vv_darcy,Res)
  use Option_module
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(Miscible_auxvar_elem_type) :: aux_var_up, aux_var_dn
  type(option_type) :: option
  PetscReal :: dd_up
  PetscReal :: aux_vars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,tor_dn
  PetscReal :: vv_darcy(:), area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec, np
  PetscReal :: fluxm(option%nflowspec),q,density_ave, v_darcy
  PetscReal :: uxmol(1:option%nflowspec),ukvr,diff,portor,Dq
  PetscReal :: gravity,dphi
  
  fluxm = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  portor = por_dn*tor_dn/dd_up*area

  ! Flow   
  do np = 1, option%nphase   ! note: nphase = 1
    select case(ibndtype(1))
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
        Dq = perm_dn / dd_up
        ! Flow term
        ukvr=0.D0
        v_darcy=0.D0

        density_ave = aux_var_up%den(np)
        
        gravity = aux_var_up%den(np)*aux_var_up%avgmw(np)*dist_gravity
       
        dphi = aux_var_up%pres - aux_var_dn%pres &
                - aux_var_up%pc(np) + aux_var_dn%pc(np) &
                + gravity
   
        ! figure out the direction of flow
        if (dphi >= 0.D0) then
          ukvr = aux_var_up%kvr(np)
        else
          ukvr = aux_var_dn%kvr(np)
        endif
     
        if (ukvr*Dq > floweps) then
          v_darcy = Dq * ukvr * dphi
        endif

      case(NEUMANN_BC) !may not work
        v_darcy = 0.D0
        if (dabs(aux_vars(1)) > floweps) then
          v_darcy = aux_vars(MIS_PRESSURE_DOF)
          if (v_darcy > 0.d0) then 
            density_ave = aux_var_up%den(np)
          else 
            density_ave = aux_var_dn%den(np)
          endif
        endif
    end select
     
    q = v_darcy*area
    vv_darcy(np) = v_darcy
    uxmol = 0.D0
     
    if (v_darcy >= 0.D0) then
      uxmol(:) = aux_var_up%xmol((np-1)*option%nflowspec+1:np*option%nflowspec)
    else
      uxmol(:) = aux_var_dn%xmol((np-1)*option%nflowspec+1:np*option%nflowspec)
    endif
    
    do ispec=1, option%nflowspec
      fluxm(ispec) = fluxm(ispec) + q*density_ave*uxmol(ispec)
    end do 
  enddo

    ! Diffusion term   
  select case(ibndtype(2))
    case(DIRICHLET_BC) 
      do np = 1, option%nphase ! note: nphase = 1
        diff = portor*aux_var_up%den(np)
        do ispec = 1, option%nflowspec
          fluxm(ispec) = fluxm(ispec) + diff* &
                   aux_var_dn%diff((np-1)*option%nflowspec+ispec)* &
                   (aux_var_up%xmol((np-1)*option%nflowspec+ispec) &
                   -aux_var_dn%xmol((np-1)*option%nflowspec+ispec))
        enddo
      enddo
  end select

  Res(1:option%nflowspec) = fluxm(:)*option%flow_dt

end subroutine MiscibleBCFlux

! ************************************************************************** !
!
! MiscibleResidual: Computes the residual equation 
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine MiscibleResidual(snes,xx,r,realization,ierr)

  use Realization_class
  use Level_module
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Grid_module 
  use Logging_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscInt :: ichange  

  field => realization%field
  grid => realization%patch%grid
  option => realization%option
  discretization => realization%discretization
  
  call PetscLogEventBegin(logging%event_r_residual,ierr)
 
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)

 ! check initial guess -----------------------------------------------
  ierr = MiscibleInitGuessCheck(realization)
  if(ierr<0)then
    !ierr = PETSC_ERR_ARG_OUTOFRANGE
    if (option%myrank==0) print *,'table out of range: ',ierr
    call SNESSetFunctionDomainError(snes,ierr) 
    return
  endif 
  ! end check ---------------------------------------------------------

  ! Communication -----------------------------------------
  ! These 3 must be called before MiscibleUpdateAuxVars()
!  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc,field%ithrm_loc,ONEDOF)

!  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)


! pass #0 prepare numerical increment  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MiscibleResidualPatch0(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

! pass #1 internal and boundary flux terms
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MiscibleResidualPatch1(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

! pass #2 for everything else
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MiscibleResidualPatch2(snes,xx,r,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if (realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rresidual.out', &
                              viewer,ierr)
    call VecView(r,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rxx.out', &
                              viewer,ierr)
    call VecView(xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call PetscLogEventEnd(logging%event_r_residual,ierr)

end subroutine MiscibleResidual

! ************************************************************************** !
!
! MiscibleResidualPatch1: Computes the Residual by Flux
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine MiscibleResidualPatch1(snes,xx,r,realization,ierr)

  use Connection_module
  use Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
  implicit none

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
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
               xx_loc_p(:), tortuosity_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
                          
               
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: iphase
  PetscInt :: icap_up, icap_dn, ithrm_up, ithrm_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: dd, f_up, f_dn, ff
  PetscReal :: perm_up, perm_dn
  PetscReal :: D_up, D_dn  ! "Diffusion" constants at upstream, downstream faces.
  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy(realization%option%nphase)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(Miscible_parameter_type), pointer :: Miscible_parameter
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: ng
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  Miscible_parameter => patch%aux%Miscible%Miscible_parameter
  aux_vars => patch%aux%Miscible%aux_vars
  aux_vars_bc => patch%aux%Miscible%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

 ! call MiscibleUpdateAuxVarsPatchNinc(realization)
  ! override flags since they will soon be out of date  
 ! patch%MiscibleAux%aux_vars_up_to_date = PETSC_FALSE 
  
  if (option%compute_mass_balance_new) then
    call MiscibleZeroMassBalDeltaPatch(realization)
  endif

! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90( r, r_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)

  r_p = 0.d0
 
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
      D_dn = Miscible_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))  
! Then need fill up increments for BCs
      do idof =1, option%nflowdof   
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(HYDROSTATIC_BC)
            xxbc(1) = boundary_condition%flow_aux_real_var(1,iconn)
            if (idof >= 2) then
              xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            endif 
          case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo

 
    call MiscibleAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),&
           global_aux_vars_bc(sum_connection),&
           realization%fluid_properties, option)

    if( associated(global_aux_vars_bc))then
      global_aux_vars_bc(sum_connection)%pres(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%pres - &
                     aux_vars(ghosted_id)%aux_var_elem(0)%pc(:)
      global_aux_vars_bc(sum_connection)%sat(:) = 1.D0
      global_aux_vars_bc(sum_connection)%den(:) = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:)
      global_aux_vars_bc(sum_connection)%den_kg = aux_vars_bc(sum_connection)%aux_var_elem(0)%den(:) &
          * aux_vars_bc(sum_connection)%aux_var_elem(0)%avgmw(:)
  !   global_aux_vars(ghosted_id)%den_kg_store
    endif

    call MiscibleBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         aux_vars_bc(sum_connection)%aux_var_elem(0), &
         aux_vars(ghosted_id)%aux_var_elem(0), &
         porosity_loc_p(ghosted_id), &
         tortuosity_loc_p(ghosted_id), &
         cur_connection_set%dist(0,iconn),perm_dn, &
         cur_connection_set%area(iconn), &
         distance_gravity,option, &
         v_darcy,Res)
    patch%boundary_velocities(:,sum_connection) = v_darcy(:)
    patch%aux%Miscible%Resold_BC(local_id,1:option%nflowdof) = &
    patch%aux%Miscible%ResOld_BC(local_id,1:option%nflowdof) - Res(1:option%nflowdof)

    if (option%compute_mass_balance_new) then
      ! contribution to boundary
      global_aux_vars_bc(sum_connection)%mass_balance_delta(:,1) = &
        global_aux_vars_bc(sum_connection)%mass_balance_delta(:,1) &
          - Res(:)/option%flow_dt 
    endif
  
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)
  enddo
  boundary_condition => boundary_condition%next
 enddo

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
                         dot_product(option%gravity, &
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

      call MiscibleFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
        tortuosity_loc_p(ghosted_id_up), &
        dd_up,perm_up, &
        aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
        tortuosity_loc_p(ghosted_id_dn), &
        dd_dn,perm_dn, &
        cur_connection_set%area(iconn),distance_gravity, &
        upweight,option,v_darcy,Res)

      patch%internal_velocities(:,sum_connection) = v_darcy(:)
      patch%aux%Miscible%Resold_FL(sum_connection,1:option%nflowdof)= Res(1:option%nflowdof)

      if (local_id_up > 0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif
   
      if (local_id_dn > 0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo    

  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90( r, r_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)

end subroutine MiscibleResidualPatch1

! ************************************************************************** !
!
! MiscibleJacobianPatch0: Computes the Residual Aux vars for numerical Jacobin
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine MiscibleResidualPatch0(snes,xx,r,realization,ierr)

  use Connection_module
  use Realization_class
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
               tortuosity_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
                          
               
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:)

  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof)


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(Miscible_parameter_type), pointer :: Miscible_parameter
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)
  PetscBool :: enthalpy_flag
  PetscInt :: ng
  PetscInt :: iconn, idof, istart, iend
  PetscReal, pointer :: delx(:)
  PetscReal :: logx2
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  Miscible_parameter => patch%aux%Miscible%Miscible_parameter
  aux_vars => patch%aux%Miscible%aux_vars
  aux_vars_bc => patch%aux%Miscible%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

 ! call MiscibleUpdateAuxVarsPatchNinc(realization)
  ! override flags since they will soon be out of date  
 ! patch%MiscibleAux%aux_vars_up_to_date = PETSC_FALSE 

! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)

  allocate(delx(option%nflowdof))

  patch%aux%Miscible%Resold_AR=0.D0
  patch%aux%Miscible%Resold_BC=0.D0
  patch%aux%Miscible%ResOld_FL=0.D0

! Multiphase flash calculation is more expensive, so calculate once per iteration

  ! Pertubations for aux terms --------------------------------
  do ng = 1, grid%ngmax
    if (grid%nG2L(ng)<0) cycle
    if (associated(patch%imat)) then
      if (patch%imat(ng) <= 0) cycle
    endif
    ghosted_id = ng   
    istart = (ng-1) * option%nflowdof + 1; iend = istart - 1 + option%nflowdof
     ! iphase =int(iphase_loc_p(ng))
    call MiscibleAuxVarCompute_Ninc(xx_loc_p(istart:iend),aux_vars(ng)%aux_var_elem(0),&
          global_aux_vars(ng),&
          realization%fluid_properties,option)
!    print *,'mis: Respatch0 ', xx_loc_p(istart:iend),aux_vars(ng)%aux_var_elem(0)%den

    if(associated(global_aux_vars)) then
      global_aux_vars(ghosted_id)%pres(:) = aux_vars(ghosted_id)%aux_var_elem(0)%pres
      global_aux_vars(ghosted_id)%sat(:) = 1D0
      global_aux_vars(ghosted_id)%den(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:)
      global_aux_vars(ghosted_id)%den_kg(:) = aux_vars(ghosted_id)%aux_var_elem(0)%den(:) &
            * aux_vars(ghosted_id)%aux_var_elem(0)%avgmw(:)
    else
      print *,'Not associated global for Miscible'
    endif

    if (option%numerical_derivatives_flow) then
      delx(1) = xx_loc_p((ng-1)*option%nflowdof+1)*dfac * 1.D-3

!     print *,'mis_res_p: ',delx(1),xx_loc_p((ng-1)*option%nflowdof+1)

!      linear formulation
#if 0
      do idof = 2, option%nflowdof
        if(xx_loc_p((ng-1)*option%nflowdof+idof) <= 0.9) then
          delx(idof) = dfac*xx_loc_p((ng-1)*option%nflowdof+idof)*1.d1 
        else
          delx(idof) = -dfac*xx_loc_p((ng-1)*option%nflowdof+idof)*1D1 
        endif

        if(delx(idof) <  1D-8 .and. delx(idof) >= 0.D0) delx(idof) = 1D-8
        if(delx(idof) > -1D-8 .and. delx(idof) <  0.D0) delx(idof) =-1D-8

        if((delx(idof)+xx_loc_p((ng-1)*option%nflowdof+idof)) > 1.D0) then
          delx(idof) = (1.D0-xx_loc_p((ng-1)*option%nflowdof+idof))*1D-4
        endif
        if((delx(idof)+xx_loc_p((ng-1)*option%nflowdof+idof)) < 0.D0) then
          delx(idof) = xx_loc_p((ng-1)*option%nflowdof+idof)*1D-4
        endif

!       print *,'mis_res_x: ',delx(idof),xx_loc_p((ng-1)*option%nflowdof+idof)
      enddo
#endif

!      log formulation
!#if 0
      do idof = 2, option%nflowdof
        logx2 = xx_loc_p((ng-1)*option%nflowdof+idof)
!         x2 = exp(logx2)
        if (logx2 <= 0.d0) then
          delx(idof) = 1.d-6 
        else
          delx(idof) = -1.d-6 
        endif

!       print *,'mis_res_x: ',delx(idof),xx_loc_p((ng-1)*option%nflowdof+idof)
      enddo
!#endif

!      store increments
      patch%aux%Miscible%delx(:,ng) = delx(:)

      call MiscibleAuxVarCompute_Winc(xx_loc_p(istart:iend),delx(:),&
            aux_vars(ng)%aux_var_elem(1:option%nflowdof),global_aux_vars(ng),&
            realization%fluid_properties,option)
    endif
  enddo

  deallocate(delx)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)

end subroutine MiscibleResidualPatch0

! ************************************************************************** !
!
! MiscibleResidualPatch2: Computes other terms in Residual
!                       (accumulation, source/sink, reaction)
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine MiscibleResidualPatch2(snes,xx,r,realization,ierr)

  use Connection_module
  use Realization_class
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
  PetscInt :: local_id, ghosted_id
  
  PetscReal, pointer ::accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:)
               
  PetscReal, pointer :: ithrm_loc_p(:)

  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: Res(realization%option%nflowdof), v_darcy(realization%option%nphase)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscReal :: psrc(1:realization%option%nphase)
  PetscViewer :: viewer
  PetscInt :: nsrcpara
  PetscReal, pointer :: msrc(:)

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(Miscible_parameter_type), pointer :: Miscible_parameter
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars_ss(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: ng
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity

! PetscReal, pointer :: iphase_loc_p(:) 

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  Miscible_parameter => patch%aux%Miscible%Miscible_parameter
  aux_vars => patch%aux%Miscible%aux_vars
  aux_vars_bc => patch%aux%Miscible%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc
  global_aux_vars_ss => patch%aux%Global%aux_vars_ss
  
 ! call MiscibleUpdateAuxVarsPatchNinc(realization)
  ! override flags since they will soon be out of date  
 ! patch%MiscibleAux%aux_vars_up_to_date = PETSC_FALSE 

! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)

! call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr); 

  ! Accumulation terms (include reaction------------------------------------
  if (.not.option%steady_state) then

    r_p = r_p - accum_p

    do local_id = 1, grid%nlmax  ! For each local node do...

      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1

      call MiscibleAccumulation(aux_vars(ghosted_id)%aux_var_elem(0),&
                            global_aux_vars(ghosted_id), &
                            porosity_loc_p(ghosted_id), &
                            volume_p(local_id), &
                            Miscible_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                            option,Res) 
      r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)

      patch%aux%Miscible%Resold_AR(local_id, :) = &
      patch%aux%Miscible%Resold_AR(local_id, :)+ Res(1:option%nflowdof)
    enddo
  endif

! call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr); 

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0 
  do 
    if (.not.associated(source_sink)) exit
    !print *, 'RES s/s begin'
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
      enthalpy_flag = PETSC_TRUE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif
      
    if (associated(source_sink%flow_condition%pressure)) then
      psrc(:) = source_sink%flow_condition%pressure%dataset%rarray(:)
    endif 
!    qsrc1 = source_sink%flow_condition%pressure%dataset%rarray(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%rarray(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%rarray(1)
!    hsrc1=0D0
!    qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
!    csrc1 = csrc1 / FMWCO2
!    msrc(1)=qsrc1; msrc(2) =csrc1
    select case(source_sink%flow_condition%itype(1))
      case(MASS_RATE_SS)
        msrc => source_sink%flow_condition%rate%dataset%rarray
        nsrcpara= 2
      case(WELL_SS)
        msrc => source_sink%flow_condition%well%dataset%rarray
        nsrcpara = 7 + option%nflowspec 
      case default
        print *, 'Flash mode does not support source/sink type: ', source_sink%flow_condition%itype(1)
        stop  
    end select

    cur_connection_set => source_sink%connection_set
       
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1 
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      
      call MiscibleSourceSink(msrc,nsrcpara,psrc,tsrc1,hsrc1,csrc1, &
                            aux_vars(ghosted_id)%aux_var_elem(0), &
                            source_sink%flow_condition%itype(1),Res, &
                            patch%ss_fluid_fluxes(:,sum_connection), &
                            enthalpy_flag,option)
      if (option%compute_mass_balance_new) then
        global_aux_vars_ss(sum_connection)%mass_balance_delta(:,1) = &
          global_aux_vars_ss(sum_connection)%mass_balance_delta(:,1) - &
          Res(:)/option%flow_dt
      endif
      r_p((local_id-1)*option%nflowdof + jh2o) = r_p((local_id-1)*option%nflowdof + jh2o) - Res(jh2o)
      r_p((local_id-1)*option%nflowdof + jglyc) = r_p((local_id-1)*option%nflowdof + jglyc) - Res(jglyc)
      patch%aux%Miscible%Resold_AR(local_id,jh2o) = patch%aux%Miscible%Resold_AR(local_id,jh2o) - Res(jh2o)    
      patch%aux%Miscible%Resold_AR(local_id,jglyc) = patch%aux%Miscible%Resold_AR(local_id,jglyc) - Res(jglyc)    
      if (enthalpy_flag)then
        r_p( local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - Res(option%nflowdof)
        patch%aux%Miscible%Resold_AR(local_id,option%nflowdof)=&
          patch%aux%Miscible%Resold_AR(local_id,option%nflowdof) - Res(option%nflowdof)
      endif 
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
  
! adjust residual to R/dt
  select case (option%idt_switch) 
    case(1) 
      r_p(:) = r_p(:)/option%flow_dt
    case(-1)
      if (option%flow_dt > 1.D0) r_p(:) = r_p(:)/option%flow_dt
  end select
  
  do local_id = 1, grid%nlmax
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    endif

!    scale residual by grid cell volume
    istart = 1 + (local_id-1)*option%nflowdof
!   if (volume_p(local_id) > 1.D0) r_p(istart:istart+2) = &
!     r_p(istart:istart+2)/volume_p(local_id)
    if (volume_p(local_id) > 1.D0) r_p(istart:istart+1) = &
      r_p(istart:istart+1)/volume_p(local_id)
!   r_p(istart:istart+1) = r_p(istart:istart+1)/volume_p(local_id)
    if(r_p(istart) > 1.E10 .or. r_p(istart) < -1.E10) then
      ghosted_id = grid%nL2G(local_id)
      print *, 'overflow in res: ', &
      local_id,ghosted_id,istart,r_p (istart:istart+1), &
      aux_vars(ghosted_id)%aux_var_elem(0)%pres, &
      aux_vars(ghosted_id)%aux_var_elem(0)%xmol(1), &
      aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2)
    endif
  enddo

! if (option%use_isothermal) then
!   do local_id = 1, grid%nlmax  ! For each local node do...
!     ghosted_id = grid%nL2G(local_id)   ! corresponding ghost index
!     if (associated(patch%imat)) then
!       if (patch%imat(ghosted_id) <= 0) cycle
!     endif
!     istart = 3 + (local_id-1)*option%nflowdof
!     r_p(istart) = 0.D0 ! xx_loc_p(2 + (ng-1)*option%nflowdof) - yy_p(p1-1)
!   enddo
! endif
 
  if (patch%aux%Miscible%inactive_cells_exist) then
    do i=1,patch%aux%Miscible%n_zero_rows
      r_p(patch%aux%Miscible%zero_rows_local(i)) = 0.d0
    enddo
  endif
 
  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
 
end subroutine MiscibleResidualPatch2


! ************************************************************************** !
!
! MiscibleJacobian: Computes the Jacobian
! author: Chuan Lu
! date: 10/10/08
!
! ************************************************************************** !
subroutine MiscibleJacobian(snes,xx,A,B,flag,realization,ierr)

  use Realization_class
  use Patch_module
  use Level_module
  use Grid_module
  use Option_module
  use Logging_module
  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B, J
  MatType :: mat_type
  type(realization_type) :: realization
  MatStructure flag
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid

  call PetscLogEventBegin(logging%event_r_jacobian,ierr)

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

 ! pass #1 for internal and boundary flux terms
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MiscibleJacobianPatch1(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

! pass #2 for everything else
 cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call MiscibleJacobianPatch2(snes,xx,J,J,flag,realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo


  if (realization%debug%matview_Jacobian) then
#if 1  
    call PetscViewerASCIIOpen(realization%option%mycomm,'Rjacobian.out', &
                              viewer,ierr)
#else
    call PetscViewerBinaryOpen(realization%option%mycomm,'Rjacobian.bin', &
                               FILE_MODE_WRITE,viewer,ierr)
#endif    
    call MatView(J,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

#if 0
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_FROBENIUS,norm,ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_INFINITY,norm,ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option) 
  endif
#endif

  call PetscLogEventEnd(logging%event_r_jacobian,ierr)

end subroutine MiscibleJacobian


! ************************************************************************** !
!
! MiscibleJacobianPatch1: Computes the Jacobian: Flux term
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
subroutine MiscibleJacobianPatch1(snes,xx,A,B,flag,realization,ierr)

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  
  implicit none

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
                          xx_loc_p(:), tortuosity_loc_p(:),&
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
  PetscInt :: natural_id_up, natural_id_dn
  
  PetscReal :: Jup(1:realization%option%nflowdof,1:realization%option%nflowdof), &
            Jdn(1:realization%option%nflowdof,1:realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: Res(realization%option%nflowdof) 
  PetscReal :: xxbc(1:realization%option%nflowdof), delxbc(1:realization%option%nflowdof)
  PetscReal :: ResInc(realization%patch%grid%nlmax,realization%option%nflowdof,&
           realization%option%nflowdof)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(Miscible_parameter_type), pointer :: Miscible_parameter
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)

  PetscReal :: vv_darcy(realization%option%nphase), voltemp
  PetscReal :: ra(1:realization%option%nflowdof,1:realization%option%nflowdof*2) 
  PetscReal :: msrc(1:realization%option%nflowspec)
  PetscReal :: psrc(1:realization%option%nphase)
  PetscReal :: dddt, dddp, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt,&
               dvdp, xphi
  PetscInt :: iphasebc                
  
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

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  Miscible_parameter => patch%aux%Miscible%Miscible_parameter
  aux_vars => patch%aux%Miscible%aux_vars
  aux_vars_bc => patch%aux%Miscible%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
!  flag = SAME_NONZERO_PATTERN

#if 0
!  call MiscibleNumericalJacobianTest(xx,realization)
#endif

 ! print *,'*********** In Jacobian ********************** '
 ! MatzeroEntries has been called in MiscibleJacobin ! clu removed on 11/04/2010 
 !  call MatZeroEntries(A,ierr)

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
!  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)

 ResInc = 0.D0

! Boundary conditions

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
      D_dn = Miscible_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      icap_dn = int(icap_loc_p(ghosted_id))

! Then need fill up increments for BCs
      delxbc=0.D0;
      do idof =1, option%nflowdof   
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
            delxbc(idof)=0.D0
          case(HYDROSTATIC_BC)
            xxbc(1) = boundary_condition%flow_aux_real_var(1,iconn)
            if(idof >= 2)then
              xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
              delxbc(idof)=patch%aux%Miscible%delx(idof,ghosted_id)
            endif 
          case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            delxbc(idof)=patch%aux%Miscible%delx(idof,ghosted_id)
        end select
      enddo
 
      call MiscibleAuxVarCompute_Ninc(xxbc,aux_vars_bc(sum_connection)%aux_var_elem(0),&
         global_aux_vars_bc(sum_connection),&
         realization%fluid_properties, option)
      call MiscibleAuxVarCompute_Winc(xxbc,delxbc,&
         aux_vars_bc(sum_connection)%aux_var_elem(1:option%nflowdof),&
         global_aux_vars_bc(sum_connection),&
         realization%fluid_properties,option)
    
      do nvar=1,option%nflowdof
        call MiscibleBCFlux(boundary_condition%flow_condition%itype, &
          boundary_condition%flow_aux_real_var(:,iconn), &
          aux_vars_bc(sum_connection)%aux_var_elem(nvar), &
          aux_vars(ghosted_id)%aux_var_elem(nvar), &
          porosity_loc_p(ghosted_id), &
          tortuosity_loc_p(ghosted_id), &
          cur_connection_set%dist(0,iconn),perm_dn, &
          cur_connection_set%area(iconn), &
          distance_gravity,option, &
          vv_darcy,Res)
        ResInc(local_id,1:option%nflowdof,nvar) = ResInc(local_id,1:option%nflowdof,nvar) - Res(1:option%nflowdof)
      enddo
    enddo
    boundary_condition => boundary_condition%next
  enddo

! Set matrix values related to single node terms: Accumulation, Source/Sink, BC
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    ra=0.D0
    do neq=1, option%nflowdof
      do nvar=1, option%nflowdof
        ra(neq,nvar) = (ResInc(local_id,neq,nvar) - patch%aux%Miscible%ResOld_BC(local_id,neq))&
          /patch%aux%Miscible%delx(nvar,ghosted_id)

!       print *,'jacobian: ',neq,nvar,local_id,ghosted_id,grid%nlmax,ra(neq,nvar),ResInc(local_id,neq,nvar), &
!         patch%aux%Miscible%ResOld_BC(local_id,neq),patch%aux%Miscible%delx(nvar,ghosted_id)
      enddo
    enddo
   
    select case(option%idt_switch)
      case(1) 
        ra(1:option%nflowdof,1:option%nflowdof) = ra(1:option%nflowdof,1:option%nflowdof) / option%flow_dt
      case(-1)
        if (option%flow_dt > 1) ra(1:option%nflowdof,1:option%nflowdof) = ra(1:option%nflowdof,1:option%nflowdof) / option%flow_dt
    end select

    Jup = ra(1:option%nflowdof,1:option%nflowdof)
    if (volume_p(local_id) > 1.D0) Jup = Jup / volume_p(local_id)
!   Jup = Jup / volume_p(local_id)
   
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)
  end do

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_bcflux.out',viewer,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  ResInc = 0.D0
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
                         dot_product(option%gravity, &
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
      D_up = Miscible_parameter%ckwet(ithrm_up)
      D_dn = Miscible_parameter%ckwet(ithrm_dn)
    
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
      
      do nvar = 1, option%nflowdof 
        call MiscibleFlux(aux_vars(ghosted_id_up)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_up), &
            tortuosity_loc_p(ghosted_id_up), &
            dd_up,perm_up, &
            aux_vars(ghosted_id_dn)%aux_var_elem(0),porosity_loc_p(ghosted_id_dn), &
            tortuosity_loc_p(ghosted_id_dn), &
            dd_dn,perm_dn, &
            cur_connection_set%area(iconn),distance_gravity, &
            upweight, option, vv_darcy, Res)
        ra(:,nvar) = (Res(:)-patch%aux%Miscible%ResOld_FL(iconn,:)) &
              /patch%aux%Miscible%delx(nvar,ghosted_id_up)
        call MiscibleFlux(aux_vars(ghosted_id_up)%aux_var_elem(0),porosity_loc_p(ghosted_id_up), &
            tortuosity_loc_p(ghosted_id_up), &
            dd_up,perm_up, &
            aux_vars(ghosted_id_dn)%aux_var_elem(nvar),porosity_loc_p(ghosted_id_dn),&
            tortuosity_loc_p(ghosted_id_dn), &
            dd_dn,perm_dn, &
            cur_connection_set%area(iconn),distance_gravity, &
            upweight, option, vv_darcy, Res)
         ra(:,nvar+option%nflowdof) = (Res(:)-patch%aux%Miscible%ResOld_FL(iconn,:))&
           /patch%aux%Miscible%delx(nvar,ghosted_id_dn)
      enddo

      select case(option%idt_switch)
        case(1)
          ra = ra / option%flow_dt
        case(-1)  
          if(option%flow_dt > 1.d0) ra = ra / option%flow_dt
      end select
    
      if (local_id_up > 0) then
        voltemp = 1.D0
        if (volume_p(local_id_up) > 1.D0) voltemp = 1.D0/volume_p(local_id_up)
        Jup(:,1:option%nflowdof) = ra(:,1:option%nflowdof)*voltemp !11
        jdn(:,1:option%nflowdof) = ra(:,1 + option%nflowdof:2*option%nflowdof)*voltemp !12

        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr)
      endif
      if (local_id_dn > 0) then
        voltemp = 1.D0
        if (volume_p(local_id_dn) > 1.D0) voltemp = 1.D0/volume_p(local_id_dn)
        Jup(:,1:option%nflowdof)= -ra(:,1:option%nflowdof)*voltemp !21
        jdn(:,1:option%nflowdof)= -ra(:,1 + option%nflowdof:2*option%nflowdof)*voltemp !22
 
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_flux.out',viewer,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)

end subroutine MiscibleJacobianPatch1

! ************************************************************************** !
!
! MiscibleJacobianPatch2: Computes the Jacobian: Accum, source, reaction
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
subroutine MiscibleJacobianPatch2(snes,xx,A,B,flag,realization,ierr)

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  
  implicit none

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
                          xx_loc_p(:), tortuosity_loc_p(:),&
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
  PetscInt ::  natural_id_up,natural_id_dn
  
  PetscReal :: Jup(1:realization%option%nflowdof,1:realization%option%nflowdof), &
            Jdn(1:realization%option%nflowdof,1:realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: Res(realization%option%nflowdof) 
  PetscReal :: xxbc(1:realization%option%nflowdof), delxbc(1:realization%option%nflowdof)
  PetscReal :: ResInc(realization%patch%grid%nlmax,realization%option%nflowdof,&
           realization%option%nflowdof)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(Miscible_parameter_type), pointer :: Miscible_parameter
  type(Miscible_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:)

  PetscReal :: vv_darcy(realization%option%nphase), voltemp
  PetscReal :: ra(1:realization%option%nflowdof,1:realization%option%nflowdof*2) 
  PetscReal, pointer :: msrc(:)
  PetscReal :: psrc(1:realization%option%nphase), ss_flow(1:realization%option%nphase)
  PetscReal :: dddt, dddp, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt,&
               dvdp, xphi
  PetscInt :: nsrcpara, flow_pc                
  
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

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  Miscible_parameter => patch%aux%Miscible%Miscible_parameter
  aux_vars => patch%aux%Miscible%aux_vars
  aux_vars_bc => patch%aux%Miscible%aux_vars_bc
  global_aux_vars => patch%aux%Global%aux_vars
  global_aux_vars_bc => patch%aux%Global%aux_vars_bc

! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s
!  flag = SAME_NONZERO_PATTERN

#if 0
!  call MiscibleNumericalJacobianTest(xx,realization)
#endif

 ! print *,'*********** In Jacobian ********************** '
!  call MatZeroEntries(A,ierr)

  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr)
! call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  ResInc = 0.D0
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
     
    do nvar =1, option%nflowdof
      call MiscibleAccumulation(aux_vars(ghosted_id)%aux_var_elem(nvar), &
             global_aux_vars(ghosted_id),& 
             porosity_loc_p(ghosted_id), &
             volume_p(local_id), &
             Miscible_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
             option,res) 
      ResInc( local_id,:,nvar) =  ResInc(local_id,:,nvar) + Res(:)
    enddo
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sinks%first
  sum_connection = 0 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
      enthalpy_flag = PETSC_TRUE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif

    if (associated(source_sink%flow_condition%pressure)) then
      psrc(:) = source_sink%flow_condition%pressure%dataset%rarray(:)
    endif
    tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%rarray(1)
 !   hsrc1=0.D0
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%rarray(1)

   ! qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
   ! csrc1 = csrc1 / FMWCO2
    select case(source_sink%flow_condition%itype(1))
      case(MASS_RATE_SS)
        msrc => source_sink%flow_condition%rate%dataset%rarray
        nsrcpara= 2
      case(WELL_SS)
        msrc => source_sink%flow_condition%well%dataset%rarray
        nsrcpara = 7 + option%nflowspec 
      case default
        print *, 'Flash mode does not support source/sink type: ', source_sink%flow_condition%itype(1)
        stop  
    end select
    cur_connection_set => source_sink%connection_set
 
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
!     if (enthalpy_flag) then
!       r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1 * option%flow_dt   
!     endif         
      do nvar =1, option%nflowdof
        call MiscibleSourceSink(msrc,nsrcpara,psrc,tsrc1,hsrc1,csrc1, aux_vars(ghosted_id)%aux_var_elem(nvar),&
                            source_sink%flow_condition%itype(1), Res,&
                            ss_flow, &
                            enthalpy_flag, option)

        ResInc(local_id,jh2o,nvar)=  ResInc(local_id,jh2o,nvar) - Res(jh2o)
        ResInc(local_id,jglyc,nvar)=  ResInc(local_id,jglyc,nvar) - Res(jglyc)
        if (enthalpy_flag) & 
          ResInc(local_id,option%nflowdof,nvar) = &
          ResInc(local_id,option%nflowdof,nvar) - Res(option%nflowdof)
      enddo 
    enddo
    source_sink => source_sink%next
  enddo
#endif

! Set matrix values related to single node terms: Accumulation, Source/Sink, BC
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    ra=0.D0
    max_dev=0.D0
    do neq=1, option%nflowdof
      do nvar=1, option%nflowdof
        ra(neq,nvar)=(ResInc(local_id,neq,nvar)-patch%aux%Miscible%ResOld_AR(local_id,neq))&
        /patch%aux%Miscible%delx(nvar,ghosted_id)
      enddo
    enddo
   
    select case(option%idt_switch)
      case(1) 
        ra(1:option%nflowdof,1:option%nflowdof) = &
          ra(1:option%nflowdof,1:option%nflowdof) / option%flow_dt
      case(-1)
        if(option%flow_dt > 1.d0) ra(1:option%nflowdof,1:option%nflowdof) = &
          ra(1:option%nflowdof,1:) / option%flow_dt
    end select

    Jup = ra(1:option%nflowdof,1:option%nflowdof)
    if (volume_p(local_id) > 1.D0) Jup = Jup / volume_p(local_id)
!   Jup = Jup / volume_p(local_id)
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES,ierr)
  end do

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_srcsink.out',viewer,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr)
! call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
 ! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
#if 0
! zero out isothermal and inactive cells
#ifdef ISOTHERMAL
  zero = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
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
#endif
#endif

  if (patch%aux%Miscible%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%Miscible%n_zero_rows, &
                          patch%aux%Miscible%zero_rows_local_ghosted,f_up, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr) 
  endif

  if (realization%debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(option%mycomm,'Rjacobian.out',viewer,ierr)
    call MatView(A,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    call MatNorm(A,NORM_1,norm,ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(A,NORM_INFINITY,norm,ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option)
!    call GridCreateVector(grid,ONEDOF,debug_vec,GLOBAL)
!    call MatGetRowMaxAbs(A,debug_vec,PETSC_NULL_INTEGER,ierr)
!    call VecMax(debug_vec,i,norm,ierr)
!    call VecDestroy(debug_vec,ierr)
  endif

end subroutine MiscibleJacobianPatch2


! ************************************************************************** !
!
! MiscibleCreateZeroArray: Computes the zeroed rows for inactive grid cells
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
subroutine MiscibleCreateZeroArray(patch,option)

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
! print *,'zero rows=', n_zero_rows
  allocate(zero_rows_local(n_zero_rows))
  allocate(zero_rows_local_ghosted(n_zero_rows))
! print *,'zero rows allocated' 
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
      endif
    enddo
  endif
!print *,'zero rows point 1'
  patch%aux%Miscible%n_zero_rows = n_zero_rows
!print *,'zero rows point 2'
  patch%aux%Miscible%zero_rows_local => zero_rows_local
!print *,'zero rows point 3'  
  patch%aux%Miscible%zero_rows_local_ghosted => zero_rows_local_ghosted
!print *,'zero rows point 4'
  call MPI_Allreduce(n_zero_rows,flag,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX, &
                     option%mycomm,ierr)
  if (flag > 0) patch%aux%Miscible%inactive_cells_exist = PETSC_TRUE

  if (ncount /= n_zero_rows) then
    print *, 'Error:  Mismatch in non-zero row count!', ncount, n_zero_rows
    stop
  endif
! print *,'zero rows', flag
end subroutine MiscibleCreateZeroArray

! ************************************************************************** !
!
! MiscibleMaxChange: Computes the maximum change in the solution vector
! author: Chuan Lu
! date: 01/15/08
!
! ************************************************************************** !
subroutine MiscibleMaxChange(realization)

  use Realization_class
  use Level_module
  use Patch_module
  use Field_module
  use Option_module
  use Field_module

  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscReal :: dcmax
  PetscInt :: idof
  PetscErrorCode :: ierr 

  option => realization%option
  field => realization%field

  cur_level => realization%level_list%first
  option%dpmax=0.D0
  option%dtmpmax=0.D0 
  option%dcmax=0.D0
  option%dsmax=0.D0
  dcmax=0.D0

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy,ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,option%dpmax,ierr)

  do idof = 1,option%nflowdof-1 
    dcmax = 0.D0
    call VecStrideNorm(field%flow_dxx,idof,NORM_INFINITY,dcmax,ierr)
    if (dcmax > option%dcmax) option%dcmax = dcmax
  enddo
  
end subroutine MiscibleMaxChange

! ************************************************************************** !
!
! MiscibleGetTecplotHeader: Returns Richards contribution to 
!                               Tecplot file header
! author: Chuan Lu
! date: 10/13/08
!
! ************************************************************************** !
function MiscibleGetTecplotHeader(realization, icolumn)

  use Realization_class
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: MiscibleGetTecplotHeader
  type(realization_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i
  
  option => realization%option
  field => realization%field
  
  string = ''

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-d(l)"'')') icolumn
  else
    write(string2,'('',"d(l)"'')')
  endif
  string = trim(string) // trim(string2)

    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-vis(l)"'')') icolumn
  else
    write(string2,'('',"vis(l)"'')')
  endif
  string = trim(string) // trim(string2)

  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xl('',i2,'')"'')') icolumn, i
    else
      write(string2,'('',"Xl('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo


  MiscibleGetTecplotHeader = string

end function MiscibleGetTecplotHeader

! ************************************************************************** !
!
! MiscibleSetPlotVariables: Adds variables to be printed to list
! author: Glenn Hammond
! date: 10/15/12
!
! ************************************************************************** !
subroutine MiscibleSetPlotVariables(realization)
  
  use Realization_class
  use Output_Aux_module
  use Variables_module

  implicit none

  type(realization_type) :: realization
  type(output_variable_type) :: output_variable
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_list_type), pointer :: list
  
  list => realization%output_option%output_variable_list
  
  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               TEMPERATURE)
  
  name = 'Liquid Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               LIQUID_PRESSURE)

  name = 'Liquid Density'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)

  name = 'Liquid Viscosity'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_VISCOSITY)

  name = 'Liquid Mole Fraction H2O'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION,ONE_INTEGER)

  name = 'Liquid Mole Fraction CO2'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION,TWO_INTEGER)

end subroutine MiscibleSetPlotVariables


end module Miscible_module
