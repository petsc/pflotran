#ifdef SURFACE_FLOW

module Surface_Flow_module

  use Global_Aux_module
  
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

  public SurfaceFlowSetup, &
         SurfaceFlowReadRequiredCardsFromInput, &
         SurfaceFlowRead, &
         SurfaceFlowResidual, &
         SurfaceFlowJacobian
  
contains

! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowSetup(surf_realization)

  use Surface_Realization_module
  
  type(surface_realization_type) :: surf_realization
  
end subroutine SurfaceFlowSetup

! ************************************************************************** !
!> This routine reads required surface flow data from the input file
!! grids.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
subroutine SurfaceFlowReadRequiredCardsFromInput(surf_realization,input,option)

  use Option_module
  use Input_module
  use String_module
  use Surface_Material_module
  use Surface_Realization_module
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module
  use Discretization_module
  use Region_module
  use Condition_module

  implicit none

  type(surface_realization_type)               :: surf_realization
  type(discretization_type),pointer            :: discretization
  type(grid_type), pointer                     :: grid
  type(input_type)                             :: input
  type(option_type)                            :: option
  type(unstructured_grid_type), pointer        :: un_str_sfgrid
  character(len=MAXWORDLENGTH)                 :: word

  discretization => surf_realization%discretization

  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  do
    call InputReadFlotranString(input,option)
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SURFACE_FLOW')
    call StringToUpper(word)
    
    select case(trim(word))
      !.........................................................................
      ! Read surface grid information
      case ('SURF_GRID')
        call InputReadFlotranString(input,option)
        if (InputCheckExit(input,option)) exit

        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','SURF_GRID')
        call StringToUpper(word)
        select case(trim(word))
          case ('TYPE')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','TYPE')
            call StringToUpper(word)

            select case(trim(word))
              case ('UNSTRUCTURED')
                discretization%itype = UNSTRUCTURED_GRID
                call InputReadNChars(input,option, &
                                     discretization%filename, &
                                     MAXSTRINGLENGTH, &
                                     PETSC_TRUE)
                call InputErrorMsg(input,option,'keyword','filename')

                grid => GridCreate()
                un_str_sfgrid => UGridCreate()
                un_str_sfgrid%grid_type = TWO_DIM_GRID
                call UGridReadSurfGrid(un_str_sfgrid, &
                                       surf_realization%subsurf_filename, &
                                       discretization%filename, &
                                       option)
                grid%unstructured_grid => un_str_sfgrid
                discretization%grid => grid
                grid%itype = discretization%itype
                grid%ctype = discretization%ctype

              case default
              option%io_buffer = 'Surface-flow supports only unstructured grid'
              call printErrMsg(option)
            end select
          case default
            option%io_buffer = 'Keyword: ' // trim(word) // &
              ' not recognized in SURF_GRID '
            call printErrMsg(option)
        end select
        call InputSkipToEND(input,option,trim(word))

    end select
  enddo

end subroutine SurfaceFlowReadRequiredCardsFromInput

! ************************************************************************** !
!> This routine reads surface flow data from the input file
!! grids.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
subroutine SurfaceFlowRead(surf_realization,input,option)

  use Option_module
  use Input_module
  use String_module
  use Surface_Material_module
  use Surface_Realization_module
  use Grid_module
  use Structured_Grid_module
  use Unstructured_Grid_module
  use Discretization_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Strata_module

  implicit none

  type(surface_realization_type)               :: surf_realization
  type(discretization_type),pointer            :: discretization
  type(grid_type), pointer                     :: grid
  type(input_type)                             :: input
  type(option_type)                            :: option
  type(unstructured_grid_type), pointer        :: un_str_sfgrid
  type(surface_material_property_type),pointer :: surf_material_property
  type(region_type), pointer                   :: region
  type(flow_condition_type), pointer           :: flow_condition
  type(coupler_type), pointer                  :: coupler
  type(strata_type), pointer                   :: strata
  character(len=MAXWORDLENGTH)                 :: word

  discretization => surf_realization%discretization

  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  do
    call InputReadFlotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SURFACE_FLOW')
    call StringToUpper(word)

    select case(trim(word))
      !.........................................................................
      ! Read surface grid information
      case ('SURF_GRID')
        call InputSkipToEND(input,option,trim(word))

      !.........................................................................
      ! Read surface material information
      case ('SURF_MATERIAL_PROPERTY')
        surf_material_property => SurfaceMaterialPropertyCreate()

        call InputReadWord(input,option,surf_material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')
        call SurfaceMaterialPropertyRead(surf_material_property,input,option)
        call SurfaceMaterialPropertyAddToList(surf_material_property, &
                                          surf_realization%surf_material_properties)
        nullify(surf_material_property)

      !.........................................................................
      case ('SURF_REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','SURF_REGION')
        call printMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,surf_realization%surf_regions)
        nullify(region)

      !.........................................................................
      case ('SURF_FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'SURF_FLOW_CONDITION','name')
        call printMsg(option,flow_condition%name)
        if (option%iflowmode == G_MODE) then
          call FlowConditionGeneralRead(flow_condition,input,option)
        else
          call FlowConditionRead(flow_condition,input,option)
        endif
        call FlowConditionAddToList(flow_condition,surf_realization%surf_flow_conditions)
        nullify(flow_condition)

      !.........................................................................
      case ('SURF_BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Boundary Condition name')
        call CouplerRead(coupler,input,option)
        call SurfaceRealizationAddCoupler(surf_realization,coupler)
        nullify(coupler)

      !.........................................................................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call SurfaceRealizationAddStrata(surf_realization,strata)
        nullify(strata)
        
      !.........................................................................
      case ('SURF_INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name') 
        call CouplerRead(coupler,input,option)
        call SurfaceRealizationAddCoupler(surf_realization,coupler)
        nullify(coupler)        

      case default
        option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
                           'not recognized'
        call printErrMsg(option)

    end select
  enddo

end subroutine SurfaceFlowRead


! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowResidual(snes,xx,r,surf_realization,ierr)

  use Surface_Realization_module
  use Surface_Field_module
  use Patch_module
  use Level_module
  use Discretization_module
  use Option_module
  use Logging_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(surface_realization_type) :: surf_realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(surface_field_type), pointer :: surf_field
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(option_type), pointer :: option
  
  call PetscLogEventBegin(logging%event_r_residual,ierr)
  
  surf_field => surf_realization%surf_field
  discretization => surf_realization%discretization
  option => surf_realization%option

  ! pass #1 for internal and boundary flux terms
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceFlowResidualPatch1(snes,xx,r,surf_realization,ierr)
      !call SurfaceFlowResidualPatch1(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfaceFlowResidual

! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowResidualPatch1(snes,xx,r,surf_realization,ierr)

  use water_eos_module

  use Connection_module
  use Surface_Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  
  implicit none

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(surface_realization_type) :: surf_realization

  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:), mannings_loc_p(:),xx_loc_p(:)

  PetscReal, pointer :: face_fluxes_p(:)
  PetscInt :: icap_up, icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity, slope
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  PetscViewer :: viewer
  PetscReal :: xm_nacl, rho, dw_kg
  
  
  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  call GridVecGetArrayF90(grid,r, r_p, ierr)
  call GridVecGetArrayF90(grid, surf_field%mannings_loc, mannings_loc_p, ierr)
  call GridVecGetArrayF90(grid, surf_field%flow_xx_loc, xx_loc_p, ierr)

  r_p = 0.d0
  
  ! Interior Flux Terms -----------------------------------
  !write(*,*),'Interior fluxes:'
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
      
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &                  ! distance_gravity = dx*g*n
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)

      slope = cur_connection_set%dist(3,iconn)/ &
                dot_product(cur_connection_set%dist(1:3,iconn), &
                            cur_connection_set%dist(1:3,iconn))

      xm_nacl = option%m_nacl * FMWNACL
      xm_nacl = xm_nacl /(1.d3 + xm_nacl)

      call nacl_den(option%reference_temperature,option%reference_pressure*1.d-6,xm_nacl,dw_kg) 
      rho = dw_kg * 1.d3

      write(*,'(3I5,4es10.2)'),sum_connection,ghosted_id_dn,ghosted_id_up, &
        !fraction_upwind,distance,distance_gravity,dd_up,dd_dn,cur_connection_set%dist(1:3,iconn)
        slope,mannings_loc_p(ghosted_id_up),&
        (xx_loc_p(ghosted_id_up)-option%reference_pressure)/option%gravity(3)/dw_kg/1.d3, dw_kg*1.d3
        
      !call SurfaceFluxKinematic()

      Res(1) = 1.0d0

      if (.not.option%use_samr) then
         
        if (local_id_up>0) then
          r_p(local_id_up) = r_p(local_id_up) + Res(1)
        endif
         
        if (local_id_dn>0) then
          r_p(local_id_dn) = r_p(local_id_dn) - Res(1)
        endif
      endif
      
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  option%io_buffer = 'stopping for debugging'
  call printErrMsgByRank(option)

  ! Boundary Flux Terms -----------------------------------
  !write(*,*),'Boundary fluxes:'
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      !write(*,'(3I5)'),sum_connection,ghosted_id,local_id
  
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call GridVecRestoreArrayF90(grid,r, r_p, ierr)
  call GridVecRestoreArrayF90(grid, surf_field%mannings_loc, mannings_loc_p, ierr)
  call GridVecRestoreArrayF90(grid, surf_field%flow_xx_loc, xx_loc_p, ierr)

  call PetscViewerASCIIOpen(option%mycomm,'r_surfflow.out',viewer,ierr)
  call VecView(r,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

end subroutine SurfaceFlowResidualPatch1  

! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine SurfaceFlowJacobian(snes,xx,A,B,flag,surf_realization,ierr)

  use Surface_Realization_module
  use Level_module
  use Patch_module
  use Grid_module
  use Option_module
  use Logging_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(surface_realization_type) :: surf_realization
  MatStructure flag
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm


  call PetscLogEventBegin(logging%event_r_jacobian,ierr)

  option => surf_realization%option

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
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch

      call SurfaceFlowJacobianPatch1(snes,xx,J,J,flag,surf_realization,ierr)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfaceFlowJacobian


! ************************************************************************** !
!
!
! ************************************************************************** !

subroutine SurfaceFlowJacobianPatch1(snes,xx,A,B,flag,surf_realization,ierr)
       
  use water_eos_module

  use Connection_module
  use Surface_Realization_module
  use Option_module
  use Patch_module
  use Grid_module
  use Coupler_module
  use Field_module
  use Debug_module
    
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(surface_realization_type) :: surf_realization
  MatStructure flag

  PetscErrorCode :: ierr

  PetscReal, pointer :: porosity_loc_p(:), &
                        perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscInt :: icap_up,icap_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: upweight
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal :: Jup(surf_realization%option%nflowdof,surf_realization%option%nflowdof), &
               Jdn(surf_realization%option%nflowdof,surf_realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  !type(field_type), pointer :: field 
  !type(richards_parameter_type), pointer :: richards_parameter
  !type(richards_auxvar_type), pointer :: rich_aux_vars(:), rich_aux_vars_bc(:) 
  !type(global_auxvar_type), pointer :: global_aux_vars(:), global_aux_vars_bc(:) 
  
  PetscViewer :: viewer

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  !field => realization%field
  !richards_parameter => patch%aux%Richards%richards_parameter
  !rich_aux_vars => patch%aux%Richards%aux_vars
  !rich_aux_vars_bc => patch%aux%Richards%aux_vars_bc
  !global_aux_vars => patch%aux%Global%aux_vars
  !global_aux_vars_bc => patch%aux%Global%aux_vars_bc

  !call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)
  !call GridVecGetArrayF90(grid,field%perm_xx_loc, perm_xx_loc_p, ierr)
  !call GridVecGetArrayF90(grid,field%perm_yy_loc, perm_yy_loc_p, ierr)
  !call GridVecGetArrayF90(grid,field%perm_zz_loc, perm_zz_loc_p, ierr)

  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_set)) exit
    
    !write(*,*),'cur_connection_set%num_connections : ',cur_connection_set%num_connections
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      
      !write(*,*),'ghosted_ids :',ghosted_id_up, ghosted_id_dn

      !if (patch%imat(ghosted_id_up) <= 0 .or. &
      !    patch%imat(ghosted_id_dn) <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      !write(*,*),'local_ids   :',ghosted_id_up, ghosted_id_dn,local_id_up,local_id_dn

      Jup = 1.0d0
      Jdn = 1.0d0

      if (local_id_up > 0) then
          call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                        Jup,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                        Jdn,ADD_VALUES,ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                              Jdn,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                               Jup,ADD_VALUES,ierr)
      endif
      
    enddo
    cur_connection_set => cur_connection_set%next

  enddo
   
  !realization%option%io_buffer = 'stopping '
  !call printErrMsg(realization%option)
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  call PetscViewerASCIIOpen(option%mycomm,'jacobian_surfflow.out',viewer,ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

end subroutine SurfaceFlowJacobianPatch1

end module Surface_Flow_module

#endif
