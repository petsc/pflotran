module Richards_LSM_module
#include "finclude/petscmat.h"
  use petscmat
  use Richards_Aux_module
  use Richards_Common_module
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  ! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: RichardsUpdateLSMAuxVarsPatch,&
            RichardsLSMFluxDerivative, &
            RichardsLSMFlux
  
contains

! ************************************************************************** !

subroutine RichardsUpdateLSMAuxVarsPatch(realization)
  ! 
  ! This routine computes the gradient of pressure using a least-square-method
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 11/24/12
  ! 

  use Realization_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Logging_module

  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(richards_auxvar_type), pointer :: rich_auxvars(:) 
  type(richards_auxvar_type), pointer :: rich_auxvars_bc(:)
  type(richards_auxvar_type), pointer :: rich_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)  
  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase, i
  PetscReal, pointer :: xx_loc_p(:), xx_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  PetscInt :: max_stencil_width
  PetscInt :: nid
  Vec :: phi,A,B
  PetscReal, pointer :: phi_p(:), b_p(:)
  PetscInt, pointer :: cell_neighbors(:,:)
  PetscReal :: distance_gravity,distance_gravity0, den_avg
  PetscReal :: dx,dy,dz
  PetscReal :: P,Pn

  max_stencil_width = 2

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  rich_auxvars => patch%aux%Richards%auxvars
  rich_auxvars_bc => patch%aux%Richards%auxvars_bc
  rich_auxvars_ss => patch%aux%Richards%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
    
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  select case(grid%itype)
    case(STRUCTURED_GRID)
      cell_neighbors => grid%structured_grid%cell_neighbors
    case(UNSTRUCTURED_GRID)
      option%io_buffer='GridComputeMinv() not implemented for unstructured grid.'
      call printErrMsg(option)
  end select

  do ghosted_id=1,grid%ngmax
    if(grid%ghosted_level(ghosted_id)<max_stencil_width) then
      call VecCreateSeq(PETSC_COMM_SELF,cell_neighbors(0,ghosted_id),phi, &
                        ierr);CHKERRQ(ierr)
      call VecGetArrayF90(phi,phi_p,ierr);CHKERRQ(ierr)
      do nid=1,cell_neighbors(0,ghosted_id)
        distance_gravity = option%gravity(1)*grid%x(cell_neighbors(nid,ghosted_id)) + &
                           option%gravity(2)*grid%y(cell_neighbors(nid,ghosted_id)) + &
                           option%gravity(3)*grid%z(cell_neighbors(nid,ghosted_id))

        Pn = xx_loc_p(cell_neighbors(nid,ghosted_id)) - &
             global_auxvars(cell_neighbors(nid,ghosted_id))%den(1)*FMWH2O*distance_gravity

        distance_gravity0 = option%gravity(1)*grid%x(ghosted_id) + &
                            option%gravity(2)*grid%y(ghosted_id) + &
                            option%gravity(3)*grid%z(ghosted_id)
        P = xx_loc_p(ghosted_id) - &
             global_auxvars(ghosted_id)%den(1)*FMWH2O*distance_gravity0

        phi_p(nid) = Pn-P
      enddo

      do nid=1,cell_neighbors(0,ghosted_id)
        if(abs(phi_p(nid))<1.D-20) phi_p(nid) = 0.d0
      enddo

      call VecRestoreArrayF90(phi,phi_p,ierr);CHKERRQ(ierr)

      call VecCreateSeq(PETSC_COMM_SELF,3,A,ierr);CHKERRQ(ierr)
      call MatMult(grid%dispT(ghosted_id),phi,A,ierr);CHKERRQ(ierr)

      call VecCreateSeq(PETSC_COMM_SELF,3,B,ierr);CHKERRQ(ierr)
      call MatMult(grid%Minv(ghosted_id),A,B,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(B,b_p,ierr);CHKERRQ(ierr)

      ! If the gradient is too small, make it zero
      if(abs(b_p(1))<1.D-20) b_p(1) = 0.d0
      if(abs(b_p(2))<1.D-20) b_p(2) = 0.d0
      if(abs(b_p(3))<1.D-20) b_p(3) = 0.d0

      global_auxvars(ghosted_id)%dphi(1,:) = b_p(:)

      call VecRestoreArrayF90(B,b_p,ierr);CHKERRQ(ierr)

      call VecDestroy(phi,ierr);CHKERRQ(ierr)
      call VecDestroy(A,ierr);CHKERRQ(ierr)
      call VecDestroy(B,ierr);CHKERRQ(ierr)
    endif
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

end subroutine RichardsUpdateLSMAuxVarsPatch

! ************************************************************************** !

subroutine RichardsLSMFluxDerivative(rich_auxvar_up,global_auxvar_up,por_up, &
                                  sir_up,dd_up,perm_up, &
                                  rich_auxvar_dn,global_auxvar_dn,por_dn, &
                                  sir_dn,dd_dn,perm_dn, &
                                  area, dist, dist_gravity,upweight, &
                                  distance, &
                                  option,sat_func_up,sat_func_dn, &
                                  jacfac,ghosted_id_up,ghosted_id_dn, &
                                  cell_neighbors, &
                                  x,y,z, &
                                  bnd_cell, &
                                  Jup,Jdn)
  ! 
  ! This routine computes the derivatives of the internal flux terms for the
  ! Jacobian, when least-squares-method is used to compute the gradients.
  ! NOTE: Each internal flux term, depends on states of 'up', 'down', and
  ! neighbors of 'up' & 'down' cells. Thus, jacobian matrix should have
  ! entries in multiple columns of the Jacobian matrix that corresponds
  ! to neigboring cells for a given internal flux term. But, the entries
  ! corresponding to neigboring cells are neglected presently.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/17/12
  ! 
  use Option_module
  use Saturation_Function_module

  implicit none

  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy, area, dist(3)
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: distance
  type(saturation_function_type) :: sat_func_up, sat_func_dn
  PetscInt :: ghosted_id_up,ghosted_id_dn
  PetscReal,pointer :: jacfac(:,:,:)
  PetscInt, pointer :: cell_neighbors(:,:)
  PetscReal, pointer :: x(:),y(:),z(:)
  PetscBool, pointer :: bnd_cell(:)
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  PetscReal :: dphi_x,dphi_y,dphi_z
  PetscReal :: dphi_x_dp_up,dphi_y_dp_up,dphi_z_dp_up
  PetscReal :: dphi_x_dp_dn,dphi_y_dp_dn,dphi_z_dp_dn

  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn
  PetscReal :: dq_dp_up, dq_dp_dn
  PetscReal :: dHup_x_dp_up, dHup_y_dp_up, dHup_z_dp_up
  PetscReal :: dHdn_x_dp_dn, dHdn_y_dp_dn, dHdn_z_dp_dn

  PetscInt :: iphase, ideriv
  PetscInt :: nid_up, nid_dn,ii
  type(richards_auxvar_type) :: rich_auxvar_pert_up, rich_auxvar_pert_dn
  type(global_auxvar_type) :: global_auxvar_pert_up, global_auxvar_pert_dn
  PetscReal :: x_up(1), x_dn(1), x_pert_up(1), x_pert_dn(1), pert_up, pert_dn, &
            res(1), res_pert_up(1), res_pert_dn(1), J_pert_up(1,1), J_pert_dn(1,1)

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  v_darcy = 0.D0
  ukvr = 0.d0

  Jup = 0.d0
  Jdn = 0.d0

  dden_ave_dp_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0

  nid_up = -1
  nid_dn = -1

  ! One of the cells is a boundary cell, so use 2-point flux
  if(bnd_cell(ghosted_id_up).or.bnd_cell(ghosted_id_dn)) then
#if 0         
    call RichardsFluxDerivative(rich_auxvar_up,global_auxvar_up,por_up, &
                                  sir_up,dd_up,perm_up, &
                                  rich_auxvar_dn,global_auxvar_dn,por_dn, &
                                  sir_dn,dd_dn,perm_dn, &
                                  area, dist, dist_gravity,upweight, &
                                  option,sat_func_up,sat_func_dn,Jup,Jdn)
#endif    
    return
  endif

  do ii = 1,cell_neighbors(0,ghosted_id_dn)
    if (cell_neighbors(ii,ghosted_id_dn) == ghosted_id_up) nid_up = ii
  enddo
  do ii = 1,cell_neighbors(0,ghosted_id_up)
    if (cell_neighbors(ii,ghosted_id_up) == ghosted_id_dn) nid_dn = ii
  enddo
  if(nid_up<0 .or. nid_dn <0 ) then
    option%io_buffer = 'Neighbors not found'
    call printErrMsg(option)
  endif

  ! Flow term
  if (global_auxvar_up%sat(1) > sir_up .or. global_auxvar_dn%sat(1) > sir_dn) then
    if (global_auxvar_up%sat(1) <eps) then
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then
      upweight=1.d0
    endif

    density_ave =        upweight*global_auxvar_up%den(1)+ &
                  (1.D0-upweight)*global_auxvar_dn%den(1)
    dden_ave_dp_up = upweight*rich_auxvar_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*rich_auxvar_dn%dden_dp

    dphi_x =       upweight *global_auxvar_up%dphi(1,1) + &
             (1.d0-upweight)*global_auxvar_dn%dphi(1,1)
    dphi_y =       upweight *global_auxvar_up%dphi(1,2) + &
             (1.d0-upweight)*global_auxvar_dn%dphi(1,2)
    dphi_z =       upweight *global_auxvar_up%dphi(1,3) + &
             (1.d0-upweight)*global_auxvar_dn%dphi(1,3)

    dphi = -(dist(1)*dphi_x + dist(2)*dphi_y + dist(3)*dphi_z)*distance

    ! H_x = P - W * rho * g *z
    dHup_x_dp_up = 1.d0 - FMWH2O*rich_auxvar_up%dden_dp*option%gravity(1)*x(ghosted_id_up)
    dHup_y_dp_up = 1.d0 - FMWH2O*rich_auxvar_up%dden_dp*option%gravity(2)*y(ghosted_id_up)
    dHup_z_dp_up = 1.d0 - FMWH2O*rich_auxvar_up%dden_dp*option%gravity(3)*z(ghosted_id_up)

    dHdn_x_dp_dn = 1.d0 - FMWH2O*rich_auxvar_dn%dden_dp*option%gravity(1)*x(ghosted_id_dn)
    dHdn_y_dp_dn = 1.d0 - FMWH2O*rich_auxvar_dn%dden_dp*option%gravity(2)*y(ghosted_id_dn)
    dHdn_z_dp_dn = 1.d0 - FMWH2O*rich_auxvar_dn%dden_dp*option%gravity(3)*z(ghosted_id_dn)

    dphi_x_dp_up =       upweight *jacfac(ghosted_id_up,0     ,1)*dHup_x_dp_up + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,nid_up,1)*dHup_x_dp_up
    dphi_y_dp_up =       upweight *jacfac(ghosted_id_up,0     ,2)*dHup_y_dp_up + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,nid_up,2)*dHup_y_dp_up
    dphi_z_dp_up =       upweight *jacfac(ghosted_id_up,0     ,3)*dHup_z_dp_up + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,nid_up,3)*dHup_z_dp_up

    dphi_x_dp_dn =       upweight *jacfac(ghosted_id_up,nid_dn,1)*dHdn_x_dp_dn + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,0     ,1)*dHdn_x_dp_dn
    dphi_y_dp_dn =       upweight *jacfac(ghosted_id_up,nid_dn,2)*dHdn_y_dp_dn + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,0     ,2)*dHdn_y_dp_dn
    dphi_z_dp_dn =       upweight *jacfac(ghosted_id_up,nid_dn,3)*dHdn_z_dp_dn + &
                   (1.d0-upweight)*jacfac(ghosted_id_dn,0     ,3)*dHdn_z_dp_dn

    dphi_dp_up = -(dist(1)*dphi_x_dp_up + dist(2)*dphi_y_dp_up + dist(3)*dphi_z_dp_up)
    dphi_dp_dn = -(dist(1)*dphi_x_dp_dn + dist(2)*dphi_y_dp_dn + dist(3)*dphi_z_dp_dn)

    if (dphi>=0.D0) then
      ukvr = rich_auxvar_up%kvr
      dukvr_dp_up = rich_auxvar_up%dkvr_dp
    else
      ukvr = rich_auxvar_dn%kvr
      dukvr_dp_dn = rich_auxvar_dn%dkvr_dp
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

end subroutine RichardsLSMFluxDerivative

! ************************************************************************** !

subroutine RichardsLSMFlux(rich_auxvar_up,global_auxvar_up, &
                           por_up,sir_up,dd_up,perm_up, &
                           rich_auxvar_dn,global_auxvar_dn, &
                           por_dn,sir_dn,dd_dn,perm_dn, &
                           area, dist, dist_gravity,upweight, distance, &
                           option, &
                           ghosted_id_up, ghosted_id_dn, cell_neighbors, &
                           bnd_cell, &
                           v_darcy,Res)
  ! 
  ! This routine computes the internal flux terms for the residual term using
  ! a least-square-method for gradient computation.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/10/12
  ! 
  use Option_module

  implicit none

  type(richards_auxvar_type) :: rich_auxvar_up, rich_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: v_darcy,area, dist(3)
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscInt  :: ghosted_id_up,ghosted_id_dn
  PetscInt, pointer :: cell_neighbors(:,:)
  PetscBool, pointer :: bnd_cell(:)

  PetscInt :: ispec
  PetscReal :: fluxm, q
  PetscReal :: ukvr,Dq
  PetscReal :: upweight,density_ave,cond,dphi
  PetscReal :: dphi_x,dphi_y,dphi_z
  PetscReal :: fraction_upwind,distance

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  fluxm = 0.d0
  v_darcy = 0.D0
  ukvr = 0.d0

 ! One of the cells is a boundary cell, so use 2-point flux
  if(bnd_cell(ghosted_id_up).or.bnd_cell(ghosted_id_dn)) then
#if 0         
    call RichardsFlux(rich_auxvar_up,global_auxvar_up, &
                      por_up,sir_up,dd_up,perm_up, &
                      rich_auxvar_dn,global_auxvar_dn, &
                      por_dn,sir_dn,dd_dn,perm_dn, &
                      area, dist, dist_gravity,upweight, &
                      option,v_darcy,Res)
#endif    
    return
  endif

  ! Flow term
  if (global_auxvar_up%sat(1) > sir_up .or. global_auxvar_dn%sat(1) > sir_dn) then
    if (global_auxvar_up%sat(1) <eps) then
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then
      upweight=1.d0
    endif

    density_ave = upweight*global_auxvar_up%den(1)+ &
                  (1.D0-upweight)*global_auxvar_dn%den(1)

    dphi_x =       upweight *global_auxvar_up%dphi(1,1) + &
             (1.d0-upweight)*global_auxvar_dn%dphi(1,1)
    dphi_y =       upweight *global_auxvar_up%dphi(1,2) + &
             (1.d0-upweight)*global_auxvar_dn%dphi(1,2)
    dphi_z =       upweight *global_auxvar_up%dphi(1,3) + &
             (1.d0-upweight)*global_auxvar_dn%dphi(1,3)

    dphi = -(dist(1)*dphi_x + dist(2)*dphi_y + dist(3)*dphi_z)*distance

    if (dphi>=0.D0) then
      ukvr = rich_auxvar_up%kvr
    else
      ukvr = rich_auxvar_dn%kvr
    endif

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
      q = v_darcy * area
      fluxm = q*density_ave
    endif
  endif

  Res(1) = fluxm
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL

end subroutine RichardsLSMFlux

end module Richards_LSM_module
