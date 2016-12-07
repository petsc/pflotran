module MFD_module
#include "finclude/petscmat.h"
  use petscmat

  use Connection_module
  use Grid_module
  use MFD_Aux_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  public :: MFDCreateJacobian, &
            MFDInitializeMassMatrices, MFDAuxGenerateStiffMatrix,&
            MFDAuxGenerateRhs, MFDAuxReconstruct, MFDAuxFluxes, MFDComputeDensity,&
            MFDAuxJacobianLocal, MFDAuxUpdateCellPressure, MFDCreateJacobianLP, &
            MFDAuxGenerateRhs_LP, MFDAuxJacobianLocal_LP

contains

! ************************************************************************** !

subroutine MFDCreateJacobian(grid, mfd_aux, mat_type, J, option)
  ! 
  ! Creates a Jacobian matrix for  the faced based dof
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 05/14/10
  ! 
#include "finclude/petscmat.h"
 use petscmat
 use Option_module
 use Grid_module
 use MFD_Aux_module

  implicit none

  type(grid_type) :: grid
  type(mfd_type) :: mfd_aux
  type(option_type) :: option
  MatType :: mat_type
  Mat :: J


  type(mfd_auxvar_type), pointer :: auxvar
  PetscInt, allocatable :: d_nnz(:), o_nnz(:)
  PetscInt :: icell, iface
  PetscInt :: ghost_face_id, local_face_id, ghost_face_id_n, local_face_id_n
  PetscInt :: icount, jcount
  PetscInt :: ndof_local
  PetscErrorCode :: ierr

  allocate(d_nnz(grid%nlmax_faces))
  allocate(o_nnz(grid%nlmax_faces))

  d_nnz = 1 ! start 1 since diagonal connection to self
  o_nnz = 0

  do icell = 1, grid%nlmax
    auxvar => MFD_aux%auxvars(icell)
    do icount = 1, auxvar%numfaces
      ghost_face_id = auxvar%face_id_gh(icount)
      local_face_id = grid%fG2L(ghost_face_id)
      if (local_face_id > 0) then
        do jcount = 1, auxvar%numfaces
          if (icount == jcount) cycle
          ghost_face_id_n = auxvar%face_id_gh(jcount)
          local_face_id_n = grid%fG2L(ghost_face_id_n)
          if (local_face_id_n > 0) then
             d_nnz(local_face_id) = d_nnz(local_face_id) + 1
          else
             o_nnz(local_face_id) = o_nnz(local_face_id) + 1
          endif
        enddo
      endif
    enddo
  enddo

  ndof_local = mfd_aux%ndof * grid%nlmax_faces

  if (option%mycommsize > 1) then
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*mfd_aux%ndof
        o_nnz = o_nnz*mfd_aux%ndof
        call MatCreateAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr);CHKERRQ(ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces, mfd_aux%mapping_ltog_faces,  &
                                        ierr);CHKERRQ(ierr)

      case(MATBAIJ)
        call MatCreateBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr);CHKERRQ(ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,mfd_aux%mapping_ltog_faces, &
                                        ierr);CHKERRQ(ierr)
        
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobian'
        call printErrMsg(option)
    end select
  else
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*mfd_aux%ndof
        call MatCreateSeqAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr);CHKERRQ(ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces, mfd_aux%mapping_ltog_faces,  &
                                        ierr);CHKERRQ(ierr)
      case(MATBAIJ)
        call MatCreateSeqBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr);CHKERRQ(ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,mfd_aux%mapping_ltog_faces, &
                                        ierr);CHKERRQ(ierr)
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobian'
        call printErrMsg(option)
    end select
  endif

  deallocate(d_nnz)
  deallocate(o_nnz)

end subroutine MFDCreateJacobian

! ************************************************************************** !

subroutine MFDCreateJacobianLP(grid, mfd_aux, mat_type, J, option)
#include "finclude/petscmat.h"
  use petscmat
  use Option_module
  use Grid_module
  use MFD_Aux_module
  use Connection_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module

  implicit none

  type(grid_type) :: grid
  type(mfd_type) :: mfd_aux
  type(option_type) :: option
  MatType :: mat_type
  Mat :: J


  type(mfd_auxvar_type), pointer :: auxvar
  PetscInt, allocatable :: d_nnz(:), o_nnz(:)
  type(connection_set_type), pointer :: conn
  PetscInt :: icell, iface, i
  PetscInt :: ghost_face_id, local_face_id, ghost_face_id_n, local_face_id_n
  PetscInt :: icount, jcount, jface, loc_id_up, loc_id_dn
  PetscInt :: ndof_local
  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  type(unstructured_grid_type),pointer :: ugrid

  allocate(d_nnz(grid%nlmax_faces + grid%nlmax))
  allocate(o_nnz(grid%nlmax_faces + grid%nlmax))

  ugrid => grid%unstructured_grid

  d_nnz = 1 ! start 1 since diagonal connection to self
  o_nnz = 0

  if (grid%itype == STRUCTURED_GRID_MIMETIC) then

    !
    ! Structured grid
    !

    do icell = 1, grid%nlmax
      auxvar => MFD_aux%auxvars(icell)
      do icount = 1, auxvar%numfaces
        ghost_face_id = auxvar%face_id_gh(icount)
        local_face_id = grid%fG2L(ghost_face_id)
        if (local_face_id > 0) then
          do jcount = 1, auxvar%numfaces
            if (icount == jcount) cycle
            ghost_face_id_n = auxvar%face_id_gh(jcount)
            local_face_id_n = grid%fG2L(ghost_face_id_n)
            if (local_face_id_n > 0) then
              d_nnz(local_face_id) = d_nnz(local_face_id) + 1
            else
              o_nnz(local_face_id) = o_nnz(local_face_id) + 1
            endif
          enddo
          d_nnz(local_face_id) = d_nnz(local_face_id) + 1   ! connection to cell-centered dof
        endif

        conn => grid%faces(ghost_face_id)%conn_set_ptr
        jface = grid%faces(ghost_face_id)%id
        if (conn%itype == INTERNAL_CONNECTION_TYPE) then
          loc_id_dn = grid%nG2L(conn%id_dn(jface))
          loc_id_up = grid%nG2L(conn%id_up(jface))
          if ((loc_id_dn > 0).and.(loc_id_up > 0)) then
            d_nnz(grid%nlmax_faces + icell) = d_nnz(grid%nlmax_faces + icell) + 1
          else
            o_nnz(grid%nlmax_faces + icell) = o_nnz(grid%nlmax_faces + icell) + 1
            if (local_face_id > 0) then
              o_nnz(local_face_id) = o_nnz(local_face_id) + 6
            endif
          endif
        else
          d_nnz(grid%nlmax_faces + icell) = d_nnz(grid%nlmax_faces + icell) + 1
        endif

        if (local_face_id > 0) then
          d_nnz(grid%nlmax_faces + icell) = d_nnz(grid%nlmax_faces + icell) + 1
        else
          o_nnz(grid%nlmax_faces + icell) = o_nnz(grid%nlmax_faces + icell) + 1
        endif
      enddo
    enddo

  else
    !
    ! Unstructured grid
    !

    ! 1) For Pressre at cell-faces
    do ghost_face_id = 1,grid%ngmax_faces

      local_face_id = grid%fG2L(ghost_face_id)
      conn => grid%faces(ghost_face_id)%conn_set_ptr
      jface = grid%faces(ghost_face_id)%id

      if (conn%itype == INTERNAL_CONNECTION_TYPE) then

        ! Interal connection

        loc_id_dn = grid%nG2L(conn%id_dn(jface))
        loc_id_up = grid%nG2L(conn%id_up(jface))

        d_nnz(local_face_id) = 1

        ! Check if downstream cell is local
        if(loc_id_dn>0) then
          d_nnz(local_face_id) = d_nnz(local_face_id) + &
                      UCellGetNFaces(ugrid%cell_type(conn%id_dn(jface)),option)
        else
          if(local_face_id>0) then
          o_nnz(local_face_id) = o_nnz(local_face_id) + &
                      UCellGetNFaces(ugrid%cell_type(conn%id_dn(jface)),option)
          endif
          o_nnz(grid%nlmax_faces+loc_id_up)=o_nnz(grid%nlmax_faces+loc_id_up)+1
        endif

        ! Check if upstream cell is local
        if(loc_id_up>0) then
          d_nnz(local_face_id) = d_nnz(local_face_id) + &
                      UCellGetNFaces(ugrid%cell_type(conn%id_up(jface)),option)
        else
          if(local_face_id>0) then
          o_nnz(local_face_id) = o_nnz(local_face_id) + &
                      UCellGetNFaces(ugrid%cell_type(conn%id_up(jface)),option)
          endif
          o_nnz(grid%nlmax_faces+loc_id_dn)=o_nnz(grid%nlmax_faces+loc_id_dn)+1
        endif

        ! GB: This is a temporary fix
        o_nnz(local_face_id) = 6*2
        d_nnz(local_face_id) = 13

      else

        loc_id_dn = grid%nG2L(conn%id_dn(jface))

        ! External connection
        d_nnz(local_face_id) = &
                  UCellGetNFaces(ugrid%cell_type(conn%id_dn(jface)),option) + 1
        o_nnz(local_face_id) = 6
        o_nnz(grid%nlmax_faces+loc_id_dn)=o_nnz(grid%nlmax_faces+loc_id_dn)+1

        ! GB: This is a temporary fix
        o_nnz(local_face_id) = 6*2
        d_nnz(local_face_id) = 13
      endif
    enddo

    ! 1) For Pressre at cell-centers
    do local_id = 1,grid%nlmax

      ghosted_id = grid%nL2G(local_id)
      d_nnz(grid%nlmax_faces+local_id) = &
                          UCellGetNFaces(ugrid%cell_type(ghosted_id),option) + 1

      do icount = 1,ugrid%cell_neighbors_local_ghosted(0,local_id)
        if(ugrid%cell_neighbors_local_ghosted(icount,local_id)>0) then
          d_nnz(grid%nlmax_faces+local_id) = d_nnz(grid%nlmax_faces+local_id) + 1
        else
          o_nnz(grid%nlmax_faces+local_id) = o_nnz(grid%nlmax_faces+local_id) + 1
        endif
      enddo

        ! GB: This is a temporary fix
        o_nnz(grid%nlmax_faces+local_id) = 6*2
        d_nnz(grid%nlmax_faces+local_id) = 13
    enddo
  endif


  ndof_local = mfd_aux%ndof * (grid%nlmax_faces + grid%nlmax)

  if (option%mycommsize > 1) then
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*mfd_aux%ndof
        o_nnz = o_nnz*mfd_aux%ndof
        call MatCreateAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr);CHKERRQ(ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_LP, mfd_aux%mapping_ltog_LP,  &
                                        ierr);CHKERRQ(ierr)

      case(MATBAIJ)
        call MatCreateBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr);CHKERRQ(ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_LP,mfd_aux%mapping_ltog_LP, &
                                        ierr);CHKERRQ(ierr)
        
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobianLP'
        call printErrMsg(option)
    end select
  else
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*mfd_aux%ndof
        call MatCreateSeqAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr);CHKERRQ(ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_LP, mfd_aux%mapping_ltog_LP,  &
                                        ierr);CHKERRQ(ierr)
      case(MATBAIJ)
        call MatCreateSeqBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr);CHKERRQ(ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_LP, mfd_aux%mapping_ltog_LP, &
                                        ierr);CHKERRQ(ierr)
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobianLP'
        call printErrMsg(option)
    end select
  endif

  deallocate(d_nnz)
  deallocate(o_nnz)


end subroutine MFDCreateJacobianLP

! ************************************************************************** !

subroutine MFDInitializeMassMatrices(grid, field, &
                                     mfd_aux, material_auxvars, option)
#include "finclude/petscmat.h"
 use petscmat
 use Option_module
 use Grid_module
 use MFD_Aux_module
 use Material_Aux_class
 use Field_module

  implicit none

  type(grid_type) :: grid
  type(field_type) :: field
  type(mfd_type) :: mfd_aux
  type(option_type) :: option
  class(material_auxvar_type), pointer :: material_auxvars(:)

  type(mfd_auxvar_type), pointer :: auxvar
  PetscInt :: ghosted_cell_id, icell, i, j
  PetscErrorCode :: ierr
  PetscReal :: PermTensor(3,3) 


  do icell = 1, grid%nlmax

    ghosted_cell_id = grid%nL2G(icell)

    PermTensor = 0.
    PermTensor(1,1) = material_auxvars(ghosted_cell_id)%permeability(perm_xx_index)
    PermTensor(2,2) = material_auxvars(ghosted_cell_id)%permeability(perm_yy_index)
    PermTensor(3,3) = material_auxvars(ghosted_cell_id)%permeability(perm_zz_index)
    PermTensor(1,3) = material_auxvars(ghosted_cell_id)%permeability(perm_xz_index)
    PermTensor(1,2) = material_auxvars(ghosted_cell_id)%permeability(perm_xy_index)
    PermTensor(2,3) = material_auxvars(ghosted_cell_id)%permeability(perm_yz_index)
    PermTensor(2,1) = PermTensor(1,2)
    PermTensor(3,1) = PermTensor(1,3)
    PermTensor(3,2) = PermTensor(2,3)

    auxvar => mfd_aux%auxvars(icell)
    call MFDAuxGenerateMassMatrixInv(grid, ghosted_cell_id, auxvar, &
                                     material_auxvars(ghosted_cell_id)%volume, &
                                     PermTensor, option)
    call MFDAuxInitResidDerivArrays(auxvar, option)
    call MFDAuxComputeGeometricValues (grid, ghosted_cell_id, auxvar, PermTensor, option)
  enddo

end subroutine MFDInitializeMassMatrices

! ************************************************************************** !

subroutine MFDAuxGenerateStiffMatrix(auxvar, rich_auxvar, global_auxvar,  &
                                       sq_faces, option)

  use Option_module
  use Richards_Aux_module
  use Global_Aux_module

  implicit none

  type(mfd_auxvar_type), pointer :: auxvar
  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscScalar, pointer :: sq_faces(:)
  type(option_type) :: option
  PetscScalar, pointer :: StiffMatrix(:,:)

  PetscInt :: iface, jface, i,j
  PetscScalar :: E
  PetscScalar, pointer :: MB(:)

  allocate(MB(auxvar%numfaces))

  E = 0
  do iface = 1, auxvar%numfaces
    MB(iface) = 0.
    do jface = 1, auxvar%numfaces
       MB(iface) = MB(iface) + auxvar%MassMatrixInv(iface,jface)*sq_faces(jface)
    enddo
    E = E + MB(iface)*sq_faces(iface)
  enddo

  do iface = 1, auxvar%numfaces
    do jface = 1, auxvar%numfaces
        auxvar%StiffMatrix(iface,jface) = sq_faces(iface)*sq_faces(jface)*    &
                                        (auxvar%MassMatrixInv(iface,jface) - &
                                        (1./E)*MB(iface)*MB(jface))
    enddo
  enddo


  deallocate(MB)

end subroutine MFDAuxGenerateStiffMatrix

! ************************************************************************** !

subroutine MFDAuxGenerateRhs(patch, grid, ghosted_cell_id, PermTensor, bc_g, source_f,  bc_h, auxvar, &
                                       rich_auxvar, global_auxvar, Accum, &
                                       porosity, volume, pres, face_pres, bnd,&
                                       sq_faces, neig_den, neig_kvr, neig_dkvr_dp, option, rhs)

  use Option_module
  use Richards_Aux_module
  use Global_Aux_module
  use Patch_module

  implicit none

  type(patch_type) :: patch
  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: auxvar
  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscScalar :: sq_faces(:), neig_den(:), neig_kvr(:), neig_dkvr_dp(:)
  type(option_type) :: option
  PetscScalar :: bc_g(:), rhs(:), bc_h(:), face_pres(:), bnd(:)
  PetscScalar :: Accum(1:option%nflowdof),source_f(1:option%nflowdof), pres(1:option%nflowdof)
  PetscScalar :: PermTensor(3,3), sat
  PetscInt :: ghosted_cell_id, bound_id

  PetscScalar :: Kg(3), dir_norm(3)
  PetscInt :: i, j, ghost_face_id, local_face_id
  type(connection_set_type), pointer :: conn
  PetscInt, parameter :: numfaces = 6

  PetscInt :: iface, jface
  PetscScalar :: E, den_gWB, den_BWB, den_BWCl , div_grav_term
  PetscScalar :: ukvr(6), dukvr_dp(6), ds_dp, porosity, volume, den_cntr, dden_cntr_dp
  PetscScalar :: PorVol_dt
  PetscScalar :: norm_len
  PetscScalar :: Wg(6), WClm(6), CWCl(6), f(6), dden_dp(6), dir, v_darcy(6)
  PetscScalar :: dbeta_dp(6), den(6), beta(6), denB(6), Clm(6), Rgr(6), BAR(6), BdARdp(6)
  PetscScalar :: BARWB, BARWClm, BdARdpWB, BdARdpWClm
  PetscScalar :: upweight, delta_p

  ukvr = 0. !rich_auxvar%kvr_x
  dukvr_dp = 0.!rich_auxvar%dkvr_x_dp
  den_cntr = global_auxvar%den(1)
  dden_cntr_dp =  rich_auxvar%dden_dp

  upweight = 0.5

  do i = 1 , numfaces
    if (bnd(i)==1) then 
        call MFDComputeDensity(global_auxvar, face_pres(i) + bc_g(i)/sq_faces(i), den(i), dden_dp(i), option)
    else  
        den(i) =  upweight*den_cntr + (1.D0-upweight)*neig_den(i)
        dden_dp(i) = dden_cntr_dp

    endif

#ifdef USE_ANISOTROPIC_MOBILITY
        if (rich_auxvar%kvr_x < 1e-10) then
            ukvr(i) = 1e-10
        else 
            ukvr(i) = rich_auxvar%kvr_x
        endif
        dukvr_dp(i) = rich_auxvar%dkvr_x_dp 
#else
        ukvr(i) = rich_auxvar%kvr
        dukvr_dp(i) = rich_auxvar%dkvr_dp 
#endif

     beta(i) = ukvr(i)*den(i)
     dbeta_dp(i) = dukvr_dp(i)*den(i) + ukvr(i)*dden_dp(i)

     denB(i) = sq_faces(i)*den(i)
     BAR(i) = sq_faces(i)*ukvr(i)*den(i)
     BdARdp(i) = sq_faces(i)*dbeta_dp(i)
  enddo

  sat = global_auxvar%sat(1)
  ds_dp = rich_auxvar%dsat_dp
  PorVol_dt = porosity*volume/option%flow_dt

  E = 0.
  f(1) = Accum(1) - source_f(1)
  Wg = 0.

  do iface = 1, numfaces
     Clm(iface) = sq_faces(iface)*face_pres(iface) + bc_g(iface)
     Rgr(iface) = den(iface)*auxvar%gr(iface)
  enddo 

  Wg = matmul(auxvar%MassMatrixInv, bc_g)
  WClm = matmul(auxvar%MassMatrixInv, Clm)

  BARWB = dot_product(auxvar%WB, BAR)
  BARWClm = dot_product(WClm, BAR)
  BdARdpWB = dot_product(auxvar%WB, BdARdp)
  BdARdpWClm = dot_product(WClm, BdARdp)
  div_grav_term = dot_product(Rgr, BAR)

   auxvar%Rp = BARWB*pres(1) - BARWClm  + f(1) + div_grav_term

   auxvar%dRp_dp = BARWB + pres(1) * BdARdpWB - BdARdpWClm
   do iface = 1, numfaces
     
     auxvar%dRp_dp = auxvar%dRp_dp + & 
       sq_faces(iface)* &
     ( dukvr_dp(iface)*den(iface)*den(iface) + 2*den(iface)*ukvr(iface)*dden_dp(iface) ) &
     * auxvar%gr(iface) 

   enddo
   auxvar%dRp_dp = auxvar%dRp_dp  + PorVol_dt*sat*dden_cntr_dp + PorVol_dt*den_cntr*ds_dp      

   do iface = 1, auxvar%numfaces
          auxvar%dRp_dl(iface) = 0.
         do jface = 1, auxvar%numfaces
            auxvar%dRp_dl(iface) = auxvar%dRp_dl(iface) - & 
              BAR(jface) * auxvar%MassMatrixInv(iface,jface) * sq_faces(iface)
         enddo
   enddo
    
    auxvar%Rl = 0
    auxvar%dRl_dp = 0
    CWCl = 0

   do iface = 1, auxvar%numfaces
     do jface = 1, auxvar%numfaces
       CWCl(iface) = CWCl(iface) + sq_faces(iface)*auxvar%MassMatrixInv(iface, jface) * &
                                   sq_faces(jface)*face_pres(jface)
     enddo
   enddo

   do iface = 1, auxvar%numfaces

     auxvar%Rl(iface) = -sq_faces(iface)*auxvar%WB(iface)*pres(1) &       
                                   + CWCl(iface)  & 
                                   + sq_faces(iface)*Wg(iface) & 
                                   - sq_faces(iface)*den(iface)*auxvar%gr(iface)&
                                   - bc_h(iface)/ukvr(iface)

     auxvar%dRl_dp(iface) = -sq_faces(iface)*auxvar%WB(iface) - sq_faces(iface)*dden_dp(iface)*auxvar%gr(iface)

   enddo

  do iface = 1, auxvar%numfaces
      rhs(iface) = -auxvar%dRl_dp(iface)*auxvar%Rp/auxvar%dRp_dp + auxvar%Rl(iface)
  enddo


end subroutine MFDAuxGenerateRhs

! ************************************************************************** !

subroutine MFDAuxGenerateRhs_LP(patch, grid, ghosted_cell_id, PermTensor, bc_g, source_f,  bc_h, auxvar, &
                                       rich_auxvar, global_auxvar, Accum, &
                                       porosity, volume, pres, face_pres, bnd,&
                                       sq_faces, neig_den, neig_kvr, neig_dkvr_dp, neig_pres, option, rhs)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module
 use Patch_module

  implicit none

  type(patch_type) :: patch
  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: auxvar
  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscScalar :: sq_faces(:), neig_den(:), neig_kvr(:), neig_dkvr_dp(:), neig_pres(:)
  type(option_type) :: option
  PetscScalar :: bc_g(:), rhs(:), bc_h(:), face_pres(:), bnd(:)
  PetscScalar :: Accum(1:option%nflowdof),source_f(1:option%nflowdof), pres(1:option%nflowdof)
  PetscScalar :: PermTensor(3,3), sat
  PetscInt :: ghosted_cell_id, bound_id



  PetscScalar :: Kg(3), dir_norm(3)
  PetscInt :: i, j, ghost_face_id, local_face_id
  type(connection_set_type), pointer :: conn
  PetscInt, parameter :: numfaces = 6



  PetscInt :: iface, jface
  PetscScalar :: E, den_gWB, den_BWB, den_BWCl , div_grav_term
  PetscScalar :: ukvr(6), dukvr_dp(6), ds_dp, porosity, volume, den_cntr, dden_cntr_dp
  PetscScalar :: PorVol_dt
  PetscScalar :: norm_len, gravity, dphi, distance_gravity
  PetscScalar :: Wg(6), WClm(6), CWCl(6), f(6), dden_dp(6), dir, v_darcy(6), dkvr_dp_neig(6)
  PetscScalar :: dbeta_dp(6), den(6), beta(6), denB(6), Clm(6), Rgr(6), BAR(6), BdARdp(6)
  PetscScalar :: BARWB, BARWClm, BdARdpWB, BdARdpWClm
  PetscScalar :: upweight, delta_p

  ukvr = 0. 
  dukvr_dp = 0.
  dkvr_dp_neig = 0.
  den_cntr = global_auxvar%den(1)
  dden_cntr_dp =  rich_auxvar%dden_dp 

  upweight = 0.5
 
  do i = 1 , numfaces
    ghost_face_id = auxvar%face_id_gh(i)
    conn => grid%faces(ghost_face_id)%conn_set_ptr
    iface = grid%faces(ghost_face_id)%id

    if (bnd(i)==1) then 
      call MFDComputeDensity(global_auxvar, face_pres(i) + bc_g(i)/sq_faces(i), den(i), dden_dp(i), option)
    else
      den(i) =  upweight*den_cntr + (1.D0-upweight)*neig_den(i)
      dden_dp(i) = dden_cntr_dp
    endif

!   UPWIND
    distance_gravity = conn%dist(0,iface) * &
                       dot_product(option%gravity, &
                                     conn%dist(1:3,iface))
    gravity = den(i) * FMWH2O * distance_gravity
    if (auxvar%gr(i) * gravity < 0) gravity = -gravity
    dphi = pres(1) - neig_pres(i) + gravity
 
    if (auxvar%gr(i) >=0 )  then
#ifdef USE_ANISOTROPIC_MOBILITY
      if (rich_auxvar%kvr_x < 1e-10 ) then
        ukvr(i) = 1e-10
      else
        ukvr(i) = rich_auxvar%kvr_x
      endif
      dukvr_dp(i) = rich_auxvar%dkvr_x_dp
#else
      if (rich_auxvar%kvr < 1e-10 ) then
        ukvr(i) = 1e-10
      else
        ukvr(i) = rich_auxvar%kvr
      endif
      dukvr_dp(i) = rich_auxvar%dkvr_dp
#endif

    else if (auxvar%gr(i) < 0)  then

      if (neig_kvr(i) < 1e-10) then
        ukvr(i) = 1e-10
      else
        ukvr(i) = neig_kvr(i)
      endif
      dkvr_dp_neig(i) = neig_dkvr_dp(i)
    endif

    beta(i) = ukvr(i)*den(i)
    dbeta_dp(i) = dukvr_dp(i)*den(i) + ukvr(i)*dden_dp(i)

    denB(i) = sq_faces(i)*den(i)
    BAR(i) = sq_faces(i)*ukvr(i)*den(i)
    BdARdp(i) = sq_faces(i)*dbeta_dp(i)
  enddo

  sat = global_auxvar%sat(1)
  ds_dp = rich_auxvar%dsat_dp
  PorVol_dt = porosity*volume/option%flow_dt

  E = 0.
  f(1) = Accum(1) - source_f(1)
  Wg = 0.

  do iface = 1, numfaces
    Clm(iface) = sq_faces(iface)*face_pres(iface) + bc_g(iface)
    Rgr(iface) = den(iface)*auxvar%gr(iface)
  enddo 

  Wg = matmul(auxvar%MassMatrixInv, bc_g)
  WClm = matmul(auxvar%MassMatrixInv, Clm)

  BARWB = dot_product(auxvar%WB, BAR)
  BARWClm = dot_product(WClm, BAR)
  BdARdpWB = dot_product(auxvar%WB, BdARdp)
  BdARdpWClm = dot_product(WClm, BdARdp)
  div_grav_term = dot_product(Rgr, BAR)

  auxvar%Rp = BARWB*pres(1) - BARWClm  + f(1) + div_grav_term
  auxvar%dRp_dp = BARWB + pres(1) * BdARdpWB - BdARdpWClm

  do iface = 1, numfaces
    auxvar%dRp_dp = auxvar%dRp_dp + &
       sq_faces(iface)* &
      (dukvr_dp(iface)*den(iface)*den(iface) + 2*den(iface)*ukvr(iface)*dden_dp(iface) ) &
      *auxvar%gr(iface)
  enddo
  auxvar%dRp_dp = auxvar%dRp_dp  + PorVol_dt*sat*dden_cntr_dp + PorVol_dt*den_cntr*ds_dp

  auxvar%dRp_dneig = 0

  do iface = 1, numfaces
    auxvar%dRp_dneig(iface)=auxvar%dRp_dneig(iface) + &
                             sq_faces(iface)*den(iface)* &
                              dkvr_dp_neig(iface)*auxvar%WB(iface)*pres(1)
    auxvar%dRp_dneig(iface)=auxvar%dRp_dneig(iface) - &
                             sq_faces(iface)*den(iface)* &
                              dkvr_dp_neig(iface)*WClm(iface)
    auxvar%dRp_dneig(iface)=auxvar%dRp_dneig(iface) + &
                             sq_faces(iface)*den(iface)*den(iface)* &
                              dkvr_dp_neig(iface)*auxvar%gr(iface)
  enddo

  do iface = 1, auxvar%numfaces
    auxvar%dRp_dl(iface) = 0.
    do jface = 1, auxvar%numfaces
      auxvar%dRp_dl(iface) = auxvar%dRp_dl(iface) - &
      BAR(jface) * auxvar%MassMatrixInv(iface,jface) * sq_faces(iface)
    enddo
  enddo
    
  auxvar%Rl = 0
  auxvar%dRl_dp = 0
  CWCl = 0

  do iface = 1, auxvar%numfaces
    do jface = 1, auxvar%numfaces
      CWCl(iface) = CWCl(iface) + sq_faces(iface)*auxvar%MassMatrixInv(iface, jface) * &
                                   sq_faces(jface)*face_pres(jface)
    enddo
  enddo

  do iface = 1, auxvar%numfaces
    auxvar%Rl(iface) = -sq_faces(iface)*auxvar%WB(iface)*pres(1) &
                        + CWCl(iface)  &
                        + sq_faces(iface)*Wg(iface) &
                        - sq_faces(iface)*den(iface)*auxvar%gr(iface)&
                        - bc_h(iface)/ukvr(iface)
    auxvar%dRl_dp(iface) = -sq_faces(iface)*auxvar%WB(iface) - sq_faces(iface)*dden_dp(iface)*auxvar%gr(iface)
  enddo

  do iface = 1, auxvar%numfaces
    rhs(iface) = auxvar%Rl(iface)
  enddo
  rhs(auxvar%numfaces + 1) = auxvar%Rp


end subroutine MFDAuxGenerateRhs_LP

! ************************************************************************** !

subroutine MFDAuxJacobianLocal( grid, auxvar, &
                                       rich_auxvar, global_auxvar,  &
                                       sq_faces, option, J)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module

  implicit none

  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: auxvar
  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscScalar, pointer :: sq_faces(:)
  type(option_type) :: option
  PetscScalar, pointer :: J(:)

  PetscInt :: iface, jface
  PetscScalar :: ukvr

#ifdef USE_ANISOTROPIC_MOBILITY
  ukvr = rich_auxvar%kvr_x
#else
  ukvr = rich_auxvar%kvr
#endif

  J = 0.

   do iface = 1, auxvar%numfaces
     do jface = 1, auxvar%numfaces

!    J(jface + (iface - 1)*auxvar%numfaces) = -ukvr*sq_faces(iface)*auxvar%MassMatrixInv(iface,jface)*sq_faces(jface) - &
!                                                 auxvar%dRl_dp(iface)*auxvar%dRp_dl(jface)/auxvar%dRp_dp

!    J(jface + (iface - 1)*auxvar%numfaces) = -ukvr*sq_faces(iface)*auxvar%MassMatrixInv(iface,jface)*sq_faces(jface) - & !TEST
!                                                 auxvar%dRl_dp(iface)*auxvar%dRp_dl(jface)/auxvar%dRp_dp                     !TEST
    J(jface + (iface - 1)*auxvar%numfaces) = -auxvar%dRl_dp(iface)*auxvar%dRp_dl(jface)/auxvar%dRp_dp   &                  !TEST
                                               + sq_faces(iface)*auxvar%MassMatrixInv(iface,jface)*sq_faces(jface)  

    enddo
  enddo

  auxvar%dRp_dneig = 0

end subroutine MFDAuxJacobianLocal

! ************************************************************************** !

subroutine MFDAuxJacobianLocal_LP( grid, auxvar, &
                                       rich_auxvar, global_auxvar,  &
                                       sq_faces, option, J)

  use Option_module
  use Richards_Aux_module
  use Global_Aux_module

  implicit none

  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: auxvar
  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscScalar, pointer :: sq_faces(:)
  type(option_type) :: option
  PetscScalar, pointer :: J(:)

  PetscInt :: iface, jface
  PetscScalar :: ukvr

#ifdef USE_ANISOTROPIC_MOBILITY
  ukvr = rich_auxvar%kvr_x
#else
  ukvr = rich_auxvar%kvr
#endif

  J = 0.
  do iface = 1, auxvar%numfaces
    do jface = 1, auxvar%numfaces
      J(iface + (jface - 1)*(auxvar%numfaces + 1)) = sq_faces(iface)*auxvar%MassMatrixInv(iface,jface)*sq_faces(jface)
    enddo
    J(auxvar%numfaces + 1 + (iface - 1)*(auxvar%numfaces + 1)) = auxvar%dRp_dl(iface)
    J(iface + auxvar%numfaces*(auxvar%numfaces + 1)) = auxvar%dRl_dp(iface)
  enddo

  J(auxvar%numfaces + 1 + auxvar%numfaces*(auxvar%numfaces + 1)) =  auxvar%dRp_dp

end subroutine MFDAuxJacobianLocal_LP

! ************************************************************************** !

subroutine MFDAuxReconstruct(face_pr, source_f, auxvar, rich_auxvar, global_auxvar, Accum, &
                                       sq_faces, option, xx)

  use Option_module
  use Richards_Aux_module
  use Global_Aux_module

  implicit none

  type(mfd_auxvar_type), pointer :: auxvar
  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscScalar, pointer :: sq_faces(:), face_pr(:)
  type(option_type) :: option
  PetscScalar :: Accum(1:option%nflowdof),source_f(1:option%nflowdof), xx(1:option%nflowdof)

  PetscInt :: iface, jface, i,j
  PetscScalar :: ukvr, E, gMB
  PetscScalar, pointer :: MB(:)

  allocate(MB(auxvar%numfaces))

  E = 0

#ifdef USE_ANISOTROPIC_MOBILITY
  ukvr = rich_auxvar%kvr_x
#else
  ukvr = rich_auxvar%kvr
#endif

  do iface = 1, auxvar%numfaces
    MB(iface) = 0.
    do jface = 1, auxvar%numfaces
       MB(iface) = MB(iface) + auxvar%MassMatrixInv(iface,jface)*sq_faces(jface)
    enddo
    E = E + MB(iface)*sq_faces(iface)
  enddo

  E = 1./E
  xx(1) = 0

  do iface = 1, auxvar%numfaces
    xx(1) = xx(1) + E*MB(iface)*sq_faces(iface)*face_pr(iface) 
  enddo 

  deallocate(MB)

end subroutine MFDAuxReconstruct

! ************************************************************************** !

subroutine  MFDAuxUpdateCellPressure(face_pres, face_DELTA_pres, mfd_auxvar,&
                                option, pressure)

  use Option_module
  use Richards_Aux_module
  use Global_Aux_module

  implicit none

  type(mfd_auxvar_type), pointer :: mfd_auxvar
  type(option_type) :: option
  PetscScalar, pointer ::  face_pres(:), face_DELTA_pres(:)
  PetscScalar :: pressure(1:option%nflowdof)
  

  PetscScalar :: DELTA_pres
  PetscInt :: iface

  DELTA_pres = mfd_auxvar%Rp

  do iface = 1, mfd_auxvar%numfaces
    DELTA_pres = DELTA_pres + mfd_auxvar%dRp_dl(iface)*face_DELTA_pres(iface)
  enddo

  DELTA_pres = -DELTA_pres/mfd_auxvar%dRp_dp
  pressure(1) = pressure(1) + DELTA_pres

end subroutine MFDAuxUpdateCellPressure

! ************************************************************************** !

subroutine MFDAuxFluxes(patch, grid, ghosted_cell_id, xx, face_pr, auxvar, PermTensor, rich_auxvar, global_auxvar, &
                                       sq_faces, bnd, neigh_den, neig_kvr,  neig_pres, option)

  use Option_module
  use Richards_Aux_module
  use Global_Aux_module
  use Grid_module
  use Connection_module
  use Patch_module

  implicit none 

  type(grid_type) :: grid
  type(patch_type) :: patch
  type(mfd_auxvar_type), pointer :: auxvar
  type(richards_auxvar_type) :: rich_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscScalar, pointer :: sq_faces(:), face_pr(:), neigh_den(:), neig_pres(:), neig_kvr(:), bnd(:)
  type(option_type) :: option
  PetscScalar :: xx(1:option%nflowdof),  PermTensor(3,3), Kg(3), real_tmp
  PetscInt :: ghosted_cell_id
  PetscScalar :: ukvr0, den0
  PetscInt, parameter :: numfaces = 6

  PetscInt :: iface, jface, i, j, ghost_face_id, local_face_id, local_id, bound_id, dir
  PetscScalar :: den(numfaces), dden_dp(numfaces) , ukvr(numfaces), den_cntr
  PetscScalar :: gravity, darcy_v, dir_norm(3), total_flux, upweight, DIV, dphi, distance_gravity
  type(connection_set_type), pointer :: conn

#ifdef USE_ANISOTROPIC_MOBILITY
  ukvr = rich_auxvar%kvr_x
#else
  ukvr = rich_auxvar%kvr
#endif
  upweight = 0.5


  do i = 1 , numfaces
    ghost_face_id = auxvar%face_id_gh(i)
    conn => grid%faces(ghost_face_id)%conn_set_ptr
    iface = grid%faces(ghost_face_id)%id
    if (bnd(i)==1) then  
      call MFDComputeDensity(global_auxvar, face_pr(i), den(i), dden_dp(i), option)
    else
      den(i) = global_auxvar%den(1)*upweight + (1 - upweight)*neigh_den(i)
    endif 

! UPWIND
    distance_gravity = conn%dist(0,iface) * &
                       dot_product(option%gravity, &
                                     conn%dist(1:3,iface))

    gravity = den(i) * FMWH2O * distance_gravity   

    if (auxvar%gr(i) * gravity < 0) gravity = -gravity 

    dphi = xx(1) - neig_pres(i) + gravity
    
    if (auxvar%gr(i) >= 0)  then

#ifdef USE_ANISOTROPIC_MOBILITY
      if (rich_auxvar%kvr_x < 1e-10 ) then
        ukvr(i) = 1e-10
      else
        ukvr(i) = rich_auxvar%kvr_x
      endif
#else
      if (rich_auxvar%kvr < 1e-10 ) then
        ukvr(i) = 1e-10
      else
        ukvr(i) = rich_auxvar%kvr
      endif
#endif

    else if (auxvar%gr(i) < 0)  then
      if (neig_kvr(i) < 1e-10 ) then
        ukvr(i) = 1e-10
      else
        ukvr(i) = neig_kvr(i)
      endif
   endif
 enddo
 
  do i = 1, auxvar%numfaces
    ghost_face_id = auxvar%face_id_gh(i)
    local_face_id = grid%fG2L(ghost_face_id)
    conn => grid%faces(ghost_face_id)%conn_set_ptr
    iface = grid%faces(ghost_face_id)%id
    if (conn%itype==INTERNAL_CONNECTION_TYPE) local_id = grid%nG2L(conn%id_up(iface))

    dir = 1
    if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(iface)==ghosted_cell_id) dir = -1

    darcy_v = 0.
    do j = 1, auxvar%numfaces
      darcy_v = darcy_v + ukvr(i)*auxvar%MassMatrixInv(i, j)* &
                                           (sq_faces(j)*(xx(1) - face_pr(j)))
    enddo

    darcy_v = darcy_v + den(i)*ukvr(i)*auxvar%gr(i)

    DIV = DIV + den(i)*darcy_v*sq_faces(i)

    if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(iface)==ghosted_cell_id.and.local_id>0) cycle
      if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        if (local_face_id > 0) then
          bound_id = grid%fL2B(local_face_id)
          if (bound_id>0) then
            patch%boundary_velocities(option%nphase, bound_id) = -darcy_v
          endif
        endif
    else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
      patch%internal_velocities(option%nphase, iface) = darcy_v * dir
    endif
  enddo

end subroutine MFDAuxFluxes

! ************************************************************************** !

subroutine MFDComputeDensity(global_auxvar, pres, den, dden_dp, option)

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module

  type(global_auxvar_type) :: global_auxvar
  PetscScalar :: pres
  type(option_type) :: option
  PetscReal, parameter :: tol = 1.d-3


  PetscInt :: numfaces, i
  PetscErrorCode :: ierr
  PetscReal :: den, dden_dp, pc, pw 
  PetscReal :: dw_kg, dw_mol, dw_dp, dw_dt, hw, hw_dp, hw_dt


  pc = option%reference_pressure - pres

  if (pc > 1.0) then
    pw = option%reference_pressure
  else 
    pw = pres
  endif

#ifndef DONT_USE_WATEOS
  call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol, &
                       dw_dp,dw_dt,ierr)
#else
  call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol,ierr)
  pert = tol*pw
  pw_pert = pw + pert
  call EOSWaterDensity(global_auxvar%temp,pw_pert,dw_kg_pert,dw_mol,ierr)
  dw_dp = (dw_kg_pert-dw_kg)/pert
  ! dw_kg = kg/m^3
  ! dw_mol = kmol/m^3
  ! FMWH2O = kg/kmol h2o
  dw_mol = dw_kg/FMWH2O
  dw_dp = dw_dp/FMWH2O
#endif

  den = dw_mol
  dden_dp = dw_dp

end subroutine MFDComputeDensity

! ************************************************************************** !

subroutine MFDAuxComputeGeometricValues (grid, ghosted_cell_id, auxvar, PermTensor, option)


 use Option_module
 use Richards_Aux_module
 use Global_Aux_module
 use Patch_module

  implicit none

  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: auxvar
  type(option_type) :: option
  PetscScalar :: PermTensor(3,3)
  PetscInt :: ghosted_cell_id

  PetscScalar :: Kg(3), dir_norm(3)
  PetscScalar :: sq_faces(6), norm_len
  PetscInt :: i, j, ghost_face_id
  type(connection_set_type), pointer :: conn
  PetscInt, parameter :: numfaces = 6
  PetscInt :: iface, jface

  Kg = matmul(PermTensor, option%gravity)

  do i=1,3
    Kg(i) = Kg(i) * FMWH2O
  enddo

  do i = 1, numfaces
    ghost_face_id = auxvar%face_id_gh(i)
    conn => grid%faces(ghost_face_id)%conn_set_ptr
    iface = grid%faces(ghost_face_id)%id
    sq_faces(i) = conn%area(iface)

    dir_norm(1) = conn%cntr(1, iface) - grid%x(ghosted_cell_id)       !direction to define outward normal
    dir_norm(2) = conn%cntr(2, iface) - grid%y(ghosted_cell_id)
    dir_norm(3) = conn%cntr(3, iface) - grid%z(ghosted_cell_id)

    norm_len = sqrt(dir_norm(1)**2 + dir_norm(2)**2 + dir_norm(3)**2)

    auxvar%gr(i) = 0
    do j = 1,3
      dir_norm(j) = dir_norm(j)/norm_len
      auxvar%gr(i) = auxvar%gr(i) + dir_norm(j) * Kg(j)
    enddo
  enddo
  auxvar%WB = matmul(auxvar%MassMatrixInv, sq_faces)

end subroutine MFDAuxComputeGeometricValues

! ************************************************************************** !

subroutine MFDAuxGenerateMassMatrixInv(grid, ghosted_cell_id,  auxvar, volume, PermTensor, option)
  ! 
  ! Create a mass matrix for cell
  ! 
  ! Author: Daniil Svyatskiy
  ! Date: 05/25/10
  ! 
#include "finclude/petscmat.h"
 use petscmat
 use Grid_module
 use Option_module
 use Connection_module

  implicit none

  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: auxvar
  PetscScalar :: PermTensor(3,3), dir_norm(3), sq_faces(6)
  PetscScalar :: N(6,3), D(6,3), H(6,3),  W1(6,6), U(3,3), rx(6), ry(6), rz(6)   ! hex only
  PetscScalar :: volume, norm, area, a1, a2, u_parm
  PetscInt :: ghost_face_id, ghosted_cell_id
  type(option_type) :: option
  type(connection_set_type), pointer :: conn


  PetscInt :: iface, i,j,k

  if (volume==0.) then
    option%io_buffer = 'Cell volume iz zero'
    call printErrMsg(option)
  endif

  N = 0.
  u_parm  = (PermTensor(1,1) + PermTensor(2,2) +PermTensor(3,3))/volume

  do i = 1, auxvar%numfaces

    ghost_face_id = auxvar%face_id_gh(i)
    conn => grid%faces(ghost_face_id)%conn_set_ptr
    iface = grid%faces(ghost_face_id)%id
    area = conn%area(iface)
    sq_faces(i) = conn%area(iface)

    dir_norm(1) = conn%cntr(1, iface) - grid%x(ghosted_cell_id)       !direction to define outward normal
    dir_norm(2) = conn%cntr(2, iface) - grid%y(ghosted_cell_id)
    dir_norm(3) = conn%cntr(3, iface) - grid%z(ghosted_cell_id)

    norm = sqrt(dir_norm(1)**2 + dir_norm(2)**2 +dir_norm(3)**2)

    do j=1,3
      dir_norm(j) = dir_norm(j)/norm
    enddo

    rx(i) = ( conn%cntr(1, iface) - grid%x(ghosted_cell_id) )*area
    ry(i) = ( conn%cntr(2, iface) - grid%y(ghosted_cell_id) )*area
    rz(i) = ( conn%cntr(3, iface) - grid%z(ghosted_cell_id) )*area

    do j=1,3
      N(i,j) = dir_norm(j)
     enddo
  enddo

#if 0

  D = 0.

  norm = dot_product(rx(1:6), rx(1:6))
  norm = sqrt(norm)
  do i = 1, auxvar%numfaces
    D(i,1) = rx(i)/norm
  enddo
  a1 = -dot_product(D(1:6,1), ry)

 ! write(*,*) "a1", a1

  do i = 1, auxvar%numfaces
    D(i,2) = ry(i) + a1*D(i,1)
  enddo

  norm = dot_product(D(1:6,2),D(1:6,2))
  norm = sqrt(norm)

  do i = 1, auxvar%numfaces
    D(i,2) = D(i,2)/norm
  enddo

  a1 = -dot_product(D(1:6,1), rz(1:6))
  a2 = -dot_product(D(1:6,2), rz(1:6))

  do i = 1, auxvar%numfaces
    D(i,3) = rz(i) + a1*D(i,1) + a2*D(i,2)
  enddo
   
  norm = dot_product(D(1:6,3),D(1:6,3))
  norm = sqrt(norm)

  do i = 1, auxvar%numfaces
    D(i,3) = D(i,3)/norm
  enddo

  do i = 1, auxvar%numfaces
    do j = 1, auxvar%numfaces
      W1(i,j) = 0.
      do k = 1, 3
         W1(i,j) = W1(i,j) + D(i,k)*D(j,k)
      enddo
    enddo
  enddo

  do i = 1, auxvar%numfaces
    W1(i,i) = W1(i,i) - 1
  enddo

  do i = 1, auxvar%numfaces
     do j = 1, 3
       H(i,j) = 0
       do k = 1, 3
         H(i,j) = H(i,j) + N(i,k)*PermTensor(k,j)
       enddo
     enddo
  enddo

  allocate(auxvar%MassMatrixInv(auxvar%numfaces, auxvar%numfaces))

  do i = 1, auxvar%numfaces
    do j = 1, auxvar%numfaces
       auxvar%MassMatrixInv(i,j) = 0
       do k =1,3
         auxvar%MassMatrixInv(i,j) = auxvar%MassMatrixInv(i,j) + H(i,k) * N(j,k)
       enddo
       auxvar%MassMatrixInv(i,j) = auxvar%MassMatrixInv(i,j) / volume
       auxvar%MassMatrixInv(i,j) = auxvar%MassMatrixInv(i,j) - u_parm*W1(i,j)
     enddo
   enddo

#else

  do i = 1, auxvar%numfaces
     do j = 1, 3
       D(i,j) = abs(N(i,j))
     enddo
  enddo

  do i = 1,3
     U(i,i) = PermTensor(i,i)/volume
     do j = i+1, 3
        U(i,j) = -abs(PermTensor(i,j))/volume
        U(j,i) = U(i,j)
     enddo
  enddo 

  do i = 1, auxvar%numfaces
     do j = 1, 3
       H(i,j) = 0
       do k = 1, 3
         H(i,j) = H(i,j) + N(i,k)*PermTensor(k,j)
       enddo
     enddo
  enddo

  allocate(auxvar%MassMatrixInv(auxvar%numfaces, auxvar%numfaces))


  do i = 1, auxvar%numfaces
    do j = 1, auxvar%numfaces
       auxvar%MassMatrixInv(i,j) = 0
       do k =1,3
         auxvar%MassMatrixInv(i,j) = auxvar%MassMatrixInv(i,j) + H(i,k) * N(j,k)
       enddo
       auxvar%MassMatrixInv(i,j) = auxvar%MassMatrixInv(i,j) / volume
     enddo
   enddo

  do i = 1, auxvar%numfaces
     do j = 1, 3
       H(i,j) = 0
       do k = 1, 3
         H(i,j) = H(i,j) + D(i,k)*U(k,j)
       enddo
    enddo
  enddo

  do i = 1, auxvar%numfaces
    do j = 1, auxvar%numfaces
       do k =1,3
         auxvar%MassMatrixInv(i,j) = auxvar%MassMatrixInv(i,j) + H(i,k) * D(j,k)
       enddo
     enddo
   enddo

#endif

end subroutine MFDAuxGenerateMassMatrixInv


end module MFD_module
