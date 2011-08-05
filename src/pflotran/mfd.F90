module MFD_module


  use Connection_module
  use Grid_module
  use MFD_Aux_module
  implicit none
  
  private 

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"


  public :: MFDCreateJacobian, &
            MFDInitializeMassMatrices, MFDAuxGenerateStiffMatrix,&
            MFDAuxGenerateRhs, MFDAuxReconstruct, MFDAuxFluxes, MFDComputeDensity,&
            MFDAuxJacobianLocal, MFDAuxUpdateCellPressure


contains

! ************************************************************************** !
!
! MFDCreateJacobian: Creates a Jacobian matrix for  the faced based dof
! author: Daniil Svyatskiy 
! date: 05/14/10
!
! ************************************************************************** !


subroutine MFDCreateJacobian(grid, mfd_aux, mat_type, J, option)

 use Option_module
 use Grid_module
 use MFD_Aux_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  type(grid_type) :: grid
  type(mfd_type) :: mfd_aux
  type(option_type) :: option
  MatType :: mat_type
  Mat :: J


  type(mfd_auxvar_type), pointer :: aux_var
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
    aux_var => MFD_aux%aux_vars(icell)
    do icount = 1, aux_var%numfaces
      ghost_face_id = aux_var%face_id_gh(icount)
      local_face_id = grid%fG2L(ghost_face_id)
      if (local_face_id > 0) then
        do jcount = 1, aux_var%numfaces
          if (icount == jcount) cycle
          ghost_face_id_n = aux_var%face_id_gh(jcount)
          local_face_id_n = grid%fG2L(ghost_face_id_n)
          if (local_face_id_n > 0) then
             d_nnz(local_face_id) = d_nnz(local_face_id) + 1
          else
             o_nnz(local_face_id) = o_nnz(local_face_id) + 1
          end if
        end do
      end if
    end do
  end do

!  do iface = 1, grid%nlmax_faces
!    write(*,*) iface, d_nnz(iface), o_nnz(iface)
!  end do 

  ndof_local = mfd_aux%ndof * grid%nlmax_faces

  if (option%mycommsize > 1) then
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*mfd_aux%ndof
        o_nnz = o_nnz*mfd_aux%ndof
        call MatCreateMPIAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces, mfd_aux%mapping_ltog_faces, ierr)
        call MatSetLocalToGlobalMappingBlock(J,mfd_aux%mapping_ltogb_faces, mfd_aux%mapping_ltogb_faces, ierr)


      case(MATBAIJ)
        call MatCreateMPIBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,mfd_aux%mapping_ltog_faces,ierr)
        call MatSetLocalToGlobalMappingBlock(J,mfd_aux%mapping_ltogb_faces, mfd_aux%mapping_ltogb_faces, ierr)
        
!       if (option%myrank==0) then
!       write(*,*) "myrank", option%myrank
!       call MatSetValuesLocal(J, 1, 100, 1, 100, &
!                                        5.0, INSERT_VALUES,ierr)
!       end if
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobian'
        call printErrMsg(option)
    end select
  else
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*mfd_aux%ndof
        call MatCreateSeqAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces, mfd_aux%mapping_ltog_faces, ierr)
        call MatSetLocalToGlobalMappingBlock(J,mfd_aux%mapping_ltogb_faces, mfd_aux%mapping_ltogb_faces, ierr)
      case(MATBAIJ)
        call MatCreateSeqBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,mfd_aux%mapping_ltog_faces,ierr)
        call MatSetLocalToGlobalMappingBlock(J,mfd_aux%mapping_ltogb_faces,mfd_aux%mapping_ltogb_faces, ierr)
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobian'
        call printErrMsg(option)
    end select
  endif



  deallocate(d_nnz)
  deallocate(o_nnz)


end subroutine MFDCreateJacobian



subroutine MFDInitializeMassMatrices(grid, field, &
                                             mfd_aux, option)

 use Option_module
 use Grid_module
 use MFD_Aux_module
 use Field_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  type(grid_type) :: grid
  type(field_type) :: field
  type(mfd_type) :: mfd_aux
  type(option_type) :: option


  type(mfd_auxvar_type), pointer :: aux_var
  PetscInt :: ghosted_cell_id, icell, ierr,i,j
  PetscReal :: PermTensor(3,3) 
  PetscReal, pointer :: volume_p(:), perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: perm_xz_loc_p(:), perm_xy_loc_p(:), perm_yz_loc_p(:)

  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecGetArrayF90(field%perm_xz_loc, perm_xz_loc_p, ierr)
  call VecGetArrayF90(field%perm_xy_loc, perm_xy_loc_p, ierr)
  call VecGetArrayF90(field%perm_yz_loc, perm_yz_loc_p, ierr)

  

  do icell = 1, grid%nlmax

    ghosted_cell_id = grid%nL2G(icell)

    PermTensor = 0.
    PermTensor(1,1) = perm_xx_loc_p(ghosted_cell_id)
    PermTensor(2,2) = perm_yy_loc_p(ghosted_cell_id)
    PermTensor(3,3) = perm_zz_loc_p(ghosted_cell_id)
    PermTensor(1,3) = perm_xz_loc_p(ghosted_cell_id)
    PermTensor(1,2) = perm_xy_loc_p(ghosted_cell_id)
    PermTensor(2,3) = perm_yz_loc_p(ghosted_cell_id)
    PermTensor(2,1) = PermTensor(1,2)
    PermTensor(3,1) = PermTensor(1,3)
    PermTensor(3,2) = PermTensor(2,3)

!    write(*,*) icell, volume_p(icell), PermTensor(1,1)

    aux_var => mfd_aux%aux_vars(icell)
    call MFDAuxGenerateMassMatrixInv(grid, ghosted_cell_id, aux_var, volume_p(icell), PermTensor, option)
    call MFDAuxInitResidDerivArrays(aux_var, option)
!    call MFDAuxInitStiffMatrix(aux_var, option)
  end do

  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xz_loc, perm_xz_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xy_loc, perm_xy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yz_loc, perm_yz_loc_p, ierr)

!  do icell = 1, grid%nlmax
!    write(*,*) "Mass Matrix ", icell
!    aux_var => mfd_aux%aux_vars(icell)
!    do i = 1, aux_var%numfaces
!       write(*,*) (aux_var%MassMatrixInv(i,j),j=1,aux_var%numfaces)
!    end do
!  end do 

!   write(*,*) "MFDInitializeMassMatrices"
!   stop
 
end subroutine MFDInitializeMassMatrices

subroutine MFDAuxGenerateStiffMatrix(aux_var, rich_aux_var, global_aux_var,  &
                                       sq_faces, option)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module

  implicit none


  type(mfd_auxvar_type), pointer :: aux_var
  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscScalar, pointer :: sq_faces(:)
  type(option_type) :: option
  PetscScalar, pointer :: StiffMatrix(:,:)


  PetscInt :: iface, jface, i,j
  PetscScalar :: E
  PetscScalar, pointer :: MB(:)

  allocate(MB(aux_var%numfaces))



  E = 0 
  do iface = 1, aux_var%numfaces
    MB(iface) = 0.
    do jface = 1, aux_var%numfaces
       MB(iface) = MB(iface) + aux_var%MassMatrixInv(iface,jface)*sq_faces(jface)
    end do
    E = E + MB(iface)*sq_faces(iface)
  end do

  do iface = 1, aux_var%numfaces
    do jface = 1, aux_var%numfaces
        aux_var%StiffMatrix(iface,jface) = sq_faces(iface)*sq_faces(jface)*    &
                                        (aux_var%MassMatrixInv(iface,jface) - &
                                        (1./E)*MB(iface)*MB(jface))
    end do
  end do


  deallocate(MB)
!    write(*,*) "StiffMatrix"
!    do i = 1, aux_var%numfaces
!       write(*,*) (aux_var%StiffMatrix(i,j),j=1,aux_var%numfaces)
!    end do
!    write(*,*)
!    write(*,*) "MassMatrix"
!    do i = 1, aux_var%numfaces
!       write(*,*) (aux_var%MassMatrixInv(i,j),j=1,aux_var%numfaces)
!    end do
!    write(*,*)

end subroutine MFDAuxGenerateStiffMatrix


!subroutine MFDAuxGenerateRhs(ghosted_cell_id, bc_g, source_f, grid,  PermTensor, aux_var, rich_aux_var, global_aux_var, Accum, &
!                                       sq_faces, option, rhs)
subroutine MFDAuxGenerateRhs(patch, grid, ghosted_cell_id, PermTensor, bc_g, source_f,  bc_h, aux_var, &
                                       rich_aux_var, global_aux_var, Accum, &
                                       porosity, volume, pres, face_pres, bnd,&
                                       sq_faces, neig_den, neig_kvr, nieg_dkvr_dp, option, rhs)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module
 use Patch_module

  implicit none

  type(patch_type) :: patch
  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: aux_var
  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscScalar :: sq_faces(:), neig_den(:), neig_kvr(:), nieg_dkvr_dp(:)
  type(option_type) :: option
  PetscScalar :: bc_g(:), rhs(:), bc_h(:), face_pres(:), bnd(:)
  PetscScalar :: Accum(1:option%nflowdof),source_f(1:option%nflowdof), pres(1:option%nflowdof)
  PetscScalar :: PermTensor(3,3), sat
!  PetscScalar, pointer :: dden_dp(:), dbeta_dp(:), den(:), beta(:)
  PetscInt :: ghosted_cell_id, bound_id



  PetscScalar :: Kg(3), dir_norm(3)
!  PetscScalar , pointer :: gr(:), f(:)
  PetscInt :: i, j, ghost_face_id, local_face_id
  type(connection_set_type), pointer :: conn
  PetscInt, parameter :: numfaces = 6



  PetscInt :: iface, jface
  PetscScalar :: E, den_gWB, den_BWB, den_BWCl , div_grav_term
  PetscScalar :: ukvr(6), dukvr_dp(6), ds_dp, porosity, volume, den_cntr, dden_cntr_dp
  PetscScalar :: PorVol_dt
  PetscScalar :: norm_len
!  PetscScalar, pointer :: WB(:), Wg(:), CWCl(:),denB(:)
  PetscScalar :: WB(6), Wg(6), WClm(6), CWCl(6), f(6), dden_dp(6), dir, v_darcy(6)
  PetscScalar :: dbeta_dp(6), den(6), beta(6), denB(6), gr(6), Clm(6), Rgr(6), BAR(6), BdARdp(6)
  PetscScalar :: BARWB, BARWClm, BdARdpWB, BdARdpWClm
  PetscScalar :: upweight, delta_p

!  allocate(WB(numfaces))
!  allocate(Wg(numfaces))
!  allocate(CWCl(numfaces))
!  allocate(f(option%nflowdof))
!  allocate(dden_dp(numfaces))
!  allocate(dbeta_dp(numfaces))
!  allocate(den(numfaces))
!  allocate(beta(numfaces))
!  allocate(denB(numfaces))
!  allocate(gr(numfaces))

  ukvr = 0. !rich_aux_var%kvr_x
  dukvr_dp = rich_aux_var%dkvr_x_dp
  den_cntr = global_aux_var%den(1)
  dden_cntr_dp =  rich_aux_var%dden_dp 
  
  Kg = matmul(PermTensor, option%gravity)

  do i=1,3
    Kg(i) = Kg(i) * FMWH2O
  end do

  do i = 1, numfaces

     ghost_face_id = aux_var%face_id_gh(i)
     local_face_id = grid%fG2L(ghost_face_id)
     conn => grid%faces(ghost_face_id)%conn_set_ptr
     iface = grid%faces(ghost_face_id)%id

     
      dir = 1

     if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(iface)==ghosted_cell_id) dir = -1

     if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        if (local_face_id > 0) then
           bound_id = grid%fL2B(local_face_id)
           if (bound_id>0) then
              v_darcy(i) = - patch%boundary_velocities(option%nphase, bound_id)
           end if
        end if
     else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
       v_darcy(i) = dir*patch%internal_velocities(option%nphase, iface)
     end if
 

     dir_norm(1) = conn%cntr(1, iface) - grid%x(ghosted_cell_id)       !direction to define outward normal
     dir_norm(2) = conn%cntr(2, iface) - grid%y(ghosted_cell_id)  
     dir_norm(3) = conn%cntr(3, iface) - grid%z(ghosted_cell_id)  

     norm_len = sqrt(dir_norm(1)**2 + dir_norm(2)**2 + dir_norm(3)**2)

     gr(i) = 0
     do j = 1,3
       dir_norm(j) = dir_norm(j)/norm_len
       gr(i) = gr(i) + dir_norm(j) * Kg(j)
     end do

  end do


  upweight = 0.5
  

  do i = 1 , numfaces
    if (bnd(i)==1) then 
        call MFDComputeDensity(global_aux_var, face_pres(i) + bc_g(i)/sq_faces(i), den(i), dden_dp(i), option)
    else  

        den(i) =  upweight*den_cntr + (1.D0-upweight)*neig_den(i)
        dden_dp(i) = dden_cntr_dp

    end if

   if (v_darcy(i) > 0)  then
        ukvr(i) = rich_aux_var%kvr_x
        dukvr_dp(i) = rich_aux_var%dkvr_x_dp 
   else 
        ukvr(i) = neig_kvr(i)
        dukvr_dp(i) = 0.!nieg_dkvr_dp(i)
   end if

!      ukvr(i) = upweight*rich_aux_var%kvr_x + (1.D0-upweight)*neig_kvr(i)
!      dukvr_dp(i) = upweight*rich_aux_var%dkvr_x_dp + (1.D0-upweight)* nieg_dkvr_dp(i)    

!     ukvr = rich_aux_var%kvr_x
!     dukvr_dp = rich_aux_var%dkvr_x_dp

     beta(i) = ukvr(i)*den(i)
     dbeta_dp(i) = dukvr_dp(i)*den(i) + ukvr(i)*dden_dp(i)


     denB(i) = sq_faces(i)*den(i)
     BAR(i) = sq_faces(i)*ukvr(i)*den(i)
     BdARdp(i) = sq_faces(i)*dbeta_dp(i)
  end do

#ifdef DASVYAT_DEBUG
   if (ghosted_cell_id==1) then
     write(*,*) "cntr pressure", pres(1)
     write(*,*) "sq ", (sq_faces(iface),iface=1,6)
     write(*,*) "p ", (face_pres(iface) + bc_g(iface)/sq_faces(iface),iface=1,6)
     write(*,*) "den ", (den(iface),iface=1,6)
     write(*,*) "den_dp ", (dden_dp(iface),iface=1,6)
     write(*,*) "dbeta_dp ", ( dbeta_dp(i),i=1,6)
     write(*,*) "Accum ", Accum(1)
  end if
#endif


  sat = global_aux_var%sat(1)
  ds_dp = rich_aux_var%dsat_dp
  PorVol_dt = porosity*volume/option%flow_dt

!  write(*,*) "Sat", global_aux_var%sat(1), Accum(1)




  E = 0.
  f(1) = Accum(1) - source_f(1)

  WB = 0.
  Wg = 0.

  do iface = 1, numfaces
     Clm(iface) = sq_faces(iface)*face_pres(iface) + bc_g(iface)
     Rgr(iface) = den(iface)*gr(iface)
  end do 

  WB = matmul(aux_var%MassMatrixInv, sq_faces)
  Wg = matmul(aux_var%MassMatrixInv, bc_g)
  WClm = matmul(aux_var%MassMatrixInv, Clm)

  BARWB = dot_product(WB, BAR)
  BARWClm = dot_product(WClm, BAR)
  BdARdpWB = dot_product(WB, BdARdp)
  BdARdpWClm = dot_product(WClm, BdARdp)
  div_grav_term = dot_product(Rgr, BAR)
  

   aux_var%Rp = BARWB*pres(1) - BARWClm  + f(1) + div_grav_term

!   write(*,*) "aux_var%Rp", aux_var%Rp, BARWB*pres(1), BARWClm, f(1),  div_grav_term

 
   aux_var%dRp_dp = BARWB + pres(1) * BdARdpWB - BdARdpWClm
   do iface = 1, numfaces
     
     aux_var%dRp_dp = aux_var%dRp_dp + & 
       sq_faces(iface)* &
     ( dukvr_dp(iface)*den(iface)*den(iface) + 2*den(iface)*ukvr(iface)*dden_dp(iface) ) &
     * gr(iface) 

   end do
   aux_var%dRp_dp = aux_var%dRp_dp  + PorVol_dt*sat*dden_cntr_dp + PorVol_dt*den_cntr*ds_dp      


   do iface = 1, aux_var%numfaces
          aux_var%dRp_dl(iface) = 0.
         do jface = 1, aux_var%numfaces
            aux_var%dRp_dl(iface) = aux_var%dRp_dl(iface) - & 
              BAR(jface) * aux_var%MassMatrixInv(iface,jface) * sq_faces(iface)
         end do
   end do
    



    aux_var%Rl = 0
    aux_var%dRl_dp = 0
    CWCl = 0

   do iface = 1, aux_var%numfaces
     do jface = 1, aux_var%numfaces
       CWCl(iface) = CWCl(iface) + sq_faces(iface)*aux_var%MassMatrixInv(iface, jface) * &
                                   sq_faces(jface)*face_pres(jface)
     end do
   end do

   do iface = 1, aux_var%numfaces

     aux_var%Rl(iface) = -sq_faces(iface)*WB(iface)*pres(1) &       
                                   + CWCl(iface)  & 
                                   + sq_faces(iface)*Wg(iface) & 
                                   - sq_faces(iface)*den(iface)*gr(iface)&
                                   - bc_h(iface)/ukvr(iface)
     
     aux_var%dRl_dp(iface) = -sq_faces(iface)*WB(iface) - sq_faces(iface)*dden_dp(iface)*gr(iface)
!     end if 
   end do

#ifdef DASVYAT
   if ((ghosted_cell_id == 1)) then
      write(*,*) "Rp ", aux_var%Rp
      write(*,*) "dRp_dp ", aux_var%dRp_dp
!      write(*,*) "Rl ", aux_var%Rl, 
      write(*,*)  "ukvr", ukvr
!      write(*,*) "xp", pres(1)
!      write(*,*) "sq_faces",(sq_faces(iface),iface=1,6)
!      write(*,*) "face_pres",(face_pres(iface),iface=1,6)
!      write(*,*) "bc_g",(bc_g(iface)/sq_faces(iface),iface=1,6)
!      write(*,*) "den", (den(iface),iface=1,6)
!      write(*,*) "Inverse"
!      do jface=1,6
!        write(*,*) (aux_var%MassMatrixInv(iface,jface),iface=1,6)
!      end do
  write(*,*) "aux_var%dRl_dp", (aux_var%dRl_dp(iface),iface=1,6)
  write(*,*) "ukvr*aux_var%dRl_dp", (ukvr(iface)*aux_var%dRl_dp(iface),iface=1,6)
  write(*,*) "aux_var%Rl", (aux_var%Rl(iface),iface=1,6)
  write(*,*) "ukvr*aux_var%Rl", (ukvr(iface)*aux_var%Rl(iface),iface=1,6)

      read(*,*)
   end if
#endif


  

  do iface = 1, aux_var%numfaces
      rhs(iface) = -aux_var%dRl_dp(iface)*aux_var%Rp/aux_var%dRp_dp + aux_var%Rl(iface)
!      write(*,*) "rhs" , rhs(iface), aux_var%face_id_gh(iface)
  end do
 !    read(*,*)

!  write(*,*) ghosted_cell_id, "rhs", (ukvr(iface)*rhs(iface),iface=1,6)

!  deallocate(WB)
!  deallocate(Wg)
!  deallocate(gr)
!  deallocate(f)
!
!  deallocate(CWCl)
!  deallocate(dden_dp)
!  deallocate(dbeta_dp)
!  deallocate(den)
!  deallocate(beta)
!  deallocate(denB)

end subroutine MFDAuxGenerateRhs


subroutine MFDAuxJacobianLocal( grid, aux_var, &
                                       rich_aux_var, global_aux_var,  &
                                       sq_faces, option, J)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module

  implicit none

  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: aux_var
  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscScalar, pointer :: sq_faces(:)
  type(option_type) :: option
  PetscScalar, pointer :: J(:)





  PetscInt :: iface, jface
  PetscScalar :: ukvr


  ukvr = rich_aux_var%kvr_x

  J = 0.

   do iface = 1, aux_var%numfaces
     do jface = 1, aux_var%numfaces

!        J(jface + (iface - 1)*aux_var%numfaces) = -ukvr*sq_faces(iface)*aux_var%MassMatrixInv(iface,jface)*sq_faces(jface) - &
!                                                 aux_var%dRl_dp(iface)*aux_var%dRp_dl(jface)/aux_var%dRp_dp

!        J(jface + (iface - 1)*aux_var%numfaces) = -ukvr*sq_faces(iface)*aux_var%MassMatrixInv(iface,jface)*sq_faces(jface) - & !TEST
!                                                 aux_var%dRl_dp(iface)*aux_var%dRp_dl(jface)/aux_var%dRp_dp                     !TEST
        J(jface + (iface - 1)*aux_var%numfaces) = -aux_var%dRl_dp(iface)*aux_var%dRp_dl(jface)/aux_var%dRp_dp   &                  !TEST
                                               + sq_faces(iface)*aux_var%MassMatrixInv(iface,jface)*sq_faces(jface)  

     end do
!        J(iface + (iface - 1)*aux_var%numfaces) = J(iface + (iface - 1)*aux_var%numfaces) + 1 !ukvr0
   end do

!   write(*,*) (J(iface + (iface - 1)*aux_var%numfaces), iface = 1,6)
!   write(*,*)

#ifdef DASVYAT_DEBUG
   do iface = 1, aux_var%numfaces
      if (J(iface + (iface - 1)*aux_var%numfaces) < 0) then 
         write(*,*) "Negative on diagonal"
         stop
      end if
   end do
#endif

end subroutine MFDAuxJacobianLocal

subroutine MFDAuxReconstruct(face_pr, source_f, aux_var, rich_aux_var, global_aux_var, Accum, &
                                       sq_faces, option, xx)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module

  implicit none


  type(mfd_auxvar_type), pointer :: aux_var
  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscScalar, pointer :: sq_faces(:), face_pr(:)
  type(option_type) :: option
  PetscScalar :: Accum(1:option%nflowdof),source_f(1:option%nflowdof), xx(1:option%nflowdof)

  PetscInt :: iface, jface, i,j
  PetscScalar :: ukvr, E, gMB
  PetscScalar, pointer :: MB(:)

  allocate(MB(aux_var%numfaces))



  E = 0

  ukvr = rich_aux_var%kvr_x


  do iface = 1, aux_var%numfaces
    MB(iface) = 0.
    do jface = 1, aux_var%numfaces
       MB(iface) = MB(iface) + aux_var%MassMatrixInv(iface,jface)*sq_faces(jface)
    end do
    E = E + MB(iface)*sq_faces(iface)
  end do


  E = 1./E


 ! write(*,*) "Source + Accum", source_f(1)+Accum(1)

!if (dabs(global_aux_var%den(1))< 1e-10) then
   xx(1) = 0 
! else
!  xx(1) = (source_f(1)+Accum(1))*E/(ukvr*global_aux_var%den(1))
! end if

  

  do iface = 1, aux_var%numfaces
    xx(1) = xx(1) + E*MB(iface)*sq_faces(iface)*face_pr(iface) 
  end do 

!   do iface = 1, aux_var%numfaces
!    write(*,*)  "MFDAuxReconstruct ", "cntr", xx(1), "face", face_pr(iface)
!   end do

  deallocate(MB)


end subroutine MFDAuxReconstruct


subroutine  MFDAuxUpdateCellPressure(face_pres, face_DELTA_pres, mfd_aux_var,&
                                option, pressure)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module

  implicit none

  type(mfd_auxvar_type), pointer :: mfd_aux_var
  type(option_type) :: option
  PetscScalar, pointer ::  face_pres(:), face_DELTA_pres(:)
  PetscScalar :: pressure(1:option%nflowdof)
  

  PetscScalar :: DELTA_pres
  PetscInt :: iface



   DELTA_pres = mfd_aux_var%Rp 


  do iface = 1, mfd_aux_var%numfaces
    DELTA_pres = DELTA_pres + mfd_aux_var%dRp_dl(iface)*face_DELTA_pres(iface)
  end do


  DELTA_pres = -DELTA_pres/mfd_aux_var%dRp_dp

! write(*,*) pres(1), DELTA_pres 

  pressure(1) = pressure(1) + DELTA_pres  

!  pres(1) = BWCl/BWB

!write(*,*) "Final pres", pres(1)

!write(*,*) "CHECK", BWB*pres(1), BWCl,( BWCl - BWB*pres(1))

!read(*,*)




end subroutine MFDAuxUpdateCellPressure


subroutine MFDAuxFluxes(patch, grid, ghosted_cell_id, xx, face_pr, aux_var, PermTensor, rich_aux_var, global_aux_var, &
                                       sq_faces, neigh_den,  option)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module
 use Grid_module
 use Connection_module
 use Patch_module

  implicit none

  type(grid_type) :: grid
  type(patch_type) :: patch
  type(mfd_auxvar_type), pointer :: aux_var
  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscScalar, pointer :: sq_faces(:), face_pr(:), neigh_den(:)
  type(option_type) :: option
  PetscScalar :: xx(1:option%nflowdof), ukvr, PermTensor(3,3), Kg(3), real_tmp
  PetscInt :: ghosted_cell_id
  PetscScalar :: ukvr0, den0
  PetscInt, parameter :: numfaces = 6
   

  PetscInt :: iface, jface, i, j, ghost_face_id, local_face_id, local_id, bound_id, dir
  PetscScalar :: gr(numfaces), den(numfaces), dden_dp(numfaces) 
  PetscScalar :: gravity, darcy_v, dir_norm(3), total_flux, upweight
  type(connection_set_type), pointer :: conn

!  allocate(gr(numfaces))
!  allocate(den(numfaces))
!  allocate(dden_dp(numfaces))



  ukvr = rich_aux_var%kvr_x
  den = global_aux_var%den(1)
  upweight = 0.5

 ! ukvr = 1123.055414382469
 ! den = 55.35245650628916  
  do i = 1 , numfaces
     ghost_face_id = aux_var%face_id_gh(i)
     conn => grid%faces(ghost_face_id)%conn_set_ptr
    if (conn%itype == BOUNDARY_CONNECTION_TYPE) then  
         call MFDComputeDensity(global_aux_var, face_pr(i), den(i), dden_dp(i), option)
    else
         den(i) = global_aux_var%den(1)*upweight + (1 - upweight)*neigh_den(i)
    end if 

 end do

  


  Kg = matmul(PermTensor, option%gravity)

  do i = 1,3 
    Kg(i) = Kg(i) * FMWH2O
!	Kg(i) = 0.
  end do



  do i = 1, aux_var%numfaces


     ghost_face_id = aux_var%face_id_gh(i)
     local_face_id = grid%fG2L(ghost_face_id)
     conn => grid%faces(ghost_face_id)%conn_set_ptr
     iface = grid%faces(ghost_face_id)%id
     if (conn%itype==INTERNAL_CONNECTION_TYPE) local_id = grid%nG2L(conn%id_up(iface)) 

     
      dir = 1

     if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(iface)==ghosted_cell_id) dir = -1
       

     dir_norm(1) = conn%cntr(1, iface) -  grid%x(ghosted_cell_id)     !direction to define outward normal
     dir_norm(2) = conn%cntr(2, iface) -  grid%y(ghosted_cell_id)
     dir_norm(3) = conn%cntr(3, iface) -  grid%z(ghosted_cell_id)

     gr(i) = dot_product(Kg, conn%dist(1:3,iface)) * den(i)
     if (dot_product(dir_norm(1:3), conn%dist(1:3,iface)).lt.0) gr(i) =  gr(i) * NEG_ONE_INTEGER

#ifdef DASVYAT
     if (ghosted_cell_id==1) then
        write(*,*) "ghosted_cell_id: ", ghosted_cell_id, "xx", xx(1), "ukvr", ukvr, "den", den(i)
        write(*,*) "gr", gr(i), "lm", face_pr(i), "dir ", dir
     end if
#endif
     if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(iface)==ghosted_cell_id.and.local_id>0) cycle

!     gr(i) = dot_product(Kg, conn%dist(1:3,iface))
!     if (conn%id_dn(iface) == ghost_id) gr(i) =  gr(i) * NEG_ONE_INTEGER


!     write(*,*) "M", aux_var%MassMatrixInv(i,i), "sq", sq_faces(i), "ukvr", ukvr
!     write(*,*) "Dq", aux_var%MassMatrixInv(i,i)*sq_faces(i)

     darcy_v = 0.
     do j = 1, aux_var%numfaces
        darcy_v = darcy_v + ukvr*aux_var%MassMatrixInv(i, j)* &
                                           (sq_faces(j)*(xx(1) - face_pr(j)))


     end do
#ifdef DASVYAT
     if (ghosted_cell_id==1) write(*,*) "flux no gr", darcy_v, "gr",  ukvr*gr(i)
#endif
     darcy_v = darcy_v + ukvr*gr(i)

     if (ghosted_cell_id==1) write(*,*) "flux ", darcy_v

     if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        if (local_face_id > 0) then
           bound_id = grid%fL2B(local_face_id)
           if (bound_id>0) then
              patch%boundary_velocities(option%nphase, bound_id) = -darcy_v
!              write(*,*) bound_id, "fl", -darcy_v, "lm", face_pr(i), "p", xx(1) 
           end if
        end if
     else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
       patch%internal_velocities(option%nphase, iface) = darcy_v * dir
     end if

  end do

!  deallocate(gr)
!  deallocate(den)
!  deallocate(dden_dp)
 
!   if (option%myrank==1) then
!       write(*,*) "End of MFDAuxFluxes"
!       write(*,*)
!       read(*,*)
!   end if

!  stop
 
end subroutine MFDAuxFluxes

subroutine MFDComputeDensity(global_aux_var, pres, den, dden_dp, option)

  use Option_module
  use Global_Aux_module
  use water_eos_module



  type(global_auxvar_type) :: global_aux_var
  PetscScalar :: pres
  type(option_type) :: option
  PetscReal, parameter :: tol = 1.d-3


  PetscInt :: numfaces, i, ierr
  PetscReal :: den, dden_dp, pc, pw 
  PetscReal :: dw_kg, dw_mol, dw_dp, dw_dt, hw, hw_dp, hw_dt


  pc = option%reference_pressure - pres

  if (pc > 1.0) then
    pw = option%reference_pressure
  else 
    pw = pres
  end if

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
#ifndef DONT_USE_WATEOS
  call wateos(global_aux_var%temp(1),pw,dw_kg,dw_mol,dw_dp,dw_dt,hw, &
              hw_dp,hw_dt,option%scale,ierr)
#else
  call density(global_aux_var%temp(1),pw,dw_kg)
  pert = tol*pw
  pw_pert = pw + pert
  call density(global_aux_var%temp(1),pw_pert,dw_kg_pert)
  dw_dp = (dw_kg_pert-dw_kg)/pert
  ! dw_kg = kg/m^3
  ! dw_mol = kmol/m^3
  ! FMWH2O = kg/kmol h2o
  dw_mol = dw_kg/FMWH2O
  dw_dp = dw_dp/FMWH2O
#endif

  den = dw_mol
   
  dden_dp = dw_dp

 
!  den = 55.3d-0 
!  dden_dp = 0
 

end subroutine MFDComputeDensity

! ************************************************************************** !
!
! MFDAuxGenerateMassMatrixInv: Create a mass matrix for cell
! author: Daniil Svyatskiy
! date: 05/25/10
!
! ************************************************************************** !
subroutine MFDAuxGenerateMassMatrixInv(grid, ghosted_cell_id,  aux_var, volume, PermTensor, option)

 use Grid_module
 use Option_module
 use Connection_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"


  type(grid_type) :: grid
  type(mfd_auxvar_type), pointer :: aux_var
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
  end if

   N = 0.
   u_parm  = (PermTensor(1,1) + PermTensor(2,2) +PermTensor(3,3))/volume


   do i = 1, aux_var%numfaces

     ghost_face_id = aux_var%face_id_gh(i)
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
     end do
  end do

#ifdef DASVYAT_DEBUG
      write(*,*) "N"
    do i=1,6
      ghost_face_id = aux_var%face_id_gh(i)
      conn => grid%faces(ghost_face_id)%conn_set_ptr
      iface = grid%faces(ghost_face_id)%id

      write(*,*) (N(i,j),j=1,3), "     ", (conn%cntr(k, iface), k=1,3)
    end do
  write(*,*) "R"

  do i=1,6
   write(*,*) rx(i), ry(i), rz(i)
  end do
#endif


#if 0

 !   write(*,*) 

 !   write(*,*) "volume", volume

  D = 0.

  norm = dot_product(rx(1:6), rx(1:6))
  norm = sqrt(norm)
  do i = 1, aux_var%numfaces
    D(i,1) = rx(i)/norm
  end do
  a1 = -dot_product(D(1:6,1), ry)

 ! write(*,*) "a1", a1

  do i = 1, aux_var%numfaces
    D(i,2) = ry(i) + a1*D(i,1)
  end do

  norm = dot_product(D(1:6,2),D(1:6,2))
  norm = sqrt(norm)

  do i = 1, aux_var%numfaces
    D(i,2) = D(i,2)/norm
  end do

  a1 = -dot_product(D(1:6,1), rz(1:6))
  a2 = -dot_product(D(1:6,2), rz(1:6))


 ! write(*,*) "a1", a1
 ! write(*,*) "a2", a2

  do i = 1, aux_var%numfaces
    D(i,3) = rz(i) + a1*D(i,1) + a2*D(i,2)
  end do
   
  norm = dot_product(D(1:6,3),D(1:6,3))
  norm = sqrt(norm)

  do i = 1, aux_var%numfaces
    D(i,3) = D(i,3)/norm
  end do

#ifdef DASVYAT_DEBUG

  write(*,*) "D"
  do i=1,6
     write(*,*) D(i,1:3)
  end do
  write(*,*)

#endif


  do i = 1, aux_var%numfaces
    do j = 1, aux_var%numfaces
      W1(i,j) = 0.
      do k = 1, 3
         W1(i,j) = W1(i,j) + D(i,k)*D(j,k)
      end do
    end do
  end do

  do i = 1, aux_var%numfaces
    W1(i,i) = W1(i,i) - 1
  end do

 !  write(*,*) "W1"
 !  do i=1,6
 !    write(*,*) W1(i, 1:6)
 !  end do
 ! write(*,*)


  do i = 1, aux_var%numfaces
     do j = 1, 3
       H(i,j) = 0
       do k = 1, 3
         H(i,j) = H(i,j) + N(i,k)*PermTensor(k,j)
       end do
     end do
  end do

  allocate(aux_var%MassMatrixInv(aux_var%numfaces, aux_var%numfaces))


  do i = 1, aux_var%numfaces
    do j = 1, aux_var%numfaces
       aux_var%MassMatrixInv(i,j) = 0
       do k =1,3
         aux_var%MassMatrixInv(i,j) = aux_var%MassMatrixInv(i,j) + H(i,k) * N(j,k)
       end do
       aux_var%MassMatrixInv(i,j) = aux_var%MassMatrixInv(i,j) / volume

       aux_var%MassMatrixInv(i,j) = aux_var%MassMatrixInv(i,j) - u_parm*W1(i,j)
     end do
   end do

#endif

#if 1

  do i = 1, aux_var%numfaces
     do j = 1, 3
       D(i,j) = abs(N(i,j))
     end do
  end do

  do i = 1,3
     U(i,i) = PermTensor(i,i)/volume
     do j = i+1, 3
        U(i,j) = -abs(PermTensor(i,j))/volume
        U(j,i) = U(i,j)
     end do
  end do 

  do i = 1, aux_var%numfaces
     do j = 1, 3
       H(i,j) = 0
       do k = 1, 3
         H(i,j) = H(i,j) + N(i,k)*PermTensor(k,j)
       end do
     end do
  end do

  allocate(aux_var%MassMatrixInv(aux_var%numfaces, aux_var%numfaces))


  do i = 1, aux_var%numfaces
    do j = 1, aux_var%numfaces
       aux_var%MassMatrixInv(i,j) = 0
       do k =1,3
         aux_var%MassMatrixInv(i,j) = aux_var%MassMatrixInv(i,j) + H(i,k) * N(j,k)
       end do
       aux_var%MassMatrixInv(i,j) = aux_var%MassMatrixInv(i,j) / volume
     end do
   end do

  do i = 1, aux_var%numfaces
     do j = 1, 3
       H(i,j) = 0
       do k = 1, 3
         H(i,j) = H(i,j) + D(i,k)*U(k,j)
       end do
    end do
  end do

  do i = 1, aux_var%numfaces
    do j = 1, aux_var%numfaces
       do k =1,3
         aux_var%MassMatrixInv(i,j) = aux_var%MassMatrixInv(i,j) + H(i,k) * D(j,k)
       end do
     end do
   end do



#endif

#ifdef DASVYAT_DEBUG

  write(*,*)
  do i=1,3
     write(*,*) PermTensor(i, 1:3)
  end do

  write(*,*) "sq"
  write(*,*) (sq_faces(i),i=1,6)

  write(*,*) "vol", volume
  write(*,*) "MassMatrix"
  do i=1,6
     write(*,*) aux_var%MassMatrixInv(i, 1:6)
  end do
  write(*,*)
  write(*,*)
  stop
#endif


!  aux_var%MassMatrixInv = 0.
!
!  do iface = 1, aux_var%numfaces
!    aux_var%MassMatrixInv(iface, iface) = (2./volume)*PermTensor(1,1)
!  end do 
 ! stop
end subroutine MFDAuxGenerateMassMatrixInv


end module MFD_module
