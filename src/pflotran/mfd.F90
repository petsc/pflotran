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
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,ierr)
        call MatSetLocalToGlobalMappingBlock(J,mfd_aux%mapping_ltogb_faces,ierr)


      case(MATBAIJ)
        call MatCreateMPIBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,ierr)
        call MatSetLocalToGlobalMappingBlock(J,mfd_aux%mapping_ltogb_faces,ierr)
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
  PetscScalar :: ukvr, den
  PetscScalar, pointer :: MB(:)

  allocate(MB(aux_var%numfaces))

  ukvr = rich_aux_var%kvr_x
  den = global_aux_var%den(1)


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
!        aux_var%StiffMatrix(iface,jface) = den*ukvr*sq_faces(iface)*sq_faces(jface)*    &
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
subroutine MFDAuxGenerateRhs(grid, ghosted_cell_id, PermTensor, bc_g, source_f,  bc_h, aux_var, &
                                       rich_aux_var, global_aux_var, Accum, &
                                       porosity, volume, pres, face_pres, bnd,&
                                       sq_faces, option, rhs)

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
  PetscScalar, pointer :: bc_g(:), rhs(:), bc_h(:), face_pres(:), bnd(:)
  PetscScalar :: Accum(1:option%nflowdof),source_f(1:option%nflowdof), pres(1:option%nflowdof)
  PetscScalar :: PermTensor(3,3), sat, dden_dp, dukvr_dp, dbeta_dp
  PetscInt :: ghosted_cell_id



  PetscScalar :: Kg(3), dir_norm(3)
  PetscScalar , pointer :: gr(:), f(:)
  PetscInt :: i, j, ghost_face_id
  type(connection_set_type), pointer :: conn


  PetscInt :: iface, jface
  PetscScalar :: E, gWB, BWB, BWCl 
  PetscScalar :: ukvr, den, ds_dp, porosity, volume
  PetscScalar :: beta, PorVol_dt
  PetscScalar :: ukvr0, den0, dden_dp0, dukvr_dp0, dbeta_dp0, beta0, norm_len
  PetscScalar, pointer :: WB(:), Wg(:), CWCl(:)

  allocate(WB(aux_var%numfaces))
  allocate(Wg(aux_var%numfaces))
  allocate(CWCl(aux_var%numfaces))
  allocate(f(option%nflowdof))

  ukvr = rich_aux_var%kvr_x
  den = global_aux_var%den(1)

  dukvr_dp = rich_aux_var%dkvr_x_dp
  dden_dp = rich_aux_var%dden_dp 


  ukvr0 = 1123.055414382469
  den0 = 55.35245650628916  
  dden_dp0 = 0.
  dukvr_dp0 = 0. 

  beta0 = ukvr0*den0
  dbeta_dp0 = dukvr_dp0*den0 + ukvr0*dden_dp0


  beta = ukvr*den
  dbeta_dp =  dukvr_dp*den + ukvr*dden_dp

  sat = global_aux_var%sat(1)
  ds_dp = rich_aux_var%dsat_dp
  PorVol_dt = porosity*volume/option%flow_dt

!  write(*,*) "Sat", global_aux_var%sat(1), Accum(1)

  Kg = matmul(PermTensor, option%gravity)



  do i=1,3
    Kg(i) = 0
    do j =1,3
       Kg(i) = Kg(i) + PermTensor(i,j)*option%gravity(j)
    end do
    Kg(i) = Kg(i) * FMWH2O
!    Kg(i) = Kg(i) * den0 * FMWH2O        !TEST
!	Kg(i) = 0
  end do


  allocate(gr(aux_var%numfaces))

  do i = 1, aux_var%numfaces


     ghost_face_id = aux_var%face_id_gh(i)
     conn => grid%faces(ghost_face_id)%conn_set_ptr
     iface = grid%faces(ghost_face_id)%id

     dir_norm(1) = conn%cntr(1, iface) - grid%x(ghosted_cell_id)       !direction to define outward normal
     dir_norm(2) = conn%cntr(2, iface) - grid%y(ghosted_cell_id)  
     dir_norm(3) = conn%cntr(3, iface) - grid%z(ghosted_cell_id)  

     norm_len = sqrt(dir_norm(1)**2 + dir_norm(2)**2 + dir_norm(3)**2)

   !  gr(i) = dot_product(Kg, conn%dist(1:3,iface)) 
   !  write(*,*) conn%dist(1:3,iface)
   !  if (dot_product(dir_norm(1:3), conn%dist(1:3,iface)).lt.0) gr(i) =  gr(i) * NEG_ONE_INTEGER

     gr(i) = 0
     do j = 1,3
       dir_norm(j) = dir_norm(j)/norm_len
       gr(i) = gr(i) + dir_norm(j) * Kg(j)
     end do

!     write(*,*) "grav" , gr(i)
!     write(*,*) "dirichlet" , bc_g(i)

  end do




  E = 0.
  f(1) = Accum(1) - source_f(1)

  BWB = 0.
  gWB = 0.
  WB = 0.
  Wg = 0.


  WB = matmul(aux_var%MassMatrixInv, sq_faces)
  Wg = matmul(aux_var%MassMatrixInv, bc_g)

  BWB = dot_product(WB, sq_faces)
  gWB = dot_product(WB, bc_g)
  
  BWCl = 0
  do iface = 1, aux_var%numfaces
     BWCl = BWCl + WB(iface)*(    sq_faces(iface)*face_pres(iface) + bc_g(iface)   )
!     write(*,*) "pres", pres(1), "lambda", face_pres(iface)
!     write(*,*) "WB(iface)", WB(iface), "sq_faces(iface)", sq_faces(iface), BWCl
  end do

     aux_var%Rp = beta*BWB*pres(1) - beta*BWCl  + f(1) 
 
   aux_var%dRp_dp = beta*BWB + PorVol_dt*den*ds_dp  & 
                     - dbeta_dp*BWCl   &         
                     + dbeta_dp*BWB*pres(1)      &  
                     + PorVol_dt*sat*dden_dp       

!       write(*,*) "aux_var%Rp", aux_var%Rp, "aux_var%dRp_dp", aux_var%dRp_dp, ukvr

!     read(*,*)

   do iface = 1, aux_var%numfaces
          aux_var%dRp_dl(iface) = -beta*sq_faces(iface)*WB(iface)
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

     aux_var%Rl(iface) = -ukvr*sq_faces(iface)*WB(iface)*pres(1) &       ! TEST
                                   + ukvr*CWCl(iface)  & 
                                   + ukvr*sq_faces(iface)*Wg(iface) & 
                                   - ukvr*sq_faces(iface)*den*gr(iface)&
                                   - bc_h(iface)

!     write(*,*) -ukvr*sq_faces(iface)*WB(iface)*pres(1) + ukvr*CWCl(iface) +  ukvr*sq_faces(iface)*Wg(iface), &
!                ukvr*sq_faces(iface)*den*gr(iface)
!   write(*,*) bc_h(iface)  

     aux_var%dRl_dp(iface) = -(ukvr + dukvr_dp*pres(1))*sq_faces(iface)*WB(iface) + dukvr_dp*CWCl(iface) & 
                                                                   + dukvr_dp*sq_faces(iface)*Wg(iface) !& ! TEST
!     if (bc_h(iface)==0) then
          aux_var%dRl_dp(iface) = aux_var%dRl_dp(iface) &
                           - (den*dukvr_dp + dden_dp*ukvr)*sq_faces(iface)*gr(iface)
 !                          - (den*dukvr_dp + dden_dp*ukvr)*sq_faces(iface)*gr(iface)
!     end if 
   end do

!    read(*,*)

  do iface = 1, aux_var%numfaces
      rhs(iface) = -aux_var%dRl_dp(iface)*aux_var%Rp/aux_var%dRp_dp + aux_var%Rl(iface)
!      write(*,*) "rhs" , rhs(iface), aux_var%face_id_gh(iface)
  end do
 !    read(*,*)


  deallocate(WB)
  deallocate(Wg)
  deallocate(gr)
  deallocate(f)

  deallocate(CWCl)

end subroutine MFDAuxGenerateRhs

subroutine MFDAuxJacobianLocal( grid, ghosted_cell_id, PermTensor,  aux_var, &
                                       rich_aux_var, global_aux_var,  &
                                       porosity, volume, pres, face_pres,&
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
  PetscScalar, pointer :: J(:), face_pres(:)
  PetscScalar :: pres(1:option%nflowdof)
  PetscScalar :: PermTensor(3,3)
  PetscInt :: ghosted_cell_id



  PetscScalar :: Kg(3), dir_norm(3)
  PetscScalar , pointer :: gr(:) 
  PetscInt :: i, ghost_face_id
  type(connection_set_type), pointer :: conn


  PetscInt :: iface, jface
  PetscScalar :: E, gWB, BWB, BWCl, Fp,  dFp_dp
  PetscScalar :: ukvr, den, ds_dp, porosity, volume
  PetscScalar :: beta, PorVol_dt, ukvr0, den0
  PetscScalar, pointer :: WB(:)

  allocate(WB(aux_var%numfaces))

  ukvr = rich_aux_var%kvr_x
  den = global_aux_var%den(1)

  ukvr0 = 1123.055414382469
  den0  = 55.35245650628916 

  beta = ukvr*den
  ds_dp = rich_aux_var%dsat_dp
  PorVol_dt = porosity*volume/option%flow_dt


  WB = 0


  WB = matmul(aux_var%MassMatrixInv, sq_faces)


  

  
  J = 0.

   do iface = 1, aux_var%numfaces
     do jface = 1, aux_var%numfaces

!        J(jface + (iface - 1)*aux_var%numfaces) = -ukvr*sq_faces(iface)*aux_var%MassMatrixInv(iface,jface)*sq_faces(jface) - &
!                                                 aux_var%dRl_dp(iface)*aux_var%dRp_dl(jface)/aux_var%dRp_dp

!        J(jface + (iface - 1)*aux_var%numfaces) = -ukvr*sq_faces(iface)*aux_var%MassMatrixInv(iface,jface)*sq_faces(jface) - & !TEST
!                                                 aux_var%dRl_dp(iface)*aux_var%dRp_dl(jface)/aux_var%dRp_dp                     !TEST
        J(jface + (iface - 1)*aux_var%numfaces) = -aux_var%dRl_dp(iface)*aux_var%dRp_dl(jface)/aux_var%dRp_dp   &                  !TEST
                                               + ukvr*sq_faces(iface)*aux_var%MassMatrixInv(iface,jface)*sq_faces(jface)  

     end do
!        J(iface + (iface - 1)*aux_var%numfaces) = J(iface + (iface - 1)*aux_var%numfaces) + 1 !ukvr0
   end do




  deallocate(WB)


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


subroutine  MFDAuxUpdateCellPressure(face_pres, face_DELTA_pres, source_f, mfd_aux_var,&
                               rich_aux_var,global_aux_var, Accum, &
                               porosity, volume, &
                               sq_faces, option, pres)

 use Option_module
 use Richards_Aux_module
 use Global_Aux_module

  implicit none

  type(mfd_auxvar_type), pointer :: mfd_aux_var
  type(richards_auxvar_type) :: rich_aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscScalar, pointer :: sq_faces(:)
  type(option_type) :: option
  PetscScalar, pointer ::  face_pres(:), face_DELTA_pres(:)
  PetscScalar :: Accum(1:option%nflowdof),source_f(1:option%nflowdof), pres(1:option%nflowdof)
  PetscScalar :: porosity, volume
  PetscScalar :: DELTA_pres

  PetscInt :: i, j

  PetscScalar :: BWB, WB(6), BWCl, ukvr



  PetscInt :: iface, jface

  ukvr = rich_aux_var%kvr_x

  BWB = 0.
  WB = 0.


  WB = matmul(mfd_aux_var%MassMatrixInv, sq_faces)

  BWB = dot_product(WB, sq_faces)

 
  BWCl = 0

  do iface = 1, mfd_aux_var%numfaces
     BWCl = BWCl + WB(iface)*sq_faces(iface)*face_pres(iface)
!     write(*,*) "WB(iface)", WB(iface), "sq_faces(iface)", sq_faces(iface), BWCl
  end do



!  do iface = 1, mfd_aux_var%numfaces
!    write(*,*) "face_pres", face_pres(iface), "face_DELTA_pres", face_DELTA_pres(iface)
!  end do



   DELTA_pres = mfd_aux_var%Rp 


  do iface = 1, mfd_aux_var%numfaces
    DELTA_pres = DELTA_pres + mfd_aux_var%dRp_dl(iface)*face_DELTA_pres(iface)
  end do


  DELTA_pres = -DELTA_pres/mfd_aux_var%dRp_dp

! write(*,*) pres(1), DELTA_pres 

  pres(1) = pres(1) + DELTA_pres  

!  pres(1) = BWCl/BWB

!write(*,*) "Final pres", pres(1)

!write(*,*) "CHECK", BWB*pres(1), BWCl,( BWCl - BWB*pres(1))

!read(*,*)




end subroutine MFDAuxUpdateCellPressure


subroutine MFDAuxFluxes(patch, grid, ghosted_cell_id, xx, face_pr, aux_var, PermTensor, rich_aux_var, global_aux_var, &
                                       sq_faces, option)

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
  PetscScalar, pointer :: sq_faces(:), face_pr(:)
  type(option_type) :: option
  PetscScalar :: xx(1:option%nflowdof), ukvr, den, PermTensor(3,3), Kg(3), real_tmp
  PetscInt :: ghosted_cell_id
  PetscScalar :: ukvr0, den0
   

  PetscInt :: iface, jface, i, j, ghost_face_id, local_face_id, local_id, bound_id, dir
  PetscScalar, pointer :: gr(:)
  PetscScalar :: gravity, darcy_v, dir_norm(3), total_flux
  type(connection_set_type), pointer :: conn



  ukvr = rich_aux_var%kvr_x
  den = global_aux_var%den(1)

 ! ukvr = 1123.055414382469
 ! den = 55.35245650628916  

  


  Kg = matmul(PermTensor, option%gravity)

  do i = 1,3 
    Kg(i) = Kg(i) * den * FMWH2O
!	Kg(i) = 0.
  end do

  allocate(gr(aux_var%numfaces))


  do i = 1, aux_var%numfaces


     ghost_face_id = aux_var%face_id_gh(i)
     local_face_id = grid%fG2L(ghost_face_id)
     conn => grid%faces(ghost_face_id)%conn_set_ptr
     iface = grid%faces(ghost_face_id)%id
     if (conn%itype==INTERNAL_CONNECTION_TYPE) local_id = grid%nG2L(conn%id_up(iface)) 

     if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(iface)==ghosted_cell_id.and.local_id>0) cycle
     
      dir = 1

     if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(iface)==ghosted_cell_id) dir = -1
       

     dir_norm(1) = conn%cntr(1, iface) -  grid%x(ghosted_cell_id)     !direction to define outward normal
     dir_norm(2) = conn%cntr(2, iface) -  grid%y(ghosted_cell_id)
     dir_norm(3) = conn%cntr(3, iface) -  grid%z(ghosted_cell_id)

     gr(i) = dot_product(Kg, conn%dist(1:3,iface))
     if (dot_product(dir_norm(1:3), conn%dist(1:3,iface)).lt.0) gr(i) =  gr(i) * NEG_ONE_INTEGER

!     gr(i) = dot_product(Kg, conn%dist(1:3,iface))
!     if (conn%id_dn(iface) == ghost_id) gr(i) =  gr(i) * NEG_ONE_INTEGER


!     write(*,*) "xx", xx(1), "lm", face_pr(i)
!     write(*,*) "M", aux_var%MassMatrixInv(i,i), "sq", sq_faces(i), "ukvr", ukvr
!     write(*,*) "Dq", aux_var%MassMatrixInv(i,i)*sq_faces(i)

     darcy_v = 0.
     do j = 1, aux_var%numfaces
        darcy_v = darcy_v + ukvr*aux_var%MassMatrixInv(i, j)* &
                                           (sq_faces(j)*(xx(1) - face_pr(j)))
     end do

     darcy_v = darcy_v + ukvr*gr(i)

     if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        if (local_face_id > 0) then
           bound_id = grid%fL2B(local_face_id)
           if (bound_id>0) then
              patch%boundary_velocities(option%nphase, bound_id) = -darcy_v
           end if
        end if
     else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
       patch%internal_velocities(option%nphase, iface) = darcy_v * dir
     end if

  end do

  deallocate(gr)
 
!   if (option%myrank==1) then
!       write(*,*) "End of MFDAuxFluxes"
!       write(*,*)
!   end if

!  stop
 
end subroutine MFDAuxFluxes

subroutine MFDComputeDensity(global_aux_var, face_pr, option)

  use Option_module
  use Global_Aux_module
  use water_eos_module



  type(global_auxvar_type) :: global_aux_var
  PetscScalar, pointer :: face_pr(:)
  type(option_type) :: option


  PetscInt :: numfaces, i
  PetscReal :: pres, pc, pw, dw_kg, dw_mol


  numfaces = 6                 !hex only
  pres = 0.
  do i=1,6
    pres = pres + face_pr(i)
  end do
  pres = pres/numfaces

  pc = option%reference_pressure - pres

  if (pc > 1.0) then
    pw = option%reference_pressure
  else 
    pw = pres
  end if

  call density(option%reference_temperature, pw, dw_kg) 

  dw_mol = dw_kg/FMWH2O

  global_aux_var%den = dw_mol
  global_aux_var%den_kg = dw_kg
  

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
  PetscScalar :: PermTensor(3,3), dir_norm(3)
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
      write(*,*) (N(i,j),j=1,3), "     ", (conn%cntr(k, iface), k=1,3)
    end do
#endif
  !write(*,*) "R"

!  do i=1,6
 !  write(*,*) rx(i), ry(i), rz(i)
!  end do




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

 ! write(*,*) "D"
!  do i=1,6
 !    write(*,*) D(i,1:3)
!  end do
 ! write(*,*)


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

#ifdef DASVYAT_DEBUG

!  write(*,*)
!  do i=1,3
!     write(*,*) PermTensor(i, 1:3)
!  end do



!  write(*,*) "vol", volume
!  write(*,*) "MassMatrix"
!  do i=1,6
!     write(*,*) aux_var%MassMatrixInv(i, 1:6)
!  end do
!  write(*,*)
!  write(*,*)
!  stop
#endif


!  aux_var%MassMatrixInv = 0.
!
!  do iface = 1, aux_var%numfaces
!    aux_var%MassMatrixInv(iface, iface) = (2./volume)*PermTensor(1,1)
!  end do 
 ! stop
end subroutine MFDAuxGenerateMassMatrixInv


end module MFD_module
