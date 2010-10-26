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

#ifdef DASVYAT

  public :: MFDCreateJacobian, &
            MFDInitializeMassMatrices, MFDAuxGenerateStiffMatrix,&
            MFDAuxGenerateRhs, MFDAuxReconstruct, MFDAuxFluxes, MFDComputeDensity


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
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,ierr)
        call MatSetLocalToGlobalMappingBlock(J,mfd_aux%mapping_ltogb_faces,ierr)
      case(MATBAIJ)
        call MatCreateSeqBAIJ(option%mycomm,mfd_aux%ndof,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,mfd_aux%mapping_ltog_faces,ierr)
        call MatSetLocalToGlobalMappingBlock(J,mfd_aux%mapping_ltogb_faces,ierr)
      case default
        option%io_buffer = 'MatType not recognized in MFDCreateJacobian'
        call printErrMsg(option)
    end select
  endif



  deallocate(d_nnz)
  deallocate(o_nnz)


end subroutine MFDCreateJacobian



subroutine MFDInitializeMassMatrices(grid, volume,  perm_xx_loc, &
                                                    perm_yy_loc, &
                                                    perm_zz_loc, &
 mfd_aux, option)

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
  Vec :: volume, perm_xx_loc, perm_yy_loc, perm_zz_loc
  type(mfd_type) :: mfd_aux
  type(option_type) :: option


  type(mfd_auxvar_type), pointer :: aux_var
  PetscInt :: icell, ierr,i,j
  PetscReal :: PermTensor(3,3) 
  PetscReal, pointer :: volume_p(:), perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)

  call VecGetArrayF90(volume, volume_p, ierr)
  call VecGetArrayF90(perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(perm_zz_loc, perm_zz_loc_p, ierr)

  

  do icell = 1, grid%nlmax

    PermTensor = 0.
    PermTensor(1,1) = perm_xx_loc_p(icell)
    PermTensor(2,2) = perm_yy_loc_p(icell)
    PermTensor(3,3) = perm_zz_loc_p(icell)

!    write(*,*) icell, volume_p(icell), PermTensor(1,1)

    aux_var => mfd_aux%aux_vars(icell)
    call MFDAuxGenerateMassMatrixInv(aux_var, volume_p(icell), PermTensor, option)
    call MFDAuxInitStiffMatrix(aux_var, option)
  end do

  call VecRestoreArrayF90(volume, volume_p, ierr)
  call VecRestoreArrayF90(perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(perm_zz_loc, perm_zz_loc_p, ierr)

!  do icell = 1, grid%nlmax
!    write(*,*) "Mass Matrix ", icell
!    aux_var => mfd_aux%aux_vars(icell)
!    do i = 1, aux_var%numfaces
!       write(*,*) (aux_var%MassMatrixInv(i,j),j=1,aux_var%numfaces)
!    end do
!  end do 

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
  PetscScalar :: ukvr
  PetscScalar, pointer :: MB(:)

  allocate(MB(aux_var%numfaces))

  ukvr = rich_aux_var%kvr_x


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
        aux_var%StiffMatrix(iface,jface) = ukvr*sq_faces(iface)*sq_faces(jface)*    &
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
subroutine MFDAuxGenerateRhs(grid, ghosted_cell_id, PermTensor, bc_g, source_f,  bc_h, aux_var, rich_aux_var, global_aux_var, Accum, &
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
  PetscScalar, pointer :: bc_g(:), rhs(:), bc_h(:)
  PetscScalar :: Accum(1:option%nflowdof),source_f(1:option%nflowdof)
  PetscScalar :: PermTensor(3,3)
  PetscInt :: ghosted_cell_id



  PetscScalar :: Kg(3), dir_norm(3)
  PetscScalar , pointer :: gr(:), f(:)
  PetscInt :: i, j, ghost_face_id
  type(connection_set_type), pointer :: conn


  PetscInt :: iface, jface
  PetscScalar :: E, gMB
  PetscScalar :: ukvr
  PetscScalar, pointer :: MB(:), Mg(:) 

  allocate(MB(aux_var%numfaces))
  allocate(Mg(aux_var%numfaces))
  allocate(f(option%nflowdof))

  ukvr = rich_aux_var%kvr_x

  allocate(gr(aux_var%numfaces))

  Kg = matmul(PermTensor, option%gravity)

  do i=1,3
    Kg(i) = Kg(i) * global_aux_var%den(1) * FMWH2O
  end do

  allocate(gr(aux_var%numfaces))

  do i = 1, aux_var%numfaces


     ghost_face_id = aux_var%face_id_gh(i)
     conn => grid%faces(ghost_face_id)%conn_set_ptr
     iface = grid%faces(ghost_face_id)%id

     dir_norm(1) = grid%x(ghosted_cell_id) - conn%cntr(1, iface)      !direction to define outward normal
     dir_norm(2) = grid%y(ghosted_cell_id) - conn%cntr(2, iface)
     dir_norm(3) = grid%z(ghosted_cell_id) - conn%cntr(3, iface)

     gr(i) = dot_product(Kg, conn%dist(1:3,iface))
     if (dot_product(dir_norm(1:3), conn%dist(1:3,iface)).lt.0) gr(i) =  gr(i) * NEG_ONE_INTEGER

  end do



!  E = Accum(1)

  E = 0.
  f(1) = (source_f(1) + Accum(1))/global_aux_var%den(1)
  
  gMB = 0.
  do iface = 1, aux_var%numfaces
    MB(iface) = 0.
    Mg(iface) = 0.
    do jface = 1, aux_var%numfaces
       MB(iface) = MB(iface) + aux_var%MassMatrixInv(iface,jface)*sq_faces(jface)
       Mg(iface) = Mg(iface) + aux_var%MassMatrixInv(iface,jface)*bc_g(jface)
    end do
    Mg(iface) = Mg(iface) !- gr(iface)
    E = E + MB(iface)*sq_faces(iface)
    gMB = gMB + MB(iface)*bc_g(iface) !- sq_faces(iface)*gr(iface)
  end do

  E = 1./E

  write(*,*) "bc_h"
  write(*,*) (bc_h(iface), iface=1,aux_var%numfaces)

  do iface = 1, aux_var%numfaces
     rhs(iface) = sq_faces(iface)*MB(iface)*E*(f(1) + ukvr*gMB) - ukvr*sq_faces(iface)*Mg(iface) - ukvr*sq_faces(iface)*gr(iface)
  end do


  deallocate(MB)
  deallocate(Mg)
  deallocate(gr)
  deallocate(f)


end subroutine MFDAuxGenerateRhs

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
  PetscScalar :: E, gMB
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


  E = 1./E


  xx(1) = (source_f(1)+Accum(1))*E/global_aux_var%den(1)


  do iface = 1, aux_var%numfaces
    xx(1) = xx(1) + E*MB(iface)*sq_faces(iface)*face_pr(iface)
  end do

!   do iface = 1, aux_var%numfaces
!    write(*,*)  "MFDAuxReconstruct ", "cntr", xx(1), "face", face_pr(iface)
!   end do

  deallocate(MB)


end subroutine MFDAuxReconstruct

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
  PetscScalar :: xx(1:option%nflowdof), ukvr, PermTensor(3,3), Kg(3)
  PetscInt :: ghosted_cell_id
   

  PetscInt :: iface, jface, i, j, ghost_face_id
  PetscScalar, pointer :: gr(:)
  PetscScalar :: gravity, darcy_v, dir_norm(3)
  type(connection_set_type), pointer :: conn



  ukvr = rich_aux_var%kvr_x


  Kg = matmul(PermTensor, option%gravity)

  do i=1,3
    Kg(i) = Kg(i) * global_aux_var%den(1) * FMWH2O
  end do

  allocate(gr(aux_var%numfaces))

  do i = 1, aux_var%numfaces


     ghost_face_id = aux_var%face_id_gh(i)
     conn => grid%faces(ghost_face_id)%conn_set_ptr
     iface = grid%faces(ghost_face_id)%id

     if (conn%itype/=BOUNDARY_CONNECTION_TYPE.and.conn%id_dn(iface)==ghosted_cell_id) cycle

     dir_norm(1) = grid%x(ghosted_cell_id) - conn%cntr(1, iface)      !direction to define outward normal
     dir_norm(2) = grid%y(ghosted_cell_id) - conn%cntr(2, iface)
     dir_norm(3) = grid%z(ghosted_cell_id) - conn%cntr(3, iface)

     gr(i) = dot_product(Kg, conn%dist(1:3,iface))
     if (dot_product(dir_norm(1:3), conn%dist(1:3,iface)).lt.0) gr(i) =  gr(i) * NEG_ONE_INTEGER

!     gr(i) = dot_product(Kg, conn%dist(1:3,iface))
!     if (conn%id_dn(iface) == ghost_id) gr(i) =  gr(i) * NEG_ONE_INTEGER


!     write(*,*) "xx", xx(1), "lm", face_pr(i)
!     write(*,*)  aux_var%MassMatrixInv(i,i), sq_faces(i), ukvr

     darcy_v = 0.
     do j = 1, aux_var%numfaces
        darcy_v = darcy_v + ukvr*aux_var%MassMatrixInv(i, j)* &
                                           (sq_faces(j)*(xx(1) - face_pr(j)))
     end do
     darcy_v = darcy_v - ukvr*gr(i)

     if (conn%itype == BOUNDARY_CONNECTION_TYPE) then
        patch%boundary_velocities(option%nphase, iface) = -darcy_v
!        write(*,*) "bound flux", iface , -darcy_v
     else if (conn%itype == INTERNAL_CONNECTION_TYPE) then
       patch%internal_velocities(option%nphase, iface) = darcy_v
!        write(*,*) "int flux", iface , darcy_v
     end if

  end do
   

  deallocate(gr)
 
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

#endif

end module MFD_module
