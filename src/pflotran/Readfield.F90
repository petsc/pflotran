module readfield

! Note:  subroutine to read in permeability field according to node natural 
!        ordering if the nature order is not appeared in the file, it will 
!        take the old value.  Repeated index will be ignored (only take the 
!        first one)

 use pflow_gridtype_module
 ! use pflow_var_module

 private 

#include "definitions.h"
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
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

   

  public Read_perm_field, Read_init_field, Read_Geom_field, Boundary_adjustment, &
         Natural2LocalIndex
  
contains

!======================================================================
   
subroutine Read_perm_field(grid)

  use Fileio_module

  implicit none

  type(pflowGrid), intent(inout) :: grid
  character(len=MAXSTRINGLENGTH) :: string 

  PetscReal, pointer :: perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),por_p(:), &
                          tor_p(:)
  PetscInt :: iln,na,nx,ny,nz,ir
  PetscErrorCode :: ierr
  PetscReal ::  px,py,pz, por, tor

  PetscInt, parameter :: fid=60
  
  open(fid, file="perm_field.in", action="read", status="old")
  read(fid,*)
     
  call VecGetArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%tor,tor_p,ierr)
  call VecGetArrayF90(grid%porosity,por_p,ierr)

  do 
    call fiReadFlotranString(fid, string, ierr)
    if (ierr /= 0) exit

    call fiReadInt(string,ir,ierr) 
    call fiDefaultMsg(grid%myrank,'na',ierr)
    ir=ir-1


    call fiReadDouble(string,px,ierr)
    call fiDefaultMsg(grid%myrank,'perm_x',ierr)

    call fiReadDouble(string,py,ierr)
    call fiDefaultMsg(grid%myrank,'perm_y',ierr)
   
    call fiReadDouble(string,pz,ierr)
    call fiDefaultMsg(grid%myrank,'perm_z',ierr)

    call fiReadDouble(string,por,ierr)
    call fiDefaultMsg(grid%myrank,'por',ierr)

    call fiReadDouble(string,tor,ierr)
    call fiDefaultMsg(grid%myrank,'tor',ierr)



!  read(fid,*)ir,px,py,pz,por,tor
!   print *, 'Read filed:',ir,nx,ny,nz,px,py,pz,por,tor
 
    do iln=1, grid%nlmax
      na = grid%nL2A(iln)
    
      if (na == ir) then
        nz = na/grid%nxy + 1
        ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
        nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
        perm_xx_p(iln)= px
        perm_yy_p(iln)= py
        perm_zz_p(iln)= pz
        if (por>=0.D0 .AND. por <=1.D0) por_p(iln)=por
        if (tor>=0.D0 .and. tor <=1.D0) tor_p(iln)=tor
         
        exit
      endif
    enddo 
  enddo
  close(fid)
 
 ! select case(card)
 ! case('LOGO' 
     
  call VecRestoreArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)              
  call VecRestoreArrayF90(grid%tor,tor_p,ierr)
  call VecRestoreArrayF90(grid%porosity,por_p,ierr)

end  subroutine Read_perm_field
  
!  subroutine pflowGridRestart(grid, fname, ntstep, kplt, iplot, iflgcut, &
 !                           ihalcnt, its)

subroutine Read_init_field(grid, kplt)

  use Fileio_module
  use Utility_module
  implicit none

  type(pflowGrid), intent(inout) :: grid
!  character(len=MAXSTRINGLENGTH) :: string 
  PetscInt :: kplt

  PetscReal, pointer :: xx_p(:), iphase_p(:)
  PetscInt :: iln,na,nx,ny,nz,ir,n, ii
  PetscErrorCode :: ierr
  PetscReal xxx(1:grid%ndof)
  PetscReal :: x,y,z,phase,pl,pg,t,sl,sg,xl,xg,vf 
  PetscInt, parameter :: fid = 60
  
  open(fid, file="pflow_init.dat", action="read", status="old")
  read(fid)
!  read(fid,'(2i10, 1p2e16.8)') kplt, grid%flowsteps, grid%t, grid%dt
  read(fid)
  read(fid)
       
  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)

  do n=1, grid%nmax 
    ir=n-1
  
!    if (grid%ny==1) then 
!     read(fid,'(1p11e12.4)')x,z,phase,pl,pg,t,sl,sg,xl,xg,vf
 !    print *, x,z,phase,pl,pg,t,sl,sg,xl,xg,vf
 !   else
      read(fid,'(i10, 1p20e16.8)') ii, x,y,z,phase, xxx
!    print *,x,y,z,phase,pg,t,sg,xl,xg,vf
 !   endif 

            call Natural2LocalIndex(ir,iln, grid%nL2A, grid%nlmax)
      if (iln>0) then
        nz = na/grid%nxy + 1
        ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
        nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
        iphase_p(iln)= phase
        xx_p(1+(iln-1)*grid%ndof: iln* grid%ndof) = xxx(:)
!        print *, na, ii,ir, xxx, phase  
   !      if ((dabs(phase-1.D0)<1e-1) .or. (dabs(phase-2.D0)<1e-1)) then  
  !        xx_p(1+(iln-1)*grid%ndof) = pg 
  !        xx_p(2+(iln-1)*grid%ndof) = t
   !       xx_p(3+(iln-1)*grid%ndof) = xg
    !      print*,'ass:', iln,ir,na,xx_p(
    !    elseif (dabs(phase-3.D0)<1e-1) then
    !      xx_p(1+(iln-1)*grid%ndof) = pg 
    !      xx_p(2+(iln-1)*grid%ndof) = t
    !      xx_p(3+(iln-1)*grid%ndof) = sg
    !    else
    !      print *, "error in phase cond:",phase
    !      stop        
    !    endif
        
      endif    
     
  enddo
  
  close(fid)
    
  call VecRestoreArrayF90(grid%xx, xx_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p,ierr)
  
  call VecCopy(grid%xx, grid%yy, ierr)
  call VecCopy(grid%iphas, grid%iphas_old, ierr)
  !propo0to034
end subroutine Read_init_field
  
 
PetscInt function nxyz2na(grid,nx,ny,nz)

  implicit none

  type(pflowGrid), intent(inout) :: grid
  PetscInt :: nx,ny,nz
  
  nxyz2na = nx -1 + (ny-1)*grid%nx + (nz-1)* grid%nxy

end function nxyz2na
 
   
PetscInt function GetLocalGhostedIdFromNaturalId(natural_id,grid)

  implicit none

  type(pflowGrid) :: grid

  PetscInt :: natural_id, local_ghosted_id
  
  do local_ghosted_id = 1, grid%ngmax
    if (natural_id == grid%nG2A(local_ghosted_id)+1) then
      GetLocalGhostedIdFromNaturalId = local_ghosted_id
      return 
    endif
  enddo
  GetLocalGhostedIdFromNaturalId = 0

end function GetLocalGhostedIdFromNaturalId
    
       
           
subroutine Read_Geom_field(grid)

  use Fileio_module
  use Utility_module
  implicit none

  type(pflowGrid), intent(inout) :: grid
  character(len=MAXSTRINGLENGTH) :: string 

 
  PetscReal, pointer :: perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),por_p(:), &
                          tor_p(:), volume_p(:)
  PetscInt :: iln,na,nx,ny,nz,ir,ng
  PetscErrorCode :: ierr
  PetscReal ::  xc,yc,zc,vc, dx,dy,dz
  PetscReal ::  px,py,pz, por, tor
  PetscReal ::  area, dist1, dist2, grav_ang
! PetscReal ::  x,y,z
  character(len=MAXCARDLENGTH) card, word
  ! character(len=MAXSTRINGLENGTH) :: string 
  PetscInt :: length
  PetscInt :: nx1, nx2, ny1, ny2, nz1,nz2, ina1, ina2,nnc, na1,na2, ncna
  PetscInt :: iln1, iln2, ng1,ng2, mconn
  PetscInt, parameter :: fid = 60
  
 ! if (grid%igeom /=-1) then
 !   print *, ' Wrong geomtry entry:: STOP !'
 !   stop
 ! Endif
  
  call VecGetArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%tor,tor_p,ierr)
  call VecGetArrayF90(grid%porosity,por_p,ierr)

  open(IUNIT1, file="geom_field.in", action="read", status="old")
 
  do
    call fiReadFlotranString(IUNIT1, string, ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    length = len_trim(word)
    call fiCharsToUpper(word,length)
    call fiReadCard(word,card,ierr)

    if (grid%myrank == 0) print *, card

    select case(card)

      case('GRID')  
        do 
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(grid%myrank,'GRID',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/') exit

          call fiReadInt(string,ir,ierr) 
          call fiDefaultMsg(grid%myrank,'na',ierr)
          ir=ir-1
     
          call fiReadDouble(string,xc,ierr)
          call fiDefaultMsg(grid%myrank,'x',ierr)
   
          call fiReadDouble(string,yc,ierr)
          call fiDefaultMsg(grid%myrank,'y',ierr)

          call fiReadDouble(string,zc,ierr)
          call fiDefaultMsg(grid%myrank,'z',ierr)
   
          call fiReadDouble(string,vc,ierr)
          call fiDefaultMsg(grid%myrank,'vol',ierr)
   
          call fiReadDouble(string,px,ierr)
          call fiDefaultMsg(grid%myrank,'perm_x',ierr)

          call fiReadDouble(string,py,ierr)
          call fiDefaultMsg(grid%myrank,'perm_y',ierr)
   
          call fiReadDouble(string,pz,ierr)
          call fiDefaultMsg(grid%myrank,'perm_z',ierr)

          call fiReadDouble(string,por,ierr)
          call fiDefaultMsg(grid%myrank,'por',ierr)

          call fiReadDouble(string,tor,ierr)
          call fiDefaultMsg(grid%myrank,'tor',ierr)

!geh          grid%x(ir+1)=xc
!geh          grid%y(ir+1)=yc
!geh          grid%z(ir+1)=zc


!         read(fid,*)ir,px,py,pz,por,tor
!         print *, 'Read geom 0:',ir,xc,yc,zc,vc, px,py, pz,por,tor
          call Natural2LocalIndex(ir,iln, grid%nL2A, grid%nlmax)
            if (iln > 0 ) then
              nz = na/grid%nxy + 1
              ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
              nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
              ng1=grid%nL2G(iln)
              grid%x(ng1)=xc
              grid%y(ng1)=yc
              grid%z(ng1)=zc

              volume_p(iln)=vc
              
              perm_xx_p(iln)= px
              perm_yy_p(iln)= py
              perm_zz_p(iln)= pz
              if (por>=0.D0 .AND. por <= 1.D0) por_p(iln)=por
              if (tor>=0.D0 .and. tor <= 1.D0) tor_p(iln)=tor
!            print *, 'Read geom 1:',ir,iln,xc,yc,zc,vc
            else
             print *, 'not found', ir,iln; stop
            endif

        enddo
!       print *, grid%x
   
   ! then need to compute connection
   
      case('CONN')  
        mconn = grid%nconn
        print *,'structured conn= ', mconn
        grid%nconn =0
        nnc=0
        do 
     
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(grid%myrank,'CONN',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/') exit
   
          call fiReadInt(string,ncna,ierr) 
          call fiDefaultMsg(grid%myrank,'na1',ierr)

          call fiReadInt(string,na1,ierr) 
          call fiDefaultMsg(grid%myrank,'na1',ierr)
          na1=na1-1
     
          call fiReadInt(string,na2,ierr) 
          call fiDefaultMsg(grid%myrank,'na2',ierr)
          na2=na2-1

          call fiReadDouble(string,dist1,ierr)
          call fiDefaultMsg(grid%myrank,'dist 1',ierr)

          call fiReadDouble(string,dist2,ierr)
          call fiDefaultMsg(grid%myrank,'dist 2',ierr)
     
          call fiReadDouble(string,area,ierr)
          call fiDefaultMsg(grid%myrank,'Area',ierr)
   
          call fiReadDouble(string,grav_ang,ierr)
          call fiDefaultMsg(grid%myrank,'cos(B)',ierr)
    
           ng1 = GetLocalGhostedIdFromNaturalId(na1,grid)
           ng2 = GetLocalGhostedIdFromNaturalId(na2,grid)
           
        if (ng1 > 0 .and. ng2 > 0) then ! both cells are local ghosted
               iln1 = grid%nG2L(ng1)
               iln2 = grid%nG2L(ng2)
      ! at least 1 cell must be local (non-ghosted) since we don't want to
      ! create a connection between two ghosted cells
          if (iln1 > 0 .or. iln2 > 0) then
               grid%nconn = grid%nconn + 1
               grid%nd1(grid%nconn) = ng1
               grid%nd2(grid%nconn) = ng2
               grid%dist1(grid%nconn) = dist1
               grid% dist2(grid%nconn) = dist2
               if (area>0.D0) grid%area(grid%nconn)= area
               grid%delz(grid%nconn) = +1.d0*abs(grid%z(ng1)-  &
                                    grid%z(ng2))
               grid%grav_ang(grid%nconn) = grav_ang


        dx = abs(grid%x(ng1)-grid%x(ng2))
        dy = abs(grid%y(ng1)-grid%y(ng2))
        dz = abs(grid%z(ng1)-grid%z(ng2))
    
        if (dx > dy .and. dx > dz) then
          grid%iperm1(grid%nconn) = 1
          grid%iperm2(grid%nconn) = 1
        else if (dy > dx .and. dy > dz) then
          grid%iperm1(grid%nconn) = 2
          grid%iperm2(grid%nconn) = 2
        else if (dz > dx .and. dz > dy) then
          grid%iperm1(grid%nconn) = 3
          grid%iperm2(grid%nconn) = 3
        endif
     !   print *, 'Read conn ',  grid%nconn,na1,na2,iln1,iln2,&
     !         grid%iperm1(grid%nconn), grid%delz(grid%nconn),grid%grav_ang(grid%nconn)
       endif
    endif

!print *, grid%myrank, local_ghosted_id_upwind, local_id_upwind, local_ghosted_id_downwind, local_id_downwind
  enddo

    end select
  enddo
 close(IUNIT1)
 
 ! select case(card)
 ! case('LOGO' 
  call VecRestoreArrayF90(grid%volume, volume_p, ierr); CHKERRQ(ierr)     
  call VecRestoreArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)              
  call VecRestoreArrayF90(grid%tor,tor_p,ierr)
  call VecRestoreArrayF90(grid%porosity,por_p,ierr)

end  subroutine Read_Geom_field
 
  
  
 subroutine Boundary_adjustment(grid)
   implicit none

   type(pflowGrid), intent(inout) :: grid
   
   PetscReal, pointer :: xx_p(:)
   PetscInt :: nc,m,ibc
   PetscErrorCode :: ierr

   call VecGetArrayF90(grid%xx,xx_p,ierr)  
     do nc=1,grid%nconnbc
         m=grid%mblkbc(nc)
        ibc = grid%ibconn(nc)
        if (grid%ibndtyp(ibc) ==3 ) grid%xxbc(1,nc)=xx_p(1+(m-1)*grid%ndof)
   enddo
   call VecRestoreArrayF90(grid%xx, xx_p, ierr)
 end subroutine Boundary_adjustment

end module readfield
