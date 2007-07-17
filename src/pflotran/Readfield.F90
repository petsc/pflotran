module readfield

! Note:  subroutine to read in permeability field according to node natural 
!        ordering if the nature order is not appeared in the file, it will 
!        take the old value.  Repeated index will be ignored (only take the 
!        first one)

 use pflow_gridtype_module
 ! use pflow_var_module

 private 

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
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

   

  public Read_perm_field, Read_init_field, Read_Geom_field, Boundary_adjustment
  
contains
  
subroutine Read_perm_field(grid)

  use fileio_module

  implicit none

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  character(len=MAXSTRINGLENGTH) :: string 

  PetscScalar, pointer :: perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),por_p(:), &
                          tor_p(:)
  integer iln,na,nx,ny,nz,ir,ierr
  real*8 ::  px,py,pz, por, tor
  
  open(60, file="perm_field.in", action="read", status="old")
  read(60,*)
     
  call VecGetArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%tor,tor_p,ierr)
  call VecGetArrayF90(grid%porosity,por_p,ierr)

  do 
    call fiReadFlotranString(60, string, ierr)
    if (ierr /= 0) exit

    call fiReadInt(string,ir,ierr) 
    call fiDefaultMsg('na',ierr)
    ir=ir-1


    call fiReadDouble(string,px,ierr)
    call fiDefaultMsg('perm_x',ierr)

    call fiReadDouble(string,py,ierr)
    call fiDefaultMsg('perm_y',ierr)
   
    call fiReadDouble(string,pz,ierr)
    call fiDefaultMsg('perm_z',ierr)

    call fiReadDouble(string,por,ierr)
    call fiDefaultMsg('por',ierr)

    call fiReadDouble(string,tor,ierr)
    call fiDefaultMsg('tor',ierr)



!  read(60,*)ir,px,py,pz,por,tor
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
  close(60)
 
 ! select case(card)
 ! case('LOGO' 
     
  call VecRestoreArrayF90(grid%perm_xx, perm_xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_yy, perm_yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%perm_zz, perm_zz_p, ierr); CHKERRQ(ierr)              
  call VecRestoreArrayF90(grid%tor,tor_p,ierr)
  call VecRestoreArrayF90(grid%porosity,por_p,ierr)

end  subroutine Read_perm_field
  
  
subroutine Read_init_field(grid)

  use fileio_module

  implicit none

!#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
!  character(len=MAXSTRINGLENGTH) :: string 

  PetscScalar, pointer :: xx_p(:), iphase_p(:)
  integer iln,na,nx,ny,nz,ir,ierr,n
  real*8 :: x,y,z,phase,pl,pg,t,sl,sg,xl,xg,vf 
  
  open(60, file="pflow_init.dat", action="read", status="old")
  read(60,*)
  read(60,*)
  read(60,*)
     
  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%iphas, iphase_p,ierr)

  do n=1, grid%nmax 
    ir=n-1
  
    if (grid%ny==1) then 
     read(60,'(1p11e12.4)')x,z,phase,pl,pg,t,sl,sg,xl,xg,vf
 !    print *, x,z,phase,pl,pg,t,sl,sg,xl,xg,vf
    else
      read(60,'(1p10e12.4)')x,y,z,phase,pg,t,sg,xl,xg,vf
!    print *,x,y,z,phase,pg,t,sg,xl,xg,vf
    endif 
    
    do iln=1, grid%nlmax
      na = grid%nL2A(iln)
    
      if (na == ir) then
        nz = na/grid%nxy + 1
        ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
        nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
        iphase_p(iln)= phase
        if ((dabs(phase-1.D0)<1e-1) .or. (dabs(phase-2.D0)<1e-1)) then  
          xx_p(1+(iln-1)*grid%ndof) = pg 
          xx_p(2+(iln-1)*grid%ndof) = t
          xx_p(3+(iln-1)*grid%ndof) = xg
    !      print*,'ass:', iln,ir,na,xx_p(
        elseif (dabs(phase-3.D0)<1e-1) then
          xx_p(1+(iln-1)*grid%ndof) = pg 
          xx_p(2+(iln-1)*grid%ndof) = t
          xx_p(3+(iln-1)*grid%ndof) = sg
        else
          print *, "error in phase cond:",phase
          stop        
        endif
        exit   
      endif    
    enddo 
  enddo
  
  close(60)
    
  call VecRestoreArrayF90(grid%xx, xx_p, ierr)
  call VecRestoreArrayF90(grid%iphas, iphase_p,ierr)
  !propo0to034
end subroutine Read_init_field
  
 
integer function nxyz2na(grid,nx,ny,nz)

  implicit none

  type(pflowGrid), intent(inout) :: grid
  integer nx,ny,nz
  
  nxyz2na = nx -1 + (ny-1)*grid%nx + (nz-1)* grid%nxy

end function nxyz2na
   
subroutine Read_Geom_field(grid)

  use fileio_module

  implicit none

#include "definitions.h"

  type(pflowGrid), intent(inout) :: grid
  character(len=MAXSTRINGLENGTH) :: string 

 
  PetscScalar, pointer :: perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),por_p(:), &
                          tor_p(:), volume_p(:)
  integer iln,na,nx,ny,nz,ir,ierr
  real*8 ::  xc,yc,zc,vc
  real*8 ::  px,py,pz, por, tor
  real*8 ::  area, dist1, dist2, grav_ang
! real*8 ::  x,y,z
  character*4 card, word
  ! character(len=MAXSTRINGLENGTH) :: string 

  integer nx1, nx2, ny1, ny2, nz1,nz2, ina1, ina2,nnc, na1,na2, ncna
  
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
    call fiCharsToUpper(word,len_trim(word))
    call fiReadCard(word,card,ierr)

    if (grid%myrank == 0) print *, card

    select case(card)

      case('GRID')  
        do 
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('GRID',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/') exit

          call fiReadInt(string,ir,ierr) 
          call fiDefaultMsg('na',ierr)
          ir=ir-1
     
          call fiReadDouble(string,xc,ierr)
          call fiDefaultMsg('x',ierr)
   
          call fiReadDouble(string,yc,ierr)
          call fiDefaultMsg('y',ierr)

          call fiReadDouble(string,zc,ierr)
          call fiDefaultMsg('z',ierr)
   
          call fiReadDouble(string,vc,ierr)
          call fiDefaultMsg('vol',ierr)
   
          call fiReadDouble(string,px,ierr)
          call fiDefaultMsg('perm_x',ierr)

          call fiReadDouble(string,py,ierr)
          call fiDefaultMsg('perm_y',ierr)
   
          call fiReadDouble(string,pz,ierr)
          call fiDefaultMsg('perm_z',ierr)

          call fiReadDouble(string,por,ierr)
          call fiDefaultMsg('por',ierr)

          call fiReadDouble(string,tor,ierr)
          call fiDefaultMsg('tor',ierr)

          grid%x(ir+1)=xc
          grid%y(ir+1)=yc
          grid%z(ir+1)=zc

!          read(60,*)ir,px,py,pz,por,tor
          print *, 'Read geom 0:',ir,xc,yc,zc,vc, px,py, pz,por,tor
 
          do iln=1, grid%nlmax
            na = grid%nL2A(iln)
            if (na == ir) then
              nz = na/grid%nxy + 1
              ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
              nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
              volume_p(iln)=vc
              perm_xx_p(iln)= px
              perm_yy_p(iln)= py
              perm_zz_p(iln)= pz
              if (por>=0.D0 .AND. por <= 1.D0) por_p(iln)=por
              if (tor>=0.D0 .and. tor <= 1.D0) tor_p(iln)=tor
              print *, 'Read geom 1:',na,xc,yc,zc,vc, px,py, pz,por,tor
              exit
            endif
          enddo 
        enddo
        print *, grid%x
   
   ! then need to compute connection
   
      case('CONN')  
        nnc=0
        do 
     
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('CONN',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/') exit
   
          call fiReadInt(string,ncna,ierr) 
          call fiDefaultMsg('na1',ierr)

          call fiReadInt(string,na1,ierr) 
          call fiDefaultMsg('na1',ierr)
          na1=na1-1
     
          call fiReadInt(string,na2,ierr) 
          call fiDefaultMsg('na2',ierr)
          na2=na2-1

          call fiReadDouble(string,dist1,ierr)
          call fiDefaultMsg('dist 1',ierr)

          call fiReadDouble(string,dist2,ierr)
          call fiDefaultMsg('dist 2',ierr)
     
          call fiReadDouble(string,area,ierr)
          call fiDefaultMsg('Area',ierr)
   
          call fiReadDouble(string,grav_ang,ierr)
          call fiDefaultMsg('cos(B)',ierr)
    
   
          do iln=1, grid%nlmax
            !iln= grid%nG2L(ign) 
            ir=-1
            if (iln>0) ir = grid%nL2A(iln)
                     
            if (na1 ==ir) then
              nnc = nnc +1 
              grid%nd1(nnc)=grid%nL2G(iln)
              grid% dist1(nnc) = dist1
              grid% dist2(nnc) = dist2 
              if (area>0D0) grid%area(nnc)= area
              grid%delz(nnc) = grid%z(na1)- grid%z(na2)
              grid%grav_ang(nnc)= grav_ang
      
           ! then need to locate na2 and determine iperm
              nz1 = na1/grid%nxy + 1
              ny1 = (na1 - (nz1-1)*grid%nxy)/grid%nx + 1
              nx1 = na1 + 1 - (ny1-1)*grid%nx - (nz1-1)*grid%nxy
              ina2=0
          ! z+1
       
              nx=nx1; ny=ny1; nz=nz1+1
              ir = nxyz2na(grid,nx,ny,nz)
              if (na2 ==ir) then
                grid%nd2(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
                                (nz-1-grid%ngzs)*grid%ngxy
                ina2=1
                grid%iperm1(nnc) = 3
                grid%iperm2(nnc) = 3
              endif
       
              nx=nx1
              ny=ny1
              nz=nz1-1
              ir = nxyz2na(grid,nx,ny,nz)
              if (na2 ==ir) then
                grid%nd2(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
                                (nz-1-grid%ngzs)*grid%ngxy
                ina2=1
                grid%iperm1(nnc) = 3
                grid%iperm2(nnc) = 3
              endif
           
              nx=nx1+1
              ny=ny1
              nz=nz1
              ir = nxyz2na(grid,nx,ny,nz)
              if (na2 ==ir) then
                grid%nd2(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
                                (nz-1-grid%ngzs)*grid%ngxy
                ina2=1
                grid%iperm1(nnc) = 1
                grid%iperm2(nnc) = 1
              endif
                    
              nx=nx1-1
              ny=ny1
              nz=nz1
              ir = nxyz2na(grid,nx,ny,nz)
              if (na2 ==ir) then
                grid%nd2(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
                                (nz-1-grid%ngzs)*grid%ngxy
                ina2=1
                grid%iperm1(nnc) = 1
                grid%iperm2(nnc) = 1
              endif
                         
            !  nx=nx1
            !  ny=ny1+1
            !   nz=nz1
            ! ir = nxyz2na(grid,nx,ny,nz)
            ! if (na2 ==ir) then
            !   grid%nd2(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
            !                   (nz-1-grid%ngzs)*grid%ngxy
            !   ina2=1
            !   grid%iperm1(nnc) = 2
            !   grid%iperm2(nnc) = 2
            ! endif

           !  nx=nx1
           !  ny=ny1-1
           !  nz=nz1
           !  ir = nxyz2na(grid,nx,ny,nz)
           !  if (na2 ==ir) then
           !    grid%nd2(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
           !    (nz-1-grid%ngzs)*grid%ngxy
           !    ina2=1
           !    grid%iperm1(nnc) = 2
           !    grid%iperm2(nnc) = 2
           !  endif
    
    
           !  if (ina2 ==0) then  !BC
               print *, 'Read conn:',ncna,nnc, na1,na2, grid%nd1(nnc), &
                         grid%nd2(nnc),area,dist1, dist2,grav_ang, &
                         grid%iperm1(nnc)
               exit
            elseif (na2==ir) then      
              nnc = nnc +1 
              grid%nd2(nnc)=grid%nL2G(iln)
              grid% dist1(nnc) = dist1
              grid% dist2(nnc) = dist2 
              if (area>0D0) grid%area(nnc)= area
              grid%delz(nnc) = grid%z(na1)- grid%z(na2)
              grid%grav_ang(nnc)= grav_ang
      
              ! then need to locate na2 and determine iperm
              nz2 = na2/grid%nxy + 1
              ny2 = (na2 - (nz2-1)*grid%nxy)/grid%nx + 1
              nx2 = na2 + 1 - (ny2-1)*grid%nx - (nz2-1)*grid%nxy
              ina1=0
              ! z+1
       
              nx=nx2
              ny=ny2
              nz=nz2+1
              ir = nxyz2na(grid,nx,ny,nz)
              if (na1 ==ir) then
                grid%nd1(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
                                (nz-1-grid%ngzs)*grid%ngxy
                ina1=1
                grid%iperm1(nnc) = 3
                grid%iperm2(nnc) = 3
              endif
       
              nx=nx2
              ny=ny2
              nz=nz2-1
              ir = nxyz2na(grid,nx,ny,nz)
              if (na1 ==ir) then
                grid%nd1(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
                                (nz-1-grid%ngzs)*grid%ngxy
                ina1=1
                grid%iperm1(nnc) = 3
                grid%iperm2(nnc) = 3
              endif
           
              nx=nx2+1
              ny=ny2
              nz=nz2
              ir = nxyz2na(grid,nx,ny,nz)
              if (na2 ==ir) then
                grid%nd1(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
                                (nz-1-grid%ngzs)*grid%ngxy
                ina1=1
                grid%iperm1(nnc) = 1
                grid%iperm2(nnc) = 1
              endif
                    
              nx=nx2-1
              ny=ny2
              nz=nz2
              ir = nxyz2na(grid,nx,ny,nz)
              if (na1 ==ir) then
                grid%nd1(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
                                (nz-1-grid%ngzs)*grid%ngxy
                ina1=1
                grid%iperm1(nnc) = 1
                grid%iperm2(nnc) = 1
              endif
              print *, 'Read conn: error !!'                     
            !  nx=nx2
            !  ny=ny2+1
            !  nz=nz2
            !  ir = nxyz2na(grid,nx,ny,nz)
            !  if (na1 ==ir) then
            !    grid%nd1(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
            !                    (nz-1-grid%ngzs)*grid%ngxy
            !    ina1=1
            !    grid%iperm1(nnc) = 2
            !    grid%iperm2(nnc) =2
            !  endif

            !  nx=nx2
            !  ny=ny2-1
            !  nz=nz2
            !  ir = nxyz2na(grid,nx,ny,nz)
            !  if (na1 ==ir) then
            !    grid%nd1(nnc) = nx- grid%ngxs + (ny-1-grid%ngys)*grid%ngx + &
            !                    (nz-1-grid%ngzs)*grid%ngxy
            !    ina1=1
            !    grid%iperm1(nnc) = 2 
            !    grid%iperm2(nnc) =2
            !  endif
    
    
            !  if (ina2 ==0) then  !BC
                                                                 
              exit
 
            endif
          enddo 
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
   
   PetscScalar, pointer :: xx_p(:)
   integer nc,m,ibc,ierr

   call VecGetArrayF90(grid%xx,xx_p,ierr)  
     do nc=1,grid%nconnbc
         m=grid%mblkbc(nc)
        ibc = grid%ibconn(nc)
        if (grid%ibndtyp(ibc) ==3 ) grid%xxbc(1,nc)=xx_p(1+(m-1)*grid%ndof)
   enddo
   call VecRestoreArrayF90(grid%xx, xx_p, ierr)
 end subroutine Boundary_adjustment

end module readfield
