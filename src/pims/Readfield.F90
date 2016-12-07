Module readfield

! Note:  subroutine to read in permeability field according to node natural ordering
!        if the nature order is not appeared in the file, it will take the old value   
!        Repeated index will be ignored (only take the first one)

#include "include/finclude/petscsnes.h"
 use petscsnes
 use pflow_gridtype_module
 ! use pflow_var_module

private 


  public Read_perm_field, Read_init_field,Boundary_adjustment
  
 contains
  
 subroutine Read_perm_field(grid, locpat)
  use fileio_module
  implicit none
  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info), intent(inout) :: locpat
#include "definitions.h"
  character(len=MAXSTRINGLENGTH) :: string 

  PetscScalar, pointer :: perm_xx_p(:),perm_yy_p(:),perm_zz_p(:),por_p(:),tor_p(:)
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
     if(ierr /= 0) exit

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
 
         do iln=1, locpat%nlmax
            na = locpat%nL2A(iln)
    
      if(na == ir) then
        nz = na/grid%nxy + 1
        ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
        nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
        perm_xx_p(iln)= px
		perm_yy_p(iln)= py
		perm_zz_p(iln)= pz
		if(por>=0.D0 .AND. por <=1.D0) por_p(iln)=por
		if(tor>=0.D0 .and. tor <=1.D0) tor_p(iln)=tor
			        
         
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
  
  
  subroutine Read_init_field(grid, locpat)
  use fileio_module
  implicit none
  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info), intent(inout) :: locpat
!#include "definitions.h"
!  character(len=MAXSTRINGLENGTH) :: string 

  PetscScalar, pointer :: xx_p(:)
  integer iln,na,nx,ny,nz,ir,ierr,n
  real*8 :: x,y,z,pl,pg,t,sl,sg,xl,xg,vf 
  
  open(60, file="pflow_init.dat", action="read", status="old")
  read(60,*)
  read(60,*)
  read(60,*)
   	
  call VecGetArrayF90(grid%xx, xx_p, ierr); CHKERRQ(ierr)

  do n=1, grid%nmax 
  ir=n-1
  
  if(grid%ny==1)then 
   read(60,'(1p11e12.4)')x,z,pl,pg,t,sl,sg,xl,xg,vf
 !  print *, x,z,phase,pl,pg,t,sl,sg,xl,xg,vf
  else
    read(60,'(1p10e12.4)')x,y,z,pg,t,sg
!	print *,x,y,z,phase,pg,t,sg,xl,xg,vf
  endif 
		
   do iln=1, locpat%nlmax
       na = locpat%nL2A(iln)
    
      if(na == ir) then
        nz = na/grid%nxy + 1
        ny = (na - (nz-1)*grid%nxy)/grid%nx + 1
        nx = na + 1 - (ny-1)*grid%nx - (nz-1)*grid%nxy
      		xx_p(1+(iln-1)*grid%ndof) = pg 
			xx_p(2+(iln-1)*grid%ndof) = t
			xx_p(3+(iln-1)*grid%ndof) = sg
	   exit
	 endif 	  
    enddo 
  enddo
  
  close(60)
    
  call VecRestoreArrayF90(grid%xx, xx_p, ierr)

  end subroutine Read_init_field
  
  
end Module readfield
