
module mphase_field_module
  
  implicit none
 
  private 
#include "definitions.h"
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
  
  type, public :: mphase_field_type
    Vec :: var_loc
    Vec :: pressure
    Vec :: temp
    Vec :: xmol
    Vec :: sat
    Vec :: conc
    Vec :: phis
    
    PetscReal, pointer :: xxphi_co2(:)
    PetscReal, pointer :: dden_co2(:)
    PetscReal, pointer :: xxphi_co2_bc(:)
    
  end type mphase_field_type
  
  type(mphase_field_type), pointer, public :: mphase_field
  
end module mphase_field_module  

module mphase_option_module
  
  implicit none
  
  private 
  
#include "definitions.h"
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
  
  type, public :: mphase_option_type
    PetscInt :: iphch
    PetscInt :: jh2o = 1
    PetscInt :: jco2 = 2
    PetscInt :: jgas = 3
    
    PetscInt :: idt_switch    

    PetscReal, pointer :: rate(:)
    PetscReal, pointer :: area_var(:)
        
!   solid reaction rate
    PetscInt :: ityprxn
    PetscReal :: rk=0.d0, phis0, areas0, pwrsrf, vbars, ceq, delHs, delEs, wfmts
    PetscReal ::qu_kin, yh2o_in_co2=0.D0
    
    PetscReal, pointer :: delx(:,:)
  end type mphase_option_type
  
  type(mphase_option_type), pointer, public :: mphase_option
  
end module mphase_option_module  

module translator_mph_module
 
  use mphase_field_module  
  use mphase_option_module 
  
  implicit none  

  private 

#include "definitions.h"
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
  
   
! **3 phase condition*************************************************
! phase              Primary Variables      index
!  e                p, T, X(e,a), X(e,c)      1
!  g                p, T, X(g,a), X(g,c)      2 
!  l                p, T, X(l,a), X(l,c)      4
!  eg               p, T,  S(g),  X(g,c)      3 
!  el               p, T,  S(l),  X(l,c)      5
!  lg               p, T,  S(g),  X(g,c)      6
!  egl              p, T,  S(g),  S(l)        7
!**********************************************************************


! phase index 1.e; 2. l; 3. g
! within each phase component index : 1. H2O; 2. CO2; 3. Air

    
  public  pri_var_trans_mph_ninc,pri_var_trans_mph_winc ,translator_mph_step_maxchange, &
          translator_mph_get_output,translator_check_phase_cond,translator_mphase_massbal, &
          Translator_MPhase_Switching
     
     
  PetscReal, private, parameter :: fmwh2o = 18.0153D0, fmwa = 28.96D0, &
                               fmwco2 = 44.0098D0, fmwnacl = 58.44277D0
  PetscReal, private, parameter :: eps=5D-7 , formeps=1D-4
  PetscReal, private, parameter ::yh2o_in_co2=1D-2   
  PetscReal, private, parameter :: rgasj   = 8.3143    ![J/K/mol]

contains

! subroutines to calculate the properties of mixture  
! will:: call other EOS mod to obtain pure fluid properties
!        apply mixing rules


subroutine translator_mphase_massbal(realization)
 
  use Realization_module
  use Grid_module
  use Option_module
  use Field_module
  
  implicit none
  
  type(realization_type) :: realization
 
  PetscInt, save :: icall
  PetscErrorCode :: ierr
  PetscInt :: nc,np,n2p,n2p0
  PetscReal :: x,y,z,nzm,nzm0, nxc,nxc0,c0, c00,nyc,nyc0,nzc,nzc0,nsm,nsm0,sm 
  PetscReal :: nxm,nxm0
  PetscInt :: index, size_var_node
  PetscInt :: dof_offset
  PetscInt :: local_id, ghosted_id
     
  PetscReal, pointer :: var_loc_p(:), porosity_loc_p(:), volume_p(:)
                           
  PetscReal, pointer :: iphase_loc_p(:)
  
  PetscReal ::  pvol,sum
  PetscReal, pointer ::  den(:),sat(:),xmol(:)
  PetscReal :: sat_avg, sat_max, sat_min, sat_var
  PetscReal :: sat_avg0, sat_max0, sat_min0, sat_var0
  PetscReal :: tot(0:realization%option%nspec,0:realization%option%nphase), &
            tot0(0:realization%option%nspec,0:realization%option%nphase)  
  data icall/0/
    
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  grid => realization%grid
  option => realization%option
  field => realization%field  

  call VecGetArrayF90(mphase_field%var_loc,var_loc_p,ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)  
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr)
 
  size_var_node=(option%ndof+1)*(2+7*option%nphase+2*option%nphase*option%nspec)
  tot=0.D0
  n2p=0
  nxc=0.; nyc=0.; nzc=0.D0; nzm=0 !grid%z(grid%nmax); 
  sm=0.; nsm=nzm; c0=0.D0; nxm =0
  sat_avg =0.D0; sat_max =0.D0; sat_min =1.D0; sat_var =0.D0
 ! pvol_tot =0.D0

  do local_id = 1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    dof_offset=(local_id-1)* option%ndof
    index=(ghosted_id-1)*size_var_node
    den=>var_loc_p(index+3+option%nphase: index+2+2*option%nphase)
    sat=>var_loc_p(index+2+1:index+2+option%nphase)
    xmol=>var_loc_p(index+2+7*option%nphase+1:index+2+7*option%nphase +&
                    option%nphase*option%nspec)    
   
    pvol=volume_p(local_id)*porosity_loc_p(ghosted_id)         
    if(dabs(iphase_loc_p(ghosted_id)- 3.D0)<.25D0)then
      n2p=n2p+1
      x=grid%x(ghosted_id)
      y=grid%y(ghosted_id)
      z=grid%z(ghosted_id)
      
      if(z>nzm) nzm=z
      if(x>nxm) nxm=x
      if(sm<sat(2))then
      !print *, n,grid%nL2A(n)+1,x,y,z,sat(2)
        sm=sat(2);nsm=z
      endif    
       sat_avg = sat_avg + sat(2) !* pvol
       sat_var = sat_var + sat(2) * sat(2)! * pvol *pvol
      ! pvol_tot =pvol_tot + pvol
       if(sat_max < sat(2))sat_max = sat(2)
       if(sat_min > sat(2))sat_min = sat(2)      

     endif

    do nc =1,option%nspec
      do np=1,option%nphase
        sum= sat(np)* xmol((np-1)*option%nspec +nc)*den(np)
        tot(nc,np)= pvol*sum + tot(nc,np)
       !tot(0,np)=tot(0,np)+tot(nc,np)
       !tot(nc,0)=tot(nc,0)+tot(nc,np)
      enddo
    enddo
    nullify(sat, den,xmol) 
!   print *,nzc,c0
  enddo
 !  call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
  call VecRestoreArrayF90(mphase_field%var_loc,var_loc_p,ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr)
 
 
  if(option%commsize >1)then
    call MPI_REDUCE(n2p, n2p0,ONE_INTEGER, MPI_INTEGER,MPI_SUM,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(nzm, nzm0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_MAX,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(nxm, nxm0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_MAX,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
!    call MPI_REDUCE(pvol_tot, pvol_tot0,1, MPI_DOUBLE_PRECISION,MPI_SUM,0, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(sat_avg, sat_avg0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_SUM,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(sat_var, sat_var0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_SUM,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(sat_min, sat_min0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_MIN,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
    call MPI_REDUCE(sat_max, sat_max0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_MAX,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
   
      
    do nc = 0,option%nspec
      do np = 0,option%nphase
        call MPI_REDUCE(tot(nc,np), tot0(nc,np),ONE_INTEGER,&
            MPI_DOUBLE_PRECISION,MPI_SUM,ZERO_INTEGER, PETSC_COMM_WORLD,ierr)
   
!       call MPI_BCAST(tot0,(option%nphase+1)*(option%nspec+1),&
!            MPI_DOUBLE_PRECISION, ZERO_INTEGER,PETSC_COMM_WORLD,ierr)
      enddo
    enddo
    if(option%myrank==0) then
      tot = tot0; n2p=n2p0;nxc=nxc0;nyc=nyc0;nzc=nzc0;nzm=nzm0;nsm=nsm0;c0=c00
      nxm =nxm0; sat_avg=sat_avg0; sat_var =sat_var0; sat_min=sat_min0; sat_max=sat_max0
    endif   
  endif 
  
  
  if(option%myrank==0)then
    if(c0<1D-6) c0=1.D0
    nxc=nxc/c0; nyc=nyc/c0;nzc=nzc/c0
    if(n2p>0)sat_avg = sat_avg/n2p
  !  sat_var = sat_var / n2p  -  sat_avg*sat_avg
  
  
    write(*,'(" Total CO2: t= ",1pe12.4," dt= ",1pe12.4," liq:",1pe13.6,&
   &" gas:",1pe13.6, " tot:", 1p2e13.6, " [kmol]",1p3e13.6)') &
    option%time/realization%output_option%tconv,option%dt/realization%output_option%tconv,tot(2,1),tot(2,2),tot(2,1)+tot(2,2) !,nzc,nzm,nsm
! & option%t/realization%output_option%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2) !,nzc,nzm,nsm
    if (icall==0) then
      open(unit=13,file='massbal.dat',status='unknown')
      write(13,'(''# time   dt   totl   totg   tot n2p'')')
      icall = 1
    endif
!   write(13,'(" Total CO2: t=",1pe13.6," liq:",1pe13.6,&
! &  " gas:",1pe13.6," tot:",1p2e13.6," [kmol]")')&
! & option%t/realization%output_option%tconv,tot(2,1),tot(2,2),tot(2,0),tot(2,1)+tot(2,2)
    write(13,'(1p19e12.4)') option%time/realization%output_option%tconv,option%dt/realization%output_option%tconv,&
    tot(2,1),tot(2,2),tot(2,1)+tot(2,2),real(n2p), nzm, nxm,&
    sat_avg, sat_min, sat_max, sat_var 
  endif    
  
end subroutine translator_mphase_massbal

PetscInt function translator_check_phase_cond(iphase, var_node,num_phase,num_spec)
 
  implicit none
   
  PetscInt :: iphase, num_phase, num_spec
  PetscReal, target:: var_node(:)
    
  PetscInt :: ibase,succ,np,nc
  PetscReal, pointer :: t,p,satu(:),den(:), avgmw(:),h(:),u(:),pc(:),&
                      kvr(:),xmol(:),diff(:)
      
  PetscReal sum     
    
  ibase=1;               t=>var_node(ibase)
  ibase=ibase+1;         p=>var_node(ibase)
  ibase=ibase+1;         satu=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; den=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; avgmw=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; h=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; u=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; pc=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; kvr=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; xmol=>var_node(ibase:ibase+num_phase*num_spec-1)
  ibase=ibase+num_phase*num_spec; diff=>var_node(ibase:ibase+num_phase*num_spec-1)
  
  succ=1
  if(iphase ==3)then
    do np =1, num_phase
      if(satu(np)>1D0-eps .or.satu(np)<eps)then
        succ=-1  
        print *, 'phase=',iphase,satu(1:2)
      endif
    enddo
  endif
 
  if(iphase ==1 .or. iphase==2)then
    do np =1, num_phase
      sum=0D0
      do nc = 1, num_spec
        sum =sum + xmol((np-1)*num_spec+nc)
      enddo 
      if(sum > 1.D0+eps)then
        succ=-1 
        print *,'phase=',iphase,sum,xmol(1:4)
      endif 
    enddo
  endif
  nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)     
  translator_check_phase_cond = succ

end function translator_check_phase_cond



subroutine translator_mph_get_output(realization)
 
  use Realization_module
  use Option_module
  use Grid_module, only : grid_type
  use Field_module 
  
  implicit none
  
  type(realization_type) :: realization

  PetscErrorCode :: ierr
      
  PetscReal, pointer :: t_p(:),p_p(:),c_p(:),s_p(:),cc_p(:),var_loc_p(:)
  PetscInt :: index_var_begin ,jn, size_var_node
! PetscReal, pointer :: p,t,satu(:),xmol(:)
  PetscInt :: local_id, ghosted_id

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  option => realization%option
  grid => realization%grid
  field => realization%field
      
  call VecGetArrayF90(mphase_field%var_loc, var_loc_p, ierr)
  call VecGetArrayF90(mphase_field%pressure, p_p, ierr)
  call VecGetArrayF90(mphase_field%temp, t_p, ierr)
  call VecGetArrayF90(mphase_field%xmol, c_p, ierr)
  call VecGetArrayF90(mphase_field%sat, s_p, ierr)
  call VecGetArrayF90(mphase_field%conc, cc_p, ierr)
  !print *,' translator_mph_get_output gotten pointers'
  
  size_var_node=(option%ndof+1)*(2+7*option%nphase +2*option%nphase*option%nspec)
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    index_var_begin= (ghosted_id-1) * size_var_node
    jn = 1 + (local_id-1)*option%nphase 
    
  p_p(jn) = var_loc_p(index_var_begin + 2) - var_loc_p(index_var_begin+5*option%nphase+3)
  p_p(jn+1) = var_loc_p(index_var_begin + 2) - var_loc_p(index_var_begin+5*option%nphase+4)
   
    t_p(local_id)=var_loc_p(index_var_begin + 1)

    c_p(jn)=var_loc_p(index_var_begin+7*option%nphase+4)
  c_p(jn+1)= var_loc_p(index_var_begin+7*option%nphase+6)
    cc_p(local_id)=c_p(jn+1)
  
  s_p(jn)=var_loc_p(index_var_begin + 3) 
    s_p(jn+1)=var_loc_p(index_var_begin + 4) 
 enddo
 
  call VecRestoreArrayF90(mphase_field%var_loc, var_loc_p, ierr)
  call VecRestoreArrayF90(mphase_field%pressure, p_p, ierr)
  call VecRestoreArrayF90(mphase_field%temp, t_p, ierr)
  call VecRestoreArrayF90(mphase_field%xmol, c_p, ierr)
  call VecRestoreArrayF90(mphase_field%sat, s_p, ierr)
  call VecRestoreArrayF90(mphase_field%conc, cc_p, ierr)
 
! work only for 2 phases
end subroutine translator_mph_get_output


subroutine translator_mph_step_maxchange(realization)

  use Realization_module
  use Option_module
  use Field_module
  use Grid_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(grid_type), pointer :: grid
  
  PetscReal, pointer :: xx_p(:), yy_p(:), iphase_loc_p(:),var_loc_p(:),iphase_old_loc_p(:)
  PetscReal :: comp1,comp,cmp  
! PetscReal :: dsm,dcm  
  PetscReal :: dsm0,dcm0  
  PetscInt :: local_id, ghosted_id, dof_offset
! PetscInt :: j
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  grid => realization%grid
  
   call VecWAXPY(field%dxx,-1.d0,field%xx,field%yy,ierr)
    call VecStrideNorm(field%dxx,0,NORM_INFINITY,option%dpmax,ierr)
    call VecStrideNorm(field%dxx,1,NORM_INFINITY,option%dtmpmax,ierr)

  call VecGetArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  call VecGetArrayF90(field%iphas_old_loc, iphase_old_loc_p,ierr)
  call VecGetArrayF90(mphase_field%var_loc, var_loc_p, ierr)
  
  comp=0.D0;comp1=0.D0
  do local_id=1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    dof_offset=(local_id-1)*option%ndof 
    if(int(iphase_loc_p(ghosted_id)) == int(iphase_old_loc_p(ghosted_id)))then
     cmp=dabs(xx_p(dof_offset+3)-yy_p(dof_offset+3))
     if(int(iphase_loc_p(ghosted_id))==1 .or.int(iphase_loc_p(ghosted_id))==2)then
       if(comp<cmp) comp=cmp
     
    endif   
       if(int(iphase_loc_p(ghosted_id))==3)then
       if(comp1<cmp) comp1=cmp
      
    endif   
    else
!  print *,'phase changed', n, iphase_loc_p(n), iphase_old_p(n)

   endif
  enddo
  !call PETSCBarrier(PETSC_NULL_OBJECT,ierr)
  call VecRestoreArrayF90(field%xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  call VecRestoreArrayF90(field%iphas_old_loc, iphase_old_loc_p,ierr)
  call VecRestoreArrayF90(mphase_field%var_loc, var_loc_p, ierr)
 
  
  if(option%commsize >1)then
    call MPI_ALLREDUCE(comp1, dsm0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    !call MPI_BCAST(dsm0,ONE_INTEGER, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(comp, dcm0,ONE_INTEGER, MPI_DOUBLE_PRECISION,MPI_MAX, PETSC_COMM_WORLD,ierr)
    !call MPI_BCAST(dcm0,ONE_INTEGER, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD,ierr)
    comp1 = dsm0
    comp = dcm0
  endif 

   option%dsmax=comp1
   option%dcmax=comp
!   print *, 'max change',grid%dpmax,option%dtmpmax,grid%dsmax,grid%dcmax

end subroutine translator_mph_step_maxchange


subroutine Translator_MPhase_Switching(xx,realization,icri,ichange)
  
  use Realization_module
  use Option_module
  use Field_module
  use Grid_module
  
  use water_eos_module
  use gas_eos_module  
  use co2eos_module
  use span_wagner_module

  implicit none
  
  type(realization_type) :: realization
  
  Vec, intent(in) :: xx
  PetscInt :: icri,ichange 

  PetscReal, pointer :: xx_p(:), yy_p(:),iphase_loc_p(:)
  PetscInt :: ipr
  PetscInt :: iipha 
  PetscErrorCode :: ierr
! PetscInt :: index,i
  PetscReal :: p2,p,tmp,t
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  PetscReal :: ug,xphi,henry,sat_pressure
  PetscReal :: xmol(realization%option%nphase*realization%option%nspec),satu(realization%option%nphase)
! PetscReal :: xla,co2_poyn
  PetscInt :: local_id, ghosted_id, dof_offset
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  
! mphase code need assemble 
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p,ierr)
 
   
  ichange = 0   
  do local_id = 1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    ipr=0
    dof_offset=(local_id-1)* option%ndof
    iipha=iphase_loc_p(ghosted_id)
    p = xx_p(dof_offset+1)
    t= xx_p(dof_offset+2)
    select case(iipha) 
    case(1) 
      xmol(2)= xx_p(dof_offset+3)
      xmol(1)=1.D0 - xmol(2)
      satu(1)=1.D0; satu(2)=0.D0
    case(2)  
      xmol(4)= xx_p(dof_offset+3)
      xmol(3)=1.D0 - xmol(4)
      satu(1)=eps; satu(2)=1.D0
    case(3) 
        satu(2)= xx_p(dof_offset+3) 
      satu(1)= 1.D0- satu(2)
      xmol(3)= yh2o_in_co2; xmol(4)=1.D0-xmol(3)
    end select
    
    p2=p!*xmol(4)
!    p2=p*xmol(4)
    
    if(p2>=5d4)then
      !call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
     
      call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
      dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
      dg = dg / fmwco2
      fg = fg * 1.D6 
      hg = hg * fmwco2
    else      
      call ideal_gaseos_noderiv(p2,t,option%scale,dg,hg,ug)
      call visco2(t,dg*fmwco2,visg)
      fg = p2
    endif

    xphi = fg/p2
!   xphi = 1.d0 ! test-pcl

!   call Henry_CO2_noderiv(xla,tmp,t,p*xmol(4),xphi,henry,co2_poyn)
!  call Henry_CO2_noderiv(xla,tmp,t,p*xmol(4),xphi,henry,co2_poyn)
       
    call PSAT(t, sat_pressure, ierr)
    sat_pressure =sat_pressure /1D5
    call Henry_duan_sun(t, p2 *1D-5, henry,xphi,option%m_nacl,option%m_nacl,sat_pressure)
    
    henry= 1D0 / (fmwh2o*1D-3) / (henry*1D-5) /xphi 
    !note: henry = H/phi
    
!    print *,'translator_mixed: ',iipha,xla,tmp,t,p,xmol(4),p*xmol(4), &
!    xphi,henry,dg,fg,hg
     
  !   tmp = (p-sat_pressure)/(henry-sat_pressure)
 !    err= dabs(tmp-xmol(4))
!   tmp=xmol(4)
 ! enddo
 ! sat_pressure =sat_pressure * 1D5
      
 !print *, 'out 2 phase solver'
    select case(iipha)     
      case(1)
        xmol(4)=xmol(2)*henry/p   
   
 ! if(iphase /= 1)then
        xmol(3)=1.D0-xmol(4)
        if(xmol(3)<0.D0)xmol(3)=0.D0
  !if(xmol(3)<0.D0) xmol(3)=0.D0
      case(2)   
  
        xmol(2)= p*xmol(4)/henry
    !xmol(2)= 1.D0/(1.D0 + 1.D0 /(henry * xmol(4) * p* 1D-5 *xphi * grid%fmwh2o ))
        
        xmol(1)=1.D0-xmol(2)
        if(xmol(1)<0.D0) xmol(1)=0.D0
      case(3)
        tmp = sat_pressure*1D5 / p 
        xmol(2)=(1.D0-tmp)/(Henry/p -tmp)
        xmol(1)= 1.D0- xmol(2)
        xmol(3)=xmol(1) * tmp
        xmol(4)= 1.D0-xmol(3)            
    end select
   

    select case(icri)
      case(0)
        tmp = sat_pressure* 1D5 /p 
        if (xmol(3) >tmp*1.05 .and. iipha==2 )then
          write(*,'('' Gas -> 2ph '',2i8,1p10e12.4)') mphase_option%iphch,local_id,xx_p(dof_offset+1:dof_offset+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_loc_p(ghosted_id) = 3
          xx_p(dof_offset+3)=1D0-formeps
          ichange = 1 ;ipr=1
!       xx_p(n0+1)= yy_p(n0+1);xx_p(n0+2)= yy_p(n0+2)
        endif
       ! gas ->  2ph 
 
    
        tmp = (1.D0-tmp)/(henry/p -tmp) *henry/p
    ! print *,n, tmp, henry,sat_pressure,p
    !if (xx_p(n0+3) > 1.025D0  .and. iipha==1) then
        if (xmol(4) > 1.05D0 *tmp  .and. iipha==1) then
          write(*,'('' Liq -> 2ph '',2i8,1p10e12.4)') mphase_option%iphch,local_id,xx_p(dof_offset+1:dof_offset+3),xmol(4), tmp
      !write(IUNIT2,'('' Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_loc_p(ghosted_id) = 3
     
     !tmp= (xmol(4)-1.D0)/(xmol(4)/xmol(2))*den(1)/den(2)
    !if(tmp>eps) tmp=eps
     ! xx_p(n0+3)=max(tmp,formeps)
          xx_p(dof_offset+3)=formeps
!     xx_p(n0+1)= yy_p(n0+1);xx_p(n0+2)= yy_p(n0+2)
    !if(dabs(xmol(4)-1.D0) < eps) xx_p(n0+1)=xx_p(n0+1)* (1.D0-eps)
          ichange = 1;ipr=1
        endif
      
        if(satu(2)>1.0D0.and. iipha==3 )then
  !if(xx_p(n0+3)> 1.D0 .and. iipha==3 )then
          write(*,'('' 2ph -> Gas '',2i8,1p10e12.4)') mphase_option%iphch,local_id,xx_p(dof_offset+1:dof_offset+3)
       ! write(IUNIT2,'('' 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_loc_p(ghosted_id) = 2
          xx_p(dof_offset + 3) = xmol(4) + (1.D0 - xmol(4)) * formeps
!    xx_p(n0+1)= yy_p(n0+1);xx_p(n0+2)= yy_p(n0+2)  
          ichange =1    ;ipr=1  
        endif
    

   ! if(sat(2)<= -formeps .and. iipha==3 )then
        if(satu(2)<= 0.0D0 .and. iipha==3 )then
          write(*,'('' 2ph -> Liq '',2i8,1p10e12.4)') mphase_option%iphch,local_id,xx_p(dof_offset+1:dof_offset+3),satu(1),satu(2)
      ! write(IUNIT2,'('' 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_loc_p(ghosted_id) = 1 ! 2ph -> Liq
          ichange = 1;ipr=1
          tmp = xmol(2) * 0.99
          xx_p(dof_offset + 3)=tmp
!    xx_p(n0+1)= yy_p(n0+1);xx_p(n0+2)= yy_p(n0+2)
    !xx_p(n0+2) =  xx_p(n0+2)*(1.D0-eps)
        endif
     
    
!  if(ichange ==1) then
!    call 
      case(1)
        tmp = sat_pressure*1D5/p    
        if (xmol(3) >tmp .and. iipha==2 )then
          write(*,'('' Gas -> 2ph '',i8,1p10e12.4)') local_id,xx_p(dof_offset+1:dof_offset+3)
      !write(IUNIT2,'('' Gas -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_loc_p(ghosted_id) = 3
          xx_p(dof_offset+3) = 1D0 - eps
          ichange = 1 ;ipr=1
        !  xx_p(n0+2)= yy_p(n0+2)
        endif
       ! gas ->  2ph 
        tmp = (1.D0-tmp)/(henry/p -tmp)*henry/p
        if (xmol(4) > 1.025*tmp  .and. iipha==1) then
          write(*,'('' Liq -> 2ph '',i8,1p10e12.4)') local_id,xx_p(dof_offset+1:dof_offset+3),xmol(4), tmp
      !write(IUNIT2,'('' Liq -> 2ph '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_loc_p(ghosted_id) = 3
     
     !tmp= (xmol(4)-1.D0)/(xmol(4)/xmol(2))*den(1)/den(2)
    !if(tmp>eps) tmp=eps
     ! xx_p(n0+3)=max(tmp,formeps)
          xx_p(dof_offset+3) = eps
    !if(dabs(xmol(4)-1.D0) < eps) xx_p(n0+1)=xx_p(n0+1)* (1.D0-eps)
          ichange = 1;ipr=1
        endif
      
        if(satu(2)>1.0D0.and. iipha==3 )then
  !if(xx_p(n0+3)> 1.D0 .and. iipha==3 )then
          write(*,'('' 2ph -> Gas '',i8,1p10e12.4)') local_id,xx_p(dof_offset+1:dof_offset+3)
       ! write(IUNIT2,'('' 2ph -> Gas '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_loc_p(ghosted_id) = 2
          xx_p(dof_offset + 3) = xmol(4)  
          ichange =1    ;ipr=1  
        endif
    

   ! if(sat(2)<= -formeps .and. iipha==3 )then
        if(satu(2)<= 0.0D0 .and. iipha==3 )then
          write(*,'('' 2ph -> Liq '',i8,1p10e12.4)') local_id,xx_p(dof_offset+1:dof_offset+3),satu(1),satu(2)
      !write(IUNIT2,'('' 2ph -> Liq '',i8,1p10e12.4)') n,xx_p(n0+1:n0+3)
          iphase_loc_p(ghosted_id) = 1 ! 2ph -> Liq
          ichange = 1;ipr=1
          tmp = xmol(4) * 0.9999995
          xx_p(dof_offset + 3)=tmp
    !xx_p(n0+2) =  xx_p(n0+2)*(1.D0-eps)
        endif
    
    end select


    if(ipr ==1)then
!        iicap=int(icap_p(n))
!    iipha = int(iphase_loc_p(n))
!    dif(1)= grid%difaq
!        dif(2)= grid%cdiff(int(ithrm_p(n)))
!        i=ithrm_p(n) 
   
!  call pri_var_trans_ninc(xx_p((n-1)*option%ndof+1:n*option%ndof),iipha,&
 !       option%scale,option%nphase,option%nspec,&
 !       iicap, grid%sir(1:option%nphase,iicap),grid%lambda(iicap),&
 !       grid%alpha(iicap),grid%pckrm(iicap),grid%pcwmax(iicap),&
 !       grid%pcbetac(iicap),grid%pwrprm(iicap),dif,&
!    var_loc_p((n-1)*size_var_node+1:(n-1)*size_var_node+size_var_use),&
!    grid%itable,ierr)

   
 !    index_var_begin=(n-1)*size_var_node+1
  !   index_var_end = index_var_begin -1 + size_var_use
     
   !  p1 = 1 + (n-1)*option%ndof    
   !  call MPHASERes_ARCont(n, var_loc_p(index_var_begin: index_var_end),&
!    porosity_loc_p(n),volume_p(n),grid%dencpr(i), grid, Res, 0,ierr)

 !  print *,res,accum_p(p1:p1-1+option%ndof)
   !accum_p(p1:p1-1+option%ndof) =res
      ipr=0
    endif

  end do

  !print *,iphase_loc_p
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p,ierr)
  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(field%yy, yy_p, ierr); CHKERRQ(ierr)

end subroutine Translator_MPhase_Switching
  

subroutine pri_var_trans_mph_ninc_3_3(x,iphase,energyscale,num_phase,num_spec,&
                    ipckrtype,pckr_sir,pckr_lambda,pckr_alpha,&
                    pckr_m,pckr_pcmax,pckr_betac,pckr_pwr,dif,&
                    var_node,itable,ierr)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
  use water_eos_module
  use gas_eos_module  
  use pckr_module
  use co2eos_module
  use span_wagner_module

  implicit none
  
  PetscInt :: num_phase,num_spec
  PetscInt :: size_var_use
  PetscReal x(1:num_spec+1),energyscale
  PetscReal,target:: var_node(:)
  PetscInt :: iphase,itable
  PetscErrorCode :: ierr
  PetscInt :: ipckrtype !, ithrmtype
!   PetscInt :: num_pricomp
  
  PetscReal  :: pckr_sir(:),pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_betac,pckr_pwr
  PetscReal  :: dif(:)
     
  PetscReal, pointer :: t,p
  PetscReal, pointer:: den(:),h(:),u(:),avgmw(:),pc(:),kvr(:)
  PetscReal, pointer :: diff(:),xmol(:),satu(:)
  PetscInt :: ibase 
! PetscReal err

  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec 
  ibase=1;               t=>var_node(ibase)
  ibase=ibase+1;           p=>var_node(ibase)
  ibase=ibase+1;           satu=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; den=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; avgmw=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; h=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; u=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; pc=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; kvr=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; xmol=>var_node(ibase:ibase+num_phase*num_spec-1)
  ibase=ibase+num_phase*num_spec; diff=>var_node(ibase:ibase+num_phase*num_spec-1)

  select case(iphase)
  
    case(1) ! only water phase
     p = x(1)
     t = x(2)
     xmol(2)= x(4)
     xmol(2)= x(3)
     
    case(2) ! only Gas Phase
    
      
    case(4) ! only Supercritical CO2 
    case(3) ! water + gas 
    case(5) ! water + CO2
    case(6) ! gas + supercritical CO2 
    case(7) ! 3 phase
  end select
  
   nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)
   
end subroutine pri_var_trans_mph_ninc_3_3
  

! **2 phase condition**************************************************
! phase                             Primary Variables      index
!   e                p, T, X(e,c)                  1
!   g                p, T, X(g,a)                  2 
!   eg                              p, T, S(g)                    3
!**********************************************************************

subroutine pri_var_trans_mph_ninc_2_2(x,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg, dif, var_node,itable,m_nacl,ierr,xphi,dco2)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
  use water_eos_module
  use gas_eos_module  
  use pckr_module
  use co2eos_module
  use span_wagner_module


  implicit none
  
  PetscInt :: num_phase,num_spec, itable
  PetscErrorCode :: ierr
  PetscInt :: size_var_use 
  PetscReal x(1:num_spec+1),energyscale
  PetscReal, target:: var_node(:)
  PetscInt ::iphase
  PetscInt :: ipckrreg !, ithrmtype
  PetscReal :: dif(:)
 
 !   PetscInt :: size_var_node = (option%ndof+1)*size_var_use

  PetscReal, pointer :: t ,p
  PetscReal, pointer :: den(:),h(:),u(:),avgmw(:),pc(:),kvr(:)
  PetscReal, pointer :: xmol(:),satu(:),diff(:)
  PetscInt :: ibase 
  
! PetscReal p1,tmp,co2_phi,co2_poyn,stea,dstea_p,dstea_t,hstea_p,hstea_t,dstea
! PetscReal pckr_swir,xla
  PetscReal p2
  PetscReal pw,dw_kg,dw_mol,hw,sat_pressure,visl,xphi,dco2
  PetscReal dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp
  PetscReal ug
  PetscReal henry,m_nacl
  PetscReal dsteamol,hstea
  PetscReal kr(num_phase)
  PetscReal err,vphi,xm_nacl,x_nacl
  PetscReal temp
  
  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
  !pckr_swir=pckr_sir(1)
  
  ibase=1;               t=>var_node(ibase)
  ibase=ibase+1;         p=>var_node(ibase)
  ibase=ibase+1;         satu=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; den=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; avgmw=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; h=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; u=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; pc=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; kvr=>var_node(ibase:ibase+num_phase-1)
  ibase=ibase+num_phase; xmol=>var_node(ibase:ibase+num_phase*num_spec-1)
  ibase=ibase+num_phase*num_spec;
  diff=>var_node(ibase:ibase+num_phase*num_spec-1)

  select case(iphase)

    case(1) ! only water phase
      p = x(1)
      t = x(2)
      xmol(2)= x(3)
  ! if(xmol(2)<0.D0) xmol(2)=0.D0
  ! if(xmol(2)>1.D0) xmol(2)=1.D0
      if(xmol(2)<0.D0) print *,'tran:',iphase, x(1:3)
      if(xmol(2)>1.D0) print *,'tran:',iphase, x(1:3)
      xmol(1)=1.D0 - xmol(2)
  ! if(xmol(3)<0.D0)xmol(3)=0.D0
      pc(1)=0.D0
      pc(2)=0.D0
      satu(1)=1.D0
      satu(2)= 0.D0
     !p2=p
   ! if( p2>5d4)then
   !  call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
     ! dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,itable)
     ! xphi=fg*1D6/p2
   !else 
      !    call ideal_gaseos_noderiv(p2,t,energyscale,dg,hg,ug)
      !   call visco2(t,dg*fmwco2,visg)
    !endif   
   
   
     ! call Henry_CO2_noderiv(xla,tmp,t,p2, xphi, henry,co2_poyn)
  ! print *,'Henry', xphi,p2, t, henry, co2_poyn  
     
    
     ! xmol(4)=henry*xmol(2)/p
   ! xmol(3)=1.D0-xmol(4)
     ! if(xmol(1)<0.D0) xmol(1)=0.D0
  
    case(2) ! only supercritical CO2 Phase
     
      p = x(1)
      t = x(2)
      xmol(4)= x(3)
      if(xmol(4)<0.D0)print *,'tran:',iphase, x(1:3)
      if(xmol(4)>1.D0) print *,'tran:',iphase, x(1:3)
     
      xmol(3)=1.D0 - xmol(4)
  !   if(xmol(3)<0.D0)xmol(3)=0.D0
   !    pc(1)=0.D0
      pc(2)=0.D0
      satu(1)=eps
      satu(2)=1.D0
  
    case(3) ! water + gas 
      p = x(1)
      t = x(2)
      satu(2)= x(3)
      if(satu(2)<0.D0) print *,'tran:',iphase, x(1:3)
      if(satu(2)>1.D0) print *,'tran:',iphase, x(1:3)
      satu(1)=1.D0  - satu(2)
      if( satu(2)<0.D0)satu(2)=0.D0
      pc(2)=0.D0
      temp = 1D-2
      xmol(1)=1.D0; xmol(2)=0.D0; xmol(3)=temp; xmol(4)=1.D0-xmol(3)

  end select  
        
    
  call PSAT(t, sat_pressure, ierr)

 !   initial guess
   
  err=1.D0

!  print *, 'in 2 phase solver'
 
   !  p2=p*xmol(4)
  !call ideal_gaseos_noderiv(p2,t,energyscale,dg,hg,ug)

  
  !do while(err>1E-8)
  
  p2=p!*xmol(4)
 !   p2=p*xmol(4)

  if(p2>=5d4)then
      !call ideal_gaseos_noderiv(pa,t,energyscale,dg,hg,ug)
     
    call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
    dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,itable)
    dg= dg / fmwco2
    fg= fg * 1.D6 
    hg= hg * fmwco2
  else      
    call ideal_gaseos_noderiv(p2,t,energyscale,dg,hg,ug)
    call visco2(t,dg*fmwco2,visg)
    fg=p2
  endif

  xphi = fg/p2
!   xphi  = 1.d0 ! test-pcl
    
!    call Henry_CO2_noderiv(xla,tmp,t,p*xmol(4),xphi,henry,co2_poyn)
  !  call Henry_duan_sun(p2*1D-5, t,  henry) ! henry = mol/kg/bars
   
  call Henry_duan_sun(t, p2 *1D-5, henry,xphi,m_nacl, m_nacl,sat_pressure*1D-5)
      
  henry= 1D0 / (fmwh2o *1D-3) / (henry*1D-5 )/xphi 
     
     !print *, "translator :; henry ", henry, xphi, p*.99/henry        
     !note: henry = H/phi
     
  !   tmp = (p-sat_pressure)/(henry-sat_pressure)
 !    err= dabs(tmp-xmol(4))
!   tmp=xmol(4)
 ! enddo
  
      
 !print *, 'out 2 phase solver'

  select case(iphase)     
    case(1)
      xmol(4)=xmol(2)*henry/p   
   
 ! if(iphase /= 1)then
      xmol(3)=1.D0-xmol(4)
      if(xmol(3)<0.D0)xmol(3)=0.D0
  !if(xmol(3)<0.D0) xmol(3)=0.D0
    case(2)   
  
      xmol(2)= p*xmol(4)/henry
    !xmol(2)= 1.D0/(1.D0 + 1.D0 /(henry * xmol(4) * p* 1D-5 *xphi * grid%fmwh2o ))
        
      xmol(1)=1.D0-xmol(2)
  !if(xmol(1)<0.D0) xmol(1)=0.D0
  !endif
    case(3)
      temp= sat_pressure /p
      xmol(2)=(1.D0-temp)/(Henry/p - temp)
      xmol(1)= 1.D0- xmol(2)
      xmol(3)=xmol(1) * temp
      xmol(4)= 1.D0-xmol(3)            
  end select


  if(p2<5d4) call visco2(t,dg*fmwco2,visg)
!***************  Liquid phase properties **************************
    avgmw(1)= xmol(1)* fmwh2o + xmol(2) * fmwco2 
    avgmw(2)= xmol(3)* fmwh2o + xmol(4) * fmwco2 
 
   ! pure water
    pw = p   
    if(num_phase>=2)then
      call pflow_pckr_noderiv(num_phase,ipckrreg,satu,pc,kr)
      pw=p !-pc(1)
     ! print *, num_phase,ipckrreg,satu,pc,kr
    endif

    
  !  call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,energyscale,ierr)
   !    call VISW(t,pw,sat_pressure,visl,tmp,tmp2,ierr)
   ! call VISW_FLO(t,dw_mol,visl)
   ! call VISW_noderiv(t,pw,sat_pressure,visl,ierr)
   !  print *,'visw  ',visl,tmp
   ! dif= 1.D-7 !effective diffusion coeff in liquid phase: check pcl

    xm_nacl = m_nacl * fmwnacl
    xm_nacl = xm_nacl /(1.D3 + xm_nacl)
    call nacl_den(t, p*1D-6, xm_nacl, dw_kg) 
    dw_kg = dw_kg * 1D3
    call nacl_vis(t,p*1D-6,xm_nacl,visl)

    diff(1:num_spec)=dif(1)
    diff(num_spec +1: 2*num_spec)=dif(2) ! m^2/s @ 25C

    !apply mixing rules
    ! should be careful with the consistance of mixing rules

!FEHM mixing
!  den(1) = xmol(2)*dg + xmol(1)*dw_mol
! ideal mixing    
  !den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 

! Garcia mixing
    x_nacl = m_nacl/(m_nacl + 1D3/fmwh2o)
   
   ! xmol(1) = xh2o + xnacl
    avgmw(1)= (xmol(1) - x_nacl) * fmwh2o + x_nacl * fmwnacl + xmol(2) * fmwco2
    vphi=1D-6*(37.51D0 + t*(-9.585D-2 + t*(8.74D-4 - t*5.044D-7)))
    den(1)=dw_kg/(1D0-(fmwco2*1D-3-dw_kg*vphi)*xmol(2)/(avgmw(1)*1D-3))
    den(1)=den(1)/avgmw(1)
  
      
 ! Hebach, J. Chem.Eng.Data 2004 (49),p950 
 !   den(1)= 949.7109D0 + p * (0.559684D-6 - 0.00097D-12 * p) &  
 !      + (t+273.15)*(0.883148 - 0.00228*(t+273.15))  
 !  den(1)=dw_kg + (den(1)-dw_kg)*xmol(2)/p*henry
 !  den(1)=den(1)/avgmw(1)
    call wateos_noderiv(t,pw,dw_kg,dw_mol,hw,energyscale,ierr)

    !h(1) = hw * xmol(1) + hg*xmol(2) 
    h(1) = hw
    u(1) = h(1) - pw /dw_mol* energyscale
    diff(1:num_spec) = dif(1)
    kvr(1) = kr(1)/visl
    !xlw = 1.D0
    !avgmw(1) = xlw*fmwh2o+(1.d0-xlw)*fmwco2

!*****************Gas phase properties 2*************************
 
    
  !  if(xmol(3)>eps)then
  !     call steameos(t,p,p2,dstea,dsteamol,dstea_p,dstea_t,&
  !          hstea,hstea_p,hstea_t,energyscale,ierr) 
!print *,'steameos', t,p,pa,p-pa,sat_pressure, dsteamol,hstea,dg,hg,dsteamol/xgw
  !     dsteamol=dsteamol
  !     hstea=hstea!/xgw
!       tmp=dsteamol
   ! else
    dsteamol=dg /xmol(4)*xmol(3)
    hstea=hg 
    !endif


!if((p*xgw)>sat_pressure)then 

 !call steameos(t,p,p-sat_pressure,dstea,dsteamol,dstea_p,dstea_t,&
 !                 hstea,hstea_p,hstea_t,energyscale,ierr) 
! dsteamol=tmp
 !   dstea_p=dg_p
 !   dstea_t=dg_t
!  hstea=tmp
!    hstea_p=hg_p
!    hstea_t=hg_t
!  Now just make the steam same as ideal gas   
 !endif
!**********************************************************************


  ! den(2)= dg*xmol(4) + dw_mol*xmol(3)
    den(2)= 1.D0/(xmol(4)/dg + xmol(3)/dw_mol)

!   call visgas_noderiv(t,pa,p,den(2),visg)
!call visgas_noderiv(t,pa,p,den(2),visg)
! call visco2(t,dg ,visg)
!    visg=8.d-6
    kvr(2)=kr(2)/visg

  
!    den(2)=1.D0/( xga/dg + xgw/dsteamol)  !*c 
  !  den(2)=1.D0/( xga/dg + 1.D0/dsteamol)
!       h(2)=  hg *xmol(4) + hw*xmol(3) 
        h(2)=  hg  
 !   h(2)= ( hg*xga  + hstea*xgw ) 
 !   h(2)= ( hg *xga + hstea ) 
!    u(2)=  h(2)-p/den(2) * energyscale
    u(2)=  h(2)-p/dg * energyscale
    pc(2)=0
   !  print *,'gas phase nonder h::',t,p,h(2),u(2),hg,hstea
    diff(num_spec+1:num_spec*2)= kr(2)* dif(2) * 1.01325D5/p*((t+273.15)/273.15)**1.8D0 ! dirty part
    dco2=den(2)/(p/rgasj/(t+273.15D0)*1D-3)
 !if(t>=100.D0) print *, p,t,xga,xgw,dg,dsteamol,hg, hstea, h(2)
  
!  print *,m_nacl,x_nacl,xm_nacl,visl,den,kr
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! remember: henry coeff is arranged in the vector as 
!    phase          1(water)          2(gas)            3(co2)
!   species      1    2    3       1    2    3        1   2    3
!   values                        1.0  1.0  1.0
!_________________________________________________________________

  
   nullify(t, p, satu, den, avgmw, h,u, pc,kvr,xmol,diff)
   
end subroutine pri_var_trans_mph_ninc_2_2


subroutine pri_var_trans_mph_ninc(x,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg,dif, var_node,itable,m_nacl,ierr,phi_co2, den_co2)
! xgw: water molar fraction in gas phase
! P/Pa, t/(Degree Centigreed), Pc/Pa, Hen(xla=Hen*xga, dimensionless)
 
  implicit none
  
  PetscInt :: num_phase,num_spec
! PetscInt :: num_pricomp
  PetscInt :: size_var_use
  PetscReal x(1:num_spec+1),energyscale
  PetscReal var_node(1:2 + 7*num_phase + 2* num_phase*num_spec)
  PetscReal :: dif(:), m_nacl
  PetscInt ::iphase, itable
  PetscErrorCode :: ierr
  PetscInt :: ipckrreg !, ithrmtype
       
  PetscReal, optional :: phi_co2, den_co2  
  PetscReal :: xphi_co2=1.D0, denco2=1.D0

  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
  if((num_phase == 2).and.( num_spec ==2)) then
    call pri_var_trans_mph_ninc_2_2( x,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg,dif, var_node,itable,m_nacl,ierr,xphi_co2, denco2)
    if(present(phi_co2)) phi_co2=xphi_co2        
    if(present(den_co2))den_co2=denco2
  else 
    print *, 'Wrong phase-specise combination. Stop.'
    stop
  endif
  ! print *, 'ninc: ',x, var_node
   
end subroutine pri_var_trans_mph_ninc   
  
  
subroutine pri_var_trans_mph_winc(x,delx,iphase,energyscale,num_phase,num_spec,&
                    ipckrreg,dif, var_node,itable,m_nacl,ierr)

  implicit none                    

  PetscInt :: num_phase,num_spec
  PetscInt :: size_var_use,size_var_node
    

  PetscReal x(1:num_spec+1),delx(1:num_spec+1),energyscale
  PetscReal var_node(:),m_nacl
  PetscReal :: dif(:)
  PetscInt ::iphase,itable
  PetscErrorCode :: ierr
  PetscInt :: ipckrreg !, ithrmtype
  PetscInt :: ispec
   
   
  PetscReal xx(1:num_spec+1)


  size_var_use = 2 + 7*num_phase + 2* num_phase*num_spec
  size_var_node = (num_spec+2)*size_var_use
  
  do ispec=1,num_spec+1
         xx=x;  xx(ispec)=x(ispec)+ delx(ispec)
  ! note: var_node here starts from 1 to option%ndof*size_var_use
    call pri_var_trans_mph_ninc(xx,iphase,energyscale,num_phase,num_spec,&
                                ipckrreg, dif,&
                                var_node((ispec-1)*size_var_use+1:ispec*size_var_use),itable,m_nacl,ierr)
  enddo


end subroutine pri_var_trans_mph_winc
 
end module   translator_mph_module

