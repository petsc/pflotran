module e4d_run

  use vars
  integer :: pf_com
  logical :: first_sol = .true.
  integer :: mcomm
  real*8 :: pf_time
contains
  
  !_____________________________________________________________________
  subroutine run_e4d
  
    implicit none
    
    if(my_rank>0) then
       call slave_run
       return
    end if

    if(.not. allocated(pf_sol)) allocate(pf_sol(pflotran_solution_vec_size))
    if(.not. allocated(sigma)) allocate(sigma(nelem))
    
100 continue

    call get_mcomm
    if(mcomm==1) then

       call get_pf_time      !!get the pflotran solution time
       call get_pf_sol       !!get the pflotran solution
       call map_pf_e4d       !!map/transform the solution to the E4D mesh
       call send_sigma       !!send the transformed solution to the slaves 
!geh
!#if 0
       call send_command(3)  !!instruct slaves to build A matrix
       call send_command(5)  !!instruct slaves to build KSP solver
       call send_command(6)  !!instruct slaves to solve  
       call get_dpred        !!assemble the simulated data
!#endif
       goto 100

    else
       call send_command(0)  !!instruct slaves to exit
       return

    end if
    
  end subroutine run_e4d
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_mcomm
    implicit none
!geh    call MPI_BCAST(mcomm,1,MPI_INTEGER,0,PFE4D_COMM,ierr)
    call MPI_BCAST(mcomm,1,MPI_INTEGER,0,PFE4D_MASTER_COMM,ierr)
  end subroutine get_mcomm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_pf_time
    implicit none
    call MPI_BCAST(pf_time,1,MPI_DOUBLE_PRECISION,0,PFE4D_MASTER_COMM, &
                   ierr)
  end subroutine get_pf_time
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_pf_sol
    implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
    integer ::  status(MPI_STATUS_SIZE)
    PetscReal, pointer :: vec_ptr(:)

    !!this code will change depending on how we send the PF solution
    !!to the E4D master process here
!geh    call MPI_RECV(pf_sol,nelem,MPI_REAL,0,0,PFE4D_COMM,status,ierr)
    call VecScatterBegin(pflotran_scatter,pflotran_solution_vec_mpi, &
                         pflotran_solution_vec_seq, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(pflotran_scatter,pflotran_solution_vec_mpi, &
                       pflotran_solution_vec_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecGetArrayF90(pflotran_solution_vec_seq,vec_ptr,ierr)
    pf_sol = vec_ptr
    call VecRestoreArrayF90(pflotran_solution_vec_seq,vec_ptr,ierr)
  end subroutine get_pf_sol
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine map_pf_e4d
    implicit none
    integer :: i 
    character*40 :: filename, word 
    !this routine will eventually do the mapping. 
    !sigma = pf_sol
    !sigma=pf_time/1e6
    sigma=0.1
    do i=1,nmap
       sigma(map_inds(i,1))=sigma(map_inds(i,1))+pf_sol(map_inds(i,2))*map(i)
    end do
   
    write(*,*) pf_time
    write(word,'(i15.15)') int(pf_time)
    filename = 'sigma_' // &
               trim(adjustl(pflotran_group_prefix)) // &
               trim(adjustl(word)) // &
               '.txt' 
    write(*,*) filename
    open(unit=86,file=trim(filename),status='replace',action='write')
    write(86,*) nelem
    do i = 1, nelem
       write(86,*) sigma(i)
    enddo
    close(86)

  end subroutine map_pf_e4d
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine send_sigma
    implicit none
    
    call send_command(4)
    call MPI_BCAST(sigma, nelem,MPI_REAL,0,E4D_COMM,ierr)
    
  end subroutine send_sigma
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_dpred
    implicit none
    integer :: opt
    integer :: i,j
    integer :: nadd
    integer :: nbuff
    integer, dimension(nm*2) :: ibuff
    real, dimension(nm) :: rbuff
    integer ::  status(MPI_STATUS_SIZE)
    character*40 :: filename, word 


    !instruct slave to assemble and send the prediceted data
    call send_command(8)
    
    !allocate dpred if not already done and zero
    if(.not. allocated(dpred)) then
       allocate(dpred(nm))
    end if
    dpred = 0
    rbuff = 0
    ibuff = 1
    do i=1,n_rank-1 
       call MPI_RECV(nbuff,1,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
       call MPI_RECV(ibuff(1:nbuff),nbuff,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
       call MPI_RECV(rbuff(1:nbuff),nbuff,MPI_REAL,i,0,E4D_COMM,status,ierr)
       
       do j=1,nbuff
          dpred(ibuff(j))=dpred(ibuff(j))+rbuff(j)
       end do
       
    end do

      
    
    write(word,'(i15.15)') int(pf_time)
    filename = 'e4d_' // &
               trim(adjustl(pflotran_group_prefix)) // &
               trim(adjustl(word)) // &
               '.dpd' 
    write(*,*) filename
    open(unit=86,file=trim(filename),status='replace',action='write')
    write(86,*) nm
    do i = 1, nm
       write(86,'(I6,4I6,F15.8)') i,s_conf(i,1:4),dpred(i)
    enddo
    close(86)
    
    !!output the predicted data to file
    !call output_dpred
  end subroutine get_dpred
  !____________________________________________________________________
  

  !______________________________________________________________
  subroutine send_command(com)
    !!Send a general command to the slaves
    integer :: com
    call MPI_BCAST(com,1,MPI_INTEGER,0,E4D_COMM,ierr)
  end subroutine send_command
  !________________________________________________________________
  



  !___________SLAVE SUBROUTINES________________________________________

  !____________________________________________________________________
  subroutine slave_run
    implicit none
    integer :: command
    
100 continue
    !Recieve a command from master
    call MPI_BCAST(command,1,MPI_INTEGER,0,E4D_COMM,ierr)
    
    !return to main
    if(command == 0) then
       return

    else if(command == 3) then     
       call build_A
       goto 100
         
    else if(command == 4) then
       call receive_sigma
       goto 100

    else if(command == 5) then
       call build_ksp
       goto 100

    else if(command == 6) then
       call forward_run
       goto 100
       
    else if(command == 8) then
       call send_dpred 
       goto 100
       
    else 
       goto 100

    end if

  end subroutine slave_run
  !____________________________________________________________________
  !__________________________________________________________________
  subroutine receive_sigma

    implicit none
    character(len=32) :: filename, word
    integer :: i
    integer, save :: num_calls = 0
    call MPI_BCAST(sigma, nelem,MPI_REAL,0,E4D_COMM,ierr)   
!geh
#if 0 
write(filename,*) num_calls
write(word,*) my_rank
filename = 'sigma_' // trim(adjustl(word)) // '_' // &
           trim(adjustl(filename)) // '.txt'
open(unit=86,file=filename)
do i = 1, nelem
  write(86,*) sigma(i)
enddo
close(86)
#endif
num_calls = num_calls + 1
!    print *, 'sigma received by slave', sigma(16)
  end subroutine receive_sigma
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine build_A
    implicit none
    
    integer, dimension(nnodes) :: ncolss
    integer, dimension(50) :: colss
    
    integer :: i,lrnl,lrnu,row,col,rbv,cbv,j,ifn
    logical :: ilow,iup
    
    integer :: rw   
    integer :: ncls 
    integer :: cls(50) 
    
    
    !zero A
    call MatZeroEntries(A,perr)
    
    do i=1,10*nelem
       row=rows(A_map(i))
       col=cols(A_map(i))
       rbv=nbounds(row)
       cbv=nbounds(col)
       
       !lower triangle
       if((rbv>=2 .and. rbv<=6) .or. (cbv>=2 .and. cbv<=6)) then
          !one or both nodes are on a boundary so set to zero for bc's
          val(1) = 0
       else
          
          val(1) = sigma(S_map(i))*delA(i)
       end if
       
       prn(1) = row-1
       pcn(1) = col-1
       
       call MatSetValues(A,1,prn,1,pcn,val,ADD_VALUES,perr)          
       !upper triangle
       if(row .ne. col) then
          call MatSetValues(A,1,pcn,1,prn,val,ADD_VALUES,perr)
       end if
       
    end do
     
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,perr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,perr)
    
    
    
  end subroutine build_A
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine build_ksp
    implicit none
#include "finclude/petscksp.h"
#include "finclude/petscksp.h90"
    real*8 :: rtol = 1e-6
    real*8 :: atol = 1e-35
    real*8 :: dtol = 500
    integer :: maxints = 10000
    !KSPSetTolerances(KSP ksp,double rtol,double atol,double dtol,int maxits);
    
    !Set up the KSP context
    call KSPCreate(PETSC_COMM_SELF,KS,perr)
    call KSPSetOperators(KS,A,A,SAME_PRECONDITIONER,perr)
    call KSPGetPC(KS,P,perr)
    !call KSPSetType(KS,KSPGMRES,perr) !use default
    !call KSPGMRESSetRestart(KS,1000,perr);
    !call KSPGetTolerances(KS,rtol,atol,dtol,maxints,perr)
    call KSPSetTolerances(KS,rtol,atol,dtol,maxints,perr)
    call KSPSetFromOptions(KS,perr)
    
  end subroutine build_ksp
  !__________________________________________________________________

  !____________________________________________________________________
  subroutine forward_run
    implicit none
    integer :: i,m,n,niter,j,enum
    integer, dimension(1) :: eindx
    real, dimension(2) :: pck
    PetscScalar :: val
    real :: tstart, tend
    
    do i=1,my_ne
       !call cpu_time(tstart)
       
       call VecGetArrayF90(psol,vloc,ierr)
       vloc(1:nnodes)=dble(poles(:,i))
       call VecRestoreArrayF90(psol,vloc,ierr) 
       enum=eind(my_rank,1)+i-1
          
       val=0.0
       call VecSet(B,val,perr)
       
       !if(i_flg) then
       !   call Add_Jpp(i)
       !end if
       
       eindx(1)=e_nods(enum)
       val=1.0
       call VecSetValues(B,1,eindx-1,val,ADD_VALUES,perr)
       !if(tank_flag) call VecSetValues(B,1,i_zpot-1,-val,ADD_VALUES,perr)
           
       call VecAssemblyBegin(B,perr)
       call VecAssemblyEnd(B,perr) 
       
       call KSPSolve(KS,B,psol,perr)
       
       call VecGetArrayF90(psol,vloc,ierr)
       poles(:,i)= real(vloc(1:nnodes))
       call VecRestoreArrayF90(psol,vloc,ierr)
       
       call KSPGetIterationNumber(KS,niter,perr)
       call cpu_time(tend)
       pck(1)=tend-tstart
       pck(2)=real(niter)
       !call MPI_SEND(pck,2,MPI_REAL,0,1,MPI_COMM_WORLD,ierr)
       write(*,*) "Slave ",my_rank," solved for pole ",eind(my_rank,1)+i-1,'in ',tend-tstart,' seconds and ',niter,' iters'
       
    end do
    
    if(first_sol) then
       call KSPSetInitialGuessNonzero(KS,PETSC_TRUE,perr)
       first_sol=.false.
    end if
    
    call KSPDestroy(KS,perr)
    
    
  end subroutine forward_run
  !____________________________________________________________________

   !__________________________________________________________________
    subroutine send_dpred
      implicit none
      integer :: flg,i
 
      call assemble_data
      call MPI_SEND(nmy_drows,1,MPI_INTEGER,0,0,E4D_COMM, ierr)
      call MPI_SEND(my_drows,nmy_drows,MPI_INTEGER,0,0,E4D_COMM,ierr)
      call MPI_SEND(my_dvals,nmy_drows,MPI_REAL,0,0,E4D_COMM,ierr)
       
    end subroutine send_dpred
    !__________________________________________________________________
     
    !___________________________________________________________
    subroutine assemble_data
      implicit none
      integer :: i,a,b,m,n,e1,e2,indx,p
      
      
      e1=eind(my_rank,1)
      e2=eind(my_rank,2) 
      if(.not.allocated(my_drows)) then
         nmy_drows=0
         do i=1,nm
            a=s_conf(i,1)
            b=s_conf(i,2)
            if((a>=e1 .and. a<=e2) .or. (b>=e1.and. b<=e2)) then
               nmy_drows=nmy_drows+1
            end if
         end do
         allocate(my_drows(nmy_drows),my_dvals(nmy_drows))
      end if
      
      indx=0
      my_dvals=0
      
     
      do i=1,nm
         a = s_conf(i,1)
         b = s_conf(i,2)
         m = s_conf(i,3)
         n = s_conf(i,4)
         
         if((a>=e1 .and. a<=e2) .or. (b>=e1.and. b<=e2)) then
            indx=indx+1
            my_drows(indx)=i
            
            do p=e1,e2
               if(p==a) then
                  if(m.ne.0) my_dvals(indx)=my_dvals(indx) + real(poles(e_nods(m),a-e1+1))
                  if(n.ne.0) my_dvals(indx)=my_dvals(indx) - real(poles(e_nods(n),a-e1+1))
               end if
               if(p==b) then
                  if(m.ne.0) my_dvals(indx)=my_dvals(indx) - real(poles(e_nods(m),b-e1+1))
                  if(n.ne.0) my_dvals(indx)=my_dvals(indx) + real(poles(e_nods(n),b-e1+1))
               end if
            end do
         end if
      end do
     
  end subroutine assemble_data
  !___________________________________________________________
end module e4d_run
