module e4d_setup

  use vars

  implicit none
  integer :: ios
  real :: xorig,yorig,zorig

contains
  !________________________________________________________________
  subroutine setup_e4d
    implicit none
    integer :: sz
    
 
    
    if(my_rank>0) then
       call slave_setup
       return
    end if
 
!#if 0
    call read_input   
    call send_info
    call setup_forward
  
!#endif
!    call send_info_geh

    call send_command(0)
  
  end subroutine setup_e4d
  !________________________________________________________________

  !________________________________________________________________
  subroutine read_input
    implicit none
    
    integer :: nchar,junk,i,check,j
    real, dimension(4) :: etmp 
    logical :: exst

    !get the name of the mesh file
    open(10,file='e4d.inp',status='old',action='read') 
    read(10,*) mshfile
    read(10,*) efile
    read(10,*) mapfile

    close(10)

    !get the number of character is in the prefix
    nchar = 0
    do i=1,40
       if(mshfile(i:i) == '.') then
          nchar = i
          exit
       end if
    end do
    
    !check to see if the files exist

    !!Allocate/read the electrode positions and survey conf
     open(10,file=efile,status='old',action='read')   
       read(10,*,IOSTAT=ios) ne
       !if(ios .ne. 0) goto 1007
       !write(51,*) "NUMBER OF ELECTRODE LOCATIONS = ",ne
       
       
       allocate(e_pos(ne,4))
       do i=1,ne
          read(10,*,IOSTAT=ios) junk,etmp
          !if(ios .ne. 0) goto 1008
          !if(junk>ne) then
          !   write(51,*) 'ELECTRODE ',junk,' IS GREATER THAN THE NUMBER OF ELLECTRODES ',ne
          !   close(51)
          !   call crash_exit
          !end if
          e_pos(junk,1:4)=etmp
       end do

       !!translate electrodes
       open(21,file='trans.txt',status='old')
       read(21,*,IOSTAT=ios) xorig,yorig,zorig
       !if(ios .ne. 0) goto 1009
       close(21)
       e_pos(:,1) = e_pos(:,1)-xorig
       e_pos(:,2) = e_pos(:,2)-yorig
       e_pos(:,3) = e_pos(:,3)-zorig

       !!Read in the survey
       read(10,*,IOSTAT=ios) nm
       allocate(s_conf(nm,4))
       do i=1,nm
          read(10,*,IOSTAT=ios) junk,s_conf(i,1:4)
       end do
       close(10) 

       open(10,file=trim(mapfile),status='old',action='read')
       read(10,*) nmap
       allocate(map_inds(nmap,2),map(nmap)) 
       do i=1,nmap
          read(10,*) map_inds(i,:),map(i)
       end do
       close(10) 
  
  end subroutine read_input
  !________________________________________________________________

  !________________________________________________________________
  subroutine send_info
    implicit none

    call send_command(7)
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    call MPI_BCAST(s_conf, 4*nm,MPI_INTEGER,0,E4D_COMM,ierr)

  end subroutine send_info
  !________________________________________________________________

  !____________________________________________________________________
  subroutine setup_forward
    implicit none
    
    integer :: npre,i,nchr,dim,bflag,ns,itmp
    real :: jnk,jnk1 
    integer :: status(MPI_STATUS_SIZE)
    
    !command slaves to do the setup
    call send_command(2)
   
    !read and send the nodes
    !!Determine the prefix of the tetgen output files
    nchr=len_trim(mshfile)
    do i=1,nchr
       if(mshfile(i:i)=='.') then
          npre=i+1;
          exit
       end if
    end do

    !!OPEN AND READ THE NODE POSITIONS 
    !!nodes(x,y,z) and nbounds
    open(10,file=mshfile(1:npre)//".node",status="old",action="read")
    read(10,*) nnodes,dim,jnk,bflag
  
    allocate(nodes(nnodes,3),nbounds(nnodes))
   
    do i=1,nnodes
       read(10,*) jnk,nodes(i,1:3),jnk1,nbounds(i)
       !write(*,*) i,nodes(i,:),nbounds(i)
       if(nbounds(i)>2) nbounds(i)=7
    end do
    close(10)
    write(*,*) 'master',my_rank,nnodes
    !!Send the nodes to the slaves
    call MPI_BCAST(nnodes, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(nbounds, nnodes, MPI_INTEGER , 0, E4D_COMM, ierr )
   
    !!Read in the elements and send each slave it's assigned elements
    open(10,file=mshfile(1:npre)//".ele",status="old",action="read")
    read(10,*) nelem,dim,jnk
    allocate(elements(nelem,4),zones(nelem))
    do i=1,nelem
       read(10,*) jnk,elements(i,1:4),zones(i)
    end do
    close(10)
    call MPI_BCAST(nelem, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,E4D_COMM,ierr)
       
    !send the assignments
    call send_dists
    
    !get electrode nodes
    call get_electrode_nodes 
 
  end subroutine setup_forward
  !____________________________________________________________________

 !____________________________________________________________________
  subroutine send_dists
    !read inputs and distribute the run info to the slave
    implicit none
    integer :: neven,nextra,ce,i,j,k

    
    !!divide the electrodes up among the slave processes
    !!note ne is the number of electrodes defined in read_inp
    tne = ne  
    neven = ne/(n_rank-1)
    nextra = ne - neven*(n_rank-1)  
    if(allocated(eind)) deallocate(eind)
    allocate(eind(n_rank-1,2))
    
    !!Build eind
    if(neven > 0) then
       ce = 1
       do i=1,n_rank-1-nextra
          eind(i,1) = ce
          eind(i,2) = ce+neven-1
          ce = ce+neven
       end do
       
       if(nextra > 0) then
          do i=n_rank-nextra,n_rank-1    
             eind(i,1) = ce
             eind(i,2) = ce+neven
             ce=ce+neven+1
          end do
       end if
    else
       !!in this case there are more processors than electrodes
       !!THIS SHOULD BE AVOIDED
       do i=1,n_rank-1
          eind(i,1) = i
          eind(i,2) = i
       end do
    end if

    !!send assignments
    call MPI_BCAST(tne, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(e_pos,4*tne ,MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(eind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )

  end subroutine send_dists
  !______________________________________________________________

  !______________________________________________________________
  subroutine get_electrode_nodes
    implicit none
    real, dimension(nnodes) :: dist
    integer :: i
    integer, dimension(1) :: indb
   
    if(.not.allocated(e_nods)) then
       allocate (e_nods(tne))
    end if
  
    do i=1,tne
       dist=sqrt( (e_pos(i,1) - nodes(:,1))**2 + &
            (e_pos(i,2) - nodes(:,2))**2 + &
            (e_pos(i,3) - nodes(:,3))**2)
       
       indb = minloc(dist)
       
       if(my_rank==1) then
          if(dist(indb(1)) > 0.1) then
             !!There should always be a node at all
             !!electrode locations. If not, something is
             !!wrong
             write(*,1001) i,dist(indb)
          end if
       end if
       e_nods(i) = indb(1)

1001   format(" WARNING: ELECTRODE ",I5," HAS BEEN MOVED ",F10.4," TO THE NEAREST NODE")	   
    end	do
    
  end subroutine get_electrode_nodes
  !____________________________________________________________________

  !______________________________________________________________
  subroutine send_command(com)
    !!Send a general command to the slaves
    integer :: com
    
    call MPI_BCAST(com,1,MPI_INTEGER,0,E4D_COMM,ierr)
    
  end subroutine send_command
  !________________________________________________________________
  
  !__________________SLAVE ROUTINES________________________________
  !________________________________________________________________
  subroutine slave_setup
    
    implicit none
    integer :: command
    
100 continue
    !Recieve a command from master
    call MPI_BCAST(command,1,MPI_INTEGER,0,E4D_COMM,ierr)
    
    
    !return to main
    if(command == 0) then
       return

    else if(command == 2) then
       call setup_frun
       goto 100

    else if(command == 7) then
       call receive_info
       goto 100
       
    else if(command == 86) then
       call receive_info_geh
       goto 100
       
    else
       goto 100

    end if
  end subroutine slave_setup
  !________________________________________________________________

  !__________________________________________________________________
  subroutine receive_info
    implicit none
    integer :: i
    
    if(allocated(s_conf)) deallocate(s_conf)
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    allocate(s_conf(nm,4))
    call MPI_BCAST(s_conf, 4*nm,MPI_INTEGER,0,E4D_COMM,ierr)
    
  end subroutine receive_info
  !__________________________________________________________________
  
  !__________________________________________________________________
  subroutine setup_frun
    implicit none
    
    integer :: status(MPI_STATUS_SIZE)
    integer :: neven, nextra, ce, i,itmp
    integer :: beg,end
    logical :: eqf = .false.
    real, dimension(:), allocatable :: dist
    integer, dimension(1) :: tmp
    real :: mx,my,mz
    
  
    !receive nodes from master
    call MPI_BCAST(nnodes, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
  

    allocate(nodes(nnodes,3),nbounds(nnodes))
    call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(nbounds, nnodes, MPI_INTEGER , 0, E4D_COMM, ierr )

    !receive elements from master
    call MPI_BCAST(nelem, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    allocate(elements(nelem,4))
    call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,E4D_COMM,ierr)

    !receive the assignments
    call receive_dists
    
    !build A_map and delA.......................................
    call build_delA

    !!Initialize the Petsc A matrix
    call MatCreateSeqAIJ(PETSC_COMM_SELF,nnodes,nnodes,d_nz,d_nnz,A,perr)
    call MatSetFromOptions(A,perr)
 
    !Get the electrode ownership indexes 
    call get_electrode_nodes
    call MatGetType(A,tp,perr)
    
    !Set up the source and solution vectors
    call VecCreate(PETSC_COMM_SELF,X,perr)
    call VecSetSizes(X,nnodes,PETSC_DECIDE,perr)
    call VecSetFromOptions(X,perr)
 
    !allocate and setup the pole solution vector
    call VecCreate(PETSC_COMM_SELF,psol,perr)
    call VecSetSizes(psol,nnodes,PETSC_DECIDE,perr)
    call VecSetFromOptions(psol,perr)
    allocate(poles(nnodes,my_ne))
    poles = 0

    call VecCreate(PETSC_COMM_SELF,B,perr)
    call VecSetSizes(B,nnodes,PETSC_DECIDE,perr)
    call VecSetFromOptions(B,perr)

    !A_map  and S_map given the coupling info so we can deallocate the 
    !nodes and elements
    deallocate(nodes)
    deallocate(elements)
    allocate(sigma(nelem))
    
  end subroutine setup_frun
  !__________________________________________________________________
  
  !__________________________________________________________________
  subroutine receive_dists
    implicit none
    
    integer :: i
   
    if(allocated(eind)) deallocate(eind)
    if(allocated(e_pos)) deallocate(e_pos)
    allocate(eind(n_rank-1,2))
    
    call MPI_BCAST(tne, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    allocate(e_pos(tne,4))
    call MPI_BCAST(e_pos,4*tne,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(eind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
    
    my_ne = eind(my_rank,2)-eind(my_rank,1)+1
  end subroutine receive_dists
  !__________________________________________________________________

  !____________________________________________________________________
  subroutine build_delA
    use e4d_mat_inv
    implicit none
    
   
    integer :: nnds,i,j,k,l,mcount,rn,cn,ncl,jj
    integer :: lrow
    real, dimension(4,4) :: atild,atildi,temp_a
    integer, dimension(4,4) :: indx
    integer, dimension(:,:),allocatable :: A_tmp
    real :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,evol
    real :: A_kij,lA_kij
    real, dimension(:), allocatable :: t1vec,t2vec
    logical :: ilow,iup
    
    !allocate the petsc matrix preallocation vectors
    allocate(d_nnz(nnodes))
    d_nnz=0
    
    allocate(rows(10*nelem),cols(10*nelem))
    allocate(A_map(10*nelem),S_map(10*nelem))
    allocate(delA(10*nelem))
    
    !!COUNT THE NUMBER OF NODE PAIRS (I.E. COUPLING MATRIX ELEMENTS) AND
    !!REALLOCATE THE COUPLING MATRIX STORAGE VECTORS
    !!These pairs will also accomodate source dependent boundaries
    rows(1)=elements(1,1)
    cols(1)=elements(1,1)
    nvals=1
    mcount = 0
    A_map = 0
    nnds=maxval(elements(:,:))
    
    !!Allocate a temporary matrix to aid in building the 
    !!mapping vector
    allocate(A_tmp(nnds,80))
    A_tmp=0
    A_tmp(rows(1),1) = 1
    A_tmp(rows(1),2) = nvals
    A_tmp(rows(1),3) = cols(1) 
    !!loop over the elements
     
    do k=1,nelem    
       
       !!loop over each node in this element
       do i=1,4
         
          !!rn is a row index of the coupling matrix
          rn=elements(k,i)    
          
          !!loop over each node in this element
          do j=1,4
             
             !!cn is a column index of the coupling matrix
             cn=elements(k,j)
        
          

             !!check to see if this pair (rn,cn) is in the lower
             !!triangle. If not, go to the next pair
             if (cn <= rn) then
                mcount=mcount+1
                !!loop over each pair found thus far (nvals pairs)
                !!and determine if the pair (rn,cn) is already 
                !!represeneted
                if(A_tmp(rn,1) == 0) then
                   nvals = nvals + 1
                   A_tmp(rn,1) = 1
                   A_tmp(rn,2) = nvals
                   A_tmp(rn,3) = cn
                   rows(nvals) = rn
                   cols(nvals) = cn
                   A_map(mcount) = nvals
                           
                else
                   
                   do l=1,A_tmp(rn,1)
                      if(cn == A_tmp(rn,2*l+1)) then
                         A_map(mcount) = A_tmp(rn,2*l)
                         goto 11
                      end if
                   end do
                   
                   !!if were here then no column index was found
                   !!for this row, so we'll add one
                   nvals = nvals+1
                   ncl = A_tmp(rn,1) + 1
                   A_tmp(rn,1) = ncl
                   A_tmp(rn,2*ncl) = nvals
                   A_tmp(rn,2*ncl+1) = cn
                   rows(nvals) = rn
                   cols(nvals) = cn
                   A_map(mcount) = nvals
                
11                 continue
     
                end if

             end if

          end do
          
       end do
       
    end do
    
    !build the petsc allocation vector
    do i=1,nnodes
       !upper part
       d_nnz(i)=A_tmp(i,1)-1;
       do j=1,A_tmp(i,1)
          cn=A_tmp(i,2*j+1)
          d_nnz(cn)=d_nnz(cn)+1
       end do
    end do   
    deallocate(A_tmp)
   
    !!now reallocate for the correct number of matrix values (i.e. pairs)
    allocate(t1vec(nvals),t2vec(nvals))
    t1vec(1:nvals)=rows(1:nvals)
    t2vec(1:nvals)=cols(1:nvals)
    deallocate(rows,cols)
    allocate(rows(nvals),cols(nvals))
    rows=t1vec
    cols=t2vec
    deallocate(t1vec,t2vec)
    allocate(trows(nvals),tcols(nvals))
    
    !!LOOP OVER THE ELEMENTS AND BUILD delA
    mcount = 0
    delA=0;
  
    do k=1,nelem
   
       !!get the 4 nodes for this element 
       x1 = nodes(elements(k,1),1)
       y1 = nodes(elements(k,1),2)
       z1 = nodes(elements(k,1),3)
       
       x2 = nodes(elements(k,2),1)
       y2 = nodes(elements(k,2),2)
       z2 = nodes(elements(k,2),3)
       
       x3 = nodes(elements(k,3),1)
       y3 = nodes(elements(k,3),2)
       z3 = nodes(elements(k,3),3) 
       
       x4 = nodes(elements(k,4),1)
       y4 = nodes(elements(k,4),2)
       z4 = nodes(elements(k,4),3)
       
       !!get the volume of this element
       call get_vol (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,evol)  
   
       if(evol==0) goto 12
       !!build the linear shape function matrix for this element
       atild(1,:) = (/1,1,1,1/)
       atild(2,:) = (/x1,x2,x3,x4/)
       atild(3,:) = (/y1,y2,y3,y4/)
       atild(4,:) = (/z1,z2,z3,z4/)
       temp_a=atild
            
       !!invert temp_a to get atildi
       !!note MIGS routine is located in mat_inv.f      
       call MIGS(temp_a,4,atildi,indx)  
      
12     continue
       !!use atildi to build the coupling coefficients
       !!for this element      
       
       !!loop over each node in this element
       do i=1,4
          !!rn is a row index of the coupling matrix
          rn=elements(k,i)
          
          !!loop over each node in this element
          do j=1,4
             
             !!cn is a column index of the coupling matrix
             cn=elements(k,j)
             
             !!check to see if this pair (rn,cn) is in the lower
             !!triangle. If not, go to the next pair
             if (cn <= rn) then
                mcount=mcount+1
                !!Compute the coupling coefficient for this element
                !!and node pair
                A_kij=0
               
                do l=2,4
                   A_kij=A_kij+(atildi(i,l)*atildi(j,l))
                end do
          
                !!delA(mcount)=lA_kij*evol;
                delA(mcount)=dble(A_kij*evol)
                S_map(mcount)=k
             end if
             
21           continue
          end do
          
       end do
       
    end do  
    
    
  end subroutine build_delA
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine get_vol( x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, volume )

    implicit none
    
    real a(4,4)
    real r4_det
    real volume
    real x1
    real x2
    real x3
    real x4
    real y1
    real y2
    real y3
    real y4
    real z1
    real z2
    real z3
    real z4
    
    a(1:4,1) = (/ x1, x2, x3, x4 /)
    a(1:4,2) = (/ y1, y2, y3, y4 /)
    a(1:4,3) = (/ z1, z2, z3, z4 /)
    a(1:4,4) = (/ 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00 /)
    r4_det = det4(a)
    volume = abs ( r4_det ) / 6.0E+00
    
    return
  end subroutine get_vol
  !____________________________________________________________________
    
  !____________________________________________________________________
  function det4 ( a1 )
    implicit none
    
    real :: a1(4,4)
    real*8 :: a(4,4)
    real :: det4
    integer ::i
    real*8 :: c1,c2,c3,c4
    
    a = dble(a1)
    c1=a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*&
         &(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) 
    
    c2=a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*&
         &(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))
    
    c3=a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*&
         &(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
    
    c4=a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(2,2)*&
         &(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
    
    det4 = real(c1 - c2 + c3 - c4)
    
    return
  end function det4
  !__________________________________________________________________________

  !________________________________________________________________
  subroutine send_info_geh
    implicit none

    call send_command(86)
    nelem = 4*4*4
    nm = nelem
    print *, 'master:', nm
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)

  end subroutine send_info_geh
  !________________________________________________________________

  !__________________________________________________________________
  subroutine receive_info_geh
    implicit none
    integer :: i
    
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    allocate(sigma(nm))
    nelem = nm
    sigma = 0.d0
    print *, 'slave:', nm
    
  end subroutine receive_info_geh
  !__________________________________________________________________
  
 
end module e4d_setup

