
module oil_pckr_module
  public

#include "definitions.h"

  PetscReal, private, parameter:: pckr_sat_water_cut = 1.D0 - 5.D-7
 ! PetscReal, private :: pckr_swir, pckr_soir, pckr_sgir
  PetscReal, private, allocatable, save :: kr_so(:,:),kr_sg(:,:)
  PetscInt, private, save :: kr_so_n, kr_sg_n
  PetscReal, private, allocatable, save :: pc_so(:,:),pc_sg(:,:)
  PetscInt, private, save :: pc_so_n, pc_sg_n
  PetscReal, private, save :: alaph_ow= 1.50D-3, alaph_ao=2.64D-4
  
  contains 

subroutine oil_pckr_init(myrank)

  use Fileio_module
   
  PetscMPIInt myrank

  character(len=MAXSTRINGLENGTH) :: string 
  character(len=MAXWORDLENGTH) :: word !, strtim
  character(len=MAXCARDLENGTH) :: card

  open(60, file="3ph_perm.in", action="read", status="old")
  do
    call fiReadFlotranString(IUNIT1, string, ierr)
    if(ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    call fiCharsToUpper(word,len_trim(word))
    call fiReadCard(word,card,ierr)

    if (myrank == 0) print *, card

    select case(card)
      case('KROw')
      call fiReadInt(string,kr_so_n,ierr) 
      call fiDefaultMsg(myrank,'rows of list for krow',ierr)
      allocate(kr_so(1:kr_so_n,3))
      i=0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg(myrank,'KROW',ierr)

        if (string(1:1) == '.' .or. string(1:1) == '/') exit
        i = i + 1
        if(i>kr_so_n)then
          print *,'Reading exceed the array , exit...'
          exit
        endif

        do j=1,3
          call fiReadDouble(string,kr_so(i,j),ierr)
          call fiDefaultMsg(myrank,'krso',ierr)
        enddo
      enddo   

    case('KROG')
      call fiReadInt(string,kr_sg_n,ierr) 
      call fiDefaultMsg(myrank,'rows of list for krog',ierr)
      allocate(kr_sg(1:kr_sg_n,3))
      i=0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg(myrank,'KROG',ierr)

        if (string(1:1) == '.' .or. string(1:1) == '/') exit
        i = i + 1
        if(i>kr_sg_n)then
        print *,'Reading exceed the array , exit...'
        exit
      endif

      do j=1,3
        call fiReadDouble(string,kr_sg(i,j),ierr)
        call fiDefaultMsg(myrank,'krsg',ierr)
      enddo
    enddo
    case default
    end select
  enddo

  close(60)
  end subroutine oil_pckr_init


  subroutine oil_pckr_noderiv(ipckrtype,pckr_swir, pckr_soir, pckr_lambda, &
                 pckr_alpha,pckr_m,pckr_pcmax,s_w,s_o,pc,kr,pckr_beta,pckr_pwr) 
    implicit none 
    PetscInt :: ipckrtype
    PetscReal :: s_w,s_o
    PetscReal :: pckr_swir,pckr_soir, pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_pwr
    PetscReal :: pc(1:3),kr(1:3)
    PetscReal :: pckr_beta

    PetscReal :: se_w,se_t,swir,sw0
    PetscReal :: pcmax
    PetscReal :: pckr_betac,pckr_un, m_r
!   PetscReal :: lam,ala,se,um,un,upc,upc_s,temp,ser,uum,betac,st

    ! if(present(pckr_beta))
      pckr_betac=pckr_beta
      swir=pckr_swir
      sw0=1.d0
      pcmax=pckr_pcmax


    pc=0.D0; kr=0.D0

    select case(ipckrtype)   
      case(10)
        se_w=(s_w-pckr_swir)/(1.D0-pckr_swir)
        se_t=(s_o-pckr_soir)/(1.D0-pckr_soir)
        if (se_w>0.D0) kr(1)=se_w
        if (se_t>0.D0) kr(3)=se_t
        kr(2)=1.D0 - s_w -s_o
        pc(2)=0.D0
        if(se_w<1D-4)then
          pc(1)=pcmax
        else
          pc(1)=(1.D0-se_w)*pcmax
        endif

       if(se_t<1D-4)then
         pc(3)=pcmax*1D-3
       else
         pc(3)=(1.D0-se_t)*pcmax*1D-3
       endif
 

     !  print *, 'oil pckr:', ipckrtype, kr
      case(11) ! Parker 1987
        m_r = 1.D0/pckr_m
        pckr_un=1.D0/(1.D0- pckr_m)

        se_w = (s_w - pckr_swir)/( 1.D0 - pckr_swir)
        se_t = (s_w + s_o - pckr_swir)/( 1.D0 - pckr_swir) 
    !  if(se_w<1D-5)then
    !    pc(1)=pcmax
    !   else
    !   pc(1)=((se_w)**pckr_m-1D0)**(1.D0/pckr_un)/alaph_ow
    !     if(pc(1)>pcmax) pc(1)=pcmax
    !   endif
    !   if(se_t<1D-5)then
    !    pc(3)=pcmax
    !   else
    !   pc(3)=((se_t)**pckr_m-1D0)**(1.D0/pckr_un)/alaph_ao
    !     if(pc(3)>pcmax) pc(1)=pcmax
    !   endif


      if (se_w<1D-5)then ! no water phase
         kr(1)=0.D0
        if( se_t<1D-5)then !no oil phase
            kr(2)=1.D0; kr(3)=0.D0
        elseif((1.D0-se_t)<1D-5)then !no gas phase
          kr(2)=0.D0; kr(3)=1.D0
        else         ! oil and gas phase
          kr(2)= (1-se_t)**.5D0 *(1.D0 - se_t**m_r)**(2.D0*pckr_m)   
          kr(3)=  (se_t-se_w)**5D-1 *(1.D0  -(1.D0 - se_t**m_r)**pckr_m)**2.D0
          endif
      else    ! have water
         kr(1) = se_w **.5D0 * (1.D0 - (1.D0 - se_w**m_r)**pckr_m)**2.D0
       if( (se_t-se_w)<1D-5)then ! no oil phase
          kr(2)= (1-se_t)**.5D0 *(1.D0 - se_t**m_r)**(2.D0*pckr_m)
                kr(3)=0.D0
       elseif((1.D0-se_t)<1D-5)then !no gas phase
          kr(2)=0.D0
        if((se_t-se_w)<1D-5)then ! no oil phase 
          kr(3)=0.D0
        else
          kr(3)= (se_t-se_w)**5D-1 *((1.D0 - se_w**m_r)**pckr_m -(1.D0 - se_t**m_r)**pckr_m)**2.D0
                endif
       else
           kr(1) = se_w **.5D0 * (1.D0 - (1.D0 - se_w**m_r)**pckr_m)**2.D0
         kr(2) = (1-se_t)**.5D0 *(1.D0 - se_t**m_r)**(2.D0*pckr_m)    
             kr(3) = (se_t-se_w)**5D-1 *((1.D0 - se_w**m_r)**pckr_m -(1.D0 - se_t**m_r)**pckr_m)**2.D0
           endif 
      endif 
       case(12) ! Stone I 
     
     
     case(13) ! Stone II           
    end select
  end subroutine oil_pckr_noderiv   
end module oil_pckr_module

