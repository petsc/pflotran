      module utilities_module

      contains

      function rnd()

!     common/rndnr/iseed
      integer*8 iseed

      iseed = iseed*125
      iseed = iseed - (iseed/2796203) * 2796203
      rnd   = iseed/2796203.0
      return
      end function rnd

      function ran1(idum)
      
      implicit none
      
      save

!-----returns a random number in the range (0,1). Set idum to neg.
!     value to initialize

      real*8 :: ran1
      real*8 :: r(97),rm1,rm2
      integer*4 :: idum,iff,ix1,ix2,ix3,j,m1,ia1,ic1,m2,ia2,ic2,m3,ia3,ic3

      parameter (M1  = 259200)
      parameter (IA1 = 7141)
      parameter (IC1 = 54773)
      parameter (RM1 = 1.0/M1)
      parameter (M2  = 134456)
      parameter (IA2 = 8121)
      parameter (IC2 = 28411)
      parameter (RM2 = 1.0/M2)
      parameter (M3  = 243000)
      parameter (IA3 = 4561)
      parameter (IC3 = 51349)

      data iff/0/

      if (idum.lt.0 .or. iff.eq.0) then
        iff=1
        ix1=mod(IC1-idum,M1)
        ix1=mod(IA1*ix1+IC1,M1)
        ix2=mod(ix1,M2)
        ix1=mod(IA1*ix1+IC1,M1)
        ix3=mod(ix1,M3)
        do j=1,97
          ix1=mod(IA1*ix1+IC1,M1)
          ix2=mod(IA2*ix2+IC2,M2)
          r(j)=(float(ix1)+float(ix2)*RM2)*RM1
        enddo
        idum=1
      endif
      ix1=mod(IA1*ix1+IC1,M1)
      ix2=mod(IA2*ix2+IC2,M2)
      ix3=mod(IA3*ix3+IC3,M3)
      j=1+(97*ix3)/M3

!     if (j.gt.97 .or. j.lt.1) pause

      ran1=r(j)
      r(j)=(float(ix1)+float(ix2)*RM2)*RM1
      
      return
      end function ran1

      function ran2(idum)
      
!     implicit none

!-----Minimal random number generator of Park and Miller 
!     in the range (0,1)

      parameter (IA = 16807)
      parameter (IM = 2147483647)
      parameter (AM = 1.0/IM)
      parameter (IQ = 127773)
      parameter (IR = 2836)
      parameter (NTAB = 32)
      parameter (NDIV = 1+(IM-1)/NTAB)
      parameter (EPS  = 1.2e-7)
      parameter (RNMX = 1.0-EPS)

      dimension iv(NTAB)

      iy = 0
      if (idum.le.0 .or. iy.eq.0) then
        if (-idum .lt. 1) then
          idum = 1
        else 
          idum = -idum
        endif
        do j = NTAB+7,0,-1
          k = idum/IQ
          idum = IA*(idum-k*IQ)-IR*k
          if (idum .lt. 0) idum = idum+IM
          if (j .lt. NTAB) iv(j) = idum
        enddo
        iy = iv(1)
      endif
      k = idum/IQ
      idum = IA*(idum-k*IQ)-IR*k
      if (idum .lt. 0) idum = idum+IM
      j= iy/NDIV
      iy = iv(j)
      iv(j) = idum
      temp = AM*iy
      if (temp .gt. RNMX) then
        ran2 = RNMX
      else 
        ran2 = temp
      endif

      return
      end function ran2
!----------------------------------------------------------------------------      
 subroutine Natural2LocalIndex(ir, nl, llist, llength)
  implicit none
  integer nl, ir,na, l_search, itt, llength
  integer llist(*)
  
  integer  nori0, nori1, nori
  
  
  nl=-1
  l_search = llength
  
  na = ir!-1 
  itt=0
  nori0 =1
  nori1 = llength
  if(na>=llist(1) .and. na <= llist(llength))then
  do while(l_search > 1 .and.itt<=50)
  
    itt=itt+1
    if(na == llist(nori0))then
      nl = nori0
      exit
    elseif(na == llist(nori1))then
       nl = nori1
      exit
    endif   
     
     ! nori = int((real(nori0 + nori1))/ 2.) + mod ( nori0 + nori1,2 )
    nori =  int(floor(real(nori0+nori1)/2D0 + .75D0))
    if( na > llist(nori)) then
      nori0 = nori
    elseif(na < llist(nori))then
      nori1 = nori
    else
      if(na == llist(nori))then
        nl = nori
        exit
      else
        print *, 'wrong index', na, nori, llist(nori); stop
      endif  
    endif
    l_search = nori1-nori0
    if (itt>=40)then
      print *, na, nori0,nori1,nori, llist(nori0), llist(nori1)
      if (itt>=50) stop
    endif
  enddo  
 endif         
          
 end subroutine Natural2LocalIndex

end module utilities_module
