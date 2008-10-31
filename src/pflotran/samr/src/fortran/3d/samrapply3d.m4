define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(PP,0)dnl
define(WW,1)dnl
define(EE,2)dnl
define(SS,3)dnl
define(NN,4)dnl
define(BB,5)dnl
define(TT,6)dnl

include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c
      subroutine samrapply7ptstencil3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & stencil,
     & ufirst0,ufirst1,ufirst2,ulast0,ulast1,ulast2,
     & u,
     & ffirst0,ffirst1,ffirst2,flast0,flast1,flast2,
     & f,
     & rfirst0,rfirst1,rfirst2,rlast0,rlast1,rlast2,
     & r)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer ufirst0,ufirst1,ufirst2
      integer ulast0,ulast1,ulast2
      integer ffirst0,ffirst1,ffirst2
      integer flast0,flast1,flast2
      integer rfirst0,rfirst1,rfirst2
      integer rlast0,rlast1,rlast2

      REAL stencil(0:6,CELL3d(ifirst,ilast,0))
      REAL u(CELL3d(ufirst,ulast,0))
      REAL f(CELL3d(ffirst,flast,0))
      REAL r(CELL3d(rfirst,rlast,0))

      integer i,j,k

c***********************************************************************
c
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               r(i,j,k)=f(i,j,k)- stencil(PP,i,j,k)*u(i,j,k)
     &                           +stencil(WW,i,j,k)*u(i-1,j,k)
     &                           +stencil(EE,i,j,k)*u(i+1,j,k)
     &                           +stencil(SS,i,j,k)*u(i,j-1,k)
     &                           +stencil(NN,i,j,k)*u(i,j+1,k)
     &                           +stencil(BB,i,j,k)*u(i,j,k-1)
     &                           +stencil(TT,i,j,k)*u(i,j,k+1)
               
            enddo
         enddo
      enddo

      return
      end
