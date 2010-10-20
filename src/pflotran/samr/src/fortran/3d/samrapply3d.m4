
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
               r(i,j,k)=f(i,j,k)-(stencil(PP,i,j,k)*u(i,j,k)
     &                           +stencil(WW,i,j,k)*u(i-1,j,k)
     &                           +stencil(EE,i,j,k)*u(i+1,j,k)
     &                           +stencil(SS,i,j,k)*u(i,j-1,k)
     &                           +stencil(NN,i,j,k)*u(i,j+1,k)
     &                           +stencil(BB,i,j,k)*u(i,j,k-1)
     &                           +stencil(TT,i,j,k)*u(i,j,k+1))
               
            enddo
         enddo
      enddo

      return
      end

c
c
      subroutine samrapply7ptblkstencil3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & ndof, 
     & stencil,
     & ugcw,
     & u,
     & fgcw,
     & f,
     & rgcw,
     & r)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer ugcw
      integer fgcw
      integer rgcw
      integer ndof
      REAL stencil(0:6,ndof,ndof,CELL3d(ifirst,ilast,0))
      REAL u(ndof,CELL3d(ifirst,ilast,ugcw))
      REAL f(ndof,CELL3d(ifirst,ilast,fgcw))
      REAL r(ndof,CELL3d(ifirst,ilast,rgcw))

      integer i,j,k,d

c***********************************************************************
c
c not coded efficiently

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
           do i = ifirst0, ilast0
	     r(:,i,j,k)=f(:,i,j,k)
	     do d = 1, ndof
               r(:,i,j,k)=r(:,i,j,k)-(stencil(PP,:,d,i,j,k)*u(d,i,j,k)
     &                           +stencil(WW,:,d,i,j,k)*u(d,i-1,j,k)
     &                           +stencil(EE,:,d,i,j,k)*u(d,i+1,j,k)
     &                           +stencil(SS,:,d,i,j,k)*u(d,i,j-1,k)
     &                           +stencil(NN,:,d,i,j,k)*u(d,i,j+1,k)
     &                           +stencil(BB,:,d,i,j,k)*u(d,i,j,k-1)
     &                           +stencil(TT,:,d,i,j,k)*u(d,i,j,k+1))
             enddo  
           enddo
         enddo
      enddo

      return
      end
      
      subroutine pflotranpcflux3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & ufirst0,ufirst1,ufirst2,ulast0,ulast1,ulast2,
     & stencilsize,
     & ndof,
     & stencil,
     & u,
     & flux0, flux1, flux2)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer ufirst0,ufirst1,ufirst2
      integer ulast0,ulast1,ulast2
      integer ndof
      integer stencilsize
      REAL stencil(stencilsize,ndof,ndof,CELL3d(ifirst,ilast,0))
      REAL u(ndof,CELL3d(ufirst,ulast,0))
      REAL flux0(ndof,SIDE3d0(ifirst,ilast,0))
      REAL flux1(ndof,SIDE3d1(ifirst,ilast,0))
      REAL flux2(ndof,SIDE3d2(ifirst,ilast,0))
      
      integer i,j,k,d

c***********************************************************************
c
c for now we assume that the stencil is a 7 pt stencil
c we also assume that the stencil operator is symmetric

      flux0=0.0
      flux1=0.0
      flux2=0.0

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
           do i = ifirst0, ilast0
	     do d=1,ndof
             flux0(:,i,j,k)=flux0(:,i,j,k)
     &       	         +stencil(2,:,d,i,j,k)*(u(d,i,j,k)-u(d,i-1,j,k))
             flux1(:,i,j,k)=flux1(:,i,j,k)
     &                   +stencil(4,:,d,i,j,k)*(u(d,i,j,k)-u(d,i,j-1,k))
             flux2(:,i,j,k)=flux2(:,i,j,k)
     &                   +stencil(6,:,d,i,j,k)*(u(d,i,j,k)-u(d,i,j,k-1))
	     enddo 
           enddo
         enddo
      enddo

c     we are forced to adjust the upper faces differently      
      i=ilast0
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
	   do d=1,ndof
             flux0(:,i+1,j,k)= flux0(:,i+1,j,k)
     &       	        +stencil(3,:,d,i,j,k)*(u(d,i+1,j,k)-u(d,i,j,k))
           enddo
         enddo
      enddo

c     we are forced to adjust the upper face differently      
      j=ilast1
      do k = ifirst2, ilast2
         do i = ifirst0, ilast0
	   do d=1,ndof
             flux1(:,i,j+1,k)= flux1(:,i,j+1,k)
     &	                 +stencil(5,:,d,i,j,k)*(u(d,i,j+1,k)-u(d,i,j,k))
           enddo
         enddo
      enddo

c     we are forced to adjust the upper face differently      
      k=ilast2
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0
	   do d=1,ndof
             flux2(:,i,j,k+1)=flux2(:,i,j,k+1)
     &	                 +stencil(7,:,d,i,j,k)*(u(d,i,j,k+1)-u(d,i,j,k))
           enddo
         enddo
      enddo

      return
      end

      subroutine pflotranpcapply3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & ufirst0,ufirst1,ufirst2,ulast0,ulast1,ulast2,
     & ffirst0,ffirst1,ffirst2,flast0,flast1,flast2,
     & sfirst0,sfirst1,sfirst2,slast0,slast1,slast2,
     & rfirst0,rfirst1,rfirst2,rlast0,rlast1,rlast2,
     & stencilsize,
     & ndof,
     & stencil,
     & u,
     & flux0, flux1, flux2,
     & f,
     & src,
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
      integer sfirst0,sfirst1,sfirst2
      integer slast0,slast1,slast2
      integer rfirst0,rfirst1,rfirst2
      integer rlast0,rlast1,rlast2
      integer ndof
      integer stencilsize
      REAL stencil(0:stencilsize-1,CELL3d(ifirst,ilast,0))
      REAL u(0:ndof-1,CELL3d(ufirst,ulast,0))
      REAL flux0(0:ndof-1,SIDE3d0(ifirst,ilast,0))
      REAL flux1(0:ndof-1,SIDE3d1(ifirst,ilast,0))
      REAL flux2(0:ndof-1,SIDE3d2(ifirst,ilast,0))
      REAL f(0:ndof-1,CELL3d(ffirst,flast,0))
      REAL src(0:ndof-1,CELL3d(sfirst,slast,0))
      REAL r(0:ndof-1,CELL3d(rfirst,rlast,0))
      integer i,j,k,d

c***********************************************************************
c
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               r(:,i,j,k)=f(:,i,j,k)-(flux0(:,i+1,j,k)-flux0(:,i,j,k)
     &                               +flux1(:,i,j+1,k)-flux1(:,i,j,k)
     &                               +flux2(:,i,j,k+1)-flux2(:,i,j,k))
     &                               -src(:,i,j,k)*u(:,i,j,k)
            enddo
         enddo
      enddo

      return
      end

      subroutine samrsetjacobiansrccoeffs3d(
     &     ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &     depth,
     &     stencil,
     &     sfirst0,sfirst1,sfirst2,slast0,slast1,slast2,
     &     src)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer sfirst0,sfirst1,sfirst2
      integer slast0,slast1,slast2
      integer depth

      REAL stencil(0:depth-1,CELL3d(ifirst,ilast,0))
      REAL src(CELL3d(sfirst,slast,0))

      integer i,j,k

c***********************************************************************
c
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               src(i,j,k)=(stencil(PP,i,j,k)
     &                    +stencil(WW,i,j,k)
     &                    +stencil(EE,i,j,k)
     &                    +stencil(SS,i,j,k)
     &                    +stencil(NN,i,j,k)
     &                    +stencil(BB,i,j,k)
     &                    +stencil(TT,i,j,k))
            enddo
         enddo
      enddo

      return
      end

