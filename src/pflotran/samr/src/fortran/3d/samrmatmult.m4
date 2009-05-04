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
      subroutine samrccellmatmult3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & stencilsize,
     & ndof,
     & stencil,
     & sgcw,
     & src,
     & dgcw,
     & dst)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer stencilsize, ndof
      integer sgcw, dgcw
      REAL dst(0:ndof-1,CELL3d(ifirst,ilast,dgcw))
      REAL src(0:ndof-1,CELL3d(ifirst,ilast,sgcw))
      REAL stencil(0:ndof*ndof*stencilsize-1, 
     &             CELL3d(ifirst,ilast,0))

      integer i,j,k
c
c***********************************************************************
c
      if(ndof.eq.1) then
         if(stencilsize.eq.7) then
            call samrcellsd7s3d(
     &               ifirst0,ifirst1,ifirst2,
     &               ilast0,ilast1,ilast2,
     &               stencil,
     &               sgcw, src,
     &               dgcw, dst)
         endif
         if(stencilsize.eq.27) then
         endif
      else
         if(stencilsize.eq.7) then
            call samrcellmd7s3d(
     &               ifirst0,ifirst1,ifirst2,
     &               ilast0,ilast1,ilast2,
     &               ndof,
     &               stencil,
     &               sgcw, src,
     &               dgcw, dst)
         endif
         if(stencilsize.eq.27) then
         endif
         
      endif
      return
      end

c
c
      subroutine samrccellmatdiagscale3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & stencilsize,
     & ndof,
     & stencil,
     & lgcw,
     & ldata,
     & rgcw,
     & rdata)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer stencilsize, ndof
      integer lgcw, rgcw
      REAL ldata(0:ndof-1,CELL3d(ifirst,ilast,lgcw))
      REAL rdata(0:ndof-1,CELL3d(ifirst,ilast,rgcw))
      REAL stencil(0:ndof*ndof*stencilsize-1, 
     &             CELL3d(ifirst,ilast,0))

      integer i,j,k
c
c***********************************************************************
c
      return
      end
c
c
      subroutine samrccellmatdiagscalelocal3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & stencilsize,
     & offsets,
     & ndof,
     & stencil,
     & dgcw,
     & ddata)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer stencilsize, ndof
      integer dgcw
      integer offsets(1:stencilsize*ndof*ndof*NDIM)
      REAL ddata(0:ndof-1,CELL3d(ifirst,ilast,dgcw))
      REAL stencil(0:ndof*ndof*stencilsize-1, 
     &             CELL3d(ifirst,ilast,0))

      integer i, j, k, io, jo,ko, be, bs, s, d
      
c
c***********************************************************************
c
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               do s = 1, stencilsize

                  io=offsets(NDIM*ndof*ndof*(s-1)+1)
                  jo=offsets(NDIM*ndof*ndof*(s-1)+2)
                  ko=offsets(NDIM*ndof*ndof*(s-1)+3)

                  do d=0,ndof-1
                  bs = ndof*ndof*(s-1)+d*ndof
                  be = bs + ndof-1
                  stencil(bs:be,i,j,k)=stencil(bs:be,i,j,k)
     &                                *ddata(d,i+io,j+jo,k+ko)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

      subroutine samrcellsd7s3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & stencil,
     & sgcw,
     & src,
     & dgcw,
     & dst)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer sgcw, dgcw
      REAL dst(CELL3d(ifirst,ilast,dgcw))
      REAL src(CELL3d(ifirst,ilast,sgcw))
      REAL stencil(0:6,CELL3d(ifirst,ilast,0))
      REAL sum
      integer i,j,k
c
c***********************************************************************
c
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               sum=0.0
               sum=sum+stencil(BB,i,j,k)*src(i,j,k-1)
               sum=sum+stencil(SS,i,j,k)*src(i,j-1,k)
               sum=sum+stencil(WW,i,j,k)*src(i-1,j,k)
               sum=sum+stencil(PP,i,j,k)*src(i,j,k)
               sum=sum+stencil(EE,i,j,k)*src(i+1,j,k)
               sum=sum+stencil(NN,i,j,k)*src(i,j+1,k)
               sum=sum+stencil(TT,i,j,k)*src(i,j,k+1)
               dst(i,j,k)=sum
c               dst(i,j,k)=stencil(PP,i,j,k)*src(i,j,k)
c     &                   +stencil(WW,i,j,k)*src(i-1,j,k)
c     &                   +stencil(EE,i,j,k)*src(i+1,j,k)
c     &                   +stencil(SS,i,j,k)*src(i,j-1,k)
c     &                   +stencil(NN,i,j,k)*src(i,j+1,k)
c     &                   +stencil(BB,i,j,k)*src(i,j,k-1)
c     &                   +stencil(TT,i,j,k)*src(i,j,k+1)
            enddo
         enddo
      enddo

      return
      end

      subroutine samrcellmd7s3d(
     & ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     & ndof,
     & stencil,
     & sgcw,
     & src,
     & dgcw,
     & dst)

c***********************************************************************
      implicit none
c
      integer ifirst0,ifirst1,ifirst2
      integer ilast0,ilast1,ilast2
      integer ndof
      integer sgcw, dgcw
      REAL dst(0:ndof-1,CELL3d(ifirst,ilast,dgcw))
      REAL src(0:ndof-1,CELL3d(ifirst,ilast,sgcw))
      REAL stencil(0:7*ndof*ndof-1,
     &             CELL3d(ifirst,ilast,0))

      integer i,j,k,l,m
      integer nd2
c
c***********************************************************************
c
      nd2 = ndof*ndof

      do k = ifirst2, ilast2
        do j = ifirst1, ilast1
          do i = ifirst0, ilast0
            do l=0,ndof-1
              dst(l,i,j,k)=0.0
              do m=0,ndof-1
               dst(l,i,j,k)=dst(l,i,j,k)
     &                    +stencil(PP*nd2+m*ndof+l,i,j,k)*src(m,i,j,k)
     &                    +stencil(WW*nd2+m*ndof+l,i,j,k)*src(m,i-1,j,k)
     &                    +stencil(EE*nd2+m*ndof+l,i,j,k)*src(m,i+1,j,k)
     &                    +stencil(SS*nd2+m*ndof+l,i,j,k)*src(m,i,j-1,k)
     &                    +stencil(NN*nd2+m*ndof+l,i,j,k)*src(m,i,j+1,k)
     &                    +stencil(BB*nd2+m*ndof+l,i,j,k)*src(m,i,j,k-1)
     &                    +stencil(TT*nd2+m*ndof+l,i,j,k)*src(m,i,j,k+1)
              enddo
            enddo
          enddo
        enddo
      enddo

      return
      end
