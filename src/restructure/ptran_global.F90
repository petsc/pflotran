!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_global.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_global.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:31:15  lichtner
! Added tconv, and tunit; isrc1 and revised maximum limits for nrgmx.
!
! Revision 1.2  2004/01/10 18:32:06  lichtner
! Began work on 2 phase capability.
!
! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
! initial entry
!

!  pFLOTRAN Version 1.0 LANL
!-----------------------------------------------------------------------
!  Date             Author(s)                Comments/Modifications
!-----------------------------------------------------------------------
!  May  2003        Peter C. Lichtner        Initial Implementation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  module ptran_global_module

  implicit none
  save
  public

#define PETSC_AVOID_DECLARATIONS
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscda.h"
#undef PETSC_AVOID_MPIF_H
#undef PETSC_AVOID_DECLARATIONS

#define PETSC_TRUE 1
#define PETSC_FALSE 0

! by default, all variables are private
! private

    ! If using_pflowGrid == PETSC_TRUE, then some parts of ptran_init 
    ! will not be executed, since they are made redundant by 
    ! pflowGrid_new() and pflowGrid_setup().
    PetscTruth :: using_pflowGrid = PETSC_FALSE

!======================================================================
!   input parameters: default values
!======================================================================
    integer :: igeom, nx, ny, nz
    integer :: kmax = 2, newton_max = 32, icut_max = 8, iwarn = 0
    integer :: loglin = 1
    integer :: idebug = 0, ibg1 = 0, ibg2 = 0
    integer :: ireg = 0, iregfld = 0
    integer :: iprint = 0, iact = 0, irestart = 0
    integer :: modetr = 2, molal = 1, isothrm = 0, iphas = 1
    real*8  :: w = 0.5d0, wlam = 1.d0, tol = 1.d-12, tolexp
    real*8, pointer :: dx0(:), dy0(:), dz0(:), rd(:)

!======================================================================

    integer :: myrank,commsize,npx,npy,npz
    integer :: ibc,iblkfmt=1,isurf
    integer :: ibug, icut, icutcum, itpetscum, nstep, &
               newton, newtcum, &
               itpetsc(100)

    integer :: i1,i2,j1,j2,k1,k2,initreg

    real*8  :: t,cb,r2norm,wm1,wlam1,rho0=1.d0
    
!...........................................................
!   time step
!...........................................................
    integer :: kplot
    real*8  :: tplot(100), dt = 1.d0, dtmax

!...........................................................
!   initial parameters
!...........................................................    
    real*8 :: temp0 = 25.d0, tempini = 25.d0, tk, pref0 = 1.d5
    integer :: icomprs = 0

!...........................................................
!   coupled flow and transport
!...........................................................
    integer :: mode = 2, iphase
    real*8, parameter :: slcutoff=1.d-4
    integer :: ipor = 0
    real*8 :: tolpor = 0.d0

!...........................................................
!   transport
!...........................................................
    real*8 :: vlx0,vly0,vlz0,vgx0,vgy0,vgz0,por0,sat0,tor0
    real*8 :: difaq,delhaq,difgas,dgexp

!...........................................................
!   dual continuum
!...........................................................
    integer :: idcdm = 0, ndloc1 = 1, ndloc2 = 0

!...........................................................
!   constants
!...........................................................
    real*8  :: zero=0.d0, quarter=0.25d0, half=0.5d0, &
               one=1.d0, two=2.d0, three=3.d0, four=4.d0, &
               five=5.d0, six=6.d0, seven=7.d0, &
               eight=8.d0, fnine=9.d0, ten=10.d0, &
               onesixth=0.1666666666667d0, onethrd=0.3333333333333d0, &
               twthrds=0.6666666666667d0

!...........................................................
! physical constants
!...........................................................
     real*8 :: pi      = 3.141592653589793d0
     real*8 :: aln10   = 2.30258509299d0
     real*8 :: aloge   = 0.434294482d0, faraday = 96485.d0
     real*8 :: aln2    = 0.6931471806d0
     real*8 :: navogadro = 6.022136736d+23
     real*8 :: rgas    = 1.98726d0 ![cal/K/mol]
     real*8 :: rgasj   = 8.3143    ![J/K/mol]
     real*8 :: rgaskj  = 8.3143d-3 ![kJ/K/mol]
     real*8 :: rgasjj  = 8.3143d-2 ![10 * rgaskj]
     real*8 :: tkelvin = 273.15d0
     real*8 :: tk0     = 298.15d0
     real*8 :: wh2o    = 0.0180153d0
     real*8 :: fmwh2o  = 18.01534d0, fmwco2 = 44.0098d0
!...........................................................
!   unit conversion
!...........................................................    
    real*8 :: yrsec = 31536000.d0, uyrsec = 3.1709792d-8

!...........................................................
!   connections, grid geometry
!...........................................................
    integer :: nconn,nconnx,nconny
    integer :: nconnbc,nconnbc_w,nconnbc_e,nconnbc_n,nconnbc_s, &
               nconnbc_t,nconnbc_b
    integer :: nlx,nly,nlz,nxs,nys,nzs,nlxy,nlxz,nlyz,nldof,nlmax
    integer :: ngx,ngy,ngz,ngxs,ngys,ngzs,ngxy,ngxz,ngyz,ngdof,ngmax
    integer :: nxe,nye,nze,ngxe,ngye,ngze
    integer :: istart,iend,jstart,jend,kstart,kend
    integer :: nxy,nmax
!...........................................................
! time stepping
!...........................................................
  integer :: iaccel = 1, imax
  real*8 :: tfac(13)
  real*8 :: dcmax0, dcdt, aaa
  data tfac/2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, &
            1.8d0, 1.6d0, 1.4d0, 1.2d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0/

  real*8 :: tconv      ! time conversion factor: s -> x
  character*2 :: tunit ! input time units
!...........................................................
! solver
!...........................................................
  real*8 ::  eps = 1.d-10, rtol_petsc, atol_petsc, dtol_petsc
  integer :: maxits_petsc

!...........................................................
! dimension parameters
!...........................................................
  integer, parameter :: &
    ncmx    = 16,   & ! primary species
    ncxmx   = 60,   & ! secondary species
    ngmx    = 6,    & ! gaseous species
    nmmx    = 50,   & ! total number of minerals
    nkmx    = 50,   & ! kinetic mineral reactions
    npmx    = 2,    & ! max no. parallel reactions???
    nxkmx   = 0,    & ! kinetic aqueous homogeneous reactions
    nexmx   = 6,    & ! exchange cations
    nsitxmx = 4,    & ! exchange sites
    nsitmx  = 6,    & ! surface complex sites
    nscxmx  = 4,    & ! surface complexes
    nrgmx   = 1200, & ! max no. of regions (COMP+BCON+SOUR, MNIR)
    nrgbcmx = 1000, & ! max no. of boundary regions
    nsrcmx  = 10,   & ! max no. of source/sinks
    nstbmx  = 100,  & ! max no. entries in source/sink
    nbrkmx  = 100,  & ! max no. breakthough points
    ndimmx  = nmmx+ncxmx+ngmx+nxkmx ! database matrix

!...........................................................
! i/o parameters
!...........................................................
  integer            :: iunit1 = 2, iunit2 = 3, ptran_iunit4 = 10, iunit5 = 11
  integer, parameter :: dbaselen = 256, strlen = 256, wordlen = 32, &
                        namlen = 20, cardlen = 4
  integer            :: nruns = 0, iread_vel = 0, iread_sat = 0, &
                        iflgco2 = 0, iflux = 0
  real*8             :: profile(ncmx+1,100)
  character(len=namlen) :: flowvx,flowvy,flowvz,fsat

!...........................................................
! flags
!...........................................................
  integer            :: iflgcut = 0

!...........................................................
! thermodynamic database
!...........................................................
  integer, parameter    :: ntmpmx = 8, ncxmx0 = 300, nmmx0 = 300, &
                           ngmx0 = 10
  integer               :: icase
  character(len=strlen) :: dbpath
  real*8                :: temptab(10),alogkeh(10), &
                           coefeh(5),coef(ndimmx,5)
  real*8 :: alnkeh

!-----Note: the following values are taken directly from the EQ3/6 database.
!     The values at temperatures 500 and 550 are taken to be equal to
!     the value at 300 degrees C.

      data alogkeh/-91.0454d0,-83.1028d0,-74.0521d0,-65.8632d0, &
          -57.8929d0,-51.6850d0,-46.7266d0,-42.6842d0,0.d0,0.d0/

      data temptab/ &
         0.d0,     25.d0,    60.d0,    100.d0, &
       150.d0,    200.d0,   250.d0,    300.d0,   500.d0,   550.d0/
     
!...........................................................
! debye-huckel parameters
!...........................................................
  real*8 :: atab(10),btab(10),bdotab(10)
  real*8 :: adebye,bdebye,bextend(ncmx),bextendx(ncxmx),sionic0
       
!-----debye huckel a
      data atab / &
        .4939d0,   .5114d0,   .5465d0,   .5995d0, &
        .6855d0,   .7994d0,   .9593d0,  1.2180d0,  1.2180d0,  1.2180d0/

!-----debye huckel b
      data btab / &
          .3253d0,  .3288d0,  .3346d0,  .3421d0, &
          .3525d0,  .3639d0,  .3766d0,  .3925d0,  .3925d0,  .3925d0/

!-----bdot
      data bdotab / &
         .0374d0,  .0410d0,  .0440d0,  .0460d0, &
         .0470d0,  .0470d0,  .0340d0,  .0000d0,  .0000d0,  .0000d0/

!...........................................................
! Pitzer model
!...........................................................

!...........................................................
! primary species
!...........................................................
  character(len=namlen) :: nam(ncmx)
  character(len=namlen) :: maspec = ' '
  integer               :: ncomp = 1, ncpri = 1, nmass, nmat
  integer               :: iloop
  integer               :: jph,joh,jo2,jo2aq,jpe,jco2,jhco3, &
                           jfe3,jfe2,jh2o,jo2g,jco2g,jh2g, &
                           iph,ife3,ife2,ihco3
  integer               :: imaster,jstep
  real*8                :: z(ncmx)
  real*8                :: a0(ncmx)
  real*8                :: wt(ncmx)
  real*8                :: atol(ncmx)
  real*8                :: guess(ncmx,nrgmx)
!...........................................................
! initial conditions
!...........................................................
  character(len=namlen) :: ncon(ncmx,nrgmx)
  integer               :: itype(ncmx,nrgmx)
  real*8                :: ctot(ncmx,nrgmx)

!...........................................................
! boundary conditions
!...........................................................
  integer               :: nblkbc=0,ibndtyp(nrgbcmx)
  integer               :: iface(nrgbcmx),invface(6),ibcreg(4)=0
  integer               :: i1bc(nrgbcmx),i2bc(nrgbcmx), &
                           j1bc(nrgbcmx),j2bc(nrgbcmx), &
                           k1bc(nrgbcmx),k2bc(nrgbcmx)
  integer               :: iregbc1(nrgbcmx),iregbc2(nrgbcmx)
  real*8                :: tempbc(nrgbcmx)
  real*8                :: psibnd(ncmx,nrgbcmx),ccbnd(ncmx,nrgbcmx)
  real*8                :: psigbnd(ncmx,nrgbcmx),pgasbnd(ngmx,nrgbcmx)
  real*8                :: cxbnd(ncxmx,nrgbcmx)
  real*8                :: xexbnd(nexmx,nrgbcmx)
  real*8                :: psisrc(ncmx,nrgbcmx)
  real*8                :: psorpbnd(ncmx,nrgbcmx)
  real*8                :: dwbc(nrgbcmx)
  
!...........................................................
! secondary species
!...........................................................
  character(len=namlen) :: namcx(ncxmx0)
  integer               :: ncmplx
  real*8                :: zx(ncxmx)
  real*8                :: ax0(ncxmx)
  real*8                :: wtx(ncxmx)
  real*8                :: eqhom(ncxmx)
  real*8                :: shom(ncmx,ncxmx0)

!...........................................................
! gas species
!...........................................................
  character(len=namlen) :: namg(ngmx)
  integer               :: ngas = 0
  real*8                :: eqgas(ngmx)
  real*8                :: sgas(ncmx,ngmx0)
  real*8                :: wtgas(ngmx)
  
!...........................................................
! total minerals
!...........................................................
  character(len=namlen) :: namrl(nmmx)
  integer               :: mnrl = 0
  real*8                :: smnrl(ncmx,nmmx0)
  real*8                :: vbar(nmmx)
  real*8                :: alnk(nmmx)
  real*8                :: ze(nmmx),ze0(nmmx)
  real*8                :: wtmin(nmmx)
  
!...........................................................
! kinetic minerals
!...........................................................
  character(len=namlen) :: namk(nkmx)
  character(len=namlen) :: namprik(ncmx,nkmx)
  character(len=namlen) :: namseck(ncxmx,nkmx)
  integer               :: nkin = 0
  integer               :: iregkin(nkmx)
  integer               :: ndxkin(nkmx)
  integer               :: npar(nkmx)
  integer               :: npar1(nkmx)
  integer               :: npar2(nkmx)
  integer               :: nkinpri(nkmx)
  integer               :: nkinsec(nkmx)
  real*8                :: qkmax = 5.d0, pwrsrf = 0.6666666666666667d0
  real*8                :: skinpri(ncmx,nkmx)
  real*8                :: skinsec(ncxmx,nkmx)
  real*8                :: tolpos(nkmx)
  real*8                :: delh(nkmx)
  real*8 :: rkfa00(nkmx),rkfb00(nkmx),rkf00(nkmx),surf00(nkmx)
  
  integer :: itypkini(100), itypkin(nkmx), &
             jpri(ncmx,nkmx),isec(ncmx,nkmx)

  real*8 ::          rkfa0(nkmx),   rkfb0(nkmx),     rkf0(nkmx), &
                     rkfa(nkmx),    rkfb(nkmx),      rkf(nkmx),  &
                     beta(nkmx),    betb(nkmx),      sigma(nkmx),&
                     eqkin(nkmx),   skin(ncmx,nkmx), fkin(nkmx), &
                     vbarkin(nkmx), wtkin(nkmx),      &
                     rlim0(nkmx),   rlim(nkmx)
     
!...........................................................
! kinetic aqueous reactions
!...........................................................
  character(len=namlen) :: namcxk(nxkmx)
  character(len=namlen) :: namrxnaq(nxkmx)
  integer               :: ncxkin = 0, nparcxk1(nxkmx), nparcxk2(nxkmx)
  integer               :: itypkiniaq(100)
  real*8                :: eqcxk(nxkmx)
  real*8                :: axk0(nxkmx)
  real*8                :: zxk(nxkmx)
  real*8                :: wtxk(nxkmx)
  real*8                :: zeaq(nxkmx)
  real*8                :: skpri(ncmx,nxkmx)
  
!...........................................................
! surface complexation
!...........................................................
  character(len=namlen) :: namscx(nsitmx)
  character(len=namlen) :: namsite(nsitmx),namsrf(nsitmx)
  character(len=namlen) :: namcoll(nsitmx)
  integer               :: nsrfmin = 0, nsrfmx = 0, nsrfsit = 0, &
                           ncolsrf = 0
  integer               :: nsite1(nkmx),nsite2(nkmx)
  integer               :: nsorp1(nsitmx),nsorp2(nsitmx)
  integer               :: msorp(nkmx)
  real*8                :: zsite(nsitmx)
  real*8                :: zsrf(nscxmx)
  real*8                :: eqsorp(nscxmx),eqsorp0(nscxmx)
  real*8                :: ssorp(ncmx,nscxmx)
  real*8                :: wtsrf(nscxmx),csorpini(nscxmx), &
                           dcsorp(ncmx,nscxmx)
  real*8                :: siteden0(nrgmx,nsitmx),sited0mf(2,nsitmx)
  real*8                :: coverage(nkmx),areamass(nkmx)
  real*8                :: rkfsrf(nscxmx),rkbsrf(nscxmx)
  
!...........................................................
! ion exchange
!...........................................................
  character(len=namlen) :: namex(nexmx),namcat(nexmx)
  integer               :: nexmax = 0, nexsolid = 0, ncollex = 0, &
                           nexsite = 0, ionex = 0
  integer               :: nsitex1(nkmx),nsitex2(nkmx), &
                           nex1(nsitxmx),nex2(nsitxmx)
  integer               :: mex(nkmx),jex(nexmx)
  real*8                :: cec0(nrgmx,nsitxmx),cec0mf(2,nsitxmx)
  real*8                :: eqiex(nexmx),alogex(nexmx)
  real*8                :: xexini(nexmx)

!...........................................................
! monod kinetics
!...........................................................
  integer :: nmonod=0
  
!...........................................................
! decay
!...........................................................
  integer :: ndecay = 0
  real*8 :: xlamhalf(ncmx)
  
!...........................................................
! colloids
!...........................................................
  integer :: ncoll = 0

!...........................................................
! compression
!...........................................................
      real*8  :: fshom(ncmx,ncxmx),ffshom(ncmx,ncmx,ncxmx), &
                 cshom(ncmx,ncxmx),sshom(ncxmx,ncmx*(ncmx+1)/2)
            
      integer :: ki(ncmx,ncxmx),kki(ncmx,ncmx,ncxmx), &
                 jcmpr(ncmx,ncxmx),ncmpr(ncxmx), &
                 jpsi(ncmx*ncmx),jpsig(ncmx*ncmx), &
                 lc(ncmx),llc(ncmx,ncmx)
                 
      real*8  :: fsgas(ncmx,ncxmx),ffsgas(ncmx,ncmx,ncxmx), &
                 csgas(ncmx,ncxmx),ssgas(ncxmx,ncmx*(ncmx+1)/2)
     
      integer :: lg(ncmx),llg(ncmx,ncmx),kg(ncmx,ngmx), &
                 kkg(ncmx,ncmx,ncxmx),njg(ngmx),jg(ncmx,ngmx)

      real*8  :: fskin(nkmx,ncmx),ffskin(nkmx,ncmx*(ncmx+1)/2)
      
      integer :: kkmin(nkmx),kminj(nkmx,ncmx*(ncmx+1)/2), &
                 kminl(ncmx,ncmx*(ncmx+1)/2),kmind(nkmx),kinj(nkmx,ncmx)
                 
!...........................................................
! region porosity, pressure, temperature
!...........................................................
  integer :: i1reg(nrgmx),i2reg(nrgmx), &
             j1reg(nrgmx),j2reg(nrgmx), &
             k1reg(nrgmx),k2reg(nrgmx)
  real*8  :: por_reg(nrgmx),tor_reg(nrgmx),pref_reg(nrgmx),temp_reg(nrgmx)
  
!...........................................................
! source/sink
!...........................................................
  character(len=namlen) :: nconsrc(ncmx,nrgmx)
  integer :: nblksrc,nsrc,isrc1
  integer :: is1(nrgmx),is2(nrgmx), &
             js1(nrgmx),js2(nrgmx), &
             ks1(nrgmx),ks2(nrgmx)
             
  integer :: itypsrc(ncmx,nrgmx)
             
  real*8  :: timesrc(nstbmx,nsrcmx),tempsrc(nstbmx,nsrcmx), &
             qsrc(nstbmx,nsrcmx),csrc(nstbmx,nsrcmx),ctotsrc(ncmx,nsrcmx)
  
!...........................................................
! solid solutions
!...........................................................
  integer :: isolidss=0
  
!...........................................................
! breakthrough curves
!...........................................................
  integer :: ibrkcrv, &
             i1brk(nbrkmx),i2brk(nbrkmx), &
             j1brk(nbrkmx),j2brk(nbrkmx), &
             k1brk(nbrkmx),k2brk(nbrkmx),ibrktyp(nbrkmx), &
             ibrkface(nbrkmx)
  
#undef PETSC_TRUE
#undef PETSC_FALSE

end module ptran_global_module
