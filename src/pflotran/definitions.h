#define MAXBCREGIONS 50000
#define MAXBCBLOCKS 50000
#define MAXSTRINGLENGTH 512
#define MAXWORDLENGTH 32
#define MAXCARDLENGTH 4
#define MAXNAMELENGTH 20
#define MAXPERMREGIONS 35000
!#define MAXINITREGIONS 80000
#define MAXINITREGIONS 8400000
#define MAXSRC 10
#define MAXSRCTIMES 100
#define IUNIT1 15
#define IUNIT2 16
#define IUNIT3 17
#define IUNIT4 18
#define HHISTORY_LENGTH 1000
! HHISTORY_LENGTH is the length of the array used to store the differencing
! values h.

! Macros that are used as 'dm_index' values.  --RTM
#define ONEDOF 1
#define NPHASEDOF 2
#define THREENPDOF 3
#define NDOF 4
#define NPHANCOMPDOF 5
#define NPHANSPECDOF 6
#define NPHANSPECNCOMPDOF 7
#define VARDOF 8

#define GLOBAL 1
#define LOCAL 2
#define NATURAL 3

! modes
#define NULL_MODE 0
#define RICHARDS_MODE 1
#define MPH_MODE 2
#define COND_MODE 3
#define TWOPH_MODE 4
#define VADOSE_MODE 5
#define LIQUID_MODE 6
#define OWG_MODE 7
#define FLASH_MODE 8
#define TH_MODE 9
#define THC_MODE 10

! grid types
#define STRUCTURED 1
#define UNSTRUCTURED 2

! boundary condition types
#define DIRICHLET_BC 1
#define NEUMANN_BC 2
