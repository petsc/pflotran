#include "include/finclude/petsc.h"

PetscInt, parameter :: MAXBCREGIONS = 50000 
PetscInt, parameter :: MAXBCBLOCKS = 50000 
PetscInt, parameter :: MAXSTRINGLENGTH = 512
PetscInt, parameter :: MAXWORDLENGTH = 32
PetscInt, parameter :: MAXCARDLENGTH = 4
PetscInt, parameter :: MAXNAMELENGTH = 20
PetscInt, parameter :: MAXPERMREGIONS = 35000
!PetscInt, parameter :: MAXINITREGIONS = 80000
PetscInt, parameter :: MAXINITREGIONS = 8400000
PetscInt, parameter :: MAXSRC = 10
PetscInt, parameter :: MAXSRCTIMES = 100
PetscInt, parameter :: IUNIT1 = 15
PetscInt, parameter :: IUNIT2 = 16
PetscInt, parameter :: IUNIT3 = 17
PetscInt, parameter :: IUNIT4 = 18
PetscInt, parameter :: HHISTORY_LENGTH = 1000
! HHISTORY_LENGTH is the length of the array used to store the differencing
! values h.

PetscInt, parameter :: ZERO_INTEGER = 0
PetscInt, parameter :: ONE_INTEGER = 1
PetscInt, parameter :: TWO_INTEGER = 2
PetscInt, parameter :: THREE_INTEGER = 3
PetscInt, parameter :: FOUR_INTEGER = 4
PetscInt, parameter :: FIVE_INTEGER = 5
PetscInt, parameter :: NEG_ONE_INTEGER = -1

PetscInt, parameter :: X_DIRECTION = 1
PetscInt, parameter :: Y_DIRECTION = 2
PetscInt, parameter :: Z_DIRECTION = 3

! Classes
PetscInt, parameter :: NULL_CLASS = 0
PetscInt, parameter :: FLOW_CLASS = 1
PetscInt, parameter :: TRANSPORT_CLASS = 2

! Macros that are used as 'dm_index' values.  --RTM
PetscInt, parameter :: ONEDOF = 1
PetscInt, parameter :: NPHASEDOF = 2
PetscInt, parameter :: THREENPDOF = 3
PetscInt, parameter :: NFLOWDOF = 4
PetscInt, parameter :: NTRANDOF = 5
PetscInt, parameter :: NPHANCOMPDOF = 6
PetscInt, parameter :: NPHANSPECDOF = 7
PetscInt, parameter :: NPHANSPECNCOMPDOF = 8
PetscInt, parameter :: VARDOF = 9

PetscInt, parameter :: GLOBAL = 1
PetscInt, parameter :: LOCAL = 2
PetscInt, parameter :: NATURAL = 3

! modes
PetscInt, parameter :: NULL_MODE = 0
PetscInt, parameter :: RICHARDS_MODE = 1
PetscInt, parameter :: MPH_MODE = 2
PetscInt, parameter :: COND_MODE = 3
PetscInt, parameter :: TWOPH_MODE = 4
PetscInt, parameter :: VADOSE_MODE = 5
PetscInt, parameter :: LIQUID_MODE = 6
PetscInt, parameter :: OWG_MODE = 7
PetscInt, parameter :: FLASH_MODE = 8
PetscInt, parameter :: TH_MODE = 9
PetscInt, parameter :: THC_MODE = 10
PetscInt, parameter :: RICHARDS_LITE_MODE = 11

! grid types
PetscInt, parameter :: STRUCTURED = 1
PetscInt, parameter :: UNSTRUCTURED = 2
PetscInt, parameter :: STRUCTURED_CARTESIAN = 10
PetscInt, parameter :: STRUCTURED_CYLINDRICAL = 11
PetscInt, parameter :: STRUCTURED_SPHERICAL = 12

! condition types
PetscInt, parameter :: DIRICHLET_BC = 1
PetscInt, parameter :: NEUMANN_BC = 2
PetscInt, parameter :: MASS_RATE = 3
PetscInt, parameter :: ZERO_GRADIENT_BC = 5
PetscInt, parameter :: HYDROSTATIC_BC = 6
PetscInt, parameter :: SEEPAGE_BC = 7

! coupler types
PetscInt, parameter :: INITIAL_COUPLER_TYPE = 1
PetscInt, parameter :: BOUNDARY_COUPLER_TYPE = 2
PetscInt, parameter :: SRC_SINK_COUPLER_TYPE = 3
PetscInt, parameter :: COUPLER_IPHASE_INDEX = 1

! connection types
PetscInt, parameter :: INTERNAL_CONNECTION_TYPE = 1
PetscInt, parameter :: BOUNDARY_CONNECTION_TYPE = 2
PetscInt, parameter :: INITIAL_CONNECTION_TYPE = 3
PetscInt, parameter :: SRC_SINK_CONNECTION_TYPE = 4

! dofs for each mode
PetscInt, parameter :: RICHARDS_PRESSURE_DOF = 1
PetscInt, parameter :: RICHARDS_CONCENTRATION_DOF = 3
PetscInt, parameter :: RICHARDS_TEMPERATURE_DOF = 2
PetscInt, parameter :: RICHARDS_ENTHALPY_DOF = 3

PetscInt, parameter :: MPH_PRESSURE_DOF = 1
PetscInt, parameter :: MPH_CONCENTRATION_DOF = 3
PetscInt, parameter :: MPH_TEMPERATURE_DOF = 2
PetscInt, parameter :: MPH_ENTHALPY_DOF = 3

PetscInt, parameter :: THC_PRESSURE_DOF = 1
PetscInt, parameter :: THC_CONCENTRATION_DOF = 3
PetscInt, parameter :: THC_TEMPERATURE_DOF = 2
PetscInt, parameter :: THC_ENTHALPY_DOF = 3

! HDF5 stuff
PetscInt, parameter :: HDF5_READ_BUFFER_SIZE = 1000000
!#define HDF5_BROADCAST

#define HASH
