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
PetscInt, parameter :: SIX_INTEGER = 6
PetscInt, parameter :: SEVERN_INTEGER = 7
PetscInt, parameter :: NEG_ONE_INTEGER = -1

PetscInt, parameter :: X_DIRECTION = 1
PetscInt, parameter :: Y_DIRECTION = 2
PetscInt, parameter :: Z_DIRECTION = 3
PetscInt, parameter :: LOWER = 1
PetscInt, parameter :: UPPER = 2

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

PetscInt, parameter :: GLOBAL = 1
PetscInt, parameter :: LOCAL = 2
PetscInt, parameter :: NATURAL = 3

! modes
PetscInt, parameter :: NULL_MODE = 0
PetscInt, parameter :: THC_MODE = 1
PetscInt, parameter :: MPH_MODE = 2
PetscInt, parameter :: RICHARDS_MODE = 3
PetscInt, parameter :: REACTIVE_TRANSPORT_MODE = 4
PetscInt, parameter :: IMS_MODE = 5

! grid types
PetscInt, parameter :: STRUCTURED_GRID = 1
PetscInt, parameter :: UNSTRUCTURED_GRID = 2
PetscInt, parameter :: AMR_GRID = 3
PetscInt, parameter :: CARTESIAN_GRID = 4
PetscInt, parameter :: CYLINDRICAL_GRID = 5
PetscInt, parameter :: SPHERICAL_GRID = 6

! condition types
PetscInt, parameter :: DIRICHLET_BC = 1
PetscInt, parameter :: PRODUCTION_WELL = -1
PetscInt, parameter :: NEUMANN_BC = 2
PetscInt, parameter :: ZERO_GRADIENT_BC = 3
PetscInt, parameter :: HYDROSTATIC_BC = 4
PetscInt, parameter :: SEEPAGE_BC = 5
PetscInt, parameter :: MASS_RATE_SS = 6
PetscInt, parameter :: VOLUMETRIC_RATE_SS = 7
PetscInt, parameter :: CONCENTRATION_SS = 9
PetscInt, parameter :: EQUILIBRIUM_SS = 10

! concentration subcondition types
PetscInt, parameter :: CONSTRAINT_NULL = 0
PetscInt, parameter :: CONSTRAINT_FREE = 1
PetscInt, parameter :: CONSTRAINT_TOTAL = 2
PetscInt, parameter :: CONSTRAINT_LOG = 3
PetscInt, parameter :: CONSTRAINT_MINERAL = 4
PetscInt, parameter :: CONSTRAINT_GAS = 5

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
PetscInt, parameter :: THC_PRESSURE_DOF = 1
PetscInt, parameter :: THC_MASS_RATE_DOF = 2
PetscInt, parameter :: THC_TEMPERATURE_DOF = 3
PetscInt, parameter :: THC_CONCENTRATION_DOF = 4
PetscInt, parameter :: THC_ENTHALPY_DOF = 5

PetscInt, parameter :: MPH_PRESSURE_DOF = 1
PetscInt, parameter :: MPH_MASS_RATE_DOF = 2
PetscInt, parameter :: MPH_TEMPERATURE_DOF = 3
PetscInt, parameter :: MPH_CONCENTRATION_DOF = 4
PetscInt, parameter :: MPH_ENTHALPY_DOF = 5

PetscInt, parameter :: RICHARDS_PRESSURE_DOF = 1
PetscInt, parameter :: RICHARDS_MASS_RATE_DOF = 2

! output definitions
PetscInt, parameter :: X_COORDINATE = 1
PetscInt, parameter :: Y_COORDINATE = 2
PetscInt, parameter :: Z_COORDINATE = 3
PetscInt, parameter :: TEMPERATURE = 4
PetscInt, parameter :: PRESSURE = 5
PetscInt, parameter :: LIQUID_SATURATION = 6
PetscInt, parameter :: GAS_SATURATION = 7
PetscInt, parameter :: LIQUID_DENSITY = 8
PetscInt, parameter :: GAS_DENSITY = 9
PetscInt, parameter :: LIQUID_ENERGY = 10
PetscInt, parameter :: GAS_ENERGY = 11
PetscInt, parameter :: LIQUID_MOLE_FRACTION = 12
PetscInt, parameter :: GAS_MOLE_FRACTION = 13
PetscInt, parameter :: PHASE = 14
PetscInt, parameter :: MATERIAL_ID = 15

PetscInt, parameter :: PRIMARY_SPEC_CONCENTRATION = 16
PetscInt, parameter :: SECONDARY_SPEC_CONCENTRATION = 17
PetscInt, parameter :: TOTAL_CONCENTRATION = 18
PetscInt, parameter :: MINERAL_VOLUME_FRACTION = 19
PetscInt, parameter :: MINERAL_SURFACE_AREA = 20

! structured grid faces
PetscInt, parameter :: WEST_FACE = 1
PetscInt, parameter :: EAST_FACE = 2
PetscInt, parameter :: SOUTH_FACE = 3
PetscInt, parameter :: NORTH_FACE = 4
PetscInt, parameter :: BOTTOM_FACE = 5
PetscInt, parameter :: TOP_FACE = 6

! mphase equation of state
PetscInt, parameter :: EOS_SPAN_WAGNER = 1
PetscInt, parameter :: EOS_MRK = 2

! HDF5 stuff
PetscInt, parameter :: HDF5_READ_BUFFER_SIZE = 1000000
!#define HDF5_BROADCAST

#define HASH
