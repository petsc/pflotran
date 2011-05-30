#include "finclude/petscsys.h"

PetscInt, parameter :: MAXHEADERLENGTH = 2048
PetscInt, parameter :: MAXSTRINGLENGTH = 512
PetscInt, parameter :: MAXWORDLENGTH = 32
PetscInt, parameter :: IUNIT1 = 15
PetscInt, parameter :: IUNIT2 = 16
PetscInt, parameter :: IUNIT3 = 17
PetscInt, parameter :: IUNIT4 = 18
PetscInt, parameter :: IUNIT_TEMP = 86
PetscInt, parameter :: HHISTORY_LENGTH = 1000
! HHISTORY_LENGTH is the length of the array used to store the differencing
! values h.

! formula weights
PetscReal, parameter :: FMWNACL = 58.44277d0
PetscReal, parameter :: FMWH2O = 18.01534d0  ! kg/kmol h2o
PetscReal, parameter :: FMWCO2 = 44.0098d0
PetscReal, parameter :: FMWAIR = 28.96d0

! conversion factors
PetscReal, parameter :: LOG_TO_LN = 2.30258509299d0
PetscReal, parameter :: LN_TO_LOG = 0.434294481904d0  

! constants
PetscReal, parameter :: IDEAL_GAS_CONST = 8.314472d0   

PetscInt, parameter :: ZERO_INTEGER = 0
PetscInt, parameter :: ONE_INTEGER = 1
PetscInt, parameter :: TWO_INTEGER = 2
PetscInt, parameter :: THREE_INTEGER = 3
PetscInt, parameter :: FOUR_INTEGER = 4
PetscInt, parameter :: FIVE_INTEGER = 5
PetscInt, parameter :: SIX_INTEGER = 6
PetscInt, parameter :: SEVEN_INTEGER = 7
PetscInt, parameter :: EIGHT_INTEGER = 8
PetscInt, parameter :: NINE_INTEGER = 9
PetscInt, parameter :: TEN_INTEGER = 10
PetscInt, parameter :: ELEVEN_INTEGER = 11
PetscInt, parameter :: TWELVE_INTEGER = 12
PetscInt, parameter :: NEG_ONE_INTEGER = -1

PetscMPIInt, parameter :: ZERO_INTEGER_MPI = ZERO_INTEGER
PetscMPIInt, parameter :: ONE_INTEGER_MPI = ONE_INTEGER
PetscMPIInt, parameter :: TWO_INTEGER_MPI = TWO_INTEGER
PetscMPIInt, parameter :: THREE_INTEGER_MPI = THREE_INTEGER
PetscMPIInt, parameter :: FOUR_INTEGER_MPI = FOUR_INTEGER
PetscMPIInt, parameter :: SEVEN_INTEGER_MPI = SEVEN_INTEGER
PetscMPIInt, parameter :: TWELVE_INTEGER_MPI = TWELVE_INTEGER
PetscMPIInt, parameter :: MAXSTRINGLENGTH_MPI = MAXSTRINGLENGTH

PetscInt, parameter :: X_DIRECTION = 1
PetscInt, parameter :: Y_DIRECTION = 2
PetscInt, parameter :: Z_DIRECTION = 3
PetscInt, parameter :: LOWER = 1
PetscInt, parameter :: UPPER = 2

PetscInt, parameter :: TIME_NULL = 0
PetscInt, parameter :: TIME_T = 1
PetscInt, parameter :: TIME_TpDT = 2

PetscInt, parameter :: SORPTION_LINEAR = 1
PetscInt, parameter :: SORPTION_LANGMUIR = 2
PetscInt, parameter :: SORPTION_FREUNDLICH  = 3

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
PetscInt, parameter :: FLASH2_MODE = 6
PetscInt, parameter :: G_MODE = 7

! grid types
PetscInt, parameter :: STRUCTURED_GRID = 1
PetscInt, parameter :: UNSTRUCTURED_GRID = 2
PetscInt, parameter :: AMR_GRID = 3
PetscInt, parameter :: CARTESIAN_GRID = 4
PetscInt, parameter :: CYLINDRICAL_GRID = 5
PetscInt, parameter :: SPHERICAL_GRID = 6
PetscInt, parameter :: STRUCTURED_GRID_MIMETIC = 7
PetscInt, parameter :: UNSTRUCTURED_GRID_MIMETIC = 8

! condition types
PetscInt, parameter :: DIRICHLET_BC = 1
PetscInt, parameter :: PRODUCTION_WELL = -1
PetscInt, parameter :: NEUMANN_BC = 2
PetscInt, parameter :: DIRICHLET_ZERO_GRADIENT_BC = 3
PetscInt, parameter :: ZERO_GRADIENT_BC = 4
PetscInt, parameter :: HYDROSTATIC_BC = 5
PetscInt, parameter :: SEEPAGE_BC = 6
PetscInt, parameter :: MASS_RATE_SS = 7
PetscInt, parameter :: VOLUMETRIC_RATE_SS = 8
PetscInt, parameter :: SCALED_MASS_RATE_SS = 9
PetscInt, parameter :: SCALED_VOLUMETRIC_RATE_SS = 10
PetscInt, parameter :: CONCENTRATION_SS = 11
PetscInt, parameter :: EQUILIBRIUM_SS = 12
PetscInt, parameter :: CONDUCTANCE_BC = 13
PetscInt, parameter :: UNIT_GRADIENT_BC = 14

! concentration subcondition types
PetscInt, parameter :: CONSTRAINT_NULL = 0
PetscInt, parameter :: CONSTRAINT_FREE = 1
PetscInt, parameter :: CONSTRAINT_TOTAL = 2
PetscInt, parameter :: CONSTRAINT_LOG = 3
PetscInt, parameter :: CONSTRAINT_PH = 4
PetscInt, parameter :: CONSTRAINT_MINERAL = 5
PetscInt, parameter :: CONSTRAINT_GAS = 6
PetscInt, parameter :: CONSTRAINT_CHARGE_BAL = 7
PetscInt, parameter :: CONSTRAINT_TOTAL_SORB = 8
PetscInt, parameter :: CONSTRAINT_SUPERCRIT_CO2 = 9

! mineral types
PetscInt, parameter :: MINERAL_REFERENCE = 1
PetscInt, parameter :: MINERAL_KINETIC = 2
PetscInt, parameter :: MINERAL_EQUILIBRIUM = 3

! surface complexation surface types
PetscInt, parameter :: NULL_SURFACE = 0
PetscInt, parameter :: COLLOID_SURFACE = 1
PetscInt, parameter :: MINERAL_SURFACE = 2

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
PetscInt, parameter :: RICHARDS_CONDUCTANCE_DOF = 2

PetscInt, parameter :: GENERAL_LIQUID_PRESSURE_DOF = 1
PetscInt, parameter :: GENERAL_GAS_PRESSURE_DOF = 1
PetscInt, parameter :: GENERAL_AIR_PRESSURE_DOF = 2
PetscInt, parameter :: GENERAL_GAS_SATURATION_DOF = 3
PetscInt, parameter :: GENERAL_LIQUID_FLUX_DOF = 1
PetscInt, parameter :: GENERAL_GAS_FLUX_DOF = 1
PetscInt, parameter :: GENERAL_TEMPERATURE_DOF = 3
PetscInt, parameter :: GENERAL_CONCENTRATION_DOF = 2
PetscInt, parameter :: GENERAL_ENTHALPY_DOF = 4
PetscInt, parameter :: GENERAL_LIQUID_CONDUCTANCE_DOF = -1
PetscInt, parameter :: GENERAL_GAS_CONDUCTANCE_DOF = -2

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
PetscInt, parameter :: POROSITY = 14
PetscInt, parameter :: PHASE = 15
PetscInt, parameter :: MATERIAL_ID = 16
PetscInt, parameter :: GAS_DENSITY_MOL = 17


PetscInt, parameter :: PRIMARY_MOLALITY = 18
PetscInt, parameter :: SECONDARY_MOLALITY = 19
PetscInt, parameter :: TOTAL_MOLALITY = 20
PetscInt, parameter :: PRIMARY_MOLARITY = 21
PetscInt, parameter :: SECONDARY_MOLARITY = 22
PetscInt, parameter :: TOTAL_MOLARITY = 23
PetscInt, parameter :: MINERAL_VOLUME_FRACTION = 24
PetscInt, parameter :: MINERAL_RATE = 25
PetscInt, parameter :: MINERAL_SURFACE_AREA = 26
PetscInt, parameter :: PH = 27
PetscInt, parameter :: SURFACE_CMPLX = 28
PetscInt, parameter :: SURFACE_CMPLX_FREE = 29
PetscInt, parameter :: KIN_SURFACE_CMPLX = 30
PetscInt, parameter :: KIN_SURFACE_CMPLX_FREE = 31
PetscInt, parameter :: PRIMARY_ACTIVITY_COEF = 32
PetscInt, parameter :: SECONDARY_ACTIVITY_COEF = 33
PetscInt, parameter :: SC_FUGA_COEFF = 34
PetscInt, parameter :: PRIMARY_KD = 35
PetscInt, parameter :: TOTAL_SORBED = 36
PetscInt, parameter :: TOTAL_SORBED_MOBILE = 37
PetscInt, parameter :: COLLOID_MOBILE = 38
PetscInt, parameter :: COLLOID_IMMOBILE = 39
PetscInt, parameter :: AGE = 40

! activity coefficients
PetscInt, parameter :: ACT_COEF_FREQUENCY_OFF = 0
PetscInt, parameter :: ACT_COEF_FREQUENCY_TIMESTEP = 1
PetscInt, parameter :: ACT_COEF_FREQUENCY_NEWTON_ITER = 2
PetscInt, parameter :: ACT_COEF_ALGORITHM_LAG = 3
PetscInt, parameter :: ACT_COEF_ALGORITHM_NEWTON = 4
PetscInt, parameter :: NO_BDOT = 5

! structured grid faces
PetscInt, parameter :: NULL_FACE = 0
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

! Tecplot stuff
PetscInt, parameter :: TECPLOT_POINT_FORMAT = 1
PetscInt, parameter :: TECPLOT_BLOCK_FORMAT = 2

PetscInt, parameter :: OBSERVATION_SCALAR = 1
PetscInt, parameter :: OBSERVATION_FLUX = 2
PetscInt, parameter :: OBSERVATION_AT_CELL_CENTER = 1
PetscInt, parameter :: OBSERVATION_AT_COORDINATE = 2

! phase ids
PetscInt, parameter :: LIQUID_PHASE = 1
PetscInt, parameter :: GAS_PHASE = 2

! thermodynamic state of fluid ids
PetscInt, parameter :: LIQUID_STATE = 1
PetscInt, parameter :: GAS_STATE = 2
PetscInt, parameter :: TWO_PHASE_STATE = 3

! variable centerings
PetscInt, parameter :: CELL_CENTERED = 0
PetscInt, parameter :: SIDE_CENTERED = 1

! approaches to coupling reactive transport
PetscInt, parameter :: GLOBAL_IMPLICIT = 0
PetscInt, parameter :: OPERATOR_SPLIT = 1

! dataset types
PetscInt, parameter :: DATASET_SCALAR = 1
PetscInt, parameter :: DATASET_VECTOR = 2
PetscInt, parameter :: DATASET_TENSOR = 3
PetscInt, parameter :: DATASET_HETEROGENEOUS = 4

! stencil type
PetscInt, parameter :: STAR_STENCIL = 1
PetscInt, parameter :: BOX_STENCIL = 2

! grid cell type
PetscInt, parameter :: HEX_TYPE          = 1
PetscInt, parameter :: WEDGE_TYPE        = 2

! grid cell properties
PetscInt, parameter :: TRI_FACE_TYPE     = 1
PetscInt, parameter :: QUAD_FACE_TYPE    = 2
PetscInt, parameter :: MAX_VERT_PER_CELL = 8
PetscInt, parameter :: MAX_DUALS         = 6
PetscInt, parameter :: MAX_VERT_PER_FACE = 4
PetscInt, parameter :: MAX_CELLS_SHARING_A_VERTEX = 16

#define HASH
