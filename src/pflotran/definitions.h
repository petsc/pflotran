#include "finclude/petscsys.h"

PetscInt, parameter :: MAXHEADERLENGTH = 4096
PetscInt, parameter :: MAXSTRINGLENGTH = 512
PetscInt, parameter :: MAXWORDLENGTH = 32
PetscInt, parameter :: OUT_UNIT = 15
PetscInt, parameter :: OUTPUT_UNIT = 16
PetscInt, parameter :: IN_UNIT = 17
! If you increase MAX_IN_UNIT, you MUST ensure that no other units #
! lie between IN_UNIT and MAX_IN_UNIT, as these units are reserved
! for embedded input files.
PetscInt, parameter :: MAX_IN_UNIT = 25
PetscInt, parameter :: IUNIT_TEMP = 86
PetscInt, parameter :: HHISTORY_LENGTH = 1000
! HHISTORY_LENGTH is the length of the array used to store the differencing
! values h.

! formula weights
PetscReal, parameter :: FMWNACL = 58.44277d0
PetscReal, parameter :: FMWH2O = 18.01534d0  ! kg/kmol h2o
PetscReal, parameter :: FMWCO2 = 44.0098d0
PetscReal, parameter :: FMWAIR = 28.96d0
PetscReal, parameter :: FMWGLYC = 76.09d0 ! propylene glycol (C3H8O2)

! conversion factors
PetscReal, parameter :: LOG_TO_LN = 2.30258509299d0
PetscReal, parameter :: LN_TO_LOG = 0.434294481904d0  

! constants
PetscReal, parameter :: IDEAL_GAS_CONST = 8.314472d0   
PetscReal, parameter :: HEAT_OF_FUSION = 3.34d5  ! J/kg
PetscReal, parameter :: PI = 3.14159265359d0

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
PetscInt, parameter :: SURF_ONEDOF = 6

PetscInt, parameter :: GLOBAL = 1
PetscInt, parameter :: LOCAL = 2
PetscInt, parameter :: NATURAL = 3

PetscInt, parameter :: NULL_MODE = 0

! flow modes
PetscInt, parameter :: THC_MODE = 1
PetscInt, parameter :: MPH_MODE = 2
PetscInt, parameter :: RICHARDS_MODE = 3
PetscInt, parameter :: REACTIVE_TRANSPORT_MODE = 4
PetscInt, parameter :: IMS_MODE = 5
PetscInt, parameter :: FLASH2_MODE = 6
PetscInt, parameter :: G_MODE = 7
PetscInt, parameter :: MIS_MODE = 8
PetscInt, parameter :: THMC_MODE = 9

! transport modes
PetscInt, parameter :: EXPLICIT_ADVECTION = 1

! grid types
!geh: moved to grid.F90
!PetscInt, parameter :: NULL_GRID = 0
!PetscInt, parameter :: STRUCTURED_GRID = 1
!PetscInt, parameter :: UNSTRUCTURED_GRID = 2
!geh: moved to structured_grid.F90 and renumbered
!PetscInt, parameter :: CARTESIAN_GRID = 3
!PetscInt, parameter :: CYLINDRICAL_GRID = 4
!PetscInt, parameter :: SPHERICAL_GRID = 5
!geh: moved to grid.F90 and renumbered
!PetscInt, parameter :: STRUCTURED_GRID_MIMETIC = 6
!PetscInt, parameter :: UNSTRUCTURED_GRID_MIMETIC = 7
!geh: moved to unstructured_grid_aux.F90 and renumbered
!PetscInt, parameter :: TWO_DIM_GRID = 8
!PetscInt, parameter :: THREE_DIM_GRID = 9

! condition types
PetscInt, parameter :: NULL_CONDITION = 0
PetscInt, parameter :: DIRICHLET_BC = 1
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
PetscInt, parameter :: SATURATION_BC = 15
PetscInt, parameter :: DISTRIBUTED_VOLUMETRIC_RATE_SS = 16
PetscInt, parameter :: DISTRIBUTED_MASS_RATE_SS = 17
PetscInt, parameter :: WELL_SS = 100

! source/sink scaling options
PetscInt, parameter :: SCALE_BY_PERM = 1
PetscInt, parameter :: SCALE_BY_NEIGHBOR_PERM = 2
PetscInt, parameter :: SCALE_BY_VOLUME = 3

! concentration subcondition types
PetscInt, parameter :: CONSTRAINT_NULL = 0
PetscInt, parameter :: CONSTRAINT_FREE = 1
PetscInt, parameter :: CONSTRAINT_TOTAL = 2
PetscInt, parameter :: CONSTRAINT_LOG = 3
PetscInt, parameter :: CONSTRAINT_PH = 4
PetscInt, parameter :: CONSTRAINT_MINERAL = 5
PetscInt, parameter :: CONSTRAINT_GAS = 6
PetscInt, parameter :: CONSTRAINT_CHARGE_BAL = 7
PetscInt, parameter :: CONSTRAINT_TOTAL_SORB_AQ_BASED = 8
PetscInt, parameter :: CONSTRAINT_TOTAL_SORB = 9
PetscInt, parameter :: CONSTRAINT_SUPERCRIT_CO2 = 10

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
!PetscInt, parameter :: THC_MASS_RATE_DOF = 2
PetscInt, parameter :: THC_MASS_RATE_DOF = 4
PetscInt, parameter :: THC_TEMPERATURE_DOF = 2
!PetscInt, parameter :: THC_CONCENTRATION_DOF = 4
PetscInt, parameter :: THC_CONCENTRATION_DOF = 3
PetscInt, parameter :: THC_ENTHALPY_DOF = 5

PetscInt, parameter :: THMC_PRESSURE_DOF = 1
PetscInt, parameter :: THMC_TEMPERATURE_DOF = 2
PetscInt, parameter :: THMC_CONCENTRATION_DOF = 3
!PetscInt, parameter :: THMC_MASS_RATE_DOF = 4
!PetscInt, parameter :: THMC_ENTHALPY_DOF = 5
PetscInt, parameter :: THMC_DISP_X_DOF = 4
PetscInt, parameter :: THMC_DISP_Y_DOF = 5
PetscInt, parameter :: THMC_DISP_Z_DOF = 6

PetscInt, parameter :: MPH_PRESSURE_DOF = 1
PetscInt, parameter :: MPH_TEMPERATURE_DOF = 2
PetscInt, parameter :: MPH_CONCENTRATION_DOF = 3

PetscInt, parameter :: RICHARDS_PRESSURE_DOF = 1
PetscInt, parameter :: RICHARDS_CONDUCTANCE_DOF = 2

PetscInt, parameter :: MIS_PRESSURE_DOF = 1
PetscInt, parameter :: MIS_CONCENTRATION_DOF = 2

PetscInt, parameter :: GENERAL_LIQUID_PRESSURE_DOF = 1
PetscInt, parameter :: GENERAL_GAS_PRESSURE_DOF = 1
PetscInt, parameter :: GENERAL_AIR_PRESSURE_DOF = 2
PetscInt, parameter :: GENERAL_GAS_SATURATION_DOF = 3
PetscInt, parameter :: GENERAL_LIQUID_FLUX_DOF = 1
PetscInt, parameter :: GENERAL_GAS_FLUX_DOF = 1
PetscInt, parameter :: GENERAL_TEMPERATURE_DOF = 3
PetscInt, parameter :: GENERAL_MOLE_FRACTION_DOF = 2
PetscInt, parameter :: GENERAL_LIQUID_CONDUCTANCE_DOF = -1
PetscInt, parameter :: GENERAL_GAS_CONDUCTANCE_DOF = -2
PetscInt, parameter :: GENERAL_FLUX_DOF = 4

! output definitions
PetscInt, parameter :: X_COORDINATE =             1
PetscInt, parameter :: Y_COORDINATE =             2
PetscInt, parameter :: Z_COORDINATE =             3
PetscInt, parameter :: TEMPERATURE =              4
PetscInt, parameter :: LIQUID_PRESSURE =          5
PetscInt, parameter :: LIQUID_SATURATION =        6
PetscInt, parameter :: GAS_SATURATION =           7
PetscInt, parameter :: LIQUID_DENSITY =           8
PetscInt, parameter :: LIQUID_DENSITY_MOL =       9
PetscInt, parameter :: GAS_DENSITY =             10
PetscInt, parameter :: GAS_DENSITY_MOL =         11
PetscInt, parameter :: LIQUID_ENERGY =           12
PetscInt, parameter :: GAS_ENERGY =              13
PetscInt, parameter :: LIQUID_VISCOSITY =        14
PetscInt, parameter :: GAS_VISCOSITY =           15
PetscInt, parameter :: LIQUID_MOBILITY =         16
PetscInt, parameter :: GAS_MOBILITY =            17
PetscInt, parameter :: LIQUID_MOLE_FRACTION =    18
PetscInt, parameter :: GAS_MOLE_FRACTION =       19
PetscInt, parameter :: POROSITY =                20
PetscInt, parameter :: PHASE =                   21
PetscInt, parameter :: MATERIAL_ID =             22

PetscInt, parameter :: PRIMARY_MOLALITY =        23
PetscInt, parameter :: SECONDARY_MOLALITY =      24
PetscInt, parameter :: TOTAL_MOLALITY =          25
PetscInt, parameter :: PRIMARY_MOLARITY =        26
PetscInt, parameter :: SECONDARY_MOLARITY =      27
PetscInt, parameter :: TOTAL_MOLARITY =          28
PetscInt, parameter :: MINERAL_VOLUME_FRACTION = 29
PetscInt, parameter :: MINERAL_RATE =            30
PetscInt, parameter :: MINERAL_SURFACE_AREA =    31
PetscInt, parameter :: MINERAL_SATURATION_INDEX =32
PetscInt, parameter :: PH =                      33
PetscInt, parameter :: SURFACE_CMPLX =           34
PetscInt, parameter :: SURFACE_CMPLX_FREE =      35
PetscInt, parameter :: SURFACE_SITE_DENSITY =    36
PetscInt, parameter :: KIN_SURFACE_CMPLX =       37
PetscInt, parameter :: KIN_SURFACE_CMPLX_FREE =  38
PetscInt, parameter :: PRIMARY_ACTIVITY_COEF =   39
PetscInt, parameter :: SECONDARY_ACTIVITY_COEF = 40
PetscInt, parameter :: SC_FUGA_COEFF =           41
PetscInt, parameter :: PRIMARY_KD =              42
PetscInt, parameter :: TOTAL_SORBED =            43
PetscInt, parameter :: TOTAL_SORBED_MOBILE =     44
PetscInt, parameter :: COLLOID_MOBILE =          45
PetscInt, parameter :: COLLOID_IMMOBILE =        46
PetscInt, parameter :: AGE =                     47
PetscInt, parameter :: STATE =                   48
PetscInt, parameter :: PROCESSOR_ID =            49
PetscInt, parameter :: ICE_SATURATION =          50
PetscInt, parameter :: TOTAL_BULK =              51
PetscInt, parameter :: ICE_DENSITY =             52
PetscInt, parameter :: GAS_PRESSURE =            53
PetscInt, parameter :: SECONDARY_TEMPERATURE =   54
PetscInt, parameter :: SECONDARY_CONCENTRATION = 55
PetscInt, parameter :: SEC_MIN_VOLFRAC =         56

PetscInt, parameter :: SURFACE_FLOW_PRESSURE =   57

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
PetscInt, parameter :: TECPLOT_FEBRICK_FORMAT = 3
PetscInt, parameter :: TECPLOT_FEQUADRILATERAL_FORMAT = 4

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

! grid cell type
PetscInt, parameter :: HEX_TYPE          = 1
PetscInt, parameter :: TET_TYPE          = 2
PetscInt, parameter :: WEDGE_TYPE        = 3
PetscInt, parameter :: PYR_TYPE          = 4
! 2D cell types:
PetscInt, parameter :: TRI_TYPE          = 5
PetscInt, parameter :: QUAD_TYPE         = 6

! grid cell properties
PetscInt, parameter :: LINE_FACE_TYPE    = 1
PetscInt, parameter :: TRI_FACE_TYPE     = 2
PetscInt, parameter :: QUAD_FACE_TYPE    = 3
PetscInt, parameter :: MAX_VERT_PER_FACE = 4
PetscInt, parameter :: MAX_FACE_PER_CELL = 6

! ids of non-petsc arrays
PetscInt, parameter :: MATERIAL_ID_ARRAY = 1
PetscInt, parameter :: SATURATION_FUNCTION_ID_ARRAY = 2

! interpolation methods
PetscInt, parameter :: INTERPOLATION_NULL = 0
PetscInt, parameter :: INTERPOLATION_STEP = 1
PetscInt, parameter :: INTERPOLATION_LINEAR = 2

! surface/subsurface flags
PetscInt, parameter :: SUBSURFACE = 0
PetscInt, parameter :: SURFACE    = 1

PetscInt, parameter :: DECOUPLED     = 0
PetscInt, parameter :: SEQ_COUPLED   = 1
PetscInt, parameter :: FULLY_COUPLED = 2

PetscInt, parameter :: KINEMATIC_WAVE = 1
PetscInt, parameter :: DIFFUSION_WAVE = 2

PetscInt, parameter :: TWO_POINT_FLUX = 0
PetscInt, parameter :: LSM_FLUX       = 1

! print secondary continuum variable ids
PetscInt, parameter :: PRINT_SEC_TEMP =           0
PetscInt, parameter :: PRINT_SEC_CONC =           1
PetscInt, parameter :: PRINT_SEC_MIN_VOLFRAC =    2
#define HASH
