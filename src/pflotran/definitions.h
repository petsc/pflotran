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
PetscInt, parameter :: XY_DIRECTION = 4
PetscInt, parameter :: XZ_DIRECTION = 5
PetscInt, parameter :: YZ_DIRECTION = 6
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
PetscInt, parameter :: TH_MODE = 10

! transport modes
PetscInt, parameter :: EXPLICIT_ADVECTION = 1

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
PetscInt, parameter :: HET_VOL_RATE_SS = 16
PetscInt, parameter :: HET_MASS_RATE_SS = 17
PetscInt, parameter :: HET_DIRICHLET = 18
PetscInt, parameter :: WELL_SS = 100

! source/sink scaling options
PetscInt, parameter :: SCALE_BY_PERM = 1
PetscInt, parameter :: SCALE_BY_NEIGHBOR_PERM = 2
PetscInt, parameter :: SCALE_BY_VOLUME = 3

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

PetscInt, parameter :: TH_PRESSURE_DOF = 1
PetscInt, parameter :: TH_MASS_RATE_DOF = 4
PetscInt, parameter :: TH_TEMPERATURE_DOF = 2
PetscInt, parameter :: TH_CONCENTRATION_DOF = 3
PetscInt, parameter :: TH_ENTHALPY_DOF = 5

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

! mphase equation of state
PetscInt, parameter :: EOS_SPAN_WAGNER = 1
PetscInt, parameter :: EOS_MRK = 2

! phase ids
PetscInt, parameter :: LIQUID_PHASE = 1
PetscInt, parameter :: GAS_PHASE = 2

! approaches to coupling reactive transport
PetscInt, parameter :: GLOBAL_IMPLICIT = 0
PetscInt, parameter :: OPERATOR_SPLIT = 1

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
