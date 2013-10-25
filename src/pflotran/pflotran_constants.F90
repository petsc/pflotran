module PFLOTRAN_Constants_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "finclude/petscsys.h"

  ! MUST INCREMENT THIS NUMBER EVERYTIME A CHECKPOINT FILE IS MODIFIED TO PREVENT
  ! COMPATIBILITY ISSUES - geh.
  PetscInt, parameter, public :: CHECKPOINT_REVISION_NUMBER = 5
  
  PetscInt, parameter, public :: MAXHEADERLENGTH = 4096
  PetscInt, parameter, public :: MAXSTRINGLENGTH = 512
  PetscInt, parameter, public :: MAXWORDLENGTH = 32
  PetscInt, parameter, public :: OUT_UNIT = 15
  PetscInt, parameter, public :: OUTPUT_UNIT = 16
  PetscInt, parameter, public :: IN_UNIT = 17
  ! If you increase MAX_IN_UNIT, you MUST ensure that no other units #
  ! lie between IN_UNIT and MAX_IN_UNIT, as these units are reserved
  ! for embedded input files.
  PetscInt, parameter, public :: MAX_IN_UNIT = 25
  PetscInt, parameter, public :: IUNIT_TEMP = 86
  PetscInt, parameter, public :: HHISTORY_LENGTH = 1000
  ! HHISTORY_LENGTH is the length of the array used to store the differencing
  ! values h.
  
  ! formula weights
  PetscReal, parameter, public :: FMWNACL = 58.44277d0
  PetscReal, parameter, public :: FMWH2O = 18.01534d0  ! kg/kmol h2o
  PetscReal, parameter, public :: FMWCO2 = 44.0098d0
  PetscReal, parameter, public :: FMWAIR = 28.96d0
  PetscReal, parameter, public :: FMWGLYC = 76.09d0 ! propylene glycol (C3H8O2)
  
  ! conversion factors
  PetscReal, parameter, public :: LOG_TO_LN = 2.30258509299d0
  PetscReal, parameter, public :: LN_TO_LOG = 0.434294481904d0  
  
  ! constants
  PetscReal, parameter, public :: IDEAL_GAS_CONST = 8.314472d0   
  PetscReal, parameter, public :: HEAT_OF_FUSION = 3.34d5  ! J/kg
  PetscReal, parameter, public :: PI = 3.14159265359d0
  PetscReal, parameter, public :: Faraday = 96485.3365d0 ! C/mol
  
  PetscInt, parameter, public :: ZERO_INTEGER = 0
  PetscInt, parameter, public :: ONE_INTEGER = 1
  PetscInt, parameter, public :: TWO_INTEGER = 2
  PetscInt, parameter, public :: THREE_INTEGER = 3
  PetscInt, parameter, public :: FOUR_INTEGER = 4
  PetscInt, parameter, public :: FIVE_INTEGER = 5
  PetscInt, parameter, public :: SIX_INTEGER = 6
  PetscInt, parameter, public :: SEVEN_INTEGER = 7
  PetscInt, parameter, public :: EIGHT_INTEGER = 8
  PetscInt, parameter, public :: NINE_INTEGER = 9
  PetscInt, parameter, public :: TEN_INTEGER = 10
  PetscInt, parameter, public :: ELEVEN_INTEGER = 11
  PetscInt, parameter, public :: TWELVE_INTEGER = 12
  PetscInt, parameter, public :: NEG_ONE_INTEGER = -1
  
  PetscMPIInt, parameter, public :: ZERO_INTEGER_MPI = ZERO_INTEGER
  PetscMPIInt, parameter, public :: ONE_INTEGER_MPI = ONE_INTEGER
  PetscMPIInt, parameter, public :: TWO_INTEGER_MPI = TWO_INTEGER
  PetscMPIInt, parameter, public :: THREE_INTEGER_MPI = THREE_INTEGER
  PetscMPIInt, parameter, public :: FOUR_INTEGER_MPI = FOUR_INTEGER
  PetscMPIInt, parameter, public :: SIX_INTEGER_MPI = SIX_INTEGER
  PetscMPIInt, parameter, public :: SEVEN_INTEGER_MPI = SEVEN_INTEGER
  PetscMPIInt, parameter, public :: TWELVE_INTEGER_MPI = TWELVE_INTEGER
  PetscMPIInt, parameter, public :: MAXSTRINGLENGTH_MPI = MAXSTRINGLENGTH
  
  PetscInt, parameter, public :: X_DIRECTION = 1
  PetscInt, parameter, public :: Y_DIRECTION = 2
  PetscInt, parameter, public :: Z_DIRECTION = 3
  PetscInt, parameter, public :: XY_DIRECTION = 4
  PetscInt, parameter, public :: XZ_DIRECTION = 5
  PetscInt, parameter, public :: YZ_DIRECTION = 6
  PetscInt, parameter, public :: LOWER = 1
  PetscInt, parameter, public :: UPPER = 2
  
  PetscInt, parameter, public :: TIME_NULL = 0
  PetscInt, parameter, public :: TIME_T = 1
  PetscInt, parameter, public :: TIME_TpDT = 2
  
  PetscInt, parameter, public :: SORPTION_LINEAR = 1
  PetscInt, parameter, public :: SORPTION_LANGMUIR = 2
  PetscInt, parameter, public :: SORPTION_FREUNDLICH  = 3
  
  ! Classes
  PetscInt, parameter, public :: NULL_CLASS = 0
  PetscInt, parameter, public :: FLOW_CLASS = 1
  PetscInt, parameter, public :: TRANSPORT_CLASS = 2
  
  ! Macros that are used as 'dm_index' values.  --RTM
  PetscInt, parameter, public :: ONEDOF = 1
  PetscInt, parameter, public :: NPHASEDOF = 2
  PetscInt, parameter, public :: THREENPDOF = 3
  PetscInt, parameter, public :: NFLOWDOF = 4
  PetscInt, parameter, public :: NTRANDOF = 5
  PetscInt, parameter, public :: SURF_ONEDOF = 6
  PetscInt, parameter, public :: NGEODOF = 7
  
  PetscInt, parameter, public :: GLOBAL = 1
  PetscInt, parameter, public :: LOCAL = 2
  PetscInt, parameter, public :: NATURAL = 3
  
  PetscInt, parameter, public :: NULL_MODE = 0
  
  ! flow modes
  PetscInt, parameter, public :: THC_MODE = 1
  PetscInt, parameter, public :: MPH_MODE = 2
  PetscInt, parameter, public :: RICHARDS_MODE = 3
  PetscInt, parameter, public :: REACTIVE_TRANSPORT_MODE = 4
  PetscInt, parameter, public :: IMS_MODE = 5
  PetscInt, parameter, public :: FLASH2_MODE = 6
  PetscInt, parameter, public :: G_MODE = 7
  PetscInt, parameter, public :: MIS_MODE = 8
  PetscInt, parameter, public :: TH_MODE = 9
  
  ! transport modes
  PetscInt, parameter, public :: EXPLICIT_ADVECTION = 1
  
  ! condition types
  PetscInt, parameter, public :: NULL_CONDITION = 0
  PetscInt, parameter, public :: DIRICHLET_BC = 1
  PetscInt, parameter, public :: NEUMANN_BC = 2
  PetscInt, parameter, public :: DIRICHLET_ZERO_GRADIENT_BC = 3
  PetscInt, parameter, public :: ZERO_GRADIENT_BC = 4
  PetscInt, parameter, public :: HYDROSTATIC_BC = 5
  PetscInt, parameter, public :: SEEPAGE_BC = 6
  PetscInt, parameter, public :: MASS_RATE_SS = 7
  PetscInt, parameter, public :: VOLUMETRIC_RATE_SS = 8
  PetscInt, parameter, public :: SCALED_MASS_RATE_SS = 9
  PetscInt, parameter, public :: SCALED_VOLUMETRIC_RATE_SS = 10
  PetscInt, parameter, public :: CONCENTRATION_SS = 11
  PetscInt, parameter, public :: EQUILIBRIUM_SS = 12
  PetscInt, parameter, public :: CONDUCTANCE_BC = 13
  PetscInt, parameter, public :: UNIT_GRADIENT_BC = 14
  PetscInt, parameter, public :: SATURATION_BC = 15
  PetscInt, parameter, public :: HET_VOL_RATE_SS = 16
  PetscInt, parameter, public :: HET_MASS_RATE_SS = 17
  PetscInt, parameter, public :: HET_DIRICHLET = 18
  PetscInt, parameter, public :: ENERGY_RATE_SS = 19
  PetscInt, parameter, public :: HET_ENERGY_RATE_SS = 20
  PetscInt, parameter, public :: HET_SURF_SEEPAGE_BC = 21
  PetscInt, parameter, public :: WELL_SS = 100
  
  ! source/sink scaling options
  PetscInt, parameter, public :: SCALE_BY_PERM = 1
  PetscInt, parameter, public :: SCALE_BY_NEIGHBOR_PERM = 2
  PetscInt, parameter, public :: SCALE_BY_VOLUME = 3
  
  ! connection types
  PetscInt, parameter, public :: INTERNAL_CONNECTION_TYPE = 1
  PetscInt, parameter, public :: BOUNDARY_CONNECTION_TYPE = 2
  PetscInt, parameter, public :: INITIAL_CONNECTION_TYPE = 3
  PetscInt, parameter, public :: SRC_SINK_CONNECTION_TYPE = 4
  
  ! dofs for each mode
  PetscInt, parameter, public :: THC_PRESSURE_DOF = 1
  PetscInt, parameter, public :: THC_TEMPERATURE_DOF = 2
  PetscInt, parameter, public :: THC_CONCENTRATION_DOF = 3
  PetscInt, parameter, public :: THC_MASS_RATE_DOF = 4
  PetscInt, parameter, public :: THC_ENTHALPY_DOF = 5
  
  PetscInt, parameter, public :: TH_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TH_TEMPERATURE_DOF = 2

  PetscInt, parameter, public :: MPH_PRESSURE_DOF = 1
  PetscInt, parameter, public :: MPH_TEMPERATURE_DOF = 2
  PetscInt, parameter, public :: MPH_CONCENTRATION_DOF = 3
  
  PetscInt, parameter, public :: RICHARDS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: RICHARDS_CONDUCTANCE_DOF = 2
  
  PetscInt, parameter, public :: MIS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: MIS_CONCENTRATION_DOF = 2
  
  ! mphase equation of state
  PetscInt, parameter, public :: EOS_SPAN_WAGNER = 1
  PetscInt, parameter, public :: EOS_MRK = 2
  
  ! phase ids
  PetscInt, parameter, public :: LIQUID_PHASE = 1
  PetscInt, parameter, public :: GAS_PHASE = 2
  
  ! approaches to coupling reactive transport
  PetscInt, parameter, public :: GLOBAL_IMPLICIT = 0
  PetscInt, parameter, public :: OPERATOR_SPLIT = 1
  
  ! ids of non-petsc arrays
  PetscInt, parameter, public :: MATERIAL_ID_ARRAY = 1
  PetscInt, parameter, public :: SATURATION_FUNCTION_ID_ARRAY = 2
  
  ! interpolation methods
  PetscInt, parameter, public :: INTERPOLATION_NULL = 0
  PetscInt, parameter, public :: INTERPOLATION_STEP = 1
  PetscInt, parameter, public :: INTERPOLATION_LINEAR = 2
  
  ! surface/subsurface flags
  PetscInt, parameter, public :: SUBSURFACE = 0
  PetscInt, parameter, public :: SURFACE    = 1
  
  PetscInt, parameter, public :: DECOUPLED     = 0
  PetscInt, parameter, public :: SEQ_COUPLED   = 1
  PetscInt, parameter, public :: FULLY_COUPLED = 2
  PetscInt, parameter, public :: SEQ_COUPLED_NEW = 3
  
  PetscInt, parameter, public :: KINEMATIC_WAVE = 1
  PetscInt, parameter, public :: DIFFUSION_WAVE = 2
  
  PetscInt, parameter, public :: TWO_POINT_FLUX = 0
  PetscInt, parameter, public :: LSM_FLUX       = 1
  
  ! print secondary continuum variable ids
  PetscInt, parameter, public :: PRINT_SEC_TEMP =           0
  PetscInt, parameter, public :: PRINT_SEC_CONC =           1
  PetscInt, parameter, public :: PRINT_SEC_MIN_VOLFRAC =    2
  
  PetscInt, parameter, public :: PROCEED = 0
  PetscInt, parameter, public :: DONE = 1
  PetscInt, parameter, public :: FAIL = 2

  ! Grid type
  PetscInt, parameter, public :: NULL_GRID = 0
  PetscInt, parameter, public :: STRUCTURED_GRID = 1
  PetscInt, parameter, public :: UNSTRUCTURED_GRID = 2
  PetscInt, parameter, public :: STRUCTURED_GRID_MIMETIC = 3
  PetscInt, parameter, public :: IMPLICIT_UNSTRUCTURED_GRID = 4
  PetscInt, parameter, public :: EXPLICIT_UNSTRUCTURED_GRID = 5
  PetscInt, parameter, public :: UNSTRUCTURED_GRID_MIMETIC = 6
  PetscInt, parameter, public :: ONE_DIM_GRID = 1
  PetscInt, parameter, public :: TWO_DIM_GRID = 2
  PetscInt, parameter, public :: THREE_DIM_GRID = 3

  ! Geomechanics
  PetscInt, parameter, public :: GEOMECH_DISP_X_DOF = 1
  PetscInt, parameter, public :: GEOMECH_DISP_Y_DOF = 2
  PetscInt, parameter, public :: GEOMECH_DISP_Z_DOF = 3
  PetscInt, parameter, public :: ONE_WAY_COUPLED = 4

  ! Macros that are used as 'vscatter_index' values
  PetscInt, parameter, public :: SURF_TO_SUBSURF = 1
  PetscInt, parameter, public :: SUBSURF_TO_SURF = 2
  PetscInt, parameter, public :: SUBSURF_TO_HYDROGEOPHY = 3

end module PFLOTRAN_Constants_module
