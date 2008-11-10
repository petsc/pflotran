module Transport_Output_module

  implicit none
  
  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"
#include "finclude/petsclog.h"

  PetscInt, parameter :: X_COORDINATE = 1
  PetscInt, parameter :: Y_COORDINATE = 2
  PetscInt, parameter :: Z_COORDINATE = 3

#if 0
  PetscInt, parameter :: TEMPERATURE = 4
  PetscInt, parameter :: PRESSURE = 5
  PetscInt, parameter :: LIQUID_SATURATION = 6
  PetscInt, parameter :: GAS_SATURATION = 7
  PetscInt, parameter :: LIQUID_ENERGY = 8
  PetscInt, parameter :: GAS_ENERGY = 9
  PetscInt, parameter :: LIQUID_MOLE_FRACTION = 10
  PetscInt, parameter :: GAS_MOLE_FRACTION = 11
  PetscInt, parameter :: VOLUME_FRACTION = 12
  PetscInt, parameter :: PHASE = 13
  PetscInt, parameter :: MATERIAL_ID = 14

  PetscInt, parameter :: TECPLOT_INTEGER = 0
  PetscInt, parameter :: TECPLOT_REAL = 1

  PetscInt, parameter :: TECPLOT_FILE = 0
  PetscInt, parameter ::  HDF5_FILE = 1

  PetscInt, parameter :: LIQUID_PHASE = 1
  PetscInt, parameter :: GAS_PHASE = 2
#endif
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  
  public :: RTOutput
  
contains

! ************************************************************************** !
!
! RTOutput: Main driver for all output subroutines
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine RTOutput(realization)

  use Transport_Realization_module
  
  implicit none
  
  type(tr_realization_type) :: realization

  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend
  
  if (realization%output_option%print_hdf5) then
    call PetscGetTime(tstart,ierr) 
!    call OutputHDF5(realization)
    call PetscGetTime(tend,ierr) 
    if (realization%option%myrank == 0) &
      print *, '      Seconds to write to HDF5 file: ', (tend-tstart)
  endif
 
  if (realization%output_option%print_tecplot) then
    call PetscGetTime(tstart,ierr) 
!    call OutputTecplot(realization)
    call PetscGetTime(tend,ierr) 
    if (realization%option%myrank == 0) &
      print *, '      Seconds to write to Tecplot file(s): ', (tend-tstart)
  endif
  
  realization%output_option%plot_number = realization%output_option%plot_number + 1

end subroutine RTOutput
#if 0
! ************************************************************************** !
!
! Output_Breakthrough: Main driver for all breakthrough output subroutines
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine OutputBreakthrough(realization)
                           ! for some flakey reason, current Intel 10.1 reports
                           ! error if 'only' statement not used.
  use Realization_module, only : realization_type 
  use Option_Module
  
  implicit none
  
  type(realization_type) :: realization

#if 0  
  if (realization%output_option%print_hdf5) then
    call OutputHDF5(realization)
  endif
#endif
 
  if (realization%output_option%print_tecplot) then
    call OutputBreakthroughTecplot(realization)
  endif

end subroutine OutputBreakthrough


! ************************************************************************** !
!
! OutputTecplot: Print to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !  
subroutine OutputTecplot(realization)

  use Realization_module
  use Grid_module
  use Structured_Grid_module
  use Option_module
  use Field_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: i
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(output_option_type), pointer :: output_option
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec
  Vec :: natural_vec
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '.tec'
  else
    if (output_option%plot_number < 10) then
      write(filename,'("pflow00",i1,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 100) then
      write(filename,'("pflow0",i2,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 1000) then
      write(filename,'("pflow",i3,".tec")') output_option%plot_number  
    else if (output_option%plot_number < 10000) then
      write(filename,'("pflow",i4,".tec")') output_option%plot_number  
    endif
  endif
  
  if (option%myrank == 0) then
    print *, '--> write tecplot output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 option%time/output_option%tconv,output_option%tunit
    ! write variables
    select case(option%imode)
      case (IMS_MODE)
        string = 'VARIABLES=' // &
                 '"X [m]",' // &
                 '"Y [m]",' // &
                 '"Z [m]",' // &
                 '"T [C]",' // &
                 '"P [Pa]",' // &
                 '"Phase",' // &
                 '"sg",' // &
                 '"Den(l)",' // &
                 '"Den(g)",' // &
                 '"Ul",' // &
                 '"Ug",'
        if (option%rk > 0.d0) then
          string = trim(string) // '"Volume Fraction"'
        endif
         if (associated(field%imat)) then
          string = trim(string) // ',"Material_ID"'
        endif
      case (MPH_MODE)
        string = 'VARIABLES=' // &
                 '"X [m]",' // &
                 '"Y [m]",' // &
                 '"Z [m]",' // &
                 '"T [C]",' // &
                 '"P [Pa]",' // &
                 '"Phase",' // &
                 '"sg",' // &
                 '"Den(l)",' // &
                 '"Den(g)",' // &
                 '"Ul",' // &
                 '"Ug",'
        do i=1,option%nspec
          write(string2,'(''"Xl('',i2,'')",'')') i
          string = trim(string) // trim(string2)
        enddo
        do i=1,option%nspec
          write(string2,'(''"Xg('',i2,'')",'')') i
          string = trim(string) // trim(string2)
        enddo
        if (option%rk > 0.d0) then
          string = trim(string) // '"Volume Fraction"'
        endif
         if (associated(field%imat)) then
          string = trim(string) // ',"Material_ID"'
        endif
      case(RICHARDS_MODE,RICHARDS_LITE_MODE)
        if (option%imode == RICHARDS_MODE) then
          string = 'VARIABLES=' // &
                   '"X [m]",' // &
                   '"Y [m]",' // &
                   '"Z [m]",' // &
                   '"T [C]",' // &
                   '"P [Pa]",' // &
                   '"sl",' // &
                   '"Ul"' 
        else
          string = 'VARIABLES=' // &
                   '"X [m]",' // &
                   '"Y [m]",' // &
                   '"Z [m]",' // &
                   '"P [Pa]",'
        endif
        if (option%imode == RICHARDS_MODE) then
          do i=1,option%nspec
            write(string2,'('',"Xl('',i2,'')"'')') i
            string = trim(string) // trim(string2)
          enddo
        endif
        if (option%rk > 0.d0) then
          string = trim(string) // ',"Volume Fraction"'
        endif
        if (associated(field%imat)) then
          string = trim(string) // ',"Material_ID"'
        endif
      case default
        string = 'VARIABLES=' // &
                 '"X [m]",' // &
                 '"Y [m]",' // &
                 '"Z [m]",' // &
                 '"T [C]",' // &
                 '"P [Pa]",' // &
                 '"sl",' // &
                 '"C [mol/L]"'
        if (option%rk > 0.d0) then
          string = trim(string) // ',"Volume Fraction"'
        endif
    end select
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                 &'', K='',i4,'','')') &
                 option%time/output_option%tconv,grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL)  
  call GridCreateVector(grid,ONEDOF,natural_vec,NATURAL)  

  ! write out coorindates
  call GetCoordinates(grid,global_vec,X_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCoordinates(grid,global_vec,Y_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCoordinates(grid,global_vec,Z_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  select case(option%imode)
    case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE,IMS_MODE)

      ! temperature
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! pressure
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,PRESSURE,ZERO_INTEGER)
          call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! liquid saturation
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select

      ! gas saturation
      select case(option%imode)
        case(MPH_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select
    
      ! liquid energy
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select
    
     ! gas energy
      select case(option%imode)
        case(MPH_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      end select
    
     ! liquid mole fraction
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE)
          ! liquid mole fractions
          do i=1,option%nspec
            call GetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i-1)
            call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
            call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
          enddo
      end select
   
     ! gas mole fraction
      select case(option%imode)
        case(MPH_MODE)
          ! gas mole fractions
          do i=1,option%nspec
            call GetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i-1)
            call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
            call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
          enddo
      end select 
      
      ! Volume Fraction
      if (option%rk > 0.d0) then
        call GetVarFromArray(realization,global_vec,VOLUME_FRACTION,ZERO_INTEGER)
        call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      endif
    
      ! phase
      select case(option%imode)
        case(MPH_MODE)
          call GetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
          call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_INTEGER)
      end select
      
      ! material id
      if (associated(field%imat)) then
        call GetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
        call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_INTEGER)
      endif
  
    case default
  
      ! temperature
      call GridGlobalToNatural(grid,field%temp,natural_vec,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

      ! pressure
      call GridGlobalToNatural(grid,field%pressure,natural_vec,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

      ! saturation
      call GridGlobalToNatural(grid,field%sat,natural_vec,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

      ! concentration
      call GridGlobalToNatural(grid,field%conc,natural_vec,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

      ! volume fraction
      if (option%rk > 0.d0) then
        call GetVarFromArray(realization,global_vec,VOLUME_FRACTION,ZERO_INTEGER)
        call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
      endif
    
  end select
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  close(IUNIT3)
  
  if (output_option%print_tecplot_velocities) then
    call OutputVelocitiesTecplot(realization)
  endif
  
  if (output_option%print_tecplot_flux_velocities) then
    if (grid%structured_grid%nx > 1) then
      call OutputFluxVelocitiesTecplot(realization,LIQUID_PHASE, &
                                       X_DIRECTION)
      select case(option%imode)
        case(MPH_MODE,IMS_MODE)
          call OutputFluxVelocitiesTecplot(realization,GAS_PHASE, &
                                           X_DIRECTION)
      end select
    endif
    if (grid%structured_grid%ny > 1) then
      call OutputFluxVelocitiesTecplot(realization,LIQUID_PHASE, &
                                       Y_DIRECTION)
      select case(option%imode)
        case(MPH_MODE,IMS_MODE)
          call OutputFluxVelocitiesTecplot(realization,GAS_PHASE, &
                                           Y_DIRECTION)
      end select
    endif
    if (grid%structured_grid%nz > 1) then
      call OutputFluxVelocitiesTecplot(realization,LIQUID_PHASE, &
                                       Z_DIRECTION)
      select case(option%imode)
        case(MPH_MODE,IMS_MODE)
          call OutputFluxVelocitiesTecplot(realization,GAS_PHASE, &
                                           Z_DIRECTION)
      end select
    endif
  endif
      
end subroutine OutputTecplot

! ************************************************************************** !
!
! OutputVelocitiesTecplot: Print velocities to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputVelocitiesTecplot(realization)
 
  use Realization_module
  use Grid_module
  use Option_module
  use Field_module
  
  implicit none

  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(output_option_type), pointer :: output_option
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  Vec :: natural_vec

  PetscReal, pointer :: vec_ptr(:)
  
  grid => realization%grid
  field => realization%field
  option => realization%option
  output_option => realization%output_option
  
  ! open file
  if (output_option%plot_number < 10) then
    write(filename,'("pflow_vel00",i1,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 100) then
    write(filename,'("pflow_vel0",i2,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 1000) then
    write(filename,'("pflow_vel",i3,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 10000) then
    write(filename,'("pflow_vel",i4,".tec")') output_option%plot_number  
  endif
  
  if (option%myrank == 0) then
    print *, '--> write tecplot velocity output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 option%time/output_option%tconv,output_option%tunit
    ! write variables
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]",' // &
             '"vlx [m/y]",' // &
             '"vly [m/y]",' // &
             '"vlz [m/y]"'
    if (option%nphase > 1) then
      string = trim(string) // &
               ',"vgx [m/y]",' // &
               '"vgy [m/y]",' // &
               '"vgz [m/y]"'
    endif
        if (associated(field%imat)) then
          string = trim(string) // ',"Material_ID"'
        endif
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                 &'', K='',i4,'','')') &
                 option%time/output_option%tconv, &
                 grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL)  
  call GridCreateVector(grid,ONEDOF,natural_vec,NATURAL)    

  ! write out coorindates
  call GetCoordinates(grid,global_vec,X_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCoordinates(grid,global_vec,Y_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCoordinates(grid,global_vec,Z_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

  if (option%nphase > 1) then
    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
    call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
    call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)

    call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
    call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_REAL)
  endif

  ! material id
  if (associated(field%imat)) then
    call GetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
    call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(IUNIT3,realization,natural_vec,TECPLOT_INTEGER)
  endif
  
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  close(IUNIT3)
  
end subroutine OutputVelocitiesTecplot

! ************************************************************************** !
!
! OutputFluxVelocitiesTecplot: Print intercellular fluxes to Tecplot file in 
!                              BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputFluxVelocitiesTecplot(realization,iphase, &
                                       direction)
!geh - specifically, the flow velocities at the interfaces between cells
 
  use Realization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  
  implicit none

  type(realization_type) :: realization
  PetscInt :: iphase
  PetscInt :: direction
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(output_option_type), pointer :: output_option
  
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscInt :: local_size, global_size
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscInt :: i, j, k
  PetscInt :: local_id, ghosted_id
  PetscInt :: adjusted_size
  PetscInt :: count, iconn
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: array(:)
  PetscInt, allocatable :: indices(:)
  Vec :: global_vec

  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
    
  nullify(array)
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  ! open file
  filename = 'pflow_'
  
  select case(iphase)
    case(LIQUID_PHASE)
      filename = trim(filename) // 'vl'
    case(GAS_PHASE)
      filename = trim(filename) // 'vg'
  end select
  
  select case(direction)
    case(X_DIRECTION)
      filename = trim(filename) // 'x'
    case(Y_DIRECTION)
      filename = trim(filename) // 'y'
    case(Z_DIRECTION)
      filename = trim(filename) // 'z'
  end select 
  
  if (output_option%plot_number < 10) then
    write(string,'("00",i1,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 100) then
    write(string,'("0",i2,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 1000) then
    write(string,'(i3,".tec")') output_option%plot_number  
  else if (output_option%plot_number < 10000) then
    write(string,'(i4,".tec")') output_option%plot_number  
  endif
  
  filename = trim(filename) // trim(string)
  
  if (option%myrank == 0) then
    print *, '--> write tecplot velocity flux output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 option%time/output_option%tconv,output_option%tunit
    ! write variables
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]",'
    select case(iphase)
      case(LIQUID_PHASE)
        string = trim(string) // '"Liquid'
      case(GAS_PHASE)
        string = trim(string) // '"Gas'
    end select
  
    select case(direction)
      case(X_DIRECTION)
        string = trim(string) // ' vlx [m/y]"'
      case(Y_DIRECTION)
        string = trim(string) // ' vly [m/y]"'
      case(Z_DIRECTION)
        string = trim(string) // ' vlz [m/y]"'
    end select 
    
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    select case(direction)
      case(X_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%time/output_option%tconv,grid%structured_grid%nx-1,grid%structured_grid%ny,grid%structured_grid%nz 
      case(Y_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%time/output_option%tconv,grid%structured_grid%nx,grid%structured_grid%ny-1,grid%structured_grid%nz 
      case(Z_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%time/output_option%tconv,grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz-1
    end select 
  
  
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  
  ! face coordinates
  local_size = grid%nlmax
  global_size = grid%nmax
!GEH - Structured Grid Dependence - Begin
  nx_local = grid%structured_grid%nlx
  ny_local = grid%structured_grid%nly
  nz_local = grid%structured_grid%nlz
  nx_global = grid%structured_grid%nx
  ny_global = grid%structured_grid%ny
  nz_global = grid%structured_grid%nz
  select case(direction)
    case(X_DIRECTION)
      global_size = grid%nmax-grid%structured_grid%ny*grid%structured_grid%nz
      nx_global = grid%structured_grid%nx-1
      if (grid%structured_grid%ngxe-grid%structured_grid%nxe == 0) then
        local_size = grid%nlmax-grid%structured_grid%nlyz
        nx_local = grid%structured_grid%nlx-1
      endif
    case(Y_DIRECTION)
      global_size = grid%nmax-grid%structured_grid%nx*grid%structured_grid%nz
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
        local_size = grid%nlmax-grid%structured_grid%nlxz
        ny_local = grid%structured_grid%nly-1
      endif
    case(Z_DIRECTION)
      global_size = grid%nmax-grid%structured_grid%nxy
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
        local_size = grid%nlmax-grid%structured_grid%nlxy
        nz_local = grid%structured_grid%nlz-1
      endif
  end select  
  allocate(indices(local_size))

  ! fill indices array with natural ids in newly sized array
  count = 0
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        indices(count) = i+grid%structured_grid%nxs+(j-1+grid%structured_grid%nys)*nx_global+ &
                         (k-1+grid%structured_grid%nzs)*nx_global*ny_global
      enddo
    enddo
  enddo
  
  ! X-coordinates
  call VecGetArrayF90(grid%structured_grid%dx,vec_ptr,ierr)
  count = 0
  allocate(array(local_size))
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        array(count) = grid%x(grid%nL2G(local_id))
        if (direction == X_DIRECTION) &
          array(count) = array(count) + 0.5d0*vec_ptr(local_id)
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(grid%structured_grid%dx,vec_ptr,ierr)
  ! warning: adjusted size will be changed in ConvertArrayToNatural
  ! thus, you cannot pass in local_size, since it is needed later
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
  ! since the array has potentially been resized, must reallocate
  deallocate(array)
  nullify(array)

  ! Y-coordinates
  count = 0
  allocate(array(local_size))
  call VecGetArrayF90(grid%structured_grid%dy,vec_ptr,ierr)
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        array(count) = grid%y(grid%nL2G(local_id))
        if (direction == Y_DIRECTION) &
          array(count) = array(count) + 0.5d0*vec_ptr(local_id)
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(grid%structured_grid%dy,vec_ptr,ierr)
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)

  ! Z-coordinates
  count = 0
  allocate(array(local_size))
  call VecGetArrayF90(grid%structured_grid%dz,vec_ptr,ierr)
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        array(count) = grid%z(grid%nL2G(local_id))
        if (direction == Z_DIRECTION) &
          array(count) = array(count) + 0.5d0*vec_ptr(local_id)
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(grid%structured_grid%dz,vec_ptr,ierr)
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)

  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL) 
  call VecZeroEntries(global_vec,ierr)
  call VecGetArrayF90(global_vec,vec_ptr,ierr)
  
  ! place interior velocities in a vector
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      ghosted_id = cur_connection_set%id_up(iconn)
      local_id = grid%nG2L(ghosted_id) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      if (local_id <= 0 .or. &
          abs(cur_connection_set%dist(direction,iconn)) < 0.99d0) cycle
      vec_ptr(local_id) = field%internal_velocities(iphase,iconn)
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! write out data set
  count = 0
  allocate(array(local_size))
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        array(count) = vec_ptr(local_id)
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(global_vec,vec_ptr,ierr)
  
  call VecDestroy(global_vec,ierr)
  
!GEH - Structured Grid Dependence - End
  
  array(1:local_size) = array(1:local_size)*output_option%tconv ! convert time units
  
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,realization,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)
  
  deallocate(indices)

  close(IUNIT3)
  
end subroutine OutputFluxVelocitiesTecplot

! ************************************************************************** !
!
! OutputVectorTecplot: Print a vector to a Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputVectorTecplot(filename,dataset_name,realization,vector)
 
  use Realization_module
  use Option_module
  use Field_module
  use Grid_module
  
  implicit none

  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXNAMELENGTH) :: dataset_name
  type(realization_type) :: realization
  Vec :: vector

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt, parameter :: fid=86

  option => realization%option
  grid => realization%grid
  field => realization%field
  
  ! open file
  if (option%myrank == 0) then
    print *, '--> write tecplot output file: ', filename
    open(unit=fid,file=filename,action="write")
  
    ! write header
    ! write title
    write(fid,'(''TITLE = "PFLOTRAN Vector"'')')
    ! write variables
    string = 'VARIABLES=' // &
             '"X [m]",' // &
             '"Y [m]",' // &
             '"Z [m]",'
    string = trim(string) // '"' // trim(dataset_name) // '"'
    if (associated(field%imat)) then
      string = trim(string) // ',"Material_ID"'
    endif
    write(fid,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',a,''",'','' I='',i4,'', J='',i4, &
                 &'', K='',i4,'','')') trim(dataset_name), &
                 grid%structured_grid%nx,grid%structured_grid%ny, &
                 grid%structured_grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(fid,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL)  
  call GridCreateVector(grid,ONEDOF,natural_vec,NATURAL)    

  ! write out coorindates
  call GetCoordinates(grid,global_vec,X_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  call GetCoordinates(grid,global_vec,Y_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  call GetCoordinates(grid,global_vec,Z_COORDINATE)
  call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  call GridGlobalToNatural(grid,vector,natural_vec,ONEDOF)
  call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_REAL)

  if (associated(field%imat)) then
    call GetVarFromArray(realization,global_vec,MATERIAL_ID,ZERO_INTEGER)
    call GridGlobalToNatural(grid,global_vec,natural_vec,ONEDOF)
    call WriteTecplotDataSetFromVec(fid,realization,natural_vec,TECPLOT_INTEGER)
  endif
  
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  close(fid)
  
end subroutine OutputVectorTecplot

! ************************************************************************** !
!
! WriteTecplotDataSetFromVec: Writes data from a Petsc Vec within a block
!                             of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSetFromVec(fid,realization,vec,datatype)

  use Realization_module
  
  implicit none

  PetscInt :: fid
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: datatype
  
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,realization,vec_ptr,datatype,ZERO_INTEGER) ! 0 implies grid%nlmax
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteTecplotDataSetFromVec

! ************************************************************************** !
!
! WriteTecplotDataSet: Writes data from an array within a block
!                      of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSet(fid,realization,array,datatype,size_flag)

  use Realization_module
  use Grid_module
  use Option_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  PetscReal :: array(:)
  PetscInt, save :: max_local_size_saved = -1
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscInt :: i
  PetscMPIInt :: iproc, recv_size
  PetscInt :: max_local_size, local_size
  PetscInt :: istart, iend, num_in_array
  PetscInt :: status(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)
  
  grid => realization%grid
  option => realization%option
  
  if (size_flag /= 0) then
    call MPI_Allreduce(size_flag,max_local_size,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                       PETSC_COMM_WORLD,ierr)
    local_size = size_flag
  else 
  ! if first time, determine the maximum size of any local array across 
  ! all procs
    if (max_local_size_saved < 0) then
      call MPI_Allreduce(grid%nlmax,max_local_size,ONE_INTEGER,MPI_INTEGER,MPI_MAX, &
                         PETSC_COMM_WORLD,ierr)
      max_local_size_saved = max_local_size
      if (option%myrank == 0) print *, 'max_local_size_saved: ', max_local_size
    endif
    max_local_size = max_local_size_saved
    local_size = grid%nlmax
  endif
  
  ! transfer the data to an integer or real array
  if (datatype == TECPLOT_INTEGER) then
    allocate(integer_data(max_local_size+10))
    allocate(integer_data_recv(max_local_size))
    do i=1,local_size
      integer_data(i) = int(array(i))
    enddo
  else
    allocate(real_data(max_local_size+10))
    allocate(real_data_recv(max_local_size))
    do i=1,local_size
      real_data(i) = array(i)
    enddo
  endif
  
  ! communicate data to processor 0, round robin style
  if (option%myrank == 0) then
    if (datatype == TECPLOT_INTEGER) then
      ! This approach makes output files identical, regardless of processor
      ! distribution.  It is necessary when diffing files.
      iend = 0
      do
        istart = iend+1
        if (iend+10 > local_size) exit
        iend = istart+9
        write(fid,'(10(i3,x))') integer_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      integer_data(1:local_size-iend) = integer_data(iend+1:local_size)
      num_in_array = local_size-iend
    else
      iend = 0
      do
        istart = iend+1
        if (iend+10 > local_size) exit
        iend = istart+9
        write(fid,'(10(es11.4,x))') real_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      real_data(1:local_size-iend) = real_data(iend+1:local_size)
      num_in_array = local_size-iend
    endif
    do iproc=1,option%commsize-1
      call MPI_Probe(iproc,MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
      recv_size = status(MPI_TAG)
      if (datatype == 0) then
        call MPI_Recv(integer_data_recv,recv_size,MPI_INTEGER,iproc, &
                      MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
        if (recv_size > 0) then
          integer_data(num_in_array+1:num_in_array+recv_size) = &
                                             integer_data_recv(1:recv_size)
          num_in_array = num_in_array+recv_size
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+10 > num_in_array) exit
          iend = istart+9
          write(fid,'(10(i3,x))') integer_data(istart:iend)
        enddo
        if (iend > 0) then
          integer_data(1:num_in_array-iend) = integer_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      else
        call MPI_Recv(real_data_recv,recv_size,MPI_DOUBLE_PRECISION,iproc, &
                      MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
        if (recv_size > 0) then
          real_data(num_in_array+1:num_in_array+recv_size) = &
                                             real_data_recv(1:recv_size)
          num_in_array = num_in_array+recv_size
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+10 > num_in_array) exit
          iend = istart+9
          write(fid,'(10(es11.4,x))') real_data(istart:iend)
        enddo
        if (iend > 0) then
          real_data(1:num_in_array-iend) = real_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      endif
    enddo
    ! Print the remaining values, if they exist
    if (datatype == 0) then
      if (num_in_array > 0) &
        write(fid,'(10(i3,x))') integer_data(1:num_in_array)
    else
      if (num_in_array > 0) &
        write(fid,'(10(es11.4,x))') real_data(1:num_in_array)
    endif
  else
    if (datatype == TECPLOT_INTEGER) then
      call MPI_Send(integer_data,local_size,MPI_INTEGER,ZERO_INTEGER,local_size, &
                    PETSC_COMM_WORLD,ierr)
    else
      call MPI_Send(real_data,local_size,MPI_DOUBLE_PRECISION,ZERO_INTEGER,local_size, &
                    PETSC_COMM_WORLD,ierr)
    endif
  endif
      
  if (datatype == TECPLOT_INTEGER) then
    deallocate(integer_data)
  else
    deallocate(real_data)
  endif

end subroutine WriteTecplotDataSet

! ************************************************************************** !
!
! OutputBreakthroughTecplot: Print to breakthrough data to TECPLOT file
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine OutputBreakthroughTecplot(realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module
  use Breakthrough_module
 
  implicit none

  type(realization_type) :: realization
  
  PetscInt :: fid, icell
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(output_option_type), pointer :: output_option
  type(breakthrough_type), pointer :: breakthrough
  logical, save :: first = .true.
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
  
  breakthrough => realization%breakthrough%first
  
  if (.not.associated(breakthrough)) return

  if (option%myrank < 10) then
    write(filename,'("breakthrough_",i1,".tec")') option%myrank  
  else if (option%myrank < 100) then
    write(filename,'("breakthrough_",i2,".tec")') option%myrank  
  else if (option%myrank < 1000) then
    write(filename,'("breakthrough_",i3,".tec")') option%myrank  
  else if (option%myrank < 10000) then
    write(filename,'("breakthrough_",i4,".tec")') option%myrank  
  else if (option%myrank < 100000) then
    write(filename,'("breakthrough_",i5,".tec")') option%myrank  
  endif
  
  ! open file
  fid = 86
  if (first) then
    first = .false.
    open(unit=fid,file=filename,action="write",status="replace")
    ! write header
    ! write title
    write(fid,'(a)',advance="no") "Time[" // trim(output_option%tunit) // "],"
    do 
      if (.not.associated(breakthrough)) exit
      do icell=1,breakthrough%region%num_cells
        call WriteBreakthroughHeaderForCell(fid,realization, &
                                            breakthrough%region%cell_ids(icell))
      enddo
      breakthrough => breakthrough%next
    enddo
    write(fid,'(a)',advance="yes") ""
    breakthrough => realization%breakthrough%first
  else
    open(unit=fid,file=filename,action="write",status="old", &
         position="append")
  endif

  write(fid,'(1es12.4)',advance="no") option%time/output_option%tconv
  do 
    if (.not.associated(breakthrough)) exit
    do icell=1,breakthrough%region%num_cells
      call WriteBreakthroughDataForCell(fid,realization, &
                                        breakthrough%region%cell_ids(icell))
    enddo
    breakthrough => breakthrough%next
  enddo
  write(fid,'(a)',advance="yes") ""
  close(fid)
      
end subroutine OutputBreakthroughTecplot

! ************************************************************************** !
!
! WriteBreakthroughHeaderForCell: Print a header for data at a cell
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine WriteBreakthroughHeaderForCell(fid,realization,local_id)

  use Realization_module
  use Grid_module
  use Field_module
  use Option_module

  implicit none
  
  PetscInt :: fid
  type(realization_type) :: realization
  PetscInt :: local_id
  
  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: cell_id_string
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  
  option => realization%option
  field => realization%field
  grid => realization%grid
  
  write(cell_id_string,*) grid%nL2A(local_id)
  cell_id_string = adjustl(cell_id_string)

  select case(option%imode)
    case (IMS_MODE)
      string = '"X [m] '// trim(cell_id_string) // '",' // &
               '"Y [m] '// trim(cell_id_string) // '",' // &
               '"Z [m] '// trim(cell_id_string) // '",' // &
               '"T [C] '// trim(cell_id_string) // '",' // &
               '"P [Pa] '// trim(cell_id_string) // '",' // &
               '"sl '// trim(cell_id_string) // '",' // &
               '"sg '// trim(cell_id_string) // '",' // &
               '"Ul '// trim(cell_id_string) // '",' // &
               '"Ug '// trim(cell_id_string) // '",'
      if (option%rk > 0.d0) then
        string = trim(string) // '"Volume Fraction '// trim(cell_id_string) // '"'
      endif
      string = trim(string) // ',"Phase '// trim(cell_id_string) // '"'
    case (MPH_MODE)
      string = '"X [m] '// trim(cell_id_string) // '",' // &
               '"Y [m] '// trim(cell_id_string) // '",' // &
               '"Z [m] '// trim(cell_id_string) // '",' // &
               '"T [C] '// trim(cell_id_string) // '",' // &
               '"P [Pa] '// trim(cell_id_string) // '",' // &
               '"sl '// trim(cell_id_string) // '",' // &
               '"sg '// trim(cell_id_string) // '",' // &
               '"Ul '// trim(cell_id_string) // '",' // &
               '"Ug '// trim(cell_id_string) // '",'
      do i=1,option%nspec
        write(string2,'(''"Xl('',i2,'') '// trim(cell_id_string) // '",'')') i
        string = trim(string) // trim(string2)
      enddo
      do i=1,option%nspec
        write(string2,'(''"Xg('',i2,'') '// trim(cell_id_string) // '",'')') i
        string = trim(string) // trim(string2)
      enddo
      if (option%rk > 0.d0) then
        string = trim(string) // '"Volume Fraction '// trim(cell_id_string) // '"'
      endif
      string = trim(string) // ',"Phase '// trim(cell_id_string) // '"'
    case(RICHARDS_MODE,RICHARDS_LITE_MODE)
      if (option%imode == RICHARDS_MODE) then
        string = '"X [m] '// trim(cell_id_string) // '",' // &
                 '"Y [m] '// trim(cell_id_string) // '",' // &
                 '"Z [m] '// trim(cell_id_string) // '",' // &
                 '"T [C] '// trim(cell_id_string) // '",' // &
                 '"P [Pa] '// trim(cell_id_string) // '",' // &
                 '"sl '// trim(cell_id_string) // '",' // &
                 '"Ul '// trim(cell_id_string) // '"' 
      else
        string = '"X [m] '// trim(cell_id_string) // '",' // &
                 '"Y [m] '// trim(cell_id_string) // '",' // &
                 '"Z [m] '// trim(cell_id_string) // '",' // &
                 '"P [Pa] '// trim(cell_id_string) // '",' // &
                 '"sl '// trim(cell_id_string) // '"'
      endif
      if (option%imode == RICHARDS_MODE) then
        do i=1,option%nspec
          write(string2,'('',"Xl('',i2,'') '// trim(cell_id_string) // '"'')') i
          string = trim(string) // trim(string2)
        enddo
      endif
      if (option%rk > 0.d0) then
        string = trim(string) // ',"Volume Fraction '// trim(cell_id_string) // '"'
      endif
    case default
      string = '"X [m]",' // &
               '"Y [m]",' // &
               '"Z [m]",' // &
               '"T [C]",' // &
               '"P [Pa]",' // &
               '"sl",' // &
               '"C [mol/L]"'
      if (option%rk > 0.d0) then
        string = trim(string) // ',"Volume Fraction"'
      endif
  end select
  write(fid,'(a)',advance="no") trim(string)

end subroutine WriteBreakthroughHeaderForCell

! ************************************************************************** !
!
! WriteBreakthroughDataForCell: Print data for data at a cell
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !  
subroutine WriteBreakthroughDataForCell(fid,realization,local_id)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none
  
  PetscInt :: fid, i
  type(realization_type) :: realization
  PetscInt :: local_id

  PetscInt :: ghosted_id
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  option => realization%option
  grid => realization%grid
  field => realization%field

100 format(es13.6)
101 format(i)
110 format(',',es13.6)
111 format(',',i)

  ghosted_id = grid%nL2G(local_id)
  ! write out coorindates
  write(fid,110,advance="no") grid%x(ghosted_id)
  write(fid,110,advance="no") grid%y(ghosted_id)
  write(fid,110,advance="no") grid%z(ghosted_id)
  
  select case(option%imode)
    case (MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE,IMS_MODE)

      ! temperature
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE,IMS_MODE)
          write(fid,110,advance="no") &
            GetVarFromArrayAtCell(realization,TEMPERATURE,ZERO_INTEGER,local_id)
      end select

      ! pressure
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE,IMS_MODE)
          write(fid,110,advance="no") &
            GetVarFromArrayAtCell(realization,PRESSURE,ZERO_INTEGER,local_id)
      end select

      ! liquid saturation
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE,IMS_MODE)
          write(fid,110,advance="no") &
            GetVarFromArrayAtCell(realization,LIQUID_SATURATION,ZERO_INTEGER,local_id)
      end select

      select case(option%imode)
        case(MPH_MODE,IMS_MODE)
          ! gas saturation
          write(fid,110,advance="no") &
            GetVarFromArrayAtCell(realization,GAS_SATURATION,ZERO_INTEGER,local_id)
      end select
    
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE)
          ! liquid energy
          write(fid,110,advance="no") &
            GetVarFromArrayAtCell(realization,LIQUID_ENERGY,ZERO_INTEGER,local_id)
      end select
    
      select case(option%imode)
        case(MPH_MODE,IMS_MODE)
          ! gas energy
          write(fid,110,advance="no") &
            GetVarFromArrayAtCell(realization,GAS_ENERGY,ZERO_INTEGER,local_id)
      end select

      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE)
          ! liquid mole fractions
          do i=1,option%nspec
            write(fid,110,advance="no") &
              GetVarFromArrayAtCell(realization,LIQUID_MOLE_FRACTION,i-1,local_id)
          enddo
      end select
  
      select case(option%imode)
        case(MPH_MODE)
          ! gas mole fractions
          do i=1,option%nspec
            write(fid,110,advance="no") &
              GetVarFromArrayAtCell(realization,GAS_MOLE_FRACTION,i-1,local_id)
          enddo
      end select 
      
      ! Volume Fraction
      if (option%rk > 0.d0) then
        write(fid,110,advance="no") &
          GetVarFromArrayAtCell(realization,VOLUME_FRACTION,ZERO_INTEGER,local_id)
      endif
    
      ! phase
      select case(option%imode)
        case(MPH_MODE,RICHARDS_MODE)
          write(fid,111,advance="no") &
            int(GetVarFromArrayAtCell(realization,PHASE,ZERO_INTEGER,local_id))
      end select
      
    case default
  
  end select

end subroutine WriteBreakthroughDataForCell

! ************************************************************************** !
!
! OutputHDF5: Print to HDF5 file in Tecplot compatible format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputHDF5(realization)

  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
  
#ifndef USE_HDF5
  implicit none
  
  type(realization_type) :: realization

  if (realization%option%myrank == 0) then
    print *
    print *, 'PFLOTRAN must be compiled with -DUSE_HDF5 to ', &
             'write to an HDF5 format.'
    print *
  endif
  stop

#else

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use hdf5
  use HDF5_module
  
  implicit none

  type(realization_type) :: realization

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  integer(HSIZE_T) :: dims(3)
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(output_option_type), pointer :: output_option
  
  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr
  
  character(len=MAXNAMELENGTH) :: filename = "pflow.h5"
  character(len=MAXSTRINGLENGTH) :: string
  logical, save :: first = .true.
  PetscReal, pointer :: array(:)
  PetscInt :: i
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option

  ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,PETSC_COMM_WORLD,MPI_INFO_NULL,hdf5_err)
#endif
  if (.not.first) call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id, &
                                      hdf5_err,prop_id)
  if (hdf5_err < 0 .or. first) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err,H5P_DEFAULT_F, &
                     prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)

  if (first) then
    if (option%myrank == 0) print *, '--> creating hdf5 output file: ', filename
  else
    if (option%myrank == 0) print *, '--> appending to hdf5 output file: ', &
                                   filename
  endif
  
  if (first) then

    ! create a group for the coordinates data set
    string = "Coordinates"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)

!GEH - Structured Grid Dependence - Begin
    ! write out coordinates in x, y, and z directions
    string = "X [m]"
    allocate(array(grid%structured_grid%nx))
    do i=1,grid%structured_grid%nx
      if (i == 1) then
        array(i) = 0.5d0*grid%structured_grid%dx0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%structured_grid%dx0(i-1)+grid%structured_grid%dx0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%nx,array,grp_id)
    deallocate(array)
  
    string = "Y [m]"
    allocate(array(grid%structured_grid%ny))
    do i=1,grid%structured_grid%ny
      if (i == 1) then
        array(i) = 0.5d0*grid%structured_grid%dy0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%structured_grid%dy0(i-1)+grid%structured_grid%dy0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%ny,array,grp_id)
    deallocate(array)
  
    string = "Z [m]"
    allocate(array(grid%structured_grid%nz))
    do i=1,grid%structured_grid%nz
      if (i == 1) then
        array(i) = 0.5d0*grid%structured_grid%dz0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%structured_grid%dz0(i-1)+grid%structured_grid%dz0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,option,grid%structured_grid%nz,array,grp_id)
    deallocate(array)
!GEH - Structured Grid Dependence - End

    call h5gclose_f(grp_id,hdf5_err)
    
  endif

  ! create a group for the data set
  write(string,'(''Time:'',es12.4,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  
  ! write out data sets 
  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL)   

  select case(option%imode)
  
    case(MPH_MODE,RICHARDS_MODE, &
         RICHARDS_LITE_MODE,IMS_MODE)

      ! temperature
      select case(option%imode)
        case (MPH_MODE,RICHARDS_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,TEMPERATURE,ZERO_INTEGER)
          string = "Temperature"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select

      ! pressure
      select case(option%imode)
        case (MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,PRESSURE,ZERO_INTEGER)
          string = "Pressure"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select

      ! liquid saturation
      select case(option%imode)
        case (MPH_MODE,RICHARDS_MODE,RICHARDS_LITE_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
          string = "Liquid Saturation"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)  
      end select

      ! gas saturation
      select case(option%imode)
        case (MPH_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,GAS_SATURATION,ZERO_INTEGER)
          string = "Gas Saturation"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      end select
      
      ! liquid energy
      select case(option%imode)
        case (MPH_MODE,RICHARDS_MODE,IMS_MODE)
          call GetVarFromArray(realization,global_vec,LIQUID_ENERGY,ZERO_INTEGER)
          string = "Liquid Energy"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
      
      ! gas energy
      select case(option%imode)
        case (MPH_MODE,IMS_MODE)    
          call GetVarFromArray(realization,global_vec,GAS_ENERGY,ZERO_INTEGER)
          string = "Gas Energy"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE) 
      end select
    
      ! liquid mole fractions
      select case(option%imode)
        case (MPH_MODE,RICHARDS_MODE)
          do i=1,option%nspec
            call GetVarFromArray(realization,global_vec,LIQUID_MOLE_FRACTION,i-1)
            write(string,'(''Liquid Mole Fraction('',i4,'')'')') i
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          enddo
      end select
      
      ! gas mole fractions
      select case(option%imode)
        case (MPH_MODE)      
          do i=1,option%nspec
            call GetVarFromArray(realization,global_vec,GAS_MOLE_FRACTION,i-1)
            write(string,'(''Gas Mole Fraction('',i4,'')'')') i
            call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
          enddo
      end select
    
      ! Volume Fraction
      if (option%rk > 0.d0) then
        call GetVarFromArray(realization,global_vec,VOLUME_FRACTION,ZERO_INTEGER)
        string = "Volume Fraction"
        call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,H5T_NATIVE_DOUBLE)
      endif
    
      ! phase
      select case(option%imode)
        case (MPH_MODE)
          call GetVarFromArray(realization,global_vec,PHASE,ZERO_INTEGER)
          string = "Phase"
          call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id,HDF_NATIVE_INTEGER) 
      end select
  
    case default
      ! temperature
      string = "Temperature"
      call HDF5WriteStructDataSetFromVec(string,realization,field%temp,grp_id, &
                                   H5T_NATIVE_DOUBLE)

      ! pressure
      string = "Pressure"
      call HDF5WriteStructDataSetFromVec(string,realization,field%pressure,grp_id, &
                                   H5T_NATIVE_DOUBLE)

      ! saturation
      string = "Saturation"
      call HDF5WriteStructDataSetFromVec(string,realization,field%sat,grp_id,H5T_NATIVE_DOUBLE)

      ! concentration
      string = "Concentration"
      call HDF5WriteStructDataSetFromVec(string,realization,field%conc,grp_id, &
                                   H5T_NATIVE_DOUBLE)

  end select
  
  if (output_option%print_hdf5_velocities) then

    ! velocities
    call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,X_DIRECTION)
    string = "Liquid X-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                 H5T_NATIVE_DOUBLE)
    call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Y_DIRECTION)
    string = "Liquid Y-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  
    call GetCellCenteredVelocities(realization,global_vec,LIQUID_PHASE,Z_DIRECTION)
    string = "Liquid Z-Velocity"
    call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                 H5T_NATIVE_DOUBLE)

    if (option%nphase > 1) then
      call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,X_DIRECTION)
      string = "Gas X-Velocity"
      call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                   H5T_NATIVE_DOUBLE)
  
      call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Y_DIRECTION)
      string = "Gas Y-Velocity"
      call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                   H5T_NATIVE_DOUBLE)
  
      call GetCellCenteredVelocities(realization,global_vec,GAS_PHASE,Z_DIRECTION)
      string = "Gas Z-Velocity"
      call HDF5WriteStructDataSetFromVec(string,realization,global_vec,grp_id, &
                                   H5T_NATIVE_DOUBLE)
    endif
  endif

  if (output_option%print_hdf5_flux_velocities) then
  
    ! internal flux velocities
    if (grid%structured_grid%nx > 1) then
      string = "Liquid X-Flux Velocities"
      call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,X_DIRECTION,grp_id)
      if (option%nphase > 1) then
        string = "Gas X-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,X_DIRECTION,grp_id)
      endif
    endif
    
    if (grid%structured_grid%ny > 1) then
      string = "Liquid Y-Flux Velocities"
      call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,Y_DIRECTION,grp_id)
      if (option%nphase > 1) then
        string = "Gas Y-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,Y_DIRECTION,grp_id)
      endif
    endif
    
    if (grid%structured_grid%nz > 1) then
      string = "Liquid Z-Flux Velocities"
      call WriteHDF5FluxVelocities(string,realization,LIQUID_PHASE,Z_DIRECTION,grp_id)
      if (option%nphase > 1) then
        string = "Gas Z-Flux Velocities"
        call WriteHDF5FluxVelocities(string,realization,GAS_PHASE,Z_DIRECTION,grp_id)
      endif
    endif
    
  endif 
  
  ! call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)

  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
  first = .false.
#endif
end subroutine OutputHDF5

#ifdef USE_HDF5
! ************************************************************************** !
!
! WriteHDF5FluxVelocities: Print flux velocities to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteHDF5FluxVelocities(name,realization,iphase,direction,file_id)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use hdf5
  use HDF5_module

  implicit none
  
  character(len=32) :: name
  type(realization_type) :: realization
  PetscInt :: iphase
  PetscInt :: direction
  integer(HID_T) :: file_id

  PetscInt :: i, j, k
  PetscInt :: count, iconn
  PetscInt :: local_id, ghosted_id
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(output_option_type), pointer :: output_option
    
  PetscReal, allocatable :: array(:)
  PetscReal, pointer :: vec_ptr(:)

  logical, save :: first = .true.
  logical, save :: trick_flux_vel_x = .false.
  logical, save :: trick_flux_vel_y = .false.
  logical, save :: trick_flux_vel_z = .false.

  Vec :: global_vec

  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
    
  grid => realization%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option  
  
  ! in a few cases (i.e. for small test problems), some processors may
  ! have no velocities to print.  This results in zero-length arrays
  ! in collective H5Dwrite().  To avoid, we switch to independent
  ! H5Dwrite() and don't write from the zero-length procs. 
!GEH - Structured Grid Dependence - Begin
  if (first) then
    nx_local = grid%structured_grid%nlx
    ny_local = grid%structured_grid%nly
    nz_local = grid%structured_grid%nlz
    if (grid%structured_grid%ngxe-grid%structured_grid%nxe == 0) then
      nx_local = grid%structured_grid%nlx-1
    endif
    call MPI_Allreduce(nx_local,i,ONE_INTEGER,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (i == 0) trick_flux_vel_x = .true.
    if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
      ny_local = grid%structured_grid%nly-1
    endif
    call MPI_Allreduce(ny_local,j,ONE_INTEGER,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (j == 0) trick_flux_vel_y = .true.
    if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
      nz_local = grid%structured_grid%nlz-1
    endif
    call MPI_Allreduce(nz_local,k,ONE_INTEGER,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (k == 0) trick_flux_vel_z = .true.
  endif

  nx_local = grid%structured_grid%nlx
  ny_local = grid%structured_grid%nly
  nz_local = grid%structured_grid%nlz
  nx_global = grid%structured_grid%nx
  ny_global = grid%structured_grid%ny
  nz_global = grid%structured_grid%nz

  select case(direction)
    case(X_DIRECTION)
      nx_global = grid%structured_grid%nx-1
      if (grid%structured_grid%ngxe-grid%structured_grid%nxe == 0) then
        nx_local = grid%structured_grid%nlx-1
      endif
      if (trick_flux_vel_x) trick_hdf5 = .true.
    case(Y_DIRECTION)
      ny_global = grid%structured_grid%ny-1
      if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
        ny_local = grid%structured_grid%nly-1
      endif
      if (trick_flux_vel_y) trick_hdf5 = .true.
    case(Z_DIRECTION)
      nz_global = grid%structured_grid%nz-1
      if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
        nz_local = grid%structured_grid%nlz-1
      endif
      if (trick_flux_vel_z) trick_hdf5 = .true.
  end select  
  allocate(array(nx_local*ny_local*nz_local))


  call GridCreateVector(grid,ONEDOF,global_vec,GLOBAL) 
  call VecZeroEntries(global_vec,ierr)
  call VecGetArrayF90(global_vec,vec_ptr,ierr)
  
  ! place interior velocities in a vector
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      ghosted_id = cur_connection_set%id_up(iconn)
      local_id = grid%nG2L(ghosted_id) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      if (local_id <= 0 .or. &
          abs(cur_connection_set%dist(direction,iconn)) < 0.99d0) cycle
      vec_ptr(local_id) = field%internal_velocities(iphase,iconn)
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  count = 0
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%structured_grid%nlx+(k-1)*grid%structured_grid%nlxy
        array(count) = vec_ptr(iphase+(direction-1)*option%nphase+ &
                               3*option%nphase*(local_id-1))
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(global_vec,vec_ptr,ierr)
  
  call VecDestroy(global_vec,ierr)
  
  array(1:nx_local*ny_local*nz_local) = &  ! convert time units
    array(1:nx_local*ny_local*nz_local) * output_option%tconv

  call HDF5WriteStructuredDataSet(name,array,file_id,H5T_NATIVE_DOUBLE, &
                        nx_global,ny_global,nz_global, &
                        nx_local,ny_local,nz_local, &
                        grid%structured_grid%nxs,grid%structured_grid%nys,grid%structured_grid%nzs)
!GEH - Structured Grid Dependence - End

  deallocate(array)
  trick_hdf5 = .false.
  first = .false.

end subroutine WriteHDF5FluxVelocities

! ************************************************************************** !
!
! WriteHDF5Coordinates: Writes structured coordinates to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteHDF5Coordinates(name,option,length,array,file_id)

  use hdf5
  use Option_module
  
  implicit none
  
  character(len=32) :: name
  type(option_type) :: option
  PetscInt :: length
  PetscReal :: array(:)
  integer(HID_T) :: file_id
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  PetscMPIInt :: rank
  
  ! write out grid structure
  rank = 1
  dims = 0
  ! x-direction
  dims(1) = length
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims)
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err) ! must be independent and only from p0
#endif
  if (option%myrank == 0) then
    call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                    hdf5_err,H5S_ALL_F,H5S_ALL_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

end subroutine WriteHDF5Coordinates
#endif

! ************************************************************************** !
!
! GetCoordinates: Extracts coordinates of cells into a PetscVec
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GetCoordinates(grid,vec,direction)

  use Grid_module
  
  implicit none
  
  type(grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  
  if (direction == X_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%x(grid%nL2G(i))
    enddo
  else if (direction == Y_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%y(grid%nL2G(i))
    enddo
  else if (direction == Z_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%z(grid%nL2G(i))
    enddo
  endif
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine GetCoordinates

! ************************************************************************** !
!
! ConvertArrayToNatural: Converts an array  to natural ordering
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine ConvertArrayToNatural(indices,array, &
                                 local_size,global_size)

  implicit none
  
  PetscInt :: local_size, global_size
  PetscInt :: indices(:)
  PetscReal, pointer :: array(:)
  
  Vec :: natural_vec
  PetscInt, allocatable :: indices_zero_based(:)
  PetscReal, pointer :: vec_ptr(:)
  
  call VecCreate(PETSC_COMM_WORLD,natural_vec,ierr)
  call VecSetSizes(natural_vec,PETSC_DECIDE,global_size,ierr)
  call VecSetType(natural_vec,VECMPI,ierr)

  allocate(indices_zero_based(local_size))
  indices_zero_based(1:local_size) = indices(1:local_size)-1

  call VecSetValues(natural_vec,local_size,indices_zero_based, &
                    array,INSERT_VALUES,ierr)

  call VecAssemblyBegin(natural_vec,ierr)
  call VecAssemblyEnd(natural_vec,ierr)

  call VecGetLocalSize(natural_vec,local_size,ierr)
  deallocate(array)
  allocate(array(local_size))
  
  call VecGetArrayF90(natural_vec,vec_ptr,ierr)
  array(1:local_size) = vec_ptr(1:local_size)
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)

  call VecDestroy(natural_vec,ierr)
  
end subroutine ConvertArrayToNatural

! ************************************************************************** !
!
! GetVarFromArrayAtCell: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function GetVarFromArrayAtCell(realization,ivar,isubvar,local_id)

  use Realization_module
  use Option_module

  use Richards_Analytical_module
  use Richards_Lite_module

  implicit none
  
  PetscReal :: GetVarFromArrayAtCell
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: local_id

  select case(realization%option%imode)
    case(RICHARDS_MODE)
      GetVarFromArrayAtCell = RichardsGetVarFromArrayAtCell(realization,ivar,isubvar,local_id)
    case(RICHARDS_LITE_MODE)
      GetVarFromArrayAtCell = RichardsLiteGetVarFromArrayAtCell(realization,ivar,isubvar,local_id)
    case(MPH_MODE)
      GetVarFromArrayAtCell = MphaseGetVarFromArrayAtCell(realization,ivar,isubvar,local_id)
    case(IMS_MODE)
      GetVarFromArrayAtCell = ImmisGetVarFromArrayAtCell(realization,ivar,isubvar,local_id)
  end select

end function GetVarFromArrayAtCell

! ************************************************************************** !
!
! GetVarFromArray: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GetVarFromArray(realization,vec,ivar,isubvar)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

#ifdef RICHARDS_ANALYTICAL
  use Richards_Analytical_module
#endif
  use Richards_Lite_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset, saturation_offset
  PetscInt :: size_var_use
  PetscInt :: size_var_node
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscReal, pointer :: var_ptr(:)
  PetscReal, pointer :: vec_ptr(:)

  option => realization%option
  grid => realization%grid
  field => realization%field

  select case(option%imode)
#ifdef RICHARDS_ANALYTICAL
    case(RICHARDS_MODE)
      call RichardsGetVarFromArray(realization,vec,ivar,isubvar)
      return
#endif  
    case(RICHARDS_LITE_MODE)
      call RichardsLiteGetVarFromArray(realization,vec,ivar,isubvar)
      return
  end select
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
      
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION)
      select case(ivar)
        case(TEMPERATURE)
          offset = 1
        case(PRESSURE)
          offset = 2
        case(LIQUID_SATURATION)
          offset = 3
        case(GAS_SATURATION)
          offset = 4
        case(LIQUID_MOLE_FRACTION)
          offset = 17+isubvar
        case(GAS_MOLE_FRACTION)
          offset = 17+option%nspec+isubvar
      end select
    
      size_var_use = 2 + 7*option%nphase + 2* option%nphase*option%nspec
      size_var_node = (option%ndof + 1) * size_var_use
        
      call VecGetArrayF90(field%var_loc,var_ptr,ierr)
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        vec_ptr(local_id) = var_ptr((ghosted_id-1)*size_var_node+offset)
      enddo
      call VecRestoreArrayF90(field%var_loc,var_ptr,ierr)

    case(LIQUID_ENERGY,GAS_ENERGY)

      select case (ivar)
        case(LIQUID_ENERGY)
          offset = 11
          saturation_offset = 3
        case(GAS_ENERGY)
          offset = 12
          saturation_offset = 4  
      end select

      size_var_use = 2 + 7*option%nphase + 2* option%nphase*option%nspec
      size_var_node = (option%ndof + 1) * size_var_use
        
      call VecGetArrayF90(field%var_loc,var_ptr,ierr)
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)      
        if (var_ptr((ghosted_id-1)*size_var_node+saturation_offset) > 1.d-30) then
          vec_ptr(local_id) = var_ptr((ghosted_id-1)*size_var_node+offset)
        else
          vec_ptr(local_id) = 0.d0
        endif
      enddo
      call VecRestoreArrayF90(field%var_loc,var_ptr,ierr)

    case(VOLUME_FRACTION)
    
      ! need to set minimum to 0.
      call VecGetArrayF90(field%phis,var_ptr,ierr)
      vec_ptr(1:grid%nlmax) = var_ptr(1:grid%nlmax)
      call VecRestoreArrayF90(field%phis,var_ptr,ierr)
     
    case(PHASE)
    
      call VecGetArrayF90(field%iphas_loc,var_ptr,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = var_ptr(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%iphas_loc,var_ptr,ierr)

    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = field%imat(grid%nL2G(local_id))
      enddo
           
  end select
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

end subroutine GetVarFromArray

! ************************************************************************** !
!
! GetCellCenteredVelocities: Computers the cell-centered velocity component 
!                            as an averages of cell face velocities
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GetCellCenteredVelocities(realization,vec,iphase,direction)

  use Realization_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: direction
  PetscInt :: iphase
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(output_option_type), pointer :: output_option
  PetscInt :: iconn
  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscReal :: velocity
  
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, allocatable :: num_additions(:)
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  output_option => realization%output_option
    
  allocate(num_additions(grid%nlmax))
  num_additions(1:grid%nlmax) = 0

  call VecSet(vec,0.d0,ierr)
  call VecGetArrayF90(vec,vec_ptr,ierr)

  ! interior velocities  
  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      velocity = field%internal_velocities(iphase,iconn)* &
                 abs(cur_connection_set%dist(direction,iconn))
      if (local_id_up > 0) then
        vec_ptr(local_id_up) = vec_ptr(local_id_up) + velocity
        num_additions(local_id_up) = num_additions(local_id_up) + 1
      endif
      if (local_id_dn > 0) then
        vec_ptr(local_id_dn) = vec_ptr(local_id_dn) + velocity
        num_additions(local_id_dn) = num_additions(local_id_dn) + 1
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! boundary velocities
  boundary_condition => realization%boundary_conditions%first
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection
    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      vec_ptr(local_id) = vec_ptr(local_id)+ &
                          field%boundary_velocities(1,iconn)* &
                          abs(cur_connection_set%dist(direction,iconn))
      num_additions(local_id) = num_additions(local_id) + 1
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! divide by number of additions
  do local_id=1,grid%nlmax
    if (num_additions(local_id) > 0) &
      vec_ptr(local_id) = vec_ptr(local_id)/real(num_additions(local_id))*output_option%tconv
  enddo
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  deallocate(num_additions)

end subroutine GetCellCenteredVelocities
#endif
end module Transport_Output_module

    
       
