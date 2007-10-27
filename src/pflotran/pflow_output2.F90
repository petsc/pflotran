module Output_module

  implicit none
  
  private

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petsclog.h"

#define X_COORDINATE 1
#define Y_COORDINATE 2
#define Z_COORDINATE 3
#define X_DIRECTION 1
#define Y_DIRECTION 2
#define Z_DIRECTION 3
#define TEMPERATURE 4
#define PRESSURE 5
#define LIQUID_SATURATION 6
#define GAS_SATURATION 7
#define LIQUID_ENERGY 8
#define GAS_ENERGY 9
#define LIQUID_MOLE_FRACTION 10
#define GAS_MOLE_FRACTION 11
#define VOLUME_FRACTION 12
#define PHASE 13

#define TECPLOT_INTEGER 0
#define TECPLOT_REAL 1

#define TECPLOT_FILE 0
#define HDF5_FILE 1

#define LIQUID_PHASE 1
#define GAS_PHASE 2

  integer :: hdf5_err
  logical :: trick_hdf5 = .false.
  PetscErrorCode :: ierr
  
  public :: Output, OutputTecplot, OutputHDF5
  
contains

! ************************************************************************** !
!
! Output: Main driver for all output subroutines
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine Output(solution,kplot,iplot)

  use Solution_module
  use Option_module
  
  implicit none
  
  type(solution_type) :: solution
  integer :: kplot
  integer :: iplot
  
  if (iplot == 1 .and. solution%option%print_hdf5) then
    call OutputHDF5(solution)
  endif
 
  if (iplot == 1 .and. solution%option%print_tecplot) then
    call OutputTecplot(solution,kplot)
  endif
  
  if (iplot == 1) then
    iplot = 0
    kplot = kplot + 1
  endif

end subroutine Output

! ************************************************************************** !
!
! OutputTecplot: Print to Tecplot file in BLOCK format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !  
subroutine OutputTecplot(solution,kplot)

  use Solution_module
  use Grid_module
  use Structured_Grid_module
  use Option_module
 
  implicit none

#include "definitions.h"

  type(solution_type) :: solution
  integer :: kplot
  
  integer :: i
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  Vec :: global
  Vec :: natural
  
  grid => solution%grid
  option => solution%option
  
  ! open file
  if (kplot < 10) then
    write(filename,'("pflow00",i1,".tec")') kplot  
  else if (kplot < 100) then
    write(filename,'("pflow0",i2,".tec")') kplot  
  else if (kplot < 1000) then
    write(filename,'("pflow",i3,".tec")') kplot  
  else if (kplot < 10000) then
    write(filename,'("pflow",i4,".tec")') kplot  
  endif
  
  if (option%myrank == 0) then
    print *, '--> write tecplot output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 option%t/option%tconv, option%tunit
    ! write variables
    select case(option%imode)
      case (TWOPH_MODE,MPH_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)
        string = 'VARIABLES=' // &
                 '"X-Coordinates",' // &
                 '"Y-Coordinates",' // &
                 '"Z-Coordinates",' // &
                 '"Temperature",' // &
                 '"Pressure",' // &
                 '"Liquid Saturation",' // &
                 '"Gas Saturation",' // &
                 '"Liquid Energy",' // &
                 '"Gas Energy",'
        do i=1,option%nspec
          write(string2,'(''"Liquid Mole Fraction('',i2,'')",'')') i
          string = trim(string) // trim(string2)
        enddo
        do i=1,option%nspec
          write(string2,'(''"Gas Mole Fraction('',i2,'')",'')') i
          string = trim(string) // trim(string2)
        enddo
        if (option%rk > 0.d0) then
          string = trim(string) // '"Volume Fraction",'
        endif
        string = trim(string) // '"Phase"'
      case default
        string = 'VARIABLES=' // &
                 '"X-Coordinates",' // &
                 '"Y-Coordinates",' // &
                 '"Z-Coordinates",' // &
                 '"Temperature",' // &
                 '"Pressure",' // &
                 '"Saturation",' // &
                 '"Concentration"'
        if (option%rk > 0.d0) then
          string = trim(string) // ',"Volume Fraction"'
        endif
    end select
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                 &'', K='',i4,'','')') &
                 option%t/option%tconv,grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call createPetscVector(grid,ONEDOF,global,GLOBAL)  
  call createPetscVector(grid,ONEDOF,natural,NATURAL)  

  ! write out coorindates
  call GetCoordinates(grid,global,X_COORDINATE)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Y_COORDINATE)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Z_COORDINATE)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  select case(option%imode)
    case (TWOPH_MODE,MPH_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)

      ! temperature
      call GetVarFromArray(solution,global,TEMPERATURE,0)
      call DMGlobalToNatural(grid,global,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

      ! pressure
      call GetVarFromArray(solution,global,PRESSURE,0)
      call DMGlobalToNatural(grid,global,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

      ! liquid saturation
      call GetVarFromArray(solution,global,LIQUID_SATURATION,0)
      call DMGlobalToNatural(grid,global,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

      ! gas saturation
      call GetVarFromArray(solution,global,GAS_SATURATION,0)
      call DMGlobalToNatural(grid,global,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)
    
      ! liquid energy
      call GetVarFromArray(solution,global,LIQUID_ENERGY,0)
      call DMGlobalToNatural(grid,global,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)
    
      ! gas energy
      call GetVarFromArray(solution,global,GAS_ENERGY,0)
      call DMGlobalToNatural(grid,global,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)
    
      ! liquid mole fractions
      do i=1,option%nspec
        call GetVarFromArray(solution,global,LIQUID_MOLE_FRACTION,i-1)
        call DMGlobalToNatural(grid,global,natural,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)
      enddo
    
      ! gas mole fractions
      do i=1,option%nspec
        call GetVarFromArray(solution,global,GAS_MOLE_FRACTION,i-1)
        call DMGlobalToNatural(grid,global,natural,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)
      enddo
    
      ! Volume Fraction
      if (option%rk > 0.d0) then
        call GetVarFromArray(solution,global,VOLUME_FRACTION,0)
        call DMGlobalToNatural(grid,global,natural,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)
      endif
    
      ! phase
      call GetVarFromArray(solution,global,PHASE,0)
      call DMGlobalToNatural(grid,global,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_INTEGER)
  
    case default
  
      ! temperature
      call DMGlobalToNatural(grid,option%temp,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

      ! pressure
      call DMGlobalToNatural(grid,option%pressure,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

      ! saturation
      call DMGlobalToNatural(grid,option%sat,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

      ! concentration
      call DMGlobalToNatural(grid,option%conc,natural,ONEDOF)
      call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

      ! volume fraction
      if (option%rk > 0.d0) then
        call GetVarFromArray(solution,global,VOLUME_FRACTION,0)
        call DMGlobalToNatural(grid,global,natural,ONEDOF)
        call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)
      endif
    
  end select
  
  call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  close(IUNIT3)
  
  if (option%print_tecplot_velocities) then
    call OutputVelocitiesTecplot(solution,kplot)
  endif
  
  if (option%print_tecplot_flux_velocities) then
    if (grid%structured_grid%nx > 1) then
      call OutputFluxVelocitiesTecplot(solution,kplot,LIQUID_PHASE,X_DIRECTION)
      call OutputFluxVelocitiesTecplot(solution,kplot,GAS_PHASE,X_DIRECTION)
    endif
    if (grid%structured_grid%ny > 1) then
      call OutputFluxVelocitiesTecplot(solution,kplot,LIQUID_PHASE,Y_DIRECTION)
      call OutputFluxVelocitiesTecplot(solution,kplot,GAS_PHASE,Y_DIRECTION)
    endif
    if (grid%structured_grid%nz > 1) then
      call OutputFluxVelocitiesTecplot(solution,kplot,LIQUID_PHASE,Z_DIRECTION)
      call OutputFluxVelocitiesTecplot(solution,kplot,GAS_PHASE,Z_DIRECTION)
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
subroutine OutputVelocitiesTecplot(solution,kplot)
 
  use Solution_module
  use Grid_module
  use Option_module
  
  implicit none

#include "definitions.h"

  type(solution_type) :: solution
  integer :: kplot
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global
  Vec :: natural

  real*8, pointer :: vec_ptr(:)
  
  grid => solution%grid
  option => solution%option
  
  ! open file
  if (kplot < 10) then
    write(filename,'("pflow_vel00",i1,".tec")') kplot  
  else if (kplot < 100) then
    write(filename,'("pflow_vel0",i2,".tec")') kplot  
  else if (kplot < 1000) then
    write(filename,'("pflow_vel",i3,".tec")') kplot  
  else if (kplot < 10000) then
    write(filename,'("pflow_vel",i4,".tec")') kplot  
  endif
  
  if (option%myrank == 0) then
    print *, '--> write tecplot velocity output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 option%t/option%tconv, option%tunit
    ! write variables
    string = 'VARIABLES=' // &
             '"X-Coordinates",' // &
             '"Y-Coordinates",' // &
             '"Z-Coordinates",' // &
             '"Liquid X-Velocity",' // &
             '"Liquid Y-Velocity",' // &
             '"Liquid Z-Velocity",' // &
             '"Gas X-Velocity",' // &
             '"Gas Y-Velocity",' // &
             '"Gas Z-Velocity"'
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                 &'', K='',i4,'','')') &
                 option%t/option%tconv,grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call createPetscVector(grid,ONEDOF,global,GLOBAL)  
  call createPetscVector(grid,ONEDOF,natural,NATURAL)    

  ! write out coorindates
  call GetCoordinates(grid,global,X_COORDINATE)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Y_COORDINATE)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Z_COORDINATE)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(solution,global,LIQUID_PHASE,X_DIRECTION)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(solution,global,LIQUID_PHASE,Y_DIRECTION)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(solution,global,LIQUID_PHASE,Z_DIRECTION)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(solution,global,GAS_PHASE,X_DIRECTION)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(solution,global,GAS_PHASE,Y_DIRECTION)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(solution,global,GAS_PHASE,Z_DIRECTION)
  call DMGlobalToNatural(grid,global,natural,ONEDOF)
  call WriteTecplotDataSetFromVec(IUNIT3,solution,natural,TECPLOT_REAL)

  call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

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
subroutine OutputFluxVelocitiesTecplot(solution,kplot,iphase,direction)
!geh - specifically, the flow velocities at the interfaces between cells
 
  use Solution_module
  use Grid_module
  use Option_module
  
  implicit none

#include "definitions.h"

  type(solution_type) :: solution
  integer :: kplot
  integer :: iphase
  integer :: direction
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  
  integer :: local_size, global_size
  integer :: nx_local, ny_local, nz_local
  integer :: nx_global, ny_global, nz_global
  integer :: i, j, k
  integer :: local_id
  integer :: adjusted_size
  integer :: count
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: array(:)
  PetscInt, allocatable :: indices(:)
  
  nullify(array)
  
  grid => solution%grid
  option => solution%option
  
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
  
  if (kplot < 10) then
    write(string,'("00",i1,".tec")') kplot  
  else if (kplot < 100) then
    write(string,'("0",i2,".tec")') kplot  
  else if (kplot < 1000) then
    write(string,'(i3,".tec")') kplot  
  else if (kplot < 10000) then
    write(string,'(i4,".tec")') kplot  
  endif
  
  filename = trim(filename) // trim(string)
  
  if (option%myrank == 0) then
    print *, '--> write tecplot velocity flux output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 option%t/option%tconv, option%tunit
    ! write variables
    string = 'VARIABLES=' // &
             '"X-Coordinates",' // &
             '"Y-Coordinates",' // &
             '"Z-Coordinates",'
    select case(iphase)
      case(LIQUID_PHASE)
        string = trim(string) // '"Liquid'
      case(GAS_PHASE)
        string = trim(string) // '"Gas'
    end select
  
    select case(direction)
      case(X_DIRECTION)
        string = trim(string) // ' X-Velocity"'
      case(Y_DIRECTION)
        string = trim(string) // ' Y-Velocity"'
      case(Z_DIRECTION)
        string = trim(string) // ' Z-Velocity"'
    end select 
    
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    select case(direction)
      case(X_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%t/option%tconv,grid%structured_grid%nx-1,grid%structured_grid%ny,grid%structured_grid%nz 
      case(Y_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%t/option%tconv,grid%structured_grid%nx,grid%structured_grid%ny-1,grid%structured_grid%nz 
      case(Z_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     option%t/option%tconv,grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz-1
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
  call WriteTecplotDataSet(IUNIT3,solution,array,TECPLOT_REAL,adjusted_size)
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
  call WriteTecplotDataSet(IUNIT3,solution,array,TECPLOT_REAL,adjusted_size)
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
  call WriteTecplotDataSet(IUNIT3,solution,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)
  
  ! write out data set
  call VecGetArrayF90(option%vl,vec_ptr,ierr)
  count = 0
  allocate(array(local_size))
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
  call VecRestoreArrayF90(option%vl,vec_ptr,ierr)
!GEH - Structured Grid Dependence - End
  
  array(1:local_size) = array(1:local_size)*option%tconv ! convert time units
  
  adjusted_size = local_size
  call ConvertArrayToNatural(indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,solution,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)
  
  deallocate(indices)

  close(IUNIT3)
  
end subroutine OutputFluxVelocitiesTecplot

! ************************************************************************** !
!
! WriteTecplotDataSetFromVec: Writes data from a Petsc Vec within a block
!                             of a Tecplot file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteTecplotDataSetFromVec(fid,solution,vec,datatype)

  use Solution_module
  
  implicit none

  integer :: fid
  type(solution_type) :: solution
  Vec :: vec
  integer :: datatype
  
  PetscScalar, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,solution,vec_ptr,datatype,0) ! 0 implies grid%nlmax
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
subroutine WriteTecplotDataSet(fid,solution,array,datatype,size_flag)

  use Solution_module

  implicit none
  
  integer :: fid
  type(solution_type) :: solution
  PetscReal :: array(:)
  integer, save :: max_local_size_saved = -1
  integer :: datatype
  integer :: size_flag ! if size_flag /= 0, use size_flag as the local size
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  integer :: i, iproc, recv_size
  integer :: max_local_size, local_size
  integer :: istart, iend, num_in_array
  integer :: status(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)
  
  grid => solution%grid
  option => solution%option
  
  if (size_flag /= 0) then
    call MPI_Allreduce(size_flag,max_local_size,1,MPI_INTEGER,MPI_MAX, &
                       PETSC_COMM_WORLD,ierr)
    local_size = size_flag
  else 
  ! if first time, determine the maximum size of any local array across 
  ! all procs
    if (max_local_size_saved < 0) then
      call MPI_Allreduce(grid%nlmax,max_local_size,1,MPI_INTEGER,MPI_MAX, &
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
        write(IUNIT3,'(10(i3,x))') integer_data(istart:iend)
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
        write(IUNIT3,'(10(es11.4,x))') real_data(istart:iend)
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
          write(IUNIT3,'(10(i3,x))') integer_data(istart:iend)
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
          write(IUNIT3,'(10(es11.4,x))') real_data(istart:iend)
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
        write(IUNIT3,'(10(i3,x))') integer_data(1:num_in_array)
    else
      if (num_in_array > 0) &
        write(IUNIT3,'(10(es11.4,x))') real_data(1:num_in_array)
    endif
  else
    if (datatype == TECPLOT_INTEGER) then
      call MPI_Send(integer_data,local_size,MPI_INTEGER,0,local_size, &
                    PETSC_COMM_WORLD,ierr)
    else
      call MPI_Send(real_data,local_size,MPI_DOUBLE_PRECISION,0,local_size, &
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
! OutputHDF5: Print to HDF5 file in Tecplot compatible format
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OutputHDF5(solution)

  use Solution_module
  use Option_module
  use Grid_module
  
#ifndef USE_HDF5
  implicit none
  
  type(solution_type) :: solution

  if (solution%option%myrank == 0) then
    print *
    print *, 'PFLOTRAN must be compiled with -DUSE_HDF5 to ', &
             'write to an HDF5 format.'
    print *
  endif
  stop

#else

  use hdf5
  
  implicit none

  type(solution_type) :: solution

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: solution_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: rank
  integer(HSIZE_T) :: dims(3)
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  Vec :: global
  Vec :: natural
  PetscScalar, pointer :: v_ptr
  
  character(len=MAXNAMELENGTH) :: filename = "pflow.h5"
  character(len=MAXSTRINGLENGTH) :: string
  logical, save :: first = .true.
  real*8, pointer :: array(:)
  integer :: i
  
  grid => solution%grid
  option => solution%option

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
    string = "X-Coordinates"
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
  
    string = "Y-Coordinates"
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
  
    string = "Z-Coordinates"
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
  write(string,'('' Time('',i4,''):'',es12.4,x,a1)') &
        option%flowsteps,option%t/option%tconv,option%tunit
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  
  ! write out data sets 
  call createPetscVector(grid,ONEDOF,global,GLOBAL)   

  select case(option%imode)
  
    case(TWOPH_MODE,MPH_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)
 
      ! temperature
      call GetVarFromArray(solution,global,TEMPERATURE,0)
      string = "Temperature"
      call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE)

      ! pressure
      call GetVarFromArray(solution,global,PRESSURE,0)
      string = "Pressure"
      call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE)

      ! liquid saturation
      call GetVarFromArray(solution,global,LIQUID_SATURATION,0)
      string = "Liquid Saturation"
      call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE)  

      ! gas saturation
      call GetVarFromArray(solution,global,GAS_SATURATION,0)
      string = "Gas Saturation"
      call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE) 
    
      ! liquid energy
      call GetVarFromArray(solution,global,LIQUID_ENERGY,0)
      string = "Liquid Energy"
      call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE) 
    
      ! gas energy
      call GetVarFromArray(solution,global,GAS_ENERGY,0)
      string = "Gas Energy"
      call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE) 
    
      ! liquid mole fractions
      do i=1,option%nspec
        call GetVarFromArray(solution,global,LIQUID_MOLE_FRACTION,i-1)
        write(string,'(''Liquid Mole Fraction('',i4,'')'')') i
        call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE)
      enddo
    
      ! gas mole fractions
      do i=1,option%nspec
        call GetVarFromArray(solution,global,GAS_MOLE_FRACTION,i-1)
        write(string,'(''Gas Mole Fraction('',i4,'')'')') i
        call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE)
      enddo
    
      ! Volume Fraction
      if (option%rk > 0.d0) then
        call GetVarFromArray(solution,global,VOLUME_FRACTION,0)
        string = "Volume Fraction"
        call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_DOUBLE)
      endif
    
      ! phase
      call GetVarFromArray(solution,global,PHASE,0)
      string = "Phase"
      call WriteHDF5DataSetFromVec(string,solution,global,grp_id,H5T_NATIVE_INTEGER) 
  
    case default
      ! temperature
      string = "Temperature"
      call WriteHDF5DataSetFromVec(string,solution,option%temp,grp_id, &
                                   H5T_NATIVE_DOUBLE)

      ! pressure
      string = "Pressure"
      call WriteHDF5DataSetFromVec(string,solution,option%pressure,grp_id, &
                                   H5T_NATIVE_DOUBLE)

      ! saturation
      string = "Saturation"
      call WriteHDF5DataSetFromVec(string,solution,option%sat,grp_id,H5T_NATIVE_DOUBLE)

      ! concentration
      string = "Concentration"
      call WriteHDF5DataSetFromVec(string,solution,option%conc,grp_id, &
                                   H5T_NATIVE_DOUBLE)

  end select
  
  if (option%print_hdf5_velocities) then

    ! velocities
    call GetCellCenteredVelocities(solution,global,LIQUID_PHASE,X_DIRECTION)
    string = "Liquid X-Velocity"
    call WriteHDF5DataSetFromVec(string,solution,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)
    call GetCellCenteredVelocities(solution,global,LIQUID_PHASE,Y_DIRECTION)
    string = "Liquid Y-Velocity"
    call WriteHDF5DataSetFromVec(string,solution,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  
    call GetCellCenteredVelocities(solution,global,LIQUID_PHASE,Z_DIRECTION)
    string = "Liquid Z-Velocity"
    call WriteHDF5DataSetFromVec(string,solution,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)

    call GetCellCenteredVelocities(solution,global,GAS_PHASE,X_DIRECTION)
    string = "Gas X-Velocity"
    call WriteHDF5DataSetFromVec(string,solution,global,grp_id, &
                               H5T_NATIVE_DOUBLE)
  
    call GetCellCenteredVelocities(solution,global,GAS_PHASE,Y_DIRECTION)
    string = "Gas Y-Velocity"
    call WriteHDF5DataSetFromVec(string,solution,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  
    call GetCellCenteredVelocities(solution,global,GAS_PHASE,Z_DIRECTION)
    string = "Gas Z-Velocity"
    call WriteHDF5DataSetFromVec(string,solution,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  endif

  if (option%print_hdf5_flux_velocities) then
  
    ! internal flux velocities
    if (grid%structured_grid%nx > 1) then
      string = "Liquid X-Flux Velocities"
      call WriteHDF5FluxVelocities(string,solution,LIQUID_PHASE,X_DIRECTION,grp_id)
      string = "Gas X-Flux Velocities"
      call WriteHDF5FluxVelocities(string,solution,GAS_PHASE,X_DIRECTION,grp_id)
    endif
    
    if (grid%structured_grid%ny > 1) then
      string = "Liquid Y-Flux Velocities"
      call WriteHDF5FluxVelocities(string,solution,LIQUID_PHASE,Y_DIRECTION,grp_id)
      string = "Gas Y-Flux Velocities"
      call WriteHDF5FluxVelocities(string,solution,GAS_PHASE,Y_DIRECTION,grp_id)
    endif
    
    if (grid%structured_grid%nz > 1) then
      string = "Liquid Z-Flux Velocities"
      call WriteHDF5FluxVelocities(string,solution,LIQUID_PHASE,Z_DIRECTION,grp_id)
      string = "Gas Z-Flux Velocities"
      call WriteHDF5FluxVelocities(string,solution,GAS_PHASE,Z_DIRECTION,grp_id)
    endif
    
  endif 
  
  ! call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
  first = .false.

end subroutine OutputHDF5

! ************************************************************************** !
!
! WriteHDF5FluxVelocities: Print flux velocities to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteHDF5FluxVelocities(name,solution,iphase,direction,file_id)

  use Solution_module
  use Grid_module
  use Option_module
  use hdf5

  implicit none
  
  character(len=32) :: name
  type(solution_type) :: solution
  integer :: iphase
  integer :: direction
  integer(HID_T) :: file_id

  integer :: i, j, k
  integer :: count
  integer :: local_id
  integer :: nx_local, ny_local, nz_local
  integer :: nx_global, ny_global, nz_global
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  real*8, allocatable :: array(:)
  PetscReal, pointer :: vec_ptr(:)

  logical, save :: first = .true.
  logical, save :: trick_flux_vel_x = .false.
  logical, save :: trick_flux_vel_y = .false.
  logical, save :: trick_flux_vel_z = .false.
  
  grid => solution%grid
  option => solution%option
  
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
    call MPI_Allreduce(nx_local,i,1,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (i == 0) trick_flux_vel_x = .true.
    if (grid%structured_grid%ngye-grid%structured_grid%nye == 0) then
      ny_local = grid%structured_grid%nly-1
    endif
    call MPI_Allreduce(ny_local,j,1,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (j == 0) trick_flux_vel_y = .true.
    if (grid%structured_grid%ngze-grid%structured_grid%nze == 0) then
      nz_local = grid%structured_grid%nlz-1
    endif
    call MPI_Allreduce(nz_local,k,1,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
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

  call VecGetArrayF90(option%vl,vec_ptr,ierr)
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
  call VecRestoreArrayF90(option%vl,vec_ptr,ierr)
  
  array(1:nx_local*ny_local*nz_local) = &  ! convert time units
    array(1:nx_local*ny_local*nz_local) * option%tconv

  call WriteHDF5DataSet(name,array,file_id,H5T_NATIVE_DOUBLE, &
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
  integer :: length
  real*8 :: array(:)
  integer(HID_T) :: file_id
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer :: rank
  
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

! ************************************************************************** !
!
! WriteHDF5DataSetFromVec: Writes data from a PetscVec to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteHDF5DataSetFromVec(name,solution,vec,file_id,data_type)

  use hdf5
  use Solution_module
  use Grid_module
  use Option_module
  
  implicit none

  character(len=32) :: name
  type(solution_type) :: solution
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscScalar, pointer :: vec_ptr(:)
  
  grid => solution%grid
  option => solution%option
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
!GEH - Structured Grid Dependence - Begin
  call WriteHDF5DataSet(name,vec_ptr,file_id,data_type, &
                        grid%structured_grid%nx,grid%structured_grid%ny,grid%structured_grid%nz, &
                        grid%structured_grid%nlx,grid%structured_grid%nly,grid%structured_grid%nlz, &
                        grid%structured_grid%nxs,grid%structured_grid%nys,grid%structured_grid%nzs)
!GEH - Structured Grid Dependence - End
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteHDF5DataSetFromVec

! ************************************************************************** !
!
! WriteHDF5DataSet: Writes data from an array into HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine WriteHDF5DataSet(name,array,file_id,data_type, &
                            nx_global,ny_global,nz_global, &
                            nx_local,ny_local,nz_local, &
                            istart_local,jstart_local,kstart_local)

  use hdf5
  
  implicit none
  
  character(len=32) :: name
  real*8 :: array(:)
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer :: nx_local, ny_local, nz_local
  integer :: nx_global, ny_global, nz_global
  integer :: istart_local, jstart_local, kstart_local
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer :: rank
  
  integer, pointer :: int_array(:)
  real*8, pointer :: double_array(:)
  integer :: i, j, k, count, id
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  integer :: ny_local_X_nz_local
  integer :: num_to_write

  ny_local_X_nz_local = ny_local*nz_local
  num_to_write = nx_local*ny_local_X_nz_local
  
  ! memory space which is a 1D vector  
  rank = 1
  dims = 0
  dims(1) = num_to_write
  if (num_to_write == 0) dims(1) = 1
  call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 3D block
  rank = 3
#define INVERT
#ifndef INVERT
    dims(1) = nx_global
    dims(2) = ny_global
    dims(3) = nz_global
#else
! have to trick hdf5 for now with inverted ordering
    dims(3) = nx_global
    dims(2) = ny_global
    dims(1) = nz_global
#endif
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims)


  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  call h5dcreate_f(file_id,name,data_type,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)
  
  ! create the hyperslab
#ifndef INVERT
  start(1) = istart_local
  start(2) = jstart_local
  start(3) = kstart_local
  length(1) =  nx_local
  length(2) =  ny_local
  length(3) =  nz_local
#else
  start(3) = istart_local
  start(2) = jstart_local
  start(1) = kstart_local
  length(3) =  nx_local
  length(2) =  ny_local
  length(1) =  nz_local
#endif
  if (num_to_write == 0) length(1) = 1
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err) 
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif
  if (num_to_write > 0) then
    if (data_type == H5T_NATIVE_INTEGER) then
      allocate(int_array(nx_local*ny_local*nz_local))
#ifdef INVERT
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            int_array(id) = int(array(count))
          enddo
        enddo
      enddo
#else
      do i=1,grid%nlmax
        int_array(i) = int(array(i))
      enddo
#endif
      call h5dwrite_f(data_set_id,data_type,int_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)
      deallocate(int_array)
    else
#ifdef INVERT
      allocate(double_array(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            double_array(id) = array(count)
          enddo
        enddo
      enddo
      call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)  
      deallocate(double_array)
#else
      call h5dwrite_f(data_set_id,data_type,array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)  
#endif
    endif
    call h5pclose_f(prop_id,hdf5_err)
  endif
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)

end subroutine WriteHDF5DataSet
!GEH - Structured Grid Dependence - End
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
  integer :: direction
  
  integer :: i
  PetscScalar, pointer :: vec_ptr(:)
  
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
  
  integer :: local_size, global_size
  integer :: indices(:)
  real*8, pointer :: array(:)
  
  Vec :: natural
  integer, allocatable :: indices_zero_based(:)
  real*8, pointer :: vec_ptr(:)
  
  call VecCreate(PETSC_COMM_WORLD,natural,ierr)
  call VecSetSizes(natural,PETSC_DECIDE,global_size,ierr)
  call VecSetType(natural,VECMPI,ierr)

  allocate(indices_zero_based(local_size))
  indices_zero_based(1:local_size) = indices(1:local_size)-1

  call VecSetValues(natural,local_size,indices_zero_based, &
                    array,INSERT_VALUES,ierr)

  call VecAssemblyBegin(natural,ierr)
  call VecAssemblyEnd(natural,ierr)

  call VecGetLocalSize(natural,local_size,ierr)
  deallocate(array)
  allocate(array(local_size))
  
  call VecGetArrayF90(natural,vec_ptr,ierr)
  array(1:local_size) = vec_ptr(1:local_size)
  call VecRestoreArrayF90(natural,vec_ptr,ierr)

  call VecDestroy(natural,ierr)
  
end subroutine ConvertArrayToNatural

! ************************************************************************** !
!
! GetVarFromArray: Extracts variables indexed by ivar from a multivar array
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GetVarFromArray(solution,vec,ivar,isubvar)

  use Solution_module
  use Grid_module
  use Option_module

  implicit none
  
  type(solution_type) :: solution
  Vec :: vec
  integer :: ivar
  integer :: isubvar

  integer :: i
  integer :: offset, saturation_offset
  integer :: size_var_use
  integer :: size_var_node
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscScalar, pointer :: var_ptr(:)
  PetscScalar, pointer :: vec_ptr(:)

  option => solution%option
  grid => solution%grid

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
        
      call VecGetArrayF90(option%var,var_ptr,ierr)
      do i=1,grid%nlmax
        vec_ptr(i) = var_ptr((i-1)*size_var_node+offset)
      enddo
      call VecRestoreArrayF90(option%var,var_ptr,ierr)

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
        
      call VecGetArrayF90(option%var,var_ptr,ierr)
      do i=1,grid%nlmax
        if (var_ptr((i-1)*size_var_node+saturation_offset) > 1.d-30) then
          vec_ptr(i) = var_ptr((i-1)*size_var_node+offset)
        else
          vec_ptr(i) = 0.d0
        endif
      enddo
      call VecRestoreArrayF90(option%var,var_ptr,ierr)

    case(VOLUME_FRACTION)
    
      ! need to set minimum to 0.
      call VecGetArrayF90(option%phis,var_ptr,ierr)
      vec_ptr(1:grid%nlmax) = var_ptr(1:grid%nlmax)
      call VecRestoreArrayF90(option%phis,var_ptr,ierr)
     
    case(PHASE)
    
      call VecGetArrayF90(option%iphas,var_ptr,ierr)
      vec_ptr(1:grid%nlmax) = var_ptr(1:grid%nlmax)
      call VecRestoreArrayF90(option%iphas,var_ptr,ierr)
     
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
subroutine GetCellCenteredVelocities(solution,vec,iphase,direction)

  use Solution_module
  use Grid_module
  use Option_module
  use Connection_module

  implicit none
  
  type(solution_type) :: solution
  Vec :: vec
  integer :: direction
  integer :: iphase
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  integer :: i, j, k, local_id, iconn
  Vec :: local_vec
  
  PetscScalar, pointer :: vec_ptr(:)
  PetscScalar, pointer :: vl_ptr(:)
  PetscScalar, pointer :: loc_vec_ptr(:)
  PetscInt, allocatable :: num_additions(:)
  
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_object
  
  grid => solution%grid
  option => solution%option
  
  allocate(num_additions(grid%nlmax))
  num_additions(1:grid%nlmax) = 0
  
  call VecSet(vec,0.d0,ierr)
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call VecGetArrayF90(option%vl,vl_ptr,ierr)
  do i=1,grid%nlmax
  ! definitely set up for a structured grid
  ! I believe that vl_ptr contains the downwind vel for all local nodes
    vec_ptr(i) = vl_ptr(iphase+(direction-1)*option%nphase+3*option%nphase*(i-1))
  enddo
  call VecRestoreArrayF90(option%vl,vl_ptr,ierr)
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
    
  call createPetscVector(grid,ONEDOF,local_vec,LOCAL)  
  call DMGlobalToLocal(grid,vec,local_vec,ONEDOF)

  call VecSet(vec,0.d0,ierr)

  call VecGetArrayF90(vec,vec_ptr,ierr)

  connection_list => grid%boundary_connection_list
  cur_connection_object => connection_list%first
  do 
    if (.not.associated(cur_connection_object)) exit
    do iconn = 1, cur_connection_object%num_connections
      if (cur_connection_object%dist(direction,iconn) < 0.99d0) cycle
      local_id = cur_connection_object%id_dn(iconn)
      if (iphase == LIQUID_PHASE) then
        vec_ptr(local_id) = vec_ptr(local_id) + option%vvlbc(iconn)
      else
        vec_ptr(local_id) = vec_ptr(local_id) + option%vvgbc(iconn)
      endif
      num_additions(local_id) = num_additions(local_id) + 1
    enddo
    cur_connection_object => cur_connection_object%next
  enddo

  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  call VecGetArrayF90(vec,vec_ptr,ierr)
  call VecGetArrayF90(local_vec,loc_vec_ptr,ierr)
  ! add in upwind portion from local vector
  
!GEH - Structured Grid Dependence - Begin
  select case(direction)
    case(X_DIRECTION)
      do local_id=1,grid%nlmax
        i= mod(local_id-1,grid%structured_grid%nlx) + 1
        if (i > 1 .or. grid%structured_grid%nxs-grid%structured_grid%ngxs > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id)-1)
          num_additions(local_id) = num_additions(local_id) + 1
        endif
        if (i < grid%structured_grid%nlx .or. grid%structured_grid%ngxe-grid%structured_grid%nxe > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id))
          num_additions(local_id) = num_additions(local_id) + 1
        endif
      enddo
    case(Y_DIRECTION)
      do local_id=1,grid%nlmax
        j= int(mod(local_id-1,grid%structured_grid%nlxy)/grid%structured_grid%nlx) + 1
        if (j > 1 .or. grid%structured_grid%nys-grid%structured_grid%ngys > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id)-grid%structured_grid%ngx)
          num_additions(local_id) = num_additions(local_id) + 1
        endif
        if (j < grid%structured_grid%nly .or. grid%structured_grid%ngye-grid%structured_grid%nye > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id))
          num_additions(local_id) = num_additions(local_id) + 1
        endif
      enddo
    case(Z_DIRECTION)
      do local_id=1,grid%nlmax
        k= int((local_id-1)/grid%structured_grid%nlxy) + 1
        if (k > 1 .or. grid%structured_grid%nzs-grid%structured_grid%ngzs > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id)-grid%structured_grid%ngxy)
          num_additions(local_id) = num_additions(local_id) + 1
        endif
        if (k < grid%structured_grid%nlz .or. grid%structured_grid%ngze-grid%structured_grid%nze > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id))
          num_additions(local_id) = num_additions(local_id) + 1
        endif
      enddo
  end select
!GEH - Structured Grid Dependence - End

  call VecRestoreArrayF90(local_vec,loc_vec_ptr,ierr)
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  call VecGetArrayF90(vec,vec_ptr,ierr)
  do i=1,grid%nlmax
    if (num_additions(i) > 0) &
      vec_ptr(i) = vec_ptr(i)/real(num_additions(i))*option%tconv
  enddo
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  call VecDestroy(local_vec,ierr)
  deallocate(num_additions)

end subroutine GetCellCenteredVelocities

end module Output_module

    
       
