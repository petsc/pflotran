module pflow_output2_module

  use pflow_gridtype_module

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
  
  public :: OutputTecplot, OutputHDF5
  
contains
  
subroutine OutputTecplot(grid,kplot)
 
  implicit none

#include "definitions.h"

  type(pflowGrid) :: grid
  integer :: kplot
  
  integer :: i
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string, string2
  Vec :: global
  Vec :: natural
  
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
  
  if (grid%myrank == 0) then
    print *, '--> write tecplot output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 grid%t/grid%tconv, grid%tunit
    ! write variables
    if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
        grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE .or. &
        grid%use_richards == PETSC_TRUE) then
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
      do i=1,grid%nspec
        write(string2,'(''"Liquid Mole Fraction('',i2,'')",'')') i
        string = trim(string) // trim(string2)
      enddo
      do i=1,grid%nspec
        write(string2,'(''"Gas Mole Fraction('',i2,'')",'')') i
        string = trim(string) // trim(string2)
      enddo
      if (grid%rk > 0.d0) then
        string = trim(string) // '"Volume Fraction",'
      endif
      string = trim(string) // '"Phase"'
    else
      string = 'VARIABLES=' // &
               '"X-Coordinates",' // &
               '"Y-Coordinates",' // &
               '"Z-Coordinates",' // &
               '"Temperature",' // &
               '"Pressure",' // &
               '"Saturation",' // &
               '"Concentration"'
      if (grid%rk > 0.d0) then
        string = trim(string) // ',"Volume Fraction"'
      endif
    endif
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                 &'', K='',i4,'','')') &
                 grid%t/grid%tconv,grid%nx,grid%ny,grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call DACreateGlobalVector(grid%da_1_dof,global,ierr)
  call DACreateNaturalVector(grid%da_1_dof,natural,ierr)

  ! write out coorindates
  call GetCoordinates(grid,global,X_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Y_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Z_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE .or. &
      grid%use_richards == PETSC_TRUE) then

    ! temperature
    call GetVarFromArray(grid,global,TEMPERATURE,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! pressure
    call GetVarFromArray(grid,global,PRESSURE,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! liquid saturation
    call GetVarFromArray(grid,global,LIQUID_SATURATION,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! gas saturation
    call GetVarFromArray(grid,global,GAS_SATURATION,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    
    ! liquid energy
    call GetVarFromArray(grid,global,LIQUID_ENERGY,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    
    ! gas energy
    call GetVarFromArray(grid,global,GAS_ENERGY,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    
    ! liquid mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,LIQUID_MOLE_FRACTION,i-1)
      call ConvertGlobalToNatural(grid,global,natural)
      call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    enddo
    
    ! gas mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,GAS_MOLE_FRACTION,i-1)
      call ConvertGlobalToNatural(grid,global,natural)
      call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    enddo
    
    ! Volume Fraction
    if (grid%rk > 0.d0) then
      call GetVarFromArray(grid,global,VOLUME_FRACTION,0)
      call ConvertGlobalToNatural(grid,global,natural)
      call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    endif
    
    ! phase
    call GetVarFromArray(grid,global,PHASE,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_INTEGER)
  
  else
  
    ! temperature
    call ConvertGlobalToNatural(grid,grid%temp,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! pressure
    call ConvertGlobalToNatural(grid,grid%pressure,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! saturation
    call ConvertGlobalToNatural(grid,grid%sat,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! concentration
    call ConvertGlobalToNatural(grid,grid%conc,natural)
    call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

    ! volume fraction
    if (grid%rk > 0.d0) then
      call GetVarFromArray(grid,global,VOLUME_FRACTION,0)
      call ConvertGlobalToNatural(grid,global,natural)
      call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
    endif
    
  endif
  
  call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  close(IUNIT3)
  
  if (grid%print_tecplot_velocities) then
    call OutputVelocitiesTecplot(grid,kplot)
  endif
  
  if (grid%print_tecplot_flux_velocities) then
    if (grid%nx > 1) then
      call OutputFluxVelocitiesTecplot(grid,kplot,LIQUID_PHASE,X_DIRECTION)
      call OutputFluxVelocitiesTecplot(grid,kplot,GAS_PHASE,X_DIRECTION)
    endif
    if (grid%ny > 1) then
      call OutputFluxVelocitiesTecplot(grid,kplot,LIQUID_PHASE,Y_DIRECTION)
      call OutputFluxVelocitiesTecplot(grid,kplot,GAS_PHASE,Y_DIRECTION)
    endif
    if (grid%nz > 1) then
      call OutputFluxVelocitiesTecplot(grid,kplot,LIQUID_PHASE,Z_DIRECTION)
      call OutputFluxVelocitiesTecplot(grid,kplot,GAS_PHASE,Z_DIRECTION)
    endif
  endif
      
end subroutine OutputTecplot

subroutine OutputVelocitiesTecplot(grid,kplot)
 
  implicit none

#include "definitions.h"

  type(pflowGrid) :: grid
  integer :: kplot
  
  character(len=MAXNAMELENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global
  Vec :: natural

  real*8, pointer :: vec_ptr(:)
  
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
  
  if (grid%myrank == 0) then
    print *, '--> write tecplot velocity output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 grid%t/grid%tconv, grid%tunit
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
                 grid%t/grid%tconv,grid%nx,grid%ny,grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  ! write out data sets  
  call DACreateGlobalVector(grid%da_1_dof,global,ierr)
  call DACreateNaturalVector(grid%da_1_dof,natural,ierr)

  ! write out coorindates
  call GetCoordinates(grid,global,X_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Y_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCoordinates(grid,global,Z_COORDINATE)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(grid,global,LIQUID_PHASE,X_DIRECTION)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(grid,global,LIQUID_PHASE,Y_DIRECTION)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(grid,global,LIQUID_PHASE,Z_DIRECTION)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(grid,global,GAS_PHASE,X_DIRECTION)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(grid,global,GAS_PHASE,Y_DIRECTION)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call GetCellCenteredVelocities(grid,global,GAS_PHASE,Z_DIRECTION)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteTecplotDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  close(IUNIT3)
  
end subroutine OutputVelocitiesTecplot

subroutine OutputFluxVelocitiesTecplot(grid,kplot,iphase,direction)
!geh - specifically, the flow velocities at the interfaces between cells
 
  implicit none

#include "definitions.h"

  type(pflowGrid) :: grid
  integer :: kplot
  integer :: iphase
  integer :: direction
  
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
  
  if (grid%myrank == 0) then
    print *, '--> write tecplot velocity flux output file: ', filename
    open(unit=IUNIT3,file=filename,action="write")
  
    ! write header
    ! write title
    write(IUNIT3,'(''TITLE = "'',1es12.4," [",a1,'']"'')') &
                 grid%t/grid%tconv, grid%tunit
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
                     grid%t/grid%tconv,grid%nx-1,grid%ny,grid%nz 
      case(Y_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     grid%t/grid%tconv,grid%nx,grid%ny-1,grid%nz 
      case(Z_DIRECTION)
        write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                     &'', K='',i4,'','')') &
                     grid%t/grid%tconv,grid%nx,grid%ny,grid%nz-1
    end select 
  
  
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

  endif
  
  ! write blocks
  
  ! face coordinates
  local_size = grid%nlmax
  global_size = grid%nmax
!GEH - Structured Grid Dependence - Begin
  nx_local = grid%nlx
  ny_local = grid%nly
  nz_local = grid%nlz
  nx_global = grid%nx
  ny_global = grid%ny
  nz_global = grid%nz
  select case(direction)
    case(X_DIRECTION)
      global_size = grid%nmax-grid%ny*grid%nz
      nx_global = grid%nx-1
      if (grid%ngxe-grid%nxe == 0) then
        local_size = grid%nlmax-grid%nlyz
        nx_local = grid%nlx-1
      endif
    case(Y_DIRECTION)
      global_size = grid%nmax-grid%nx*grid%nz
      ny_global = grid%ny-1
      if (grid%ngye-grid%nye == 0) then
        local_size = grid%nlmax-grid%nlxz
        ny_local = grid%nly-1
      endif
    case(Z_DIRECTION)
      global_size = grid%nmax-grid%nxy
      nz_global = grid%nz-1
      if (grid%ngze-grid%nze == 0) then
        local_size = grid%nlmax-grid%nlxy
        nz_local = grid%nlz-1
      endif
  end select  
  allocate(indices(local_size))

  ! fill indices array with natural ids in newly sized array
  count = 0
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        indices(count) = i+grid%nxs+(j-1+grid%nys)*nx_global+ &
                         (k-1+grid%nzs)*nx_global*ny_global
      enddo
    enddo
  enddo
  
  ! X-coordinates
  call VecGetArrayF90(grid%dx,vec_ptr,ierr)
  count = 0
  allocate(array(local_size))
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%nlx+(k-1)*grid%nlxy
        array(count) = grid%x(grid%nL2G(local_id))
        if (direction == X_DIRECTION) &
          array(count) = array(count) + 0.5d0*vec_ptr(local_id)
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(grid%dx,vec_ptr,ierr)
  ! warning: adjusted size will be changed in ConvertArrayToNatural
  ! thus, you cannot pass in local_size, since it is needed later
  adjusted_size = local_size
  call ConvertArrayToNatural(grid,indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,grid,array,TECPLOT_REAL,adjusted_size)
  ! since the array has potentially been resized, must reallocate
  deallocate(array)
  nullify(array)

  ! Y-coordinates
  count = 0
  allocate(array(local_size))
  call VecGetArrayF90(grid%dy,vec_ptr,ierr)
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%nlx+(k-1)*grid%nlxy
        array(count) = grid%y(grid%nL2G(local_id))
        if (direction == Y_DIRECTION) &
          array(count) = array(count) + 0.5d0*vec_ptr(local_id)
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(grid%dy,vec_ptr,ierr)
  adjusted_size = local_size
  call ConvertArrayToNatural(grid,indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,grid,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)

  ! Z-coordinates
  count = 0
  allocate(array(local_size))
  call VecGetArrayF90(grid%dz,vec_ptr,ierr)
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%nlx+(k-1)*grid%nlxy
        array(count) = grid%z(grid%nL2G(local_id))
        if (direction == Z_DIRECTION) &
          array(count) = array(count) + 0.5d0*vec_ptr(local_id)
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(grid%dz,vec_ptr,ierr)
  adjusted_size = local_size
  call ConvertArrayToNatural(grid,indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,grid,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)
  
  ! write out data set
  call VecGetArrayF90(grid%vl,vec_ptr,ierr)
  count = 0
  allocate(array(local_size))
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%nlx+(k-1)*grid%nlxy
        array(count) = vec_ptr(iphase+(direction-1)*grid%nphase+ &
                               3*grid%nphase*(local_id-1))
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(grid%vl,vec_ptr,ierr)
!GEH - Structured Grid Dependence - End
  
  array(1:local_size) = array(1:local_size)*grid%tconv ! convert time units
  
  adjusted_size = local_size
  call ConvertArrayToNatural(grid,indices,array,adjusted_size,global_size)
  call WriteTecplotDataSet(IUNIT3,grid,array,TECPLOT_REAL,adjusted_size)
  deallocate(array)
  nullify(array)
  
  deallocate(indices)

  close(IUNIT3)
  
end subroutine OutputFluxVelocitiesTecplot

subroutine WriteTecplotDataSetFromVec(fid,grid,vec,datatype)

  implicit none

  integer :: fid
  type(pflowGrid) :: grid
  Vec :: vec
  integer :: datatype
  
  PetscScalar, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteTecplotDataSet(fid,grid,vec_ptr,datatype,0) ! 0 implies grid%nlmax
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteTecplotDataSetFromVec

subroutine WriteTecplotDataSet(fid,grid,array,datatype,size_flag)

  implicit none
  
  integer :: fid
  type(pflowGrid) :: grid
  PetscReal :: array(:)
  integer, save :: max_local_size_saved = -1
  integer :: datatype
  integer :: size_flag ! if size_flag /= 0, use size_flag as the local size
  
  integer :: i, iproc, recv_size
  integer :: max_local_size, local_size
  integer :: istart, iend, num_in_array
  integer :: status(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)
  
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
      if (grid%myrank == 0) print *, 'max_local_size_saved: ', max_local_size
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
  if (grid%myrank == 0) then
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
    do iproc=1,grid%commsize-1
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

#ifndef USE_HDF5
subroutine OutputHDF5(grid)

  implicit none
  
  type(pflowGrid) :: grid

  if (grid%myrank == 0) then
    print *
    print *, 'PFLOTRAN must be compiled with -DUSE_HDF5 to ', &
             'write to an HDF5 format.'
    print *
  endif
  stop

end subroutine OutputHDF5

#else
subroutine OutputHDF5(grid)

  use hdf5
  
  implicit none

  type(pflowGrid) :: grid

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: rank
  integer(HSIZE_T) :: dims(3)
  
  Vec :: global
  Vec :: natural
  PetscScalar, pointer :: v_ptr
  
  character(len=MAXNAMELENGTH) :: filename = "pflow.h5"
  character(len=MAXSTRINGLENGTH) :: string
  logical, save :: first = .true.
  real*8, pointer :: array(:)
  integer :: i

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
    if (grid%myrank == 0) print *, '--> creating hdf5 output file: ', filename
  else
    if (grid%myrank == 0) print *, '--> appending to hdf5 output file: ', &
                                   filename
  endif
  
  if (first) then

    ! create a group for the coordinates data set
    string = "Coordinates"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)

!GEH - Structured Grid Dependence - Begin
    ! write out coordinates in x, y, and z directions
    string = "X-Coordinates"
    allocate(array(grid%nx))
    do i=1,grid%nx
      if (i == 1) then
        array(i) = 0.5d0*grid%dx0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%dx0(i-1)+grid%dx0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,grid,grid%nx,array,grp_id)
    deallocate(array)
  
    string = "Y-Coordinates"
    allocate(array(grid%ny))
    do i=1,grid%ny
      if (i == 1) then
        array(i) = 0.5d0*grid%dy0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%dy0(i-1)+grid%dy0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,grid,grid%ny,array,grp_id)
    deallocate(array)
  
    string = "Z-Coordinates"
    allocate(array(grid%nz))
    do i=1,grid%nz
      if (i == 1) then
        array(i) = 0.5d0*grid%dz0(1)
      else
        array(i) = array(i-1) + 0.5d0*(grid%dz0(i-1)+grid%dz0(i))
      endif
    enddo
    call WriteHDF5Coordinates(string,grid,grid%nz,array,grp_id)
    deallocate(array)
!GEH - Structured Grid Dependence - End

    call h5gclose_f(grp_id,hdf5_err)
    
  endif

  ! create a group for the data set
  write(string,'('' Time('',i4,''):'',es12.4,x,a1)') &
        grid%flowsteps,grid%t/grid%tconv,grid%tunit
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  
  ! write out data sets  
  call DACreateGlobalVector(grid%da_1_dof,global,ierr)
  !call DACreateNaturalVector(grid%da_1_dof,natural,ierr)

  if (grid%use_2ph == PETSC_TRUE .or. grid%use_mph == PETSC_TRUE .or. &
      grid%use_vadose == PETSC_TRUE .or. grid%use_flash == PETSC_TRUE .or. &
      grid%use_richards == PETSC_TRUE) then
  
    ! temperature
    call GetVarFromArray(grid,global,TEMPERATURE,0)
    !call DAGlobalToNaturalBegin(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)
    !call DAGlobalToNaturalEnd(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)
    string = "Temperature"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)

    ! pressure
    call GetVarFromArray(grid,global,PRESSURE,0)
    string = "Pressure"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)

    ! liquid saturation
    call GetVarFromArray(grid,global,LIQUID_SATURATION,0)
    string = "Liquid Saturation"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)  

    ! gas saturation
    call GetVarFromArray(grid,global,GAS_SATURATION,0)
    string = "Gas Saturation"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! liquid energy
    call GetVarFromArray(grid,global,LIQUID_ENERGY,0)
    string = "Liquid Energy"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! gas energy
    call GetVarFromArray(grid,global,GAS_ENERGY,0)
    string = "Gas Energy"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE) 
    
    ! liquid mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,LIQUID_MOLE_FRACTION,i-1)
      write(string,'(''Liquid Mole Fraction('',i4,'')'')') i
      call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)
    enddo
    
    ! gas mole fractions
    do i=1,grid%nspec
      call GetVarFromArray(grid,global,GAS_MOLE_FRACTION,i-1)
      write(string,'(''Gas Mole Fraction('',i4,'')'')') i
      call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)
    enddo
    
    ! Volume Fraction
    if (grid%rk > 0.d0) then
      call GetVarFromArray(grid,global,VOLUME_FRACTION,0)
      string = "Volume Fraction"
      call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_DOUBLE)
    endif
    
    ! phase
    call GetVarFromArray(grid,global,PHASE,0)
    string = "Phase"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id,H5T_NATIVE_INTEGER) 
  
  else
  
    ! temperature
    string = "Temperature"
    call WriteHDF5DataSetFromVec(string,grid,grid%temp,grp_id, &
                                 H5T_NATIVE_DOUBLE)

    ! pressure
    string = "Pressure"
    call WriteHDF5DataSetFromVec(string,grid,grid%pressure,grp_id, &
                                 H5T_NATIVE_DOUBLE)

    ! saturation
    string = "Saturation"
    call WriteHDF5DataSetFromVec(string,grid,grid%sat,grp_id,H5T_NATIVE_DOUBLE)

    ! concentration
    string = "Concentration"
    call WriteHDF5DataSetFromVec(string,grid,grid%conc,grp_id, &
                                 H5T_NATIVE_DOUBLE)

  endif

  if (grid%print_hdf5_velocities) then

    ! velocities
    call GetCellCenteredVelocities(grid,global,LIQUID_PHASE,X_DIRECTION)
    string = "Liquid X-Velocity"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)
    call GetCellCenteredVelocities(grid,global,LIQUID_PHASE,Y_DIRECTION)
    string = "Liquid Y-Velocity"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  
    call GetCellCenteredVelocities(grid,global,LIQUID_PHASE,Z_DIRECTION)
    string = "Liquid Z-Velocity"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)

    call GetCellCenteredVelocities(grid,global,GAS_PHASE,X_DIRECTION)
    string = "Gas X-Velocity"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id, &
                               H5T_NATIVE_DOUBLE)
  
    call GetCellCenteredVelocities(grid,global,GAS_PHASE,Y_DIRECTION)
    string = "Gas Y-Velocity"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  
    call GetCellCenteredVelocities(grid,global,GAS_PHASE,Z_DIRECTION)
    string = "Gas Z-Velocity"
    call WriteHDF5DataSetFromVec(string,grid,global,grp_id, &
                                 H5T_NATIVE_DOUBLE)
  endif

  if (grid%print_hdf5_flux_velocities) then
  
    ! internal flux velocities
    if (grid%nx > 1) then
      string = "Liquid X-Flux Velocities"
      call WriteHDF5FluxVelocities(string,grid,LIQUID_PHASE,X_DIRECTION,grp_id)
      string = "Gas X-Flux Velocities"
      call WriteHDF5FluxVelocities(string,grid,GAS_PHASE,X_DIRECTION,grp_id)
    endif
    
    if (grid%ny > 1) then
      string = "Liquid Y-Flux Velocities"
      call WriteHDF5FluxVelocities(string,grid,LIQUID_PHASE,Y_DIRECTION,grp_id)
      string = "Gas Y-Flux Velocities"
      call WriteHDF5FluxVelocities(string,grid,GAS_PHASE,Y_DIRECTION,grp_id)
    endif
    
    if (grid%nz > 1) then
      string = "Liquid Z-Flux Velocities"
      call WriteHDF5FluxVelocities(string,grid,LIQUID_PHASE,Z_DIRECTION,grp_id)
      string = "Gas Z-Flux Velocities"
      call WriteHDF5FluxVelocities(string,grid,GAS_PHASE,Z_DIRECTION,grp_id)
    endif
    
  endif 
  
  ! call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  call h5gclose_f(grp_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
  first = .false.

end subroutine OutputHDF5

subroutine WriteHDF5FluxVelocities(name,grid,iphase,direction,file_id)

  use hdf5

  implicit none
  
  character(len=32) :: name
  type(pflowGrid) :: grid
  integer :: iphase
  integer :: direction
  integer(HID_T) :: file_id

  integer :: i, j, k
  integer :: count
  integer :: local_id
  integer :: nx_local, ny_local, nz_local
  integer :: nx_global, ny_global, nz_global
  
  real*8, allocatable :: array(:)
  PetscReal, pointer :: vec_ptr(:)

  logical, save :: first = .true.
  logical, save :: trick_flux_vel_x = .false.
  logical, save :: trick_flux_vel_y = .false.
  logical, save :: trick_flux_vel_z = .false.
  
  ! in a few cases (i.e. for small test problems), some processors may
  ! have no velocities to print.  This results in zero-length arrays
  ! in collective H5Dwrite().  To avoid, we switch to independent
  ! H5Dwrite() and don't write from the zero-length procs. 
!GEH - Structured Grid Dependence - Begin
  if (first) then
    nx_local = grid%nlx
    ny_local = grid%nly
    nz_local = grid%nlz
    if (grid%ngxe-grid%nxe == 0) then
      nx_local = grid%nlx-1
    endif
    call MPI_Allreduce(nx_local,i,1,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (i == 0) trick_flux_vel_x = .true.
    if (grid%ngye-grid%nye == 0) then
      ny_local = grid%nly-1
    endif
    call MPI_Allreduce(ny_local,j,1,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (j == 0) trick_flux_vel_y = .true.
    if (grid%ngze-grid%nze == 0) then
      nz_local = grid%nlz-1
    endif
    call MPI_Allreduce(nz_local,k,1,MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
    if (k == 0) trick_flux_vel_z = .true.
  endif

  nx_local = grid%nlx
  ny_local = grid%nly
  nz_local = grid%nlz
  nx_global = grid%nx
  ny_global = grid%ny
  nz_global = grid%nz

  select case(direction)
    case(X_DIRECTION)
      nx_global = grid%nx-1
      if (grid%ngxe-grid%nxe == 0) then
        nx_local = grid%nlx-1
      endif
      if (trick_flux_vel_x) trick_hdf5 = .true.
    case(Y_DIRECTION)
      ny_global = grid%ny-1
      if (grid%ngye-grid%nye == 0) then
        ny_local = grid%nly-1
      endif
      if (trick_flux_vel_y) trick_hdf5 = .true.
    case(Z_DIRECTION)
      nz_global = grid%nz-1
      if (grid%ngze-grid%nze == 0) then
        nz_local = grid%nlz-1
      endif
      if (trick_flux_vel_z) trick_hdf5 = .true.
  end select  
  allocate(array(nx_local*ny_local*nz_local))

  call VecGetArrayF90(grid%vl,vec_ptr,ierr)
  count = 0
  do k=1,nz_local
    do j=1,ny_local
      do i=1,nx_local
        count = count + 1
        local_id = i+(j-1)*grid%nlx+(k-1)*grid%nlxy
        array(count) = vec_ptr(iphase+(direction-1)*grid%nphase+ &
                               3*grid%nphase*(local_id-1))
      enddo
    enddo
  enddo
  call VecRestoreArrayF90(grid%vl,vec_ptr,ierr)
  
  array(1:nx_local*ny_local*nz_local) = &  ! convert time units
    array(1:nx_local*ny_local*nz_local) * grid%tconv

  call WriteHDF5DataSet(name,array,file_id,H5T_NATIVE_DOUBLE, &
                        nx_global,ny_global,nz_global, &
                        nx_local,ny_local,nz_local, &
                        grid%nxs,grid%nys,grid%nzs)
!GEH - Structured Grid Dependence - End

  deallocate(array)
  trick_hdf5 = .false.
  first = .false.

end subroutine WriteHDF5FluxVelocities

subroutine WriteHDF5Coordinates(name,grid,length,array,file_id)

  use hdf5
  
  implicit none
  
  character(len=32) :: name
  type(pflowGrid) :: grid
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
  if (grid%myrank == 0) then
    call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,array,dims, &
                    hdf5_err,H5S_ALL_F,H5S_ALL_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

end subroutine WriteHDF5Coordinates

subroutine WriteHDF5DataSetFromVec(name,grid,vec,file_id,data_type)

  use hdf5
  
  implicit none

  character(len=32) :: name
  type(pflowGrid) :: grid
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  PetscScalar, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
!GEH - Structured Grid Dependence - Begin
  call WriteHDF5DataSet(name,vec_ptr,file_id,data_type, &
                        grid%nx,grid%ny,grid%nz, &
                        grid%nlx,grid%nly,grid%nlz, &
                        grid%nxs,grid%nys,grid%nzs)
!GEH - Structured Grid Dependence - End
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteHDF5DataSetFromVec

!GEH - Structured Grid Dependence - Begin
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
                            hdf5_err) ! must be independent and only from p0
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err) ! must be independent and only from p0
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

subroutine GetCoordinates(grid,vec,direction)

  implicit none
  
  type(pflowGrid) :: grid
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

subroutine ConvertGlobalToNatural(grid,global,natural)

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: global
  Vec :: natural
  
  call DAGlobalToNaturalBegin(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)
  call DAGlobalToNaturalEnd(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)

end subroutine 

subroutine ConvertArrayToNatural(grid,indices,array, &
                                 local_size,global_size)

  implicit none
  
  type(pflowGrid) :: grid
  
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

subroutine GetVarFromArray(grid,vec,ivar,isubvar)

  use pflow_gridtype_module

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: vec
  integer :: ivar
  integer :: isubvar

  integer :: i
  integer :: offset, saturation_offset
  integer :: size_var_use
  integer :: size_var_node
  PetscScalar, pointer :: var_ptr(:)
  PetscScalar, pointer :: vec_ptr(:)

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
          offset = 17+grid%nspec+isubvar
      end select
    
      size_var_use = 2 + 7*grid%nphase + 2* grid%nphase*grid%nspec
      size_var_node = (grid%ndof + 1) * size_var_use
        
      call VecGetArrayF90(grid%var,var_ptr,ierr)
      do i=1,grid%nlmax
        vec_ptr(i) = var_ptr((i-1)*size_var_node+offset)
      enddo
      call VecRestoreArrayF90(grid%var,var_ptr,ierr)

    case(LIQUID_ENERGY,GAS_ENERGY)

      select case (ivar)
        case(LIQUID_ENERGY)
          offset = 11
          saturation_offset = 3
        case(GAS_ENERGY)
          offset = 12
          saturation_offset = 4  
      end select

      size_var_use = 2 + 7*grid%nphase + 2* grid%nphase*grid%nspec
      size_var_node = (grid%ndof + 1) * size_var_use
        
      call VecGetArrayF90(grid%var,var_ptr,ierr)
      do i=1,grid%nlmax
        if (var_ptr((i-1)*size_var_node+saturation_offset) > 1.d-30) then
          vec_ptr(i) = var_ptr((i-1)*size_var_node+offset)
        else
          vec_ptr(i) = 0.d0
        endif
      enddo
      call VecRestoreArrayF90(grid%var,var_ptr,ierr)

    case(VOLUME_FRACTION)
    
      ! need to set minimum to 0.
      call VecGetArrayF90(grid%phis,var_ptr,ierr)
      vec_ptr(1:grid%nlmax) = var_ptr(1:grid%nlmax)
      call VecRestoreArrayF90(grid%phis,var_ptr,ierr)
     
    case(PHASE)
    
      call VecGetArrayF90(grid%iphas,var_ptr,ierr)
      vec_ptr(1:grid%nlmax) = var_ptr(1:grid%nlmax)
      call VecRestoreArrayF90(grid%iphas,var_ptr,ierr)
     
  end select
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

end subroutine GetVarFromArray

subroutine GetCellCenteredVelocities(grid,vec,iphase,direction)

  use pflow_gridtype_module

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: vec
  integer :: direction
  integer :: iphase
  
  integer :: i, j, k, local_id, iconn
  Vec :: local_vec
  
  PetscScalar, pointer :: vec_ptr(:)
  PetscScalar, pointer :: vl_ptr(:)
  PetscScalar, pointer :: loc_vec_ptr(:)
  PetscInt, allocatable :: num_additions(:)
  
  allocate(num_additions(grid%nlmax))
  num_additions(1:grid%nlmax) = 0
  
  call VecSet(vec,0.d0,ierr)
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call VecGetArrayF90(grid%vl,vl_ptr,ierr)
  do i=1,grid%nlmax
  ! I believe that vl_ptr contains the downwind vel for all local nodes
    vec_ptr(i) = vl_ptr(iphase+(direction-1)*grid%nphase+3*grid%nphase*(i-1))
  enddo
  call VecRestoreArrayF90(grid%vl,vl_ptr,ierr)
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
    
  call DACreateLocalVector(grid%da_1_dof,local_vec,ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof,vec,INSERT_VALUES,local_vec,ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof,vec,INSERT_VALUES,local_vec,ierr)

  call VecSet(vec,0.d0,ierr)

  call VecGetArrayF90(vec,vec_ptr,ierr)
  do iconn = 1, grid%nconnbc
    if (grid%ipermbc(iconn) == direction) then ! direction 1=x,2=y,3=z
      local_id = grid%mblkbc(iconn)  ! node of bc m = local id
      if (iphase == LIQUID_PHASE) then
        vec_ptr(local_id) = vec_ptr(local_id) + grid%vvlbc(iconn)
      else
        vec_ptr(local_id) = vec_ptr(local_id) + grid%vvgbc(iconn)
      endif
      num_additions(local_id) = num_additions(local_id) + 1
    endif
  enddo
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  call VecGetArrayF90(vec,vec_ptr,ierr)
  call VecGetArrayF90(local_vec,loc_vec_ptr,ierr)
  ! add in upwind portion from local vector
  
!GEH - Structured Grid Dependence - Begin
  select case(direction)
    case(X_DIRECTION)
      do local_id=1,grid%nlmax
        i= mod(local_id-1,grid%nlx) + 1
        if (i > 1 .or. grid%nxs-grid%ngxs > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id)-1)
          num_additions(local_id) = num_additions(local_id) + 1
        endif
        if (i < grid%nlx .or. grid%ngxe-grid%nxe > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id))
          num_additions(local_id) = num_additions(local_id) + 1
        endif
      enddo
    case(Y_DIRECTION)
      do local_id=1,grid%nlmax
        j= int(mod(local_id-1,grid%nlxy)/grid%nlx) + 1
        if (j > 1 .or. grid%nys-grid%ngys > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id)-grid%ngx)
          num_additions(local_id) = num_additions(local_id) + 1
        endif
        if (j < grid%nly .or. grid%ngye-grid%nye > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id))
          num_additions(local_id) = num_additions(local_id) + 1
        endif
      enddo
    case(Z_DIRECTION)
      do local_id=1,grid%nlmax
        k= int((local_id-1)/grid%nlxy) + 1
        if (k > 1 .or. grid%nzs-grid%ngzs > 0) then
          vec_ptr(local_id) = vec_ptr(local_id) + &
                              loc_vec_ptr(grid%nL2G(local_id)-grid%ngxy)
          num_additions(local_id) = num_additions(local_id) + 1
        endif
        if (k < grid%nlz .or. grid%ngze-grid%nze > 0) then
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
      vec_ptr(i) = vec_ptr(i)/real(num_additions(i))*grid%tconv
  enddo
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  call VecDestroy(local_vec,ierr)
  deallocate(num_additions)

end subroutine GetCellCenteredVelocities

end module pflow_output2_module

    
       
