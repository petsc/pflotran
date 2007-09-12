module tecplot_output_module

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

  PetscErrorCode :: ierr
  
  public :: OutputTecplot
  
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
    print *, '--> write output file: ', filename
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
               '"X-coordinates",' // &
               '"Y-coordinates",' // &
               '"Z-coordinates",' // &
               '"Temperature",' // &
               '"Pressure",' // &
               '"Liquid Saturation",' // &
               '"Gas Saturation",' // &
               '"Liquid Energy",' // &
               '"Gas Energy",'
      do i=1,grid%nspec
        write(string2,'(''"Liquid Mole Fraction('',i2,'')",'')') i
        string = trim(string) // trim(string2)
        write(string2,'(''"Gas Mole Fraction('',i2,'')",'')') i
        string = trim(string) // trim(string2)
      enddo
      if (grid%rk > 0.d0) then
        string = trim(string) // '"Volume Fraction",'
      endif
      string = trim(string) // '"Phase"'
    else
      string = '"X-coordinates",' // &
               '"Y-coordinates",' // &
               '"Z-coordinates",' // &
               '"Temperature",' // &
               '"Pressure",' // &
               '"Saturation",' // &
               '"Concentration",'
      if (grid%rk > 0.d0) then
        string = trim(string) // '"Volume Fraction",'
      endif
      string = trim(string) // '"Phase"'
    endif
    write(IUNIT3,'(a)') trim(string)
  
    ! write zone header
    write(string,'(''ZONE T= "'',1es12.4,''",'','' I='',i4,'', J='',i4, &
                 &'', K='',i4,'','')') &
                 grid%t/grid%tconv,grid%nx,grid%ny,grid%nz 
    string = trim(string) // ' DATAPACKING=BLOCK'
    write(IUNIT3,'(a)') trim(string)

    ! write blocks
  
    ! write out coorindates
    write(IUNIT3,'(10(es11.4,x))') (grid%x(i),i=1,grid%nmax)
    write(IUNIT3,'(10(es11.4,x))') (grid%y(i),i=1,grid%nmax)
    write(IUNIT3,'(10(es11.4,x))') (grid%z(i),i=1,grid%nmax)
    
  endif
    
  ! write out data sets  
  call DACreateGlobalVector(grid%da_1_dof,global,ierr)
  call DACreateNaturalVector(grid%da_1_dof,natural,ierr)

  ! temperature
  call GetVarFromArray(grid,global,TEMPERATURE,0)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  ! pressure
  call GetVarFromArray(grid,global,PRESSURE,0)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  ! liquid saturation
  call GetVarFromArray(grid,global,LIQUID_SATURATION,0)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)

  ! gas saturation
  call GetVarFromArray(grid,global,GAS_SATURATION,0)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
  
  ! liquid energy
  call GetVarFromArray(grid,global,LIQUID_ENERGY,0)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
  
  ! gas energy
  call GetVarFromArray(grid,global,GAS_ENERGY,0)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
  
  ! liquid mole fractions
  do i=1,grid%nspec
    call GetVarFromArray(grid,global,LIQUID_MOLE_FRACTION,i-1)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
  enddo
  
  ! gas mole fractions
  do i=1,grid%nspec
    call GetVarFromArray(grid,global,GAS_MOLE_FRACTION,i-1)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
  enddo
  
  ! Volume Fraction
  if (grid%rk > 0.d0) then
    call GetVarFromArray(grid,global,VOLUME_FRACTION,0)
    call ConvertGlobalToNatural(grid,global,natural)
    call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_REAL)
  endif
  
  ! phase
  call GetVarFromArray(grid,global,PHASE,0)
  call ConvertGlobalToNatural(grid,global,natural)
  call WriteDataSetFromVec(IUNIT3,grid,natural,TECPLOT_INTEGER)
  
  call VecDestroy(natural,ierr)
  call VecDestroy(global,ierr)

  close(IUNIT3)
      
end subroutine OutputTecplot

subroutine GetVarFromArray(grid,vec,ivar,isubvar)

  use pflow_gridtype_module

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: vec
  integer :: ivar
  integer :: isubvar

  integer :: i
  integer :: offset
  integer :: size_var_use
  integer :: size_var_node
  PetscScalar, pointer :: var_ptr(:)
  PetscScalar, pointer :: vec_ptr(:)

  call VecGetArrayF90(vec,vec_ptr,ierr)
      
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_ENERGY,GAS_ENERGY,LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION)
      select case(ivar)
        case(TEMPERATURE)
          offset = 1
        case(PRESSURE)
          offset = 2
        case(LIQUID_SATURATION)
          offset = 3
        case(GAS_SATURATION)
          offset = 4
        case(LIQUID_ENERGY)
          offset = 11
        case(GAS_ENERGY)
          offset = 12    
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

subroutine ConvertGlobalToNatural(grid,global,natural)

  implicit none
  
  type(pflowGrid) :: grid
  Vec :: global
  Vec :: natural
  
  call DAGlobalToNaturalBegin(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)
  call DAGlobalToNaturalEnd(grid%da_1_dof,global,INSERT_VALUES,natural,ierr)

end subroutine 
  
subroutine WriteDataSetFromVec(fid,grid,vec,datatype)

  implicit none

  integer :: fid
  type(pflowGrid) :: grid
  Vec :: vec
  integer :: datatype
  
  PetscScalar, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
  call WriteDataSet(fid,grid,vec_ptr,datatype)
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine WriteDataSetFromVec

subroutine WriteDataSet(fid,grid,vec_ptr,datatype)

  implicit none
  
  integer :: fid
  type(pflowGrid) :: grid
  PetscReal, pointer :: vec_ptr(:)
  integer, save :: max_local_size = -1
  integer :: datatype
  
  integer :: i, iproc, recv_size
  integer :: status(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:)
  PetscReal, allocatable :: real_data(:)
  
  ! if first time, determine the maximum size of any local array across all procs
  if (max_local_size < 0) then
    call MPI_Allreduce(grid%nlmax,max_local_size,1,MPI_INTEGER,MPI_MAX, &
                       PETSC_COMM_WORLD,ierr)
    if (grid%myrank == 0) print *, 'max_local_size: ', max_local_size
  endif
  
  ! transfer the data to an integer or real array
  if (datatype == TECPLOT_INTEGER) then
    allocate(integer_data(max_local_size))
    do i=1,grid%nlmax
      integer_data(i) = int(vec_ptr(i))
    enddo
  else
    allocate(real_data(max_local_size))
    do i=1,grid%nlmax
      real_data(i) = vec_ptr(i)
    enddo
  endif
  
  ! communicate data to processor 0, round robin style
  if (grid%myrank == 0) then
    if (datatype == TECPLOT_INTEGER) then
      write(IUNIT3,'(10(i3,x))') (integer_data(i),i=1,grid%nlmax)
    else
      write(IUNIT3,'(10(es11.4,x))') (real_data(i),i=1,grid%nlmax)
    endif
    do iproc=1,grid%commsize-1
      call MPI_Probe(iproc,MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
      recv_size = status(MPI_TAG)
      if (datatype == 0) then
        call MPI_Recv(integer_data,recv_size,MPI_INTEGER,iproc,MPI_ANY_TAG, &
                      PETSC_COMM_WORLD,status,ierr)
        write(IUNIT3,'(10(i3,x))') (integer_data(i),i=1,recv_size)
      else
        call MPI_Recv(real_data,recv_size,MPI_DOUBLE_PRECISION,iproc, &
                      MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
        write(IUNIT3,'(10(es11.4,x))') (real_data(i),i=1,recv_size)
      endif
    enddo
  else
    if (datatype == TECPLOT_INTEGER) then
      call MPI_Send(integer_data,grid%nlmax,MPI_INTEGER,0,grid%nlmax, &
                    PETSC_COMM_WORLD,ierr)
    else
      call MPI_Send(real_data,grid%nlmax,MPI_DOUBLE_PRECISION,0,grid%nlmax, &
                    PETSC_COMM_WORLD,ierr)
    endif
  endif
      
  if (datatype == TECPLOT_INTEGER) then
    deallocate(integer_data)
  else
    deallocate(real_data)
  endif

end subroutine WriteDataSet
  
end module tecplot_output_module

    
       
