module output

  use globals, only: wordlen

  implicit none
  save

#define DEBUG 0
#define USE_OLD

#define COMPS 1
#define MNRLS 2
#define NAPLS 3

#define X_START        1
#define X_END          2
#define Y_START        3
#define Y_END          4
#define Z_START        5
#define Z_END          6
#define N_SPATIAL_DOF  7

#define X_DIM 1
#define Y_DIM 2
#define Z_DIM 3
#define XY_DIM 4
#define YZ_DIM 5
#define XZ_DIM 6
#define XYZ_DIM 7




! by default, all variables are private
  private

! specify public types
!  public ::

! specify public variables
!  public :: 

! specify public functions and subroutines

!  public :: ouPrintResults, ouPrintGMSResults, ouPrintMineralResults, &
!            ouPrintVector, ouPrintDomains, ouPrintGMSGrid

  public :: ouPostProcessResults

  public :: ouPrintResults, ouInitOutput, ouFinalizeOutput

  public :: ouPrintPETScMatrix, ouPrintPETScVector, ouPrintGMSGRid, &
            ouPrintL_GConversionArrays

  public :: ouPrintBoundaryConnections, ouPrintConnections, &
            ouPrintArray, ouPrintVectorBinary_old, ouPrintGMSVector, &
            ouPrintVelocities, ouDeallocateOutput, &
            ouPrintDispersionCoefficients, &
            ouPrintConcentrationsBinary

! type declarations

! file globals - variables
! static
  logical :: petsc = .false.
  real*8, public :: ou_porosity

  integer :: timeunit       ! unit # for time output file
  integer :: tempb_unit  ! unit # for binary output file
  integer :: tempb_priunit  ! unit # for binary primary species output file
  integer :: tempb_mnrlunit ! unit # for binary mineral output file
  integer :: tempb_naplunit ! unit # for binary NAPL output file

  character(len=wordlen) :: time_filename
  character(len=wordlen) :: binary_filename
  character(len=wordlen) :: pri_filename
  character(len=wordlen) :: mnrl_filename
  character(len=wordlen) :: napl_filename

  integer, allocatable :: indices(:,:)
  real*8, allocatable :: output_buffer(:)

! dynamic

! __________________________________________________________________________ !
! __________________________________________________________________________ !

contains

! ************************************************************************** !
!
! ouInitOutput: Initialize output files
! author: Glenn Hammond
! date: 10/27/03
!
! ************************************************************************** !
subroutine ouInitOutput()

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscviewer.h"

  integer :: ierr, iproc
  integer :: local_indices(7), status(MPI_STATUS_SIZE)

  time_filename = 'time.tmp'
  binary_filename = 'data.tmp'

! Temporary files for binary output during model run (much faster)
  if (myrank == 0) then
    timeunit = 29
    open(unit=timeunit,file=time_filename)
    call PetscBinaryOpen(binary_filename,PETSC_BINARY_CREATE, &
                         tempb_unit,ierr)
  endif

  local_indices(X_START) = lxs+1
  local_indices(X_END) = lxe
  local_indices(Y_START) = lys+1
  local_indices(Y_END) = lye
  local_indices(Z_START) = lzs+1
  local_indices(Z_END) = lze
  local_indices(N_SPATIAL_DOF) = lnx*lny*lnz

  ! obtain indices of other processors
  if (myrank == 0) then

    allocate(indices(N_SPATIAL_DOF,commsize))
    indices(1:N_SPATIAL_DOF,1) = local_indices(1:N_SPATIAL_DOF)     

    do iproc=2,commsize
      call MPI_Recv(indices(1:N_SPATIAL_DOF,iproc),N_SPATIAL_DOF, &
                    MPI_INTEGER,iproc-1,MPI_ANY_TAG,PETSC_COMM_WORLD, &
                    status,ierr)
    enddo

  else
    call MPI_Send(local_indices,N_SPATIAL_DOF,MPI_INTEGER, &
                  0,myrank,PETSC_COMM_WORLD,ierr)
  endif

end subroutine ouInitOutput

! ************************************************************************** !
!
! ouInitOutput: Initialize output files
! author: Glenn Hammond
! date: 10/27/03
!
! ************************************************************************** !
subroutine ouInitOutput_old()

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscviewer.h"

  integer :: ierr, iproc
  integer :: local_indices(7), status(MPI_STATUS_SIZE)

  time_filename = 'time.tmp'
  pri_filename = 'pri.tmp'
  mnrl_filename = 'mnrl.tmp'
  napl_filename = 'napl.tmp'  

! Temporary files for binary output during model run (much faster)
  if (myrank == 0) then
    timeunit = 29
    open(unit=timeunit,file=time_filename)
    call PetscBinaryOpen(pri_filename,PETSC_BINARY_CREATE, &
                         tempb_priunit,ierr)
    if (nmnrl > 0) &
      call PetscBinaryOpen(mnrl_filename,PETSC_BINARY_CREATE, &
                           tempb_mnrlunit,ierr)
    if (nnapl > 0) &
      call PetscBinaryOpen(napl_filename,PETSC_BINARY_CREATE, &
                           tempb_naplunit,ierr)
  endif

  local_indices(X_START) = lxs+1
  local_indices(X_END) = lxe
  local_indices(Y_START) = lys+1
  local_indices(Y_END) = lye
  local_indices(Z_START) = lzs+1
  local_indices(Z_END) = lze
  local_indices(N_SPATIAL_DOF) = lnx*lny*lnz

  ! obtain indices of other processors
  if (myrank == 0) then

    allocate(indices(N_SPATIAL_DOF,commsize))
    indices(1:N_SPATIAL_DOF,1) = local_indices(1:N_SPATIAL_DOF)     

    do iproc=2,commsize
      call MPI_Recv(indices(1:N_SPATIAL_DOF,iproc),N_SPATIAL_DOF, &
                    MPI_INTEGER,iproc-1,MPI_ANY_TAG,PETSC_COMM_WORLD, &
                    status,ierr)
    enddo

  else
    call MPI_Send(local_indices,N_SPATIAL_DOF,MPI_INTEGER, &
                  0,myrank,PETSC_COMM_WORLD,ierr)
  endif

end subroutine ouInitOutput_old

! ************************************************************************** !
!
! ouPostProcessOutput: Simplifies partran.F90
! author: Glenn Hammond
! date: 12/12/02
!
! ************************************************************************** !
subroutine ouPostProcessResults(da_sdof)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"

  DA :: da_sdof

  character(len=wordlen) :: filename
  integer :: ierr
  real*8, allocatable :: global_del_x(:), global_del_y(:), global_del_z(:)


  call ouComputeGlobalGridSpacing(da_sdof,global_del_x,global_del_y, &
                                  global_del_z)

  if ((ny == 1 .and. nz == 1) .or. (nx == 1 .and. nz == 1) .or. &
      (nx == 1 .and. ny == 1)) then

    if (myrank == 0) print *, 'Printing EXCEL Dataset'

    filename = 'excel_new.dat'
    call ouPrintExcelResults_new('data.tmp',filename,global_del_x, &
                             global_del_y,global_del_z)

#if 0
    filename = 'excel_pri.dat'
    call ouPrintExcelResults(pri_filename,filename,COMPS,global_del_x, &
                             global_del_y,global_del_z)
    if (nmnrl > 0) then
      filename = 'excel_mnrl.dat'
      call ouPrintExcelResults(pri_filename,filename,MNRLS,global_del_x, &
                               global_del_y,global_del_z)
    endif
    if (nnapl > 0) then
      filename = 'excel_napl.dat'
      call ouPrintExcelResults(pri_filename,filename,NAPLS,global_del_x, &
                               global_del_y,global_del_z)
    endif
#endif
  else

    if (print_gms) then
      if (myrank == 0) then
        print *, 'Printing GMS grid file.'
        call ouPrintGMSGrid(global_del_x,global_del_y,global_del_z)
      endif

      if (myrank == 0) print *, 'Printing GMS Dataset'

      filename = 'new_pri.dat'
      call ouPrintGMSResults_new('data.tmp',filename)

#if 0

      filename = 'gms_pri.dat'
      call ouPrintGMSResults(pri_filename,filename,COMPS)

      if (nmnrl > 0) then
        filename = 'gms_mnrl.dat'
        call ouPrintGMSResults(mnrl_filename,filename,MNRLS)
      endif
      if (nnapl > 0) then
        filename = 'gms_napl.dat'
        call ouPrintGMSResults(napl_filename,filename,NAPLS)
      endif
#endif
    endif
  endif

  if (allocated(global_del_x)) &
    deallocate(global_del_x,global_del_y,global_del_z)

  if (myrank == 0) print *, 'Total nodes', nx*ny*nz

end subroutine ouPostProcessResults

! ************************************************************************** !
!
! ouPrintResults: Prints vectors to associated files
! author: Glenn Hammond
! date: 10/27/03
!
! ************************************************************************** !
subroutine ouPrintResults(da_mdof, da_sdof, Gconc, Gmnrl, Gnapl, time, ts)

  use globals, only : myrank, nmnrl, nnapl

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  DA :: da_mdof, da_sdof
  Vec :: Gconc, Gmnrl, Gnapl, temp
  integer :: ts
  real*8 :: time

  integer :: dummy

  call ouPrintConcentrationsBinary(time,da_sdof,temp,temp,temp, &
                                   temp,Gconc,Gmnrl,temp,temp,temp)
  return

  call ouPrintVectorBinary_old(tempb_priunit,da_mdof,da_sdof,Gconc)
  if (nmnrl > 0) call ouPrintVectorBinary_old(tempb_mnrlunit,da_mdof, &
                                          da_sdof,Gmnrl)
  if (nnapl > 0) call ouPrintVectorBinary_old(tempb_naplunit,da_mdof, &
                                          da_sdof,Gnapl)

100   format('Time: ', es20.8, '  Time Step: ', i5)

end subroutine ouPrintResults

! ************************************************************************** !
!
!
! ************************************************************************** !
subroutine ouPrintConcentrationsBinary(time,da_sdof, &
                                       Gflow, Gpermeability, Gporosity, &
                                       Gsaturation, Gconc, Gmnrl_conc, &
                                       Gmnrl_rates, Gnapl_conc, Gnapl_rates)

  use globals
  use fileio
  use eqrxn

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"


  DA :: da_sdof
  Vec :: Gflow, Gpermeability, Gporosity, Gsaturation
  Vec :: Gconc
  Vec :: Gmnrl_conc, Gmnrl_rates
  Vec :: Gnapl_conc, Gnapl_rates
  integer :: binaryunit
  real*8 :: time

  Vec :: global_sdof_vec, seq_vec
  character(len=namlen) :: allchar
  real*8, pointer :: vec_ptr(:)
  integer :: i, idof, ierr, num_spatial_dof, num_dof
  integer :: icomp, icplx, imnrl
#ifndef MPIUNI
  integer :: status(MPI_STATUS_SIZE)
#endif
  type(list_type), pointer :: cur_component

  num_spatial_dof = nx*ny*nz
  allchar = 'all'

  if (myrank == 0 .and. commsize > 1) then
    call VecCreateSeq(PETSC_COMM_SELF,num_spatial_dof,seq_vec,ierr)
  endif
  call DACreateGlobalVector(da_sdof,global_sdof_vec,ierr)

! secondary components
  call VecGetBlockSize(Gconc,num_dof,ierr)

  cur_component => primary_print_list
  if (associated(cur_component)) then
    if (fiStringCompare(cur_component%name,allchar,namlen)) then
      do icomp=1,ncomp
        call VecStrideGather(Gconc,icomp-1,global_sdof_vec, &
                             INSERT_VALUES,ierr)
        call ouPrintSDOFVector(global_sdof_vec,seq_vec,num_spatial_dof, &
                               time,pri_names(icomp))
      enddo
    else
      do

        if (.not.associated(cur_component)) exit
 
        do icomp=1,ncomp
          if(fiStringCompare(pri_names(icomp),cur_component%name,namlen)) then
            call VecStrideGather(Gconc,icomp-1,global_sdof_vec, &
                                 INSERT_VALUES,ierr)
            call ouPrintSDOFVector(global_sdof_vec,seq_vec,num_spatial_dof, &
                                   time,cur_component%name)
            exit
          endif
        enddo

        cur_component => primary_print_list%next

      enddo
    endif
  endif

  cur_component => secondary_print_list
  if (associated(cur_component)) then 
    if (fiStringCompare(cur_component%name,allchar,namlen)) then
      do icplx=1,ncplx
        call eqComputeComplexConcentration(Gconc,global_sdof_vec,icplx)
        call ouPrintSDOFVector(global_sdof_vec,seq_vec,num_spatial_dof, &
                               time,sec_names(icplx))
      enddo
    else
      do

        if (.not.associated(cur_component)) exit

        do icplx=1,ncplx
          if(fiStringCompare(sec_names(icplx),cur_component%name,namlen)) then
            call eqComputeComplexConcentration(Gconc,global_sdof_vec,icplx)
            call ouPrintSDOFVector(global_sdof_vec,seq_vec,num_spatial_dof, &
                                   time,cur_component%name)
            exit
          endif
        enddo

        cur_component => cur_component%next

      enddo
    endif
  endif

  cur_component => mineral_print_list
  if (associated(cur_component)) then
    if (fiStringCompare(cur_component%name,allchar,namlen)) then
      do imnrl=1,nmnrl
        call VecStrideGather(Gmnrl_conc,imnrl-1,global_sdof_vec, &
                             INSERT_VALUES,ierr)
        call ouPrintSDOFVector(global_sdof_vec,seq_vec,num_spatial_dof, &
                               time,mnrl_names(imnrl))
      enddo
    else
      do

        if (.not.associated(cur_component)) exit
 
        do imnrl=1,nmnrl
          if(fiStringCompare(mnrl_names(imnrl),cur_component%name,namlen)) then
            call VecStrideGather(Gmnrl_conc,imnrl-1,global_sdof_vec, &
                                 INSERT_VALUES,ierr)
            call ouPrintSDOFVector(global_sdof_vec,seq_vec,num_spatial_dof, &
                                   time,cur_component%name)
            exit
          endif
        enddo

        cur_component => cur_component%next

      enddo
    endif
  endif

  if (myrank == 0 .and. commsize > 1) then
    call VecDestroy(seq_vec,ierr)
  endif
  call VecDestroy(global_sdof_vec,ierr)

end subroutine ouPrintConcentrationsBinary

! ************************************************************************** !
!
! ouPrintSDOFVector: Converts a global single dof vector to natural format 
!                    and prints it in binary format using PETSc binary.
! author: Glenn Hammond
! date: 4/29/04
!
! ************************************************************************** !
subroutine ouPrintSDOFVector(global_sdof_vec,seq_vec,n,time,name)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  Vec :: global_sdof_vec, seq_vec
  integer :: n
  real*8 :: time
  character(len=namlen) :: name

  integer :: ierr
  real*8, pointer :: vec_ptr(:)

  if (myrank == 0) then
    write(timeunit,111) name, time, nx, ny, nz
  endif

  if (commsize > 1) then
    call ouVecMPIToSeqZero(global_sdof_vec,seq_vec)
    if (myrank == 0) then
      call VecGetArrayF90(seq_vec,vec_ptr,ierr)   
      call PetscBinaryWrite(tempb_unit,vec_ptr,n,PETSC_DOUBLE,0,ierr)
      call VecRestoreArrayF90(seq_vec,vec_ptr,ierr)
    endif
  else
    call VecGetArrayF90(global_sdof_vec,vec_ptr,ierr)   
    call PetscBinaryWrite(tempb_unit,vec_ptr,n,PETSC_DOUBLE,0,ierr)
    call VecRestoreArrayF90(global_sdof_vec,vec_ptr,ierr)
  endif

111   format(a20,es20.8,i5,i5,i5)

end subroutine ouPrintSDOFVector

! ************************************************************************** !
!
! ouPrintVectorBinary: Converts a global vector to natural format and prints
!                      it in binary format using PETSc binary.
! author: Glenn Hammond
! date: 10/27/03
!
! ************************************************************************** !
subroutine ouPrintVectorBinary_old(binaryunit, da_mdof, da_sdof, global_vec)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"


  DA :: da_mdof, da_sdof
  Vec :: global_vec
  integer :: binaryunit, debugunit

  Vec :: global_sdof_vec, seq_vec
  real*8, pointer :: vec_ptr(:)
  integer :: i, idof, ierr, num_spatial_dof, num_dof
#ifndef MPIUNI
  integer :: status(MPI_STATUS_SIZE)
#endif

  num_spatial_dof = nx*ny*nz
  call VecGetBlockSize(global_vec,num_dof,ierr)


  if (myrank == 0 .and. commsize > 1) then
    call VecCreateSeq(PETSC_COMM_SELF,num_spatial_dof,seq_vec,ierr)
  endif

  if (num_dof > 1) then
    call DACreateGlobalVector(da_sdof,global_sdof_vec,ierr)
    do idof=1,num_dof

      call VecStrideGather(global_vec,idof-1,global_sdof_vec, &
                           INSERT_VALUES,ierr)
      if (commsize > 1) then
        call ouVecMPIToSeqZero(global_sdof_vec,seq_vec)
        if (myrank == 0) then
          call VecGetArrayF90(seq_vec,vec_ptr,ierr)   
          call PetscBinaryWrite(binaryunit,vec_ptr,num_spatial_dof, &
                                PETSC_DOUBLE,0,ierr)
          call VecRestoreArrayF90(seq_vec,vec_ptr,ierr)
        endif
      else
        call VecGetArrayF90(global_sdof_vec,vec_ptr,ierr)   
        call PetscBinaryWrite(binaryunit,vec_ptr,num_spatial_dof, &
                              PETSC_DOUBLE,0,ierr)
        call VecRestoreArrayF90(global_sdof_vec,vec_ptr,ierr)
      endif
    enddo
    call VecDestroy(global_sdof_vec,ierr)
  else
    if (commsize > 1) then
      call ouVecMPIToSeqZero(global_vec,seq_vec)
      if (myrank == 0) then
        call VecGetArrayF90(seq_vec,vec_ptr,ierr)   
        call PetscBinaryWrite(binaryunit,vec_ptr,num_spatial_dof, &
                              PETSC_DOUBLE,0,ierr)
        call VecRestoreArrayF90(seq_vec,vec_ptr,ierr)
      endif
    else
      call VecGetArrayF90(global_vec,vec_ptr,ierr)   
      call PetscBinaryWrite(binaryunit,vec_ptr,num_spatial_dof, &
                            PETSC_DOUBLE,0,ierr)
      call VecRestoreArrayF90(global_vec,vec_ptr,ierr)
    endif
  endif

  if (myrank == 0 .and. commsize > 1) then
    call VecDestroy(seq_vec,ierr)
  endif

end subroutine ouPrintVectorBinary_old

! ************************************************************************** !
!
! ouPrintGMSVector: Converts a 1 dof global vector to natural format and 
!                   prints it in ascii format using PETSc binary.
! author: Glenn Hammond
! date: 10/27/03
!
! ************************************************************************** !
subroutine ouPrintGMSVector(global_vec,filename,dataset_name)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"


  Vec :: global_vec
  character(len=namlen) :: filename, dataset_name

  Vec :: seq_vec
  real*8, pointer :: vec_ptr(:)
  integer :: outunit, num_spatial_dof, ierr

  num_spatial_dof = nx*ny*nz

  outunit = 86

  if (myrank == 0) open(unit=outunit,file=filename)

  if (myrank == 0 .and. commsize > 1) then
    call VecCreateSeq(PETSC_COMM_SELF,num_spatial_dof,seq_vec,ierr)
  endif

  if (commsize > 1) then
    call ouVecMPIToSeqZero(global_vec,seq_vec)
    if (myrank == 0) then
      call VecGetArrayF90(seq_vec,vec_ptr,ierr)   
      call ouPrintGMSDataSetHeader(outunit,dataset_name,num_spatial_dof)
      call ouPrintGMSTimeStep(outunit,0.d0,vec_ptr,num_spatial_dof)
      call ouPrintGMSDataSetEnd(outunit)
      call VecRestoreArrayF90(seq_vec,vec_ptr,ierr)
    endif
  else
    call VecGetArrayF90(global_vec,vec_ptr,ierr)   
    call ouPrintGMSDataSetHeader(outunit,dataset_name,num_spatial_dof)
    call ouPrintGMSTimeStep(outunit,0.d0,vec_ptr,num_spatial_dof)
    call ouPrintGMSDataSetEnd(outunit)
    call VecRestoreArrayF90(global_vec,vec_ptr,ierr)
  endif

  if (myrank == 0 .and. commsize > 1) then
    call VecDestroy(seq_vec,ierr)
  endif

  if (myrank == 0) close(outunit)

end subroutine ouPrintGMSVector

! ************************************************************************** !
!
! ouFinalizeOutput: Close output files
! author: Glenn Hammond
! date: 10/27/03
!
! ************************************************************************** !
subroutine ouFinalizeOutput()

  use globals, only : myrank, nmnrl, nnapl

  implicit none

#include "include/finclude/petsc.h"

  integer :: ierr

  if (myrank == 0) then
    close(timeunit)
    call PetscBinaryClose(tempb_unit,ierr)
    call PetscBinaryClose(tempb_priunit,ierr)
    if (nmnrl > 0) call PetscBinaryClose(tempb_mnrlunit,ierr)
    if (nnapl > 0) call PetscBinaryClose(tempb_naplunit,ierr)
  endif

end subroutine ouFinalizeOutput

! ************************************************************************** !
!
! ouDeallocateOutput: Deallocate internal file globals
! author: Glenn Hammond
! date: 12/02/03
!
! ************************************************************************** !
subroutine ouDeallocateOutput()

  implicit none

  if (allocated(output_buffer)) deallocate(output_buffer)
  if (allocated(indices)) deallocate(indices)

end subroutine ouDeallocateOutput

! ************************************************************************** !
!
! ouPrintPETScMatrix
! author: Glenn Hammond
! date: 3/8/01
!
! ************************************************************************** !
subroutine ouPrintPETScMatrix(mat,filename)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscviewer.h"

  Mat      mat
  PetscViewer   matviewer
  character(len=wordlen) :: filename

  integer :: ierr

  call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)

  if (len_trim(filename) == 0) then
    call MatView(mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
  else
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,matviewer,ierr)
    call PetscViewerSetFormat(matviewer,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call MatView(mat,matviewer,ierr)
    call PetscViewerDestroy(matviewer,ierr)
  endif

end subroutine ouPrintPETScMatrix

! ************************************************************************** !
!
! ouPrintPETScVector
! author: Glenn Hammond
! date: 3/8/01
!
! ************************************************************************** !
subroutine ouPrintPETScVector(vec,filename)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscviewer.h"

  Vec      vec
  PetscViewer   vecviewer
  character(len=wordlen) :: filename

  integer :: ierr

  if (len_trim(filename) == 0) then
    call VecView(vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
  else
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,vecviewer,ierr)
    call PetscViewerSetFormat(vecviewer,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call VecView(vec,vecviewer,ierr)
    call PetscViewerDestroy(vecviewer,ierr)
  endif

end subroutine ouPrintPETScVector

! ************************************************************************** !
!
! ouPrintL_GConversionArrays
! author: Glenn Hammond
! date: 1/23/01
!
! ************************************************************************** !
subroutine ouPrintL_GConversionArrays(outunit)

  use globals

  implicit none

  integer :: outunit

  integer :: i

  write(outunit,*)
  write(outunit,*) 'local to ghosted'
  do i=1,lnum_nodes
    write(outunit,*) 'local: ', i, '   ghosted: ', node_id_ltog(i)
  enddo
  write(outunit,*)
  write(outunit,*)

  write(outunit,*) 'ghosted to local'
  do i=1,gnum_nodes
    write(outunit,*) 'ghosted: ', i, '   local: ', node_id_gtol(i)
  enddo
  write(outunit,*)

end subroutine ouPrintL_GConversionArrays

! ************************************************************************** !
!
! ouPrintDomains: Prints a GMS grid with the domain decomposition mapped
!                      using different materials
! author: Glenn Hammond
! date: 7/19/00
!
! ************************************************************************** !
subroutine ouPrintGMSGrid(del_x,del_y,del_z)

  use globals

  implicit none

  real*8 :: del_x(*), del_y(*), del_z(*)

  integer :: i, outunit, ierr, n
  real*8 :: tempreal

  outunit = 42

  if (nz == 1) then

    n = nx*ny
    open(outunit, file='partran.2dg')
  
    write(outunit,10) !'GRID2D'
    write(outunit,20) !'TYPE 1'
    write(outunit,30) !'IJ +y +x'
    write(outunit,40) 0.d0, 0.d0, 0.d0 !'ORIGIN '
    write(outunit,50) !'ROTZ 0'
    write(outunit,60) nx+1, ny+1 !'DIM'

    tempreal = 0.d0
    write(outunit,90) tempreal
    do i=1,nx
      tempreal = tempreal + del_x(i)
      write(outunit,90) tempreal
    enddo

    tempreal = 0.d0
    write(outunit,90) tempreal
    do i=1,ny
      tempreal = tempreal + del_y(i)
      write(outunit,90) tempreal
    enddo

  else

    n = nx*ny*nz
    open(outunit, file='partran.3dg')
  
    write(outunit,15) !'GRID3D'
    write(outunit,20) !'TYPE 1'
    write(outunit,35) !'IJ +y +x'
    write(outunit,40) 0.d0, 0.d0, 0.d0 !'ORIGIN'
    write(outunit,50) !'ROTZ 0'
    write(outunit,65) nx+1, ny+1, nz+1 !'DIM'

    tempreal = 0.d0
    write(outunit,90) tempreal
    do i=1,nx
      tempreal = tempreal + del_x(i)
      write(outunit,90) tempreal
    enddo

    tempreal = 0.d0
    write(outunit,90) tempreal
    do i=1,ny
      tempreal = tempreal + del_y(i)
      write(outunit,90) tempreal
    enddo

    tempreal = 0.d0
    write(outunit,90) tempreal
    do i=1,nz
      tempreal = tempreal + del_z(i)
      write(outunit,90) tempreal
    enddo

  endif

  write(outunit,95) !'DMAT 1'

  close(outunit)

10 format('GRID2D')
15 format('GRID3D')
20 format('TYPE 1')
30 format('IJ +y +x')
35 format('IJK +y +x +z')
40 format('ORIGIN',f12.4,f12.4,f12.4)
50 format('ROTZ 0')
60 format('DIM',i7,i7)
65 format('DIM',i7,i7,i7)
90 format(f12.4)
95 format('DMAT 1')
100 format('Writing of GMS grid completed!')

end subroutine ouPrintGMSGrid

! ************************************************************************** !
!
! ouReadTimes: Reads in output times from the time.tmp file
! author: Glenn Hammond
! date: 1/31/01, 10/28/03
!
! ************************************************************************** !
subroutine ouReadTimes(num_times,print_times)

  use globals, only: strlen, wordlen
  use fileio

  implicit none

  integer :: num_times
  real*8, allocatable :: print_times(:)

  integer :: inunit, ierr, cur_time
  character(len=strlen) :: string
  character(len=wordlen) :: word

  inunit = 41
  open(inunit, file='time.tmp')
  num_times = 0
  do 
    call fiReadString(inunit,string,ierr)
    call fiReadWord(string,word,.true.,ierr)
    if (ierr /= 0) exit
    num_times = num_times + 1
  enddo
  allocate(print_times(num_times))
  rewind(inunit)
  do cur_time=1,num_times
    call fiReadString(inunit,string,ierr)
    call fiReadWord(string,word,.true.,ierr)
    call fiReadDouble(string,print_times(cur_time),ierr)
  enddo
  close(inunit)

end subroutine ouReadTimes

! ************************************************************************** !
!
! ouPrintGMSResults: Prints vector in a standard GMS format
! author: Glenn Hammond
! date: 1/31/01, 10/28/03
!
! ************************************************************************** !
subroutine ouPrintGMSResults(binary_filename,output_filename,data_type)

  use fileio
  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

  character(len=wordlen) :: binary_filename, output_filename
  integer :: data_type

  integer :: outunit, num_times, ierr, cur_time, cur_node, n, &
             binaryinunit, dummyint, offset, idof, ndof
  character(len=namlen) :: dataset_name
  real*8, allocatable :: values(:), print_times(:)
  real*8 :: value

  if (myrank /= 0) return

  outunit = 42

  call ouReadTimes(num_times,print_times)

  n = nx*ny*nz
  allocate(values(n))

  open(outunit,file=output_filename)
  call PetscBinaryOpen(binary_filename,PETSC_BINARY_RDONLY,binaryinunit,ierr)

  if (data_type == COMPS) then
    ndof = ncomp
  elseif(data_type == MNRLS) then
    ndof = nmnrl
  else
    ndof = nnapl
  endif

  do idof=1,ndof

    if (data_type == COMPS) then
      dataset_name = pri_names(idof)
    elseif(data_type == MNRLS) then
      dataset_name = mnrl_names(idof)
    else
      dataset_name = napl_names(idof)
    endif

    call ouPrintGMSDataSetHeader(outunit,dataset_name,n)

    do cur_time=1,num_times
    
      offset = ((cur_time-1)*ndof*n+(idof-1)*n)*PETSC_BINARY_SCALAR_SIZE
      call PetscBinarySeek(binaryinunit,offset,PETSC_BINARY_SEEK_SET, &
                           dummyint,ierr)
      call PetscBinaryRead(binaryinunit,values,n,PETSC_DOUBLE,ierr)

      call ouPrintGMSTimeStep(outunit,print_times(cur_time),values,n)

    enddo
     
    call ouPrintGMSDataSetEnd(outunit)

  enddo
  
  if (data_type == COMPS) then
    write(*,101)
  elseif(data_type == MNRLS) then
    write(*,102)
  else
    write(*,103)
  endif

  call PetscBinaryClose(binaryinunit,ierr)
  close(outunit)

  if (allocated(print_times)) deallocate(print_times)
  if (allocated(values)) deallocate(values)

101 format('Writing of components completed!')
102 format('Writing of minerals completed!')
103 format('Writing of napls completed!')

end subroutine ouPrintGMSResults

! ************************************************************************** !
!
! ouPrintGMSResults: Prints vector in a standard GMS format
! author: Glenn Hammond
! date: 1/31/01, 10/28/03
!
! ************************************************************************** !
subroutine ouPrintGMSResults_new(binary_filename,output_filename)

  use fileio
  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

  character(len=wordlen) :: binary_filename, output_filename

  integer :: inunit, outunit, num_times, ierr, cur_time, cur_node, n, &
             binaryinunit, dummyint, offset, idof, ndof, i, tmp_nx, tmp_ny, &
             tmp_nz
  character(len=namlen) :: dataset_name, first_name
  character(len=strlen) :: string
  real*8, allocatable :: values(:)
  real*8 :: time

  if (myrank /= 0) return

  outunit = 42
  open(outunit,file=output_filename)
  call PetscBinaryOpen(binary_filename,PETSC_BINARY_RDONLY,binaryinunit,ierr)

  inunit = 41
  open(inunit, file='time.tmp')
  call fiReadString(inunit,string,ierr)
  call fiReadNChars(string,first_name,namlen,.true.,ierr)
  ndof = 1
  do
    call fiReadString(inunit,string,ierr)
    call fiReadNChars(string,dataset_name,namlen,.true.,ierr)
    if (ierr /= 0) exit
    if (fiStringCompare(first_name,dataset_name,namlen)) then
      exit
    else
      ndof = ndof + 1
    endif
  enddo


  do idof=1,ndof

    cur_time = 1

    rewind(inunit)
    do i=1,idof-1  ! remove lines before current dof
      call fiReadString(inunit,string,ierr)
    enddo

    do 
      call fiReadString(inunit,string,ierr)
      call fiReadNChars(string,dataset_name,namlen,.true.,ierr)
      if (ierr /= 0) exit
      call fiReadDouble(string,time,ierr)
      call fiReadInt(string,tmp_nx,ierr)
      call fiReadInt(string,tmp_ny,ierr)
      call fiReadInt(string,tmp_nz,ierr)

      n = tmp_nx*tmp_ny*tmp_nz

      if (cur_time == 1) &
        call ouPrintGMSDataSetHeader(outunit,dataset_name,n)

      offset = ((cur_time-1)*ndof*n+(idof-1)*n)*PETSC_BINARY_SCALAR_SIZE
      call PetscBinarySeek(binaryinunit,offset,PETSC_BINARY_SEEK_SET, &
                           dummyint,ierr)
      allocate(values(n))
      call PetscBinaryRead(binaryinunit,values,n,PETSC_DOUBLE,ierr)

      call ouPrintGMSTimeStep(outunit,time,values,n)
      deallocate(values)

      do i=1,ndof-1
        call fiReadString(inunit,string,ierr)
        if (ierr /= 0) exit
      enddo
      if (ierr /= 0) exit
    
      cur_time = cur_time + 1

    enddo
    call ouPrintGMSDataSetEnd(outunit)

  enddo
  
  call PetscBinaryClose(binaryinunit,ierr)
  close(outunit)
  close(inunit)

end subroutine ouPrintGMSResults_new

! ************************************************************************** !
!
! ouPrintGMSDataSetHeader: Prints a data set header in GMS format
! author: Glenn Hammond
! date: 11/18/03
!
! ************************************************************************** !
subroutine ouPrintGMSDataSetHeader(gms_unit,dataset_name,num_nodes)

  use fileio
  use globals

  implicit none

  integer :: gms_unit, num_nodes
  character(len=namlen) :: dataset_name

  write(gms_unit,10) !'DATASET'
  if (nz == 1) then
    write(gms_unit,20) !'OBJTYPE grid2d'
  else
    write(gms_unit,25) !'OBJTYPE grid3d'
  endif

  write(gms_unit,30) !'BEGSCL'  
  write(gms_unit,40) num_nodes
  write(gms_unit,50) num_nodes

  write(gms_unit,60) dataset_name

10 format('DATASET')
20 format('OBJTYPE grid2d')
25 format('OBJTYPE grid3d')
30 format('BEGSCL')
40 format('ND ',i7)
50 format('NC ',i7)
60 format('NAME "',a20,'"')

end subroutine ouPrintGMSDataSetHeader

! ************************************************************************** !
!
! ouPrintGMSTimeStep: Prints timestep portion of GMS data set
! author: Glenn Hammond
! date: 11/18/03
!
! ************************************************************************** !
subroutine ouPrintGMSTimeStep(gms_unit,time,values,num_values)

  use fileio
  use globals

  implicit none

  integer :: gms_unit, num_values
  real*8 :: time, values(num_values)

  integer :: cur_node

  write(gms_unit,70) 0, time

  do cur_node = 1,num_values
    write(gms_unit,75) max(values(cur_node),9.9999d-99)
  enddo

70 format('TS ', i3,es15.4)
75 format(es15.4)

end subroutine ouPrintGMSTimeStep

! ************************************************************************** !
!
! ouPrintGMSDataSetEnd: Prints 'ENDDS' at the end of a data set
! author: Glenn Hammond
! date: 11/18/03
!
! ************************************************************************** !
subroutine ouPrintGMSDataSetEnd(gms_unit)

  implicit none

  integer :: gms_unit

  write(gms_unit,80) 

80 format('ENDDS')

end subroutine ouPrintGMSDataSetEnd

! ************************************************************************** !
!
! ouPrintExcelResults: Prints vector in Excel format
! author: Glenn Hammond
! date: 1/31/01, 10/28/03
!
! ************************************************************************** !
subroutine ouPrintExcelResults(binary_filename,output_filename,data_type, &
                               del_x,del_y,del_z)

  use fileio
  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

  character(len=wordlen) :: binary_filename, output_filename
  integer :: data_type
  real*8 :: del_x(*), del_y(*), del_z(*)

  integer :: inunit, outunit, num_times, ierr, cur_time, cur_node, n, &
             binaryinunit, dummyint, offset, idof, ndof, xyz
  character(len=namlen) :: dataset_name
  real*8, allocatable :: values(:), print_times(:)
  real*8 :: value

  if (myrank /= 0) return

  outunit = 42

  call ouReadTimes(num_times,print_times)
  call ouCalculateXYZ(xyz)

  n = nx*ny*nz

  open(outunit,file=output_filename)
  call PetscBinaryOpen(binary_filename,PETSC_BINARY_RDONLY,binaryinunit,ierr)

  if (data_type == COMPS) then
    ndof = ncomp
  elseif(data_type == MNRLS) then
    ndof = nmnrl
  else
    ndof = nnapl
  endif

  allocate(values(n*ndof))

  do cur_time=1,num_times

    if (data_type == COMPS) then
      call ouPrintTabularDataHeader(outunit,print_times(cur_time),ncomp, &
                                    pri_names,xyz)
    elseif(data_type == MNRLS) then
      call ouPrintTabularDataHeader(outunit,print_times(cur_time),nmnrl, &
                                    mnrl_names,xyz)
    else
      call ouPrintTabularDataHeader(outunit,print_times(cur_time),nnapl, &
                                    napl_names,xyz)
    endif

!    offset = ((cur_time-1)*ndof*n*PETSC_BINARY_SCALAR_SIZE
!    call PetscBinarySeek(binaryinunit,offset,PETSC_BINARY_SEEK_SET, &
!                         dummyint,ierr)
    call PetscBinaryRead(binaryinunit,values,n,PETSC_DOUBLE,ierr)

    call ouPrintTabularData(outunit,ndof,values,del_x,del_y,del_z,xyz)

    write(outunit,*)
    write(outunit,*)

  enddo
  
  if (data_type == COMPS) then
    write(*,201)
  elseif(data_type == MNRLS) then
    write(*,202)
  else
    write(*,203)
  endif

  call PetscBinaryClose(binaryinunit,ierr)
  close(outunit)

  if (allocated(print_times)) deallocate(print_times)
  if (allocated(values)) deallocate(values)

201 format('Writing of components completed!')
202 format('Writing of minerals completed!')
203 format('Writing of napls completed!')

end subroutine ouPrintExcelResults


! ************************************************************************** !
!
! ouPrintExcelResults: Prints vector in Excel format
! author: Glenn Hammond
! date: 1/31/01, 10/28/03
!
! ************************************************************************** !
subroutine ouPrintExcelResults_new(binary_filename,output_filename, &
                                   del_x,del_y,del_z)

  use fileio
  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

  character(len=wordlen) :: binary_filename, output_filename
  integer :: data_type
  real*8 :: del_x(*), del_y(*), del_z(*)

  integer :: inunit, outunit, num_times, dummyint, ierr, n, &
             binaryinunit, offset, idof, ndof, xyz
  integer :: itime, ntimes, tmp_nx, tmp_ny, tmp_nz
  character(len=strlen) :: string
  character(len=namlen) :: dataset_name, first_name
  character(len=namlen), allocatable :: names(:)
  real*8, allocatable :: values(:), print_times(:)

  if (myrank /= 0) return

  call ouCalculateXYZ(xyz)

  outunit = 42
  open(outunit,file=output_filename)
  call PetscBinaryOpen(binary_filename,PETSC_BINARY_RDONLY,binaryinunit,ierr)
  if (ierr /= 0) return

  inunit = 41
  open(inunit, file='time.tmp')
  call fiReadString(inunit,string,ierr)
  call fiReadNChars(string,first_name,namlen,.true.,ierr)
  ndof = 1
  ntimes = 1
  do
    call fiReadString(inunit,string,ierr)
    call fiReadNChars(string,dataset_name,namlen,.true.,ierr)
    if (ierr /= 0) exit
    if (fiStringCompare(first_name,dataset_name,namlen)) then
      ntimes = ntimes + 1
    else
      if (ntimes == 1) ndof = ndof + 1
    endif
  enddo

  allocate(names(ndof))
  allocate(print_times(ntimes))

  rewind(inunit)
  itime = 0
  idof = 0
  do
    call fiReadString(inunit,string,ierr)
    call fiReadNChars(string,dataset_name,namlen,.true.,ierr)
    if (ierr /= 0) exit
    if (fiStringCompare(first_name,dataset_name,namlen)) then
      itime = itime + 1
      call fiReadDouble(string,print_times(itime),ierr)
      if (itime == 1) then
        call fiReadInt(string,tmp_nx,ierr)
        call fiReadInt(string,tmp_ny,ierr)
        call fiReadInt(string,tmp_nz,ierr)
        n = tmp_nx*tmp_ny*tmp_nz
      endif
    endif
    if (itime == 1) then
      idof = idof + 1
      names(idof) = dataset_name
    endif
  enddo

  allocate(values(n*ndof))

  do itime=1,ntimes

    call ouPrintTabularDataHeader(outunit,print_times(itime),ndof, &
                                  names,xyz)

    offset = ((itime-1)*ndof*n)*PETSC_BINARY_SCALAR_SIZE
    call PetscBinarySeek(binaryinunit,offset,PETSC_BINARY_SEEK_SET, &
                         dummyint,ierr)
    call PetscBinaryRead(binaryinunit,values,ndof*n,PETSC_DOUBLE,ierr)
    call ouPrintTabularData_new(outunit,ndof,values,del_x,del_y,del_z,xyz, &
                            tmp_nx,tmp_ny,tmp_nz)

  enddo
  
  call PetscBinaryClose(binaryinunit,ierr)
  close(outunit)
  close(inunit)

  if (allocated(names)) deallocate(names)
  if (allocated(print_times)) deallocate(print_times)
  if (allocated(values)) deallocate(values)

end subroutine ouPrintExcelResults_new

! ************************************************************************** !
!
! ouPrintTabularDataHeader: Prints time and component names in tabular format
! author: Glenn Hammond
! date: 12/02/03
!
! ************************************************************************** !
subroutine ouPrintTabularDataHeader(outunit,time,ndof,names,xyz)

  use globals, only: namlen
  
  implicit none

  integer :: ndof, outunit, xyz
  character(len=namlen) :: names(ndof)
  character(len=3) :: tab
  real*8 :: time

  integer :: i

  tab = achar(9)

  write(outunit,100) tab, time
!  write(outunit,110) (names(i),i=1,ndof)

  select case(xyz)
    case(X_DIM)
      write(outunit,110) 'X', tab, (names(i),tab,i=1,ndof)
    case(Y_DIM)
      write(outunit,110) 'Y, tab', (names(i),tab,i=1,ndof)
    case(Z_DIM)
      write(outunit,110) 'Z, tab', (names(i),tab,i=1,ndof)
    case(XY_DIM)
      write(outunit,120) 'X, tab', 'Y', tab, (names(i),tab,i=1,ndof)
    case(YZ_DIM)
      write(outunit,120) 'Y, tab', 'Z', tab, (names(i),tab,i=1,ndof)
    case(XZ_DIM)
      write(outunit,120) 'X, tab', 'Z', tab, (names(i),tab,i=1,ndof)
    case(XYZ_DIM)
      write(outunit,130) 'X, tab', 'Y', tab, 'Z', tab, (names(i),tab,i=1,ndof)
  end select


100 format(' Time: ',a2,es15.6)
110 format(a3,a2,100(a20,a2))
120 format(2(a3,a2),100(a20,a2))
130 format(3(a3,a2),100(a20,a2))

end subroutine ouPrintTabularDataHeader

! ************************************************************************** !
!
! ouPrintTabularData: Prints 1, 2 or 3D data in a tabular format
! author: Glenn Hammond
! date: 12/02/03
!
! ************************************************************************** !
subroutine ouPrintTabularData(outunit,ndof,values,del_x,del_y,del_z,xyz)

  use globals, only : nx, ny, nz
  implicit none

  integer :: outunit, xyz, ndof
  real*8 :: values(*), del_x(*), del_y(*), del_z(*)

  integer :: count, i, j, k
  real*8 :: location_x, location_y, location_z

  count = 1
  do k=1,nz
    if (k==1) then
      location_z = 0.5d0*del_z(count)
    else
      location_z = location_z+0.5d0*(del_z(count-1)+del_z(count))
    endif
    do j=1,ny
      if (j==1) then
        location_y = 0.5d0*del_y(count)
      else
        location_y = location_y+0.5d0*(del_y(count-1)+del_y(count))
      endif
      do i=1,nx
        if (i==1) then
          location_x = 0.5d0*del_x(count)
        else
          location_x = location_x+0.5d0*(del_x(count-1)+del_x(count))
        endif

        select case(xyz)
          case(X_DIM)
            write(outunit,100) location_x, &
                               values((count-1)*ndof+1:count*ndof)
          case(Y_DIM)
            write(outunit,100) location_y, &
                               values((count-1)*ndof+1:count*ndof)
          case(Z_DIM)
            write(outunit,100) location_z, &
                               values((count-1)*ndof+1:count*ndof)
          case(XY_DIM)
            write(outunit,110) location_x, location_y, &
                               values((count-1)*ndof+1:count*ndof)
          case(YZ_DIM)
            write(outunit,110) location_y, location_z, &
                               values((count-1)*ndof+1:count*ndof)
          case(XZ_DIM)
            write(outunit,110) location_x, location_z, &
                               values((count-1)*ndof+1:count*ndof)
          case(XYZ_DIM)
            write(outunit,120) location_x, location_y, location_z, &
                               values((count-1)*ndof+1:count*ndof)
        end select 

        count = count + 1
      enddo
    enddo
  enddo

100 format(es12.5,100(es18.6))
110 format(2(es12.5),100(es18.6))
120 format(3(es12.5),100(es18.6))

end subroutine ouPrintTabularData

! ************************************************************************** !
!
! ouPrintTabularData: Prints 1, 2 or 3D data in a tabular format
! author: Glenn Hammond
! date: 12/02/03
!
! ************************************************************************** !
subroutine ouPrintTabularData_new(outunit,ndof,values,del_x,del_y,del_z, &
                                  xyz,nx,ny,nz)

  implicit none

  integer :: outunit, xyz, ndof, nx, ny, nz, n
  real*8 :: values(*), del_x(*), del_y(*), del_z(*)

  integer :: count, i, j, k, idof
  real*8 :: location_x, location_y, location_z
  character(len=1) :: tab

  tab = achar(9)

  n = nx*ny*nz

  count = 1
  do k=1,nz
    if (k==1) then
      location_z = 0.5d0*del_z(count)
    else
      location_z = location_z+0.5d0*(del_z(count-1)+del_z(count))
    endif
    do j=1,ny
      if (j==1) then
        location_y = 0.5d0*del_y(count)
      else
        location_y = location_y+0.5d0*(del_y(count-1)+del_y(count))
      endif
      do i=1,nx
        if (i==1) then
          location_x = 0.5d0*del_x(count)
        else
          location_x = location_x+0.5d0*(del_x(count-1)+del_x(count))
        endif

        select case(xyz)
          case(X_DIM)
!            write(outunit,100) location_x, &
!                               values((count-1)*ndof+1:count*ndof)
            write(outunit,100) location_x, tab, &
                               (values(count+idof-1),tab,idof=1,ndof*n,n)
          case(Y_DIM)
            write(outunit,100) location_y, tab, &
                               (values(count+idof-1),tab,idof=1,ndof*n,n)
          case(Z_DIM)
            write(outunit,100) location_z, tab, &
                               (values(count+idof-1),tab,idof=1,ndof*n,n)
          case(XY_DIM)
            write(outunit,110) location_x, tab, location_y, tab, &
                               (values(count+idof-1),tab,idof=1,ndof*n,n)
          case(YZ_DIM)
            write(outunit,110) location_y, tab, location_z, tab, &
                               (values(count+idof-1),tab,idof=1,ndof*n,n)
          case(XZ_DIM)
            write(outunit,110) location_x, tab, location_z, tab, &
                               (values(count+idof-1),tab,idof=1,ndof*n,n)
          case(XYZ_DIM)
            write(outunit,120) location_x, tab, location_y, tab, location_z, tab, &
                               (values(count+idof-1),tab,idof=1,ndof*n,n)
        end select 

        count = count + 1
      enddo
    enddo
  enddo

  write(outunit,*) ! put blank line between data sets

100 format(es12.5,a2,100(es18.6,a2))
110 format(2(es12.5,a2),100(es18.6,a2))
120 format(3(es12.5,a2),100(es18.6,a2))

end subroutine ouPrintTabularData_new

! ************************************************************************** !
!
! ouComputeGlobalGridSpacing: Computes coordinates of grid nodes/cells
! author: Glenn Hammond
! date: 12/02/03
!
! ************************************************************************** !
subroutine ouComputeGlobalGridSpacing(da_sdof,del_x,del_y,del_z)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  real*8, allocatable :: del_x(:), del_y(:), del_z(:)

  DA :: da_sdof
  Vec :: global_sdof_vec, seq_vec
  integer :: n, ierr, size, cur_lnode, gnode_id

  PetscScalar, pointer :: vec_ptr(:)

  call PetscBarrier(da_sdof,ierr)

  if (.not.allocated(indices).and.myrank == 0) then
    print *, 'ouInitOutput needs to be called before ', &
             'ouComputeGlobalGridSpacing'
    stop
  endif

  n = nx*ny*nz
  if (myrank == 0) allocate(del_x(n),del_y(n),del_z(n))

  if (.not. uniform_grid) then
    call DACreateGlobalVector(da_sdof,global_sdof_vec,ierr)

    call VecGetSize(global_sdof_vec,size,ierr)

    if (myrank == 0 .and. commsize > 1) then
      call VecCreateSeq(PETSC_COMM_SELF,size,seq_vec,ierr)
    endif

! dist x
    call VecGetArrayF90(global_sdof_vec,vec_ptr,ierr)
    do cur_lnode = 1,lnum_nodes
      gnode_id = node_id_ltog(cur_lnode)
      vec_ptr(cur_lnode) = dx(gnode_id)*0.1d0
    enddo
    call VecRestoreArrayF90(global_sdof_vec,vec_ptr,ierr)

    call ouVecMPIToSeqZero(global_sdof_vec,seq_vec)
    
    if (myrank == 0) then
      call VecGetArrayF90(seq_vec,vec_ptr,ierr)
      del_x(1:n) = vec_ptr(1:n)
      call VecRestoreArrayF90(seq_vec,vec_ptr,ierr)
    endif

! dist y
    call VecGetArrayF90(global_sdof_vec,vec_ptr,ierr)
    do cur_lnode = 1,lnum_nodes
      gnode_id = node_id_ltog(cur_lnode)
      vec_ptr(cur_lnode) = dy(gnode_id)*0.1d0
    enddo
    call VecRestoreArrayF90(global_sdof_vec,vec_ptr,ierr)

    call ouVecMPIToSeqZero(global_sdof_vec,seq_vec)
    
    if (myrank == 0) then
      call VecGetArrayF90(seq_vec,vec_ptr,ierr)
      del_y(1:n) = vec_ptr(1:n)
      call VecRestoreArrayF90(seq_vec,vec_ptr,ierr)
    endif

! dist z
    call VecGetArrayF90(global_sdof_vec,vec_ptr,ierr)
    do cur_lnode = 1,lnum_nodes
      gnode_id = node_id_ltog(cur_lnode)
      vec_ptr(cur_lnode) = dz(gnode_id)*0.1d0
    enddo
    call VecRestoreArrayF90(global_sdof_vec,vec_ptr,ierr)

    call ouVecMPIToSeqZero(global_sdof_vec,seq_vec)
    
    if (myrank == 0) then
      call VecGetArrayF90(seq_vec,vec_ptr,ierr)
      del_z(1:n) = vec_ptr(1:n)
      call VecRestoreArrayF90(seq_vec,vec_ptr,ierr)
    endif

    call VecDestroy(global_sdof_vec,ierr)
    call VecDestroy(seq_vec,ierr)
  
  else ! uniform grid

    if (myrank == 0) then
      del_x(1:n) = uni_dx*0.1d0
      del_y(1:n) = uni_dy*0.1d0
      del_z(1:n) = uni_dz*0.1d0
    endif

  endif

end subroutine ouComputeGlobalGridSpacing

! ************************************************************************** !
!
! ouPrintConnections
! author: Glenn Hammond
! date: 1/23/01
!
! ************************************************************************** !
subroutine ouPrintConnections(outunit)

  use globals

  implicit none

  integer :: outunit

  integer :: i

  write(outunit,*) 'Local Processor Bounds'
  write(outunit,*) 'Remember that these are based on C numbering.'
  write(outunit,*) 'lnx: ', lnx, '   -> ', lxs, ' -', lxe
  write(outunit,*) 'lny: ', lny, '   -> ', lys, ' -', lye
  write(outunit,*) 'lnz: ', lnz, '   -> ', lzs, ' -', lze

  write(outunit,*)
  write(outunit,*) 'gnx: ', gnx, '   -> ', gxs, ' -', gxe
  write(outunit,*) 'gny: ', gny, '   -> ', gys, ' -', gye
  write(outunit,*) 'gnz: ', gnz, '   -> ', gzs, ' -', gze

  write(outunit,*)
  write(outunit,*) 'number of connections: ', num_connections
  write(outunit,*)

  write(outunit,*)
  write(outunit,*) 'Connections (ghosted indexing)'
  do i=1,size(connection,2)
!    print *, i, connection(upwind,i), connection(downwind,i)
    write(outunit,*) i, connection(upwind,i), connection(downwind,i)
  enddo
  write(outunit,*)

end subroutine ouPrintConnections

! ************************************************************************** !
!
! ouPrintBoundaryConnections
! author: Glenn Hammond
! date: 1/23/01
!
! ************************************************************************** !
subroutine ouPrintBoundaryConnections(outunit)

  use globals

  implicit none

  integer :: outunit

  integer :: i
  type(bcon_type), pointer :: cur_bc

  write(outunit,*)
  write(outunit,*) 'Boundary Connections (local indexing)'
  write(outunit,*)
  cur_bc => bcs
  do 
    if (.not.associated(cur_bc)) exit

    write(outunit,*)
    write(outunit,*) cur_bc%location, '  -    local indexing   global indexing'
    write(outunit,*)
    do i=1,cur_bc%num_connections
      write(outunit,*) i, cur_bc%connection(i), '      ', &
                       node_id_ltog(cur_bc%connection(i))
    enddo
    cur_bc => cur_bc%next
  enddo
  write(outunit,*)

end subroutine ouPrintBoundaryConnections

! ************************************************************************** !
!
! ouPrintVelocities
! author: Glenn Hammond
! date: 1/23/01
!
! ************************************************************************** !
subroutine ouPrintVelocities(Lporosity)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

#ifndef MPIUNI
  integer :: status(MPI_STATUS_SIZE)
#endif

  Vec :: Lporosity

  integer, parameter :: xstart = 1, xend = 2, numx = 3, ystart = 4, &
                        yend = 5, numy = 6, zstart = 7, zend = 8, numz = 9
  integer :: ierr, myid, myid2, i, j, k, & !, x_conn, y_conn, z_conn
             gnode_id, count, &
             lnode_id, global_node_id, size, max_size, global_nxXny, n, &
             cur_proc, ii, jj, kk, cur_connection
  integer, allocatable :: indices(:,:)
  real*8 :: velx, vely, velz, vel_sum
  real*8, allocatable :: global_velx(:), global_vely(:), global_velz(:), &
                         global_sum(:), vel_info(:)
  type(bcon_type), pointer :: cur_bc
  real*8, pointer :: lvec_ptr(:)

  if (.not.allocated(vel_x)) return

  allocate(vel_info(4*lnum_nodes))
  vel_info = 0.d0
  velx = 0.d0
  vely = 0.d0
  velz = 0.d0

  call VecGetArrayF90(Lporosity,lvec_ptr,ierr)
  ou_porosity = lvec_ptr(1)
  call VecRestoreArrayF90(Lporosity,lvec_ptr,ierr)

  count = 0
!  do k=lzs-gzs+1,lzs-gxs+lnz
!    do j=lys-gys+1,lys-gys+lny
!      do i=lxs-gxs+1,lxs-gxs+lnx

  do k=1,gnz
    do j=1,gny
      do i=1,gnx

        if (k > lzs-gzs .and. k <= lzs-gzs+lnz .and. &
            j > lys-gys .and. j <= lys-gys+lny .and. &
            i > lxs-gxs .and. i <= lxs-gxs+lnx) then

          count = count + 1
          vel_sum = 0.d0

          gnode_id = i + (j-1)*gnx + (k-1)*gnxXny

          if (nx > 1) then
            if (i == 1) then
              velx = vel_x(gnode_id)/2.d0
              vel_sum = vel_sum - vel_x(gnode_id)
            else    
              velx = (vel_x(gnode_id) + vel_x(gnode_id-1))/2.d0
              vel_sum = vel_sum - vel_x(gnode_id) + vel_x(gnode_id-1)
            endif
          endif
  
          if (ny > 1) then
            if (j == 1) then
              vely = vel_y(gnode_id)/2.d0
              vel_sum = vel_sum - vel_y(gnode_id)
            else
              vely = (vel_y(gnode_id) + vel_y(gnode_id-gnx))/2
              vel_sum = vel_sum - vel_y(gnode_id) + vel_y(gnode_id-gnx)
            endif
          endif
  
          if (nz > 1) then
            ii = i-(lxs-gxs)
            jj = j-(lys-gys)
            if (k == 1) then
              velz = vel_z(gnode_id)/2.d0
              vel_sum = vel_sum - vel_z(gnode_id)
            else
              velz = (vel_z(gnode_id) + vel_z(gnode_id-gnxXny))/2
              vel_sum = vel_sum - vel_z(gnode_id) + vel_z(gnode_id-gnxXny)
            endif
          endif

          vel_info(count) = vel_sum
          vel_info(lnum_nodes+count) = velx
          vel_info(2*lnum_nodes+count) = vely
          vel_info(3*lnum_nodes+count) = velz
!          if (myrank == 1) print *, count, velx, vely
        endif

      enddo
    enddo
  enddo

  cur_bc => bcs

! Boundary Connections
  do

    if (.not.associated(cur_bc)) exit

    if (cur_bc%location == 'x_upwind') then
      ii = lnum_nodes
      do cur_connection=1,cur_bc%num_connections
        lnode_id = cur_bc%connection(cur_connection)
        velx = cur_bc%vel(cur_connection)/2.d0
        vel_info(ii+lnode_id) = vel_info(ii+lnode_id) + velx
        vel_info(lnode_id) = vel_info(lnode_id) + velx
      enddo
    elseif (cur_bc%location == 'x_downwind') then
      ii = lnum_nodes
      do cur_connection=1,cur_bc%num_connections
        lnode_id = cur_bc%connection(cur_connection)
        velx = cur_bc%vel(cur_connection)/2.d0
        vel_info(ii+lnode_id) = vel_info(ii+lnode_id) - velx
        vel_info(lnode_id) = vel_info(lnode_id) - velx
      enddo
    elseif (cur_bc%location == 'y_upwind  ') then
      ii = 2*lnum_nodes
      do cur_connection=1,cur_bc%num_connections
        lnode_id = cur_bc%connection(cur_connection)
        vely = cur_bc%vel(cur_connection)/2.d0
        vel_info(ii+lnode_id) = vel_info(ii+lnode_id) + vely
        vel_info(lnode_id) = vel_info(lnode_id) + vely
      enddo
    elseif (cur_bc%location == 'y_downwind') then
      ii = 2*lnum_nodes
      do cur_connection=1,cur_bc%num_connections
        lnode_id = cur_bc%connection(cur_connection)
        vely = cur_bc%vel(cur_connection)/2.d0
        vel_info(ii+lnode_id) = vel_info(ii+lnode_id) - vely
        vel_info(lnode_id) = vel_info(lnode_id) - vely
      enddo
    elseif (cur_bc%location == 'z_upwind  ') then
      ii = 3*lnum_nodes
      do cur_connection=1,cur_bc%num_connections
        lnode_id = cur_bc%connection(cur_connection)
        velz = cur_bc%vel(cur_connection)/2.d0
        vel_info(ii+lnode_id) = vel_info(ii+lnode_id) + velz
        vel_info(lnode_id) = vel_info(lnode_id) + velz
      enddo
    else
      ii = 3*lnum_nodes
      do cur_connection=1,cur_bc%num_connections
        lnode_id = cur_bc%connection(cur_connection)
        velz = cur_bc%vel(cur_connection)/2.d0
        vel_info(ii+lnode_id) = vel_info(ii+lnode_id) - velz
        vel_info(lnode_id) = vel_info(lnode_id) - velz
      enddo
    endif

    cur_bc => cur_bc%next
  
  enddo


#if 0
  do k=1,gnz
    do j=1,gny
      do i=1,gnx

        if (k > lzs-gzs .and. k <= lzs-gzs+lnz .and. &
            j > lys-gys .and. j <= lys-gys+lny .and. &
            i > lxs-gxs .and. i <= lxs-gxs+lnx) then

          count = count + 1
          vel_sum = 0.d0

          if (nx > 1) then
            jj = j-(lys-gys)
            kk = k-(lzs-gzs)
            if (i == 1) then
              x_conn = 1 + (jj-1)*(gnx-1) + (kk-1)*(gnx-1)*lny
              velx = vel_x(x_conn)/2.d0
              vel_sum = vel_sum - vel_x(x_conn)
            elseif (i == gnx) then
              x_conn = jj*(gnx-1) + (kk-1)*(gnx-1)*lny
              velx = vel_x(x_conn)/2.d0
              vel_sum = vel_sum + vel_x(x_conn)
            else    
              x_conn = (i-1) + (jj-1)*(gnx-1) + (kk-1)*(gnx-1)*lny
              velx = (vel_x(x_conn) + vel_x(x_conn+1))/2.d0
              vel_sum = vel_sum - vel_x(x_conn) + vel_x(x_conn+1)
            endif
          endif
  
          if (ny > 1) then
            ii = i-(lxs-gxs)
            kk = k-(lzs-gzs)
            if (j == 1) then
              y_conn = max_x_connection + 1 + (ii-1)*(gny-1) + (kk-1)*lnx*(gny-1)
              vely = vel_y(y_conn)/2.d0
              vel_sum = vel_sum - vel_y(y_conn)
            elseif (j == gny) then
              y_conn = max_x_connection + ii*(gny-1) + (kk-1)*lnx*(gny-1)
              vely = vel_y(y_conn)/2.d0
              vel_sum = vel_sum + vel_y(y_conn)
            else
              y_conn = max_x_connection + (j-1) + (ii-1)*(gny-1) + (kk-1)*lnx*(gny-1)
              vely = (vel_y(y_conn) + vel_y(y_conn+1))/2
              vel_sum = vel_sum - vel_y(y_conn) + vel_y(y_conn+1)
            endif
          endif
  
          if (nz > 1) then
            ii = i-(lxs-gxs)
            jj = j-(lys-gys)
            if (k == 1) then
              z_conn = max_y_connection + 1 + (ii-1)*(gnz-1) + (jj-1)*lnx*(gnz-1)
              velz = vel_z(z_conn)/2.d0
              vel_sum = vel_sum - vel_z(z_conn)
            elseif (k == gnz) then
              z_conn = max_y_connection + ii*(gnz-1) + (jj-1)*lnx*(gnz-1)
              velz = vel_z(z_conn)/2.d0
              vel_sum = vel_sum + vel_z(z_conn)
            else
              z_conn = max_y_connection + (k-1) + (ii-1)*(gnz-1) + (jj-1)*lnx*(gnz-1)
              velz = (vel_z(z_conn) + vel_z(z_conn+1))/2
              vel_sum = vel_sum - vel_z(z_conn) + vel_z(z_conn+1)
            endif
          endif

          vel_info(count) = vel_sum
          vel_info(lnum_nodes+count) = velx
          vel_info(2*lnum_nodes+count) = vely
          vel_info(3*lnum_nodes+count) = velz
!          if (myrank == 1) print *, count, velx, vely
        endif

      enddo
    enddo
  enddo
#endif

  if (myrank /= 0) then

!    do i=1,lnum_nodes
!      print *, i, vel_info(i), vel_info(i+lnum_nodes), vel_info(i+2*lnum_nodes)
!    enddo

    allocate(indices(numz,1))
    indices(xstart,1) = lxs
    indices(xend,1) = lxe
    indices(numx,1) = lnx
    indices(ystart,1) = lys
    indices(yend,1) = lye
    indices(numy,1) = lny
    indices(zstart,1) = lzs
    indices(zend,1) = lze 
    indices(numz,1) = lnz 

#ifndef MPIUNI 
!    print *, 'Beginning Send', myrank
    call MPI_Send(indices,numz,MPI_INTEGER, &
                  0,myrank,PETSC_COMM_WORLD,ierr)
!    print *, 'Ending Send', myrank

    call MPI_Send(vel_info,4*lnum_nodes,MPI_DOUBLE_PRECISION, &
                  0,myrank,PETSC_COMM_WORLD,ierr)
#endif

  else

    allocate(indices(numz,commsize))

    n = nx*ny*nz
    allocate(global_velx(n),global_vely(n),global_velz(n),global_sum(n))
    global_nxXny = nx*ny

    global_velx = 0.d0
    global_vely = 0.d0
    global_velz = 0.d0
    global_velx = 0.d0

    count = 0
    do k=lzs,lze-1
      do j=lys,lye-1
        do i=lxs+1,lxe
          count = count + 1
          global_node_id = i + j*nx + k*global_nxXny
          global_sum(global_node_id) = vel_info(count)
          global_velx(global_node_id) = vel_info(lnum_nodes+count)
          global_vely(global_node_id) = vel_info(2*lnum_nodes+count)
          global_velz(global_node_id) = vel_info(3*lnum_nodes+count)
        enddo
      enddo
    enddo

    indices(xstart,1) = lxs
    indices(xend,1) = lxe
    indices(numx,1) = lnx
    indices(ystart,1) = lys
    indices(yend,1) = lye
    indices(numy,1) = lny
    indices(zstart,1) = lzs
    indices(zend,1) = lze 
    indices(numz,1) = lnz 

    max_size = (lxe-lxs)*(lye-lys)*(lze-lzs)

#ifndef MPIUNI 
    do cur_proc=2,commsize

!      print *, 'Beginning Recv', cur_proc
      call MPI_Recv(indices(:,cur_proc),numz,MPI_INTEGER, &
                    cur_proc-1,MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
!      print *, 'Ending Recv', cur_proc
      size = indices(numx,cur_proc)*indices(numy,cur_proc)* &
             indices(numz,cur_proc)
      if (max_size < size) max_size = size

    enddo
#endif

    deallocate(vel_info)
    allocate(vel_info(4*max_size))
    vel_info = 0.d0

#ifndef MPIUNI 
    do cur_proc=2,commsize

      size = indices(numx,cur_proc)*indices(numy,cur_proc)* &
             indices(numz,cur_proc)
!      print *, 'Beginning Recv', cur_proc
      call MPI_Recv(vel_info,4*size,MPI_DOUBLE_PRECISION, &
                    cur_proc-1,MPI_ANY_TAG,PETSC_COMM_WORLD,status,ierr)
!      print *, 'Ending Recv', cur_proc

      count = 0
      do k=indices(zstart,cur_proc),indices(zend,cur_proc)-1
        do j=indices(ystart,cur_proc),indices(yend,cur_proc)-1
          do i=indices(xstart,cur_proc)+1,indices(xend,cur_proc)
            count = count + 1
            global_node_id = i + j*nx + k*global_nxXny
            global_sum(global_node_id) = vel_info(count)
            global_velx(global_node_id) = vel_info(size+count)
            global_vely(global_node_id) = vel_info(2*size+count)
            global_velz(global_node_id) = vel_info(3*size+count)
!            print *, count, global_node_id, vel_info(size+count), vel_info(2*size+count)
          enddo
        enddo
      enddo
    enddo
#endif

    myid=67
    myid2=68
    open(unit=myid,file='flow_vec.dat',iostat=ierr)
    open(unit=myid2,file='flow_mass.dat',iostat=ierr)

    ! sec -> yr
    ! dm -> m
    vel_sum = 3600.d0*24.d0*365.d0/10.d0

    global_velx = global_velx*vel_sum
    global_vely = global_vely*vel_sum
    global_velz = global_velz*vel_sum

    write(myid,10)
    if (nz == 1) then
      write(myid,15)
    else
      write(myid,16)
    endif
    write(myid,30)
    write(myid,35) n
    write(myid,40) n
    write(myid,50)
    write(myid,60)

    write(myid2,10)
    if (nz == 1) then
      write(myid2,15)
    else
      write(myid2,16)
    endif
    write(myid2,31)
    write(myid2,35) n
    write(myid2,40) n
    write(myid2,51)
    write(myid2,60)

    do i=1,n
      ! dm -> m
      write(myid,65) global_velx(i)/ou_porosity, &
                     global_vely(i)/ou_porosity, &
                     global_velz(i)/ou_porosity
      write(myid2,66) global_sum(i)
    enddo

    write(myid,95)
    write(myid2,95)

    close(myid)
    close(myid2)

  endif

  if (allocated(indices)) deallocate(indices)
  if (allocated(vel_info)) deallocate(vel_info)
  if (allocated(global_velx)) deallocate(global_velx,global_vely,global_velz, &
                                         global_sum)

10 format('DATASET')
15 format('OBJTYPE "grid2d"')
16 format('OBJTYPE "grid3d"')
30 format('BEGVEC')
31 format('BEGSCL')
35 format('ND',i10)
40 format('NC',i10)
50 format('NAME "velocity"')
51 format('NAME "mass balance"')
60 format('TS 0 0')
65 format(es12.4,es12.4,es12.4)
66 format(es12.4)
95 format('ENDDS')

end subroutine ouPrintVelocities

! ************************************************************************** !
!
! ouVecMPIToSeqZero: converts a MPI global vec to a sequential on p0
! author: Glenn Hammond
! date: 11/17/03
!
! ************************************************************************** !
subroutine ouVecMPIToSeqZero(global,seq)

  use globals

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  Vec :: global, seq

  integer :: i, j, k, iproc, count, i_start, j_start, ierr
  integer :: status(MPI_STATUS_SIZE)
  real*8, pointer :: gvec_ptr(:), svec_ptr(:)

  if (.not.allocated(indices).and.myrank == 0) then
    print *, 'ouInitOutput needs to be called before ouVecMPIToSeqZero'
    stop
  endif

  call VecGetArrayF90(global,gvec_ptr,ierr)
  if (myrank == 0) then
    call VecGetArrayF90(seq,svec_ptr,ierr)

    ! reorder p0 data
    count = 1
    iproc = 1
    do k=indices(Z_START,iproc),indices(Z_END,iproc)
      j_start = (k-1)*nx*ny
      do j=indices(Y_START,iproc),indices(Y_END,iproc)
        i_start = j_start + (j-1)*nx
        do i=indices(X_START,iproc),indices(X_END,iproc)
          svec_ptr(i_start+i) = gvec_ptr(count)
!          print *, i_start+i, svec_ptr(i_start+i), count, gvec_ptr(count)
          count = count + 1
        enddo
      enddo
    enddo

    allocate(output_buffer(maxval(indices(N_SPATIAL_DOF,1:commsize))))

    do iproc=2,commsize

      call MPI_Recv(output_buffer,indices(N_SPATIAL_DOF,iproc), &
                    MPI_DOUBLE_PRECISION,iproc-1,MPI_ANY_TAG, &
                    PETSC_COMM_WORLD,status,ierr)

      ! reorder
      count = 1
      do k=indices(Z_START,iproc),indices(Z_END,iproc)
        j_start = (k-1)*nx*ny
        do j=indices(Y_START,iproc),indices(Y_END,iproc)
          i_start = j_start + (j-1)*nx
          do i=indices(X_START,iproc),indices(X_END,iproc)
            svec_ptr(i_start+i) = output_buffer(count)
!            print *, i_start+i, svec_ptr(i_start+i), count, output_buffer(count)
            count = count + 1
          enddo
        enddo
      enddo
    enddo

    deallocate(output_buffer)

    call VecRestoreArrayF90(seq,svec_ptr,ierr)
  else

    call MPI_Send(gvec_ptr,lnum_nodes,MPI_DOUBLE_PRECISION, &
                  0,myrank,PETSC_COMM_WORLD,ierr)

  endif
  call VecRestoreArrayF90(global,gvec_ptr,ierr)

end subroutine ouVecMPIToSeqZero

! ************************************************************************** !
!
! ouPrintArray
! author: Glenn Hammond
! date: 1/24/01
!
! ************************************************************************** !
subroutine ouPrintArray(filename,vector,length,label)

  use globals

  implicit none

  character(len=*) :: filename
  character(len=*) :: label
  integer :: length
  real*8 :: vector(length)

  integer :: outunit
  integer :: i

  outunit = 99
  open(unit=outunit,file=filename)

  length = size(vector)
  write(outunit,*)
  write(outunit,*) label
  write(outunit,*)
  do i=1,length
    write(outunit,*) i, vector(i)
  enddo
  write(outunit,*)
  close(outunit)

end subroutine ouPrintArray

! ************************************************************************** !
!
! ouPrintDispersionCoefficients
! author: Glenn Hammond
! date: 9/5/01
!
! ************************************************************************** !
subroutine ouPrintDispersionCoefficients(filename)

  use globals

  implicit none


  character(len=*) :: filename

  integer :: outunit
  integer :: i


  outunit = 99
  open(unit=outunit,file=filename)

  if (uniform_flow) then
    write(outunit,*) 'vx vy vz disp_xx disp_yy ', &
                  'disp_zz disp_xy disp_yx disp_xz disp_zx disp_yz disp_zy'
    write(outunit,100) uni_vel_x, uni_vel_y, uni_vel_z, &
                       uni_disp_xx, uni_disp_yy, uni_disp_zz, &
                       uni_disp_xy, uni_disp_xy, &
                       uni_disp_xz, uni_disp_xz, &
                       uni_disp_yz, uni_disp_yz
  else
    if (gnz > 1) then
      write(outunit,*) 'cur_node vx vy vz disp_xx disp_yy ', &
                  'disp_zz disp_xy disp_yx disp_xz disp_zx disp_yz disp_zy'
      do i=1,gnum_nodes
        write(outunit,200) i, vel_x(i), vel_y(i), vel_z(i), &
                           disp_xx(i), disp_yy(i), disp_zz(i), &
                           disp_xy(i), disp_yx(i), &
                           disp_xz(i), disp_zx(i), &
                           disp_yz(i), disp_zy(i)
      enddo
    else
      write(outunit,*) 'cur_node vx vy disp_xx disp_yy ', &
                  'disp_xy disp_yx'
      do i=1,gnum_nodes
        write(outunit,200) i, vel_x(i), vel_y(i), &
                           disp_xx(i), disp_yy(i), &
                           disp_xy(i), disp_yx(i)
      enddo
    endif
  endif

  close(outunit)

100 format(20(es17.8))
200 format(i5,20(es17.8))

end subroutine ouPrintDispersionCoefficients

! ************************************************************************** !
!
! ouCalculateXYZ: Returns the dimensions of the output (e.g. X, XY, XYZ)
! author: Glenn Hammond
! date: 12/02/03
!
! ************************************************************************** !
subroutine ouCalculateXYZ(xyz)

  use globals, only: nx, ny, nz

  implicit none

  integer :: xyz

  xyz = 0
  if (nx > 1 .and. ny > 1 .and. nz > 1) then
    xyz = XYZ_DIM
  elseif (nx > 1 .and. ny > 1) then
    xyz = XY_DIM
  elseif (ny > 1 .and. nz > 1) then
    xyz = YZ_DIM
  elseif (nx > 1 .and. nz > 1) then
    xyz = XZ_DIM
  elseif (nx > 1) then
    xyz = X_DIM
  elseif (ny > 1) then
    xyz = Y_DIM
  elseif (nz > 1) then
    xyz = Z_DIM
  endif

end subroutine ouCalculateXYZ

end module output
