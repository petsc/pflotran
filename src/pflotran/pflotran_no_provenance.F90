
module PFLOTRAN_Provenance_module

!
!-_-! write_warning_comment !-_-!
! 
! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
!
! Purpose: This file is a template / dummy file for
! pflotran_provenance.F90. pflotran-provenance.py automatically
! generantes pflotran_provenance.F90 from pflotran_no_provenance.F90.
! If you can not run pflotran-provenance.py to generate the provenance
! information, then modify your build system to compile and link this
! file instead of pflotran_provenance.F90
!
 
  implicit none

  public

  ! Size of provenance information
  integer(kind=8), parameter :: provenance_max_str_len = 7


  ! PFLOTRAN provenance information
  character(len=*), parameter :: pflotran_compile_date_time = "unknown"
  character(len=*), parameter :: pflotran_changeset = "unknown"
  character(len=*), parameter :: pflotran_status = "unknown"

  integer, parameter :: detail_pflotran_fflags_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_pflotran_fflags(detail_pflotran_fflags_len) = (/ "unknown" /)

  integer, parameter :: detail_pflotran_status_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_pflotran_status(detail_pflotran_status_len) = (/ "unknown" /)

  integer, parameter :: detail_pflotran_parent_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_pflotran_parent(detail_pflotran_parent_len) = (/ "unknown" /)

  ! FIXME(bja, 2013-11-25): break gcc when diffs are present
  integer, parameter :: detail_pflotran_diff_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_pflotran_diff(detail_pflotran_diff_len) = (/ "unknown" /)

  ! PETSc provenance information
  character(len=*), parameter :: petsc_status = "unknown"
  character(len=*), parameter :: petsc_changeset = "unknown"

  integer, parameter :: detail_petsc_status_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_petsc_status(detail_petsc_status_len) = (/ "unknown" /)

  integer, parameter :: detail_petsc_parent_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_petsc_parent(detail_petsc_parent_len) = (/ "unknown" /)

  integer, parameter :: detail_petsc_config_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_petsc_config(detail_petsc_config_len) = (/ "unknown" /)

  public :: PrintProvenanceToScreen, &
       WriteProvenanceToHDF5

contains

subroutine PrintProvenanceToScreen()

  implicit none

  integer i

  write(*, '(''------------------------------ Provenance --------------------------------------'')')

  write(*, '(''pflotran_compile_date_time = '', a)') pflotran_compile_date_time
  write(*, '(''pflotran_changeset = '', a)') pflotran_changeset
  write(*, '(''pflotran_status = '', a)') pflotran_status
  write(*, '(''petsc_changeset = '', a)') petsc_changeset
  write(*, '(''petsc_status = '', a)') petsc_status

!-_-! write_provenance_details !-_-!

  write(*, '(''--------------------------------------------------------------------------------'')')

end subroutine PrintProvenanceToScreen

! ************************************************************************** !

subroutine WriteProvenanceToHDF5(option, output_option, var_list_type)
  !
  ! if hdf5 is compiled in, open the hdf5 output file and write
  ! pflotran and petsc provenance information
  !

  use Option_module, only : option_type
  use Output_Aux_module, only : output_option_type

#include "finclude/petscsysdef.h"

  implicit none

  type(option_type), intent(in) :: option
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: var_list_type  

#ifdef PETSC_HAVE_HDF5
  call ProvenanceToHDF5(option, output_option, var_list_type)
#endif

end subroutine WriteProvenanceToHDF5

! ************************************************************************** !

#ifdef PETSC_HAVE_HDF5
subroutine ProvenanceToHDF5(option, output_option, var_list_type)
  !
  ! open the hdf5 output file and write pflotran and petsc provenance
  ! information
  !

  use Option_module, only : option_type
  use Output_Aux_module, only : output_option_type
  use Output_HDF5_module, only : OutputHDF5OpenFile, OutputHDF5CloseFile

#include "finclude/petscsysdef.h"

  use hdf5

  implicit none

  type(option_type), intent(in) :: option
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: var_list_type  

  character(len=32) :: filename, name
  integer(HID_T) :: file_id, prop_id, provenance_id, string_type
  integer :: hdf5_err
  PetscBool :: first

  call OutputHDF5OpenFile(option, output_option, var_list_type, file_id, first)

  ! NOTE(bja, 2014-03) don't use 'first' because this function is only
  ! called once from main.

  ! create the provenance group
  name = "Provenance"
  call h5gcreate_f(file_id, name, provenance_id, hdf5_err, OBJECT_NAMELEN_DEFAULT_F)

  ! create fixed length string datatype
  call h5tcopy_f(H5T_FORTRAN_S1, string_type, hdf5_err)
  call h5tset_size_f(string_type, provenance_max_str_len, hdf5_err)

  call ProvenanceToHDF5_PFLOTRAN(option, provenance_id, string_type)
  call ProvenanceToHDF5_PETSc(provenance_id, string_type)

  ! close the provenance group
  call h5tclose_f(string_type, hdf5_err)
  call h5gclose_f(provenance_id, hdf5_err)

  call OutputHDF5CloseFile(option, file_id)

end subroutine ProvenanceToHDF5

! ************************************************************************** !

subroutine ProvenanceToHDF5_PFLOTRAN(option, provenance_id, string_type)
  !
  ! write the pflotran provenance data as attributes (small) or
  ! datasets (big details)
  !

  use Option_module, only : option_type
  use Output_HDF5_module, only : OutputHDF5DatasetStringArray, OutputHDF5AttributeStringArray

  use hdf5

  implicit none

  type(option_type), intent(in) :: option
  integer(HID_T), intent(in) :: provenance_id
  integer(HID_T), intent(in) :: string_type

  character(len=32) :: name
  integer(HID_T) :: pflotran_id
  integer :: hdf5_err

  ! Create the pflotran group under provenance
  name = "PFLOTRAN"
  call h5gcreate_f(provenance_id, name, pflotran_id, hdf5_err, OBJECT_NAMELEN_DEFAULT_F)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, "pflotran_compile_date_time", &
       1, pflotran_compile_date_time)

  call OutputHDF5AttributeStringArray(pflotran_id, string_type, "pflotran_status", &
       1, pflotran_status)

  call OutputHDF5AttributeStringArray(pflotran_id, string_type, "pflotran_changeset", &
       1, pflotran_changeset)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, "detail_pflotran_fflags", &
       detail_pflotran_fflags_len, detail_pflotran_fflags)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, "detail_pflotran_status", &
       detail_pflotran_status_len, detail_pflotran_status)

  call OutputHDF5DatasetStringArray(pflotran_id, string_type, "detail_pflotran_parent", &
       detail_pflotran_parent_len, detail_pflotran_parent)

  ! FIXME(bja, 2013-11-25): break gcc when diffs are present  
  call OutputHDF5DatasetStringArray(pflotran_id, string_type, "detail_pflotran_diff", &
       detail_pflotran_diff_len, detail_pflotran_diff)

  call ProvenanceToHDF5_input(option, pflotran_id)

  ! close pflotran group
  call h5gclose_f(pflotran_id, hdf5_err)

end subroutine ProvenanceToHDF5_PFLOTRAN

! ************************************************************************** !

subroutine ProvenanceToHDF5_input(option, pflotran_id)
  !
  ! open the pflotran input file, figure out how long it is, read it
  ! into a buffer, then write the buffer as a pflotran provenance
  ! group dataset.
  !
  use hdf5
  use Input_Aux_module, only : input_type, InputCreate, InputDestroy, &
       InputGetLineCount, InputReadToBuffer
  use Option_module, only : option_type
  use Output_HDF5_module, only : OutputHDF5DatasetStringArray
  use PFLOTRAN_Constants_module, only : IN_UNIT, MAXSTRINGLENGTH

  implicit none

  type(option_type), intent(in) :: option
  integer(HID_T), intent(in) :: pflotran_id

  integer(HID_T) :: input_string_type
  type(input_type), pointer :: input
  integer :: i, input_line_count
  character(len=MAXSTRINGLENGTH), allocatable :: input_buffer(:)
  integer :: hdf5_err

  input => InputCreate(IN_UNIT, option%input_filename, option)
  input_line_count = InputGetLineCount(input)
  allocate(input_buffer(input_line_count))
  call InputReadToBuffer(input, input_buffer)
  call h5tcopy_f(H5T_FORTRAN_S1, input_string_type, hdf5_err)
  call h5tset_size_f(input_string_type, int(MAXSTRINGLENGTH, kind=8), hdf5_err)
  call OutputHDF5DatasetStringArray(pflotran_id, input_string_type, "pflotran_input_file", &
       input_line_count, input_buffer)
  call h5tclose_f(input_string_type, hdf5_err)
  deallocate(input_buffer)
  call InputDestroy(input)

end subroutine ProvenanceToHDF5_input

! ************************************************************************** !

subroutine ProvenanceToHDF5_PETSc(provenance_id, string_type)
  !
  ! write the petsc provenance data as attributes (small) or datasets
  ! (big details)
  !

  use hdf5
  use Output_HDF5_module, only : OutputHDF5DatasetStringArray, OutputHDF5AttributeStringArray

  implicit none

  integer(HID_T), intent(in) :: provenance_id
  integer(HID_T), intent(in) :: string_type

  character(len=32) :: name
  integer(HID_T) :: petsc_id
  integer :: hdf5_err

  ! create the petsc group under provenance
  name = "PETSc"
  call h5gcreate_f(provenance_id, name, petsc_id, hdf5_err, OBJECT_NAMELEN_DEFAULT_F)

  call OutputHDF5AttributeStringArray(petsc_id, string_type, "petsc_status", &
       1, petsc_status)

  call OutputHDF5AttributeStringArray(petsc_id, string_type, "petsc_changeset", &
       1, petsc_changeset)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, "detail_petsc_status", &
       detail_petsc_status_len, detail_petsc_status)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, "detail_petsc_parent", &
       detail_petsc_parent_len, detail_petsc_parent)

  call OutputHDF5DatasetStringArray(petsc_id, string_type, "detail_petsc_config", &
       detail_petsc_config_len, detail_petsc_config)

  ! close the petsc group
  call h5gclose_f(petsc_id, hdf5_err)

end subroutine ProvenanceToHDF5_PETSc

! ************************************************************************** !

#endif
! PETSC_HAVE_HDF
end module PFLOTRAN_Provenance_module
