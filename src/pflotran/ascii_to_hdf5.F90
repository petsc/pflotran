program ascii_to_hdf5

  use fileio_module
  use dataset_write_module
  use hdf5

  implicit none

  integer fid, iostatus, ierr
  character(len=512) :: string 
  character(len=4) :: card
  character(len=32) :: filename
  
  integer, allocatable :: cell_id(:)
  integer, allocatable :: material_id(:)
  real*8, allocatable :: x(:)
  real*8, allocatable :: y(:)
  real*8, allocatable :: z(:)
  real*8, allocatable :: volume(:)
  real*8, allocatable :: permx(:)
  real*8, allocatable :: permy(:)
  real*8, allocatable :: permz(:)
  real*8, allocatable :: porosity(:)
  real*8, allocatable :: tortuosity(:)
  integer, allocatable :: connection_id(:)
  integer, allocatable :: id_upwind(:)
  integer, allocatable :: id_downwind(:)
  real*8, allocatable :: distance_upwind(:)
  real*8, allocatable :: distance_downwind(:)
  real*8, allocatable :: area(:)
  real*8, allocatable :: cosB(:)
  
  real*8, pointer :: dummy_real(:)
  integer, pointer :: dummy_int(:)
    
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
  integer :: hdf5_err
  integer :: i
  integer :: num_cells, num_connections, count1
  integer :: one_tenth

  nullify(dummy_real)
  nullify(dummy_int)
  
  print *, 'Enter filename:'
  read *, filename

  fid = 86
  open(fid, file=filename, action="read", status="old", iostat=iostatus)
  if (iostatus /= 0) then
    print *, 'Error opening file: ', trim(filename)
    stop
  endif
  
! open the hdf5 file
  ! initialize fortran interface to hdf5
  call h5open_f(hdf5_err)
  do i=1,len_trim(filename)
    if (filename(i:i) == '.') exit
  enddo
  filename = filename(1:min(i-1,len_trim(filename))) // '.h5'
  print *, 'Opening HDF5 file: ', trim(filename)
  ! create hdf5 file
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
  call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err,H5P_DEFAULT_F, &
                   prop_id)
  call h5pclose_f(prop_id,hdf5_err)

! read the GRID cell data
! GRID information
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadNChars(string,card,4,.true., ierr)
    if (fiStringCompare(card,"GRID",4) .or. ierr /= 0) exit
  enddo

  ! report error if card does not exist
  if (ierr /= 0) then
    print *, 'ERROR: Card (',card, ') not found in file'
    stop
  endif
  
  print *, 'Reading grid cells'
  
  call fiReadInt(string,num_cells,ierr)
  if (ierr /= 0) then
    print *, 'ERROR: The number of grid cells must be listed after ', &
             ' the keyword GRID (e.g. GRID 2000).'
    stop
  endif

  one_tenth = num_cells/10
   
  allocate(cell_id(num_cells))
  allocate(material_id(num_cells))
  allocate(x(num_cells))
  allocate(y(num_cells))
  allocate(z(num_cells))
  allocate(volume(num_cells))
  allocate(permx(num_cells))
  allocate(permy(num_cells))
  allocate(permz(num_cells))
  allocate(porosity(num_cells))
  allocate(tortuosity(num_cells))
   
  count1 = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    count1 = count1 + 1
    call fiReadInt(string,cell_id(count1),ierr) 
    call fiErrorMsg('cell_id',card,ierr)

    call fiReadInt(string,material_id(count1),ierr) 
    call fiErrorMsg('material_id',card,ierr)

    call fiReadDouble(string,x(count1),ierr)
    call fiErrorMsg('xcoord',card,ierr)
   
    call fiReadDouble(string,y(count1),ierr)
    call fiErrorMsg('ycoord',card,ierr)

    call fiReadDouble(string,z(count1),ierr)
    call fiErrorMsg('zcoord',card,ierr)
   
    call fiReadDouble(string,volume(count1),ierr)
    call fiErrorMsg('volume',card,ierr)
   
    call fiReadDouble(string,permx(count1),ierr)
    call fiErrorMsg('permx',card,ierr)

    call fiReadDouble(string,permy(count1),ierr)
    call fiErrorMsg('permy',card,ierr)
   
    call fiReadDouble(string,permz(count1),ierr)
    call fiErrorMsg('permz',card,ierr)

    call fiReadDouble(string,porosity(count1),ierr)
    call fiErrorMsg('porosity',card,ierr)

    call fiReadDouble(string,tortuosity(count1),ierr)
    call fiErrorMsg('tortuosity',card,ierr)

    if (mod(count1,one_tenth) == 0) print *, count1, ' out of ', &
      num_cells, ' cells read.'

  enddo


  print *, count1, ' cells read'
  
  string = "Grid Cells"
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  
  !write the grid cell data to hdf5
  string = "Cell Id"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,cell_id,dummy_real,num_cells,grp_id, &
                        H5T_NATIVE_INTEGER)
  string = "Material Id"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,material_id,dummy_real,num_cells,grp_id, &
                        H5T_NATIVE_INTEGER)
  string = "X-Coordinate"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,x,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "Y-Coordinate"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,y,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "Z-Coordinate"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,z,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "Volume"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,volume,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "X-Permeability"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,permx,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "Y-Permeability"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,permy,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "Z-Permeability"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,permz,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "Porosity"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,porosity,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "Tortuosity"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,tortuosity,num_cells,grp_id, &
                        H5T_NATIVE_DOUBLE)
  call h5gclose_f(grp_id,hdf5_err)
  
  ! deallocate the grid cell data to leave more room for connection info
  deallocate(cell_id)
  deallocate(material_id)
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(volume)
  deallocate(permx)
  deallocate(permy)
  deallocate(permz)
  deallocate(porosity)
  deallocate(tortuosity)
  
  print *, 'Reading grid connections'

! CONNection information
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadNChars(string,card,4,.true., ierr)
    if (fiStringCompare(card,"CONN",4) .or. ierr /= 0) exit
  enddo

  ! report error if card does not exist
  if (ierr /= 0) then
    print *, 'ERROR: Card (',card, ') not found in file'
    stop
  endif
  
  call fiReadInt(string,num_connections,ierr)
  if (ierr /= 0) then
    print *, 'ERROR: The number of connections must be listed after ', &
             ' the keyword CONN (e.g. CONN 6000).'
    stop
  endif

  one_tenth = num_connections/10

  allocate(connection_id(num_connections))
  allocate(id_upwind(num_connections))
  allocate(id_downwind(num_connections))
  allocate(distance_upwind(num_connections))
  allocate(distance_downwind(num_connections))
  allocate(area(num_connections))
  allocate(cosB(num_connections))

  count1 = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(card,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/') exit

    count1 = count1 + 1
    call fiReadInt(string,connection_id(count1),ierr) 
    call fiErrorMsg('connection_id',card,ierr)

    call fiReadInt(string,id_upwind(count1),ierr) 
    call fiErrorMsg('id_upwind',card,ierr)

    call fiReadInt(string,id_downwind(count1),ierr) 
    call fiErrorMsg('id_downwind',card,ierr)

    call fiReadDouble(string,distance_upwind(count1),ierr) 
    call fiErrorMsg('distance_upwind',card,ierr)

    call fiReadDouble(string,distance_downwind(count1),ierr) 
    call fiErrorMsg('distance_downwind',card,ierr)

    call fiReadDouble(string,area(count1),ierr) 
    call fiErrorMsg('area',card,ierr)

    call fiReadDouble(string,cosB(count1),ierr) 
    call fiErrorMsg('cosB',card,ierr)

    if (mod(count1,one_tenth) == 0) print *, count1, ' out of ', &
      num_connections, ' connections read.'

  enddo

  print *, count1, ' connections read'

  close(fid)

  !write the connection data to hdf5
  string = "Connections"
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  string = "Connection Id"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,connection_id,dummy_real,num_connections, &
                        grp_id,H5T_NATIVE_INTEGER)
  string = "Id Upwind"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,id_upwind,dummy_real,num_connections,grp_id, &
                        H5T_NATIVE_INTEGER)
  string = "Id Downwind"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,id_downwind,dummy_real,num_connections, &
                        grp_id,H5T_NATIVE_INTEGER)
  string = "Distance Upwind"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,distance_upwind,num_connections, &
                        grp_id,H5T_NATIVE_DOUBLE)
  string = "Distance Downwind"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,distance_downwind,num_connections, &
                        grp_id,H5T_NATIVE_DOUBLE)
  string = "Area"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,area,num_connections,grp_id, &
                        H5T_NATIVE_DOUBLE)
  string = "CosB"
  print *, 'Writing ', trim(string), ' to HDF5 file.'
  call WriteHDF5DataSet(string,dummy_int,cosB,num_connections,grp_id, &
                        H5T_NATIVE_DOUBLE)
  call h5gclose_f(grp_id,hdf5_err)

  deallocate(connection_id)
  deallocate(id_upwind)
  deallocate(id_downwind)
  deallocate(distance_upwind)
  deallocate(distance_downwind)
  deallocate(area)
  deallocate(cosB)
  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

  print *, 'Done'
 
end program ascii_to_hdf5
