program test

  include 'mpif.h'

  integer :: ierr
  integer :: rank
  integer :: commsize
  character(len=32) :: word
  character(len=32) :: word2

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,commsize,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)

  write(word,*) rank
  write(word2,*) commsize
  print *, ' ' // trim(adjustl(word)) // ' of ' // trim(adjustl(word2)) // ' processors.'
  
  call mpi_finalize(ierr)

end program test

