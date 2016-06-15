#!/files0/software/python2.7/bin/python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--CFLAGS=-g -O0',
    '--CXXFLAGS=-g -O0',
    '--FFLAGS=-g -O0',
    '--download-fblaslapack=yes',
    '--download-hypre=yes',
    '--download-metis=yes',
    '--download-parmetis=yes',
    '--with-cc=/files0/software/mpich2-1.4.1p1/gcc-4.7.2/bin/mpicc',
    '--with-clanguage=c',
    '--with-cxx=/files0/software/mpich2-1.4.1p1/gcc-4.7.2/bin/mpicxx',
    '--with-debugging=1',
    '--with-fc=/files0/software/mpich2-1.4.1p1/gcc-4.7.2/bin/mpif90',
    '--with-hdf5-dir=/files0/software/hdf5-1.8.9/gcc-4.7.2',
    '--with-hdf5=1',
    '--with-shared-libraries=0',
    '--with-valgrind=1',
    'PETSC_ARCH=fuji-gnu',
  ]
  configure.petsc_configure(configure_options)
