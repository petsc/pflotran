#!/usr/bin/env python2.7
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--CFLAGS=-g -O0',
    '--CXXFLAGS=-g -O0',
    '--FFLAGS=-g -O0',
    '--download-cmake=yes',
    '--download-fblaslapack=yes',
    '--download-metis=yes',
    '--download-parmetis=yes',
    '--with-hdf5-dir=/data/software/sl-6.x86_64/modules/gcc/4.7.3/hdf5/1.8.11-gcc-p',
    '--with-mpi-dir=/data/software/sl-6.x86_64/modules/gcc/4.7.3/openmpi/1.6.5-gcc',
    '--with-clanguage=c',
    '--with-fortran=1',
    '--with-debugging=1',
    '--with-shared-libraries=0',
    'PETSC_ARCH=aurelien-gnu',
  ]
  configure.petsc_configure(configure_options)
