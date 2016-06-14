#!/usr/bin/env python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--CFLAGS=-g -O0',
    '--CXXFLAGS=-g -O0',
    '--FFLAGS=-g -O0',
    '--download-hdf5=1',
    '--with-blas-lapack-lib=/System/Library/Frameworks/Accelerate.framework/Versions/Current/Accelerate',
    '--download-parmetis=yes',
    '--download-metis=yes',
    '--with-cc=/opt/local/bin/mpicc',
    '--with-cxx=/opt/local/bin/mpic++',
    '--with-fc=/opt/local/bin/mpif90',
    '--with-mpiexec=/opt/local/bin/mpiexec',
    '--with-shared-libraries=0',
    'PETSC_ARCH=margaux-gnu',
  ]
  configure.petsc_configure(configure_options)
