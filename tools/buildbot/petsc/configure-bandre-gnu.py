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
    '--with-cc=/opt/local/bin/openmpicc',
    '--with-cxx=/opt/local/bin/openmpic++',
    '--with-fc=/opt/local/bin/openmpif90',
    '--with-mpiexec=/opt/local/bin/openmpiexec',
    '--with-shared-libraries=0',
    'PETSC_ARCH=bandre-gnu',
  ]
  configure.petsc_configure(configure_options)
