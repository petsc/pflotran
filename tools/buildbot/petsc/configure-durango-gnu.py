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
    '--with-hdf5-dir=/opt/local',
    '--with-blas-lapack-lib=/System/Library/Frameworks/Accelerate.framework/Versions/Current/Accelerate',
    '--download-parmetis=yes',
    '--download-metis=yes',
    '--with-cc=/opt/local/bin/mpicc-mp',
    '--with-cxx=/opt/local/bin/mpicxx-mp',
    '--with-fc=/opt/local/bin/mpif90-mp',
    '--with-mpiexec=/opt/local/bin/mpiexec-mp',
    '--with-shared-libraries=0',
    'PETSC_ARCH=durango-gnu',
  ]
  configure.petsc_configure(configure_options)
