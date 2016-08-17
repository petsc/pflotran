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
    '--download-fblaslapack=yes',
    '--download-metis=yes',
    '--download-parmetis=yes',
    '--download-hdf5=yes',
    '--download-mpich=yes',
    '--download-cmake=yes',
    '--with-cc=/global/software/sl-6.x86_64/modules/langs/gcc/4.8.5/bin/gcc',
    '--with-cxx=/global/software/sl-6.x86_64/modules/langs/gcc/4.8.5/bin/g++',
    '--with-fc=/global/software/sl-6.x86_64/modules/langs/gcc/4.8.5/bin/gfortran',
    '--with-clanguage=c',
    '--with-debugging=1',
    '--with-shared-libraries=0',
    'PETSC_ARCH=bryce-gnu',
  ]
  configure.petsc_configure(configure_options)
