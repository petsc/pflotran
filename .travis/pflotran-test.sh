#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=petsc-arch

cd src/pflotran; 

# Run unit tests
make utest

# Run regression tests
cd ../../regression_tests

make test

cat *.testlog

