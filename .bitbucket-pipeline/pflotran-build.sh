#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=linux-gnu

cd src/pflotran; 
make pflotran; 

