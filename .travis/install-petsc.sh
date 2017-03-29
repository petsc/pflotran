#!/bin/sh

git clone https://bitbucket.org/petsc/petsc.git

PETSC_GIT_HASH=`cat tools/buildbot/petsc/petsc-git-version.txt`

cd petsc

git checkout ${PETSC_GIT_HASH}

export PETSC_DIR=$PWD
export PETSC_ARCH=petsc-arch
export PETSC_WITH_HDF5=1

./configure PETSC_ARCH=petsc-arch --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=0 --download-hdf5 --download-metis --download-parmetis

make

