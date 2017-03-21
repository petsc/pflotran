#!/bin/sh

sudo apt-get update -qq
sudo apt-get install -y cmake gcc libopenmpi-dev openmpi-bin liblapack-dev gfortran mercurial git 

git clone https://bitbucket.org/petsc/petsc.git

#PETSC_GIT_HASH=`cat tools/buildbot/petsc/petsc-git-version.txt`

cd petsc

# pflotran-xsdk follows the tip (master branch) of PETSc
#git checkout ${PETSC_GIT_HASH}

export PETSC_DIR=$PWD
export PETSC_ARCH=linux-gnu
export PETSC_WITH_HDF5=1

#if [ ${PETSC_WITH_HDF5} == 1 ]; then
./configure PETSC_ARCH=linux-gnu --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=1 --download-hdf5 --download-metis --download-parmetis --download-zlib
#elif [ ${PETSC_WITH_HDF5} == 0 ]; then
#  ./configure PETSC_ARCH=linux-gnu --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=1 --download-metis --download-parmetis
#else
#  echo "Unknown value for PETSC_WITH_HDF5 = ${PETSC_WITH_HDF5}"
#  exit
#fi

make

