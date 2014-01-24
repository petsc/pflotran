#!/usr/bin/env bash
#
# script to build pfunit
#

PFUNIT_DIR=
COMPILER=
################################################################################
function build-pfunit() {

    pushd ${PFUNIT_DIR}
    make F90=${COMPILER}
    ./tests/tests.x
    popd
}


################################################################################
#
# main program
#
################################################################################
function usage() {
     echo "
Usage: $0 [options]
    -c FORTRAN_COMPILER   path to fortran compiler
    -d PFUNIT_DIR         path to pfunit root directory
    -h                print this help message

Notes:

  * eventually add flags for mpi

"
}

# setup based on commandline args
while getopts "c:d:h" FLAG
do
  case ${FLAG} in
    c) COMPILER=${OPTARG};;
    d) PFUNIT_DIR=${OPTARG};;
    h) usage;;
  esac
done

# verify all required info is set
if [ -z "${PFUNIT_DIR}" ]; then
    echo "ERROR: The pFUnit root directory must be provided on the command line."
    exit 1
fi

if [ -z "${COMPILER}" ]; then
    echo "ERROR: The compiler must be provided on the command line."
    exit 1
fi


echo "PFUNIT_DIR: ${PFUNIT_DIR}"
echo "FC: ${COMPILER}"

build-pfunit

exit ${BUILD_STATUS}
