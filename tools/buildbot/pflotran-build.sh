#!/usr/bin/env bash

################################################################################
#
# Build script for the pflotran-dev buildbot system.
#
# Author: Ben Andre <bandre@lbl.gov>
#
# All the logic for buildbot is moved out of the buildbot config and
# placed here.
#
# Requirements:
#
#  - Assumes that the script is being called from the root of the
#    pflotran repo. If not, then the pflotran directory must be set on
#    the command line (-p).
#
# Builder requirements:
#
#  - git, mercurial, working mpi compilers
#
#  - a working petsc configure script that will build petsc and any
#    additional TPLs needed for that builder, i.e. hdf5, metis and
#    parmetis for unstructured mesh.
#
#  - NOTE: petsc must be built with '--with-shared-libraries=0'
#    because we check for the presence of libpetsc.a to figure out if
#    a petsc build is acceptable.
#
#
################################################################################

################################################################################
#
# Global variables
#
################################################################################
BUILDER_ID=
PFLOTRAN_DIR=`pwd`
PETSC_DIR=
PETSC_ARCH=
BUILD_STATUS=0
PETSC_REPO=https://bitbucket.org/petsc/petsc.git

################################################################################
#
# work routines
#
################################################################################

function set-builder-info() {
    BUILDER_ID=`hostname -s`
    echo "pflotran builder id : ${BUILDER_ID}"

    _petsc_version_file=${PFLOTRAN_DIR}/tools/buildbot/petsc/petsc-git-version.txt
    if [ ! -f ${_petsc_version_file} ]; then
        echo "ERROR: could not find petsc version file : ${_petsc_version_file}"
        exit 1
    fi
    PETSC_REQUIRED_VERSION=`cat ${_petsc_version_file}`
    echo "PFLOTRAN requires PETSc git reversion ${PETSC_REQUIRED_VERSION}"

    PETSC_DIR=${PFLOTRAN_DIR}/../petsc.git.${PETSC_REQUIRED_VERSION:0:8}
    PETSC_ARCH=${BUILDER_ID}-${COMPILER}
    echo "Requiring petsc env: "
    echo "  PETSC_DIR=${PETSC_DIR}"
    echo "  PETSC_ARCH=${PETSC_ARCH}"
    echo ""
    echo ""
}

function stage-petsc() {
    echo "PETSc stage"

    _lib_petsc=${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a

    if [ -d ${PETSC_DIR} ]; then
        echo "Found existing PETSc directory."
        cd ${PETSC_DIR}
        _petsc_install_version=`git log --pretty="%H" -1 HEAD`
        echo "Found PETSc version: ${_petsc_install_version}"
        if [ ${_petsc_install_version} != ${PETSC_REQUIRED_VERSION} ]; then
            # petsc dir name != petsc version, shouldn't happen...
            echo "Rebuilding PETSc version: ${PETSC_REQUIRED_VERSION}"
            petsc-build ${PETSC_REQUIRED_VERSION}
        elif [ ! -f ${_lib_petsc} ]; then
            echo "PETSc : could not find libpetsc.a for this PETSC_ARCH. Rebuilding."
            echo "    ${_lib_petsc}"
            petsc-build ${PETSC_REQUIRED_VERSION}
        fi
    else
        echo "Building PETSc from scratch"
        git clone ${PETSC_REPO} ${PETSC_DIR}
        petsc-build ${PETSC_REQUIRED_VERSION}
    fi

    if [ ! -f ${_lib_petsc} ]; then
        echo "PETSc : build failed? : libpetsc.a missing after attempted build."
    else
        echo "PETSc appears to be installed at the correct version"
        if [ -z "${BUILD_STATUS}" ]; then
            BUILD_STATUS=0
        fi
    fi
}

function petsc-build() {
    cd ${PETSC_DIR}
    git checkout $1

    _petsc_config_file="${PFLOTRAN_DIR}/tools/buildbot/petsc/configure-${PETSC_ARCH}.py"
    if [ ! -f ${_petsc_config_file} ]; then
        echo "ERROR: pflotran repository does not contain a petsc configure script for this builder. Expected: '${_petsc_config_file}'"
        exit 1
    else
        echo "Linking PETSc config file from pflotran repo: ${_petsc_config_file}"
    fi
    ln -s  ${_petsc_config_file} ${PETSC_DIR}/configure-${PETSC_ARCH}.py
    echo "Configuring PETSc..."
    python configure-${PETSC_ARCH}.py
    echo "Building PETSc..."
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all
    echo "Testing PETSc..."
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} test
    BUILD_STATUS=$?
}

function stage-pflotran-build() {
    echo "Building PFLOTRAN with :"
    echo "  PETSC_DIR=${PETSC_DIR}"
    echo "  PETSC_ARCH=${PETSC_ARCH}"
    export PETSC_DIR PETSC_ARCH
    _pflotran_flags=
    _flags_file=${PFLOTRAN_DIR}/tools/buildbot/build-flags/${BUILD_FLAGS}.txt
    if [ -f ${_flags_file} ]; then
        _pflotran_flags=`cat ${_flags_file}`
        echo "  pflotran build flags=${_pflotran_flags}"
    else
        echo "Could not find build flags file: ${_flags_file}. Building with 'make pflotran'."
    fi
    
    cd ${PFLOTRAN_DIR}/src/pflotran
    make ${_pflotran_flags} clean
    make ${_pflotran_flags} pflotran
    BUILD_STATUS=$?
}

function stage-pflotran-test() {
    echo "Running PFLOTRAN regression and unit tests with:"
    echo "  PETSC_DIR=${PETSC_DIR}"
    echo "  PETSC_ARCH=${PETSC_ARCH}"
    export PETSC_DIR PETSC_ARCH
    _test_dir=${PFLOTRAN_DIR}/regression_tests
    _pflotran_flags=
    _flags_file=${PFLOTRAN_DIR}/tools/buildbot/build-flags/${BUILD_FLAGS}.txt
    if [ -f ${_flags_file} ]; then
        _pflotran_flags=`cat ${_flags_file}`
        echo "  pflotran build flags=${_pflotran_flags}"
    else
        echo "Could not find build flags file: ${_flags_file}. Testing with 'make pflotran'."
    fi

    cd ${PFLOTRAN_DIR}/src/pflotran
    make ${_pflotran_flags} clean-tests &> /dev/null
    make ${_pflotran_flags} test
    BUILD_STATUS=$?
    cat ${_test_dir}/*.testlog
}


################################################################################
#
# main program
#
################################################################################
function usage() {
     echo "
Usage: $0 [options]
    -b BUILD_FLAGS    group of build flags to use.
    -c COMPILER       compiler name: gnu, pgi, intel
    -h                print this help message
    -p PFLOTRAN_DIR   root directory for the build (default: '.')
    -s BUILD_STAGE    build stage must be one of:
                        all petsc pflotran-build pflotran-test

Notes:

  - The build flags group must have a corresponding file in the
    tools/buildbot/build-flags/ directory.

"
}

# setup based on commandline args
BUILD_FLAGS="__NONE__"
BUILD_STAGE=
COMPILER=
while getopts "b:c:hp:s:" FLAG
do
  case ${FLAG} in
    b) BUILD_FLAGS=${OPTARG};;
    c) COMPILER=${OPTARG};;
    h) usage; exit 0;;
    p) PFLOTRAN_DIR=${OPTARG};;
    s) BUILD_STAGE=${OPTARG};;
  esac
done

# verify all required info is set
if [ -z "${BUILD_STAGE}" ]; then
    echo "ERROR: The build stage name must be provided on the command line."
    exit 1
fi

if [ -z "${COMPILER}" ]; then
    echo "ERROR: The compiler name must be provided on the command line."
    exit 1
fi

set-builder-info

echo "PFLOTRAN_DIR: ${PFLOTRAN_DIR}"

case ${BUILD_STAGE} in
    all) stage-petsc; stage-pflotran-build; stage-pflotran-test;;
    petsc) stage-petsc;;
    pflotran-build) stage-pflotran-build;;
    pflotran-test) stage-pflotran-test;;
    *) echo "ERROR: The requested build stage '${BUILD_STAGE}' is invalid."; exit 1;;
esac

exit ${BUILD_STATUS}


