# Copyright 2005, The Regents of the University 
# of California. This software was produced under
# a U.S. Government contract (W-7405-ENG-36) 
# by Los Alamos National Laboratory, which is
# operated by the University of California for the
# U.S. Department of Energy. The U.S.
# Government is licensed to use, reproduce, and
# distribute this software. Permission is granted
# to the public to copy and use this software
# without charge, provided that this Notice and
# any statement of authorship are reproduced on
# all copies. Neither the Government nor the
# University makes any warranty, express or
# implied, or assumes any liability or
# responsibility for the use of this software.
#
##
## File:        Makefile
## Package:     SAMRAI applications
## Copyright:   (c) 1997-1999 The Regents of the University of California
## Release:     $Name$
## Revision:    $Revision: 2854 $
## Modified:    $Date: 2006-08-23 15:10:50 -0600 (Wed, 23 Aug 2006) $
## Description: makefile for application to test implicit/nonlinear solver
##
BUILDTYPE :=$(shell echo "`uname -p`-`uname -s`")

# for now use the SAMRAI definition for F77 which should default to pgf95 or fortran
F90=${F77}
FLINKER=${F77}

ALL: SAMRAIDriver

TOPDIR  = ../../../../..

OBJECT      = ${SAMRAI}

include $(OBJECT)/config/Makefile.config

M4DIRS        = -DSAMRAI_FORTDIR=$(SAMRAI)/\`include\'
PDIM = 3

ifeq "$(BUILDTYPE)" "x86_64-Linux"
   CXXFLAGS_EXTRA += -DNDIM=$(PDIM)
   CPPFLAGS_EXTRA += -DNDIM=$(PDIM)
endif
ifeq "$(BUILDTYPE)" "i386-Darwin"
   CXXFLAGS_EXTRA += -DNDIM=$(PDIM)
   CPPFLAGS_EXTRA += -DNDIM=$(PDIM)
endif

CPPFLAGS_EXTRA += -I.
CPPFLAGS_EXTRA += -I${AMRUTILITIES_HOME}/include
CXXFLAGS_EXTRA += -I${AMRUTILITIES_HOME}/include
CPPFLAGS_EXTRA += -I${PRECONDITIONER_HOME}/include
CXXFLAGS_EXTRA += -I${PRECONDITIONER_HOME}/include

#CPPFLAGS_EXTRA += -I${TIMEINTEGRATOR_HOME}/include
#CXXFLAGS_EXTRA += -I${TIMEINTEGRATOR_HOME}/include
#LDFLAGS_EXTRA += -L${AMRUTILITIES_HOME}/lib   
#LDLIBS_EXTRA += ${AMRUTILITIES_HOME}/lib/libAMRUtils3d.a
#LDFLAGS_EXTRA += -L${TIMEINTEGRATOR_HOME}/lib   
#LDLIBS_EXTRA += ${TIMEINTEGRATOR_HOME}/lib/libtimeintegration3d.a

SRCDIR = .
VPATH = .

PFLOW_DIR = ../..

CSRC=${wildcard ${SRCDIR}/*.c}
COBJS=${CSRC:.c=.o}

CXXSRC=${wildcard ${SRCDIR}/*.C}
CXXOBJS=${CXXSRC:.C=.o}

F90SRC  = ${wildcard ./*.F90}
F90OBJS = ${F90SRC:.F90=.o}
#F90OBJS = ../pflow_gridtype.o grid_create.o

FSRC  = ${wildcard ./*.f}
#FOBJS = ${FSRC:.f=.o}

FFLAGS		 = -g

FPPFLAGS         = -g

PFLOW_LD_FLAGS = -L${PFLOW_DIR} -lpflotranamr

LDFLAGS += ${PFLOW_LD_FLAGS}

MYFLAGS = ${PETSCFLAGS} -I${PFLOW_DIR} -g -p ${PFLOW_DIR}

ifdef petscrel
 include ${PETSC_DIR}/bmake/common/base
else
 include ${PETSC_DIR}/conf/base
endif

.F90.mod :
	${FLINKER} $< -c ${PETSC_INCLUDE} 

.F90.o :
	${F90} ${FOPTFLAGS} $< ${MYFLAGS} -c ${PETSC_INCLUDE}

.f.o :
	${F90} ${FOPTFLAGS} $< ${MYFLAGS} -c ${PETSC_INCLUDE}


SAMRAIDriver: SAMRAIDriver.o $(COBJS) $(CXXOBJS) ${F90OBJS}
	$(CXX) -g $(CXXOBJS) $(COBJS) ${F90OBJS} ${FOBJS} ${PETSC_FORTRAN_LIB} $(LDFLAGS) $(LIBSAMRAI3D) $(LIBSAMRAI) -lhdf5_fortran  $(LDLIBS) -lU77 -o $@

clean:
	$(RM) SAMRAIDriver *.f core *~ *.mod
	$(RM) *.o *.ii *.int.c
	$(RM) -r ti_files ii_files

veryclean:
	$(MAKE) clean
	$(RM) -fr .deps/*.d
