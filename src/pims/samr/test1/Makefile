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

F90=pgf90
FLINKER=pgf90

ALL: tests

TOPDIR  = ../../../../..

OBJECT      = ${SAMRAI}

include $(OBJECT)/config/Makefile.config

M4DIRS        = -DSAMRAI_FORTDIR=$(SAMRAI)/\`include\'
PDIM = 3

ifeq "$(BUILDTYPE)" "x86_64-Linux"
   CXXFLAGS_EXTRA += -DNDIM=$(PDIM)
   CPPFLAGS_EXTRA += -DNDIM=$(PDIM)
endif

SRCDIR = .
VPATH = .

TESTS := $(patsubst test_%.C,test_%,$(filter test_%, $(wildcard *.C)))
CXXSRC=${wildcard ${SRCDIR}/*.C}
CXXOBJS=${CXXSRC:.C=.o}

F90SRC  = ${wildcard ./*.F90}
F90OBJS = ${F90SRC:.F90=.o}
#F90OBJS = ../pflow_gridtype.o grid_create.o

FSRC  = ${wildcard ./*.f}
#FOBJS = ${FSRC:.f=.o}

FFLAGS		 = -g

FPPFLAGS         = -g


MYFLAGS = ${PETSCFLAGS} -g

include ${PETSC_DIR}/bmake/common/base

.F90.mod :
	${FLINKER} $< -c ${PETSC_INCLUDE} 

.F90.o :
	${F90} ${FOPTFLAGS} $< ${MYFLAGS} -c ${PETSC_INCLUDE}

.f.o :
	${F90} ${FOPTFLAGS} $< ${MYFLAGS} -c ${PETSC_INCLUDE}

CPPFLAGS_EXTRA += -I.

LDFLAGS_EXTRA +=

tests:  $(TESTS)

test_%: test_%.o $(CXXOBJS) ${F90OBJS}
	$(CXX) -g $(CXXOBJS) ${F90OBJS} ${FOBJS} ${MYFLAGS} ${PETSC_FORTRAN_LIB} $(LDFLAGS) $(LIBSAMRAI3D) $(LIBSAMRAI)  $(LDLIBS) -o $@

check: $(TESTS)

clean:
	$(RM) $(TESTS) *.f core *~
	$(RM) *.o *.ii *.int.c
	$(RM) -r ti_files ii_files

veryclean:
	$(MAKE) clean
	$(RM) -fr .deps/*.d

test_3.o: BogusTagAndInitStrategy.o BoundaryConditionStrategy.o