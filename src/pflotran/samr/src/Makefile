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
BUILDTYPE :=$(shell echo "`uname -m`-`uname -s`")
 .SUFFIXES: .F90
# for now use the SAMRAI definition for F77 which should default to pgf95 or fortran
F90=${F77}
FLINKER=${F77}

ALL: SAMRAIDriver

TOPDIR  = ../../../../..

OBJECT      = ${SAMRAI}

include $(OBJECT)/config/Makefile.config

M4DIRS        = -DSAMRAI_FORTDIR=$(SAMRAI)/\`include\'
PDIM = 3

FORTDIR  = fortran/${PDIM}d
FORTM4SRC  = ${wildcard ${FORTDIR}/*.m4}
FORTOBJS = ${FORTM4SRC:.m4=.o}

ifeq "$(BUILDTYPE)" "x86_64-Linux"
   CXXFLAGS_EXTRA += -DNDIM=$(PDIM)
   CPPFLAGS_EXTRA += -DNDIM=$(PDIM)
endif
ifeq "$(BUILDTYPE)" "i386-Darwin"
   CXXFLAGS_EXTRA += -DNDIM=$(PDIM)
   CPPFLAGS_EXTRA += -DNDIM=$(PDIM)
endif

ifndef jaguar
HDF5_FLAGS = -I$(HDF5_INCLUDE) -I$(HDF5_LIB) -DUSE_HDF5
HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran
else
HDF5_FFLAGS = $(HDF5_FLIB) -DUSE_HDF5
HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz -lsz
endif

CPPFLAGS_EXTRA += ${HDF5_FLAGS}

CPPFLAGS_EXTRA += -I.

CPPFLAGS_EXTRA += -I${AMRUTILITIES_HOME}/include
CXXFLAGS_EXTRA += -I${AMRUTILITIES_HOME}/include

CPPFLAGS_EXTRA += -I${PRECONDITIONER_HOME}/include
CXXFLAGS_EXTRA += -I${PRECONDITIONER_HOME}/include


PFLOTRAN_LDFLAGS = -L${AMRUTILITIES_HOME}/lib   
PFLOTRAN_LDFLAGS += -L${PRECONDITIONER_HOME}/lib   

PFLOTRAN_LDLIBS = ${PRECONDITIONER_HOME}/lib/libpreconditionerbase3d.a
PFLOTRAN_LDLIBS += ${PRECONDITIONER_HOME}/lib/libmlsolvers3d.a
PFLOTRAN_LDLIBS += ${PRECONDITIONER_HOME}/lib/liblinearops3d.a
PFLOTRAN_LDLIBS += ${AMRUTILITIES_HOME}/lib/libAMRUtils3d.a
PFLOTRAN_LDLIBS += -lpflotranamr
PFLOTRAN_LDLIBS += ${HDF5_LIBS}


SRCDIR = .
VPATH = .

vpath %.m4 fortran/${PDIM}d

PFLOTRAN_DIR = ../..

CSRC=${wildcard ${SRCDIR}/*.c}
COBJS=${CSRC:.c=.o}

CXXSRC=${wildcard ${SRCDIR}/*.C}
CXXOBJS=${CXXSRC:.C=.o}

F90SRC  = ${wildcard ./*.F90}
F90OBJS = ${F90SRC:.F90=.o}

ifndef jaguar
FFLAGS		 = -g
else
FFLAGS += ${HDF5_FFLAGS}
endif

PFLOTRAN_LDFLAGS += -L${PFLOTRAN_DIR}

#LDLIBS_EXTRA += ${PFLOTRAN_LDLIBS}
#LDFLAGS_EXTRA += ${PFLOTRAN_LDFLAGS}
ifdef ABSOFT
MYFLAGS = ${PETSCFLAGS} -I${PFLOTRAN_DIR} -g -p ${PFLOTRAN_DIR}
LDLIBS_EXTRA += -lU77
else
MYFLAGS = ${PETSCFLAGS} -I${PFLOTRAN_DIR} -g
endif

ifdef petscrel
 include ${PETSC_DIR}/bmake/common/base
else
# include ${PETSC_DIR}/conf/base
endif

ifndef jaguar
CXX=mpicxx
endif

ifndef jaguar
FPPFLAGS         = -g ${MYFLAGS} -I${PETSC_DIR}/include
else
FPPFLAGS         = -g ${MYFLAGS} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
endif

$(FORTOBJS):$(FORTM4SRC)
	cd $(FORTDIR); ${MAKE}

OBJS =${CXXOBJS} ${COBJS} ${F90OBJS} ${FORTOBJS}

%.o: %.F90
	${F90} -c ${FFLAGS} ${FC_FLAGS} ${FPPFLAGS} -o $@ $<

SAMRAIDriver: SAMRAIDriver.o ${OBJS} ${PFLOTRAN_DIR}/libpflotranamr.a
	$(CXX) -g ${OBJS} ${PFLOTRAN_LDLIBS} ${PFLOTRAN_LDFLAGS} $(LIBSAMRAI3D) $(LIBSAMRAI) ${LDFLAGS} ${LDLIBS} -o $@

clean:
	$(RM) SAMRAIDriver *.f core *~ *.mod
	$(RM) *.o *.ii *.int.c

veryclean:
	$(MAKE) clean
	$(RM) -fr .deps/*.d
