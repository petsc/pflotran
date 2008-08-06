BUILDTYPE :=$(shell echo "`uname -p`-`uname -s`")

SRCDIR  = .
VPATH   = .
TOPDIR=../../../../../..

# for now use the SAMRAI definition for F77 which should default to pgf95 or fortran
F90=${F77}
FLINKER=${F77}
OBJECT      = ${SAMRAI}

include $(OBJECT)/config/Makefile.config

M4DIRS        = -DSAMRAI_FORTDIR=$(SAMRAI)/\`include\'
PDIM = 3

FFLAGS		 = -g

FPPFLAGS         = -g

FORTM4SRC  = ${wildcard ${SRCDIR}/*.m4}
FORTOBJS = ${FORTM4SRC:.m4=.o}

all: $(FORTOBJS)

FORTOBJS: $(FORTM4SRC)

# preserve the intermediate.f files else fortran source cannot be viewed in totalview
.PRECIOUS:%.f

%.f:%.m4
	$(M4) $(M4DIRS) $< > $@

clean:
	rm -f *.o *~
