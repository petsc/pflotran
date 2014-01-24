# -*- mode: makefile -*-
#
# makefile controling generation and execution of unit tests with pFUnit
#

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

SRC_DIR = ..
PFUNIT_DIR = ../../../tpls/pfunit

include $(PFUNIT_DIR)/include/extensions.mk

TEST_SRCS = $(wildcard *.pf)
TEST_OBJS = $(TEST_SRCS:.pf=$(OBJ_EXT))

# this list should only include the modules in the test suite!
PFLOTRAN_OBJS = \
	$(SRC_DIR)/eos_water.o \
	$(SRC_DIR)/input_aux.o \
	$(SRC_DIR)/logging.o \
	$(SRC_DIR)/option.o \
	$(SRC_DIR)/string.o \
	$(SRC_DIR)/utility.o

FFLAGS = -I. -I$(SRC_DIR) -I$(PFUNIT_DIR)/source -I$(PFUNIT_DIR)/include

all : test

# uncomment the PRECIOUS rule to keep the intermediate F90 files
#.PRECIOUS : %.F90

testSuites.inc : $(TEST_SRCS)

test : pfunit unittests
	@echo "----------------------------------------------------------------------"
	@echo "Running pflotran unit tests :"
	@./unittests
	@echo "----------------------------------------------------------------------"

unittests : testSuites.inc $(TEST_OBJS)
	$(FC) -o $@ $(PFUNIT_DIR)/include/driver.F90 \
	$(FFLAGS) \
	$(TEST_OBJS) $(PFLOTRAN_OBJS) -L$(PFUNIT_DIR)/source -lpfunit  ${PETSC_LIB} \


%.F90 : %.pf
	$(PFUNIT_DIR)/bin/pFUnitParser.py $< $@

%$(OBJ_EXT) : %.F90
	$(FC) -c $(FC_FLAGS) $(FFLAGS) ${FCPPFLAGS} $<

pfunit : $(PFUNIT_DIR)/source/libpfunit.a

$(PFUNIT_DIR)/source/libpfunit.a :
	./build-pfunit.sh -c $(FC) -d $(PFUNIT_DIR)

# note, this removes intermediate F90 files!
clean-tests :
	-rm unittests *~ *.F90 *.o *.mod
