common_src = ./

ifdef gcc
CXX              = g++
CXXFLAGS         = -O3
CPPFLAGS         =
else
CXX              = icpc
CXXFLAGS         = -fast
CPPFLAGS         =
endif

CLEANFILES       = performance_test

obj = ${common_src}Species.o \
      ${common_src}Reaction.o \
      ${common_src}Beaker.o \
      ${common_src}Driver.o \

performance_test : ${obj}
	${CXX} -o performance_test ${obj}

%.o : %.cpp
	${CXX} ${CXXFLAGS} $< -c

#Dependencies
Driver.o : Beaker.o
Beaker.o : Reaction.o Species.o
Reaction.o : Species.o
