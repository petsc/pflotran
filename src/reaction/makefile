common_src   = ./

CXX		= g++
CXXFLAGS	= -g
CPPFLAGS	= 
CLEANFILES	= reaction


reaction_obj = ${common_src}Block.o \
	${common_src}FileIO.o \
	${common_src}LU.o \
	${common_src}Speciation.o \
	${common_src}Reaction.o

reaction: $(reaction_obj)
	${CXX} -o reaction ${reaction_obj}

%.o : %.cpp
	${CXX} ${CXXFLAGS} $< -c

clean:
	rm -f reaction ${reaction_obj}

#Dependencies
Reaction.o : Speciation.o
Speciation.o : Block.o LU.o
