common_src   = ./

#hdf5 = 1
ifdef hdf5
  HDF5_FLAGS = -I${HDF5_INCLUDE} -DUSE_HDF5
  HDF5_LIBS = -L${HDF5_LIB} -lhdf5 -lz
endif

LIBS = ${HDF5_LIBS}

CFLAGS		 =  ${HDF5_FLAGS}
FFLAGS		 = 
CPPFLAGS         = ${HDF5_FLAGS}
FPPFLAGS         = 
CLEANFILES       = mdt


ifdef dev
 include ${PETSC_DIR}/conf/base
else
 include ${PETSC_DIR}/bmake/common/base
endif

mdt_obj = ${common_src}AsciiGrid.o \
	${common_src}BoundaryConnection.o \
	${common_src}BoundarySet.o \
	${common_src}Condition.o \
	${common_src}Connection.o \
	${common_src}FileIO.o \
	${common_src}Grid.o \
	${common_src}GridCell.o \
	${common_src}GridVertex.o \
	${common_src}Hanford300.o \
	${common_src}HDF.o \
	${common_src}Material.o \
	${common_src}ModelDevelopmentToolkit.o \
	${common_src}Output.o \
	${common_src}Polygon.o \
	${common_src}Source.o \
	${common_src}StructuredGrid.o \
	${common_src}TestCase.o

mdt: $(mdt_obj)
	${CLINKER} -o mdt ${mdt_obj} \
	${PETSC_LIB} ${LIBS}

#Dependencies
AsciiGrid.o : FileIO.o
BoundaryConnection.o : Condition.o
BoundarySet.o : Condition.o Connection.o
Condition.o : FileIO.o
Connection.o : 
FileIO.o : 
GridCell.o :
Grid.o : BoundarySet.o GridCell.o GridVertex.o Source.o StructuredGrid.o
GridVertex.o :
Hanford300.o : AsciiGrid.o BoundarySet.o Connection.o Grid.o Polygon.o
HDF.o :
Material.o :
ModelDevelopmentToolkit.o : Grid.o Hanford300.o Output.o
Output.o : BoundarySet.o Grid.o HDF.o
Polygon.o :
Source.o :
StructuredGrid.o : BoundaryConnection.o GridCell.o GridVertex.o Source.o
TestCase.o : BoundarySet.o Connection.o Grid.o 
