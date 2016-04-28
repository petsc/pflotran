common_src   = ./

CFLAGS		 = 
FFLAGS		 = 
CPPFLAGS         = 
FPPFLAGS         = 
CLEANFILES       = mdt


ifdef rel
 include ${PETSC_DIR}/bmake/common/base
else
 include ${PETSC_DIR}/conf/variables
 include ${PETSC_DIR}/conf/rules
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
	${common_src}Hanford300v2.o \
	${common_src}HDF.o \
	${common_src}IFC.o \
	${common_src}IFC_Scope3.o \
	${common_src}IFC_Plot_120x120x20.o \
	${common_src}IFC_Plot_120x120x15.o \
	${common_src}IFC_2D.o \
	${common_src}MAD.o \
	${common_src}MidIFC.o \
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
Hanford300v2.o : AsciiGrid.o BoundarySet.o Connection.o Grid.o Polygon.o
IFC.o : AsciiGrid.o BoundarySet.o Connection.o Grid.o Polygon.o
IFC_Plot_120x120x20.o : AsciiGrid.o BoundarySet.o Connection.o Grid.o Polygon.o
IFC_Plot_120x120x15.o : AsciiGrid.o BoundarySet.o Connection.o Grid.o Polygon.o
IFC_2D.o : AsciiGrid.o BoundarySet.o Connection.o Grid.o Polygon.o
MAD.o : AsciiGrid.o BoundarySet.o Connection.o Grid.o Polygon.o
HDF.o :
Material.o :
ModelDevelopmentToolkit.o : Grid.o Hanford300.o Hanford300v2.o Output.o IFC.o \
                            IFC_Scope3.o IFC_2D.o MAD.o MidIFC.o
Output.o : BoundarySet.o Grid.o HDF.o
Polygon.o :
Source.o :
StructuredGrid.o : BoundaryConnection.o GridCell.o GridVertex.o Source.o
TestCase.o : BoundarySet.o Connection.o Grid.o 
