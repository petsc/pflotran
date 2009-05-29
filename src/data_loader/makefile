common_src   = ./

ifndef jaguar
  HDF5_FLAGS = -I${HDF5_INCLUDE} -DUSE_HDF5
  HDF5_LIBS = -L${HDF5_LIB} -lhdf5 -lz
else
# the below is for jaguar
  HDF5_FLAGS = $(HDF5_CLIB) -DUSE_HDF5
  HDF5_LIBS = $(HDF5_CLIB) -lhdf5 -lz
endif

LIBS = ${HDF5_LIBS}

CFLAGS		 =  ${HDF5_FLAGS}
FFLAGS		 = 
CPPFLAGS         = ${HDF5_FLAGS}
FPPFLAGS         = 
CLEANFILES       = dataloader


dataloader_obj = \
	${common_src}FileIO.o \
	${common_src}HDF.o \
	${common_src}DataLoader.o

dataloader: $(dataloader_obj)
	${CC} -o dataloader ${dataloader_obj} ${LIBS}

#Dependencies
FileIO.o : 
HDF.o : 
DataLoader.o : FileIO.o HDF.o 
