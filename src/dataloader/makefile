common_src   = ./

HDF5_FLAGS = -I/ccs/home/ghammond/software/lib/hdf5-1.6.5/include
HDF5_LIBS = -L/ccs/home/ghammond/software/lib/hdf5-1.6.5/lib -lhdf5

LIBS = ${HDF5_LIBS}

CC = g++

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
	$(CC) -o dataloader ${dataloader_obj} ${LIBS}

#Dependencies
FileIO.o : 
HDF.o : 
DataLoader.o : FileIO.o HDF.o 
