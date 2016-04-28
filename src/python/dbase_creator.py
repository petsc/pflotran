import sys
from h5py import *
import numpy

filename = '543_dbase.h5'
h5file = File(filename,mode='w')

num_realizations = 10

variables = []
variables.append('ALPHA-sf3')
variables.append('NORTH_MAX_Z')
variables.append('PERMEABILITY_X::SOIL4')
variables.append('POROSITY::SOIL1')
variables.append('RIVER_TRACER2')

values = []
values.append(1.0211e-4)
values.append(60.)
values.append(1.e-9)
values.append(0.25)
values.append(1.e-7)

rarray = numpy.zeros((num_realizations),'=f8')

irealization = 5

for i in range(len(variables)):
  rarray[irealization-1] = values[i]
  h5dset = h5file.create_dataset(variables[i], data=rarray)

h5file.close()

print('done')
