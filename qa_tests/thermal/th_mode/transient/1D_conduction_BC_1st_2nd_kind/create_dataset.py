import sys
from h5py import *
import numpy as np

filename = 'dataset.h5'
h5file = File(filename,mode='w')

# 1d line in x
# Temperature initial condition
L = 100.
h5grp = h5file.create_group('initial')
h5grp.attrs['Dimension'] = np.string_('X')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [1.]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
nx = 101
rarray = np.zeros(nx,'=f8')
for i in range(nx):
  if (0. <= i < (L/10.)):
    rarray[i] = 0.
  if ((L/10.) <= i < (4.*L/10.)):
    rarray[i] = (10./(3.*L))*float(i) - (1./3.)
  if ((4.*L/10.) <= i < (6.*L/10.)):
    rarray[i] = 1.
  if ((6.*L/10.) <= i < (9.*L/10.)):
    rarray[i] = 3. - (10./(3.*L))*float(i)
  if ((9.*L/10.) <= i < L):
    rarray[i] = 0.
h5dset = h5grp.create_dataset('Data', data=rarray)