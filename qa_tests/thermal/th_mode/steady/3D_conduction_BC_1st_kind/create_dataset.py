import sys
from h5py import *
import numpy as np

filename = 'dataset.h5'
h5file = File(filename,mode='w')

T0 = 1.  # [C]
L = 1.   # [m]

# 2D Surface:
# -----------------------------------------------------------------------------
# Temperature boundary condition x=0 face
h5grp = h5file.create_group('node_centered_surf_west')
h5grp.attrs['Dimension'] = np.string_('YZ')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1,0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.,0.]
# Load the dataset values
ny = 11
nz = 11
rarray = np.zeros((ny,nz),'=f8')
for j in range(ny):
  for k in range(nz):
    y = (float(j)*L)/(ny-1)
    z = (float(k)*L)/(nz-1)
    rarray[j][k] = T0*(0+(y/L)+(z/L))
h5dset = h5grp.create_dataset('Data', data=rarray)

# 2D Surface:
# -----------------------------------------------------------------------------
# Temperature boundary condition x=L face
h5grp = h5file.create_group('node_centered_surf_east')
h5grp.attrs['Dimension'] = np.string_('YZ')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1,0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.,0.]
# Load the dataset values
ny = 11
nz = 11
rarray = np.zeros((ny,nz),'=f8')
for j in range(ny):
  for k in range(nz):
    y = (float(j)*L)/(ny-1)
    z = (float(k)*L)/(nz-1)
    rarray[j][k] = T0*(1+(y/L)+(z/L))
h5dset = h5grp.create_dataset('Data', data=rarray)

# 2D Surface:
# -----------------------------------------------------------------------------
# Temperature boundary condition z=0 face
h5grp = h5file.create_group('node_centered_surf_bottom')
h5grp.attrs['Dimension'] = np.string_('XY')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1,0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.,0.]
# Load the dataset values
nx = 11
ny = 11
rarray = np.zeros((nx,ny),'=f8')
for i in range(nx):
  for j in range(ny):
    x = (float(i)*L)/(nx-1)
    y = (float(j)*L)/(ny-1)
    rarray[i][j] = T0*((x/L)+(y/L)+0)
h5dset = h5grp.create_dataset('Data', data=rarray)

# 2D Surface:
# -----------------------------------------------------------------------------
# Temperature boundary condition z=L face
h5grp = h5file.create_group('node_centered_surf_top')
h5grp.attrs['Dimension'] = np.string_('XY')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1,0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.,0.]
# Load the dataset values
nx = 11
ny = 11
rarray = np.zeros((nx,ny),'=f8')
for i in range(nx):
  for j in range(ny):
    x = (float(i)*L)/(nx-1)
    y = (float(j)*L)/(ny-1)
    rarray[i][j] = T0*((x/L)+(y/L)+1)
h5dset = h5grp.create_dataset('Data', data=rarray)

# 2D Surface:
# -----------------------------------------------------------------------------
# Temperature boundary condition y=0 face
h5grp = h5file.create_group('node_centered_surf_south')
h5grp.attrs['Dimension'] = np.string_('XZ')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1,0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.,0.]
# Load the dataset values
nx = 11
nz = 11
rarray = np.zeros((nx,nz),'=f8')
for i in range(nx):
  for k in range(nz):
    x = (float(i)*L)/(nx-1)
    z = (float(k)*L)/(nz-1)
    rarray[i][k] = T0*((x/L)+0+(z/L))
h5dset = h5grp.create_dataset('Data', data=rarray)

# 2D Surface:
# -----------------------------------------------------------------------------
# Temperature boundary condition y=L face
h5grp = h5file.create_group('node_centered_surf_north')
h5grp.attrs['Dimension'] = np.string_('XZ')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [0.1,0.1]
# Location of origin
h5grp.attrs['Origin'] = [0.,0.]
# Load the dataset values
nx = 11
nz = 11
rarray = np.zeros((nx,nz),'=f8')
for i in range(nx):
  for k in range(nz):
    x = (float(i)*L)/(nx-1)
    z = (float(k)*L)/(nz-1)
    rarray[i][k] = T0*((x/L)+1+(z/L))
h5dset = h5grp.create_dataset('Data', data=rarray)
