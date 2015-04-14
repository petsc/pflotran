from h5py import *
import numpy

filename = 'mesh.h5'
h5file = File(filename,mode='w')

nx = 5
ny = 4
nz = 3

# for irregular grid spacing, dx,dy,dz must be specified.
dx = [10.,11.,12.,13.,14.]
dy = [13.,12.,11.,10.]
dz = [15.,20.,25.]
# for uniform grid spacing, specify one value in each direction
dx = [10.]
dy = [10.]
dz = [1.]
# note that you can mix and match irregular spacing along the different axes

nxp1 = nx+1
nyp1 = ny+1
nzp1 = nz+1

x_origin = 0.
y_origin = 0.
z_origin = 0.

element_array = numpy.zeros((nx*ny*nz,9),'=i4')

for k in range(nz):
  for j in range(ny):
    for i in range(nx):
      cell_id = i + j*nx + k*nx*ny
      element_array[cell_id,0] = 8
      element_array[cell_id,1] = i + j*nxp1 + k*nxp1*nyp1
      element_array[cell_id,2] = i + 1 + j*nxp1 + k*nxp1*nyp1
      element_array[cell_id,3] = i + 1 + (j+1)*nxp1 + k*nxp1*nyp1
      element_array[cell_id,4] = i + (j+1)*nxp1 + k*nxp1*nyp1
      element_array[cell_id,5] = i + j*nxp1 + (k+1)*nxp1*nyp1
      element_array[cell_id,6] = i + 1 + j*nxp1 + (k+1)*nxp1*nyp1
      element_array[cell_id,7] = i + 1 + (j+1)*nxp1 + (k+1)*nxp1*nyp1
      element_array[cell_id,8] = i + (j+1)*nxp1 + (k+1)*nxp1*nyp1
h5dset = h5file.create_dataset('Domain/Cells', data=element_array)

vertex_array = numpy.zeros((nxp1*nyp1*nzp1,3),'=f8')

x = 0.
y = 0.
z = 0.
for k in range(nzp1):
  if len(dz) == 1:
    z = k*dz[0] + z_origin
  else:
    if k == 0:
      z = z_origin
    else:
      z += dz[k-1]
  for j in range(nyp1):
    if len(dy) == 1:
      y = j*dy[0] + y_origin
    else:
      if j == 0:
        y = y_origin
      else:
        y += dy[j-1]
    for i in range(nxp1):
      vertex_id = i + j*nxp1 + k*nxp1*nyp1
      if len(dx) == 1:
        x = i*dx[0] + x_origin
      else:
        if i == 0:
          x = x_origin
        else:
          x += dx[i-1]
      vertex_array[vertex_id,0] = x
      vertex_array[vertex_id,1] = y
      vertex_array[vertex_id,2] = z

h5dset = h5file.create_dataset('Domain/Vertices', data=vertex_array)

h5file.close()
print('done')
