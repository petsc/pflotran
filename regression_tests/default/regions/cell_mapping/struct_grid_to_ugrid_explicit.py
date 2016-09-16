from h5py import *
import numpy

nx = 5
ny = 4
nz = 3

# for irregular grid spacing, dx,dy,dz must be specified.
dx = [10.,11.,12.,13.,14.]
dy = [13.,12.,11.,10.]
dz = [15.,20.,25.]
# for uniform grid spacing, specify one value in each direction
#dx = [1.]
#dy = [1.]
#dz = [1.]
# note that you can mix and match irregular spacing along the different axes

nxp1 = nx+1
nyp1 = ny+1
nzp1 = nz+1

x_origin = 0.
y_origin = 0.
z_origin = 0.

f = open('543.uge','w')

f.write('CELLS %d\n' % (nx*ny*nz))
z = z_origin + 0.5*dz[0]
for k in range(nz):
  if len(dz) == 1:
    lenz = dz[0]
  else:
    lenz = dz[k]
  y = y_origin + 0.5*dy[0]
  for j in range(ny):
    if len(dy) == 1:
      leny = dy[0]
    else:
      leny = dy[j]
    x = x_origin + 0.5*dx[0]
    for i in range(nx):
      cell_id = i + j*nx + k*nx*ny
      if len(dx) == 1:
        lenx = dx[0]
      else:
        lenx = dx[i]
      volume = lenx*leny*lenz
      f.write('%d %f %f %f %f\n' % (cell_id + 1,x,y,z,volume))
      if len(dx) == 1:
        x += dx[0]
      else:
        x += 0.5*dx[i]
        if i+1 < nx:
          x += 0.5*dx[i+1]
    if len(dy) == 1:
      y += dy[0]
    else:
      y += 0.5*dy[j]
      if j+1 < ny:
        y += 0.5*dy[j+1]
  if len(dz) == 1:
    z += dz[0]
  else:
    z += 0.5*dz[k]
    if k+1 < nz:
      z += 0.5*dz[k+1]

f.write('CONNECTIONS %d\n' % ((nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1)))

# x connections
z = z_origin + 0.5*dz[0]
for k in range(nz):
  if len(dz) == 1:
    lenz = dz[0]
  else:
    lenz = dz[k]
  y = y_origin + 0.5*dy[0]
  for j in range(ny):
    if len(dy) == 1:
      leny = dy[0]
    else:
      leny = dy[j]
    x = x_origin + dx[0]
    for i in range(nx-1):
      cell_id = i + j*nx + k*nx*ny
      area = leny*lenz
      f.write('%d %d %f %f %f %f\n' % (cell_id + 1,
                                       cell_id + 2,
                                       x,y,z,area))
      if len(dx) == 1:
        x += dx[0]
      else:
        x += dx[i]
    if len(dy) == 1:
      y += dy[0]
    else:
      y += 0.5*dy[j]
      if j+1 < ny:
        y += 0.5*dy[j+1]
  if len(dz) == 1:
    z += dz[0]
  else:
    z += 0.5*dz[k]
    if k+1 < nz:
      z += 0.5*dz[k+1]

# y connections
z = z_origin + 0.5*dz[0]
for k in range(nz):
  if len(dz) == 1:
    lenz = dz[0]
  else:
    lenz = dz[k]
  x = x_origin + 0.5*dx[0]
  for i in range(nx):
    if len(dx) == 1:
      lenx = dx[0]
    else:
      lenx = dx[i]
    y = y_origin + dy[0]
    for j in range(ny-1):
      cell_id = i + j*nx + k*nx*ny
      area = lenx*lenz
      f.write('%d %d %f %f %f %f\n' % (cell_id + 1,
                                       cell_id + 1 + nx,
                                       x,y,z,area))
      if len(dy) == 1:
        y += dy[0]
      else:
        y += dy[j]
    if len(dx) == 1:
      x += dx[0]
    else:
      x += 0.5*dx[i]
      if i+1 < nx:
        x += 0.5*dx[i+1]
  if len(dz) == 1:
    z += dz[0]
  else:
    z += 0.5*dz[k]
    if k+1 < nz:
      z += 0.5*dz[k+1]

# z connections
y = y_origin + 0.5*dy[0]
for j in range(ny):
  if len(dy) == 1:
    leny = dy[0]
  else:
    leny = dy[j]
  x = x_origin + 0.5*dx[0]
  for i in range(nx):
    if len(dx) == 1:
      lenx = dx[0]
    else:
      lenx = dx[i]
    z = z_origin + dz[0]
    for k in range(nz-1):
      cell_id = i + j*nx + k*nx*ny
      area = leny*lenx
      f.write('%d %d %f %f %f %f\n' % (cell_id + 1,
                                       cell_id + 1 + nx*ny,
                                       x,y,z,area))
      if len(dz) == 1:
        z += dz[0]
      else:
        z += dz[k]
    if len(dx) == 1:
      x += dx[0]
    else:
      x += 0.5*dx[i]
      if i+1 < nx:
        x += 0.5*dx[i+1]
  if len(dy) == 1:
    y += dy[0]
  else:
    y += 0.5*dy[j]
    if j+1 < ny:
      y += 0.5*dy[j+1]

f.close()

print('done')
