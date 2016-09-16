from h5py import *
import numpy

def explicit_face_to_implicit(iface):
  pft_to_exo_face = [-999,0,2,1,3,4,5] # pft is 1-based
  # this is rotated 90 degrees counter clockwise around Z
  indices = -999
  if pft_to_exo_face[iface] == 0: # west
    indices = [3,0,4,7]
  elif pft_to_exo_face[iface] == 1: # south
    indices = [0,1,5,4]
  elif pft_to_exo_face[iface] == 2: # east
    indices = [1,2,6,5]
  elif pft_to_exo_face[iface] == 3: # north
    indices = [2,3,7,6]
  elif pft_to_exo_face[iface] == 4: # bottom
    indices = [3,2,1,0]
  elif pft_to_exo_face[iface] == 5: # top
    indices = [4,5,6,7]
  else:
    print('Error in explicit_face_to_implicit')
  return indices

def map_sideset_from_explicit(filename,elements):
  faces = []
  for line in open(filename,'r'):
    w = line.split()
    cell_id = int(w[0])-1
    vertex_ids = explicit_face_to_implicit(int(w[1]))
#    print(vertex_ids)
    vertices = []
#    print(elements[cell_id,1:9])
#    print(cell_id,elements[cell_id])
    for i in range(len(vertex_ids)):
      vertices.append(elements[cell_id,vertex_ids[i]+1])
#    print(vertices)
    faces.append(vertices)
  return faces

# get material ids from ../543/543.h5
filename = '../543/543.h5'
h5file1 = File(filename,mode='r')
material_ids = h5file1['Materials/Material Ids']
#print(material_ids[:])

filename = '../543/parameters-543.h5'
h5file2 = File(filename,mode='r')
permeability = h5file2['Permeability']
porosity = h5file2['Porosity']

filename = '543_ugi.h5'
h5file = File(filename,mode='w')
filename = '543.ugi'
f_ascii = open(filename,'w')

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

element_array = numpy.zeros((nx*ny*nz,9),'=i4')

for k in range(nz):
  for j in range(ny):
    for i in range(nx):
      element_id = i + j*nx + k*nx*ny
      element_array[element_id,0] = 8
      element_array[element_id,1] = i + j*nxp1 + k*nxp1*nyp1
      element_array[element_id,2] = i + 1 + j*nxp1 + k*nxp1*nyp1
      element_array[element_id,3] = i + 1 + (j+1)*nxp1 + k*nxp1*nyp1
      element_array[element_id,4] = i + (j+1)*nxp1 + k*nxp1*nyp1
      element_array[element_id,5] = i + j*nxp1 + (k+1)*nxp1*nyp1
      element_array[element_id,6] = i + 1 + j*nxp1 + (k+1)*nxp1*nyp1
      element_array[element_id,7] = i + 1 + (j+1)*nxp1 + (k+1)*nxp1*nyp1
      element_array[element_id,8] = i + (j+1)*nxp1 + (k+1)*nxp1*nyp1

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

# create side sets based on original vertex numbering and then convert to new
ss_filenames = []
ss_filenames.append('west_inactive.txt')
ss_filenames.append('east_inactive.txt')
ss_filenames.append('south_inactive.txt')
ss_filenames.append('north_inactive.txt')
ss_filenames.append('bottom_inactive.txt')
ss_filenames.append('top_inactive.txt')
ss = []
for filename in ss_filenames:
 ss.append(map_sideset_from_explicit(filename,element_array))

# re-map elements and vertex ids accounting for inactive cells
iactive_element_array = numpy.ones((nx*ny*nz),'=i4')
# inactivate following elements
iactive_element_array[28:31] = 0
iactive_element_array[33:37] = 0
iactive_element_array[38:40] = 0
iactive_element_array[41:43] = 0
iactive_element_array[45:48] = 0
iactive_element_array[50:] = 0
#print(iactive_element_array)

iactive_vertex_array = numpy.zeros((nxp1*nyp1*nzp1),'=i4')
vertex_map = numpy.zeros((nxp1*nyp1*nzp1),'=i4')
vertex_map[:] = -999
active_element_count = 0
for element_id in range(nx*ny*nz):
  if iactive_element_array[element_id] > 0:
    active_element_count += 1
    for i in range(element_array[element_id,0]):
      iactive_vertex_array[element_array[element_id,i+1]] = 1
#print(iactive_vertex_array)
active_vertex_count = 0
for vertex_id in range(nxp1*nyp1*nzp1):
  if iactive_vertex_array[vertex_id] > 0:
    vertex_map[vertex_id] = active_vertex_count
    active_vertex_count += 1
#print(vertex_map)

print(active_element_count)
print(active_vertex_count)
new_element_array = numpy.zeros((active_element_count,9),'=i4')
new_vertex_array = numpy.zeros((active_vertex_count,3),'=f8')
new_material_ids = numpy.zeros(active_element_count,'=i4')
cell_ids = numpy.arange(1,active_element_count+1,dtype='=i4')
new_permeability = numpy.zeros(active_element_count,'=f8')
new_porosity = numpy.zeros(active_element_count,'=f8')

ielement = 0
for element_id in range(nx*ny*nz):
  if iactive_element_array[element_id] > 0:
    new_material_ids[ielement] = material_ids[element_id]
    new_permeability[ielement] = permeability[element_id]
    new_porosity[ielement] = porosity[element_id]
    new_element_array[ielement,0] = element_array[element_id,0]
    for i in range(1,9):
      new_element_array[ielement,i] = vertex_map[element_array[element_id,i]]
    ielement += 1
ivertex = 0
for vertex_id in range(nxp1*nyp1*nzp1):
  if iactive_vertex_array[vertex_id] > 0:
    new_vertex_array[ivertex,:] = vertex_array[vertex_id,:]
    ivertex += 1

# map sideset vertices
for ifile in range(len(ss_filenames)):
  f = open(ss_filenames[ifile].split('.')[0]+'_ugi_before.ss','w')
  array = ss[ifile]
  f.write('%d\n' % len(array))
  for i in range(len(array)):
    f.write('Q')
    for j in range(4):
      f.write(' %d' % (array[i][j]+1))
    f.write('\n')
  f.close()
  for i in range(len(array)):
    for j in range(len(array[i])):
      array[i][j] = vertex_map[array[i][j]]

element_array = new_element_array
vertex_array = new_vertex_array 
material_ids = new_material_ids
permeability = new_permeability
porosity = new_porosity

h5dset = h5file.create_dataset('Domain/Cells', data=element_array)
h5dset = h5file.create_dataset('Domain/Vertices', data=vertex_array)
h5dset = h5file.create_dataset('Materials/Cell Ids', data=cell_ids)
h5dset = h5file.create_dataset('Cell Ids', data=cell_ids)
h5dset = h5file.create_dataset('Materials/Material Ids', data=material_ids)
h5dset = h5file.create_dataset('Permeability', data=permeability)
h5dset = h5file.create_dataset('Porosity', data=porosity)

f_ascii.write('%d %d\n' % (active_element_count,active_vertex_count))
for element_id in range(active_element_count):
  f_ascii.write('H')
  for iv in range(element_array[element_id][0]):
    f_ascii.write(' %d' % (element_array[element_id,iv+1]+1))
  f_ascii.write('\n')
for i in range(active_vertex_count):
  f_ascii.write('%12.6e %12.6e %12.6e\n' % 
                (vertex_array[i,0],vertex_array[i,1],vertex_array[i,2]))

# write sidesets
for ifile in range(len(ss_filenames)):
  f = open(ss_filenames[ifile].split('.')[0]+'_ugi.ss','w')
  array = ss[ifile]
  f.write('%d\n' % len(array))
  for i in range(len(array)):
    f.write('Q')
    for j in range(4):
      f.write(' %d' % (array[i][j]+1))
    f.write('\n')
  f.close()

h5file.close()
f_ascii.close()
h5file1.close()
h5file2.close()
print('done')
