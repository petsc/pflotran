import sys
import math
from h5py import *
import numpy

h5filename = '543_regions.h5'
h5file = File(h5filename,mode='w')

filename = '543_cell_by_cell_regions_hdf5.txt'
f_region_block_h5 = open(filename,'w')
filename = '543_cell_by_cell_regions_ascii.txt'
f_region_block_ascii = open(filename,'w')
filename = '543_transport_block.txt'
f_transport_block = open(filename,'w')
filename = '543_condition_cell_by_cell_block.txt'
f_condition_cbc_block = open(filename,'w')
filename = '543_condition_eighths_block.txt'
f_condition_eighths_block = open(filename,'w')

nx = 5
ny = 4
nz = 3
n = nx*ny*nz

region_group = h5file.create_group("Regions")

region_template = '''REGION %s
  FILE ./%s
END
'''

transport_condition_template = '''TRANSPORT_CONDITION %s
  TYPE DIRICHLET
  CONSTRAINT %s
    CONCENTRATIONS
      Tracer  %9.3e  F
    /
  /
END
'''

condition_template = '''INITIAL_CONDITION
  TRANSPORT_CONDITION %s
  REGION %s
END
'''

for k in range(nz):
  for j in range(ny):
    for i in range(nx):
      index = i + j*nx + k*nx*ny + 1
      # hdf5
      region_name = '%d' % index
      region_string = region_name
      new_group = region_group.create_group(region_name)
      iarray = numpy.zeros(1,'=i4')
      iarray[0] = index
      dataset_name = 'Cell Ids'
      new_group.create_dataset(dataset_name, data=iarray)
      # ascii
      region_name = region_name + '.txt'
      f_region = open(region_name,'w')
      f_region.write('%d\n' % index)
      f_region.close()
      # region block
      f_region_block_h5.write(region_template % (region_string,h5filename))
      f_region_block_ascii.write(region_template % (region_string,region_name))
      concentration = 1.e-5 * float(index)
      f_transport_block.write(transport_condition_template % 
                              (region_string,region_string,concentration))
      f_condition_cbc_block.write(condition_template % 
                                  (region_string,region_string))

# divide into eighths
istart = [0,2]
iend = [2,5]
jstart = [0,2]
jend = [2,4]
kstart = [0,1]
kend = [1,3]

region_names = ['bottom_south_west','bottom_south_east',
                'bottom_north_west','bottom_north_east',
                'top_south_west','top_south_east',
                'top_north_west','top_north_east']

iregion = 0
for kk in range(2):
  for jj in range(2):
    for ii in range(2):
      n = (iend[ii]-istart[ii])*(jend[jj]-jstart[jj])*(kend[kk]-kstart[kk])
      iarray = numpy.zeros(n,'=i4')
      count = 0
      for k in range(kstart[kk],kend[kk]):
        for j in range(jstart[jj],jend[jj]):
          for i in range(istart[ii],iend[ii]):
            index = i + j*nx + k*nx*ny + 1
            iarray[count] = index
            count += 1
      # hdf5
      region_name = region_names[iregion]
      new_group = region_group.create_group(region_name)
      dataset_name = 'Cell Ids'
      new_group.create_dataset(dataset_name, data=iarray)
      # ascii
      region_name = region_name + '.txt'
      f_region = open(region_name,'w')
      for iii in range(n):
        f_region.write('%d\n' % iarray[iii])
      f_region.close()
      # region block
      f_region_block_h5.write(region_template % 
                              (region_names[iregion],h5filename))
      f_region_block_ascii.write(region_template % 
                                 (region_names[iregion],region_name))
      concentration = 1.e-5 * float(iregion+1)
      f_transport_block.write(transport_condition_template % 
                              (region_names[iregion],region_names[iregion],
                               concentration))
      f_condition_eighths_block.write(condition_template % 
                                      (region_names[iregion],
                                       region_names[iregion]))
      iregion += 1 
      
h5file.close()
f_region_block_h5.close()
f_region_block_ascii.close()
f_transport_block.close()
f_condition_cbc_block.close()
f_condition_eighths_block.close()

print 'done with everything'
  
