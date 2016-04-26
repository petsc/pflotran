import sys
import math
from h5py import *
import numpy

# Region object
class Region:
  def __init__(self,region_group,region_name,filename):
    self.region_name = region_name
    self.filename = filename
    self.group = region_group.create_group(self.region_name)
  def writeRegion(self,n):
    # in text file, region data should be specified as follows
    # cell_id face_id
    # one connection per line
    cell_id_array = numpy.zeros(n,'=i4')
    face_id_array = numpy.zeros(n,'=i4')
    f = open(self.filename)
    count = 0
    while (1):
      s = f.readline()
      if len(s) < 2:
        break
      w = s.split()
      cell_id_array[count] = int(w[0])
      face_id_array[count] = int(w[1])
      count += 1
    iarray = numpy.zeros(count,'=i4')
    iarray[0:count] = cell_id_array[0:count]
    dataset_name = 'Cell Ids'
    self.group.create_dataset(dataset_name, data=iarray)
    cell_id_array = 0
    iarray[0:count] = face_id_array[0:count]
    dataset_name = 'Face Ids'
    self.group.create_dataset(dataset_name, data=iarray)
    face_id_array = 0
    iarray = 0
    print 'done with Region:', self.region_name


filename = 'material_and_regions.h5'
h5file = File(filename,mode='w')

nx = 2
ny = 2
nz = 2
n = nx*ny*nz

# material ids
# create the Materials group
materials_group = h5file.create_group("Materials")

# create integer array for specifying cell and material ids
iarray = numpy.zeros((nx*ny*nz),'=i4')

# add cell ids to Materials group
for i in range(n):
  iarray[i] = i+1
dataset_name = 'Cell Ids'
h5dset = materials_group.create_dataset(dataset_name, data=iarray)

# add material ids to Materials group
filename = 'material_ids.txt'
f = open(filename,'r')
for i in range(n):
  s = f.readline()
  iarray[i] = int(s)
dataset_name = 'Material Ids'
h5dset = materials_group.create_dataset(dataset_name, data=iarray)
f.close()
iarray = 0 # hopefully frees memory

region_group = h5file.create_group("Regions")
# For each region call:

# Region 1
filename = 'region1.txt'
region_name = filename.split('.')[0]
region = Region(region_group,region_name,filename)
region.writeRegion(n)

# Region 2
filename = 'region2.txt'
region_name = filename.split('.')[0]
region = Region(region_group,region_name,filename)
region.writeRegion(n)

# ....

h5file.close()
print 'done with everything'
  
