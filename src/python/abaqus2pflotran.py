#-------------------------------------------------------------------------------
# Tool to convert an Abaqus formatted mesh to PFLOTRAN HDF5 mesh
# Usage: python abaqus2pflotran.py abaqus_mesh.txt
# 
# Author: Glenn Hammond (gehammo@sandia.gov)
# Applied Systems Analysis and Research, 6224
# Sandia National Laboratories
# Date: 04/08/15
#-------------------------------------------------------------------------------
from h5py import * 
import numpy 
import os
import sys
import time

# ASSUMPTIONS:
# 1. Single set of nodes and elements in file denoted by '*NODE' and '*ELEMENT'
# 2. Hex elements (8 nodes per elements)
# 3. Node and element IDs are in ascending order

def abaqus_to_pflotran_mesh():
  start_time = time.time()
#  infilename = sys.argv[1]
  infilename = 'filename.inp'
  outfilename = infilename.split('.')[0] + '_usg.h5'
#  if len(sys.argv) != 2:
#    print("ERROR: Command line arguments not provided.\n")
#    print("Usage: python abaqus2pflotran.py abaqus_mesh.txt")
#    sys.exit(0)

  print("Infile: %s Outfile: %s" %(infilename, outfilename))

  nodes = []
  elements = []
  # open Abaqus file and read in nodes and elements
  abaqus_file = open(infilename,'r')
  while(1):
    line = abaqus_file.readline()
    if line.startswith('*NODE'):
      while(True):
        line = abaqus_file.readline()
        if line.startswith('*'):
          break
        w = line.split(',')
        node_id = int(w[0])
        node_coordinates = [0.]*3
        for icoord in range(3):
          node_coordinates[icoord] = float(w[icoord+1])
        nodes.append(node_coordinates)
        if node_id != len(nodes):
          sys.exit('Non-contiguous node ids.')
        if len(nodes) % 10000 == 0:
          print('Node %d' % len(nodes))
    elif line.startswith('*ELEMENT'):
      w = line.split('ELSET=')
      material_id = int(w[1].lstrip('EB'))
      while(True):
        line = abaqus_file.readline()
        if line.startswith('*ELEMENT'):
          print(line)
          w = line.split('ELSET=')
          material_id = int(w[1].lstrip('EB'))
          continue
        elif line.startswith('**'):
          break
        w = line.split(',')
        element_id = int(w[0])
        node_ids = [0.]*8
        for inode in range(8):
          node_ids[inode] = int(w[inode+1])
        elements.append([[element_id,material_id],node_ids])
        if len(elements) % 10000 == 0:
          print('Element %d' % len(elements))
    else:
      if len(nodes) > 0 and len(elements) > 0:
        break
  num_materials = material_id #save largest number, ers 10.20.15
  abaqus_file.close()

  # Create mapping to output elements  and material ids in ascending order
  mapping = numpy.zeros((len(elements)),'i4')
  for i in range(len(elements)):
    element_id = elements[i][0][0]
    mapping[element_id-1] = i

  # Open PFLOTRAN HDF5 Mesh file
  h5_file = File(outfilename, mode='w')

  # write node coordinates
  print('Creating Vertices data set.')
  float_array = numpy.zeros((len(nodes),3),'f8')
  for inode in range(len(nodes)):
    float_array[inode][:] = nodes[inode][:]
  # create node data set
  h5_file.create_dataset('Domain/Vertices',data=float_array)

  # write element nodes  
  print('Creating Cells data set.')
  int_array = numpy.zeros((len(elements),9),'i4')
  # 8 vertices per cell
  int_array[:][:] = 8
  for i in range(len(elements)):
    element_nodes = elements[mapping[i]][1]
    int_array[i][1:9] = element_nodes
    if i+1 != elements[mapping[i]][0][0]:
      sys.exit('Non-contiguous element ids.')

  if numpy.amax(int_array) > len(nodes):
    sys.exit('Node ids do not match elements.')
  # create element data set
  h5_file.create_dataset('Domain/Cells',data=int_array)

  print('Creating material IDs.')
  num_elements = len(elements)
  int_array = numpy.arange(1,num_elements+1,1)
  h5dataset = h5_file.create_dataset('Materials/Cell Ids', data=int_array)
  int_array = numpy.zeros((num_elements),'int')
  for i in range(num_elements):
    material_id = elements[mapping[i]][0][1]
    int_array[i] = material_id
  h5dataset = h5_file.create_dataset('Materials/Material Ids', data=int_array)

  # write regions, label them with material numbers, ers 10.20.15 
  #geh: added mapping[i] to ensure ascending order
  #geh: removed num_elements_in_mat[] as it is not needed
  print ('Creating regions from materials.')
  for i in range(num_materials):
    elements_in_mat = []  
    for j in range(len(elements)):
      if elements[mapping[j]][0][1] == i+1: # added mapping 
        elements_in_mat.append(elements[mapping[j]][0][0]) # added mapping
    int_array = numpy.zeros(len(elements_in_mat),'int')
    for j in range(len(elements_in_mat)): # use len() to determine size
      int_array[j] = elements_in_mat[j] #suspect this is really inefficient 
    dataset_name = 'Regions/Region%i/Cell Ids' %(i+1)
    h5_file.create_dataset(dataset_name,data=int_array)
  

  h5_file.close()
  print('done')
  
if __name__ == "__main__":
  abaqus_to_pflotran_mesh()
