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
  infilename = 'DBH_Mesh_C.inp'
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
        node_coordinates = [0.]*3
        for icoord in range(3):
          node_coordinates[icoord] = float(w[icoord+1])
        nodes.append(node_coordinates)
        if len(nodes) % 10000 == 0:
          print('Node %d' % len(nodes))
    elif line.startswith('*ELEMENT'):
      while(True):
        line = abaqus_file.readline()
        if line.startswith('*'):
          break
        w = line.split(',')
        node_ids = [0.]*8
        for inode in range(8):
          node_ids[inode] = int(w[inode+1])
        elements.append(node_ids)
        if len(elements) % 10000 == 0:
          print('Element %d' % len(elements))
    if len(nodes) > 0 and len(elements) > 0:
      break
  abaqus_file.close()

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
  for ielem in range(len(elements)):
    int_array[ielem][1:9] = elements[ielem][:]
  # create element data set
  h5_file.create_dataset('Domain/Cells',data=int_array)

  print('Creating material IDs.')
  num_elements = len(elements)
  int_array = numpy.arange(1,num_elements+1,1)
  h5dataset = h5_file.create_dataset('Materials/Cell Ids', data=int_array)
  int_array = numpy.zeros((num_elements),'int')
  int_array[:] = 1
  h5dataset = h5_file.create_dataset('Materials/Material Ids', data=int_array)

  h5_file.close()
  print('done')
  
if __name__ == "__main__":
  abaqus_to_pflotran_mesh()
