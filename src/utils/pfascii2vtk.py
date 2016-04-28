#!/usr/bin/env python

import sys
import array
import visit_writer

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if len(argv) > 1:
    pflowdat = open(argv[1], 'r')
  if len(argv) > 2:
    vlxflag = True
    pflowvlx = open(argv[2], 'r')
  if len(argv) > 3:
    vlyflag = True
    pflowvly = open(argv[3], 'r')
    
  # Read the header information from the pflow.out file.
  line = pflowdat.readline()
  title = line.split('"')[1]
  line = pflowdat.readline()
  line = line.rstrip()
  line = line.replace("VARIABLES=", "")
  line = line.replace('"', '')
  varnames = line.split(",")
  line = pflowdat.readline()
  t = float(line.split('"')[1])
  nx = int(line.split()[5])
  ny = int(line.split()[8])
  nz = int(line.split()[10])
  
  #x = array.array('f')
  x = []
  #y = array.array('f')
  y = []
  #z = array.array('f')
  z = []
  var = [None] * (len(varnames) - 3)
  for i in range(len(var)):
    #var[i] = array.array('f')
    var[i] = []

  for k in range(nz):
    for j in range(ny):
      for i in range(nx):
        line = pflowdat.readline()
        if line == "":
          print "Warning!  Encountered EOF early!"
          break
        fields = line.split()
        x.append(float(fields[0]))
        y.append(float(fields[1]))
        z.append(float(fields[2]))
        for i in range(len(var)):
          var[i].append(float(fields[i+3]))

  visitvars = [None] * len(var)
  for i in range(len(var)):
    visitvars[i] = (varnames[i], 1, 1, var[i])
  visit_writer.WriteRectilinearMesh(argv[1]+".vtk", 0, x, y, z, visitvars)  

if __name__ == "__main__":
  sys.exit(main())
