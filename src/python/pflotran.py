# pflotran.py
import numpy as np
import sys

class FileType:
  NULL = -1
  OBSERVATION = 0
  TECPLOT_POINT = 1
  TECPLOT_BLOCK = 2

def get_full_paths(paths,filenames):
  full_paths = []
  for ipath in range(len(paths)):
    for ifilename in range(len(filenames)):
      full_path = []
      full_path.append(paths[ipath])
      full_path.append(filenames[ifilename])
      full_path = '/'.join(full_path)
      full_paths.append(full_path)
  return full_paths

def get_tec_filenames(ids):
  filenames = []
  for i in range(len(ids)):
    ifile = ids[i]
    if ifile < 10:
      filename = 'pflotran-00%d.tec' % ifile
    elif ifile < 100:
      filename = 'pflotran-0%d.tec' % ifile
    elif ifile < 1000:
      filename = 'pflotran-%d.tec' % ifile
    else:
      print('File numbers greater than 999 not handled')
      exit(0)
    filenames.append(filename)
  return filenames
      

class Dataset:
  def __init__(self,filename,xcol,ycol):
    self.filetype = FileType.NULL
    self.dictionary = {}
    self.title = ''
    self.variables = []
    self.arrays = []

    self.f = open(filename,'r')
    self.read_file_header()
    # xcol and ycol are 1-based
    self.dictionary['xname'] = self.variables[xcol-1]
    self.dictionary['yname'] = self.variables[ycol-1]
    self.read_dataset_from_columns(xcol,ycol)

  def get_array(self,dictionary_entry):
    iarray = self.dictionary[dictionary_entry]
    return self.arrays[iarray]

  def get_name(self,dictionary_entry):
    return self.dictionary[dictionary_entry]

  def read_file_header(self):
    for line in self.f:
      if line.strip().startswith('"'):
        # observation file
        w = line.strip().split(',')
        for i in range(len(w)):
          self.variables.append(w[i].strip('"'))
        self.filetype = FileType.OBSERVATION
        return
      elif line.strip().startswith('VARIABLES'):
        # read variable names
        line = line.strip().split('S=')[1]
        w = line.strip().split(',')
        for i in range(len(w)):
          self.variables.append(w[i].strip('"'))
      elif line.strip().startswith('ZONE'):
        w = line.strip().split(',')
        zone_name = w[0].split('=')[1]
        self.title = zone_name.strip('"')
        for i in range(len(w)):
          if w[i].strip().startswith('DATAPACKING'):
            s = w[i].split('=')[1]
            if s.endswith('POINT'):
              self.filetype = FileType.TECPLOT_POINT
            elif s.endswith('BLOCK'):
              self.filetype = FileType.TECPLOT_BLOCK
              print('Block datapacking not yet supported')
              exit(0)
            else:
              print('Datapacking method undefined')
              exit(0)
        return

  def read_dataset_from_columns(self,xcol,ycol):
    size = 100
    arrays = []
    array1 = np.zeros(size,'=f8')
    array2 = np.zeros(size,'=f8')
    count = 0
    for line in self.f:
      if count >= size:
        size *= 2
        array1.resize(size)
        array2.resize(size)
      w = line.split()
      # xcol and ycol are 1-based
      array1[count] = float(w[xcol-1])
      array2[count] = float(w[ycol-1])
      count += 1
    array1.resize(count)
    array2.resize(count)
    self.arrays.append(array1)
    self.arrays.append(array2)
    self.dictionary['x'] = 0
    self.dictionary['y'] = 1
