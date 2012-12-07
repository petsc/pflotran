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
  return get_tec_filenames('pflotran',ids)
      
def get_tec_filenames_with_prefix(prefix,ids):
  filenames = []
  for i in range(len(ids)):
    ifile = ids[i]
    if ifile < 10:
      filename = '%s-00%d.tec' % (prefix,ifile)
    elif ifile < 100:
      filename = '%s-0%d.tec' % (prefix,ifile)
    elif ifile < 1000:
      filename = '%s-%d.tec' % (prefix,ifile)
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
    self.var_dict = {}
    self.arrays = []

    self.f = open(filename,'r')
    self.read_file_header()

    # if variables are strings, index to ids
    if isinstance(xcol,str):
      xcol_id = self.var_dict[xcol.strip()]
    else:
      # xcol and ycol are 1-based
      xcol_id = xcol-1
    if isinstance(ycol,str):
      ycol_id = self.var_dict[xcol.strip()]
    else:
      ycol_id = ycol-1

    if self.filetype == FileType.TECPLOT_BLOCK:
      self.dictionary['zname'] = self.variables[xcol_id]
      self.read_dataset_from_block(xcol_id)
      if isinstance(xcol,str):
        self.dictionary[xcol] = self.dictionary['z']
      self.dictionary['X'] = self.dictionary['x']
      self.dictionary['Y'] = self.dictionary['y']
    else:
      self.dictionary['xname'] = self.variables[xcol_id]
      self.dictionary['yname'] = self.variables[ycol_id]
      self.read_dataset_from_columns(xcol_id,ycol_id)
      if isinstance(xcol,str):
        self.dictionary[xcol] = self.dictionary['x']
      if isinstance(ycol,str):
        self.dictionary[ycol] = self.dictionary['y']
        
    self.f.close()

  def get_array(self,dictionary_entry):
    try:
      iarray = self.dictionary[dictionary_entry]
    except KeyError:
      print('Array %s not stored in Dataset object' % dictionary_entry)
      sys.exit(1)
    return self.arrays[iarray]

  def get_name(self,dictionary_entry):
    try:
      iarray = self.dictionary[dictionary_entry]
    except KeyError:
      print('Array %s not stored in Dataset object' % dictionary_entry)
      sys.exit(1)
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
          variable = w[i].strip('"')
          self.variables.append(variable)
          self.var_dict[variable] = i
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
            else:
              print('Datapacking method undefined')
              exit(0)
        if self.filetype == FileType.TECPLOT_BLOCK:
          nx = int(w[1].split("=")[1])-1
          ny = int(w[2].split("=")[1])-1
          nz = int(w[3].split("=")[1])-1
          if nz > 1:
            print('Tecplot BLOCK format is only supported for 2D problems in X-Y.')
            sys.exit()
          self.read_discretization(nx,ny,nz)
        return

  def read_discretization(self,nx,ny,nz):
    x = np.zeros(nx,'=f8')
    y = np.zeros(ny,'=f8')
    n = (nx+1)*(ny+1)*(nz+1)
    temp_array = np.zeros(n,'=f8')
    # X
    count = 0
    for line in self.f:
      w = line.split()
      for i in range(len(w)):
        temp_array[count] = float(w[i])
        count += 1
      if count >= n:
        break
    for i in range(nx):
      x[i] = 0.5*(temp_array[i]+temp_array[i+1])
    # Y
    count = 0
    for line in self.f:
      w = line.split()
      for i in range(len(w)):
        temp_array[count] = float(w[i])
        count += 1
      if count >= n:
        break
    for i in range(ny):
      y[i] = 0.5*(temp_array[i*ny]+temp_array[(i+1)*ny])
    # Z
    count = 0
    for line in self.f:
      w = line.split()
      for i in range(len(w)):
        temp_array[count] = float(w[i])
        count += 1
      if count >= n:
        break
    self.dictionary['x'] = len(self.arrays)
    self.arrays.append(x)
    self.dictionary['y'] = len(self.arrays)
    self.arrays.append(y)
  
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
      array1[count] = float(w[xcol])
      array2[count] = float(w[ycol])
      count += 1
    array1.resize(count)
    array2.resize(count)
    self.dictionary['x'] = len(self.arrays)
    self.arrays.append(array1)
    self.dictionary['y'] = len(self.arrays)
    self.arrays.append(array2)

  def read_dataset_from_block(self,ivar):
    nx = len(self.arrays[self.dictionary['x']])
    ny = len(self.arrays[self.dictionary['y']])
    n = nx*ny
    array = np.zeros(n,'=f8')

    lines_per_var = n/10
    if not n == lines_per_var*10:
      lines_per_var += 1

    # skip to variable of interest
    for i in range(ivar-3):
      count = 0
      for line in self.f:
        count += 1
        if count >= lines_per_var:
          break

    count = 0
    for line in self.f:
      w = line.split()
      for i in range(len(w)):
        array[count] = float(w[i])
        count += 1
        if count >= n:
          break
      if count >= n:
        break

    self.dictionary['z'] = len(self.arrays)
    self.arrays.append(array)
    
