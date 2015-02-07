import sys

# Author: Glenn Hammond
# Date: 06/03/13
# This python script calculates the dependencies between PFLOTRAN source files
# and writes the dependencies to pflotran_dependencies.txt.  The developer
# must then insert the contents of this file within the PFLOTRAN makefile.

def get_filename(root,suffix):
  filename = []
  filename.append(root)
  filename.append(suffix)
  filename = '.'.join(filename)
  return filename

def skip_if0(line,f):
  if line.lstrip().startswith('#if 0'):
    if_count = 1
    while 1:
      success, line = get_line(f)
      if not success:
        break
      line = line.strip()
      if line.startswith('#if'):
        if_count += 1
      if line.startswith('#endif'):
        if_count -= 1
        if if_count == 0:
          break

def get_line(f):
  string = ''
  strings = []
  success = True
  for line in f:
    skip_if0(line,f)
    line = line.split('!')[0].strip()
    if line.endswith('&'):
      line = line.rstrip('&').strip()
      strings.append(line)
    else:
      strings.append(line)
      string = ''.join(strings)
      break
  else:
    success = False
#  print(string)
  return success,string
    

class DerivedType:
  def __init__(self,f,name):
    self.name = name
    print('Derived type %s read' % self.name)
    self.read(f)
  def read(self,f):
    end_found = False
    while 1:
      success,line = get_line(f)
      if not success:
        break
      if line.startswith('end type') and line.endswith(self.name):
        end_found = True
        break
    if not end_found:
      print('End of derived type %s not found.' % self.name)
      exit(1)

class Subroutine:
  def __init__(self,f,name):
    self.name = name
    self.modules_used = []
    self.derived_types_used = []
    self.subroutines_used = []
    print('Reading subroutine %s' % self.name)
    self.read(f)
    print('Subroutine %s read' % self.name)
  def add_module(self,name):
    self.modules_used.append(name)
  def add_subroutine(self,name):
    self.subroutines_used.append(name)
  def add_derived_type(self,name):
    self.derived_types_used.append(name)
  def read(self,f):
    while 1:
      success,line = get_line(f)
      if not success:
        break
      elif line.startswith('end subroutine') and line.endswith(self.name):
        end_found = True
        break
    if not end_found:
      print('End of module %s not found.' % self.name)
      exit(1)
#      elif line.startswith('call'):

class Function:
  def __init__(self,f,name):
    self.name = name
    self.modules_used = []
    self.derived_types_used = []
    self.subroutines_used = []
    print('Reading function %s' % self.name)
    self.read(f)
    print('Function %s read' % self.name)
  def add_module(self,name):
    self.modules_used.append(name)
  def add_subroutine(self,name):
    self.subroutines_used.append(name)
  def add_derived_type(self,name):
    self.derived_types_used.append(name)
  def read(self,f):
    end_found = False
    while 1:
      success,line = get_line(f)
#      print(line)
      if not success:
        break
      if line.startswith('end function') and line.endswith(self.name):
        end_found = True
        break
    if not end_found:
      print('End of function %s not found.' % self.name)
      exit(1)
#      elif line.startswith('call'):

class SourceFile:
  def __init__(self,root):
    self.name = get_filename(root,'F90')
    self.modules = []
    self.read()
  def read(self):
#    print('Opening %s' % self.name)
    f = open(self.name)
    print('%s opened' % self.name)
    line = ''
    while 1:
      success,line = get_line(f)
      if not success:
        break
      print(line)
      if line.startswith('module'):
        w = line.split()
        # skip 'module procedure'
        if not w[1].startswith('procedure'):
          self.modules.append(Module(f,w[1]))
    f.close()

class Module:
  def __init__(self,f,name):
    self.name = name
    self.derived_types = []
    self.subroutines = []
    self.functions = []
    print('Reading module %s' % self.name)
    self.read(f)
    print('Module %s read' % self.name)
  def read(self,f):
    line = ''
    end_found = False
    while 1:
      success,line = get_line(f)
      if not success:
        break
#      if line.startswith('end'):
#        w = line.split()
#        if len(w) > 1:
#          card = w[1].strip()
#          if card.startswith('type'):
#          elif card.startswith('function'):
#          elif card.startswith('subroutine'):
      
      if line.startswith('subroutine'):
        temp = line.lstrip('subroutine')
        w = temp.split('(') # remove comments
        subroutine_name = w[0].strip()
        self.subroutines.append(Subroutine(f,subroutine_name))
      elif line.startswith('type'):
        print('---%s' % line)
        temp = line.lstrip('type')
        if not temp.lstrip().startswith('('):
          w = temp.split('!') # remove comments
          w = w[0].split('::') # split at ::
          w = w[1].strip().split() 
          derived_type_name = w[0] # take first word after ::
          self.derived_types.append(DerivedType(f,derived_type_name))
      elif line.find('function ') > -1:
        temp = line.split('!')[0]
        if temp.find('function ') > -1:
          temp = temp.split('function ')[1]
          function_name = temp.split('(')[0].strip()
          self.functions.append(Function(f,function_name))
      elif line.startswith('end module') and line.endswith(self.name):
        end_found = True
        break
    if not end_found:
      print('End of module %s not found.' % self.name)
      exit(1)
    

#print('here')
source_file = SourceFile('waypoint')
print('done')
exit(0)

# Obtain list of source files
source_file_roots = []
source_block = False
for line in open('makefile','r'):
  if line.startswith('# Begin Source Block'):
    source_block = True
  if line.startswith('# End Source Block'):
    source_block = False
  if source_block:
    # find .o file
    # coule use re.split() here, but too complicated.
    w = line.split('}')
    if len(w) == 2:
      w2 = w[1].split('.o')
#      print(w2[0])
      source_file_roots.append(w2[0])
source_file_roots.append('pflotran')

# Alphabetize
source_files = []
source_file_roots.sort()
for root in source_file_roots:
  source_files = SourceFile(root)

print('done')
exit(0)


# for each source file generate a list of derived type and subroutines

# Obtain a list of modules
module_list = []
for root in source_file_roots:
#  print(root)
  for line in open(get_filename(root,'F90')):
    if line.lstrip().startswith('module'):
      w = line.split()
      # skip 'module procedure'
      if not w[1].startswith('procedure'):
#        print('  '+w[1])
        module_link = []
        module_link.append(w[1])
        module_link.append(get_filename(root,'o'))
        module_list.append(module_link)
module_dictionary = dict(module_list)
#print(module_dictionary.keys())

f = open('pflotran_dependencies.txt','w')
# now loop over all source files and create the dependency list
for root in source_file_roots:
#  print(root)
  try:
    modules_to_remove = differing_pflotran_rxn_dependencies[root]
    num_times_to_print = 2
    f.write('ifdef PFLOTRAN_RXN_FLAG\n')
  except KeyError:
    num_times_to_print = 1
  for iprint in range(num_times_to_print):
    module_list = []
    for line in open(get_filename(root,'F90')):
      if line.lstrip().startswith('use '):
        w = line.split()
        module_list.append(w[1].strip(','))
    # remove duplicate modules
    module_list = set(module_list)
    file_list = []
    for module in module_list:
      try:
        key = module_dictionary[module]
      except KeyError:
        # need to skip hdf5
        if not module.startswith('hdf5') and not module.startswith('h5lt'):
          print('Module "%s" not found in dictionary.\n' % module)
          print(root, module)
          sys.exit()
#      print(key)
      file_list.append(key)
    # remove duplicates first as it will destroy an sorting
    file_set = set(file_list)
    # convert back to list
    file_list = list(file_set)
    # sort
    sorted_file_list = sorted(file_list)
    # remove the root file if it is listed
    filename = get_filename(root,'o')
    if filename in sorted_file_list:
      sorted_file_list.remove(filename)
      print('Removing %s from own dependency.\n' % filename)
#    print(sorted_file_list)
    # remove files that are not needed for pflotran rxn
    if num_times_to_print == 2 and iprint == 0:
      for root2 in modules_to_remove:
        sorted_file_list.remove(get_filename(root2,'o'))
#      print(sorted_file_list)
    string = '%s :' % get_filename(root,'o')
    f.write('%s' % string)
    string_len = len(string)
    indentation = string_len
    for file in sorted_file_list:
      string = ' %s' % file
      word_len = len(string)
      string_len += word_len
      if string_len > 78:
        f.write(' \\\n')
        for i in range(indentation):
          f.write(' ')
        string_len  = indentation + word_len
      f.write('%s' % string)
    f.write('\n')
    if num_times_to_print == 2 and iprint == 0:
      f.write('else\n')
    elif num_times_to_print == 2: 
      f.write('endif\n')
f.close()  
print('done!')
