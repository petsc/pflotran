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

pflotran_rxn_list = []
remove_file_list = []
remove_file_list.append('logging')
pflotran_rxn_list.append(('constraint',remove_file_list))
remove_file_list = []
remove_file_list.append('co2_span_wagner')
remove_file_list.append('co2eos')
remove_file_list.append('water_eos')
pflotran_rxn_list.append(('reaction',remove_file_list))
differing_pflotran_rxn_dependencies = dict(pflotran_rxn_list)

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
source_file_roots.sort()
#print(source_file_roots)

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
        if not module.startswith('hdf5'):
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
    # remove files that are not needed for pfltoran rxn
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
