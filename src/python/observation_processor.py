import sys
import math

paths = []
paths.append('.')

names = []
names.append('XXX-obs-0.tec')

filenames = []
for ipath in range(len(paths)):
  for iname in range(len(names)):
    filename = []
    filename.append(paths[ipath])
    filename.append(names[iname])
    filename = '/'.join(filename)
    filenames.append(filename)

files = []
# open files and determine column ids
for ifilename in range(len(filenames)):
  files.append(open(filenames[ifilename],'r'))
  header = files[ifilename].readline()
  if ifilename == 0: # only read header from first file
    while(1):
      headings = header.split(',')
      for icol in range(len(headings)):
        headings[icol] = headings[icol].replace('"','')
        print '%d: %s' % (icol+1,headings[icol])
      print 'Enter the column ids for the desired data, delimiting with spaces:'
      s = raw_input('-> ')
      w = s.split()
      columns = []
      columns.append(int(0))
      for i in range(len(w)):
        if not w[i].isdigit():
          print 'Entry %s not a recognized integer' % w[i]
        else:
          icol = int(w[i])-1
          if icol > 0:
            columns.append(int(w[i])-1)
      print 'Desired columns include:'
      for i in range(len(columns)):
        print '%d: %s' % (i+1,headings[columns[i]])
      yes = 0
      while(1):
        s = raw_input('Are these correct [y/n]?: ').lstrip().lower()
        if (s.startswith('y') and len(s) == 1) or \
           (s.startswith('yes') and len(s) == 3):
          yes = 1
          break
      if yes == 1:
        break
    
filename = 'processor_output.txt'
fout = open(filename,'w')

# write out filename header
for ifile in range(len(files)):
  for icol in range(len(columns)):
    if icol == 0:
      fout.write(filenames[ifile].strip())
    fout.write('\t')
  fout.write('\t')
fout.write('\n')
# write out variable header
for ifile in range(len(files)):
  for icol in range(len(columns)):
    fout.write(headings[columns[icol]])
    fout.write('\t')
  fout.write('\t')
fout.write('\n')
while(1):
  end_of_all_files_found = 1
  for ifile in range(len(files)):
    s = files[ifile].readline()
    w = s.split()
    if len(w) >= len(columns): # ensure that a valid line was read
      end_of_all_files_found = 0
      for icol in range(len(columns)): # skip time (only printed once)
        fout.write(w[columns[icol]])
        fout.write('\t')
    else:
      for icol in range(len(columns)):
        fout.write(' \t')
    fout.write('\t')
  fout.write('\n')
  if end_of_all_files_found == 1:
    break
      
for ifile in range(len(files)):
  files[ifile].close()
fout.close()

print 'done'
  
