# Jacobian_diff.py
import sys

filename1 = 'Gjacobian-1.out'
file1 = open(filename1,'r')

filename2 = 'numerical_jacobian-1.out'
file2 = open(filename2,'r')

file1.readline()
file1.readline()
file2.readline()
file2.readline()

row_count = 0
for line1 in file1:
  line2 = file2.readline()
  row1 = line1.split(':')
  row2 = line2.split(':')
  row1 = row1[1].strip(' ()\n').split(')  (')
  row2 = row2[1].strip(' ()\n').split(')  (')
  icount1 = 0
  icount2 = 0
  increment1 = True
  increment2 = True
  while icount1 < len(row1) and icount2 < len(row2):
    if increment1:
      w1 = row1[icount1].split(',')
      i1 = int(w1[0])
    if increment2:
      w2 = row2[icount2].split(',')
      i2 = int(w2[0])
    if i1 < i2:
      increment1 = True
      increment2 = False
      icount1 += 1
      value1 = float(w1[1])
      if value1 > 1.e-20:
        string = 'Missing value (%e) in "%s" at row %d column %d\n' % \
                 (value1,filename2,row_count,i1)
        print(string)
    elif i1 > i2:
      increment1 = False
      increment2 = True
      icount2 += 1
      value2 = float(w2[1])
      if value2 > 1.e-20:
        string = 'Missing value (%e) in "%s" at row %d column %d\n' % \
                 (value2,filename1,row_count,i2)
        print(string)
    else:
      increment1 = True
      increment2 = True
      value1 = float(w1[1])
      value2 = float(w2[1])
      if value1 > 0. and (value1-value2)/value1 > 1.e-3:
        string = 'Large difference on row %d, column %d\n' % \
                 (row_count,i1)
        print(string)
      icount1 += 1
      icount2 += 1
  row_count += 1

print('done')
        

