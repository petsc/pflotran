# swap_comment_character.py
import sys
import shutil
import os
import fnmatch

dir1 = '1'
dir2 = '4'
suffix = '*.tec'

for root, dirnames, filenames in os.walk('1'):
  for filename in fnmatch.filter(filenames,suffix):
    print(filename)
    f1 = open(dir1+'/'+filename,'r')
    f2 = open(dir2+'/'+filename,'r')
    count = 0
    for line1 in f1:
      line2 = f2.readline()
      count += 1
      if line1 != line2:
        print('\nDifference on line %d ------------------\n' % count)
    f1.close()
    f2.close()
      
print('\nDone.')
