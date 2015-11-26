import sys
from numpy import *
import numpy as np

import datetime
import time

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

print(st)
print('======================================')

"""

Author: P.C. Lichtner
Date: 2-12-2014

-----------------------------
database structures         |
-----------------------------
thermXu5.dat    hanford.dat |
-----------------------------
basis           basis       |
aqueous         aqueous     |
minerals        gases       |
gases           minerals    |
srf cmplx       srf cmplex  |
-----------------------------

Note: it is necessary to copy the gaseous species block by hand 
above the mineral block. Up to 8 primary species in a reaction
is allowed.

Files needed: convert_tough.py (this file), database thermXu5.dat, 
and dictionary: species_dict.dat

Outfiles: pflotran.dat and species.dat (not needed except 
for checking species names)

"""

print 'start'

filename = 'thermXu5.dat'
f = open(filename,'r')

filename = 'pflotran_rev.dat'
g = open(filename,'w')

filename = 'species_rev.dat'
species = open(filename,'w')

filename = 'species_dict.dat'
spec = open(filename,'r')

t2p = {} # initialize dictionary
for line in spec:
  w = line.strip().split(':')
  t2p[w[0].strip()] = w[1].strip()

#skip first 11 lines of database
for i in range(1,11):
  s = f.readline()

#print temperatures
s = f.readline()
print s.title()
g.write(s.title())

name_exit = "'null'"

#group 1: primary species

g.write(":basis_species  a0  valence  formula weight [g]\n")

print '-> Processing primary species'

npri = 0
for i in range(1,100):
  s = f.readline()
  w = s.split()
  name_tough = w[0]
  if name_tough == name_exit:
  	break

  name = t2p[name_tough]
  npri = npri + 1

  list = name+' '+w[1]+' '+w[2]+' '+w[3]
  g.writelines("%s\n" % list)
  list = name_tough+' ' + ':' + ' '+name
  species.writelines("%s\n" % list)

g.write("'null' 0 0 0\n")

#group 2: aqueous complexes

g.write(":species name  num (n_i A_i, i=1,num)  log K (1:8)  a0  valence  formula weight [g]\n")

print '-> Processing aqueous species'

nn2 = 1 # initialize species counter

naq = 0
for i in range(1,1000):
  s1 = f.readline()

  w = s1.split()
  name_tough = w[0]
  if name_tough == name_exit:
  	  break

  s2 = f.readline()
  s3 = f.readline()
  logk = s2.split()

  n = int(w[4])
  nn2 = max(nn2,n)

  name = t2p[name_tough]
  naq = naq + 1

  stoi = name+' '+w[4]
  for k in range(n):
    stoi = stoi +' '+w[2*k+5]+' '+t2p[w[2*k+6]]

  list_logk = ' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+\
  logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]

  list_ptr = stoi+list_logk
  g.writelines("%s\n" % list_ptr)
  list = name_tough+' ' + ':' + ' '+name
  species.write("%s\n" % list)

g.write("'null' 1 0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.\n")

#group 3: minerals

g.write(":mineral  name molar_vol  num (n_i A_i, i=1,num) log K (1:8)  formula weight [g]\n")

print '-> Processing minerals'

nn3 = 1
nmin = 0
for i in range(1,1000):
  s1 = f.readline()
  w = s1.split()
  name_tough = w[0]
  if name_tough == name_exit:
  	  break
  s2 = f.readline()
  s3 = f.readline()
  logk = s2.split()

  n = int(w[3])
  nn3 = max(nn3,n)

  name = t2p[name_tough]
  nmin = nmin + 1

  stoi = name+' '+w[2]+' '+w[3]
  for k in range(n):
    stoi = stoi +' '+w[2*k+4]+' '+t2p[w[2*k+5]]

  list_logk = ' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+\
  logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]

  list_ptr = stoi+list_logk
  g.writelines("%s\n" % list_ptr)
  list = name_tough+' ' + ':' + ' '+name
  species.write("%s\n" % list)

g.write("'null'  0. 1 0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0.\n")

#group 4: gases

g.write(":gas  name molar_vol  num (n_i A_i, i=1,num) log K (1:8)  formula weight [g]\n")

print '-> Processing gases'

nn4 = 1
ngas = 0
for i in range(1,1000):
  s1 = f.readline()
  w = s1.split()
  name_tough = w[0]
  if name_tough == name_exit:
  	  break
  s2 = f.readline()
  s3 = f.readline()
  logk = s2.split()

  n = int(w[3])
  nn4 = max(nn4,n)

  name = t2p[name_tough]
  ngas = ngas + 1

  stoi = name+' '+w[2]+' '+w[3]
  for k in range(n):
    stoi = stoi +' '+w[2*k+4]+' '+t2p[w[2*k+5]]

  list_logk = ' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+\
    logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]

  list_ptr = stoi+list_logk
  g.writelines("%s\n" % list_ptr)
  list = name_tough+' ' + ':' + ' '+name
  species.write("%s\n" % list)

g.write("'null' 0. 1 1. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0.\n")

string='maximum number of primary species in rxns: '
print string.title(),nn2,nn3,nn4
print 'npri = ',npri,' naq = ',naq,' nmin = ',nmin,' ngas = ',ngas

f.close()
g.close()
species.close()
spec.close()

print 'end'
