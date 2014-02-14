import sys
from numpy import *
import numpy as np

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

filename = 'pflotran.dat'
g = open(filename,'w')

filename = 'species.dat'
species = open(filename,'w')

filename = 'species_dict.dat'
spec = open(filename,'r')

t2p = {}
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
# print name,w[1],w[2],w[3]

  list = name+' '+w[1]+' '+w[2]+' '+w[3]
  g.writelines("%s\n" % list)
  list = name_tough+' ' + ':' + ' '+name
  species.writelines("%s\n" % list)

g.write("'null' 0 0 0\n")

#group 2: aqueous complexes

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
# print 'n=',n,logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8]

  name = t2p[name_tough]
  naq = naq + 1

  if n == 1:
#   print name,n,w[5],w[6],logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[4]+' '+w[5]+' '+t2p[w[6]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+\
    logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]
  elif n == 2:
#   print name,n,w[5],w[6],w[7],w[8],logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[4]+' '+w[5]+' '+t2p[w[6]]+' '+w[7]+' '+t2p[w[8]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+\
    logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]
  elif n == 3:
#   print name,n,w[5],w[6],w[7],w[8],w[9],w[10],\
    logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[4]+' '+w[5]+' '+t2p[w[6]]+' '+w[7]+' '+t2p[w[8]]+' '+w[9]+' '+t2p[w[10]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]

  elif n == 4:
#   print name,n,w[5],w[6],w[7],w[8],w[9],w[10],w[11],w[12],\
    logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[4]+' '+w[5]+' '+t2p[w[6]]+' '+w[7]+' '+t2p[w[8]]+' '+w[9]+' '+t2p[w[10]]+' '+w[11]+' '+t2p[w[12]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]

  elif n == 5:
#   print name,n,w[5],w[6],w[7],w[8],w[9],w[10],w[11],w[12],w[13],w[14],\
    logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[4]+' '+w[5]+' '+t2p[w[6]]+' '+w[7]+' '+t2p[w[8]]+' '+w[9]+' '+t2p[w[10]]+' '+w[11]+' '+t2p[w[12]]+' '+w[13]+' '+t2p[w[14]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]

  elif n == 6:
#   print name,n,w[5],w[6],w[7],w[8],w[9],w[10],w[11],w[12],w[13],w[14],w[15],w[16],\
    logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[4]+' '+w[5]+' '+t2p[w[6]]+' '+w[7]+' '+t2p[w[8]]+' '+w[9]+' '+t2p[w[10]]+' '+w[11]+' '+t2p[w[12]]+' '+\
    w[13]+' '+w[14]+' '+w[15]+' '+w[16]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]

  elif n == 7:
#   print name,n,w[5],w[6],w[7],w[8],w[9],w[10],w[11],w[12],w[13],w[14],w[15],w[16],w[17],w[18],\
    logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[4]+' '+w[5]+' '+t2p[w[6]]+' '+w[7]+' '+t2p[w[8]]+' '+w[9]+' '+t2p[w[10]]+' '+w[11]+' '+t2p[w[12]]+' '+\
    w[13]+' '+t2p[w[14]]+' '+w[15]+' '+t2p[w[16]]+' '+w[17]+' '+t2p[w[18]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]

  elif n == 8:
#   print name,n,w[5],w[6],w[7],w[8],w[9],w[10],\
    w[11],w[12],w[13],w[14],w[15],w[16],w[17],w[18],w[19],w[20],\
    logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[4]+' '+w[5]+' '+t2p[w[6]]+' '+w[7]+' '+t2p[w[8]]+' '+w[9]+' '+t2p[w[10]]+' '+\
    w[11]+' '+t2p[w[12]]+' '+w[13]+' '+t2p[w[14]]+' '+w[15]+' '+t2p[w[16]]+' '+w[17]+' '+t2p[w[18]]+' '+w[19]+' '+t2p[w[20]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]
  else:
    print 'n>8'

  g.writelines("%s\n" % list)
  list = name_tough+' ' + ':' + ' '+name
  species.write("%s\n" % list)

g.write("'null' 1 0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.\n")

#group 3: minerals

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
# print logk

  n = int(w[3])
  nn3 = max(nn3,n)
# print 'n=',n,logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8]

  name = t2p[name_tough]
  nmin = nmin + 1

  if n == 1:
#   print name,w[2],n,w[4],w[5],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 2:
#   print name,w[2],n,w[4],w[5],w[6],w[7],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 3:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8],\
    w[2],w[3],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+\
    w[2]+' '+w[3]+' '+w[1]
  elif n == 4:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],w[10],w[11],\
    logk[1],logk[2],logk[3],logk[4],\
    logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+w[10]+' '+t2p[w[11]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+\
    logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 5:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    w[10],w[11],w[12],w[13],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    w[10]+' '+t2p[w[11]]+' '+w[12]+' '+t2p[w[13]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 6:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    w[10],w[11],w[12],w[13],w[14],w[15],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    w[10]+' '+t2p[w[11]]+' '+w[12]+' '+t2p[w[13]]+' '+w[14]+' '+t2p[w[15]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 7:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    w[10],w[11],w[12],w[13],w[14],w[15],w[16],w[17],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    w[10]+' '+t2p[w[11]]+' '+w[12]+' '+t2p[w[13]]+' '+w[14]+' '+t2p[w[15]]+' '+w[16]+' '+t2p[w[17]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 8:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    w[10],w[11],w[12],w[13],w[14],w[15],w[16],w[17],w[18],w[19],\
    logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    w[10]+' '+t2p[w[11]]+' '+w[12]+' '+t2p[w[13]]+' '+w[14]+' '+t2p[w[15]]+' '+w[16]+' '+t2p[w[17]]+' '+w[18]+' '+t2p[w[19]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  else:
    print 'n>8'
    stop

  g.writelines("%s\n" % list)
  list = name_tough+' ' + ':' + ' '+name
  species.write("%s\n" % list)

g.write("'null'  0. 1 0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0.\n")

#group 4: gases

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
# print logk

  n = int(w[3])
  nn4 = max(nn4,n)
# print 'n=',n,logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8]
  logkk=logk[1],logk[2],logk[3],logk[4],logk[5],logk[6],logk[7],logk[8]

  name = t2p[name_tough]
  ngas = ngas + 1

  if n == 1:
#   print name,w[2],n,w[4],w[5],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 2:
#   print name,w[2],n,w[4],w[5],w[6],w[7],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 3:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[2],w[3],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[2]+' '+w[3]+' '+w[1]
  elif n == 4:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],w[10],w[11],\
    logk[1],logk[2],logk[3],logk[4],\
    logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+w[10]+' '+t2p[w[11]]+' '+\
    logk[1]+' '+logk[2]+' '+logk[3]+' '+logk[4]+' '+\
    logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 5:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    w[10],w[11],w[12],w[13],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    w[10]+' '+t2p[w[11]]+' '+w[12]+' '+t2p[w[13]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 6:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    w[10],w[11],w[12],w[13],w[14],w[15],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    w[10]+' '+t2p[w[11]]+' '+w[12]+' '+t2p[w[13]]+' '+w[14]+' '+t2p[w[15]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 7:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    w[10],w[11],w[12],w[13],w[14],w[15],w[16],w[17],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    w[10]+' '+t2p[w[11]]+' '+w[12]+' '+t2p[w[13]]+' '+w[14]+' '+t2p[w[15]]+' '+w[16]+' '+t2p[w[17]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  elif n == 8:
#   print name,w[2],n,w[4],w[5],w[6],w[7],w[8],w[9],\
    w[10],w[11],w[12],w[13],w[14],w[15],w[16],w[17],w[18],w[19],logk[1],logk[2],logk[3],\
    logk[4],logk[5],logk[6],logk[7],logk[8],w[1]

    list = name+' '+w[2]+' '+w[3]+' '+w[4]+' '+t2p[w[5]]+' '+w[6]+' '+t2p[w[7]]+' '+w[8]+' '+t2p[w[9]]+' '+\
    w[10]+' '+t2p[w[11]]+' '+w[12]+' '+t2p[w[13]]+' '+w[14]+' '+t2p[w[15]]+' '+w[16]+' '+t2p[w[17]]+' '+w[18]+' '+t2p[w[19]]+' '+logk[1]+' '+logk[2]+' '+logk[3]+' '+\
    logk[4]+' '+logk[5]+' '+logk[6]+' '+logk[7]+' '+logk[8]+' '+w[1]
  else:
    print 'n>8'
    stop

  g.writelines("%s\n" % list)
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
