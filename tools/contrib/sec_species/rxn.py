import sys
#from numpy import *
#import numpy as np
#import numpy
#import scipy
#import matplotlib
import datetime
import time

# NB: first run "find_values_no_temp_data.py" to generate text input files
# before setting this flag to true
exclude_species_without_temp_data = False

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y.%m.%d %H:%M:%S')

print('Author: PCL ',st)
print('=============================================')

"""

Author: P.C. Lichtner
Date: 10-07-2015

-----------------------------
database structures         |
-----------------------------
hanford.dat                 |
-----------------------------
basis           basis       |
aqueous         aqueous     |
minerals        gases       |
gases           minerals    |
srf cmplx       srf cmplex  |
-----------------------------

Input files: aq_spec.dat, gases.dat, minerals.dat
Skip files: aq_skip.dat, gas_skip.dat, min_skip.dat

Outfiles: chem.out

It is necessary to consider both standard (normal) and nonstandard ordering of the primary species.
Standard ordering is in the same order as the database; nonstandard implies swaping of
basis species (e.g. OH-/H+, CO2(aq)/HCO3-, Al(OH)4-/Al+++, etc.). If a database contains
reactions with both types of ordering then standard refers to form of the majority of reactions.

"""

#Example Lists of Primary Species: for MPHASE use CO2(aq) as primary species

#Enter list of primary species

#pri=['Al+++','Fe++','Mg++','Mn++','Na+','K+','Li+','H+','SiO2(aq)','Cl-','SO4--','O2(aq)','H2O']

"""
pri=[  
  'H+',
  'Al+++',
  'Ca++',
  'Mg++',
  'Fe++',
  'Li+',
  'Mn++',
  'UO2++',
  'VO++',
  'H2AsO4-',
  'K+',
  'Na+',
  'B(OH)3(aq)',
  'HCO3-',
  'Cl-',
  'F-',
  'HPO4--',
  'NH3(aq)',
  'SO4--',
  'SiO2(aq)',
  'O2(aq)',
  'Tracer',
  'Tracer2',
  'H2O'
]
"""
"""
pri=[
  'Ca++',
  'Mg++',
  'Fe++',
  'Na+',
  'K+',
  'Al+++',
  'Cl-',
  'SiO2(aq)',
  'O2(aq)',
  'H+',
  'CO2(aq)',
  'H2O'
  ]
"""

#pri=['H+','HCO3-','Na+','Cl-','H2O']
#pri=['H+','CO2(aq)','Na+','Cl-','H2O']
#pri=['OH-','HCO3-','Na+','Cl-','H2O']
#pri=['OH-','CO2(aq)','Na+','Cl-','H2O']

#pri=['H+','Al+++','HCO3-','Na+','Cl-','H2O']
#pri=['H+','Al(OH)4-','HCO3-','Na+','Cl-','H2O']
#pri=['OH-','Al(OH)4-','HCO3-','Na+','Cl-','H2O']

#pri=['OH-','Al(OH)4-','CO2(aq)','NaCl(aq)','Cl-','H2O']
#pri=['OH-','Al(OH)4-','CO2(aq)','NaCl(aq)','Cl-']

#pri=['H+','Al+++','H2O']
#pri=['OH-','Al(OH)4-','H2O']
#pri=['H+','Al(OH)4-','H2O']

#pri=['H+','H2O']
#pri=['OH-','H2O']

#pri=['Fe++','H+','SO4--','O2(aq)','H2O']

#pri=['CrO4--','O2(aq)','H+','H2O']

#pri=['Fe++','O2(aq)','H+','H2O']
#pri=['Fe+++','O2(aq)','H+','H2O']
#pri=['Fe++','Fe+++','H+','H2O']

#"""
pri=[
'H+',
'Na+',
'K+',
'Mg++',
'Al+++',
'CO2(aq)',
'Cl-',
'SiO2(aq)',
'H2O'
]
#"""

"""
pri=[
  'H+',
  'Na+',
  'Al+++',
  'CO2(aq)',
  'SiO2(aq)',
  'Cl-',
  'H2O'
  ]
"""

#pri=['Cd++','H+','NO3-','HCO3-','O2(aq)','H2O']
#pri=['UO2++','H+','Ca++','SO4--','HCO3-','O2(aq)','H2O']
#pri=['H+','Na+','Al+++','CO2(aq)','SiO2(aq)','Cl-','H2O']

#rare earth elements
"""
pri=[
'H2O','H+','HCO3-','HPO4--','SO4--','Cl-','F-','NO3-','O2(aq)',
'La+++','Sc+++','Y+++','Ce+++','Pr+++','Nd+++','Pm+++','Sm+++',
'Eu+++','Gd+++','Tb+++','Dy+++','Ho+++','Er+++','Tm+++','Yb+++','Lu+++'
]
"""

#=============== end of examples ==================

faq = open('aq_spec.dat','r')
fgas = open('gases.dat','r')
fmin = open('minerals.dat','r')
fchem = open('chem.out','w')

faq_skip = open('aq_skip.dat','r')
fgas_skip = open('gas_skip.dat','r')
fmin_skip = open('min_skip.dat','r')

if exclude_species_without_temp_data:
  fh = open('sec_species_no_temp_data.txt','r')
  sec_species_no_temp = sorted([x.strip() for x in fh.readlines()])
  fh.close()

  fh = open('gas_species_no_temp_data.txt','r')
  gas_species_no_temp = sorted([x.strip() for x in fh.readlines()])
  fh.close()

  fh = open('minerals_no_temp_data.txt','r')
  minerals_no_temp = sorted([x.strip() for x in fh.readlines()])
  fh.close()
 

#skip aqueous species listed in faq_skip
list_aq = {}
i=0
for record in faq:
  faq_skip.seek(0)

  w = record.strip().split(' ')
  name = w[0]
  skip = 0
  for line in faq_skip:
    line = line.rstrip('\n')
#   line = line[1:-1]
#   print 'line: ',name,line
    if name == line:
#     print 'line-found: ',name #line.rstrip('\n')
      skip = 1
      break
  if exclude_species_without_temp_data:
    if name in sec_species_no_temp:
      skip = 1
  if skip == 0:
#   print record
# list.append(record)
    i = i + 1
    list_aq[i] = record.rstrip('\n')
#   print list_aq[i]

#for i in list_aq:
#  print i,list_aq[i]

#check for duplicate primary species
for pri_spec in pri:
  ndup = 0
  ispec = 0
  for pri_spec1 in pri:
    ispec = ispec + 1
    if pri_spec == pri_spec1:
      ndup = ndup + 1
      if ndup > 1:
        print 'Error: duplicate primary species: ',ispec,pri_spec
        exit("Run Stopped")

#check for presence of H2O in list of primary species
for pri_spec in pri:
  if pri_spec == 'H2O':
    break
else:
  print 'Error: H2O missing from list of primary species!'
  exit("Run Stopped")

print 'PRIMARY_SPECIES'
fchem.write('PRIMARY_SPECIES\n')
npri=0
iflgo2 = 0
for pri_spec in pri:
  npri=npri+1

  if pri_spec == 'O2(aq)': # set O2 flag
    iflgo2 = 1

  print pri_spec
  fchem.write("  %s\n" % pri_spec) 
print 'END'
fchem.write("END\n") 

print 'SECONDARY_SPECIES'
fchem.write('SECONDARY_SPECIES\n')

#initialize dictionaries (probably could use numpy arrays instead)
iflg   = {}
sec    = {}
accept = {}

#read in reactions: first pass
nsec = 0 # counter for secondary species
for rxn in list_aq: # loop over reactions in database
  accept[rxn] = 0
  w = list_aq[rxn].strip().split(' ')
  name = w[0].strip()
  n = int(w[1].strip())

# loop over species in reaction: assume standard ordering
  for k in range(n):
    iflg[k] = 0 
    for pri_spec in pri:
      if pri_spec == w[2*(k+1)+1].strip("'"):
        iflg[k] = 1 # species found in list of primary species
        break
      if 'O2(g)' == w[2*(k+1)+1].strip("'"):
        if iflgo2 == 1:
          iflg[k] = 1

  indx = 1
  for k in range(n):
    if iflg[k] == 0:
      indx = 0 # reject reaction at first occurence of primary species not found
      break

  if indx == 1:
    nsec = nsec+1
    sec[nsec] = name.strip("'")
    accept[rxn] = 1 # accept reaction 

#   print 'pass1: ',rxn
#   print 'pass1: ',nsec,name.strip("'")
    print name.strip("'")
    fchem.write("  %s\n" % name.strip("'"))

# check for additional secondary species: assume nonstandard ordering
  iflg   = {}
  for pri_spec in pri: # check if primary species occurs as first species in reaction

    if pri_spec == name.strip("'"):

      kzero = 0

      for k in range(n):
        iflg[k] = 0
        for pri_spec1 in pri:
          if pri_spec1 == w[2*(k+1)+1].strip("'"):
            iflg[k] = 1

      kzero = 0
      for k in range(n):
        if iflg[k] == 0:
          kzero = kzero + 1

      for k in range(n):

        if iflg[k] == 0:
          reject = 0

          for new_sec in range(nsec):

            if w[2*(k+1)+1] == 'O2(g)':
              reject = 0
              break

            if sec[new_sec+1] == w[2*(k+1)+1].strip("'"):
              reject = 1
              break

          if reject == 0:

            nsec = nsec + 1
            accept[rxn] = 1
            if kzero == 1:
              if w[2*(k+1)+1].strip("'") == 'O2(g)':
                sec[nsec] = 'O2(aq)'.strip("'")
                iflgo2 = 1
              else:
                sec[nsec] = w[2*(k+1)+1].strip("'")
            elif kzero > 1:
              sec[nsec] = w[2*(k+1)+1].strip("'")

              if 'O2(g)' == w[2*(k+1)+1].strip("'"):
                nsec = nsec - 1 
                continue

#           print 'pass2: ',nsec,sec[nsec].strip("'")
            print sec[nsec].strip("'")
            fchem.write("  %s\n" % sec[nsec].strip("'"))

#read in reactions: second pass
#faq.seek(0)

irxn = 0
for rxn in list_aq:
  irxn = irxn + 1

  if accept[rxn] == 1: # skip reaction if already accepted
    continue

  w = list_aq[rxn].strip().split(' ')
  name = w[0].strip()
  n = int(w[1].strip())

# loop over species in reaction
  for k in range(n):
    iflg[k] = 0
    for pri_spec in pri:
      if pri_spec == w[2*(k+1)+1].strip("'"):
        iflg[k] = 1
        break

    for i in sec:
      if sec[i] == w[2*(k+1)+1].strip("'"):
        iflg[k] = 1
        break
  indx = 1
  for k in range(n):
    if iflg[k] == 0:
      indx = 0

  if indx == 1:
    reject = 0 # drop if all species in reaction already included: need only check first species
    for pri_spec in pri:
      if pri_spec == w[0].strip("'"):
        reject = 1
        break
    if (reject == 0):

#     if ('O2(aq)*' == name.strip("'")):
#       break
      nsec = nsec + 1
      sec[nsec] = name.strip("'")
#     print 'pass3: ',rxn
#     print 'pass3: ','nsec= ',nsec,name.strip("'")
      print name.strip("'")
      fchem.write("  %s\n" % name.strip("'"))

print 'END'
fchem.write("END\n") 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#skip gas species listed in fgas_skip
list_gas = {}
i=0
for record in fgas:
  fgas_skip.seek(0)

  w = record.strip().split(' ')
  name = w[0]
  skip = 0
  for line in fgas_skip:
    line = line.rstrip('\n')
    if name == line:
      skip = 1
      break
  if exclude_species_without_temp_data:
    if name in gas_species_no_temp:
      skip = 1
  if skip == 0:
    i = i + 1
    list_gas[i] = record.rstrip('\n')
#   print list_gas[i]

ngas=0
iflg = {}

fchem.write('GAS_SPECIES\n') 
print 'GAS_SPECIES'

#read in reactions
ngas = 0
for rxn in list_gas:
  w = list_gas[rxn].strip().split(' ')
  name = w[0].strip()
  n = int(w[2].strip())

# loop over species in reaction
  for k in range(n):
    iflg[k] = 0
    for pri_spec in pri:
      if pri_spec == w[2*(k+1)+2].strip("'"):
        iflg[k] = 1

    if 'O2(g)' == w[2*(k+1)+2].strip("'"):
        if iflgo2 == 1:
          iflg[k] = 1

    for i in sec:
      if sec[i] == w[2*(k+1)+2].strip("'"):
        iflg[k] = 1
        break

  indx = 1
  for k in range(n):
    if iflg[k] == 0:
      indx = 0

  if indx == 1:
    ngas=ngas+1
    print name.strip("'")
    fchem.write("  %s\n" % name.strip("'"))

print 'END'
fchem.write("END\n") 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#skip minerals listed in fmin_skip
list_min = {}
i=0
for record in fmin:
  fmin_skip.seek(0)

  w = record.strip().split(' ')
  name = w[0]
  skip = 0
  for line in fmin_skip:
    line = line.rstrip('\n')
    if name == line:
      skip = 1
      break
  if exclude_species_without_temp_data:
    if name in minerals_no_temp:
      skip = 1
  if skip == 0:
    i = i + 1
    list_min[i] = record.rstrip('\n')

nmin=0
iflg = {}

fchem.write("MINERALS\n") 
print 'MINERALS'

#read in reactions
nmin = 0
for rxn in list_min:
  w = list_min[rxn].strip().split(' ')
  name = w[0].strip()
  n = int(w[2].strip())

# loop over species in reaction
  for k in range(n):
    iflg[k] = 0
    for pri_spec in pri:
      if pri_spec == w[2*(k+1)+2].strip("'"):
        iflg[k] = 1

      if 'O2(g)' == w[2*(k+1)+2].strip("'"):
        if iflgo2 == 1:
          iflg[k] = 1

    for i in sec:
      if sec[i] == w[2*(k+1)+2].strip("'"):
        iflg[k] = 1
        break

  indx = 1
  for k in range(n):
    if iflg[k] == 0:
      indx = 0

  if indx == 1:
    nmin=nmin+1
    print name.strip("'")
    fchem.write("  %s\n" % name.strip("'"))

print 'END'
fchem.write("END\n") 
print('================================================')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print 'npri = ',npri,' nsec = ',nsec,' ngas = ',ngas,' nmin = ',nmin
print '\nFinished!'
fchem.write("\nnpri = %s\n" % npri)
fchem.write("nsec = %s\n" % nsec)
fchem.write("ngas = %s\n" % ngas)
fchem.write("nmin = %s\n" % nmin)


faq.close()
fgas.close()
fmin.close()
fchem.close()

faq_skip.close()
fgas_skip.close()
fmin_skip.close()
