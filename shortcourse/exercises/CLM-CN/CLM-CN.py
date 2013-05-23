import sys
import os
try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
  sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft

filename = 'CLM-CN-obs-0.tec'

f = plt.figure(figsize=(14,10))
f.suptitle("CLM-CN",fontsize=16)
plt.subplot(2,2,1)
plt.xlabel('Time [y]')
plt.ylabel('Concentration [mol/m^3]')

#plt.xlim(0.,1.)
#plt.ylim(4.8,8.2)
#plt.yscale('log')
#plt.grid(True)

data = pft.Dataset(filename,1,4)

# 12 g C / 1 g N
CN_ratio_12 = 12./12. / (1. / 14.)
# 10 g C / 1 g N
CN_ratio_10 = 10./12. / (1. / 14.)

# if is a dummy array 
new_array = pft.Dataset(filename,1,2).get_array('y')

arrayN = pft.Dataset(filename,1,3).get_array('y')
arrayC = pft.Dataset(filename,1,4).get_array('y')
arraySOM1 = pft.Dataset(filename,1,5).get_array('y')
arraySOM2 = pft.Dataset(filename,1,6).get_array('y')
arraySOM3 = pft.Dataset(filename,1,7).get_array('y')
arraySOM4 = pft.Dataset(filename,1,8).get_array('y')
arrayLabileC = pft.Dataset(filename,1,9).get_array('y')
arrayCelluloseC = pft.Dataset(filename,1,10).get_array('y')
arrayLigninC = pft.Dataset(filename,1,11).get_array('y')
arrayLabileN = pft.Dataset(filename,1,12).get_array('y')
arrayCelluloseN = pft.Dataset(filename,1,13).get_array('y')
arrayLigninN = pft.Dataset(filename,1,14).get_array('y')

for i in range(len(new_array)):
  new_array[i] = arrayLabileC[i] + arrayLabileN[i]
plt.plot(data.get_array('x'),new_array,label='Labile')
for i in range(len(new_array)):
  new_array[i] = arrayCelluloseC[i] + arrayCelluloseN[i]
plt.plot(data.get_array('x'),new_array,label='Cellulose')
for i in range(len(new_array)):
  new_array[i] = arrayLigninC[i] + arrayLigninN[i]
plt.plot(data.get_array('x'),new_array,label='Lignin')

#'best'         : 0, (only implemented for axis legends)
#'upper right'  : 1,
#'upper left'   : 2,
#'lower left'   : 3,
#'lower right'  : 4,
#'right'        : 5,
#'center left'  : 6,
#'center right' : 7,
#'lower center' : 8,
#'upper center' : 9,
#'center'       : 10,
plt.legend(loc=1)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

# SOM
plt.subplot(2,2,2)
plt.xlabel('Time [y]')
plt.ylabel('Concentration [mol/m^3]')
plt.plot(data.get_array('x'),arraySOM1,label='SOM1')
plt.plot(data.get_array('x'),arraySOM2,label='SOM2')
plt.plot(data.get_array('x'),arraySOM3,label='SOM3')
plt.plot(data.get_array('x'),arraySOM4,label='SOM4')
plt.legend(loc=1)
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))

# C
plt.subplot(2,2,3)
plt.xlabel('Time [y]')
plt.ylabel('Concentration [mol/m^3]')
for i in range(len(arrayC)):
  new_array[i] = arrayC[i] + arrayLabileC[i]  + \
                 arrayCelluloseC[i] + arrayLigninC[i] + \
                 (arraySOM1[i] + arraySOM2[i] + \
                  arraySOM3[i] + arraySOM4[i])
plt.plot(data.get_array('x'),new_array,label='Total C')
for i in range(len(arrayC)):
  new_array[i] = arrayLabileC[i] + \
                 arrayCelluloseC[i] + arrayLigninC[i] + \
                 (arraySOM1[i] + arraySOM2[i] + \
                  arraySOM3[i] + arraySOM4[i])
plt.plot(data.get_array('x'),new_array,label='Immobile C')
plt.plot(data.get_array('x'),arrayC,label='Mineral C')
plt.legend(loc=1)
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))

# N
plt.subplot(2,2,4)
plt.xlabel('Time [y]')
plt.ylabel('Concentration [mol/m^3]')
for i in range(len(arrayN)):
  new_array[i] = arrayN[i] + arrayLabileN[i] + \
                 arrayCelluloseN[i] + arrayLigninN[i] + \
                 (arraySOM1[i] + arraySOM2[i]) / CN_ratio_12 + \
                 (arraySOM3[i] + arraySOM4[i]) / CN_ratio_10
plt.plot(data.get_array('x'),new_array,label='Total N')
for i in range(len(arrayN)):
  new_array[i] = arrayLabileN[i] + \
                 arrayCelluloseN[i] + arrayLigninN[i] + \
                 (arraySOM1[i] + arraySOM2[i]) / CN_ratio_12 + \
                 (arraySOM3[i] + arraySOM4[i]) / CN_ratio_10
plt.plot(data.get_array('x'),new_array,label='Immobile N')
#plt.plot(data.get_array('x'),arrayLabileN,label='Labile N')
#plt.plot(data.get_array('x'),arraySOM1*CN_fraction_12_N,label='SOM1 N')
plt.plot(data.get_array('x'),arrayN,label='Mineral N')
plt.legend(loc=1)
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.12,right=.9)

plt.show()
