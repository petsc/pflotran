import sys
import os
import matplotlib
import matplotlib.pyplot as plt

f = plt.figure(figsize=(10,10))
f.suptitle("PFLOTRAN Build Speedup",fontsize=16)

plt.xlim(0.9,10)
plt.ylim(30,500)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('# Cores (make -j # pflotran)')
plt.ylabel('Time [s]')
cores = [1, 2, 4, 8]
intel_O_mint = [437,228,157,136]
gnu_O_mint = [188,95,62,50]
ideal = [300,150,75,37.5]
plt.plot(cores,intel_O_mint,marker='o',markersize=10,label='Intel Opt Mint')
plt.plot(cores,gnu_O_mint,marker='o',markersize=10,label='GNU Opt Mint')
plt.plot(cores,ideal,label='Ideal')

plt.legend(loc=1,title='Platform')
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.12,right=.9)

plt.show()
