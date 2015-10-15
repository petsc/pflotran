import sys
import os
import matplotlib
import matplotlib.pyplot as plt

f = plt.figure(figsize=(10,10))
f.suptitle("PFLOTRAN Build Speedup with Intel on Linux Mint",fontsize=16)

plt.xlim(0.9,9)
plt.ylim(40,500)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('# Cores (make -j # pflotran)')
plt.ylabel('Time [s]')
cores = [1, 2, 4, 8]
times = [437,228,157,136]
ideal = [400,200,100,50]
plt.plot(cores,times,marker='o',markersize=10)
plt.plot(cores,ideal)

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.12,right=.9)

plt.show()
