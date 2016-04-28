# water_saturation_pressure.py
import math
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

class WaterSaturationPressureEOS:
  def __init__(self):
    self.A = [-7.691234564,-2.608023696e1,-1.681706546e2, \
              6.423285504e1,-1.189646225e2,4.167117320, \
              2.097506760e1,1.e9,6.]
    self.TC = 0.
    self.one_m_tc = 0.
    self.one_m_tc_sq = 0.
    self.SC = 0.
    self.E1 = 0.
    self.E2 = 0.
    self.E2_bottom = 0.
    self.PCAP = 0.
  def getSaturationPressure(self,T):
    if T > 500.:
      print('Temperatures over 500 degrees C not supported')
      return
    self.getVariables(T)
    PS = self.PCAP*2.212e7
    return PS            
  def getSaturationPressureDerivative(self,T):
    if T > 500.:
      print('Temperatures over 500 degrees C not supported')
      return
    self.getVariables(T)
    dTC_dT = 1./647.3
    dSC_dTC = -self.A[0]-2.*self.A[1]*self.one_m_tc-3.*self.A[2]*self.one_m_tc_sq- \
              4.*self.A[3]*self.one_m_tc**3.-5.*self.A[4]*self.one_m_tc**4.
    dE1_dTC = (1.+self.A[5]*self.one_m_tc+self.A[6]*self.one_m_tc_sq)+ \
              self.TC*(-self.A[5]-2.*self.A[6]*self.one_m_tc)
    dE2_dTC = -1./self.E2_bottom+self.one_m_tc/(self.E2_bottom*self.E2_bottom)*2.*self.one_m_tc
    dPC_dTC = (-self.SC/(self.E1*self.E1)*dE1_dTC-dE2_dTC)*self.PCAP
    dPC_dSC = 1./self.E1*self.PCAP
    dPS_dT = (dPC_dSC*dSC_dTC+dPC_dTC)*dTC_dT*2.212e7
    return dPS_dT
  def getVariables(self,T):
    self.TC = (T+273.15)/647.3
    self.one_m_tc = 1.-self.TC
    self.one_m_tc_sq = self.one_m_tc*self.one_m_tc
    self.SC = self.A[0]*self.one_m_tc+self.A[1]*self.one_m_tc_sq+self.A[2]*self.one_m_tc**3.+ \
              self.A[3]*self.one_m_tc**4.+self.A[4]*self.one_m_tc**5.
    self.E1 = self.TC*(1.+self.A[5]*self.one_m_tc+self.A[6]*self.one_m_tc_sq)
    self.E2_bottom = self.A[7]*self.one_m_tc_sq+self.A[8]
    self.E2 = self.one_m_tc/self.E2_bottom
    self.PCAP = math.exp(self.SC/self.E1-self.E2)

test = WaterSaturationPressureEOS()
print(test.getSaturationPressure(25.))
print(test.getSaturationPressureDerivative(25.))

f = plt.figure(figsize=(8,6))
plt.subplot(1,1,1)
plt.xscale('log')
plt.yscale('log')

for t in range(0,101,10):
  string = ('%f' % t)
  permutations = []
  derivatives = []
  for i in range(-100,20,1):
    permutation = 10**(float(i/10))
    permutations.append(permutation)
    derivatives.append((test.getSaturationPressure(t+permutation)- \
                        test.getSaturationPressure(t)) / permutation)
  plt.plot(permutations,derivatives,label=string)

plt.legend(loc=2,title='Temperature')
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.14,right=.9)
plt.show()
            
