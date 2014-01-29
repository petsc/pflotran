# pflotran.py
import math
import sys

class AnalyticalSolution:
  def __init__(self,initial_concentration,final_concentration,
               Darcy_velocity,dispersion_coefficient,retardation,
               decay_rate,porosity):
    self.c0 = initial_concentration
    self.c1 = final_concentration
    self.U = Darcy_velocity
    self.D = dispersion_coefficient
    self.R = retardation
    self.lam = decay_rate
    self.porosity = porosity
  def ogata_banks(self,x,t):
    v = self.U/self.porosity
    temp = 0.5* \
           (math.erfc((x-v/self.R*t)/(2.*math.sqrt(self.D/self.R*t))) + \
            math.exp(v*x/self.D) * \
            math.erfc((x+v/self.R*t)/(2.*math.sqrt(self.D/self.R*t))))
    value = temp*(self.c1-self.c0) + self.c0
    return value
  def de_marsily_no_reaction(self,x,t):
    U_over_porR = self.U/(self.porosity*self.R)
    two_sqrt_Dt_over_porR = 2.*math.sqrt(self.D*t/(self.porosity*self.R))
    temp = 0.5* \
           (math.erfc((x-t*U_over_porR)/two_sqrt_Dt_over_porR) + \
            math.exp(self.U*x/self.D) * \
            math.erfc((x+t*U_over_porR)/two_sqrt_Dt_over_porR))
    value = temp*(self.c1-self.c0) + self.c0
    return value
  def de_marsily(self,x,t):
    Ux_over_2D = self.U/(2.*self.D)
    temp = self.U/(self.porosity*self.R)
    sqrt_term = math.sqrt(temp*temp + 4.*self.lam*self.D / \
                                        (self.porosity*self.R))
    two_sqrt_Dt_over_porR = 2.*math.sqrt(self.D*t/(self.porosity*self.R))
    beta = math.sqrt(Ux_over_2D*Ux_over_2D + \
                     self.lam*self.porosity*self.R/self.D)
    temp = 0.5*math.exp(self.U*x/(2.*self.D)) * \
               (math.exp(-1.*beta*x) * \
                math.erfc((x-t*sqrt_term)/two_sqrt_Dt_over_porR) + \
                math.exp(beta*x) * \
                math.erfc((x+t*sqrt_term)/two_sqrt_Dt_over_porR))
    value = temp*(self.c1-self.c0) + self.c0
    return value
  # incorrect implementation in NUTS
  def nuts(self,x,t):
    v = self.U/self.porosity
    beta = math.sqrt(v*v/(4.*self.D*self.D))+ self.lam/self.D
    temp = 0.5 * \
           math.exp(v*x/(2.*self.D)) * \
           (math.exp(-1.*beta*x) * \
            math.erfc((x-math.sqrt(v*v+4.*self.lam*self.D)*t) / \
                      (2.*math.sqrt(self.D*t))) + \
            math.exp(beta*x) * \
            math.erfc((x+math.sqrt(v*v+4.*self.lam*self.D)*t) / \
                      (2.*math.sqrt(self.D*t))))
    value = temp*(self.c1-self.c0) + self.c0
    return value
