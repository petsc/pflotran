# pflotran.py
import math
import sys

class AnalyticalSolution:
  def __init__(self,initial_concentration,final_concentration,
               Darcy_velocity,diffusion_coefficient,dispersivity,
               tortuosity,porosity,retardation,half_life):
    self.c0 = initial_concentration
    self.c1 = final_concentration
    self.R = retardation
    self.lam = -1.*math.log(0.5)/half_life
    self.porosity = porosity
    self.U = Darcy_velocity
    self.D = dispersivity*Darcy_velocity + \
             porosity*tortuosity*diffusion_coefficient
    print(self.D)
  def ogata_banks(self,x,t):
    # Based on "Solution of the Differential Equation of Longitudinal
    # Dispersion in Porous Media", USGS Professional Paper 411-A by
    # Akio Ogata and R.B. Banks, 1961.
    D_ = self.D/self.porosity
    v = self.U/self.porosity
    temp = 0.5* \
           (math.erfc((x-v/self.R*t)/(2.*math.sqrt(D_/self.R*t))) + \
            math.exp(v*x/D_) * \
            math.erfc((x+v/self.R*t)/(2.*math.sqrt(D_/self.R*t))))
    value = temp*(self.c1-self.c0) + self.c0
    return value
  def de_marsily_no_reaction(self,x,t):
    # Based on Equation 10.3.2 in Quantitative Hydrogeology, Ghislain de
    # Marsily, 1986.
    D_ = self.D
    U_over_porR = self.U/(self.porosity*self.R)
    two_sqrt_Dt_over_porR = 2.*math.sqrt(D_*t/(self.porosity*self.R))
    temp = 0.5* \
           (math.erfc((x-t*U_over_porR)/two_sqrt_Dt_over_porR) + \
            math.exp(self.U*x/D_) * \
            math.erfc((x+t*U_over_porR)/two_sqrt_Dt_over_porR))
    value = temp*(self.c1-self.c0) + self.c0
    return value
  def de_marsily(self,x,t):
    # Based on Equation 10.3.4 in Quantitative Hydrogeology, Ghislain de
    # Marsily, 1986.
    D_ = self.D
    Ux_over_2D = self.U/(2.*D_)
    temp = self.U/(self.porosity*self.R)
    sqrt_term = math.sqrt(temp*temp + 4.*self.lam*D_ / \
                                        (self.porosity*self.R))
    two_sqrt_Dt_over_porR = 2.*math.sqrt(D_*t/(self.porosity*self.R))
    beta = math.sqrt(Ux_over_2D*Ux_over_2D + \
                     self.lam*self.porosity*self.R/D_)
    temp = 0.5*math.exp(self.U*x/(2.*D_)) * \
               (math.exp(-1.*beta*x) * \
                math.erfc((x-t*sqrt_term)/two_sqrt_Dt_over_porR) + \
                math.exp(beta*x) * \
                math.erfc((x+t*sqrt_term)/two_sqrt_Dt_over_porR))
    value = temp*(self.c1-self.c0) + self.c0
    return value
  def bear(self,x,t):
    # Based on Equation 10.6.22 in "Dynamics of Fluids in Porous Media", 
    # Jacob Bear, 1988.
    D_ = self.D/self.porosity
    v = self.U/self.porosity
    beta = math.sqrt(v*v/(4.*D_*D_)+ self.lam/D_)
    temp = 0.5 * \
           math.exp(v*x/(2.*D_)) * \
           (math.exp(-1.*beta*x) * \
            math.erfc((x-math.sqrt(v*v+4.*self.lam*D_)*t) / \
                      (2.*math.sqrt(D_*t))) + \
            math.exp(beta*x) * \
            math.erfc((x+math.sqrt(v*v+4.*self.lam*D_)*t) / \
                      (2.*math.sqrt(D_*t))))
    value = temp*(self.c1-self.c0) + self.c0
    return value

