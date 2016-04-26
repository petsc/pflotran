# Plane -- Defines a plain in space based on 3 points
class Plane:
  def __init__(self,pt1,pt2,pt3):
    self.point1 = pt1
    self.point2 = pt2 
    self.point3 = pt3
    self.A = 0
    self.B = 0
    self.C = 0
    self.D = 0
  def setZ(self,ipoint,z):
    if ipoint == 1:
      self.point1.setZ(z)
    elif ipoint == 2:
      self.point2.setZ(z)
    elif ipoint == 3:
      self.point3.setZ(z)
  def getZ(self,ipoint):
    if ipoint == 1:
      return self.point1.getZ()
    elif ipoint == 2:
      return self.point2.getZ()
    elif ipoint == 3:
      return self.point3.getZ()
    else: 
      return -999.
  def getPoint(self,name):
    if self.point1.sameName(name) == 1:
      return 1
    elif self.point2.sameName(name) == 1:
      return 2
    elif self.point3.sameName(name) == 1:
      return 3
    else:
      return -1
  def computePlaneCoefficients(self):
    x1 = self.point1.getX()
    y1 = self.point1.getY()
    z1 = self.point1.getZ()
    x2 = self.point2.getX()
    y2 = self.point2.getY()
    z2 = self.point2.getZ()
    x3 = self.point3.getX()
    y3 = self.point3.getY()
    z3 = self.point3.getZ()
#    print x1,y1,z1
#    print x2,y2,z2
#    print x3,y3,z3
    self.A = y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2)
    self.B = z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2)
    self.C = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)
    self.D = -1.*(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1))
  def printCoefficients(self):
    print self.A, self.B, self.C, self.D
  def computeX(self,y,z):
    # X = (-D-By-Cz)/A
    return -1.*(self.D+self.B*t+self.C*z)/self.A
  def computeY(self,x,z):
    # Y = (-D-Ax-Cz)/B
    return (-self.D-self.A*x-self.C*z)/self.B
  def computeZ(self,x,y):
    # Z = (-D-Ax-By)/C
    return (-self.D-self.A*x-self.B*y)/self.C
  def computeGradientZ_X(self,x,y):
    # dZ_dx = -A/C
    return -self.A/self.C
  def computeGradientZ_Y(self,x,y):
    # dZ_dy = -B/C
    return -self.B/self.C
  def printPointNames(self):
    self.point1.printName()
    self.point2.printName()
    self.point3.printName()
    
