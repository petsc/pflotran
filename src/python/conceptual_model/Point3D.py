import math
# Point3D -- Defines a point within space
class Point3D:
  def __init__(self,x,y,z):
    self.x = x
    self.y = y 
    self.z = z
    self.name = ""
  def setName(self,name):
    self.name = name
  def setX(self,x):
    self.x = x
  def setY(self,y):
    self.y = y
  def setZ(self,z):
    self.z = z
  def getName(self):
    return self.name
  def sameName(self,name):
    if self.name.strip().startswith(name.strip()) \
       and self.name.strip().endswith(name.strip()):
#      print name, self.name, 'yes'
      return 1
    else:
#      print name, self.name, 'no'
      return 0
  def getX(self):
    return self.x
  def getY(self):
    return self.y
  def getZ(self):
    return self.z
  def translate(self,dx,dy,dz):
    self.x += dx
    self.y += dy
    self.z += dz
  def rotate(self,rotation):
    rotation_radians = rotation/180*math.pi
    new_x = math.cos(rotation_radians)*self.x-math.sin(rotation_radians)*self.y
    new_y = math.sin(rotation_radians)*self.x+math.cos(rotation_radians)*self.y
    self.x = new_x
    self.y = new_y
  def translate_and_rotate(self,point,rotation):
    self.translate(point.x,point.y,0.)
    self.rotate(rotation)
  def rotate_around_point(self,point,rotation):
    self.translate(-point.x,-point.y,0.)
    self.rotate(rotation)
    self.translate(point.x,point.y,0.)
  def printName(self):
    print self.name
  def printInfo(self):
    print self.name, self.x, self.y, self.z