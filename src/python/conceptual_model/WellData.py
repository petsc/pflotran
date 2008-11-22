import sys
import math
from datetime import *
from Point3D import *

# Point3D -- Defines a point within space
class WellData:
  def __init__(self,filename):
    self.filename = filename
    f = open(filename,'r')
    self.f = f
    s = self.f.readline()
    w = s.split("\t")
    self.well_name = w[1]
    s = self.f.readline()
    w = s.split("\t")
    if len(w) < 4:
      self.location = Point3D(float(w[1]),float(w[2]),0.)
    else:
      self.location = Point3D(float(w[1]),float(w[2]),float(w[3]))
    # remove date/time header
    s = self.f.readline()
    # get bounds on date/time
    self.startdate = self.readDate()
    while 1:
      prev_line = s
      s = self.f.readline().strip()
      if len(s) == 0:
        break
    self.enddate = readDateFromString(prev_line)
    self.rewind()
    self.time1 = self.startdate
    self.time2 = self.startdate
    self.value1 = 0.
    self.value2 = 0.
    self.readDataValue()
    self.incrementDataValue()
  def getValueAtTime(self,time):
    if time < self.time1:
      self.rewind()
    while 1:
      if self.time1 <= time and self.time2 >= time:
        break
      self.incrementDataValue()
    delta1 = time-self.time1
    delta2 = self.time2-self.time1
    weight = float(delta1.days*3600*24+delta1.seconds)/ \
             float(delta2.days*3600*24+delta2.seconds)
    value = (1.-weight)*self.value1+weight*self.value2
    return value
  def incrementDataValue(self):
    self.time1 = self.time2
    self.value1 = self.value2
    self.readDataValue()
  def readDataValue(self):
    s = self.f.readline()
    if len(s.strip()) < 2:
      print 'Error blank line read from file: ', self.filename
      quit()
    w = s.split("\t")
    self.time2 = readDateFromString(s)
    self.value2 = float(w[1])
  def setLocation(self,x,y,z):
    self.location = Point3D(x,y,z)
  def isWithinTimeBounds(self,date_time):
    return date_time > self.startdate and date_time < self.enddate
  def rewind(self):
    self.f.seek(0)
    # skip headers
    s = self.f.readline()
    s = self.f.readline()
    s = self.f.readline()
  def readDate(self):
    s = self.f.readline()
    return readDateFromString(s)
  def printInfo(self):
    print 'Well Name: ', self.well_name
    print '  Location: ', self.location.x, self.location.y, self.location.z
    print '  Start Date: ', self.startdate.ctime() 
    print '  End Date: ', self.enddate.ctime()
  def finalize(self):
    self.f.close()
def readDateFromString(s):
  w = s.split("\t")
  w = w[0].split("/")
  month = int(w[0])
  day = int(w[1])
  w = w[2].split()
  year = int(w[0])
  if year < 50:
    year += 2000
  elif year < 100:
    year += 1900
  hour = 0
  minute = 0
  if len(w) > 1:
    w2 = w[1].split(":")
    hour = int(w2[0])
    if len(w) > 2:
      w3 = w[2].strip().upper()
      if w3.startswith("AM"):
        if hour == 12:
            hour = 0
      elif w3.startswith("PM"):
        if hour < 12:
          hour += 12
    minute = int(w2[1].strip())
  date_time = datetime(year,month,day,hour,minute)
  return date_time
