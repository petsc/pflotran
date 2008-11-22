import sys
import math
from WellData import *
from Plane import *
from datetime import *

new_origin = Point3D(0.,0.,0.)
translation_angle = -14.

print 'Enter input file name:'
filename = sys.stdin.readline().strip()
f = open(filename,'r')
s = f.readline()
w = s.split(':')
w = w[1].split()
origin = Point3D(float(w[0]),float(w[1]),0.)
delta = Point3D(new_origin.x-origin.x,new_origin.y-origin.y,0.)
origin.translate_and_rotate(delta,translation_angle)
origin.setName('Origin')
origin.printInfo()


s = f.readline()
w = s.split(':')
w = w[1].split()
datum = Point3D(float(w[0]),float(w[1]),0.)
datum.translate_and_rotate(delta,translation_angle)
datum.setName('Datum')
datum.printInfo()

well = []
for i in range(3):
  s = f.readline()
  w = s.split(":")
  filename = w[1].strip()
  well.append(WellData(filename))
  well[i].location.translate_and_rotate(delta,translation_angle)
  well[i].location.setName(w[0].strip())
  well[i].printInfo()

well_plane = Plane(well[0].location,well[1].location,well[2].location)

s = f.readline()
w = s.split(": ")
start_date = readDateFromString(w[1])

s = f.readline()
w = s.split(": ")
end_date = readDateFromString(w[1])

time = start_date
time_increment = timedelta(seconds=3600)
end_time = end_date

for i in range(3):
  if not well[i].isWithinTimeBounds(start_date):
    print 'Start date before beginning of well data:'
    well[i].printInfo()
    quit()
  if not well[i].isWithinTimeBounds(end_date):
    print 'Start date after end of well data:'
    well[i].printInfo()
    quit()
    
s = f.readline()
w = s.split(":")
filename_base = w[1].strip()
datum_file_out = open(filename_base+".datum",'w')
gradient_file_out = open(filename_base+".gradient",'w')
count = 0
average = [0.]*3
max = [-1.e20]*3
min = [1.e20]*3

while time <= end_time:
  for i in range(3):
    well_plane.setZ(i+1,well[i].getValueAtTime(time))
  well_plane.computePlaneCoefficients()
  datum.setZ(well_plane.computeZ(datum.getX(),datum.getY()))            
  datum_file_out.write("%.7e" % (float(count)*3600.)) # convert to seconds
  datum_file_out.write(" %.7e %.7e %.7e\n" % (datum.getX(), \
                                   datum.getY(),datum.getZ()))
  gradz_x = well_plane.computeGradientZ_X(datum.getX(),datum.getY())
  gradz_y = well_plane.computeGradientZ_Y(datum.getX(),datum.getY())
  gradient_file_out.write("%.7e" % (float(count)*3600.))
  gradient_file_out.write(" %.7e %.7e %.7e\n" % (gradz_x,gradz_y,0.))
  average[0] += datum.getZ()
  average[1] += gradz_x
  average[2] += gradz_y
  if datum.getZ() < min[0]:
    min[0] = datum.getZ()
  if datum.getZ() > max[0]:
    max[0] = datum.getZ()
  if gradz_x < min[1]:
    min[1] = gradz_x
  if gradz_x > max[1]:
    max[1] = gradz_x
  if gradz_y < min[2]:
    min[2] = gradz_y
  if gradz_y > max[2]:
    max[2] = gradz_y
  print count, time.ctime() 
  time += time_increment
  count += 1

#  print count
for i in range(3):
  average[i] /= float(count)
print 'average:', average
print 'max:', max
print 'min:', min
for i in range(3):
  well[i].finalize()
f.close()  
datum_file_out.close()
gradient_file_out.close()
print 'done'  