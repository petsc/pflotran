import sys
import os
import math
from WellData import *
from Plane import *
from datetime import *

new_origin = Point3D(0.,0.,0.)
rotation_angle = 0. # positive counter clock-wise

lines = []

print 'Enter input file name:'
filename = sys.stdin.readline().strip()
source_filename = filename
source_path = os.getcwd()

f = open(filename,'r')
s = f.readline()
lines.append(s)
w = s.split(':')
w = w[1].split()
origin = Point3D(float(w[0]),float(w[1]),0.)
delta = Point3D(new_origin.x-origin.x,new_origin.y-origin.y,0.)

s = f.readline()
lines.append(s)
w = s.split(':')
rotation_angle = float(w[1])
print 'Rotation angle: ', rotation_angle

origin.translate_and_rotate(delta,rotation_angle)
origin.setName('Origin')
origin.printInfo()

s = f.readline()
lines.append(s)
w = s.split(':')
w = w[1].split()
datum = Point3D(float(w[0]),float(w[1]),0.)
datum.translate_and_rotate(delta,rotation_angle)
datum.setName('Datum')
datum.printInfo()

well = []
for i in range(3):
  s = f.readline()
  lines.append(s)
  w = s.split(":")
  filename = w[1].strip()
  well.append(WellData(filename))
  well[i].location.translate_and_rotate(delta,rotation_angle)
  well[i].location.setName(w[0].strip())
  well[i].printInfo()

well_plane = Plane(well[0].location,well[1].location,well[2].location)

s = f.readline()
lines.append(s)
w = s.split(": ")
start_date = readDateFromString(w[1])

s = f.readline()
lines.append(s)
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
lines.append(s)
w = s.split(":")
filename_base = w[1].strip()
datum_file_out = open(filename_base+".datum",'w')
gradient_file_out = open(filename_base+".gradient",'w')
count = 0

sum_datum = [0.]*3
sum_gradient = [0.]*3

max_datum = [-1.e20]*3
min_datum = [1.e20]*3

max_gradient = [-1.e20]*3
min_gradient = [1.e20]*3

datum_file_out.write(": Source File - %s/%s\n:\n" % (source_path,source_filename))
gradient_file_out.write(": Source File - %s/%s\n:\n" % (source_path,source_filename))

datum_file_out.write(": Datafile generation instructions\n")
gradient_file_out.write(": Datafile generation instructions\n")

for i in range(len(lines)):
  datum_file_out.write(": %s" % lines[i])
  gradient_file_out.write(": %s" % lines[i])

datum_file_out.write(":\n: Data - Time(seconds) X(meters) Y(meters) Z(meters)\n")
gradient_file_out.write(":\n: Data - Time(seconds) GradX(-) GradY(-) GradZ(-)\n")

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
  gradz_z = 0.
  gradient_file_out.write("%.7e" % (float(count)*3600.))
  gradient_file_out.write(" %.7e %.7e %.7e\n" % (gradz_x,gradz_y,0.))
  
  
  sum_datum[0] += datum.getX()
  sum_datum[1] += datum.getY()
  sum_datum[2] += datum.getZ()
  
  sum_gradient[0] += gradz_x
  sum_gradient[1] += gradz_y
  sum_gradient[2] += gradz_z
  
  
  if datum.getX() < min_datum[0]:
    min_datum[0] = datum.getX()
  if datum.getY() < min_datum[1]:
    min_datum[1] = datum.getY()
  if datum.getZ() < min_datum[2]:
    min_datum[2] = datum.getZ()
    
  if datum.getX() > max_datum[0]:
    max_datum[0] = datum.getX()
  if datum.getY() > max_datum[1]:
    max_datum[1] = datum.getY()
  if datum.getZ() > max_datum[2]:
    max_datum[2] = datum.getZ()
    
    
  if gradz_x < min_gradient[0]:
    min_gradient[0] = gradz_x
  if gradz_y < min_gradient[1]:
    min_gradient[1] = gradz_y
  if gradz_z < min_gradient[2]:
    min_gradient[2] = gradz_z
    
  if gradz_x > max_gradient[0]:
    max_gradient[0] = gradz_x
  if gradz_y > max_gradient[1]:
    max_gradient[1] = gradz_y
  if gradz_z > max_gradient[2]:
    max_gradient[2] = gradz_z
    
  print count, time.ctime() 
  time += time_increment
  count += 1

#  print count
for i in range(3):
  sum_gradient[i] /= float(count)
  sum_datum[i] /= float(count)
print 'datum average:', sum_datum[0], sum_datum[1], sum_datum[2]
print 'datum max:', max_datum[0], max_datum[1], max_datum[2]
print 'datum min:', min_datum[0], min_datum[1], min_datum[2]
print 'gradient average:', sum_gradient[0], sum_gradient[1], sum_gradient[2]
print 'gradient max:', max_gradient[0], max_gradient[1], max_gradient[2]
print 'gradient min:', min_gradient[0], min_gradient[1], min_gradient[2]
for i in range(3):
  well[i].finalize()
f.close()  
datum_file_out.close()
gradient_file_out.close()
print 'done'  