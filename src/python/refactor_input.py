# swap_comment_character.py
import sys
import shutil
import os
import fnmatch

simulation = '''\nSIMULATION
  SIMULATION_TYPE %s
  PROCESS_MODELS\n'''

flow_string = '''    SUBSURFACE_FLOW flow
      MODE %s
    /\n'''

transport_string = '''    SUBSURFACE_TRANSPORT transport
      GLOBAL_IMPLICIT
    /\n'''

simulation_end = '''  /
END\n\n'''

subsurface_string = 'SUBSURFACE\n\n'
end_subsurface_string = 'END_SUBSURFACE\n'

surface_string = '\nSURFACE\n\n'
end_surface_string = '\nEND_SURFACE\n'

def refactor_file(filename):
  # search through file for keywords
  flow = False
  transport = False
  for line in open(filename):
    if line.strip().startswith('MODE'):
      w = line.strip().split()
      flow_mode = w[1]
      flow = True
    elif line.strip().startswith('CHEMISTRY'):
      transport = True
    if flow and transport:
      break
  f = open(filename,'r')
  f2 = open(filename+'.tmp','w')
  # copy comment lines
  for line in f:
    if not line.strip().startswith('#'):
      break
    f2.write(line)
  # add new simulation cards
  if flow and transport:
    f2.write(simulation % 'SUBSURFACE_FLOW_AND_TRAN')
  elif flow:
    f2.write(simulation % 'SUBSURFACE_FLOW')
  elif transport:
    f2.write(simulation % 'SUBSURFACE_TRANSPORT')
  if flow:
    f2.write(flow_string % flow_mode)
  if transport:
    f2.write(transport_string)
  f2.write(simulation_end)
  f2.write(subsurface_string)
  # write remainder of file
  surface_found = False
  for line in f:
    if line.strip().startswith('SURFACE_FLOW') and not surface_found:
      surface_found = True
      f2.write(end_subsurface_string)
      f2.write(surface_string)
    f2.write(line)
  # append end string
  if surface_found:
    f2.write(end_surface_string)
  else:
    f2.write(end_subsurface_string)
  f.close()
  f2.close()
  # using shutil.move adds ^M to end of lines.
  os.remove(filename)
  shutil.copy(filename+'.tmp',filename)
  os.remove(filename+'.tmp')

suffix = '*.in'
for root, dirnames, filenames in os.walk('.'):
  for filename in fnmatch.filter(filenames,suffix):
    filename = os.path.join(root,filename)
    print(filename)
    refactor_file(filename)

print('done')
