# swap_comment_character.py
import sys
import shutil
import os
import fnmatch

simulation = '''\nSIMULATION
  SIMULATION_TYPE %s
  PROCESS_MODELS\n'''

flow_string = '''    SUBSURFACE_FLOW flow
      MODE %s\n'''

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
  surface_flow = False
  flow_option_strings = []
  skip_count = 0
  f = open(filename,'r')
  while 1:
    line = f.readline()
    if not line:
      break
    string = line.strip().upper()
    if string.startswith('SKIP'):
      skip_count += 1
    elif string.startswith('NOSKIP'):
      skip_count -= 1
    elif skip_count == 0:
      if string.startswith('MODE'):
        w = string.strip().split()
        flow_mode = w[1]
        flow = True
        if len(w) > 2 and w[2].startswith('MORE'):
          while 1:
            string = f.readline()
            if string.strip().startswith('/') or \
               string.strip().startswith('END'):
              break
            flow_option_strings.append(string)
        elif string.strip().endswith(' FREEZING'):
          flow_option_strings.append('FREEZING')
      elif string.startswith('ICE_MODEL'):
        flow_option_strings.append(string)
      elif string.startswith('CHEMISTRY'):
        transport = True
      elif string.startswith('SURFACE_FLOW'):
        surface_flow = True
    if flow and transport and surface_flow:
      break
  f.close()
  if surface_flow:
    # skip surface flow files.
    print('\n%s skipped due to surface flow\n' % filename)
    return

#  print('Options:')
#  for i in range(len(flow_option_strings)):
#    print(flow_option_strings[i])
#  print('after')

  f = open(filename,'r')
  f2 = open(filename+'.tmp','w')
  # copy comment lines
  for line in f:
    if not (line.strip().startswith('#') or line.strip().startswith('!')):
      # if we hit mode on the first non-comment line, we have to remove
      # the entire block - argh
      if line.strip().startswith('MODE'):
        if line.strip().endswith('MORE'):
          for line in f:
            if line.strip().startswith('/') or line.strip().startswith('END'):
              break
      break
    f2.write('%s\n' % line.rstrip())
  # add new simulation cards
#  if flow and transport:
#    f2.write(simulation % 'SUBSURFACE_FLOW_AND_TRAN')
#  elif flow:
#    f2.write(simulation % 'SUBSURFACE_FLOW')
#  elif transport:
#    f2.write(simulation % 'SUBSURFACE_TRANSPORT')
  f2.write(simulation % 'SUBSURFACE')
  if flow:
    f2.write(flow_string % flow_mode)
    if len(flow_option_strings) > 0:
      f2.write('      OPTIONS\n')
      for i in range(len(flow_option_strings)):
        f2.write('        %s\n' % flow_option_strings[i].strip())
      f2.write('      /\n')
    f2.write('    /\n')
  if transport:
    f2.write(transport_string)
  f2.write(simulation_end)
  f2.write(subsurface_string)
  # write remainder of file
  surface_found = False
  for line in f:
    line = line.rstrip()
    if line.strip().startswith('MODE'):
      if line.strip().endswith('MORE'):
        for line in f:
          if line.strip().startswith('/') or line.strip().startswith('END'):
            break
      continue
    elif line.strip().startswith('ICE_MODEL'):
      continue
    if line.strip().startswith('SURFACE_FLOW') and not surface_found:
      surface_found = True
      f2.write(end_subsurface_string)
      f2.write(surface_string)
    f2.write('%s\n' % line)
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
