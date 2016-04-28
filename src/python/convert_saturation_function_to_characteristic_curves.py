# Saturation function to Characteristic Curve Converter for PLOTRAN ----------!
# This short script converts an old saturation function input files to
# the new characteristic curves file (ONLY in GENERAL MODE).
# Go to the bottom of this file for main script input.
#
# - Heeho Park 10/14/2014

import os
import sys

#----------- functions (main script input is located at the bottom) ----------!

def scanfolder(parent_dir,extension):
    file_list = []
    # find all files with an extension under parent dir and its sub dirs
    for path, subdirs, files in os.walk(parent_dir):
        for target_file in files:
            # probably not the best way to find file extensions but it works
            # to some extent
            if target_file.endswith(extension):
                file_list.extend([os.path.join(path,target_file)])
    #return file_path strings in lists
    return file_list

def convert_sf2cc(filename,delete_files,old_extension):
	"""converts saturation function card into chracteristic curve card
	filename = filename including the path and extension
	delete_files = True or False"""

	try:
		old = open(filename,'r')
	except:
		print 'File: ' + filename + ' does not exist. skip.\n'
		return

#	#if the file is not in general mode, it skips the file
#	for line in old:
#		words = line.split()
#		if words == []:
#			continue
#		elif words[0].upper() == 'MODE':
#			if words[1].upper() != 'GENERAL':
#				print 'File below is not in GENERAL MODE thus skips.'
#				print filename+'\n'
#				return
	old.close()

	#renames the old file and keep the new files as original filename
	#but does not rename something that has been renamed.
	#old files stays old; new file stays new
	if filename[-len(old_extension):] == old_extension:
		old_filename = filename
		filename = filename[:-len(old_extension)]+'.in'
		new = open(filename,'w+')
	else:
		old_filename = filename[:-3]+old_extension
		try:
			os.rename(filename,old_filename)
		except:
			pass
		new = open(filename,'w+')

	try:
		old = open(old_filename,'r')
	except:
		print 'File: ' + old_filename + ' does not exist. skip.'
		return

	# 	defaults
	flag = ''
	ccname = ''
	rgas = '-999.0'
	rliq = '-999.0'
	lamb = '-999.0'
	alph = '-999.0'
	mcap = '1.d+08'
	satf = 'BROOKS_COREY'
	perf = 'MUALEM'
	for line in old:
		words = line.split()
		
		if words == []:
			new.write(line)
			continue

		for i in range(len(words)):
			# convert to uppercase
			words[i] = words[i].upper()
		
		if flag == 'mp':
			#replace keywords under material properties
			if words[0] == 'SATURATION_FUNCTION':
				new.write('  CHARACTERISTIC_CURVES '+words[1]+'\n')
				ccname = words[1]
				flag = ''
			elif words[0] == 'END':
				new.write(line)
				flag = ''
			else:
				new.write(line)
		elif flag == 'sf':
			if words[0] == 'SATURATION_FUNCTION_TYPE':
				new.write('  SATURATION_FUNCTION '+words[1]+'\n')
				satf = words[1]
			elif words[0] == 'RESIDUAL_SATURATION':
				if words[1] == 'LIQUID' or words[1] == 'WATER' or \
				words[1] == 'LIQUID_PHASE' or words[1] == 'WATER_PHASE':
					rliq = words[2]
				elif words[1] == 'GAS' or words[1] == 'CO2' or \
				words[1] == 'GAS_PHASE' or words[1] == 'CO2_PHASE':
					rgas = words[2]
				else:
					rliq = words[1]
			elif words[0] == 'LAMBDA':
				lamb = words[1]
			elif words[0] == 'ALPHA':
				alph = words[1]
			elif words[0] == 'MAX_CAPILLARY_PRESSURE':
				mcap = words[1]
			elif words[0] == 'PERMEABILITY_FUNCTION_TYPE':
				perf = words[1]
			elif words[0] == 'END' or words[0] == '/':
				flag = ''
				if rgas == '-999.0' or rliq == '-999.0' or \
				lamb == '-999.0' or alph == '-999.0':
				    print 'SATURATION FUNCTION' + ccname + \
				    'is missing some parameters.'
				    print 'Please check the converted file:' + filename
				if satf == 'BROOKS_COREY':
					new.write('    LAMBDA '+lamb+'\n')
				elif satf == 'VAN_GENUCHTEN':
					new.write('    M '+lamb+'\n')
				new.write('    ALPHA  '+alph+'\n')
				new.write('    LIQUID_RESIDUAL_SATURATION '+rliq+'\n')
				new.write('    MAX_CAPILLARY_PRESSURE '+mcap+'\n')
				new.write('  /\n')
				if perf == 'MUALEM' and satf == 'BROOKS_COREY':
					new.write('  PERMEABILITY_FUNCTION MUALEM_BC_LIQ\n')
				elif perf == 'BURDINE' and satf == 'BROOKS_COREY':
					new.write('  PERMEABILITY_FUNCTION BURDINE_BC_LIQ\n')
				elif perf == 'MUALEM' and satf == 'VAN_GENUCHTEN':
					new.write('  PERMEABILITY_FUNCTION MUALEM_VG_LIQ\n')
				elif perf == 'BURDINE' and satf == 'VAN_GENUCHTEN':
					new.write('  PERMEABILITY_FUNCTION BURDINE_VG_LIQ\n')
				if satf == 'BROOKS_COREY':
					new.write('    LAMBDA '+lamb+'\n')
				elif satf == 'VAN_GENUCHTEN':
					new.write('    M '+lamb+'\n')
				new.write('    LIQUID_RESIDUAL_SATURATION '+rliq+'\n')
				new.write('  /\n')
				if perf == 'MUALEM' and satf == 'BROOKS_COREY':
					new.write('  PERMEABILITY_FUNCTION MUALEM_BC_GAS\n')
				elif perf == 'BURDINE' and satf == 'BROOKS_COREY':
					new.write('  PERMEABILITY_FUNCTION BURDINE_BC_GAS\n')
				elif perf == 'MUALEM' and satf == 'VAN_GENUCHTEN':
					new.write('  PERMEABILITY_FUNCTION MUALEM_VG_GAS\n')
				elif perf == 'BURDINE' and satf == 'VAN_GENUCHTEN':
					new.write('  PERMEABILITY_FUNCTION BURDINE_VG_GAS\n')
				if satf == 'BROOKS_COREY':
					new.write('    LAMBDA '+lamb+'\n')
				elif satf == 'VAN_GENUCHTEN':
					new.write('    M '+lamb+'\n')
				new.write('    LIQUID_RESIDUAL_SATURATION '+rliq+'\n')
				new.write('    GAS_RESIDUAL_SATURATION    '+rgas+'\n')
				new.write('  /\n')
				new.write('END\n\n')
				ccname = ''
				rgas = '-999.0'
				rliq = '-999.0'
				lamb = '-999.0'
				alph = '-999.0'
				mcap = '-999.0'
				satf = 'BROOKS_COREY'
				perf = 'BURDINE'
		elif flag == 'cf':
			if words[0] == 'END':
				flag = ''
				new.write(line)
			elif words[0] == 'CHARACTERISTIC_CURVES':
				new.write('SATURATION_FUNCTION '+words[1]+'\n')
			else:
				new.write(line)
		elif words[0] == 'MATERIAL_PROPERTY':
			flag = 'mp'
			new.write(line)
		elif words[0] == 'CHARACTERISTIC_CURVES':
			flag = 'cf'
			new.write(line)
		elif words[0] == 'SATURATION_FUNCTION':
			flag = 'sf'
			new.write('CHARACTERISTIC_CURVES '+words[1]+'\n')
		else:
			new.write(line)
	old.close()
	new.close()
	print filename+' Converted'
	if delete_files:
		os.remove(old_filename)
		print old_filename + ' Deleted\n'


############################## MAIN PROGRAM ###################################
# This code will change parent directory files and all subdirectory files.
# set one_file to '' if you want to have all files under a directory 
# to be converted
one_file = ''
parent_dir = 'C:/Sandia/98806/PFLOTRAN_test_cases/case1/pflotran'
extension = '.in'
# False if you want to keep old files which will have old_file_extension 
delete_old_files = False
# extension you'd like to add to your old file
old_file_extension = '_orig.in'

if one_file != '':
	fn = one_file + extension
	convert_sf2cc(fn,delete_old_files,old_file_extension)
else:
	all_file_path = scanfolder(parent_dir,extension)
	for i in xrange(len(all_file_path)):
		filename = all_file_path[i]
		convert_sf2cc(filename,delete_old_files,old_file_extension)
print 'Your original input Permeability Function is used for GAS/LIQ phase. Change as needed.'
