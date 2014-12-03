# This python script will generate the <Files></Files> portion of the *.vfproj 
# file for PFLOTRAN. It will read the MakeFile to create the <Files> portion
# here is the instruction
# 1. define source_path
# 2. define printout
# 3. define makefile
# 4. open the textfile written by printout
# 5. ctrl+a, ctrl+c: copy the whole text
# 6. go to your visual studio path 
#    e.g. C:\Users\heepark\Documents\Visual Studio 2013\Projects\PROJECT_NAME\PROJECT_NAME
# 7. open PROJECT_NAME.vfproj
# 8. paste it over <Files> ... </Files>
# 9. clean build. 

# Define your *.vfproj relative path to your source code 'src\pflotran'
# in Windows backslash format. no backslash at the end.
source_path = '..\..\..\..\..\..\..\software\pflotran-dev\src\pflotran'

# where to save the printout with the file name
printout = 'C:\Users\heepark\Desktop\copypaste.txt'

# where the Make file is relative to this file location
makefile = '..\pflotran\makefile'

mf = open(makefile,'r')
po = open(printout,'w+')

po.write('  <Files>\n'
'    <Filter Name="Header Files" Filter="fi;fd"/>\n'
'    <Filter Name="Resource Files" Filter="rc;ico;cur;bmp;'
'dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"/>\n'
'    <Filter Name="Source Files" Filter="f90;for;f;fpp;ftn;def;odl;idl">\n')
for line in mf:
	if line.strip().startswith('# Begin Source Block'):
		for line in mf:
			if line.strip().startswith('${common_src}'):
				filename = line.split('}')[1]
				filename = filename.split('.')[0]
				filename = '\\' + filename + '.F90'
				po.write('    <File RelativePath="')
				po.write(source_path)
				po.write(filename)
				po.write('"/>\n')
			elif line.strip().startswith('# End Source Block'):
				break

po.write('  </Filter></Files> ')

mf.close()
po.close()




