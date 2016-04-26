import os

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
    
def replace_oldchar_with_newchar(filename,oldchar,newchar):
    #may take your file size of RAM but works simply.
    f = open(filename,'r')
    oldfile = f.read()
    f.close()

    newfile = oldfile.replace(oldchar,newchar)

    f = open(filename,'w+')
    f.write(newfile)
    f.close()
    print filename + ' done!'
    

# This code will change parent directory files and all subdirectory files.
parent_dir = 'C:/sandia/98806/python'
extension = '.in'
all_file_path = scanfolder(parent_dir,extension)
# appends all the files names that you have changed
output = open('list_of_changed_files','a+')

for i in xrange(len(all_file_path)):
    filename = all_file_path[i]
    replace_oldchar_with_newchar(filename, ':', '!')
    output.write(filename+'\n')

output.close()