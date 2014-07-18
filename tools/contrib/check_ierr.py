import os,re

def IsCommentedLine(line):
    try:
        if (line.strip())[0] == "!": return True
    except:
        pass
    return False

class FortranCall:
    def __init__(self,row,lines):

        # Initialize
        self.row       = row
        self.call_tab  = ""
        self.call_name = ""
        self.call_args = []

        # Take possible multiline calls and create a single string with the call in it
        call   = lines[row]
        popen  = lines[row].count("(")
        pclose = lines[row].count(")")
        self.shift  = 0
        while popen > pclose:
            self.shift  += 1
            if IsCommentedLine(lines[row+self.shift]): continue
            call   += lines[row+self.shift]
            popen  += lines[row+self.shift].count("(")
            pclose += lines[row+self.shift].count(")")
        call = call.replace("\n","").replace("&","")

        # Now search for important information in that string
        single_call = re.search(r"(\s*)call\s(\w*)\((.*)\).*",call)
        if single_call:
            self.call_tab  = single_call.group(1)
            self.call_name = single_call.group(2)
            self.call_args = single_call.group(3)
            self.call_args = self.call_args.split(",")
            for i,arg in enumerate(self.call_args):
                self.call_args[i] = arg.strip()
        return
    
    def IsPETScCall(self):
        sub = self.call_name.lower()
        prefix = False
        if (sub[:5] == "petsc" or sub[:3] == "mat"  or sub[:3] == "vec" or
            sub[:3] == "ksp"   or sub[:4] == "snes" or sub[:2] == "ts"  or
            sub[:2] == "ao"    or sub[:2] == "dm"   or sub[:2] == "pc"  or
            sub[:2] == "is"): prefix = True
        if prefix and (self.call_args[-1]).strip() == "ierr": return True
        return False

    def IsErrorChecked(self,lines):
        try:
            if (lines[self.row+self.shift+1].count("CHKERRQ(ierr)")>0): return True
        except:
            pass
        return False

    def AddErrorCheck(self,lines):
        if self.IsErrorChecked(lines): return
        row   = self.row
        shift = self.shift

        # is the call part of a single-line if statement?
        single_if = re.search("(\s*)if.*(\s*)call\s.*",lines[row])
        if single_if:
            lines[row] = lines[row].replace("call","then\n%s  call" % single_if.group(1))
            lines[row+shift] = lines[row+shift] + "%s  CHKERRQ(ierr)\n%sendif\n" % (single_if.group(1),single_if.group(1))
            return

        # the line above the call could be a single-line if with continuation
        contin_if = re.search("if.*&",lines[row-1])
        if contin_if:
            lines[row-1]     = lines[row-1].replace("&","then")
            lines[row+shift] = lines[row+shift].replace("\n","\n%sendif\n" % self.call_tab[:-2])

        # standalone call replacement
        lines[row+shift] = lines[row+shift].replace("ierr)","ierr)\n%sCHKERRQ(ierr)" % self.call_tab)

added = []
for (dirpath, dirnames, filenames) in os.walk("./"):
    for filename in filenames:
        if filename[-4:] != ".F90": continue
        lines = file(filename).readlines()
        for row,line in enumerate(lines):
            if IsCommentedLine(line): continue
            if line.count("call") > 0:
                fc = FortranCall(row,lines)
                if fc.IsPETScCall():
                    if fc.IsErrorChecked(lines) == False:
                        added.append("%s:%d %s" % (filename,row+1,fc.call_name))
                        fc.AddErrorCheck(lines)
        file(filename,"w").writelines(lines)

for add in added:
    print add
if len(added) == 0:
    print "All PETSc calls have ierr checked"
else:
    print "Added CHKERRQ(ierr) to the above %d PETSc call(s)" % (len(added))
