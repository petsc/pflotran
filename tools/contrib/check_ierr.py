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

        # strip off possible existing error checks
        call = call.replace("CHKERRQ(ierr)","").replace(";","")

        # Now search for important information in that string
        single_call = re.search(r"(\s*)call\s(\w*)\((.*)\).*",call)
        if single_call:
            self.call_tab  = single_call.group(1)
            self.call_name = single_call.group(2)
            self.call_args = single_call.group(3)
            self.call_args = self.call_args.split(",")
            for i,arg in enumerate(self.call_args):
                self.call_args[i] = arg.strip()

        self.call = self.call_name + "(" + ",".join(self.call_args) + ")"
        return
    
    def __str__(self):
        return self.call

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
        if (lines[self.row+self.shift  ].count("CHKERRQ") == 1): return True        
        if (lines[self.row+self.shift+1].count("CHKERRQ") == 1 and
            lines[self.row+self.shift+1].count(" call ") == 0): return True        
        return False

    def AddErrorCheck(self,lines,column_limit=80):
        if self.IsErrorChecked(lines): return
        row    = self.row
        shift  = self.shift
        lenchk = len(";CHKERRQ(ierr)")

        # If the call is part of a single-line 'if' statement, break it up
        single_if = re.search("(\s*)if.*(\s*)call\s.*",lines[row])
        if single_if:
            lines[row]       = lines[row].replace("call","then\n%s  call" % single_if.group(1))
            lines[row+shift] = lines[row+shift].replace("\n","\n%sendif\n" % self.call_tab[:-2])

        # If the line above the call is a single-line 'if' with continuation, break it up
        contin_if = re.search("if.*&",lines[row-1])
        if contin_if:
            lines[row-1]     = lines[row-1].replace("&","then")
            lines[row+shift] = lines[row+shift].replace("\n","\n%sendif\n" % self.call_tab[:-2])

        # stand-alone call replacement if the line isn't too long
        if len(lines[row+shift]) < column_limit-lenchk:
            lines[row+shift] = lines[row+shift].replace("ierr)","ierr);CHKERRQ(ierr)")
            return False
        
        # wrap the last argument and add check
        lines[row+shift] = lines[row+shift].replace("ierr)"," &\n%s%sierr);CHKERRQ(ierr)" % (self.call_tab," "*(len(self.call_name)+6)))
        return True

added  = []
nwraps = 0
for (dirpath, dirnames, filenames) in os.walk("./"):
    for filename in filenames:
        if filename[-4:] != ".F90": continue
        lines = file(filename).readlines()
        for row,line in enumerate(lines):
            if IsCommentedLine(line): continue
            if line.count(" call ") > 0:
                fc = FortranCall(row,lines)
                if fc.IsPETScCall():
                    if fc.IsErrorChecked(lines) == False:
                        added.append("%s:%d %s" % (filename,row+1,fc.call_name))
                        nwraps += fc.AddErrorCheck(lines)
        file(filename,"w").writelines(lines)

if len(added) == 0:
    print "All PETSc calls have ierr checked"
else:
    print "Added CHKERRQ(ierr) to %d PETSc call(s)" % (len(added))
    print "Had to multiline %d call(s) due to a 80 character column limit violation" % (nwraps)
