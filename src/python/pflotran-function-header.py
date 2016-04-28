#!/usr/bin/env python

"""PFLOTRAN function header refactor

* move function header information just below the function definition
  line.

* cleanup header comments, keep all of the old "author / creation
  date" info but remove doxygen style comments, etc

* Insert a new separator between functions.

Author: Ben Andre <bandre@lbl.gov>

"""

from __future__ import print_function

from copy import deepcopy
import fileinput
import os
import re
import sys
import traceback



def process_file(filename):
    print("Processing '{0}' ... ".format(filename), end='')
    
    with open(filename, 'r') as orig_file:
        src = orig_file.readlines()

    new_src = deepcopy(src)

    functions = process_functions(src)
    for f in functions:
        if f.has_key('comment'):
            cleanup_comments(f)
        update_function_header(f, new_src)


    with open(filename, 'w') as src:
        for l in new_src:
            src.write(l)

    print("done")

def process_functions(src):
    """loop through the file backwards and record the start of each
     function.  assume that any line following the begining of a
     function that is not a comment or empty is the end of the
     function header.

    """

    func_sub_re = re.compile("^[\s\w]*(function|subroutine)[\s]+([\w]+)[\s]*\(",
                             re.IGNORECASE)
    # catch the end of commented out functions
    end_func_sub_re = re.compile("[\s!]*end[\s]*(function|subroutine)", re.IGNORECASE)

    contains_re = re.compile("^[\s]*contains[\s]*$")

    functions = []
    start_func = None
    end_func = None
    for num, line in reversed(list(enumerate(src))):
        #print("{0} : {1}".format(num, line), end='')

# Move to bottom of loop - geh
#        if contains_re.match(line):
#            # stop if we find a contains statement for a module
#            break

        if start_func > 0:
            temp = line.strip()
            found_end_func = end_func_sub_re.match(temp)
            if len(temp) > 0 and (temp[0] != '!' or found_end_func):
                # NOTE: num + 1 because we want to move back down a
                # line to the last comment/blank
                end_func = num + 1
                functions[-1]['end'] = end_func
                comment = src[end_func:start_func+1]
                functions[-1]['comment'] = comment
                start_func = None
                end_func = None

        match = func_sub_re.search(line)
        if match:
            # NOTE: num - 1 because we want to move up a line before
            # the start of the function!
            start_func = num - 1
            if start_func < 0:
                # special case of only functions (no module, no leading comments)
                start_func = 0
            functions.append({'type': match.group(1), 'name': match.group(2), 'start': start_func})

            #print("file : {0} --> found {1}".format(filename, match.group(0)))
            #print(match.groups())
 
        # placed at bottom so that subroutine/function immediately below 
        # contains statement is updated too. - geh
        if contains_re.match(line):
            # stop if we find a contains statement for a module
            break

    # for files with just free functions instead of modules
    if start_func != None and end_func == None:
        functions[-1]['end'] = start_func

    return functions


def cleanup_comments(func):
    """go through the header comments we identified for this function and
     remove every thing we don't want. Reformat what is left and save
     the final comment.

    """
    separator_1_re = re.compile("^![\s]*[*=cC\-@]+[\s]*!*$") # '! ******* !' style comments
    separator_2_re = re.compile("^!$") # empty comments '!'
    separator_3_re = re.compile("^![\s]*[*\-]{6,}[*a-zA-Z0-9 ]+[*\-]+[\s!]*$") # '! *** words **' or '!***1***2' style comments
    author_re = re.compile("^![\s]*author[\s]*:[\s]*(.+)", re.IGNORECASE)
    written_by_re = re.compile("^![\s]*written by[\s:]*(.+)", re.IGNORECASE)
    dox_author_re = re.compile("^![\s]+@author[\s]*", re.IGNORECASE)
    date_re = re.compile("^![\s]*date[\s]*:[\s]*(.+)", re.IGNORECASE)

    author_str = "Author: {0}"
    date_str = "Date: {0}"
    comment_str = "  ! {0}\n"

    # NOTE: since we traversed the file backwords when we created the
    # function list, the list starts at the bottom of the file. This
    # is what we want because we can modify lines based on indices
    # without invalidating the indices for other functions!
    author_next = False
    new_comment = []
    author = None
    date = None
    for c in func['comment']:
        c = c.strip()
        # convert various comment styles to the same format
        c = c.replace("!!", '!')
        c = c.replace("!>", '!')
        c = c.replace("!c", '!')
        c = c.replace("!/*", '!')
        c = c.replace("!*/", '!')
        c = c.replace("!*", '!')

        if author_next is True:
            author_next = False
            c = c.replace('!', '')
            c = c.strip()
            author = author_str.format(c)
        elif len(c) == 0:
            #print("--> skipping : '{0}'".format(c))
            pass
        elif separator_2_re.search(c):
            #print("--> skipping : '{0}'".format(c))
            pass
        elif separator_1_re.search(c):
            #print("--> skipping : '{0}'".format(c))
            pass
        elif separator_3_re.search(c):
            #print("--> skipping : '{0}'".format(c))
            pass
        elif author_re.search(c):
             #print("--> matched author_re : '{0}'".format(c))
            match = author_re.search(c)
            c = match.group(1).strip()
            author = author_str.format(c)
        elif written_by_re.search(c):
            #print("--> matched written_by_re : '{0}'".format(c))
            match = written_by_re.search(c)
            c = match.group(1).strip()
            author = author_str.format(c)
        elif dox_author_re.search(c):
            #print("--> skipping : '{0}'".format(c))
            author_next = True
        elif date_re.search(c):
            #print("--> {0} mathched date_re : '{1}'".format(func['name'], c))
            match = date_re.search(c)
            c = match.group(1).strip()
            date = date_str.format(c)
        else:
            cmt = c.find('!')
            if cmt > 0:
                c = c[cmt:]
            new_comment.append(c)

    for i, c in enumerate(new_comment):
        # match the function name followed by some text
        func_name_re = "^![\s]*{0}[\s]*:(.*)".format(func['name'])
        match = re.search(func_name_re, c, re.IGNORECASE)
        if match:
            c = match.group(1).strip()
        else:
            #print("no match for {0} : '{1}'".format(func_name_re, c))
            pass
        # match any other comment text
        match = re.search("^![\s]*([\S]+.*)", c) # remove whitspace between ! and text
        if match:
            c = match.group(1).strip()
        new_comment[i] = c

    # remove empty comments
    new_comment = [x for x in new_comment if x != '']
    # reinsert the comment characters
    for i, c in enumerate(new_comment):
        if len(c) > 0:
            new_comment[i] = comment_str.format(c)
    # insert author and date if found
    if (len(new_comment) > 0 and 
        (author is not None or date is not None)):
        # inserate a blank line separator
        new_comment.append(comment_str.format(''))
    if author is not None:
        new_comment.append( comment_str.format(author))
    if date is not None:
        new_comment.append( comment_str.format(date))
    # pad the new comment with blank lines
    if len(new_comment) > 0:
        new_comment.insert(0, comment_str.format(''))
        new_comment.insert(len(new_comment), comment_str.format(''))

    func['comment'] = new_comment


def update_function_header(func, new_src):

    function_separator = """
! ************************************************************************** !

"""

    #print_func_dict(func)

    # new we insert the new comment in the correct location. we do
    # this first because it occurs after the text we are going to
    # remove and this will not invalidate our indicies!
    if func.has_key('comment'):
        start = func['start'] + 1
        while new_src[start].find('&') > 0:
            start += 1
        for i, c in enumerate(func['comment']):
            new_src.insert(start + 1 + i, c)

    if func['end'] != 0:
        # now we remove the old comment. NOTE: This will invalidate
        # the indices for other functions that occur below this one in
        # the file!
        for i in reversed(range(func['end'], func['start']+1)):
            #print("{0} : {1}".format(i, new_src[i]), end='')
            new_src.pop(i)

        # now we insert the function separator
        new_src.insert(func['end'], function_separator)


def print_func_dict(func):
    print("{0} :".format(func['name']))
    print("  type = {0}".format(func['type']))
    print("  range = ({0}, {1})".format(func['start'], func['end']))
    print("  comment :")
    if func.has_key('comment'):
        for c in func['comment']:
            print("    {0}".format(func['comment']))


def main():
    if len(sys.argv) > 1:
        file_list = sys.argv[1:]
    else:
        raise RuntimeError("must provide a list of files to process on the command line.")

    for f in file_list:
        if not os.path.isfile(f):
            raise RuntimeError("{0} is not a valid file!")
        process_file(f)

    return 0

if __name__ == "__main__":
    try:
        status = main()
        sys.exit(status)
    except Exception as error:
        print(str(error))
        traceback.print_exc()
        sys.exit(1)
