
# ========================================================== #
#                   Python utility module                    #
#                                                            #
#                       version 1.0                          #
#                                                            #
#                   authors: Ing. Matteo Pini                #
# ========================================================== #

import numpy as np

# ================================ #
# Types conversion
# ================================ #

# convert an integer into a boolean
def int_to_bool (inval):
    if   inval == 0:
        outval = False
    elif inval == 1:
        outval = True
    else:
        print('Error: input value has to be 0 or 1')
        quit()
    return outval   
    
# convert a string into a boolean
def str_to_bool (instr):
    if   instr == 'False' or instr == 'false':
        outval = False
    elif instr == 'True'  or instr == 'true':
        outval = True
    else:
        print('Error: only (T)(t)rue/(F)(f)alse input strings are admitted')
        quit()
    return outval   
    
# ================================ #
# Files management
# ================================ #

# create file
def create_file(filename):
    file_name = open(filename,'w')
    file_name.close()
    
# copy file into a given folder
def copy_file_to_folder(filename,folder):
    import shutil as sh
    sh.copy(filename,folder)
    
# copy file from a folder into the current directory
def copy_file_from_folder(filename,path_to_folder):
    import os
    import shutil as sh
    path = path_to_folder+filename
    folder = os.getcwd()
    sh.copy(path,folder)

# copy file (if exists)
def copy_file(src,dst):
    import os
    import shutil as sh
    if os.path.exists( src ):
	sh.copyfile(src,dst)    
    
# delete file (if exists)
def delete_file(filename):
    import os
    if os.path.exists( filename ):
	os.remove(filename)    
    
# rename file by adding a label at the end to the file name
def rename_file(old_name,label):
    import os
    import os.path as pth
    root, ext = pth.splitext(old_name)
    new_name = root+str(label)+ext
    os.rename(old_name,new_name)
    return new_name

# read the last line of a file
def read_last_line(filename):
    import os
    file_name = open(filename,'r')
    if os.path.exists( filename ):
        string = file_name.readlines()[-1]
    return string

# read all the lines of a file and return a 
# user requested one
def read_user_line(filename,linenumber=1):
    import os
    file_name = open(filename,'r')
    if os.path.exists( filename ):
        string = file_name.readlines()
    return string[linenumber].rstrip('\n')     
    
# ================================ #
# Folders management
# ================================ #

# create folder
def create_folder(folder_name):
    import os
    os.mkdir(folder_name)
    
# delete folder
def delete_folder(folder_name):
    import os
    import shutil as sh
    if os.path.exists( folder_name ):
	if os.listdir(folder_name) == []:
	    os.rmdir(folder_name)
	else:
	    sh.rmtree(folder_name)

# ================================ #
# Strings management
# ================================ #
	
# Write a string in a file 
def write_string (filein,string):
    file_name = open(filein,'a')
    file_name.write(string)
    file_name.close()

# Split a selected string in a file according to input character, 
# being parts the number of parts the string has to be splitted
def split_string (filein,string,character,parts):
    file_name = open(filein,'r+')
    strings = list(range(parts))
    for line in file_name: 
	if line.startswith(string):
	    strings[:] = line.split(character,parts-1)
    file_name.close()
    return strings[:]

# Extract a single string (identified by a given label) from a file filein
def extract_sg_string (filein,label,character,parts,strnum):
    strings    = list( range(parts) )
    strings[:] = split_string (filein,label,character,parts)
    string     = strings[strnum-1]
    return string

# Extract multiple (nstrings) strings (identified by given labels) from a file filein
def extract_mp_string (filein,nstrings,labels,character,parts,strnum):
    string     = list( range(nstrings) )
    for i in range( 0, nstrings ):
        string[i] = extract_sg_string (filein,labels[i],character,parts,strnum)
    return string[:]    
    
# Create a new fileout replacing the string oldstr of the given file
# filein by a string newstr
def replace_string (filein,fileout,instr,oldstr,newstr):
    file_old = open(filein, 'r')
    file_new = open(fileout,'w')
    for line in file_old:
        if line.startswith(instr):
	    new_line = line.replace(oldstr,newstr)
	    file_new.write(new_line)
	else:
	    file_new.write(line)
    file_old.close()
    file_new.close()
        
# Check a string in filein and append it in fileout 
def copy_string (filein,fileout,string):
    file_name_in  = open(filein,'r')
    file_name_out = open(fileout,'a')	
    for line in file_name_in:
	if line.startswith(string):
	    file_name_out.write(str(line))
    file_name_in.close()
    file_name_out.close()
    
# Find a character in a string
def find_char (string, character):
    Index = 0
    while Index < len(string):
        if string[Index] == character:
            return 0
        Index = Index + 1
    return -1

# Cut blank spaces present in a string
def cut_blank_space (string):
    string = string.replace(' ', '')
    return string
    
# ================================ #
# Data files management
# ================================ #

# Read a data file containing floating/integer arrays
# NOTE: 
# Each row of the textfile must have the 
# same number of values 

def read_data (fname,typ=float,delim='',skiph=0,skipf=0):
    data = np.genfromtxt(fname, dtype=typ, delimiter=delim,  skip_header=skiph, \
                         skip_footer=skipf, defaultfmt=None )
    return data
    
def write_data (fname,data,frmt='%.12e',delim=' ',newl='\n',skiph=0,skipf=0,header='\n'):
    file_new = open(fname,'w')
    for i in range( 0, skiph):
        write_string (fname,header)
    f = open(fname,'a')
    np.savetxt(f, data, fmt=frmt, delimiter=delim, newline=newl)  # v. 1.6
    #np.savetxt(fname, data, fmt=frmt, delimiter=delim, newline='\n', header=head, footer=foot)  # v. 1.7
    
