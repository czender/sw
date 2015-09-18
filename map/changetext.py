#!/usr/bin/env python
#
# chg_txt.py
#
# Replacing all chosen text entries in all files
# in a directory
#
# Copyright Department of Geophysics
# University of Oslo
#
#

import sys
import os
import re
import string

def fixup(ifile, ofile, old_text, new_text):

   #print ifile, ofile
   rfile = open(ifile,'r')
   wfile = open(ofile,'w')

   while 1: 
      line = rfile.readline()
      if not line:
         break
      ofs = line.rfind(old_text)
      if ofs != -1:
         out_line = line.replace (old_text, new_text)
         print line
         print 'Outline: ', out_line
         wfile.write(out_line)
      else:
         wfile.write(line)


   # Close files
   rfile.close()
   wfile.close()


def find_files(directory, old_text, new_text):

   if directory != '':
      os.chdir(directory)

   os.system('ls -l > /tmp/thisdir')
   hlist = open('/tmp/thisdir','r')
   contents = hlist.readlines()
   for i in contents:
      #print i[:-1]
      if i[:3] == '-rw':
         str = string.split(i)
         l = len(str)
         #print str[l-1]
         ifile = str[l-1]+'.ORG'
         ofile = str[l-1]
         if ofile[-4:] == '.F90' or ofile[-3:] == '.sh':
            cmd = 'cp '+ofile+' '+ifile
            print 'Copying ', cmd
            os.system(cmd)
            fixup(ifile, ofile, old_text, new_text)



if __name__ == '__main__':
   l = len(sys.argv)
   #print l
   if l < 3:
      print 'Usage: chg_txt directory[optional] old_text new_text '
      sys.exit()

   directory = ''

   if l == 4:
      find_files(sys.argv[1], sys.argv[2], sys.argv[3])
   else:
      find_files(directory, sys.argv[1], sys.argv[2])
