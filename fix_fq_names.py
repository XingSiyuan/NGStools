# Hendrik-Jan Megens, 29-06-2013
# Animal Breeding & Genomics Centre
# Wageningen University
# NOTE: currently this is a python2.7 script
# need to modify for python3 (utf-8 encoding)
import sys
import gzip
import re
def next_sequence_gzip(filename):
   try:
      file = gzip.open(filename)
      line = file.readline()
	
      while line:
         if line and line[0] == '@':
            line1 = line
            line2 = file.readline()
            line3 = file.readline()
            line4 = file.readline()
         yield (line1,line2,line3,line4)
         line = file.readline()
   finally:
      file.close()

def fix_fq_name(file):
   seqs = next_sequence_gzip(file)
   for seq in seqs:
       fqname=seq[0][0:-1]
       fqname=re.subn(' |\/','#',fqname)[0]
       qs = seq[3][0:-1]
       print fqname+'\n'+seq[1]+seq[2]+qs

if __name__=="__main__":
   fix_fq_name(sys.argv[1])

