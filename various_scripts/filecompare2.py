import argparse
import sys
import os
import re
import mysql.connector
import gzip

# python3 filecompare2.py -1 /lustre/nobackup/WUR/ABGC/shared/Pig/Mapping_results/ -2 /lustre/backup/WUR/ABGC/shared/ABGC_datastore/SequenceData/Pig/Mapping_results/


parser = argparse.ArgumentParser( description='Check if files in path1 are redundantly present in path2, and remove path1 files if True')
parser.add_argument("-1", "--path_to_first", help="/path/to/", nargs=1)
parser.add_argument("-2", "--path_to_second", help="/path/to/", nargs=1)

def get_filesize(filenm):
   filesize=os.path.getsize(filenm)
   return filesize

def do_md5check(filenm):
   md5 = os.popen('md5sum '+filenm).read().split()[0]
   return md5

def return_flist(path):
   myfiles=os.walk(path)
   finallist=[]
   for myfile in myfiles:
      #print(myfile)
      if myfile[1]==[]:
         for mfile in myfile[2]:
            #print(mfile)
            finallist.append(['/'.join([myfile[0],mfile]),'/'.join([myfile[0],mfile]).replace(path,'')])
   return finallist

if __name__=="__main__":
   args = parser.parse_args()
   path1 = args.path_to_first[0]
   path2 = args.path_to_second[0]

   firstlist=return_flist(path1)
   secondlist=return_flist(path2)
   seconddict={mfile[1]:mfile[0] for mfile in secondlist}
   for mfile in firstlist:
      print("++++++++++++++++++++++")
      if mfile[1] in seconddict.keys():
         
         if get_filesize(mfile[0]) == get_filesize(seconddict[mfile[1]]):
            print("filesizes match: ",mfile[0],get_filesize(mfile[0]),seconddict[mfile[1]],get_filesize(seconddict[mfile[1]]))
            md5sum1=do_md5check(mfile[0])
            md5sum2=do_md5check(seconddict[mfile[1]])
            if md5sum1 == md5sum2:
               print("md5sums match: ",mfile[0],md5sum1,seconddict[mfile[1]],md5sum2)
               print("***REMOVING***: ",mfile[0])
               os.remove(mfile[0])
            else:
               print("WARNING: md5sums do not match! ",mfile[0],md5sum1,seconddict[mfile[1]],md5sum2)
            
         else:
            print("WARNING: file sizes do not match! ",mfile[0],get_filesize(mfile[0]),seconddict[mfile[1]],get_filesize(seconddict[mfile[1]]))
         
      else:
         print("WARNING: ",mfile[1], "not present")


