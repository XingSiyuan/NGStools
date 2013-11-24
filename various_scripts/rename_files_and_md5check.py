import argparse
import sys
import os
import re

# python3 rename_files.py -t transtable.txt
# where 'transtable.txt' contains the following (first 5lines):
# 53258b1a7ca265bc04810d5fb73da5c4        bov-wur-58      run0241 bov-wur-58_AGTCAA_L005_R1_001.fastq.gz
#404ad16f1cbc6607e55fed98f16c971f        bov-wur-58      run0241 bov-wur-58_AGTCAA_L005_R2_001.fastq.gz
#6989f5a3aea8ff4f0ac4b4cb8c3d62dd        bov-wur-58      run0247 bov-wur-58_AGTCAA_L001_R1_001.fastq.gz
#f3ceedc008087db21527c8e8ecde4577        bov-wur-58      run0247 bov-wur-58_AGTCAA_L001_R2_001.fastq.gz

parser = argparse.ArgumentParser( description='Renames a bunch of files, creates dirs, and reorganizes files')
parser.add_argument("-t", "--translatetable", help="translatetable", nargs=1)

def rename():
  translist = translatetable(transfile)
  newlist = fq_filelist(translist)
  for fq in newlist:
    oldname = fq[1]+'/'+fq[2]+'/'+fq[3]
    newname = fq[1]+'/'+fq[4]
    print(oldname+' ---> '+newname)
    os.rename(oldname,newname)
    do_md5check(fq[0],newname)
    
def do_md5check(md5original,filenm):
      newmd5= os.popen('md5sum '+filenm).read().split()[0]
      if md5original != newmd5:
         print('WARNING! - conflict in md5sums: original --> '+md5original+' :: newmd5 --> '+newmd5)
      else:
         print('md5 ok: original --> '+md5original+' :: newmd5 --> '+newmd5)

def fq_filelist(flist):
  newlist=[]
  for item in flist:
     parts = item[3].split('_',1)
     newname = parts[0]+'_'+item[2]+'_'+parts[1]
     newlist.append(item+[newname])
  return newlist

def translatetable(trans_file):
  flist=[]
  with open(trans_file) as trans:
    for l in trans.readlines():
      l=l.rstrip().split()
      flist.append([l[0],l[1],l[2],l[3]])
  return flist

args = parser.parse_args()
transfile=args.translatetable[0]

if __name__=="__main__":
  rename()
