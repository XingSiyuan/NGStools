import argparse
import sys
import os
import re
import mysql.connector
import gzip

parser = argparse.ArgumentParser( description='create a bunch of parameters of primary data files')
parser.add_argument("-a", "--archive_name", help="sequence archive name", nargs=1)
parser.add_argument("-p", "--path_to_abgsa", help="/path/to/sequence/archive", nargs=1)

def next_sequence_gzip(filename):
    try:
        fileh = gzip.open(filename)
        line = fileh.readline()[:-1].decode('utf-8')
        lines=[]
        while line:
           if line and line[0] == '@':
              lines.append(line)
              lines.append(fileh.readline().decode('utf-8'))
              lines.append(fileh.readline().decode('utf-8'))
              lines.append(fileh.readline().decode('utf-8'))

           yield lines
           line = fileh.readline()[:-1].decode('utf-8')
    finally:
        fileh.close()

def check_illumina_or_sanger(file_name):
    maxQ=0
    seqs = next_sequence_gzip(file_name)
    offset='sanger'
    maxlength=0;
    for i in range(10000):
      seq = next(seqs)
      qs = seq[3][0:-1]
      if len(qs)>maxlength:
         maxlength=len(qs)
      for q in qs:
        if (ord(q)-33)>maxQ:
           maxQ=ord(q)-33
    if maxQ>41:
        offset='illumina'
    return offset,maxlength

def check_numbases(file_name):
    allbases=0
    allreads=0
    firstreadname='noname'
    lastreadname='noname'
    seqs = next_sequence_gzip(file_name)
    for seq in seqs:
       lastreadname = seq[0][0:-1]
       if allreads==0:
          firstreadname=lastreadname
       sq = seq[1][0:-1]
       allbases+=len(sq)
       allreads+=1
    return allbases,allreads,firstreadname,lastreadname

def do_md5check(md5check,md5original,filenm):
   if md5check:
      if md5original != os.popen('md5sum '+filenm).read().split()[0]:
         raise SystemExit
      else:
         print('md5 ok: '+md5original)

def get_numlines_in_zipfile(filenm):
   numlines = os.popen('gunzip -c '+filenm+' | wc -l').read().split()[0]
   return numlines
    
def get_filesize(filenm):
   filesize=os.path.getsize(filenm)
   return filesize

def get_info_from_db(archive_name):
   output=[]
   stmt_select = "select ABG_individual_id, archive_name, lane_names_orig,md5sum_gzip from ABGSAschema_main where archive_name = '"+archive_name+"' order by lane_names_orig"
   cursor.execute(stmt_select)
   for row in cursor.fetchall():
      output.append([row[0],row[1],row[2],row[3]])
   for archive in output:
      yield archive

def get_filesize_from_db(fqname):
   output=[]
   stmt_select = "select filesize from filesize where lane_names_orig = '"+fqname+"'"
   #print(stmt_select)
   cursor.execute(stmt_select)
   for row in cursor.fetchall():
      output.append([row[0]])
   return output[0]

def insert(values):
   stmt_insert  = ("INSERT INTO fqfile_attributes (lane_names_orig, md5sum, filesize, numlines_in_gzip, qval_offset, maxlength_seq, totnumbases, totnumreads, firstreadname, lastreadname) VALUES (%s, %s, %s, %s, %s)")
   emp_no = cursor.lastrowid
   cursor.execute(stmt_select)
   db.commit()
# example: http://dev.mysql.com/doc/refman/5.6/en/myconnpy_example_cursor_transaction.html

if __name__=="__main__":
   # initialize db cursor
   db = mysql.connector.Connect(user='XXXX',host='localhost',database='ABGSAschema', password='XXXXXXX')
   cursor = db.cursor()
   
   # get command line options
   args = parser.parse_args()
   archive_name=args.archive_name[0]
   abgsapath = args.path_to_abgsa[0]
   archives=get_info_from_db(archive_name)
   for archive in archives: 
      filenm=abgsapath+archive[1]+'/'+archive[2] 
      filesize = get_filesize(filenm)
      filesize_indb = get_filesize_from_db(archive[2])
      numlines_in_gzip = get_numlines_in_zipfile(filenm)
      allbases,allreads,firstreadname,lastreadname = check_numbases(filenm)
      offset,maxlength = check_illumina_or_sanger(filenm)
      print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(archive[0],archive[1],archive[2],archive[3],filesize,filesize_indb,numlines_in_gzip,offset,maxlength,allbases,allreads,firstreadname,lastreadname))
   cursor.close()
   db.close()
