import argparse
import sys
import os
import re
import mysql.connector
import gzip

parser = argparse.ArgumentParser( description='create a bunch of parameters of primary data files')
#parser.add_argument("-a", "--archive_name", help="sequence archive name", nargs=1)
parser.add_argument("-p", "--path_to_abgsa", help="/path/to/sequence/archive", nargs=1)


def do_gc(filenm):
   gc= float(os.popen('gunzip -c '+filenm+' | head -4000000 | perl GCcontfq.pl ').read().rstrip().split(' ')[1])
   #print(gc)
   if gc > 0:
      return gc
   else:
      return 0

def get_info_from_db():
   output=[]
   stmt_select = "select ABG_individual_id, archive_name, lane_names_orig,md5sum_gzip from ABGSAschema_main where lane_names_orig not in ( select lane_names_orig from gc_cont ) order by archive_name,lane_names_orig"
   #stmt_select = "select ABG_individual_id, archive_name, lane_names_orig,md5sum_gzip from ABGSAschema_main where archive_name = '"+archive_name+"' order by lane_names_orig"
   cursor.execute(stmt_select)
   for row in cursor.fetchall():
      output.append([row[0],row[1],row[2],row[3]])
   for archive in output:
      yield archive

def insert(values):
   stmt_insert  = ("INSERT INTO gc_cont (lane_names_orig, gc_cont) VALUES (%s, %s)")
   cursor.execute(stmt_insert,values)
   db.commit()
# example: http://dev.mysql.com/doc/refman/5.6/en/myconnpy_example_cursor_transaction.html
#CREATE TABLE `fqfile_attributes` (`id` INT PRIMARY KEY AUTO_INCREMENT,`lane_names_orig` VARCHAR(120) default null, `md5check` varchar(120) default null, filesize bigint(15) default null, numlines_in_gzip bigint(15) default null, qval_offset varchar(10) default null, maxlength_seq int(6) default null, totnumbases bigint(15) default null, totnumreads int(12) default null, firstreadname varchar(120) default null, lastreadname varchar(120) default null, `ts_create` TIMESTAMP DEFAULT CURRENT_TIMESTAMP)


if __name__=="__main__":
   # initialize db cursor
   db = mysql.connector.Connect(user='ABGSAroot',host='localhost',database='ABGSAschema', password='abgsaroot@1234')
   cursor = db.cursor()
   # abgsaroot@1234 = 39970aea463b78d7
   # get command line options
   args = parser.parse_args()
#  archive_name=args.archive_name[0]
   abgsapath = args.path_to_abgsa[0]
#  archives=get_info_from_db(archive_name)
   archives=get_info_from_db()
   for archive in archives: 
      filenm=abgsapath+archive[1]+'/'+archive[2] 
      gc = do_gc(filenm)
      print('{}\t{}\t{}\t{}\t{}'.format(archive[0],archive[1],archive[2],gc,archive[3]))
      if gc > 0:
         insert((archive[2],gc))
      else:
         print('warning for '+archive[2]+': nothing inserted')
   cursor.close()
   db.close()
