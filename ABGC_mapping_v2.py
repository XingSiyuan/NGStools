import argparse
import sys
import os
import re
import mysql.connector
import gzip

parser = argparse.ArgumentParser( description='creates run file for automated trimming, mapping of fastq archives')
parser.add_argument("-i", "--individual_name", help="name of individual to be mapped", nargs=1)
parser.add_argument("-o", "--outfile", help="output file", nargs=1)

   #db = mysql.connector.Connect(user='hapmapuser',password='hapmapuser@1234',host='localhost',database='pig_hapmap2')

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
           line = fileh.readline()
    finally:
        fileh.close()

def check_illumina_or_sanger(file_name):
    maxQ=0
    seqs = next_sequence_gzip(file_name)
    offset='sanger'
    for i in range(1000):
      seq = next(seqs)
      qs = seq[3][0:-1]
      for q in qs:
        if (ord(q)-33)>maxQ:
           maxQ=ord(q)-33
    if maxQ>41:
        offset='illumina'
    return offset

def get_info_from_db(individual):
   output=[]
   #print(snp_name)
   stmt_select = "select ABG_individual_id, archive_name, lane_names_orig from ABGSAschema_main where ABG_individual_id = '"+individual+"' order by lane_names_orig"
   cursor.execute(stmt_select)
   for row in cursor.fetchall():
      output.append([row[1],row[2]])
   for archive in output:
      yield archive
def qsub_headers():
   print('#!/bin/bash')
   print('#$ -cwd')
   print('#$ -S /bin/sh')
   print('#$ -l h_vmem=10G')

def prepare_temp_fq_files(abgsa,archive_dir,filenm,tempdir):
   print("gunzip -c "+abgsa+archive_dir+'/'+filenm+" | sed 's/ /#/' | pigz >"+tempdir+filenm)

def trim(tempdir,file1,file2,offset):
      print('# quality trimming of reads by sickle')
      stub1=file1.replace('.gz','')
      stub2=file2.replace('.gz','')
      print('sickle pe -f '+tempdir+file1+' -r '+tempdir+file2+' -o '+tempdir+stub1+'.tr -p '+tempdir+stub2+'.tr -s '+tempdir+stub1+'.singles.tr -t '+offset)
      print('pigz '+tempdir+stub1+'.tr')
      print('pigz '+tempdir+stub2+'.tr')
      if offset == 'illumina':
         print('# since sequences have offset +64 we need to convert to sanger (offset +33)')
         print('python convert_ill_to_sang.py '+tempdir+stub1+'.tr.gz | gzip -c >'+tempdir+'/'+stub1+'.tr.sa.gz')
         print('python convert_ill_to_sang.py '+tempdir+stub2+'.tr.gz | gzip -c >'+tempdir+'/'+stub2+'.tr.sa.gz')
         print('rm '+tempdir+stub1+'.tr.gz')
         print('rm '+tempdir+stub2+'.tr.gz')
         print('mv '+tempdir+stub1+'.tr.sa.gz '+tempdir+stub1+'.tr.gz')
         print('mv '+tempdir+'/'+stub2+'.tr.sa.gz '+tempdir+stub2+'.tr.gz')

def map_bwa_mem(bwapath, samtoolspath,archive_dir,index,ref,tempdir,file1,file2,sample):
   print('# maping using the bwa-mem algorithm, including sorting of bam')
   stub1=file1.replace('.gz','')
   stub2=file2.replace('.gz','')
   #print('/opt/bwa/bwa-0.7.5a/bwa mem -t 4 -R '+"'"+r'@RG\tID:'+archive_dir+'_'+index+r'\tSM:1 '+ref+' '+tempdir+'/'+stub1+'.tr.gz '+tempdir+'/'+stub2+'.tr.gz | gzip -3 > aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sam.gz')
   print(bwapath+'bwa mem -t 4 -R '+"'"+r'@RG\tID:'+archive_dir+'_'+index+r'\tSM:'+sample+r"' "+ref+' '+tempdir+stub1+'.tr.gz '+tempdir+stub2+'.tr.gz | '+samtoolspath+'samtools view -Shb -q 10 > '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.bam')
   print(samtoolspath+'samtools sort '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.bam '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sorted')
   print('rm '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.bam')
   return tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sorted.bam'
 
def merge_bams(samtoolspath, bams,sample,tempdir):
   if len(bams)>1:
      print('#multiple bam files --> do merge')
      make_new_header(bams,samtoolspath)
      stub = samtoolspath+'samtools merge '+tempdir+'tmpmerged'+tempdir+sample+'.bam '
      for bam in bams:
        stub = stub+bam+' '
      print(stub)
      print(samtoolspath+'samtools reheader '+tempdir+'newheader.txt '+tempdir+'tmpmerged'+sample+'.bam >'+tempdir+sample+'_rh.bam')
      print('rm tmpmerged'+sample+'.bam')
   else:
      print('#only one bam file, no need for merging')
      print('mv '+tempdir+bams[0]+' '+tempdir+sample+'_rh.bam')

def make_new_header(bams, samtoolspath):
    counter=1
    for bam in bams:
       if counter==1:
          print(samtoolspath+'samtools view -H '+bam+' >newheader.txt')
       else:
          print(samtoolspath+'samtools view -H '+bam+' | grep @RG >>newheader.txt')
       counter+=1

def variant_calling_pileup(samtoolspath_v12,tempdir,sample,filterdepth, ref):
   print('# old-school variant calling using the pileup algortithm')
   print(samtoolspath_v12+'samtools view -u '+tempdir+sample+'_rh.bam | '+samtoolspath_v12+'samtools pileup -vcf '+ref+' - >vars-raw_'+tempdir+sample+'.txt')

   print(samtoolspath_v12+'samtools.pl varFilter -D'+filterdepth+' vars-raw_'+tempdir+sample+r".txt | awk '($3=="+'"*"&&$6>=50)||($3!="*"&&$6>=20)'+r"' >vars-flt_"+tempdir+sample+'-final.txt')

def create_shell_script(sample):
   qsub_headers()
   #abgsa = '/srv/mds01/shared/Sus/ABGSA/';
   abgsa = '/media/InternBkp1/repos/ABGSA/';
   ref = '/media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa';
   #allbam = 'samtools merge merged.bam';
   tempdir = 'tmp'+sample
   bwapath='/opt/bwa/bwa-0.7.5a/'
   samtoolspath='/opt/samtools/samtools-0.1.19/'
   samtoolspath_v12='/opt/samtools/samtools-0.1.12/'
   filterdepth=20
   print('mkdir '+tempdir)
   tempdir=tempdir+'/'
   archives=get_info_from_db(sample)
   count=0
   bams=[]
   for archive in archives:
      count+=1
      print("# archive number "+str(count)+": "+archive[0])
      archive_dir=archive[0]
      file1=archive[1]
      prepare_temp_fq_files(abgsa,archive_dir,file1,tempdir)
      archive = next(archives)
      file2=archive[1]
      prepare_temp_fq_files(abgsa,archive_dir,file2,tempdir)
      offset=check_illumina_or_sanger(abgsa+archive_dir+'/'+file1)
      trim(tempdir,file1,file2,offset)
      bams.append(map_bwa_mem(bwapath, samtoolspath,archive_dir,str(count),ref,tempdir,file1,file2,sample))
   print("#number of bams: "+str(len(bams)))   
   merge_bams(samtoolspath, bams,sample,tempdir)
   variant_calling_pileup(samtoolspath_v12,tempdir,sample,str(filterdepth),ref)


if __name__=="__main__":
  db = mysql.connector.Connect(user='anonymous',host='localhost',database='ABGSAschema', password='anonymous')
  cursor = db.cursor()
  args = parser.parse_args()
  individual=args.individual_name[0]
  alignedfile=args.outfile[0]
  create_shell_script(individual)
#sickle pe -f ABGSA0189/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz -r ABGSA0189/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz -o ABGSA0189_trim/test_R1.tr.fq -p ABGSA0189_trim/test_R2.tr.fq -s ABGSA0189_trim/singles/test.singles.tr.fq -t illumina

