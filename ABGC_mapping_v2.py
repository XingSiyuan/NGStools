# single-sample oriented version of the new ABGC mapping pipeline.
# very early, rough draft - handle with care: many components are still to be validated.
# Note that access to sample database (MySQL) and access to primary data is required for the script to run! 
# Hendrik-Jan Megens, 29-06-2013
# Animal Breeding & Genomics Centre
# Wageningen University
# example usage:
# % python3 ABGC_mapping_v2.py -i LW22F08 -a /media/InternBkp1/repos/ABGSA/ -r /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa
# This will create a file 'runLW22F08.sh', which can be submitted using SGE (qsub)
# for more information:
# % python3 ABGC_mapping_v2.py -h

import argparse
import sys
import os
import re
import mysql.connector
import gzip

parser = argparse.ArgumentParser( description='creates run file for automated trimming, mapping of fastq archives')
parser.add_argument("-i", "--individual_name", help="name of individual to be mapped", nargs=1)
parser.add_argument("-a", "--path_to_abgsa", help="/path/to/abgsa/", nargs=1)
parser.add_argument("-r", "--path_to_reference_fasta", help="/path/to/reference/ref.fa", nargs=1)
parser.add_argument("-t", "--number_of_threads", help="number of threads to be used by aligner", nargs=1)
parser.add_argument("-m", "--mapper", help="mapping method < bwa-mem | bwa-aln | mosaik >", nargs=1)
parser.add_argument("-d", "--dedup_method", help="dedup method < samtools | picard >", nargs=1)

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
   stmt_select = "select ABG_individual_id, archive_name, lane_names_orig from ABGSAschema_main where ABG_individual_id = '"+individual+"' order by lane_names_orig"
   cursor.execute(stmt_select)
   for row in cursor.fetchall():
      output.append([row[1],row[2]])
   for archive in output:
      yield archive

def qsub_headers():
   qf.write('#!/bin/bash'+'\n')
   qf.write('#$ -cwd'+'\n')
   qf.write('#$ -S /bin/sh'+'\n')
   qf.write('#$ -l h_vmem=10G'+'\n')

def prepare_temp_fq_files(abgsa,archive_dir,filenm,tempdir):
   qf.write("gunzip -c "+abgsa+archive_dir+'/'+filenm+" | sed 's/ /#/' | pigz >"+tempdir+filenm+'\n')
   return tempdir+filenm

def trim(tempdir,seqfiles,offset):
   qf.write('# quality trimming of reads by sickle'+'\n')
   stub1=seqfiles[1].replace('.gz','')
   stub2=seqfiles[2].replace('.gz','')
   qf.write('sickle pe -f '+seqfiles[1]+' -r '+seqfiles[2]+' -o '+stub1+'.tr -p '+stub2+'.tr -s '+stub1+'.singles.tr -l 45 -t '+offset+'\n')
   qf.write('pigz '+stub1+'.tr'+'\n')
   qf.write('pigz '+stub2+'.tr'+'\n')
   if offset == 'illumina':
      qf.write('# since sequences have offset +64 we need to convert to sanger (offset +33)'+'\n')
      qf.write('python convert_ill_to_sang.py '+stub1+'.tr.gz | gzip -c >'+stub1+'.tr.sa.gz'+'\n')
      qf.write('python convert_ill_to_sang.py '+stub2+'.tr.gz | gzip -c >'+stub2+'.tr.sa.gz'+'\n')
      qf.write('rm '+stub1+'.tr.gz'+'\n')
      qf.write('rm '+stub2+'.tr.gz'+'\n')
      qf.write('mv '+stub1+'.tr.sa.gz '+stub1+'.tr.gz'+'\n')
      qf.write('mv '+stub2+'.tr.sa.gz '+stub2+'.tr.gz'+'\n')
   seqfiles={}
   seqfiles[1]=stub1+'.tr.gz'
   seqfiles[2]=stub2+'.tr.gz'
   return seqfiles

def map_bwa_mem(bwapath, samtoolspath,archive_dir,index,ref,tempdir,seqfiles,sample,numthreads):
   # BWA-mem is a new algorithm, we need to consider if this is suitable
   qf.write('# maping using the bwa-mem algorithm, including sorting of bam'+'\n')
   qf.write("echo 'start mapping using BWA-mem algorithm'"+'\n')
   qf.write(bwapath+'bwa mem -t '+str(numthreads)+' -R '+"'"+r'@RG\tID:'+archive_dir+'_'+index+r'\tSM:'+sample+r"' "+ref+' '+seqfiles[1]+' '+seqfiles[2]+' >'+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sam'+'\n')
   qf.write(samtoolspath+'samtools view -Shb -q 10 '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sam'+' > '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.bam'+'\n')
   qf.write('rm '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sam'+'\n')
   qf.write("echo 'start sorting'"+'\n')
   qf.write(samtoolspath+'samtools sort '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.bam '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sorted'+'\n')
   qf.write('rm '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.bam'+'\n')
   return tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sorted.bam'
 
def map_bwa_aln(bwapath, samtoolspath,archive_dir,index,ref,tempdir,seqfiles,sample,numthreads):
   #taken from old BWA aligning pipeline from HJM
   stub1=seqfiles[1].replace('.gz','')
   stub2=seqfiles[2].replace('.gz','')
   qf.write('# maping using the bwa-aln algorithm, including sorting of bam'+'\n')
   qf.write("echo 'start mapping using BWA-aln algorithm'"+'\n')
   qf.write(bwapath+'bwa aln -n 0.07 -t '+str(numthreads)+' '+ref+' '+seqfiles[1]+'  >'+stub1+'.sai'+'\n')
   qf.write(bwapath+'bwa aln -n 0.07 -t '+str(numthreads)+' '+ref+' '+seqfiles[2]+'  >'+stub2+'.sai'+'\n')
   qf.write(bwapath+'bwa sampe -P '+ref+r" -r '@RG\tID:"+archive_dir+'_'+index+r'\tSM:'+sample+r"' "+stub1+'.sai '+stub2+'.sai '+seqfiles[1]+' '+seqfiles[2]+' | '+samtoolspath+'samtools view -q 20 -Suh - | '+samtoolspath+'samtools sort -m 5000000000 - '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'PE2.sorted'+'\n')
   return tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'PE2.sorted'

def map_Mosaik(mosaikref, mosaikjump, archive_dir, index,tempdir,seqfiles,sample,numthreads):
   #taken from old Mosaik alignment pipeline from HJM
   qf.write('# maping using Mosaik, including sorting of bam'+'\n')
   qf.write("echo 'start mapping using Mosaik'"+'\n')
   qf.write('MosaikBuild -q '+seqfiles[1]+' -q2 '+seqfiles[2]+' -out '+tempdir+sample+'-'+index+'.dat -st sanger'+'\n');
   qf.write('MosaikAligner -in '+tempdir+sample+'-'+index+'.dat -out '+tempdir+'aln-'+sample+'-'+index+'_build10.dat -ia '+mosaikref+' -hs 15 -mmp 0.07 -m all -mhp 10 -p '+str(numthreads)+' -act 20 -j '+mosaikjump+'\n')

   qf.write('MosaikSort -in '+tempdir+'aln-'+sample+'-'+index+'_build10.dat -out '+tempdir+'aln-'+sample+'-'+index+'_build10-sorted.dat'+'\n')

   qf.write('MosaikText -in '+tempdir+'aln-'+sample+'-'+index+'_build10-sorted.dat -bam '+tempdir+'aln-'+sample+'-'+index+'_build10-sorted.bam'+'\n')
   qf.write('rm '+tempdir+'*.dat'+'\n')
   return tempdir+'aln-'+sample+'-'+index+'_build10-sorted.bam'

def merge_bams(samtoolspath, bams,sample,tempdir):
   # based on old Mosaik/BWA alignment pipeline from HJM
   if len(bams)>1:
      qf.write('#multiple bam files --> do merge'+'\n')
      make_new_header(bams,samtoolspath,tempdir,sample)
      stub = samtoolspath+'samtools merge '+tempdir+'tmpmerged'+sample+'.bam '
      for bam in bams:
        stub = stub+bam+' '
      qf.write(stub+'\n')
      qf.write(samtoolspath+'samtools reheader '+tempdir+'newheader.txt '+tempdir+'tmpmerged'+sample+'.bam >'+tempdir+sample+'_rh.bam'+'\n')
      qf.write('rm '+tempdir+'tmpmerged'+sample+'.bam'+'\n')
   else:
      qf.write('#only one bam file, no need for merging'+'\n')
      qf.write('mv '+bams[0]+' '+tempdir+sample+'_rh.bam'+'\n')
   return tempdir+sample+'_rh.bam'

def make_new_header(bams, samtoolspath,tempdir,sample):
   counter=1
   for bam in bams:
      if counter==1:
         qf.write(samtoolspath+'samtools view -H '+bam+r" | sed 's/SM:unknown/SM:"+sample+r"/' >"+tempdir+'newheader.txt'+'\n')
      else:
         qf.write(samtoolspath+'samtools view -H '+bam+r" | sed 's/SM:unknown/SM:"+sample+r"/' | grep @RG >>"+tempdir+'newheader.txt'+'\n')
      counter+=1

def dedup_picard(samtoolspath,picardpath,bam):
   # based on Qingyuan's pipleline
   qf.write("# dedup using Picard" +'\n')
   bamstub=bam.replace('.bam','')
   qf.write("echo 'dedupping using picard MarkDuplicates'"+'\n')
   qf.write('java7 -Xmx4g -jar '+picardpath+'MarkDuplicates.jar ASSUME_SORTED=true REMOVE_DUPLICATES=true INPUT='+bam+' OUTPUT='+bamstub+'.dedup_pi.bam METRICS_FILE='+bamstub+'.dedup.metrics'+'\n')
   qf.write(samtoolspath+'samtools sort '+bamstub+'.dedup_pi.bam '+bamstub+'.dedup_pi.sorted'+'\n')
   qf.write('rm '+bamstub+'.dedup_pi.bam'+'\n')
   qf.write('mv '+bamstub+'.dedup_pi.sorted.bam '+bamstub+'.dedup_pi.bam'+'\n')
   # consider removing original bam file
   return bamstub+'.dedup_pi.bam'

def dedup_samtools(samtoolspath,bam):
   bamstub=bam.replace('.bam','')
   qf.write("# dedup using samtools"+'\n')
   qf.write("echo 'dedupping using samtools'"+'\n')
   qf.write(samtoolspath+'samtools dedup '+bamstub+'.bam '+bamstub+'.dedup_st.bam'+'\n')
   # consider removing original bam file
   return bamstub+'.dedup_st.bam'
 
def re_align(tempdir,sample,bam,ref,GATKpath):
   # based on Qingyuan's pipeline
   bamstub=bam.replace('.bam','')
   qf.write("java7 -jar "+GATKpath+"GenomeAnalysisTK.jar -T RealignerTargetCreator -R "+ref+" -I "+bam+" -o "+bamstub+".reA.intervals"+'\n')

   qf.write("java7 -jar "+GATKpath+'GenomeAnalysisTK.jar -T IndelRealigner -R '+ref+' -I '+bam+' -targetIntervals '+bamstub+'.reA.intervals -o '+bamstub+'.reA.bam' +'\n')
   qf.write(samtoolspath+'samtools sort '+bamstub+'.reA.bam '+bamstub+'.reA.sorted'+'\n')
   qf.write('rm '+bamstub+'.reA.bam'+'\n')
   qf.write('mv '+bamstub+'.reA.sorted.bam '+bamstub+'.reA.bam'+'\n')
   # consider removing original bam file
   return bamstub+'.reA.bam'

def recalibrate():
   # based on Qingyuan's pipeline
   bamstub=bam.replace('.bam','')
   qf.write('java -jar '+GATKpath+'GenomeAnalysisTK.jar -T BaseRecalibrator -R '+ref+' -I '+bam+' -knownSites '+dbSNPfile+' -o '+bamstub+'.recal.grp'+'\n')

   qf.write('java -jar '+GATKpath+'GenomeAnalysisTK.jar -T PrintReads -R '+ref+' -I '+bam+' -BQSR '+bamstub+'.recal.grp -o '+bamstub+'.recal.bam'+'\n')
   # consider removing original bam file
   return bamstub+'.recal.bam'

def variant_calling_GATK():
   bamstub=bam.replace('.bam','')
   # in progress   

def variant_calling_mpileup(tempdir,bam,sample,ref,samtoolspath,maxfilterdepth,minfilterdepth):
   # based on Qingyuan's pipeline
   qf.write("# variant calling using the mpileup function of samtools"+'\n')
   bamstub=bam.replace('.bam','')
   bcftoolspath=samtoolspath+'bcftools/'
   qf.write(samtoolspath+"samtools mpileup -C50 -ugf "+ref+' '+bam+' | '+bcftoolspath+'bcftools view -bvcg -| '+bcftoolspath+'bcftools view - | perl '+bcftoolspath+'bcftools/vcfutils.pl varFilter -D '+maxfilterdepth+' -d '+minfilterdepth+' >'+bamstub+'.var.flt.vcf'+'\n')

def variant_calling_pileup(samtoolspath_v12,tempdir,sample, ref,bam):
   # based on HJM's old Mosaik/BWA pipeline
   bamstub=bam.replace('.bam','')
   qf.write('# old-school variant calling using the pileup algorithm'+'\n')
   qf.write("echo 'old-school variant calling using the pileup algorithm'"+'\n')
   qf.write(samtoolspath_v12+'samtools view -u '+bam+' | '+samtoolspath_v12+'samtools pileup -vcf '+ref+' - >'+bamstub+'_vars-raw.txt'+'\n')
   qf.write(r'VAR=`cat '+bamstub+'_vars-raw.txt'+r" | cut -f8 | head -100000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`"+'\n')
   qf.write('let VAR=2*VAR'+'\n')
   qf.write('echo "max depth is $VAR"'+'\n')
   qf.write(samtoolspath_v12+'misc/samtools.pl varFilter -D$VAR '+bamstub+r"_vars-raw.txt | awk '($3=="+'"*"&&$6>=50)||($3!="*"&&$6>=20)'+r"' >"+bamstub+'_vars-flt_final.txt'+'\n')

def create_shell_script(sample,abgsa,ref,mapper,numthreads):
   # print qsub header lines
   qsub_headers()

   # set a bunch of variables and paths - consider doing by config-file
   tempdir = 'tmp'+sample
   bwapath='/opt/bwa/bwa-0.7.5a/'
   samtoolspath='/opt/samtools/samtools-0.1.19/'
   samtoolspath_v12='/opt/samtools/samtools-0.1.12a/'
   picardpath='/opt/picard/picard-tools-1.93/'
   GATKpath='/opt/GATK/GATK2.6/'
   mosaikref='/path/to/mosaik/ref.dat'
   mosaikjump='/path/to/mosaikjump/ref.j15'
   maxfilterdepth=20
   minfilterdepth=4

   qf.write('mkdir '+tempdir+'\n')
   tempdir=tempdir+'/'
    
   # get sequence info from database
   archives=get_info_from_db(sample)

   count=0
   bams=[]

   # preparing fq, trimming, mapping, in a loop,
   # per two gzipped fq files
   for archive in archives:
      seqfiles={}
      count+=1
      qf.write("# archive number "+str(count)+": "+archive[0]+'\n')
      archive_dir=archive[0]
      seqfiles[1]=archive[1]
      offset=check_illumina_or_sanger(abgsa+archive_dir+'/'+seqfiles[1])
      seqfiles[1]=prepare_temp_fq_files(abgsa,archive_dir,seqfiles[1],tempdir)
      archive = next(archives)
      seqfiles[2]=archive[1]
      seqfiles[2]=prepare_temp_fq_files(abgsa,archive_dir,seqfiles[2],tempdir)
      seqfiles=trim(tempdir,seqfiles,offset)
      if mapper == 'bwa-mem':
         bams.append(map_bwa_mem(bwapath, samtoolspath,archive_dir,str(count),ref,tempdir,seqfiles,sample,numthreads))
      elif mapper == 'bwa-aln':
         bams.append(map_bwa_aln(bwapath, samtoolspath,archive_dir,str(count),ref,tempdir,seqfiles,sample,numthreads))
      elif mapper == 'mosaik':  
         bams.append(map_Mosaik(mosaikref, mosaikjump, archive_dir,str(count),tempdir,seqfiles,sample,numthreads))
   qf.write("#number of bams: "+str(len(bams))+'\n') 

   # merge and reheader bam files  
   bam=merge_bams(samtoolspath, bams,sample,tempdir)
   
   # further optimization of bam files
   if dedup == 'samtools':
      bam=dedup_samtools(samtoolspath,bam)
   elif dedup == 'picard':
      bam=dedup_picard(samtoolspath,picardpath,bam)

   # variant calling
   variant_calling_pileup(samtoolspath_v12,tempdir,sample,ref,bam)
   variant_calling_mpileup(tempdir,bam,sample,ref,samtoolspath,str(maxfilterdepth),str(minfilterdepth))

if __name__=="__main__":
   # initialize db cursor
   db = mysql.connector.Connect(user='anonymous',host='localhost',database='ABGSAschema', password='anonymous')
   cursor = db.cursor()
   
   # get command line options
   args = parser.parse_args()
   individual=args.individual_name[0]
   abgsa = args.path_to_abgsa[0]
   mapper=args.mapper[0]
   numthreads=args.number_of_threads[0]
   ref = args.path_to_reference_fasta[0]
   dedup=args.dedup_method[0]
   
   # open qsub-file (qf)
   qf=open('run'+individual+'.sh','w')
   # invoke master subroutine
   create_shell_script(individual,abgsa,ref,mapper,numthreads)
   qf.close()

   # optional submitting job
   #os.sys('qsub -q all.q run'+individual+'.sh') 
