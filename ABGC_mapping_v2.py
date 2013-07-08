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
parser.add_argument("-t", "--number_of_threads", help="number of threads to be used by aligner",type=int, default=1)
parser.add_argument("-m", "--mapper", help="mapping method", type=str, choices=['bwa-mem','bwa-aln','mosaik'], default='bwa-mem')
parser.add_argument("-d", "--dedup_method", help="dedup method", type=str, choices=['samtools','picard'], default='samtools')
parser.add_argument("-c", "--domd5check", help="check md5 integrity of sequence archive against database", action="store_true")

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
   stmt_select = "select ABG_individual_id, archive_name, lane_names_orig,md5sum_gzip from ABGSAschema_main where ABG_individual_id = '"+individual+"' order by lane_names_orig"
   cursor.execute(stmt_select)
   for row in cursor.fetchall():
      output.append([row[1],row[2],row[3]])
   for archive in output:
      yield archive

def do_md5check(md5check,md5original,filenm):
   if md5check:
      if md5original != os.popen('md5sum '+filenm).read().split()[0]:
         raise SystemExit
      else:
         print('md5 ok: '+md5original)

def qsub_headers():
   qf.write('#!/bin/bash'+'\n')
   qf.write('#$ -cwd'+'\n')
   qf.write('#$ -S /bin/bash'+'\n')
   qf.write('#$ -l h_vmem=10G'+'\n')

def prepare_temp_fq_files(abgsamapping_toolpath, abgsa,archive_dir,filenm,tempdir):
   qf.write("python "+abgsamapping_toolpath+"fix_fq_names.py "+abgsa+archive_dir+'/'+filenm+" | pigz >"+tempdir+filenm+'\n')
   return tempdir+filenm

def trim(abgsamapping_toolpath, tempdir,seqfiles,offset):
   qf.write('# quality trimming of reads by sickle'+'\n')
   stub1=seqfiles[1].replace('.gz','')
   stub2=seqfiles[2].replace('.gz','')
   qf.write('sickle pe -f '+seqfiles[1]+' -r '+seqfiles[2]+' -o '+stub1+'.tr -p '+stub2+'.tr -s '+stub1+'.singles.tr -l 45 -t '+offset+'\n')
   qf.write('pigz '+stub1+'.tr'+'\n')
   qf.write('pigz '+stub2+'.tr'+'\n')
   if offset == 'illumina':
      qf.write('# since sequences have offset +64 we need to convert to sanger (offset +33)'+'\n')
      qf.write('python '+abgsamapping_toolpath+'convert_ill_to_sang.py '+stub1+'.tr.gz | gzip -c >'+stub1+'.tr.sa.gz'+'\n')
      qf.write('python '+abgsamapping_toolpath+'convert_ill_to_sang.py '+stub2+'.tr.gz | gzip -c >'+stub2+'.tr.sa.gz'+'\n')
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
   qf.write(bwapath+'bwa mem -t '+str(numthreads)+' -M -R '+"'"+r'@RG\tID:'+archive_dir+'_'+index+r'\tSM:'+sample+r"\tPL:ILLUMINA' "+ref+' '+seqfiles[1]+' '+seqfiles[2]+' >'+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'-pe.sam'+'\n')
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
   qf.write(bwapath+'bwa sampe -P '+ref+r" -r '@RG\tID:"+archive_dir+'_'+index+r'\tSM:'+sample+r"\tPL:ILLUMINA' "+stub1+'.sai '+stub2+'.sai '+seqfiles[1]+' '+seqfiles[2]+' | '+samtoolspath+'samtools view -q 20 -Suh - | '+samtoolspath+'samtools sort -m 5000000000 - '+tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'PE2.sorted'+'\n')
   return tempdir+'aln-'+archive_dir+'-'+sample+'-'+index+'PE2.sorted.bam'

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
      qf.write(samtoolspath+'samtools reheader '+tempdir+'newheader.txt '+bam[0]+' >'+tempdir+sample+'_rh.bam'+'\n')
   qf.write(samtoolspath+'samtools index '+tempdir+sample+'_rh.bam'+'\n')
   return tempdir+sample+'_rh.bam'

def make_new_header(bams, samtoolspath,tempdir,sample):
   counter=1
   for bam in bams:
      if counter==1:
         qf.write(samtoolspath+'samtools view -H '+bam+r" | sed 's/SM:unknown/SM:"+sample+r"/'  | sed 's/PL:sanger/PL:ILLUMINA/' >"+tempdir+'newheader.txt'+'\n')
      else:
         qf.write(samtoolspath+'samtools view -H '+bam+r" | sed 's/SM:unknown/SM:"+sample+r"/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>"+tempdir+'newheader.txt'+'\n')
      counter+=1

def dedup_picard(samtoolspath,picardpath,bam):
   # based on Qingyuan's pipleline
   # currently does not work - too many issues, mostly regarding memory use
   qf.write("# dedup using Picard" +'\n')
   bamstub=bam.replace('.bam','')
   qf.write("echo 'dedupping using picard MarkDuplicates'"+'\n')
   qf.write('java7 -Xmx4g -jar '+picardpath+'MarkDuplicates.jar MAX_RECORDS_IN_RAM=5000000 ASSUME_SORTED=true REMOVE_DUPLICATES=true INPUT='+bam+' OUTPUT='+bamstub+'.dedup_pi.bam METRICS_FILE='+bamstub+'.dedup.metrics'+'\n')
   qf.write(samtoolspath+'samtools sort '+bamstub+'.dedup_pi.bam '+bamstub+'.dedup_pi.sorted'+'\n')
   qf.write('rm '+bamstub+'.dedup_pi.bam'+'\n')
   qf.write('mv '+bamstub+'.dedup_pi.sorted.bam '+bamstub+'.dedup_pi.bam'+'\n')
   qf.write(samtoolspath+'samtools index '+bamstub+'.dedup_pi.bam'+'\n')
   # consider removing original bam file
   # Question: is it really necesary to re-sort? Couldn't find info. Investigate
   return bamstub+'.dedup_pi.bam'

def dedup_samtools(samtoolspath,bam):
   bamstub=bam.replace('.bam','')
   qf.write("# dedup using samtools"+'\n')
   qf.write("echo 'dedupping using samtools'"+'\n')
   qf.write(samtoolspath+'samtools rmdup '+bamstub+'.bam '+bamstub+'.dedup_st.bam'+'\n')
   qf.write(samtoolspath+'samtools index '+bamstub+'.dedup_st.bam'+'\n')
   # consider removing original bam file
   qf.write("cp "+bamstub+'.dedup_st.bam.bai '+bamstub+'.dedup_st.bai'+'\n')
   return bamstub+'.dedup_st.bam'
 
def re_align(samtoolspath,bam,ref,GATKpath):
   # based on Qingyuan's pipeline
   qf.write('# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner'+'\n')
   bamstub=bam.replace('.bam','')
   qf.write("java7 -jar "+GATKpath+"GenomeAnalysisTK.jar -nt "+str(numthreads)+" -T RealignerTargetCreator -R "+ref+" -I "+bam+" -o "+bamstub+".reA.intervals"+'\n')

   qf.write("java7 -jar "+GATKpath+'GenomeAnalysisTK.jar -T IndelRealigner -R '+ref+' -I '+bam+' -targetIntervals '+bamstub+'.reA.intervals -o '+bamstub+'.reA.bam' +'\n')
   #qf.write(samtoolspath+'samtools sort '+bamstub+'.reA.bam '+bamstub+'.reA.sorted'+'\n') # re-sorting does not seem needed?
   #qf.write('rm '+bamstub+'.reA.bam'+'\n')
   #qf.write('mv '+bamstub+'.reA.sorted.bam '+bamstub+'.reA.bam'+'\n')
   #qf.write(samtoolspath+'samtools index '+bamstub+'.reA.bam'+'\n') # re-indexing needed?
   # consider removing original bam file
   # Question 1: is mate information retained?
   # Do we need to add Picard's FixMateInformation?.jar ?
   # Question 2: is het really necesary to re-sort? Maybe only re-index? Couldn't find info. Investigate
   qf.write("cp "+bamstub+'.reA.bai '+bamstub+'.reA.bam.bai'+'\n')
   return bamstub+'.reA.bam'

def recalibrate(GATKpath,samtoolspath,dbSNPfile,ref,bam):
   # based on Qingyuan's pipeline
   qf.write('# Recalibration of BAM using GATK-BaseRecalibrator+PrintReads'+'\n')
   bamstub=bam.replace('.bam','')
   qf.write('java7 -jar '+GATKpath+'GenomeAnalysisTK.jar -nct '+str(numthreads)+' -T BaseRecalibrator -R '+ref+' -I '+bam+' -knownSites '+dbSNPfile+' -o '+bamstub+'.recal.grp'+'\n')
   qf.write('java7 -jar '+GATKpath+'GenomeAnalysisTK.jar -nct '+str(numthreads)+' -T PrintReads -R '+ref+' -I '+bam+' -BQSR '+bamstub+'.recal.grp -o '+bamstub+'.recal.bam'+'\n')
   # consider removing original bam file
   #qf.write(samtoolspath+'samtools index '+bamstub+'.recal.bam'+'\n') # re-indexing needed?
   qf.write("cp "+bamstub+'.recal.bai '+bamstub+'.recal.bam.bai'+'\n')
   return bamstub+'.recal.bam'

def variant_calling_GATK(GATKpath,dbSNPfile,bam,numthreads):
   # GATK variant calling, including annotation of dbSNP rs numbers
   qf.write('# Variant calling using GATK UnifiedGenotyper - parameters need tweaking'+'\n')
   bamstub=bam.replace('.bam','')
   # in progress   
   qf.write('java7 -jar '+GATKpath+'GenomeAnalysisTK.jar -nt '+str(numthreads)+' -R '+ref+' -T UnifiedGenotyper -I '+bam+' --dbsnp '+dbSNPfile+' --genotype_likelihoods_model BOTH -o '+bamstub+'.UG.raw.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0  -dcov 200'+'\n')  
   return bamstub+'.UG.raw.vcf'

def create_gVCF(gatk_gvcf_path, gvcftools_path, ref, bam, numthreads):
   # create gVCF file
   qf.write('# Create gVCF file using modified GATK UnifiedGenotyper - parameters need tweaking'+'\n')
   bamstub=bam.replace('.bam','')
   qf.write(gvcftools_path+'bin/getBamAvgChromDepth.pl '+bam+' >'+bamstub+'.avgdepth.txt'+'\n')
   qf.write('java7 -jar '+gatk_gvcf_path+'GenomeAnalysisTK.jar -nt '+str(numthreads)+' -R '+ref+' -T UnifiedGenotyper -I '+bam+' -glm BOTH -l OFF -stand_call_conf 20.0 -stand_emit_conf 10.0  -dcov 200  -out_mode EMIT_ALL_SITES | '+gvcftools_path+'bin/gatk_to_gvcf --chrom-depth-file '+bamstub+'.avgdepth.txt | bgzip -c >'+bamstub+'.gvcf.gz'+'\n')
   return bamstub+'.gvcf.gz'

def variant_calling_mpileup(bam,ref,samtoolspath,maxfilterdepth,minfilterdepth):
   # based on Qingyuan's pipeline
   qf.write("# variant calling using the mpileup function of samtools"+'\n')
   bamstub=bam.replace('.bam','')
   bcftoolspath=samtoolspath+'bcftools/'
   qf.write(samtoolspath+"samtools mpileup -C50 -ugf "+ref+' '+bam+' | '+bcftoolspath+'bcftools view -bvcg -| '+bcftoolspath+'bcftools view - | perl '+bcftoolspath+'vcfutils.pl varFilter -D '+maxfilterdepth+' -d '+minfilterdepth+' >'+bamstub+'.var.mpileup.flt.vcf'+'\n')
   return bamstub+'.var.mpileup.flt.vcf'

def variant_calling_pileup(samtoolspath_v12, ref,bam):
   # based on HJM's old Mosaik/BWA pipeline
   # currently quits without error message after having completed chrom 6 - no idea why...
   bamstub=bam.replace('.bam','')
   qf.write('# old-school variant calling using the pileup algorithm'+'\n')
   qf.write("echo 'old-school variant calling using the pileup algorithm'"+'\n')
   qf.write(samtoolspath_v12+'samtools view -u '+bam+' | '+samtoolspath_v12+'samtools pileup -vcf '+ref+' - >'+bamstub+'_vars-raw.txt'+'\n')
   qf.write(r'VAR=`cat '+bamstub+'_vars-raw.txt'+r" | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`"+'\n')
   qf.write('let VAR=2*VAR'+'\n')
   qf.write('echo "max depth is $VAR"'+'\n')
   qf.write(samtoolspath_v12+'misc/samtools.pl varFilter -D$VAR '+bamstub+r"_vars-raw.txt | awk '($3=="+'"*"&&$6>=50)||($3!="*"&&$6>=20)'+r"' >"+bamstub+'.vars-flt_final.txt'+'\n')
   return bamstub+'.vars-flt_final.txt'

def course_nucdiv(abgsamapping_toolpath,samtoolspath_v12,varfile,bam,ref):
   # doing some course nucleotide diversity calculations
   qf.write('# course nucleotide diversity stat generator'+'\n')
   bamstub=bam.replace('.bam','')
   qf.write('VAR=`cat '+varfile+r" | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`"+'\n')
   qf.write('let VAR=2*VAR'+'\n')
   qf.write('echo "max depth is $VAR"'+'\n')
   qf.write(samtoolspath_v12+'samtools view -u '+bam+' | '+samtoolspath_v12+'samtools pileup -f '+ref+r" -c - | awk '$8>4' | awk -v VAR=$VAR '$8<VAR' | perl "+abgsamapping_toolpath+"extract_stats-pileup-bins_allchroms.pl -f "+bamstub+'\n')

def create_shell_script(sample,abgsa,ref,mapper,numthreads,md5check):
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
   dbSNPfile='/media/InternBkp1/repos/dbSNP/Ssc_dbSNP138.vcf'
   gatk_gvcf_path='/opt/GATK/GATK_gVCFmod/'
   gvcftools_path='/opt/gvcftools/v0.13-2-gd92e721/'
   abgsamapping_toolpath='/opt/abgsascripts/'
   maxfilterdepth=20
   minfilterdepth=4

   qf.write('mkdir '+tempdir+'\n')
   tempdir=tempdir+'/'
    
   # get sequence info from database
   archives=get_info_from_db(sample)

   count=0
   bams=[]
   varfiles=[]
   logfile=tempdir+sample+'.log' 
   firstline=r'DATE=`date`; echo "++++++++++++++++++++++++++++" >>'+logfile+'; '
   # preparing fq, trimming, mapping, in a loop,
   # per two gzipped fq files
   # do some logging
   qf.write(firstline+"echo 'starting time: '$DATE  >>"+logfile+'\n')

   for archive in archives:
      seqfiles={}
      count+=1
      qf.write("# archive number "+str(count)+": "+archive[0]+'\n')
      archive_dir=archive[0]
      seqfiles[1]=archive[1]
      offset=check_illumina_or_sanger(abgsa+archive_dir+'/'+seqfiles[1])
      do_md5check(md5check,archive[2], abgsa+archive_dir+'/'+seqfiles[1])
      seqfiles[1]=prepare_temp_fq_files(abgsamapping_toolpath,abgsa,archive_dir,seqfiles[1],tempdir)
      archive = next(archives)
      seqfiles[2]=archive[1]
      do_md5check(md5check,archive[2], abgsa+archive_dir+'/'+seqfiles[2])
      seqfiles[2]=prepare_temp_fq_files(abgsamapping_toolpath,abgsa,archive_dir,seqfiles[2],tempdir)
      seqfiles=trim(abgsamapping_toolpath,tempdir,seqfiles,offset)
      # report when mapping starts
      qf.write(firstline+"echo 'starting "+mapper+" mapping of "+sample+" archive "+str(count)+": '$DATE  >>"+logfile+'\n')
      if mapper == 'bwa-mem':
         bams.append(map_bwa_mem(bwapath, samtoolspath,archive_dir,str(count),ref,tempdir,seqfiles,sample,numthreads))
      elif mapper == 'bwa-aln':
         bams.append(map_bwa_aln(bwapath, samtoolspath,archive_dir,str(count),ref,tempdir,seqfiles,sample,numthreads))
      elif mapper == 'mosaik':  
         bams.append(map_Mosaik(mosaikref, mosaikjump, archive_dir,str(count),tempdir,seqfiles,sample,numthreads))

      # report back when finished mapping, and how large file size is
      qf.write(firstline+"echo 'finished, produced BAM file "+bams[count-1]+" archive "+str(count)+": '$DATE  >>"+logfile+'\n')
      qf.write(r'FSIZE=`stat --printf="%s" '+bams[count-1]+'`; echo "size of file '+bams[count-1]+' is "$FSIZE  >>'+logfile+'\n')

   qf.write("#number of bams: "+str(len(bams))+'\n') 

   # merge and reheader bam files  
   bam=merge_bams(samtoolspath, bams,sample,tempdir)
   print(bam)
   # report back when finished merging, and how large file size is
   qf.write(firstline+"echo 'finished merging, produced BAM file "+bam+": '$DATE  >>"+logfile+'\n')
   qf.write(r'FSIZE=`stat --printf="%s" '+bam+'`; echo "size of file '+bam+' is "$FSIZE  >>'+logfile+'\n')

   # further optimization of bam files
   if dedup == 'samtools':
      bam=dedup_samtools(samtoolspath,bam)
   elif dedup == 'picard':
      bam=dedup_picard(samtoolspath,picardpath,bam)
   print(bam)
   # report back when finished de-dupping, and how large file size is
   qf.write(firstline+"echo 'finished dedupping using "+dedup+", produced BAM file "+bam+": '$DATE  >>"+logfile+'\n')
   qf.write(r'FSIZE=`stat --printf="%s" '+bam+'`; echo "size of file '+bam+' is "$FSIZE  >>'+logfile+'\n')

   bam=re_align(samtoolspath,bam,ref,GATKpath)
   print(bam)

   # report back when finished re-aligning, and how large file size is
   qf.write(firstline+"echo 'finished re-aligning, produced BAM file "+bam+": '$DATE  >>"+logfile+'\n')
   qf.write(r'FSIZE=`stat --printf="%s" '+bam+'`; echo "size of file '+bam+' is "$FSIZE  >>'+logfile+'\n')

   bam=recalibrate(GATKpath,samtoolspath,dbSNPfile,ref,bam)
   print(bam)

   # report back when finished recalibrating, and how large file size is
   qf.write(firstline+"echo 'finished re-aligning, produced BAM file "+bam+": '$DATE  >>"+logfile+'\n')
   qf.write(r'FSIZE=`stat --printf="%s" '+bam+'`; echo "size of file '+bam+' is "$FSIZE  >>'+logfile+'\n')

   # variant calling
   varfiles.append(variant_calling_pileup(samtoolspath_v12,ref,bam))
   varfiles.append(variant_calling_mpileup(bam,ref,samtoolspath,str(maxfilterdepth),str(minfilterdepth)))
   varfiles.append(variant_calling_GATK(GATKpath,dbSNPfile,bam,numthreads))
   varfiles.append(create_gVCF(gatk_gvcf_path, gvcftools_path, ref, bam, numthreads))

   # report back when finished variant calling
   qf.write(firstline+"echo 'finished variant calling: '$DATE  >>"+logfile+'\n')
   print(varfiles)
   
   course_nucdiv(abgsamapping_toolpath,samtoolspath_v12,varfiles[0],bam,ref)

if __name__=="__main__":
   # initialize db cursor
   db = mysql.connector.Connect(user='anonymous',host='localhost',database='ABGSAschema', password='anonymous')
   cursor = db.cursor()
   
   # get command line options
   args = parser.parse_args()
   individual=args.individual_name[0]
   abgsa = args.path_to_abgsa[0]
   mapper=args.mapper
   numthreads=args.number_of_threads
   ref = args.path_to_reference_fasta[0]
   dedup=args.dedup_method
   md5check=args.domd5check
   print('md5check: '+str(md5check))
   print('mapper: '+mapper)
   print('dedupper: '+dedup) 

   # open qsub-file (qf)
   qf=open('run'+individual+'.sh','w')
   # invoke master subroutine
   create_shell_script(individual,abgsa,ref,mapper,numthreads,md5check)
   qf.close()

   cursor.close()
   db.close()

   # optional submitting job
   #os.sys('qsub -q all.q run'+individual+'.sh') 
