#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h_vmem=20G
mkdir tmpBOV-WUR-4

DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'starting time: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
# archive number 1: SZAIPI019133-23
python /opt/abgsascripts/fix_fq_names.py /srv/mds01/shared/Bulls1000/F12FPCEUHK0755_alq121122/cleandata/SZAIPI019133-23/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /srv/mds01/shared/Bulls1000/F12FPCEUHK0755_alq121122/cleandata/SZAIPI019133-23/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.gz
## quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.gz -r tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.gz -o tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr -p tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr -s tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.singles.tr -l 50 -t illumina
pigz tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr
pigz tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.gz | gzip -c >tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.sa.gz
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.gz | gzip -c >tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.sa.gz
rm tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.gz
rm tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.gz
mv tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.sa.gz tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.gz
mv tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.sa.gz tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'starting bwa-aln mapping of BOV-WUR-4 archive 1: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.5.9/bwa aln -t 12 /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.gz  >tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa aln -t 12 /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.gz  >tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa sampe -P /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa -r '@RG\tID:SZAIPI019133-23_1\tSM:HOLNLDM000292559732\tPL:ILLUMINA' tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.sai tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.sai tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_1.fq.tr.gz tmpBOV-WUR-4/121205_I288_FCD1JLDACXX_L3_SZAIPI019133-23_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-1PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'finished, produced BAM file tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-1PE2.sorted.bam archive 1: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-1PE2.sorted.bam`; echo "size of file tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-1PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-4/BOV-WUR-4.log
# archive number 2: SZAIPI019133-23
python /opt/abgsascripts/fix_fq_names.py /srv/mds01/shared/Bulls1000/F12FPCEUHK0755_alq121122/cleandata/SZAIPI019133-23/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /srv/mds01/shared/Bulls1000/F12FPCEUHK0755_alq121122/cleandata/SZAIPI019133-23/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.gz -r tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.gz -o tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr -p tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr -s tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.singles.tr -l 50 -t illumina
pigz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr
pigz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.gz | gzip -c >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.sa.gz
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.gz | gzip -c >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.sa.gz
rm tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.gz
rm tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.gz
mv tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.sa.gz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.gz
mv tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.sa.gz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'starting bwa-aln mapping of BOV-WUR-4 archive 2: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.5.9/bwa aln -t 12 /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.gz  >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa aln -t 12 /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.gz  >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa sampe -P /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa -r '@RG\tID:SZAIPI019133-23_2\tSM:HOLNLDM000292559732\tPL:ILLUMINA' tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.sai tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.sai tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_1.fq.tr.gz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L1_SZAIPI019133-23_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-2PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'finished, produced BAM file tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-2PE2.sorted.bam archive 2: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-2PE2.sorted.bam`; echo "size of file tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-2PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-4/BOV-WUR-4.log
# archive number 3: SZAIPI019133-23
python /opt/abgsascripts/fix_fq_names.py /srv/mds01/shared/Bulls1000/F12FPCEUHK0755_alq121122/cleandata/SZAIPI019133-23/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /srv/mds01/shared/Bulls1000/F12FPCEUHK0755_alq121122/cleandata/SZAIPI019133-23/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.gz -r tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.gz -o tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr -p tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr -s tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.singles.tr -l 50 -t illumina
pigz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr
pigz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.gz | gzip -c >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.sa.gz
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.gz | gzip -c >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.sa.gz
rm tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.gz
rm tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.gz
mv tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.sa.gz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.gz
mv tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.sa.gz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'starting bwa-aln mapping of BOV-WUR-4 archive 3: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.5.9/bwa aln -t 12 /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.gz  >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa aln -t 12 /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.gz  >tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa sampe -P /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa -r '@RG\tID:SZAIPI019133-23_3\tSM:HOLNLDM000292559732\tPL:ILLUMINA' tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.sai tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.sai tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_1.fq.tr.gz tmpBOV-WUR-4/121205_I812_FCD1JHHACXX_L2_SZAIPI019133-23_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-3PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'finished, produced BAM file tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-3PE2.sorted.bam archive 3: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-3PE2.sorted.bam`; echo "size of file tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-3PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-4/BOV-WUR-4.log
#number of bams: 3
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-1PE2.sorted.bam | sed 's/SM:unknown/SM:HOLNLDM000292559732/'  | sed 's/PL:sanger/PL:ILLUMINA/' >tmpBOV-WUR-4/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-2PE2.sorted.bam | sed 's/SM:unknown/SM:HOLNLDM000292559732/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBOV-WUR-4/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-3PE2.sorted.bam | sed 's/SM:unknown/SM:HOLNLDM000292559732/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBOV-WUR-4/newheader.txt
/opt/samtools/samtools-0.1.19/samtools merge tmpBOV-WUR-4/tmpmergedBOV-WUR-4.bam tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-1PE2.sorted.bam tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-2PE2.sorted.bam tmpBOV-WUR-4/aln-SZAIPI019133-23-BOV-WUR-4-3PE2.sorted.bam 
/opt/samtools/samtools-0.1.19/samtools reheader tmpBOV-WUR-4/newheader.txt tmpBOV-WUR-4/tmpmergedBOV-WUR-4.bam >tmpBOV-WUR-4/BOV-WUR-4_rh.bam
rm tmpBOV-WUR-4/tmpmergedBOV-WUR-4.bam
/opt/samtools/samtools-0.1.19/samtools index tmpBOV-WUR-4/BOV-WUR-4_rh.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'finished merging, produced BAM file tmpBOV-WUR-4/BOV-WUR-4_rh.bam: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-4/BOV-WUR-4_rh.bam`; echo "size of file tmpBOV-WUR-4/BOV-WUR-4_rh.bam is "$FSIZE  >>tmpBOV-WUR-4/BOV-WUR-4.log
# dedup using samtools
echo 'dedupping using samtools'
/opt/samtools/samtools-0.1.19/samtools rmdup tmpBOV-WUR-4/BOV-WUR-4_rh.bam tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bam
/opt/samtools/samtools-0.1.19/samtools index tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bam
cp tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bam.bai tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'finished dedupping using samtools, produced BAM file tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bam: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bam`; echo "size of file tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bam is "$FSIZE  >>tmpBOV-WUR-4/BOV-WUR-4.log
# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 12 -T RealignerTargetCreator -R /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa -I tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bam -o tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.intervals
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /srv/mds01/shared/Bulls1000/UMD31/umd_3_1.fa -I tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.bam -targetIntervals tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.intervals -o tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bam
cp tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bai tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'finished re-aligning, produced BAM file tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bam: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bam`; echo "size of file tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bam is "$FSIZE  >>tmpBOV-WUR-4/BOV-WUR-4.log
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-4/BOV-WUR-4.log; echo 'finished re-aligning, produced BAM file tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bam: '$DATE  >>tmpBOV-WUR-4/BOV-WUR-4.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bam`; echo "size of file tmpBOV-WUR-4/BOV-WUR-4_rh.dedup_st.reA.bam is "$FSIZE  >>tmpBOV-WUR-4/BOV-WUR-4.log
