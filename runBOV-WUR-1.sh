#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h_vmem=20G
mkdir tmpBOV-WUR-1
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'starting time: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# archive number 1: SZAIPI019130-16
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.gz -r tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.gz -o tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr -p tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr -s tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.singles.tr -l 50 -t illumina
pigz tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr
pigz tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.gz | gzip -c >tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.sa.gz
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.gz | gzip -c >tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.sa.gz
rm tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.gz
rm tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.gz
mv tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.sa.gz tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.gz
mv tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.sa.gz tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'starting bwa-aln mapping of BOV-WUR-1 archive 1: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.gz  >tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.gz  >tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa sampe -P /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -r '@RG\tID:SZAIPI019130-16_1\tSM:HOLNLDM000120873995\tPL:ILLUMINA' tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.sai tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.sai tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.gz tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished, produced BAM file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam archive 1: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam`; echo "size of file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# archive number 2: SZAIPI019130-16
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.gz -r tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.gz -o tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr -p tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr -s tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.singles.tr -l 50 -t illumina
pigz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr
pigz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.gz | gzip -c >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.sa.gz
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.gz | gzip -c >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.sa.gz
rm tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.gz
rm tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.gz
mv tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.sa.gz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.gz
mv tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.sa.gz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'starting bwa-aln mapping of BOV-WUR-1 archive 2: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.gz  >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.gz  >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa sampe -P /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -r '@RG\tID:SZAIPI019130-16_2\tSM:HOLNLDM000120873995\tPL:ILLUMINA' tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.sai tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.sai tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.gz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished, produced BAM file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam archive 2: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam`; echo "size of file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# archive number 3: SZAIPI019130-16
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.gz -r tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.gz -o tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr -p tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr -s tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.singles.tr -l 50 -t illumina
pigz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr
pigz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.gz | gzip -c >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.sa.gz
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.gz | gzip -c >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.sa.gz
rm tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.gz
rm tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.gz
mv tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.sa.gz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.gz
mv tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.sa.gz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'starting bwa-aln mapping of BOV-WUR-1 archive 3: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.gz  >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.gz  >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa sampe -P /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -r '@RG\tID:SZAIPI019130-16_3\tSM:HOLNLDM000120873995\tPL:ILLUMINA' tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.sai tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.sai tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.gz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished, produced BAM file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam archive 3: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam`; echo "size of file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
#number of bams: 3
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam | sed 's/SM:unknown/SM:HOLNLDM000120873995/'  | sed 's/PL:sanger/PL:ILLUMINA/' >tmpBOV-WUR-1/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam | sed 's/SM:unknown/SM:HOLNLDM000120873995/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBOV-WUR-1/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam | sed 's/SM:unknown/SM:HOLNLDM000120873995/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBOV-WUR-1/newheader.txt
/opt/samtools/samtools-0.1.19/samtools merge tmpBOV-WUR-1/tmpmergedBOV-WUR-1.bam tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam 
/opt/samtools/samtools-0.1.19/samtools reheader tmpBOV-WUR-1/newheader.txt tmpBOV-WUR-1/tmpmergedBOV-WUR-1.bam >tmpBOV-WUR-1/BOV-WUR-1_rh.bam
rm tmpBOV-WUR-1/tmpmergedBOV-WUR-1.bam
/opt/samtools/samtools-0.1.19/samtools index tmpBOV-WUR-1/BOV-WUR-1_rh.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished merging, produced BAM file tmpBOV-WUR-1/BOV-WUR-1_rh.bam: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/BOV-WUR-1_rh.bam`; echo "size of file tmpBOV-WUR-1/BOV-WUR-1_rh.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# dedup using samtools
echo 'dedupping using samtools'
/opt/samtools/samtools-0.1.19/samtools rmdup tmpBOV-WUR-1/BOV-WUR-1_rh.bam tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam
/opt/samtools/samtools-0.1.19/samtools index tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam
cp tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam.bai tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished dedupping using samtools, produced BAM file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam`; echo "size of file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -I tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam -o tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.intervals
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -I tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam -targetIntervals tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.intervals -o tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam
cp tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bai tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished re-aligning, produced BAM file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam`; echo "size of file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished re-aligning, produced BAM file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam`; echo "size of file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
