#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h_vmem=20G
mkdir tmpBOV-WUR-2
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'starting time: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
# archive number 1: SZAIPI019131-17
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019131-17/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019131-17/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.gz -r tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.gz -o tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr -p tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr -s tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.singles.tr -l 45 -t illumina
pigz tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr
pigz tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.gz | gzip -c >tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.sa.gz
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.gz | gzip -c >tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.sa.gz
rm tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.gz
rm tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.gz
mv tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.sa.gz tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.gz
mv tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.sa.gz tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'starting bwa-aln mapping of BOV-WUR-2 archive 1: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.gz  >tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.gz  >tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa sampe -P /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -r '@RG\tID:SZAIPI019131-17_1\tSM:BOV-WUR-2\tPL:ILLUMINA' tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.sai tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.sai tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_1.fq.tr.gz tmpBOV-WUR-2/121205_I288_FCD1JLDACXX_L1_SZAIPI019131-17_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-1PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'finished, produced BAM file tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-1PE2.sorted.bam archive 1: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-1PE2.sorted.bam`; echo "size of file tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-1PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-2/BOV-WUR-2.log
# archive number 2: SZAIPI019131-17
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019131-17/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019131-17/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.gz -r tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.gz -o tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr -p tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr -s tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.singles.tr -l 45 -t illumina
pigz tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr
pigz tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.gz | gzip -c >tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.sa.gz
python /opt/abgsascripts/convert_ill_to_sang.py tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.gz | gzip -c >tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.sa.gz
rm tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.gz
rm tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.gz
mv tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.sa.gz tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.gz
mv tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.sa.gz tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'starting bwa-aln mapping of BOV-WUR-2 archive 2: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.gz  >tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa aln -t 4 /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.gz  >tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.sai
/opt/bwa/bwa-0.5.9/bwa sampe -P /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -r '@RG\tID:SZAIPI019131-17_2\tSM:BOV-WUR-2\tPL:ILLUMINA' tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.sai tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.sai tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_1.fq.tr.gz tmpBOV-WUR-2/121205_I812_FCD1JHHACXX_L7_SZAIPI019131-17_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-2PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'finished, produced BAM file tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-2PE2.sorted.bam archive 2: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-2PE2.sorted.bam`; echo "size of file tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-2PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-2/BOV-WUR-2.log
#number of bams: 2
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-1PE2.sorted.bam | sed 's/SM:unknown/SM:BOV-WUR-2/'  | sed 's/PL:sanger/PL:ILLUMINA/' >tmpBOV-WUR-2/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-2PE2.sorted.bam | sed 's/SM:unknown/SM:BOV-WUR-2/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBOV-WUR-2/newheader.txt
/opt/samtools/samtools-0.1.19/samtools merge tmpBOV-WUR-2/tmpmergedBOV-WUR-2.bam tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-1PE2.sorted.bam tmpBOV-WUR-2/aln-SZAIPI019131-17-BOV-WUR-2-2PE2.sorted.bam 
/opt/samtools/samtools-0.1.19/samtools reheader tmpBOV-WUR-2/newheader.txt tmpBOV-WUR-2/tmpmergedBOV-WUR-2.bam >tmpBOV-WUR-2/BOV-WUR-2_rh.bam
rm tmpBOV-WUR-2/tmpmergedBOV-WUR-2.bam
/opt/samtools/samtools-0.1.19/samtools index tmpBOV-WUR-2/BOV-WUR-2_rh.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'finished merging, produced BAM file tmpBOV-WUR-2/BOV-WUR-2_rh.bam: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-2/BOV-WUR-2_rh.bam`; echo "size of file tmpBOV-WUR-2/BOV-WUR-2_rh.bam is "$FSIZE  >>tmpBOV-WUR-2/BOV-WUR-2.log
# dedup using samtools
echo 'dedupping using samtools'
/opt/samtools/samtools-0.1.19/samtools rmdup tmpBOV-WUR-2/BOV-WUR-2_rh.bam tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bam
/opt/samtools/samtools-0.1.19/samtools index tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bam
cp tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bam.bai tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'finished dedupping using samtools, produced BAM file tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bam: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bam`; echo "size of file tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bam is "$FSIZE  >>tmpBOV-WUR-2/BOV-WUR-2.log
# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -I tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bam -o tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.intervals
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /media/InternBkp1/repos/cowrepo_test/UMD31/umd_3_1 -I tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.bam -targetIntervals tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.intervals -o tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bam
cp tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bai tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'finished re-aligning, produced BAM file tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bam: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bam`; echo "size of file tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bam is "$FSIZE  >>tmpBOV-WUR-2/BOV-WUR-2.log
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-2/BOV-WUR-2.log; echo 'finished re-aligning, produced BAM file tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bam: '$DATE  >>tmpBOV-WUR-2/BOV-WUR-2.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bam`; echo "size of file tmpBOV-WUR-2/BOV-WUR-2_rh.dedup_st.reA.bam is "$FSIZE  >>tmpBOV-WUR-2/BOV-WUR-2.log
