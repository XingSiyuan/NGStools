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
sickle pe -f tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.gz -r tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.gz -o tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr -p tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr -s tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.singles.tr -l 45 -t illumina
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
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.gz  >tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.gz  >tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -r '@RG\tID:SZAIPI019130-16_1\tSM:BOV-WUR-1\tPL:ILLUMINA' tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.sai tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.sai tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_1.fq.tr.gz tmpBOV-WUR-1/121202_I598_FCD1JRLACXX_L6_SZAIPI019130-16_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished, produced BAM file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam archive 1: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam`; echo "size of file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# archive number 2: SZAIPI019130-16
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.gz -r tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.gz -o tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr -p tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr -s tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.singles.tr -l 45 -t illumina
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
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.gz  >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.gz  >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -r '@RG\tID:SZAIPI019130-16_2\tSM:BOV-WUR-1\tPL:ILLUMINA' tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.sai tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.sai tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_1.fq.tr.gz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L1_SZAIPI019130-16_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished, produced BAM file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam archive 2: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam`; echo "size of file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# archive number 3: SZAIPI019130-16
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/cowrepo_test/SZAIPI019130-16/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.gz.clean.dup.clean.gz | pigz >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.gz -r tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.gz -o tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr -p tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr -s tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.singles.tr -l 45 -t illumina
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
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.gz  >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.gz  >tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -r '@RG\tID:SZAIPI019130-16_3\tSM:BOV-WUR-1\tPL:ILLUMINA' tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.sai tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.sai tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_1.fq.tr.gz tmpBOV-WUR-1/121205_I812_FCD1JHHACXX_L2_SZAIPI019130-16_2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished, produced BAM file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam archive 3: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam`; echo "size of file tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
#number of bams: 3
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-1PE2.sorted.bam | sed 's/SM:unknown/SM:BOV-WUR-1/'  | sed 's/PL:sanger/PL:ILLUMINA/' >tmpBOV-WUR-1/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-2PE2.sorted.bam | sed 's/SM:unknown/SM:BOV-WUR-1/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBOV-WUR-1/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBOV-WUR-1/aln-SZAIPI019130-16-BOV-WUR-1-3PE2.sorted.bam | sed 's/SM:unknown/SM:BOV-WUR-1/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBOV-WUR-1/newheader.txt
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
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -I tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam -o tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.intervals
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -I tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.bam -targetIntervals tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.intervals -o tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam
cp tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bai tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished re-aligning, produced BAM file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam`; echo "size of file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# Recalibration of BAM using GATK-BaseRecalibrator+PrintReads
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nct 4 -T BaseRecalibrator -R /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -I tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam -knownSites /media/InternBkp1/repos/cowrepo_test/refs/dbSNP/dbSNP.vcf -o tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.grp
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nct 4 -T PrintReads -R /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -I tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.bam -BQSR tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.grp -o tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam
cp tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bai tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished re-aligning, produced BAM file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
FSIZE=`stat --printf="%s" tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam`; echo "size of file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam is "$FSIZE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# old-school variant calling using the pileup algorithm
echo 'old-school variant calling using the pileup algorithm'
/opt/samtools/samtools-0.1.12a/samtools view -u tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.12a/samtools pileup -vcf /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa - >tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal_vars-raw.txt
VAR=`cat tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal_vars-raw.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/opt/samtools/samtools-0.1.12a/misc/samtools.pl varFilter -D$VAR tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal_vars-raw.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.vars-flt_final.txt
# variant calling using the mpileup function of samtools
/opt/samtools/samtools-0.1.19/samtools mpileup -C50 -ugf /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.19/bcftools/bcftools view -bvcg -| /opt/samtools/samtools-0.1.19/bcftools/bcftools view - | perl /opt/samtools/samtools-0.1.19/bcftools/vcfutils.pl varFilter -D 20 -d 4 >tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.var.mpileup.flt.vcf
# Variant calling using GATK UnifiedGenotyper - parameters need tweaking
java7 -Xmx8g -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -R /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -T UnifiedGenotyper -I tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam --dbsnp /media/InternBkp1/repos/cowrepo_test/refs/dbSNP/dbSNP.vcf --genotype_likelihoods_model BOTH -o tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.UG.raw.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0  -dcov 200
bgzip tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.UG.raw.vcf
tabix -p vcf tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.UG.raw.vcf.gz
# Create gVCF file using modified GATK UnifiedGenotyper - parameters need tweaking
/opt/gvcftools/v0.13-2-gd92e721/bin/getBamAvgChromDepth.pl tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam >tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.avgdepth.txt
java7 -Xmx8g -jar /opt/GATK/GATK_gVCFmod/GenomeAnalysisTK.jar -nt 4 -R /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -T UnifiedGenotyper -I tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam -glm BOTH -l OFF -stand_call_conf 20.0 -stand_emit_conf 10.0  -dcov 200  -out_mode EMIT_ALL_SITES | /opt/gvcftools/v0.13-2-gd92e721/bin/gatk_to_gvcf --chrom-depth-file tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.avgdepth.txt | bgzip -c >tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.gvcf.gz
tabix -p vcf tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.gvcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBOV-WUR-1/BOV-WUR-1.log; echo 'finished variant calling: '$DATE  >>tmpBOV-WUR-1/BOV-WUR-1.log
# predicting function using VEP
gunzip -c tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.UG.raw.vcf.gz | perl /opt/VEP/variant_effect_predictor.pl --dir /opt/VEP/ --species sus_scrofa -o tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.UG.raw.vep.txt --fork 4 --canonical --sift b --coding_only --no_intergenic --offline --force_overwrite
# course nucleotide diversity stat generator
VAR=`cat tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.vars-flt_final.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/opt/samtools/samtools-0.1.12a/samtools view -u tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.12a/samtools pileup -f /media/InternBkp1/repos/cowrepo_test/refs/umd_3_1_reference_1000_bull_genomes.fa -c - | awk '$8>4' | awk -v VAR=$VAR '$8<VAR' | perl /opt/abgsascripts/extract_stats-pileup-bins_allchroms.pl -f tmpBOV-WUR-1/BOV-WUR-1_rh.dedup_st.reA.recal
