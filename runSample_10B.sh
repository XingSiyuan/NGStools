#!/bin/bash
#SBATCH --time=10000
#SBATCH --mem=16000
#SBATCH --ntasks=12
#SBATCH --output=outputSample_10B_%j.txt
#SBATCH --error=error_outputSample_10B_%j.txt
#SBATCH --job-name=Sample_10B
#SBATCH --partition=ABGC_Research
mkdir tmpSample_10B
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting time: '$DATE  >>tmpSample_10B/Sample_10B.log
# archive number 1: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_001_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_001_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_001_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_001_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_001_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_001_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 1: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_001_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_001_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_1\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_001_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_001_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_001_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-1PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-1PE2.sorted.bam archive 1: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-1PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-1PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 2: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_002_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_002_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_002_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_002_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_002_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_002_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 2: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_002_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_002_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_2\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_002_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_002_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_002_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-2PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-2PE2.sorted.bam archive 2: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-2PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-2PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 3: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_003_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_003_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_003_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_003_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_003_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_003_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 3: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_003_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_003_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_3\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_003_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_003_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_003_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-3PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-3PE2.sorted.bam archive 3: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-3PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-3PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 4: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_004_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_004_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_004_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_004_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_004_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_004_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 4: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_004_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_004_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_4\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_004_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_004_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_004_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-4PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-4PE2.sorted.bam archive 4: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-4PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-4PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 5: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_005_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_005_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_005_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_005_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_005_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_005_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 5: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_005_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_005_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_5\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_005_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_005_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_005_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-5PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-5PE2.sorted.bam archive 5: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-5PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-5PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 6: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_006_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_006_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_006_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_006_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_006_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_006_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 6: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_006_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_006_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_6\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_006_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_006_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_006_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-6PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-6PE2.sorted.bam archive 6: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-6PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-6PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 7: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_007_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_007_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_007_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_007_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_007_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_007_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 7: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_007_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_007_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_7\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_007_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_007_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_007_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-7PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-7PE2.sorted.bam archive 7: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-7PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-7PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 8: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_008_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_008_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_008_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_008_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_008_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_008_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 8: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_008_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_008_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_8\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_008_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_008_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_008_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-8PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-8PE2.sorted.bam archive 8: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-8PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-8PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 9: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_009_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_009_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_009_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_009_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_009_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_009_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 9: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_009_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_009_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_9\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_009_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_009_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_009_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-9PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-9PE2.sorted.bam archive 9: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-9PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-9PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 10: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_010_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_010_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_010_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_010_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_010_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_010_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 10: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_010_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_010_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_10\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_010_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_010_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_010_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-10PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-10PE2.sorted.bam archive 10: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-10PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-10PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 11: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_011_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_011_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_011_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_011_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_011_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_011_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 11: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_011_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_011_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_11\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_011_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_011_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_011_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-11PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-11PE2.sorted.bam archive 11: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-11PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-11PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 12: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_012_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_012_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_012_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_012_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_012_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_012_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 12: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_012_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_012_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_12\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_012_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_012_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_012_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-12PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-12PE2.sorted.bam archive 12: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-12PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-12PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 13: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_013_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_013_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_013_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_013_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_013_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_013_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 13: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_013_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_013_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_13\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_013_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_013_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_013_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-13PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-13PE2.sorted.bam archive 13: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-13PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-13PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 14: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_014_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_014_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_014_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_014_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_014_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_014_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 14: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_014_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_014_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_14\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_014_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_014_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_014_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-14PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-14PE2.sorted.bam archive 14: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-14PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-14PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 15: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_015_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_015_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_015_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_015_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_015_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_015_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 15: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_015_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_015_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_15\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_015_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_015_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_015_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-15PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-15PE2.sorted.bam archive 15: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-15PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-15PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 16: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_016_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_016_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_016_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_016_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_016_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_016_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 16: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_016_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_016_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_16\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_016_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_016_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_016_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-16PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-16PE2.sorted.bam archive 16: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-16PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-16PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
# archive number 17: Sample_10B
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_017_R1.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Turkey/DataBeltsville/Sample_10B/10B_ACAGTG_L004_017_R2.fastq.gz | pigz >tmpSample_10B/10B_ACAGTG_L004_017_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.gz -r tmpSample_10B/10B_ACAGTG_L004_017_R2.fastq.gz -o tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.tr -p tmpSample_10B/10B_ACAGTG_L004_017_R2.fastq.tr -s tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.tr
pigz tmpSample_10B/10B_ACAGTG_L004_017_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'starting bwa-aln mapping of Sample_10B archive 17: '$DATE  >>tmpSample_10B/Sample_10B.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa aln -t 12 /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/10B_ACAGTG_L004_017_R2.fastq.tr.gz  >tmpSample_10B/10B_ACAGTG_L004_017_R2.fastq.tr.sai
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa sampe -P /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -r '@RG\tID:Sample_10B_17\tSM:Sample_10B\tPL:ILLUMINA' tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_017_R2.fastq.tr.sai tmpSample_10B/10B_ACAGTG_L004_017_R1.fastq.tr.gz tmpSample_10B/10B_ACAGTG_L004_017_R2.fastq.tr.gz | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Suh - | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpSample_10B/aln-Sample_10B-Sample_10B-17PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished, produced BAM file tmpSample_10B/aln-Sample_10B-Sample_10B-17PE2.sorted.bam archive 17: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/aln-Sample_10B-Sample_10B-17PE2.sorted.bam`; echo "size of file tmpSample_10B/aln-Sample_10B-Sample_10B-17PE2.sorted.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
#number of bams: 17
#multiple bam files --> do merge
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-1PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' >tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-2PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-3PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-4PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-5PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-6PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-7PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-8PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-9PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-10PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-11PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-12PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-13PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-14PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-15PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-16PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpSample_10B/aln-Sample_10B-Sample_10B-17PE2.sorted.bam | sed 's/SM:unknown/SM:Sample_10B/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpSample_10B/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools merge tmpSample_10B/tmpmergedSample_10B.bam tmpSample_10B/aln-Sample_10B-Sample_10B-1PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-2PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-3PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-4PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-5PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-6PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-7PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-8PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-9PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-10PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-11PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-12PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-13PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-14PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-15PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-16PE2.sorted.bam tmpSample_10B/aln-Sample_10B-Sample_10B-17PE2.sorted.bam 
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools reheader tmpSample_10B/newheader.txt tmpSample_10B/tmpmergedSample_10B.bam >tmpSample_10B/Sample_10B_rh.bam
rm tmpSample_10B/tmpmergedSample_10B.bam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools index tmpSample_10B/Sample_10B_rh.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'Produced BAM file tmpSample_10B/Sample_10B_rh.bam: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/Sample_10B_rh.bam`; echo "size of file tmpSample_10B/Sample_10B_rh.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
MD5BAM=`md5sum tmpSample_10B/Sample_10B_rh.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpSample_10B/Sample_10B_rh.bam is "$MD5BAM  >>tmpSample_10B/Sample_10B.log
# dedup using samtools
echo 'dedupping using samtools'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools rmdup tmpSample_10B/Sample_10B_rh.bam tmpSample_10B/Sample_10B_rh.dedup_st.bam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools index tmpSample_10B/Sample_10B_rh.dedup_st.bam
cp tmpSample_10B/Sample_10B_rh.dedup_st.bam.bai tmpSample_10B/Sample_10B_rh.dedup_st.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'Produced BAM file tmpSample_10B/Sample_10B_rh.dedup_st.bam: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/Sample_10B_rh.dedup_st.bam`; echo "size of file tmpSample_10B/Sample_10B_rh.dedup_st.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
MD5BAM=`md5sum tmpSample_10B/Sample_10B_rh.dedup_st.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpSample_10B/Sample_10B_rh.dedup_st.bam is "$MD5BAM  >>tmpSample_10B/Sample_10B.log
# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner
java7 -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 12 -T RealignerTargetCreator -R /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -I tmpSample_10B/Sample_10B_rh.dedup_st.bam -o tmpSample_10B/Sample_10B_rh.dedup_st.reA.intervals
java7 -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -I tmpSample_10B/Sample_10B_rh.dedup_st.bam -targetIntervals tmpSample_10B/Sample_10B_rh.dedup_st.reA.intervals -o tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam
cp tmpSample_10B/Sample_10B_rh.dedup_st.reA.bai tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'Produced BAM file tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam`; echo "size of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam is "$FSIZE  >>tmpSample_10B/Sample_10B.log
MD5BAM=`md5sum tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam is "$MD5BAM  >>tmpSample_10B/Sample_10B.log
# Calculate coverage statistics
java7 -Xmx8G -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -T DepthOfCoverage -R /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -I tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam.coverage
# old-school variant calling using the pileup algorithm
echo 'old-school variant calling using the pileup algorithm'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools view -u tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools pileup -vcf /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa - >tmpSample_10B/Sample_10B_rh.dedup_st.reA_vars-raw.txt
VAR=`cat tmpSample_10B/Sample_10B_rh.dedup_st.reA_vars-raw.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/misc/samtools.pl varFilter -D$VAR tmpSample_10B/Sample_10B_rh.dedup_st.reA_vars-raw.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >tmpSample_10B/Sample_10B_rh.dedup_st.reA.vars-flt_final.txt
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'Produced variant file tmpSample_10B/Sample_10B_rh.dedup_st.reA.vars-flt_final.txt: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/Sample_10B_rh.dedup_st.reA.vars-flt_final.txt`; echo "size of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.vars-flt_final.txt is "$FSIZE  >>tmpSample_10B/Sample_10B.log
MD5VAR=`md5sum tmpSample_10B/Sample_10B_rh.dedup_st.reA.vars-flt_final.txt | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.vars-flt_final.txt is "$MD5VAR  >>tmpSample_10B/Sample_10B.log
# variant calling using the mpileup function of samtools
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools mpileup -C50 -ugf /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/bcftools view -bvcg -| /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/bcftools view - | perl /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/vcfutils.pl varFilter -D 20 -d 4 >tmpSample_10B/Sample_10B_rh.dedup_st.reA.var.mpileup.flt.vcf
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'Produced variant file tmpSample_10B/Sample_10B_rh.dedup_st.reA.var.mpileup.flt.vcf: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/Sample_10B_rh.dedup_st.reA.var.mpileup.flt.vcf`; echo "size of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.var.mpileup.flt.vcf is "$FSIZE  >>tmpSample_10B/Sample_10B.log
MD5VAR=`md5sum tmpSample_10B/Sample_10B_rh.dedup_st.reA.var.mpileup.flt.vcf | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.var.mpileup.flt.vcf is "$MD5VAR  >>tmpSample_10B/Sample_10B.log
# Variant calling using GATK UnifiedGenotyper - parameters need tweaking
java7 -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 12 -R /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -T UnifiedGenotyper -I tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam --genotype_likelihoods_model BOTH -o tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0  -dcov 200
bgzip tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf
tabix -p vcf tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'Produced variant file tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf.gz: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf.gz`; echo "size of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf.gz is "$FSIZE  >>tmpSample_10B/Sample_10B.log
MD5VAR=`md5sum tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf.gz | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf.gz is "$MD5VAR  >>tmpSample_10B/Sample_10B.log
# Create gVCF file using modified GATK UnifiedGenotyper - parameters need tweaking
/cm/shared/apps/WUR/ABGC/gvcftools/gvcftools-0.16/bin/getBamAvgChromDepth.pl tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam >tmpSample_10B/Sample_10B_rh.dedup_st.reA.avgdepth.txt
java7 -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK_gVCFmod/GenomeAnalysisTK.jar -nt 12 -R /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -T UnifiedGenotyper -I tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam -glm BOTH -l OFF -stand_call_conf 20.0 -stand_emit_conf 10.0  -dcov 200  -out_mode EMIT_ALL_SITES | /cm/shared/apps/WUR/ABGC/gvcftools/gvcftools-0.16/bin/gatk_to_gvcf --chrom-depth-file tmpSample_10B/Sample_10B_rh.dedup_st.reA.avgdepth.txt | bgzip -c >tmpSample_10B/Sample_10B_rh.dedup_st.reA.gvcf.gz
tabix -p vcf tmpSample_10B/Sample_10B_rh.dedup_st.reA.gvcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'Produced variant file tmpSample_10B/Sample_10B_rh.dedup_st.reA.gvcf.gz: '$DATE  >>tmpSample_10B/Sample_10B.log
FSIZE=`stat --printf="%s" tmpSample_10B/Sample_10B_rh.dedup_st.reA.gvcf.gz`; echo "size of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.gvcf.gz is "$FSIZE  >>tmpSample_10B/Sample_10B.log
MD5VAR=`md5sum tmpSample_10B/Sample_10B_rh.dedup_st.reA.gvcf.gz | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpSample_10B/Sample_10B_rh.dedup_st.reA.gvcf.gz is "$MD5VAR  >>tmpSample_10B/Sample_10B.log
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpSample_10B/Sample_10B.log; echo 'finished variant calling: '$DATE  >>tmpSample_10B/Sample_10B.log
# predicting function using VEP
gunzip -c tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vcf.gz | perl /cm/shared/apps/WUR/ABGC/variant_effect_predictor/VEP231213/variant_effect_predictor.pl --dir /cm/shared/apps/WUR/ABGC/variant_effect_predictor/VEP231213//cache --species meleagris_gallopavo -o tmpSample_10B/Sample_10B_rh.dedup_st.reA.UG.raw.vep.txt --fork 12 --canonical --coding_only --no_intergenic --offline --force_overwrite
# course nucleotide diversity stat generator
VAR=`cat tmpSample_10B/Sample_10B_rh.dedup_st.reA.vars-flt_final.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools view -u tmpSample_10B/Sample_10B_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools pileup -f /lustre/nobackup/WUR/ABGC/shared/Turkey/UMD5/UMD_5_genome_nospaces_60.fa -c - | awk '$8>4' | awk -v VAR=$VAR '$8<VAR' | perl /cm/shared/apps/WUR/ABGC/abgsascripts/extract_stats-pileup-bins_allchroms.pl -f tmpSample_10B/Sample_10B_rh.dedup_st.reA
