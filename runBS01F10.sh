#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h_vmem=10G
mkdir tmpBS01F10
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'starting time: '$DATE  >>tmpBS01F10/BS01F10.log
# archive number 1: ABGSA0215
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0215/ABGSA0215_BS01F10_r131_R1.fq.gz | pigz >tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0215/ABGSA0215_BS01F10_r131_R2.fq.gz | pigz >tmpBS01F10/ABGSA0215_BS01F10_r131_R2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.gz -r tmpBS01F10/ABGSA0215_BS01F10_r131_R2.fq.gz -o tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.tr -p tmpBS01F10/ABGSA0215_BS01F10_r131_R2.fq.tr -s tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.singles.tr -l 45 -t sanger
pigz tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.tr
pigz tmpBS01F10/ABGSA0215_BS01F10_r131_R2.fq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'starting bwa-aln mapping of BS01F10 archive 1: '$DATE  >>tmpBS01F10/BS01F10.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.tr.gz  >tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpBS01F10/ABGSA0215_BS01F10_r131_R2.fq.tr.gz  >tmpBS01F10/ABGSA0215_BS01F10_r131_R2.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -r '@RG\tID:ABGSA0215_1\tSM:BS01F10\tPL:ILLUMINA' tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.tr.sai tmpBS01F10/ABGSA0215_BS01F10_r131_R2.fq.tr.sai tmpBS01F10/ABGSA0215_BS01F10_r131_R1.fq.tr.gz tmpBS01F10/ABGSA0215_BS01F10_r131_R2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBS01F10/aln-ABGSA0215-BS01F10-1PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'finished, produced BAM file tmpBS01F10/aln-ABGSA0215-BS01F10-1PE2.sorted.bam archive 1: '$DATE  >>tmpBS01F10/BS01F10.log
FSIZE=`stat --printf="%s" tmpBS01F10/aln-ABGSA0215-BS01F10-1PE2.sorted.bam`; echo "size of file tmpBS01F10/aln-ABGSA0215-BS01F10-1PE2.sorted.bam is "$FSIZE  >>tmpBS01F10/BS01F10.log
# archive number 2: ABGSA0215
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0215/ABGSA0215_BS01F10_r134_R1.fq.gz | pigz >tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0215/ABGSA0215_BS01F10_r134_R2.fq.gz | pigz >tmpBS01F10/ABGSA0215_BS01F10_r134_R2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.gz -r tmpBS01F10/ABGSA0215_BS01F10_r134_R2.fq.gz -o tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.tr -p tmpBS01F10/ABGSA0215_BS01F10_r134_R2.fq.tr -s tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.singles.tr -l 45 -t sanger
pigz tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.tr
pigz tmpBS01F10/ABGSA0215_BS01F10_r134_R2.fq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'starting bwa-aln mapping of BS01F10 archive 2: '$DATE  >>tmpBS01F10/BS01F10.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.tr.gz  >tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpBS01F10/ABGSA0215_BS01F10_r134_R2.fq.tr.gz  >tmpBS01F10/ABGSA0215_BS01F10_r134_R2.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -r '@RG\tID:ABGSA0215_2\tSM:BS01F10\tPL:ILLUMINA' tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.tr.sai tmpBS01F10/ABGSA0215_BS01F10_r134_R2.fq.tr.sai tmpBS01F10/ABGSA0215_BS01F10_r134_R1.fq.tr.gz tmpBS01F10/ABGSA0215_BS01F10_r134_R2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBS01F10/aln-ABGSA0215-BS01F10-2PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'finished, produced BAM file tmpBS01F10/aln-ABGSA0215-BS01F10-2PE2.sorted.bam archive 2: '$DATE  >>tmpBS01F10/BS01F10.log
FSIZE=`stat --printf="%s" tmpBS01F10/aln-ABGSA0215-BS01F10-2PE2.sorted.bam`; echo "size of file tmpBS01F10/aln-ABGSA0215-BS01F10-2PE2.sorted.bam is "$FSIZE  >>tmpBS01F10/BS01F10.log
# archive number 3: ABGSA0215
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0215/ABGSA0215_BS01F10_r137_R1.fq.gz | pigz >tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0215/ABGSA0215_BS01F10_r137_R2.fq.gz | pigz >tmpBS01F10/ABGSA0215_BS01F10_r137_R2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.gz -r tmpBS01F10/ABGSA0215_BS01F10_r137_R2.fq.gz -o tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.tr -p tmpBS01F10/ABGSA0215_BS01F10_r137_R2.fq.tr -s tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.singles.tr -l 45 -t sanger
pigz tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.tr
pigz tmpBS01F10/ABGSA0215_BS01F10_r137_R2.fq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'starting bwa-aln mapping of BS01F10 archive 3: '$DATE  >>tmpBS01F10/BS01F10.log
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.tr.gz  >tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpBS01F10/ABGSA0215_BS01F10_r137_R2.fq.tr.gz  >tmpBS01F10/ABGSA0215_BS01F10_r137_R2.fq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -r '@RG\tID:ABGSA0215_3\tSM:BS01F10\tPL:ILLUMINA' tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.tr.sai tmpBS01F10/ABGSA0215_BS01F10_r137_R2.fq.tr.sai tmpBS01F10/ABGSA0215_BS01F10_r137_R1.fq.tr.gz tmpBS01F10/ABGSA0215_BS01F10_r137_R2.fq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpBS01F10/aln-ABGSA0215-BS01F10-3PE2.sorted
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'finished, produced BAM file tmpBS01F10/aln-ABGSA0215-BS01F10-3PE2.sorted.bam archive 3: '$DATE  >>tmpBS01F10/BS01F10.log
FSIZE=`stat --printf="%s" tmpBS01F10/aln-ABGSA0215-BS01F10-3PE2.sorted.bam`; echo "size of file tmpBS01F10/aln-ABGSA0215-BS01F10-3PE2.sorted.bam is "$FSIZE  >>tmpBS01F10/BS01F10.log
#number of bams: 3
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpBS01F10/aln-ABGSA0215-BS01F10-1PE2.sorted.bam | sed 's/SM:unknown/SM:BS01F10/'  | sed 's/PL:sanger/PL:ILLUMINA/' >tmpBS01F10/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBS01F10/aln-ABGSA0215-BS01F10-2PE2.sorted.bam | sed 's/SM:unknown/SM:BS01F10/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBS01F10/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpBS01F10/aln-ABGSA0215-BS01F10-3PE2.sorted.bam | sed 's/SM:unknown/SM:BS01F10/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpBS01F10/newheader.txt
/opt/samtools/samtools-0.1.19/samtools merge tmpBS01F10/tmpmergedBS01F10.bam tmpBS01F10/aln-ABGSA0215-BS01F10-1PE2.sorted.bam tmpBS01F10/aln-ABGSA0215-BS01F10-2PE2.sorted.bam tmpBS01F10/aln-ABGSA0215-BS01F10-3PE2.sorted.bam 
/opt/samtools/samtools-0.1.19/samtools reheader tmpBS01F10/newheader.txt tmpBS01F10/tmpmergedBS01F10.bam >tmpBS01F10/BS01F10_rh.bam
rm tmpBS01F10/tmpmergedBS01F10.bam
/opt/samtools/samtools-0.1.19/samtools index tmpBS01F10/BS01F10_rh.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'finished merging, produced BAM file tmpBS01F10/BS01F10_rh.bam: '$DATE  >>tmpBS01F10/BS01F10.log
FSIZE=`stat --printf="%s" tmpBS01F10/BS01F10_rh.bam`; echo "size of file tmpBS01F10/BS01F10_rh.bam is "$FSIZE  >>tmpBS01F10/BS01F10.log
# dedup using samtools
echo 'dedupping using samtools'
/opt/samtools/samtools-0.1.19/samtools rmdup tmpBS01F10/BS01F10_rh.bam tmpBS01F10/BS01F10_rh.dedup_st.bam
/opt/samtools/samtools-0.1.19/samtools index tmpBS01F10/BS01F10_rh.dedup_st.bam
ln -s tmpBS01F10/BS01F10_rh.dedup_st.bam.bai tmpBS01F10/BS01F10_rh.dedup_st.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'finished dedupping using samtools, produced BAM file tmpBS01F10/BS01F10_rh.dedup_st.bam: '$DATE  >>tmpBS01F10/BS01F10.log
FSIZE=`stat --printf="%s" tmpBS01F10/BS01F10_rh.dedup_st.bam`; echo "size of file tmpBS01F10/BS01F10_rh.dedup_st.bam is "$FSIZE  >>tmpBS01F10/BS01F10.log
# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpBS01F10/BS01F10_rh.dedup_st.bam -o tmpBS01F10/BS01F10_rh.dedup_st.reA.intervals
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpBS01F10/BS01F10_rh.dedup_st.bam -targetIntervals tmpBS01F10/BS01F10_rh.dedup_st.reA.intervals -o tmpBS01F10/BS01F10_rh.dedup_st.reA.bam
ln -s tmpBS01F10/BS01F10_rh.dedup_st.reA.bai tmpBS01F10/BS01F10_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'finished re-aligning, produced BAM file tmpBS01F10/BS01F10_rh.dedup_st.reA.bam: '$DATE  >>tmpBS01F10/BS01F10.log
FSIZE=`stat --printf="%s" tmpBS01F10/BS01F10_rh.dedup_st.reA.bam`; echo "size of file tmpBS01F10/BS01F10_rh.dedup_st.reA.bam is "$FSIZE  >>tmpBS01F10/BS01F10.log
# Recalibration of BAM using GATK-BaseRecalibrator+PrintReads
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nct 4 -T BaseRecalibrator -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpBS01F10/BS01F10_rh.dedup_st.reA.bam -knownSites /media/InternBkp1/repos/dbSNP/Ssc_dbSNP138.vcf -o tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.grp
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nct 4 -T PrintReads -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpBS01F10/BS01F10_rh.dedup_st.reA.bam -BQSR tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.grp -o tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam
ln -s tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bai tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'finished re-aligning, produced BAM file tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam: '$DATE  >>tmpBS01F10/BS01F10.log
FSIZE=`stat --printf="%s" tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam`; echo "size of file tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam is "$FSIZE  >>tmpBS01F10/BS01F10.log
# old-school variant calling using the pileup algorithm
echo 'old-school variant calling using the pileup algorithm'
/opt/samtools/samtools-0.1.12a/samtools view -u tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.12a/samtools pileup -vcf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa - >tmpBS01F10/BS01F10_rh.dedup_st.reA.recal_vars-raw.txt
VAR=`cat tmpBS01F10/BS01F10_rh.dedup_st.reA.recal_vars-raw.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/opt/samtools/samtools-0.1.12a/misc/samtools.pl varFilter -D$VAR tmpBS01F10/BS01F10_rh.dedup_st.reA.recal_vars-raw.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.vars-flt_final.txt
# variant calling using the mpileup function of samtools
/opt/samtools/samtools-0.1.19/samtools mpileup -C50 -ugf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.19/bcftools/bcftools view -bvcg -| /opt/samtools/samtools-0.1.19/bcftools/bcftools view - | perl /opt/samtools/samtools-0.1.19/bcftools/vcfutils.pl varFilter -D 20 -d 4 >tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.var.mpileup.flt.vcf
# Variant calling using GATK UnifiedGenotyper - parameters need tweaking
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -T UnifiedGenotyper -I tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam --dbsnp /media/InternBkp1/repos/dbSNP/Ssc_dbSNP138.vcf --genotype_likelihoods_model BOTH -o tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.UG.raw.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0  -dcov 200
# Create gVCF file using modified GATK UnifiedGenotyper - parameters need tweaking
/opt/gvcftools/v0.13-2-gd92e721/bin/getBamAvgChromDepth.pl tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam >tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.avgdepth.txt
java7 -jar /opt/GATK/GATK_gVCFmod/GenomeAnalysisTK.jar -nt 4 -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -T UnifiedGenotyper -I tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam -glm BOTH -l OFF -stand_call_conf 20.0 -stand_emit_conf 10.0  -dcov 200  -out_mode EMIT_ALL_SITES | /opt/gvcftools/v0.13-2-gd92e721/bin/gatk_to_gvcf --chrom-depth-file tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.avgdepth.txt | bgzip -c >tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.gvcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpBS01F10/BS01F10.log; echo 'finished variant calling: '$DATE  >>tmpBS01F10/BS01F10.log
# course nucleotide diversity stat generator
VAR=`cat tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.vars-flt_final.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/opt/samtools/samtools-0.1.12a/samtools view -u tmpBS01F10/BS01F10_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.12a/samtools pileup -f /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -c - | awk '$8>4' | awk -v VAR=$VAR '$8<VAR' | perl /opt/abgsascripts/extract_stats-pileup-bins_allchroms.pl -f tmpBS01F10/BS01F10_rh.dedup_st.reA.recal
