#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h_vmem=10G
mkdir tmpTA01F19
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'starting time: '$DATE  >>tmpTA01F19/TA01F19.log
# archive number 1: ABGSA0063
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0063/ABGSA0063_TA01F19_R1.PF.fastq.gz | pigz >tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0063/ABGSA0063_TA01F19_R2.PF.fastq.gz | pigz >tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.gz -r tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.gz -o tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr -p tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr -s tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.singles.tr -l 45 -t sanger
pigz tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr
pigz tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'starting bwa-mem mapping of TA01F19 archive 1: '$DATE  >>tmpTA01F19/TA01F19.log
# maping using the bwa-mem algorithm, including sorting of bam
echo 'start mapping using BWA-mem algorithm'
/opt/bwa/bwa-0.7.5a/bwa mem -t 4 -M -R '@RG\tID:ABGSA0063_1\tSM:TA01F19\tPL:ILLUMINA' /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr.gz tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr.gz >tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sam
/opt/samtools/samtools-0.1.19/samtools view -Shb -q 10 tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sam > tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.bam
rm tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sam
echo 'start sorting'
/opt/samtools/samtools-0.1.19/samtools sort tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.bam tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sorted
rm tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'finished, produced BAM file tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sorted.bam archive 1: '$DATE  >>tmpTA01F19/TA01F19.log
FSIZE=`stat --printf="%s" tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sorted.bam`; echo "size of file tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sorted.bam is "$FSIZE  >>tmpTA01F19/TA01F19.log
# archive number 2: ABGSA0073
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0073/ABGSA0073_TA01F19_R1.PF.fastq.gz | pigz >tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0073/ABGSA0073_TA01F19_R2.PF.fastq.gz | pigz >tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.gz -r tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.gz -o tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr -p tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr -s tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.singles.tr -l 45 -t sanger
pigz tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr
pigz tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'starting bwa-mem mapping of TA01F19 archive 2: '$DATE  >>tmpTA01F19/TA01F19.log
# maping using the bwa-mem algorithm, including sorting of bam
echo 'start mapping using BWA-mem algorithm'
/opt/bwa/bwa-0.7.5a/bwa mem -t 4 -M -R '@RG\tID:ABGSA0073_2\tSM:TA01F19\tPL:ILLUMINA' /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr.gz tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr.gz >tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sam
/opt/samtools/samtools-0.1.19/samtools view -Shb -q 10 tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sam > tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.bam
rm tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sam
echo 'start sorting'
/opt/samtools/samtools-0.1.19/samtools sort tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.bam tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sorted
rm tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'finished, produced BAM file tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sorted.bam archive 2: '$DATE  >>tmpTA01F19/TA01F19.log
FSIZE=`stat --printf="%s" tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sorted.bam`; echo "size of file tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sorted.bam is "$FSIZE  >>tmpTA01F19/TA01F19.log
# archive number 3: ABGSA0123
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0123/ABGSA0123_TA01F19_R1.PF.fastq.gz | pigz >tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.gz
python /opt/abgsascripts/fix_fq_names.py /media/InternBkp1/repos/ABGSA/ABGSA0123/ABGSA0123_TA01F19_R2.PF.fastq.gz | pigz >tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.gz -r tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.gz -o tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr -p tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr -s tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.singles.tr -l 45 -t sanger
pigz tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr
pigz tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'starting bwa-mem mapping of TA01F19 archive 3: '$DATE  >>tmpTA01F19/TA01F19.log
# maping using the bwa-mem algorithm, including sorting of bam
echo 'start mapping using BWA-mem algorithm'
/opt/bwa/bwa-0.7.5a/bwa mem -t 4 -M -R '@RG\tID:ABGSA0123_3\tSM:TA01F19\tPL:ILLUMINA' /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr.gz tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr.gz >tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sam
/opt/samtools/samtools-0.1.19/samtools view -Shb -q 10 tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sam > tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.bam
rm tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sam
echo 'start sorting'
/opt/samtools/samtools-0.1.19/samtools sort tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.bam tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sorted
rm tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'finished, produced BAM file tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sorted.bam archive 3: '$DATE  >>tmpTA01F19/TA01F19.log
FSIZE=`stat --printf="%s" tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sorted.bam`; echo "size of file tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sorted.bam is "$FSIZE  >>tmpTA01F19/TA01F19.log
#number of bams: 3
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sorted.bam | sed 's/SM:unknown/SM:TA01F19/'  | sed 's/PL:sanger/PL:ILLUMINA/' >tmpTA01F19/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sorted.bam | sed 's/SM:unknown/SM:TA01F19/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpTA01F19/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sorted.bam | sed 's/SM:unknown/SM:TA01F19/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpTA01F19/newheader.txt
/opt/samtools/samtools-0.1.19/samtools merge tmpTA01F19/tmpmergedTA01F19.bam tmpTA01F19/aln-ABGSA0063-TA01F19-1-pe.sorted.bam tmpTA01F19/aln-ABGSA0073-TA01F19-2-pe.sorted.bam tmpTA01F19/aln-ABGSA0123-TA01F19-3-pe.sorted.bam 
/opt/samtools/samtools-0.1.19/samtools reheader tmpTA01F19/newheader.txt tmpTA01F19/tmpmergedTA01F19.bam >tmpTA01F19/TA01F19_rh.bam
rm tmpTA01F19/tmpmergedTA01F19.bam
/opt/samtools/samtools-0.1.19/samtools index tmpTA01F19/TA01F19_rh.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'finished merging, produced BAM file tmpTA01F19/TA01F19_rh.bam: '$DATE  >>tmpTA01F19/TA01F19.log
FSIZE=`stat --printf="%s" tmpTA01F19/TA01F19_rh.bam`; echo "size of file tmpTA01F19/TA01F19_rh.bam is "$FSIZE  >>tmpTA01F19/TA01F19.log
# dedup using samtools
echo 'dedupping using samtools'
/opt/samtools/samtools-0.1.19/samtools rmdup tmpTA01F19/TA01F19_rh.bam tmpTA01F19/TA01F19_rh.dedup_st.bam
/opt/samtools/samtools-0.1.19/samtools index tmpTA01F19/TA01F19_rh.dedup_st.bam
cp tmpTA01F19/TA01F19_rh.dedup_st.bam.bai tmpTA01F19/TA01F19_rh.dedup_st.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'finished dedupping using samtools, produced BAM file tmpTA01F19/TA01F19_rh.dedup_st.bam: '$DATE  >>tmpTA01F19/TA01F19.log
FSIZE=`stat --printf="%s" tmpTA01F19/TA01F19_rh.dedup_st.bam`; echo "size of file tmpTA01F19/TA01F19_rh.dedup_st.bam is "$FSIZE  >>tmpTA01F19/TA01F19.log
# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpTA01F19/TA01F19_rh.dedup_st.bam -o tmpTA01F19/TA01F19_rh.dedup_st.reA.intervals
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpTA01F19/TA01F19_rh.dedup_st.bam -targetIntervals tmpTA01F19/TA01F19_rh.dedup_st.reA.intervals -o tmpTA01F19/TA01F19_rh.dedup_st.reA.bam
cp tmpTA01F19/TA01F19_rh.dedup_st.reA.bai tmpTA01F19/TA01F19_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'finished re-aligning, produced BAM file tmpTA01F19/TA01F19_rh.dedup_st.reA.bam: '$DATE  >>tmpTA01F19/TA01F19.log
FSIZE=`stat --printf="%s" tmpTA01F19/TA01F19_rh.dedup_st.reA.bam`; echo "size of file tmpTA01F19/TA01F19_rh.dedup_st.reA.bam is "$FSIZE  >>tmpTA01F19/TA01F19.log
# Recalibration of BAM using GATK-BaseRecalibrator+PrintReads
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nct 4 -T BaseRecalibrator -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpTA01F19/TA01F19_rh.dedup_st.reA.bam -knownSites /media/InternBkp1/repos/dbSNP/Ssc_dbSNP138.vcf -o tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.grp
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nct 4 -T PrintReads -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpTA01F19/TA01F19_rh.dedup_st.reA.bam -BQSR tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.grp -o tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam
cp tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bai tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'finished re-aligning, produced BAM file tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam: '$DATE  >>tmpTA01F19/TA01F19.log
FSIZE=`stat --printf="%s" tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam`; echo "size of file tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam is "$FSIZE  >>tmpTA01F19/TA01F19.log
# old-school variant calling using the pileup algorithm
echo 'old-school variant calling using the pileup algorithm'
/opt/samtools/samtools-0.1.12a/samtools view -u tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.12a/samtools pileup -vcf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa - >tmpTA01F19/TA01F19_rh.dedup_st.reA.recal_vars-raw.txt
VAR=`cat tmpTA01F19/TA01F19_rh.dedup_st.reA.recal_vars-raw.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/opt/samtools/samtools-0.1.12a/misc/samtools.pl varFilter -D$VAR tmpTA01F19/TA01F19_rh.dedup_st.reA.recal_vars-raw.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.vars-flt_final.txt
# variant calling using the mpileup function of samtools
/opt/samtools/samtools-0.1.19/samtools mpileup -C50 -ugf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.19/bcftools/bcftools view -bvcg -| /opt/samtools/samtools-0.1.19/bcftools/bcftools view - | perl /opt/samtools/samtools-0.1.19/bcftools/vcfutils.pl varFilter -D 20 -d 4 >tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.var.mpileup.flt.vcf
# Variant calling using GATK UnifiedGenotyper - parameters need tweaking
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -T UnifiedGenotyper -I tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam --dbsnp /media/InternBkp1/repos/dbSNP/Ssc_dbSNP138.vcf --genotype_likelihoods_model BOTH -o tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.UG.raw.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0  -dcov 200
bgzip tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.UG.raw.vcf
tabix -p vcf tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.UG.raw.vcf.gz
# Create gVCF file using modified GATK UnifiedGenotyper - parameters need tweaking
/opt/gvcftools/v0.13-2-gd92e721/bin/getBamAvgChromDepth.pl tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam >tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.avgdepth.txt
java7 -jar /opt/GATK/GATK_gVCFmod/GenomeAnalysisTK.jar -nt 4 -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -T UnifiedGenotyper -I tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam -glm BOTH -l OFF -stand_call_conf 20.0 -stand_emit_conf 10.0  -dcov 200  -out_mode EMIT_ALL_SITES | /opt/gvcftools/v0.13-2-gd92e721/bin/gatk_to_gvcf --chrom-depth-file tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.avgdepth.txt | bgzip -c >tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.gvcf.gz
tabix -p vcf tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.gvcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpTA01F19/TA01F19.log; echo 'finished variant calling: '$DATE  >>tmpTA01F19/TA01F19.log
# predicting function using VEP
gunzip -c tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.UG.raw.vcf.gz | perl /opt/VEP/variant_effect_predictor.pl --dir /opt/VEP/ --species sus_scrofa -o tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.UG.raw.vep.txt --fork 4 --canonical --sift b --coding_only --no_intergenic --offline --force_overwrite
# course nucleotide diversity stat generator
VAR=`cat tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.vars-flt_final.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/opt/samtools/samtools-0.1.12a/samtools view -u tmpTA01F19/TA01F19_rh.dedup_st.reA.recal.bam | /opt/samtools/samtools-0.1.12a/samtools pileup -f /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -c - | awk '$8>4' | awk -v VAR=$VAR '$8<VAR' | perl /opt/abgsascripts/extract_stats-pileup-bins_allchroms.pl -f tmpTA01F19/TA01F19_rh.dedup_st.reA.recal
