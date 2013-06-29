#!/bin/bash
#$ -cwd
#$ -S /bin/sh
#$ -l h_vmem=10G
mkdir tmpTA01F19
# archive number 1: ABGSA0063
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0063/ABGSA0063_TA01F19_R1.PF.fastq.gz | sed 's/ /#/' | pigz >tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.gz
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0063/ABGSA0063_TA01F19_R2.PF.fastq.gz | sed 's/ /#/' | pigz >tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.gz -r tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.gz -o tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr -p tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr -s tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.singles.tr -l 45 -t sanger
pigz tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr
pigz tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr.gz  >tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr.gz  >tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -r '@RG\tID:ABGSA0063_1\tSM:TA01F19' tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr.sai tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr.sai tmpTA01F19/ABGSA0063_TA01F19_R1.PF.fastq.tr.gz tmpTA01F19/ABGSA0063_TA01F19_R2.PF.fastq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpTA01F19/aln-ABGSA0063-TA01F19-1PE2.sorted
# archive number 2: ABGSA0073
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0073/ABGSA0073_TA01F19_R1.PF.fastq.gz | sed 's/ /#/' | pigz >tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.gz
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0073/ABGSA0073_TA01F19_R2.PF.fastq.gz | sed 's/ /#/' | pigz >tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.gz -r tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.gz -o tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr -p tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr -s tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.singles.tr -l 45 -t sanger
pigz tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr
pigz tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr.gz  >tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr.gz  >tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -r '@RG\tID:ABGSA0073_2\tSM:TA01F19' tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr.sai tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr.sai tmpTA01F19/ABGSA0073_TA01F19_R1.PF.fastq.tr.gz tmpTA01F19/ABGSA0073_TA01F19_R2.PF.fastq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpTA01F19/aln-ABGSA0073-TA01F19-2PE2.sorted
# archive number 3: ABGSA0123
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0123/ABGSA0123_TA01F19_R1.PF.fastq.gz | sed 's/ /#/' | pigz >tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.gz
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0123/ABGSA0123_TA01F19_R2.PF.fastq.gz | sed 's/ /#/' | pigz >tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.gz -r tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.gz -o tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr -p tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr -s tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.singles.tr -l 45 -t sanger
pigz tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr
pigz tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr.gz  >tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr.gz  >tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -r '@RG\tID:ABGSA0123_3\tSM:TA01F19' tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr.sai tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr.sai tmpTA01F19/ABGSA0123_TA01F19_R1.PF.fastq.tr.gz tmpTA01F19/ABGSA0123_TA01F19_R2.PF.fastq.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpTA01F19/aln-ABGSA0123-TA01F19-3PE2.sorted
#number of bams: 3
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpTA01F19/aln-ABGSA0063-TA01F19-1PE2.sorted.bam | sed 's/SM:unknown/SM:TA01F19/' >tmpTA01F19/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpTA01F19/aln-ABGSA0073-TA01F19-2PE2.sorted.bam | sed 's/SM:unknown/SM:TA01F19/' | grep @RG >>tmpTA01F19/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpTA01F19/aln-ABGSA0123-TA01F19-3PE2.sorted.bam | sed 's/SM:unknown/SM:TA01F19/' | grep @RG >>tmpTA01F19/newheader.txt
/opt/samtools/samtools-0.1.19/samtools merge tmpTA01F19/tmpmergedTA01F19.bam tmpTA01F19/aln-ABGSA0063-TA01F19-1PE2.sorted.bam tmpTA01F19/aln-ABGSA0073-TA01F19-2PE2.sorted.bam tmpTA01F19/aln-ABGSA0123-TA01F19-3PE2.sorted.bam 
/opt/samtools/samtools-0.1.19/samtools reheader tmpTA01F19/newheader.txt tmpTA01F19/tmpmergedTA01F19.bam >tmpTA01F19/TA01F19_rh.bam
rm tmpTA01F19/tmpmergedTA01F19.bam
# dedup using Picard
echo 'dedupping using picard MarkDuplicates'
java7 -Xmx4g -jar /opt/picard/picard-tools-1.93/MarkDuplicates.jar ASSUME_SORTED=true REMOVE_DUPLICATES=true INPUT=tmpTA01F19/TA01F19_rh.bam OUTPUT=tmpTA01F19/TA01F19_rh.dedup_pi.bam METRICS_FILE=tmpTA01F19/TA01F19_rh.dedup.metrics
/opt/samtools/samtools-0.1.19/samtools sort tmpTA01F19/TA01F19_rh.dedup_pi.bam tmpTA01F19/TA01F19_rh.dedup_pi.sorted
rm tmpTA01F19/TA01F19_rh.dedup_pi.bam
mv tmpTA01F19/TA01F19_rh.dedup_pi.sorted.bam tmpTA01F19/TA01F19_rh.dedup_pi.bam
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpTA01F19/TA01F19_rh.dedup_pi.bam -o tmpTA01F19/TA01F19_rh.dedup_pi.reA.intervals
java7 -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpTA01F19/TA01F19_rh.dedup_pi.bam -targetIntervals tmpTA01F19/TA01F19_rh.dedup_pi.reA.intervals -o tmpTA01F19/TA01F19_rh.dedup_pi.reA.bam
/opt/samtools/samtools-0.1.19/samtools sort tmpTA01F19/TA01F19_rh.dedup_pi.reA.bam tmpTA01F19/TA01F19_rh.dedup_pi.reA.sorted
rm tmpTA01F19/TA01F19_rh.dedup_pi.reA.bam
mv tmpTA01F19/TA01F19_rh.dedup_pi.reA.sorted.bam tmpTA01F19/TA01F19_rh.dedup_pi.reA.bam
java -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpTA01F19/TA01F19_rh.dedup_pi.reA.bam -knownSites /path/to/dbsnp/file.vcf -o tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.grp
java -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -T PrintReads -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpTA01F19/TA01F19_rh.dedup_pi.reA.bam -BQSR tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.grp -o tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.bam
# old-school variant calling using the pileup algorithm
echo 'old-school variant calling using the pileup algorithm'
/opt/samtools/samtools-0.1.12a/samtools view -u tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.bam | /opt/samtools/samtools-0.1.12a/samtools pileup -vcf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa - >tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal_vars-raw.txt
VAR=`cat tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal_vars-raw.txt | cut -f8 | head -100000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/opt/samtools/samtools-0.1.12a/misc/samtools.pl varFilter -D$VAR tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal_vars-raw.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.vars-flt_final.txt
# variant calling using the mpileup function of samtools
/opt/samtools/samtools-0.1.19/samtools mpileup -C50 -ugf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.bam | /opt/samtools/samtools-0.1.19/bcftools/bcftools view -bvcg -| /opt/samtools/samtools-0.1.19/bcftools/bcftools view - | perl /opt/samtools/samtools-0.1.19/bcftools/bcftools/vcfutils.pl varFilter -D 20 -d 4 >tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.var.mpileup.flt.vcf
java -jar /opt/GATK/GATK2.6/GenomeAnalysisTK.jar -R /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -T UnifiedGenotyper -I tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.bam --dbsnp /path/to/dbsnp/file.vcf --genotype_likelihoods_model BOTH -o tmpTA01F19/TA01F19_rh.dedup_pi.reA.recal.UG.raw.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0  -dcov 50
