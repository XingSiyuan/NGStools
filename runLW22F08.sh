#!/bin/bash
#SBATCH --time=10000
#SBATCH --mem=16000
#SBATCH --ntasks=4
#SBATCH --output=outputLW22F08_%j.txt
#SBATCH --error=error_outputLW22F08_%j.txt
#SBATCH --job-name=LW22F08
#SBATCH --partition=ABGC_Research
mkdir tmpLW22F08
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'starting time: '$DATE  >>tmpLW22F08/LW22F08.log
# archive number 1: ABGSA0189
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Pig/ABGSA/ABGSA0189/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz | pigz >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Pig/ABGSA/ABGSA0189/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz | pigz >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz -r tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz -o tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.tr -p tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.tr -s tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.singles.tr -l 50 -t illumina
pigz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.tr
pigz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/convert_ill_to_sang.py tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.tr.gz | gzip -c >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.tr.sa.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/convert_ill_to_sang.py tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.tr.gz | gzip -c >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.tr.sa.gz
rm tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.tr.gz
rm tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.tr.gz
mv tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.tr.sa.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.tr.gz
mv tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.tr.sa.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'starting bwa-mem mapping of LW22F08 archive 1: '$DATE  >>tmpLW22F08/LW22F08.log
# maping using the bwa-mem algorithm, including sorting of bam
echo 'start mapping using BWA-mem algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa mem -t 4 -M -R '@RG\tID:ABGSA0189_1\tSM:LW22F08\tPL:ILLUMINA' /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.tr.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.tr.gz >tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Shb tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sam > tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.bam
rm tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sam
echo 'start sorting'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.bam tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted
rm tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished, produced BAM file tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted.bam archive 1: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted.bam`; echo "size of file tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted.bam is "$FSIZE  >>tmpLW22F08/LW22F08.log
# archive number 2: ABGSA0189
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Pig/ABGSA/ABGSA0189/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz | pigz >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/nobackup/WUR/ABGC/shared/Pig/ABGSA/ABGSA0189/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz | pigz >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz
# quality trimming of reads by sickle
sickle pe -f tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz -r tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz -o tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.tr -p tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.tr -s tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.singles.tr -l 50 -t illumina
pigz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.tr
pigz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/convert_ill_to_sang.py tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.tr.gz | gzip -c >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.tr.sa.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/convert_ill_to_sang.py tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.tr.gz | gzip -c >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.tr.sa.gz
rm tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.tr.gz
rm tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.tr.gz
mv tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.tr.sa.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.tr.gz
mv tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.tr.sa.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.tr.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'starting bwa-mem mapping of LW22F08 archive 2: '$DATE  >>tmpLW22F08/LW22F08.log
# maping using the bwa-mem algorithm, including sorting of bam
echo 'start mapping using BWA-mem algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa mem -t 4 -M -R '@RG\tID:ABGSA0189_2\tSM:LW22F08\tPL:ILLUMINA' /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.tr.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.tr.gz >tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Shb tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sam > tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.bam
rm tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sam
echo 'start sorting'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.bam tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted
rm tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished, produced BAM file tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted.bam archive 2: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted.bam`; echo "size of file tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted.bam is "$FSIZE  >>tmpLW22F08/LW22F08.log
#number of bams: 2
#multiple bam files --> do merge
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted.bam | sed 's/SM:unknown/SM:LW22F08/'  | sed 's/PL:sanger/PL:ILLUMINA/' >tmpLW22F08/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted.bam | sed 's/SM:unknown/SM:LW22F08/'  | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpLW22F08/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools merge tmpLW22F08/tmpmergedLW22F08.bam tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted.bam tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted.bam 
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools reheader tmpLW22F08/newheader.txt tmpLW22F08/tmpmergedLW22F08.bam >tmpLW22F08/LW22F08_rh.bam
rm tmpLW22F08/tmpmergedLW22F08.bam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools index tmpLW22F08/LW22F08_rh.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished re-aligning, produced BAM file tmpLW22F08/LW22F08_rh.bam: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/LW22F08_rh.bam`; echo "size of file tmpLW22F08/LW22F08_rh.bam is "$FSIZE  >>tmpLW22F08/LW22F08.log
MD5BAM=`md5sum tmpLW22F08/LW22F08_rh.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpLW22F08/LW22F08_rh.bam is "$MD5BAM  >>tmpLW22F08/LW22F08.log
# dedup using samtools
echo 'dedupping using samtools'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools rmdup tmpLW22F08/LW22F08_rh.bam tmpLW22F08/LW22F08_rh.dedup_st.bam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools index tmpLW22F08/LW22F08_rh.dedup_st.bam
cp tmpLW22F08/LW22F08_rh.dedup_st.bam.bai tmpLW22F08/LW22F08_rh.dedup_st.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished re-aligning, produced BAM file tmpLW22F08/LW22F08_rh.dedup_st.bam: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/LW22F08_rh.dedup_st.bam`; echo "size of file tmpLW22F08/LW22F08_rh.dedup_st.bam is "$FSIZE  >>tmpLW22F08/LW22F08.log
MD5BAM=`md5sum tmpLW22F08/LW22F08_rh.dedup_st.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpLW22F08/LW22F08_rh.dedup_st.bam is "$MD5BAM  >>tmpLW22F08/LW22F08.log
# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner
java7 -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpLW22F08/LW22F08_rh.dedup_st.bam -o tmpLW22F08/LW22F08_rh.dedup_st.reA.intervals
java7 -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpLW22F08/LW22F08_rh.dedup_st.bam -targetIntervals tmpLW22F08/LW22F08_rh.dedup_st.reA.intervals -o tmpLW22F08/LW22F08_rh.dedup_st.reA.bam
cp tmpLW22F08/LW22F08_rh.dedup_st.reA.bai tmpLW22F08/LW22F08_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished re-aligning, produced BAM file tmpLW22F08/LW22F08_rh.dedup_st.reA.bam: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/LW22F08_rh.dedup_st.reA.bam`; echo "size of file tmpLW22F08/LW22F08_rh.dedup_st.reA.bam is "$FSIZE  >>tmpLW22F08/LW22F08.log
MD5BAM=`md5sum tmpLW22F08/LW22F08_rh.dedup_st.reA.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpLW22F08/LW22F08_rh.dedup_st.reA.bam is "$MD5BAM  >>tmpLW22F08/LW22F08.log
# Calculate coverage statistics
java7 -Xmx8G -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -T DepthOfCoverage -R /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -I tmpLW22F08/LW22F08_rh.dedup_st.reA.bam --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o tmpLW22F08/LW22F08_rh.dedup_st.reA.bam.coverage
# old-school variant calling using the pileup algorithm
echo 'old-school variant calling using the pileup algorithm'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools view -u tmpLW22F08/LW22F08_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools pileup -vcf /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa - >tmpLW22F08/LW22F08_rh.dedup_st.reA_vars-raw.txt
VAR=`cat tmpLW22F08/LW22F08_rh.dedup_st.reA_vars-raw.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/misc/samtools.pl varFilter -D$VAR tmpLW22F08/LW22F08_rh.dedup_st.reA_vars-raw.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >tmpLW22F08/LW22F08_rh.dedup_st.reA.vars-flt_final.txt
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished re-aligning, produced BAM file tmpLW22F08/LW22F08_rh.dedup_st.reA.vars-flt_final.txt: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/LW22F08_rh.dedup_st.reA.vars-flt_final.txt`; echo "size of file tmpLW22F08/LW22F08_rh.dedup_st.reA.vars-flt_final.txt is "$FSIZE  >>tmpLW22F08/LW22F08.log
MD5VAR=`md5sum tmpLW22F08/LW22F08_rh.dedup_st.reA.vars-flt_final.txt | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpLW22F08/LW22F08_rh.dedup_st.reA.vars-flt_final.txt is "$MD5VAR  >>tmpLW22F08/LW22F08.log
# variant calling using the mpileup function of samtools
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools mpileup -C50 -ugf /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/LW22F08_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/bcftools view -bvcg -| /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/bcftools view - | perl /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/vcfutils.pl varFilter -D 20 -d 4 >tmpLW22F08/LW22F08_rh.dedup_st.reA.var.mpileup.flt.vcf
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished re-aligning, produced BAM file tmpLW22F08/LW22F08_rh.dedup_st.reA.var.mpileup.flt.vcf: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/LW22F08_rh.dedup_st.reA.var.mpileup.flt.vcf`; echo "size of file tmpLW22F08/LW22F08_rh.dedup_st.reA.var.mpileup.flt.vcf is "$FSIZE  >>tmpLW22F08/LW22F08.log
MD5VAR=`md5sum tmpLW22F08/LW22F08_rh.dedup_st.reA.var.mpileup.flt.vcf | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpLW22F08/LW22F08_rh.dedup_st.reA.var.mpileup.flt.vcf is "$MD5VAR  >>tmpLW22F08/LW22F08.log
# Variant calling using GATK UnifiedGenotyper - parameters need tweaking
java7 -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -R /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -T UnifiedGenotyper -I tmpLW22F08/LW22F08_rh.dedup_st.reA.bam --dbsnp /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/dbSNP/dbSNP.vcf --genotype_likelihoods_model BOTH -o tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0  -dcov 200
bgzip tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf
tabix -p vcf tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished re-aligning, produced BAM file tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf.gz: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf.gz`; echo "size of file tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf.gz is "$FSIZE  >>tmpLW22F08/LW22F08.log
MD5VAR=`md5sum tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf.gz | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf.gz is "$MD5VAR  >>tmpLW22F08/LW22F08.log
# Create gVCF file using modified GATK UnifiedGenotyper - parameters need tweaking
/cm/shared/apps/WUR/ABGC/gvcftools/gvcftools-0.16/bin/getBamAvgChromDepth.pl tmpLW22F08/LW22F08_rh.dedup_st.reA.bam >tmpLW22F08/LW22F08_rh.dedup_st.reA.avgdepth.txt
java7 -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK_gVCFmod/GenomeAnalysisTK.jar -nt 4 -R /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -T UnifiedGenotyper -I tmpLW22F08/LW22F08_rh.dedup_st.reA.bam -glm BOTH -l OFF -stand_call_conf 20.0 -stand_emit_conf 10.0  -dcov 200  -out_mode EMIT_ALL_SITES | /cm/shared/apps/WUR/ABGC/gvcftools/gvcftools-0.16/bin/gatk_to_gvcf --chrom-depth-file tmpLW22F08/LW22F08_rh.dedup_st.reA.avgdepth.txt | bgzip -c >tmpLW22F08/LW22F08_rh.dedup_st.reA.gvcf.gz
tabix -p vcf tmpLW22F08/LW22F08_rh.dedup_st.reA.gvcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished re-aligning, produced BAM file tmpLW22F08/LW22F08_rh.dedup_st.reA.gvcf.gz: '$DATE  >>tmpLW22F08/LW22F08.log
FSIZE=`stat --printf="%s" tmpLW22F08/LW22F08_rh.dedup_st.reA.gvcf.gz`; echo "size of file tmpLW22F08/LW22F08_rh.dedup_st.reA.gvcf.gz is "$FSIZE  >>tmpLW22F08/LW22F08.log
MD5VAR=`md5sum tmpLW22F08/LW22F08_rh.dedup_st.reA.gvcf.gz | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpLW22F08/LW22F08_rh.dedup_st.reA.gvcf.gz is "$MD5VAR  >>tmpLW22F08/LW22F08.log
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpLW22F08/LW22F08.log; echo 'finished variant calling: '$DATE  >>tmpLW22F08/LW22F08.log
# predicting function using VEP
gunzip -c tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vcf.gz | perl /cm/shared/apps/WUR/ABGC/variant_effect_predictor/VEP231213/variant_effect_predictor.pl --dir /cm/shared/apps/WUR/ABGC/variant_effect_predictor/VEP231213//cache --species sus_scrofa -o tmpLW22F08/LW22F08_rh.dedup_st.reA.UG.raw.vep.txt --fork 4 --canonical --sift b --coding_only --no_intergenic --offline --force_overwrite
# course nucleotide diversity stat generator
VAR=`cat tmpLW22F08/LW22F08_rh.dedup_st.reA.vars-flt_final.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools view -u tmpLW22F08/LW22F08_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools pileup -f /lustre/nobackup/WUR/ABGC/shared/Pig/Sscrofa_build10_2/Ensembl72/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -c - | awk '$8>4' | awk -v VAR=$VAR '$8<VAR' | perl /cm/shared/apps/WUR/ABGC/abgsascripts/extract_stats-pileup-bins_allchroms.pl -f tmpLW22F08/LW22F08_rh.dedup_st.reA
