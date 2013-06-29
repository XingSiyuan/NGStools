#!/bin/bash
#$ -cwd
#$ -S /bin/sh
#$ -l h_vmem=10G
mkdir tmpLW22F08
# archive number 1: ABGSA0189
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0189/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz | sed 's/ /#/' | pigz >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0189/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz | sed 's/ /#/' | pigz >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz
sickle pe -f tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz -r tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz -o tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr -p tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr -s tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.singles.tr -l 45 -t illumina
pigz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr
pigz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python convert_ill_to_sang.py tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz | gzip -c >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sa.gz
python convert_ill_to_sang.py tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz | gzip -c >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sa.gz
rm tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz
rm tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz
mv tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sa.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz
mv tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sa.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz  >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz  >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -r '@RG\tID:ABGSA0189_1\tSM:LW22F08' tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sai tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sai tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpLW22F08/aln-ABGSA0189-LW22F08-1PE2.sorted
# archive number 2: ABGSA0189
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0189/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz | sed 's/ /#/' | pigz >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0189/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz | sed 's/ /#/' | pigz >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz
sickle pe -f tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz -r tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz -o tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr -p tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr -s tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.singles.tr -l 45 -t illumina
pigz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr
pigz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python convert_ill_to_sang.py tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz | gzip -c >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sa.gz
python convert_ill_to_sang.py tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz | gzip -c >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sa.gz
rm tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz
rm tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz
mv tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sa.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz
mv tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sa.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz
# maping using the bwa-aln algorithm, including sorting of bam
echo 'start mapping using BWA-aln algorithm'
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz  >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sai
/opt/bwa/bwa-0.7.5a/bwa aln -n 0.07 -t 4 /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz  >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sai
/opt/bwa/bwa-0.7.5a/bwa sampe -P /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa -r '@RG\tID:ABGSA0189_2\tSM:LW22F08' tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sai tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sai tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -q 20 -Suh - | /opt/samtools/samtools-0.1.19/samtools sort -m 5000000000 - tmpLW22F08/aln-ABGSA0189-LW22F08-2PE2.sorted
#number of bams: 2
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpLW22F08/aln-ABGSA0189-LW22F08-1PE2.sorted | sed 's/SM:unknown/SM:LW22F08/' >tmpLW22F08/newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpLW22F08/aln-ABGSA0189-LW22F08-2PE2.sorted | sed 's/SM:unknown/SM:LW22F08/' | grep @RG >>tmpLW22F08/newheader.txt
/opt/samtools/samtools-0.1.19/samtools merge tmpLW22F08/tmpmergedLW22F08.bam tmpLW22F08/aln-ABGSA0189-LW22F08-1PE2.sorted tmpLW22F08/aln-ABGSA0189-LW22F08-2PE2.sorted 
/opt/samtools/samtools-0.1.19/samtools reheader tmpLW22F08/newheader.txt tmpLW22F08/tmpmergedLW22F08.bam >tmpLW22F08/LW22F08_rh.bam
rm tmpLW22F08/tmpmergedLW22F08.bam
# dedup using Picard
echo 'dedupping using picard MarkDuplicates
java7 -Xmx4g -jar /opt/picard/picard-tools-1.93/MarkDuplicates.jar ASSUME_SORTED=true REMOVE_DUPLICATES=true INPUT=tmpLW22F08/LW22F08_rh.bam OUTPUT=tmpLW22F08/LW22F08_rh.dedup_pi.bam METRICS_FILE=tmpLW22F08/LW22F08_rh.dedup.metrics
/opt/samtools/samtools-0.1.19/samtools sort tmpLW22F08/LW22F08_rh.dedup_pi.bam tmpLW22F08/LW22F08_rh.dedup_pi.sorted
rm tmpLW22F08/LW22F08_rh.dedup_pi.bam
mv tmpLW22F08/LW22F08_rh.dedup_pi.sorted.bam tmpLW22F08/LW22F08_rh.dedup.bam
# old-school variant calling using the pileup algorithm
echo 'old-school variant calling using the pileup algorithm'
/opt/samtools/samtools-0.1.12a/samtools view -u tmpLW22F08/LW22F08_rh.dedup_pi.bam | /opt/samtools/samtools-0.1.12a/samtools pileup -vcf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa - >tmpLW22F08/LW22F08_rh.dedup_pi_vars-raw.txt
VAR=`cat tmpLW22F08/LW22F08_rh.dedup_pi_vars-raw.txt | cut -f8 | head -100000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/opt/samtools/samtools-0.1.12a/misc/samtools.pl varFilter -D$VAR tmpLW22F08/LW22F08_rh.dedup_pi_vars-raw.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >tmpLW22F08/LW22F08_rh.dedup_pi_vars-flt_final.txt
# variant calling using the mpileup function of samtools
/opt/samtools/samtools-0.1.19/samtools mpileup -C50 -ugf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/LW22F08_rh.dedup_pi.bam | /opt/samtools/samtools-0.1.19/bcftools/bcftools view -bvcg -| /opt/samtools/samtools-0.1.19/bcftools/bcftools view - | perl /opt/samtools/samtools-0.1.19/bcftools/bcftools/vcfutils.pl varFilter -D 20 -d 4 >tmpLW22F08/LW22F08_rh.dedup_pi.var.flt.vcf
