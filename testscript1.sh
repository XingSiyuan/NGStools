#!/bin/bash
#$ -cwd
#$ -S /bin/sh
#$ -l h_vmem=10G
mkdir tmpLW22F08
# archive number 1: ABGSA0189
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0189/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz | sed 's/ /#/' | pigz >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0189/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz | sed 's/ /#/' | pigz >tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz
# quality trimming of reads by sickle
sickle pe -f tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz -r tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz -o tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr -p tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr -s tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.singles.tr -t illumina
pigz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr
pigz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python convert_ill_to_sang.py tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz | gzip -c >tmpLW22F08//120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sa.gz
python convert_ill_to_sang.py tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz | gzip -c >tmpLW22F08//120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sa.gz
rm tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz
rm tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz
mv tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sa.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz
mv tmpLW22F08//120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sa.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz
# maping using the bwa-mem algorithm, including sorting of bam
/opt/bwa/bwa-0.7.5a/bwa mem -t 4 -R '@RG\tID:ABGSA0189_1\tSM:LW22F08' /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz tmpLW22F08/120423_I652_FCC0E33ACXX_L3_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -Shb -q 10 > tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.bam
/opt/samtools/samtools-0.1.19/samtools sort tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.bam tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted
rm tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.bam
# archive number 2: ABGSA0189
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0189/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz | sed 's/ /#/' | pigz >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz
gunzip -c /media/InternBkp1/repos/ABGSA/ABGSA0189/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz | sed 's/ /#/' | pigz >tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz
# quality trimming of reads by sickle
sickle pe -f tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.gz.clean.dup.clean.gz -r tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.gz.clean.dup.clean.gz -o tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr -p tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr -s tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.singles.tr -t illumina
pigz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr
pigz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr
# since sequences have offset +64 we need to convert to sanger (offset +33)
python convert_ill_to_sang.py tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz | gzip -c >tmpLW22F08//120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sa.gz
python convert_ill_to_sang.py tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz | gzip -c >tmpLW22F08//120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sa.gz
rm tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz
rm tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz
mv tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.sa.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz
mv tmpLW22F08//120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.sa.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz
# maping using the bwa-mem algorithm, including sorting of bam
/opt/bwa/bwa-0.7.5a/bwa mem -t 4 -R '@RG\tID:ABGSA0189_2\tSM:LW22F08' /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_1.fq.clean.dup.clean.tr.gz tmpLW22F08/120506_I224_FCC0RYYACXX_L4_SZAIPI008160-111_2.fq.clean.dup.clean.tr.gz | /opt/samtools/samtools-0.1.19/samtools view -Shb -q 10 > tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.bam
/opt/samtools/samtools-0.1.19/samtools sort tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.bam tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted
rm tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.bam
#number of bams: 2
#multiple bam files --> do merge
/opt/samtools/samtools-0.1.19/samtools view -H tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted.bam >newheader.txt
/opt/samtools/samtools-0.1.19/samtools view -H tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted.bam | grep @RG >>newheader.txt
/opt/samtools/samtools-0.1.19/samtools merge tmpLW22F08/tmpmergedtmpLW22F08/LW22F08.bam tmpLW22F08/aln-ABGSA0189-LW22F08-1-pe.sorted.bam tmpLW22F08/aln-ABGSA0189-LW22F08-2-pe.sorted.bam 
/opt/samtools/samtools-0.1.19/samtools reheader tmpLW22F08/newheader.txt tmpLW22F08/tmpmergedLW22F08.bam >tmpLW22F08/LW22F08_rh.bam
rm tmpmergedLW22F08.bam
# old-school variant calling using the pileup algortithm
/opt/samtools/samtools-0.1.12/samtools view -u tmpLW22F08/LW22F08_rh.bam | /opt/samtools/samtools-0.1.12/samtools pileup -vcf /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa - >vars-raw_tmpLW22F08/LW22F08.txt
/opt/samtools/samtools-0.1.12/samtools.pl varFilter -D20 vars-raw_tmpLW22F08/LW22F08.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >vars-flt_tmpLW22F08/LW22F08-final.txt
