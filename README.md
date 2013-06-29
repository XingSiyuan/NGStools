Single-sample oriented version of the new ABGC mapping pipeline.
Very early, rough draft - handle with care: many components are still to be validated.
Note that access to sample database (MySQL) and access to primary data is required for the script to run!

The python script will create a shell script that contains various steps to process and map fastq data. The python scripts takes an individual id as argument (-1). It will also check if data is in Illumina or Sanger (offset +64 or +32) quality coding, and add translation to sanger if needed. For an example of the resulting shell script, see 'testscript1.sh'. 

example usage:
% python3 ABGC_mapping_v2.py -i LW22F08 -a /media/InternBkp1/repos/ABGSA/ -r /media/InternBkp1/repos/refs/Sus_scrofa.Sscrofa10.2.72.dna.toplevel.fa

This will create a file 'runLW22F08.sh', which can be submitted using SGE (qsub)

for more information:
% python3 ABGC_mapping_v2.py -h

The script is designed to:
- pull data from the repositories
- prepare data for mapping:
  * create temporary repository
  * trimming by sickle
  * convert from Illumina fq to Sanger fq, if needed (automated detection)
- mapping using BWA mem (currently only option, more mapping options will follow)
- sort BAM files
- merge all BAM files of same individual, and reheader
- do old-school variant calling using samtools pileup

Still to validate:
  * filtering of BAM files based on mapping Q
  * local re-alignment using GATK
  * dedup (Picard)
  * dedup (samtools)
  * recalibration
  * variant calling (samtools mpileup)
  * variant calling (GATK UnifiedGenotyper)

To do:
  * do intersection of variants by Samtools and GATK (maybe)
  * do Variant Effect Predictor (VEP)


