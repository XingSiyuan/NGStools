#!/usr/bin/perl -w
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Std;
my %opts = ();
$opts{f}= 'final';
$opts{d}= 'dir';
$opts{c}= 'chrom';
$opts{m}= 2;
getopt('dfcm', \%opts);
my $filesfile = $opts{f};
my $dir = $opts{d};
my $chrom = $opts{c};
my $maxdepth = $opts{m};
open (IN, "$filesfile") or die "can not open file $filesfile $!\n";
my @files;
while (<IN>){
        my $file = $_;
        chomp $file;
        push(@files,$file);
}
close(IN);
for my $file (@files){
        my $stub = $file;
        $stub =~ s/_rh\.bam//;
        open (OUT, ">run_".$stub."_create".$chrom."_fq.txt");
        print OUT '#!/bin/bash'."\n";
        print OUT '#$ -cwd'."\n";
        print OUT '#$ -S /bin/sh'."\n";
        print OUT ' VAR=`cat /srv/nfs02/shared/Sus/vars_hjm_newbuild10_2/vars_pileup/vars-flt_'.$stub.'-final.txt | grep '.$chrom.' | cut -f8 | head -100000 | sort | uniq -c | sed '."'".'s/^ \+//'."'".' | sed '."'".'s/ \+/\t/'."'".' | sort -k1 -nr | head -1 | cut -f2`'."\n";
        print OUT 'let VAR='.$maxdepth.'*VAR'."\n";
        print OUT 'echo "max depth is $VAR"'."\n";
        print OUT 'samtools view -u /srv/nfs02/shared/Sus/BAM_files_hjm_newbuild10_2/'.$stub.'_rh.bam '.$chrom.' | samtools pileup -cf /srv/nfs02/shared/Sus/Sscrofa_build10_2/FASTA/Ssc10_2_wUn.fa -c - | awk '."'".'$6!~/[ID]/'."'".' | samtools.pl pileup2fq -D$VAR -d2 >'.$stub."_".$chrom.".fq\n";
        close(OUT);
        my $qsubfile = 'run_'.$stub.'_create'.$chrom.'_fq.txt';
        `qsub -q all.q $qsubfile`;

}
exit;

