#!/usr/bin/perl -w
# Hendrik-Jan Megens, 
# Animal Breeding & Genomics Centre
# Wageningen University
# hendrik-jan.megens -at- wur.nl
# starts nucleotide diversity estimates using the
# 'extract_stats-pileup-bins_allchroms.pl' based on
# a file that contains list of bam files

use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Std;
my %opts = ();
$opts{f}= 'final';
$opts{i}= 'ind';
getopt('dfi', \%opts);
my $filesfile = $opts{f};
my $ind = $opts{i};
open (IN, "$filesfile") or die "can not open file $filesfile $!\n";
my @files;
while (<IN>){
        my $file = $_;
        chomp $file;
        push(@files,$file);
}
close(IN);

my $bamfile='pathtonowhere';

for my $file (@files){
        my $query = $ind.'_rh.bam';
        if ($file =~ m/$query/){

                $bamfile = $file;
        }
}
if ($bamfile eq 'pathtonowhere'){
        print "corresponding BAM file not found for $ind\n";
}

open (OUT, ">run_".$ind."_nucdiv.txt");
print OUT '#!/bin/bash'."\n";
print OUT '#$ -cwd'."\n";
print OUT '#$ -S /bin/sh'."\n";
print OUT '#$ -l h_vmem=5G'."\n";
print OUT ' VAR=`cat /srv/mds01/shared/Sus/vars_hjm_newbuild10_2/vars-flt_'.$ind.'-final.txt | cut -f8 | head -100000 | sort | uniq -c | sed '."'".'s/^ \+//'."'".' | sed '."'".'s/ \+/\t/'."'".' | sort -k1 -nr | head -1 | cut -f2`'."\n";
print OUT 'let VAR=2*VAR'."\n";
print OUT 'echo "max depth is $VAR"'."\n";
print OUT 'samtools view -u '.$bamfile.' | samtools pileup -f /srv/mds01/shared/Sus/Sscrofa_build10_2/FASTA/Ssc10_2_wUn.fa -c - | awk '."'".'$8>4'."'".' | awk -v VAR=$VAR '."'".'$8<VAR'."'".' | perl extract_stats-pileup-bins_allchroms.pl -f bins_'.$ind."\n";
close(OUT);
my $qsubfile = 'run_'.$ind.'_nucdiv.txt';
`qsub -q all.q $qsubfile`;

exit;

