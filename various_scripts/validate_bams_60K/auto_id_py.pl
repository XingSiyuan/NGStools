#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
my %opts = ();

$opts{f}='final_sample_name';
$opts{d}=200000;
getopt('fd', \%opts);
print '#!/bin/bash'."\n";
print '#$ -cwd'."\n";
print '#$ -S /bin/sh'."\n";
print '#$ -l h_vmem=50G'."\n";
my @files = `ls | grep vars-flt`;
my $teller=0;
foreach my $file (@files){
        ++$teller;
        chomp $file;
	print "echo '-----------------------'\necho 'QUERY: '\necho $file\necho 'RESULT: '\ncat $file | awk "."'".'$4~/[RYSWKM]/'."'"." | awk '".'$5>20&&$6>20&&$8<20'."' | head -$opts{d} |  perl validate_pop_all.py -p snppos10_2.txt -o out\n";
} 
exit;

  
