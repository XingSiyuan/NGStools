#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
my %opts = ();

$opts{f}='final_sample_name';
$opts{a}='mydir';
getopt('fa', \%opts);
my $sample = $opts{f};
my $archivestring=$opts{a};
print '#!/bin/bash'."\n";
print '#$ -cwd'."\n";
print '#$ -S /bin/sh'."\n";
print '#$ -l h_vmem=50G'."\n";
my $abgsa = '/srv/mds01/shared/Sus/ABGSA/';
my $allbam = 'samtools merge merged.bam';
print "mkdir $sample\n";
my @archives=split(' ',$archivestring);
my $sampledir=$sample.'/';
my $orphdir = $sampledir.'orph';
print "mkdir $orphdir\n";
foreach my $archive (@archives){
	my $fqdir=$abgsa.$archive.'/';
	my @files = `ls $fqdir | grep .gz`;
	while (@files){
	        my $file1 = shift @files;
	        my $file2 = shift @files;
	        chomp $file1;
	        chomp $file2;
		my $trimext='_qtrim';
		print "gunzip -c $fqdir$file1 | sed 's/ /#/' | pigz >$sampledir$file1\n";;
		print "gunzip -c $fqdir$file2 | sed 's/ /#/' | pigz >$sampledir$file2\n";
		print "perl qualitytrim_dv.pl -f  $sampledir$file1 -r $sampledir$file2 -m 45 -q s >$sampledir$file1$trimext\n";
		print "rm $sampledir$file1\n";
		print "rm $sampledir$file2\n";
		print "mv $sampledir*QT_orph* $orphdir/\n";
	}
}
exit;

