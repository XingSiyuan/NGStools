#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
my %opts = ();

$opts{f}='final_sample_name';
$opts{d}=25;
getopt('fd', \%opts);
my $sample = $opts{f};
my $filterdepth = $opts{d};
print '#!/bin/bash'."\n";
print '#$ -cwd'."\n";
print '#$ -S /bin/sh'."\n";
print '#$ -l h_vmem=50G'."\n";
print '#$ -v MOSAIK_TMP=/srv/mds01/users/megen002/'."\n";
my $allbam = 'samtools merge merged.bam';
my @files = `ls *.gz`;
my $teller=0;
my $firstbam = '';
while (@files){
	++$teller;
	my $file1 = shift @files;
	my $file2 = shift @files;
	chomp $file1;
	chomp $file2;
	#print "$file1 - $file2\n";
	my $datfile = $file1;
	$datfile =~ s/(.+_\d)_1.+\.gz/$1/;
	print 'MosaikBuild -q '.$file1.' -q2 '.$file2.' -out '.$datfile.'.dat -st sanger'."\n";

	print 'MosaikAligner -in '.$datfile.'.dat -out aligned_'.$datfile.'_build10.dat -ia /srv/mds01/shared/Sus/Sscrofa_build10_2/Mosaiklibs/Ssc10_2_wUn.dat -hs 15 -mmp 0.07 -m all -mhp 10 -p 12 -act 20 -j /srv/mds01/shared/Sus/Sscrofa_build10_2/Mosaiklibs/Ssc10_2_wUn.jump15'."\n";

	print 'MosaikSort -in aligned_'.$datfile.'_build10.dat -out aligned_'.$datfile.'_build10-sorted.dat'."\n";

	print 'MosaikText -in aligned_'.$datfile.'_build10-sorted.dat -bam aligned_'.$datfile.'_build10-sorted.bam'."\n";
	if ($teller == 1){
		print 'samtools view -H aligned_'.$datfile.'_build10-sorted.bam | sed '."\'".'s/SM:unknown/SM:'.$sample.'/'."\'".' | sed '."\'".'s/PL:sanger/PL:ILLUMINA/'."\'".' >newheader.txt'."\n";
		$firstbam = 'aligned_'.$datfile.'_build10-sorted.bam';
	}
	else {
		print 'samtools view -H aligned_'.$datfile.'_build10-sorted.bam | sed '."\'".'s/SM:unknown/SM:'.$sample.'/'."\'".' | sed '."\'".'s/PL:sanger/PL:ILLUMINA/'."\'".' | grep @RG >>newheader.txt'."\n\n";
	}
	$allbam = $allbam.' aligned_'.$datfile.'_build10-sorted.bam';		


}
if ($teller == 1){
	print "mv $firstbam merged.bam\n\n";	
}
else {
	print "$allbam\n\n";
}

print 'samtools reheader newheader.txt merged.bam >'.$sample.'_rh.bam'."\n"; 
print 'samtools view -u '.$sample.'_rh.bam | samtools pileup -vcf /srv/mds01/shared/Sus/Sscrofa_build10_2/FASTA/Ssc10_2_wUn.fa - >vars-raw_'.$sample.'.txt'."\n";

print 'samtools.pl varFilter -D'.$filterdepth.' vars-raw_'.$sample.'.txt | awk '."\'".'($3=="*"&&$6>=50)||($3!="*"&&$6>=20)'."\'".' >vars-flt_'.$sample.'-final.txt'."\n";
print 'rm *.dat'."\n";
exit;

