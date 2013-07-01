#!/usr/bin/perl -w
# Hendrik-Jan Megens, 
# Animal Breeding & Genomics Centre
# Wageningen University
# hendrik-jan.megens -at- wur.nl

use strict;
use warnings;
use Getopt::Std;
my %opts = ();

$opts{f}='infile';
$opts{b}='100000';
getopt('fb', \%opts);
my $infile = $opts{f};
my $binsize=$opts{b};
my $oldchrom = 'NA';
#open (PRO, "$infile") or die "No such file! $!\n";
open (VAR, ">sv-binvars_$infile");
#open(EXON,"genefile_$infile.csv");
open(EXONINFO,">bininfofile_$infile");
my($refcount,$prev1,$prev2,$prev3,$count,$countgoodread,$bin,$countexonpos,$countexonvar,$countexondif,$totalreaddepth,$alltotalreaddepth,$allpos)=initialize($binsize);
#print FAS ">Sus_verrucosus_MT\n";
while (<>){
	my $line = $_;
	chomp $line;
	my @elements = split("\t",$line);
	my $chrom = $elements[0];
	my $refpos = $elements[1];
	my $refbase = $elements[2];
	my $consensusbase = $elements[3];
	my $consensusqual = $elements[4];
	my $snpqual = $elements[5];
	my $mappingqual = $elements[6];
	my $readdepth = $elements[7];
	my $readbases = $elements[8];
	my $basequals = $elements[9];
	if ($oldchrom eq 'NA'){
		$oldchrom = $chrom;
	}
	if ($oldchrom ne $chrom){
		print EXONINFO "$oldchrom\t$bin\t$countexonpos\t$countexondif\t$countexonvar\t$totalreaddepth\t$allpos\t$alltotalreaddepth\n";
		($refcount,$prev1,$prev2,$prev3,$count,$countgoodread,$bin,$countexonpos,$countexonvar,$countexondif,$totalreaddepth,$alltotalreaddepth,$allpos)=initialize($binsize);
		$oldchrom = $chrom;
	}
	if ($refpos > $bin){
		print EXONINFO "$chrom\t$bin\t$countexonpos\t$countexondif\t$countexonvar\t$totalreaddepth\t$allpos\t$alltotalreaddepth\n";
		$bin = $bin + $binsize;
		$countexonpos=0;
		$countexonvar=0;
		$countexondif=0;
		$totalreaddepth=0;
		$alltotalreaddepth=0;
		$allpos=0;

	}
	++$allpos;
	++$count;
	$alltotalreaddepth=$alltotalreaddepth+$readdepth;
	if (($prev1+$prev2+$prev3<100) && ($consensusqual+$snpqual)>20 && $mappingqual>20){
		if ($refpos <= $bin && $refpos >= ($bin-$binsize)){
			++$countexonpos;
			$totalreaddepth=$totalreaddepth+$readdepth;
			if ($refbase =~ m/[ACGT]/ && $consensusbase =~m/[KMSWRY]/ && $snpqual> 20 ){
					#print VAR "$line\t$bin\n";
					++$countexonvar;
			}
			elsif ($refbase =~ m/[ACGT]/ && $consensusbase =~m/[ACGT]/ && $refbase ne $consensusbase && $snpqual> 20 ){
					#print VAR "$line\t$bin\n";
					++$countexondif;
			}
			
		}
	}
	$prev3=$prev2;
	$prev2=$prev1;
	$prev1=$readdepth;
			
}
#print "$count\t$countgoodread\n";
#close(PRO);
close(VAR);
#close(EXON);
close(EXONINFO);
exit;

sub initialize {
	my ($binsize) = @_;
	my $refcount = 0;
	my $prev1=0;
	my $prev2=0;
	my $prev3=0;
	my $count=0;
	my $countgoodread=0;
	my $bin=$binsize;
	my $countexonpos=0;
	my $countexonvar=0;
	my $countexondif=0;
	my $totalreaddepth=0;
	my $alltotalreaddepth=0;
	my $allpos=0;
	return ($refcount,$prev1,$prev2,$prev3,$count,$countgoodread,$bin,$countexonpos,$countexonvar,$countexondif,$totalreaddepth,$alltotalreaddepth,$allpos)
}

