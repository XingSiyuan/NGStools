#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
my %opts = ();

use DBI;
my $database = 'pig_hapmap2';
my $server = 'localhost';
my $user = 'xxxxx';
my $pwd = 'xxxxxx';
my $pig_hapmap = DBI->connect("dbi:mysql:$database:$server", $user, $pwd);
$opts{f}='infile';
getopt('fg', \%opts);
my $posfile = $opts{f};
#my $binsize=$opts{b};

open (SNPPOS, $posfile);
my %positions10=();
while (<SNPPOS>){
        my $line = $_;
        chomp $line;
        my @elements = split("\t",$line);
        my $snp = $elements[0];
        my $chrom = $elements[1];
        my $pos = $elements[2];
		$positions10{$chrom."_".$pos}=$snp;
}
close(SNPPOS);
my %idhash=();
my $query = "select dna_name from sample_sheet8 where callrate > 0.7";
my $sql = $pig_hapmap->prepare($query);
$sql->execute();
while (my $row = $sql->fetchrow_arrayref) {
      my ($sample) = @$row;
      my $key = $sample;
      my $value=0;
      #print  "$key\t=>\t$value\n";
      $idhash{$key} = $value;

}

my $foundsnp=0;
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
	if ($snpqual > 20 && $consensusqual > 20 && $consensusbase =~ /[RYSWKM]/){
		if (exists($positions10{$chrom."_".$refpos})){
			#print "$line\n"; # prints lines when found, disabled
			my $snp_name = $positions10{$chrom."_".$refpos};
			#$snp_name = "##".$snp_name."##";
			my $snpexists=0;

			my $query = "select dna_name,fwallele1,fwallele2 from allgenotypes16 where SNP = '$snp_name'";
			my $sql = $pig_hapmap->prepare($query);
			$sql->execute();
			my @genos = ();
			while (my $row = $sql->fetchrow_arrayref) {
				my ($sample,$allele1,$allele2) = @$row;
				my $key = $sample;
				my $value=$allele1.$allele2;
				my $string = $key."##".$value;
				#print  "$key\t=>\t$value\n";
				push(@genos,$string);

			}
			foreach my $element (@genos){
				++$snpexists;
				my @int = split("##",$element);
				if (exists($idhash{$int[0]})){
					my ($al1,$al2)=split('',$int[1]);
					if ($al1 ne $al2){
						$idhash{$int[0]}=$idhash{$int[0]}+1;
					}				
					#print "\tgeno:$element";
				}
			}
			if($snpexists>0){
				++$foundsnp;
			}
			#print "\n";
		}
	
			
	}

}
$pig_hapmap->disconnect;
foreach my $key (keys(%idhash)){
	my $good = $idhash{$key};
	my $pct = $good/$foundsnp;
	if ($pct > 0.7){
		print "$key\t$pct\n";
	}
}
exit;
