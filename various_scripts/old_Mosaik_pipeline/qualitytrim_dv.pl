#!/usr/bin/perl -w
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Std;
my %opts = ();
$opts{q}= 'i';
$opts{r}= '';
$opts{f}= '';
$opts{t}= '13';
$opts{m}= '40';
#grab comandline options
getopt('qrftm', \%opts);
my $ascivalue = get_ascivalue($opts{q});
my $minlength = $opts{m};
my $infile_forward = $opts{f};
my $infile_reverse = $opts{r};
my $treshold = $opts{t};
my $outfile_forward = $infile_forward;
$outfile_forward =~ s/\./_QT\./;
my $outfile_reverse = $infile_reverse;
$outfile_reverse =~ s/\./_QT\./;
my $outfile_reverse_orphaned = $outfile_reverse;
$outfile_reverse_orphaned =~ s/QT\./QT_orph\./;
my $outfile_forward_orphaned = $outfile_forward;
$outfile_forward_orphaned =~ s/QT\./QT_orph\./;

open(my $in1, "<:gzip",$infile_forward) or die "could not open file $outfile_forward for writing! $! \n";
open(my $in2, "<:gzip",$infile_reverse) or die "could not open file $outfile_forward for writing! $! \n";
open(my $out1, ">:gzip",$outfile_forward) or die "could not open file $outfile_forward for writing! $! \n";
open(my $out2, ">:gzip",$outfile_forward_orphaned) or die "could not open file $outfile_forward_orphaned for writing! $! \n";
open(my $out3, ">:gzip",$outfile_reverse) or die "could not open file $outfile_reverse for writing! $! \n";
open(my $out4, ">:gzip",$outfile_reverse_orphaned) or die "could not open file $outfile_reverse_orphaned for writing! $! \n";

my $forward_totalbases_original=0;
my $reverse_totalbases_original=0;
my $forward_totalreads_original=0;
my $reverse_totalreads_original=0;
	
my $forward_totalbases=0;
my $reverse_totalbases=0;
my $forward_totalbases_orphaned=0;
my $reverse_totalbases_orphaned=0;

my $forward_totalreads=0;
my $reverse_totalreads=0;
my $forward_totalreads_orphaned=0;
my $reverse_totalreads_orphaned=0;

while(<$in1>){
	my @forward_read = get_read1($_,$in1);
	#print "first read: $_";
	my @reverse_read = ();
	@reverse_read = get_read2($in2);
	
	#print "second read: $reverse_read[0]\n";
	if (get_readname($forward_read[0]) eq get_readname($reverse_read[0])){
		++$forward_totalreads_original;
		++$reverse_totalreads_original;
		$forward_totalbases_original += length $forward_read[1];
		$reverse_totalbases_original += length $reverse_read[1];

		#print "$forward_read[0]\t$reverse_read[0]\n";
		my @forward_read_QT = get_QT_read(\@forward_read,$treshold,$minlength,$ascivalue);
		my @reverse_read_QT = get_QT_read(\@reverse_read,$treshold,$minlength,$ascivalue);
		#if (length $forward_read_QT[1] >= $minlength){
		#	while ($forward_read_QT[1] =~ m/^N/ && $forward_read_QT[3] =~ m/^#/){
		#		$forward_read_QT[1] =~ s/^N//;
		#		$forward_read_QT[3] =~ s/^#//;
		#	}
		#	if (length $forward_read_QT[1] != length $forward_read_QT[3]){
		#		die "length of seq and qual no longer equal! FATAL ERROR\n";
		#	}
		#}
		#if (length $reverse_read_QT[1] >= $minlength){
                #        while ($reverse_read_QT[1] =~ m/^N/ && $reverse_read_QT[3] =~ m/^B/){
                #                $reverse_read_QT[1] =~ s/^N//;
                #                $reverse_read_QT[3] =~ s/^B//;
                #        }
                #        if (length $reverse_read_QT[1] != length $reverse_read_QT[3]){
                #                die "length of seq and qual no longer equal! FATAL ERROR\n";
                #        }
                #}

		
				
				
		if (length $forward_read_QT[1] >= $minlength && length $reverse_read_QT[1] >= $minlength){
			print $out1 "$forward_read_QT[0]\n$forward_read_QT[1]\n$forward_read_QT[2]\n$forward_read_QT[3]\n";
			print $out3 "$reverse_read_QT[0]\n$reverse_read_QT[1]\n$reverse_read_QT[2]\n$reverse_read_QT[3]\n";
			$forward_totalbases += length $forward_read_QT[1];
			++$forward_totalreads;
			$reverse_totalbases += length $reverse_read_QT[1];
			++$reverse_totalreads;
			
		}
		else {
			if (length $forward_read_QT[1] >= $minlength){
				print $out2 "$forward_read_QT[0]\n$forward_read_QT[1]\n$forward_read_QT[2]\n$forward_read_QT[3]\n";
				$forward_totalbases_orphaned += length $forward_read_QT[1];
				++$forward_totalreads_orphaned;

			}
			if (length $reverse_read_QT[1] >= $minlength){
				print $out4 "$reverse_read_QT[0]\n$reverse_read_QT[1]\n$reverse_read_QT[2]\n$reverse_read_QT[3]\n";
				$reverse_totalbases_orphaned += length $forward_read_QT[1];
				++$reverse_totalreads_orphaned;


			}
		}
		
	}
	else {
		die "the forward and reverse file are not symmetrical/sorted!!!\n";
	}

}

print "total reads forward original file: $forward_totalreads_original\n";
print "total reads reverse original file: $reverse_totalreads_original\n";
print "total bases forward original file: $forward_totalbases_original\n";
print "total bases reverse original file: $reverse_totalbases_original\n";
print "total reads forward: $forward_totalreads\n";
print "total reads reverse: $reverse_totalreads\n";	
print "total bases forward: $forward_totalbases\n";	
print "total bases reverse: $reverse_totalbases\n";	
print "total reads reverse orphaned: $forward_totalreads_orphaned\n";
print "total reads reverse orphaned: $reverse_totalreads_orphaned\n";	
print "total bases forward orphaned: $forward_totalbases_orphaned\n";	
print "total bases reverse orphaned: $reverse_totalbases_orphaned\n";	
my $recovered_forward_reads = $forward_totalreads/$forward_totalreads_original;
my $recovered_reverse_reads = $reverse_totalreads/$reverse_totalreads_original;
my $recovered_forward_bases = $forward_totalbases/$forward_totalbases_original;
my $recovered_reverse_bases = $reverse_totalbases/$reverse_totalbases_original;
print "recovered proportion reads forward: $recovered_forward_reads\n";
print "recovered proportion reads reverse: $recovered_reverse_reads\n";
print "recovered proportion bases forward: $recovered_forward_bases\n";
print "recovered proportion bases reverse: $recovered_reverse_bases\n";


exit;

sub get_ascivalue {
	my ($sanger_or_illumina) = @_;
	my $ascivalue = -9;
	if ($sanger_or_illumina eq 'i'){
		$ascivalue = 64;
	}
	elsif ($sanger_or_illumina eq 's'){
		$ascivalue = 33;
	}
	else {
		die "unsupported technology; only Illumina fq ('i') or Sanger fq ('s') are supported!\n";
	}
	return $ascivalue;
}

sub get_readname {
	my ($readname) = @_;;
	my @int = split('#',$readname);
	$readname = $int[0];
	return $readname;
}
sub get_read1 {
	my ($firstline,$fh)=@_;
	my $line1 = $firstline;
	chomp $line1;
	if ($line1 and $line1 =~ m/^@/){

		my $line2 = <$fh>;
		chomp $line2;
		#print $line2;
		my $line3 = <$fh>;
		chomp $line3;
		#print $line3;
		my $line4 = <$fh>;
		chomp $line4;
		return ($line1,$line2,$line3,$line4);
	}
}
sub get_read2 {
	my ($fh)=@_;
	my $line1 = <$fh>;
	chomp $line1;
	if ($line1 and $line1 =~ m/^@/){

		my $line2 = <$fh>;
		chomp $line2;
		#print $line2;
		my $line3 = <$fh>;
		chomp $line3;
		#print $line3;
		my $line4 = <$fh>;
		chomp $line4;
		return ($line1,$line2,$line3,$line4);
	}
}

sub get_QT_read {
	my ($readref,$treshold,$minlength,$ascivalue)=@_;
	my @read = @$readref;
	if (length $read[3]>= $minlength){
		my @illqs = split('',$read[3]);
		if ($ascivalue == 64){
			$read[3]='';
			$read[3] .= chr( ord($_) - 31) for @illqs;
		}

		my $index = $minlength -4;
		my $prev2 = ord($illqs[$index])-$ascivalue; 
		++$index;
		my $prev = ord($illqs[$index])-$ascivalue;
		++$index;
		my $here = ord($illqs[$index])-$ascivalue;
		++$index;
		while ((($prev2+$prev+$here)/3 >= $treshold) && $index < scalar @illqs){
		#while (($prev2 > $treshold || $prev > $treshold || $here > $treshold) && $index < scalar @illqs){
			$prev2 = $prev;
			$prev = $here;
			$here = ord($illqs[$index])-$ascivalue;
			#print scalar @illqs;
			#print " $index\n";
			++$index;
		}
		#print "$index\n";
		#print "untrimmed: $read[1]\n$read[3]\nlength: ". length($read[1]) ."\n";
		if ($index < length $read[1]){
			$index = $index -3;
		}
		elsif (($prev2+$prev+$here)/3<$treshold){	
			$index = $index -3;
		}
		$read[0] =~ s/-(\d+)$/:$1/;
		$read[2] =~ s/-(\d+)$/:$1/;
		$read[1] = substr($read[1],0,$index);
		$read[3] = substr($read[3],0,$index);
		#print "trimmed: $read[1]\n$read[3]\nlength: ". length($read[1]) ."\n";
	}
	return @read;
}
	
