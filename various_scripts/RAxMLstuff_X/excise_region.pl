#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Std;
my %opts = ();

#grab comandline options
getopt('fsbe', \%opts);
my $file = $opts{f};
my $seqname = $opts{s};
my $begin = $opts{b};
my $end = $opts{e};

my $seq = '';
my $id = '';

open (FA, $file) or die $!;
while (<FA>){
	my $line = $_;
	chomp $line;
	if ($line =~ m/^>/){
		if ($id){
			do_something($id,uc $seq,$seqname,$begin,$end,$file);
		}
		
		$id = $file;
		#$id = $line;
		$id =~ s/\.fa//;
		#$id =~ s/^>//;
		$seq = '';
	}
	else {
		$seq = $seq . $line;
	}
}
close(FA);
do_something($id,uc $seq,$seqname,$begin,$end,$file);
exit;

sub do_something {
	my($id,$seq,$seqname,$begin,$end,$file) = @_;
	if ($seqname =~ m/$id/){
		my @int = split(/\./,$file);
		my $label = $int[0];
		my $length = $end-$begin+1;
		my $returnseq = substr($seq,($begin-1),$length);
		print ">$label-$begin-$end\n";
		format_seq($returnseq,60);
	}
}
sub format_seq {
	my($seq,$offset)=@_;
	my $index = 0;
	while ($index < length $seq){
	        my $intseq;
	        if ($index+$offset < length $seq){
	                $intseq = substr($seq,$index,$offset);
	        }
	        else {
	                my $remaining = (length $seq) - $index;
	                $intseq = substr($seq,$index,$remaining);
	        }
	        print "$intseq\n";
	        $index += $offset;
	}
}


			




