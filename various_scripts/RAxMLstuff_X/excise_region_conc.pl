#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Std;
my %opts = ();

#grab comandline options
getopt('fsbzed', \%opts);
my $file = $opts{f};
my $seqname = $opts{s};
my $begin = $opts{b};
my $end = $opts{e};
my $binsize = $opts{z};
my $dir=$opts{d};
my $seq = '';
my $id = '';
my $fafile = $dir.$file;
print "$dir\t$file\t$fafile\n";
open (FA, $fafile) or die $!;
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
		for (my $beginfrag = $begin; $beginfrag<$end ; $beginfrag =$beginfrag+$binsize){
		
			my $returnseq = substr($seq,($beginfrag-1),$binsize);
			my $endfrag = $beginfrag+$binsize;
			open(FRAG, ">>frag_$beginfrag-$endfrag") or die "can not open file $!\n";
			print FRAG ">$label-$beginfrag-$endfrag\n";
			format_seq($returnseq,60);
			close(FRAG);
		}
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
	        print FRAG "$intseq\n";
	        $index += $offset;
	}
}


			




