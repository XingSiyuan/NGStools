#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Std;
my %opts = ();

#grab comandline options
getopt('fsbzed', \%opts);
my $begin = $opts{b};
my $end = $opts{e};
my $binsize = $opts{z};
my $seq = '';
my $id = '';
for (my $beginfrag = $begin; $beginfrag<$end ; $beginfrag =$beginfrag+$binsize){
	my $endfrag = $beginfrag+$binsize;
	`./convert.sh frag_$beginfrag-$endfrag >frag_$beginfrag-$endfrag.phy`;
	`./raxmlHPC-PTHREADS -T 20 -s frag_$beginfrag-$endfrag.phy -o OM001_Warthog_SSCX -m GTRGAMMA -n frag_$beginfrag-$endfrag -p 1 | grep -v 'undetermined values'`; 
	`mv RAxML_bestTree.frag_$beginfrag-$endfrag bestree.frag_$beginfrag-$endfrag`;
	`rm RAxML_*`;
	`rm frag_$beginfrag-$endfrag.phy*`;
}
