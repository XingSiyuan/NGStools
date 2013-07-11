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
	open (TMP, '>tmp.txt');
	print TMP "./convert.sh frag_$beginfrag-$endfrag >frag_$beginfrag-$endfrag.phy\n";
	print TMP "./raxmlHPC-PTHREADS -T 2 -s frag_$beginfrag-$endfrag.phy -o OM001_Warthog_SSCX -m GTRGAMMA -n frag_$beginfrag-$endfrag -p 1 -y | grep -v 'undetermined values'\n";
	my $command = 'VAR=`wc -l RAxML_info.frag_'.$beginfrag.'-'.$endfrag.' | sed '."'".'s/ \+/\t/'."'".' | cut -f1`';
	print TMP "$command\n";  
	$command = "echo bestree.frag_$beginfrag-$endfrag ".'$VAR '.'>>goodtrees1b.txt';
	print TMP "$command\n";
	$command = 'if [ $VAR -ge 6000 ]; then rm bestree.frag_'.$beginfrag.'-'.$endfrag.' ; rm frag_'.$beginfrag.'-'.$endfrag.'; fi';
	print TMP "$command\n";
	print TMP "rm RAxML_*.frag_$beginfrag-$endfrag\n";
	print TMP "rm frag_$beginfrag-$endfrag.phy*\n";
	close(TMP);	
	`sh tmp.txt`;
}
