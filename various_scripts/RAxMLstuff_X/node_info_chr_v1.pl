#!/usr/bin/perl -w
# Modified from Laurent's script
use strict;
use warnings;
#use File::Basename;
use Bio::TreeIO;
use Getopt::Std;
#my @files = @ARGV;
my %nodes;
my %clades;
my $count;
my %clades_chr;
my %count_chr;
my %opts = ();

$opts{f}='infile';
getopt('fb', \%opts);
my $files = $opts{f};
my $binsize=$opts{b};
print STDERR "$files\n";
#foreach my $files (@files) {
#	($name,$dir,$ext) = fileparse($files,'\..*');
open (TREE, "<$files");
my @alltrees=<TREE>;
close(TREE);
my $numtrees = scalar @alltrees;
my $countrees=0;
print STDERR "numtrees: $numtrees\n";
my @allbins;
while ($countrees < ($numtrees - $binsize)){	
	open (OUT,">intrees.txt");
	for (my $i=$countrees; $i<($binsize+$countrees); ++$i){
		print OUT $alltrees[$i];
	}
	++$countrees;
	close(OUT);
	#@n = split (/_/, $name);
	#$chr = $n[3];
	my $chr=$countrees;
	push(@allbins,$chr);
	print STDERR "bin: $chr\n";
	my $treeio = Bio::TreeIO->new(-format => 'newick',-file => "intrees.txt");
	while( my $tree = $treeio->next_tree ) {
		++$count;
		++$count_chr{$chr};
		my @nodes = $tree->get_nodes;
		foreach my $n (@nodes) {
		    my @desc = $n->get_all_Descendents;
			if ($n->is_Leaf) {
 	       			# degenerate clades...
    				next;
			}
    			else {
				my @lvs = grep { $_->is_Leaf } @desc;
				$clades{ join(',',sort map {$_->id} @lvs) }++;
				$clades_chr{$chr}{join(',',sort map {$_->id} @lvs) }++;
    			}
		}
	}
}
print "clade\tall";
#foreach my $c1 (sort keys %clades_chr) {
#	print "\t$c1";
#}
foreach my $c1 (@allbins){
	print "\t$c1";
}
print "\n";
my @sorteds = sort {$clades{$b} <=> $clades{$a} } keys %clades;
foreach my $sorted (@sorteds) {
	my $res = $clades{$sorted} / $count;
	print "$sorted\t$res";
	#foreach my $c (sort keys %clades_chr) {
	foreach my $c (@allbins){
		my $cf = $clades_chr{$c}{$sorted} / $count_chr{$c};
		print "\t$cf";
	}
	print "\n";
}

