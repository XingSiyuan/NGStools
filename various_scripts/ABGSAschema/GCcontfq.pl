#!/usr/bin/perl -w
use strict;
use warnings;

my $bases=0;
my $a=0;
my $c=0;
my $g=0;
my $t=0;
my $count=0;
while (<>){
        my $line = $_;
        chomp $line;
	++$count;
	if ($count==2){
		
		my @array = split("",$line);
		foreach my $base (@array){
			if ($base eq 'G'){
				++$g;
			} 
			if ($base eq 'C'){
				++$c;
			}
			if ($base eq 'T'){  
        	                ++$t;
        	        }
        	        if ($base eq 'A'){
        	                ++$a;
        	        }
			if ($base ne 'N'){
				++$bases;
			}
			unless ($base =~ /[ACGTN]/){
				die "wrong line: $line\n";
			}
		}
	}
	if ($count==4){
		$count=0;
	}
}
#print "total number of bases: $bases\n";
#print "A: $a\nC: $c\nG: $g\nT: $t\n";
my $gccontent=($g+$c)/$bases;
print "GCcontent: $gccontent\n";
exit;



