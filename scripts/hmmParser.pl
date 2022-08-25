#!/usr/bin/perl

use warnings;
use strict;
use English;

my $INPUT=$ARGV[0];

my %HammerOutput;

open(IN, $INPUT) or die "can't open $INPUT";
while (<IN>){
	if ($_ !~ '#'){
        chomp;
        my @HMMrow = (split ' ', $_);
        my $GeneName = $HMMrow[0];
        # add all domain in one row, comma separated
        if(exists($HammerOutput{"$GeneName"})) {
            $HammerOutput{"$GeneName"} = $HammerOutput{"$GeneName"} . ', ' . $HMMrow[2];
        }else{
            $HammerOutput{"$GeneName"} = $HMMrow[2];
            }
        }
    }
close(IN);
   
   
while ( my ($k,$v) = each %HammerOutput ) {
    print "$k\t$v\n";
}
