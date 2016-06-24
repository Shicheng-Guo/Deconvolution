#!/usr/bin/env perl
use strict;
my @f=qw /0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0/;
foreach my $f(@f){
my $f1=int(1386115*$f);
my $f2=int(24850020*(1-$f));
my $cmd1="perl ../randomSampleFromHaploInfo.pl $f1 6-T-Merge.hapinfo.txt > 6T";
my $cmd2="perl ../randomSampleFromHaploInfo.pl $f2 WB.hapinfo.merge > WB";
my $cmd3="cat WB 6T > Wcmix.$f.hapinfo";
print "$cmd1\n$cmd2\n$cmd3\n";
}

