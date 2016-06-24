#!/usr/bin/env perl -w
use strict;
use Cwd;
chdir getcwd;
my $HapinfoDIR="/home/shg047/oasis/monod/hapinfo/hapinfo";
my $config="/home/shg047/oasis/monod/saminfo.txt";
my $OutDIR="/home/shg047/oasis/monod/hapinfo/MergeHapinfo";
my @file=glob("*hapInfo.txt");

open F,$config;
my %sam;
while(<F>){
chomp;
my ($sample,$group)=split/\t/;
push(@{$sam{$group}},$sample);
}
close F;

foreach my $group(sort keys %sam){
        print "$group\n";
        system("touch $OutDIR/$group.hapinfo.merge");
        foreach my $label(sort @{$sam{$group}}){
                foreach my $file(@file){
                        if($file=~/$label/){
                                system("cat $file >> $OutDIR/$group.hapinfo.merge");
                        }
                }
        }
}
