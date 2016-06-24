
#!/usr/bin/env perl
use strict;
my @f=qw /0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0/;
use Cwd;
my $dir=getcwd;
foreach my $f(@f){
    open(JOB_FILE, ">Wcmix.$f.hapinfo.job") || die("Error in opening file Wcmix.$f.hapinfo\n");
    print JOB_FILE "#!/bin/csh\n";
    print JOB_FILE "#PBS -q hotel\n";
    print JOB_FILE "#PBS -l nodes=1:ppn=1\n";
    print JOB_FILE "#PBS -l walltime=8:00:00\n";
    print JOB_FILE "#PBS -o Wcmix.$f.log\n";
    print JOB_FILE "#PBS -e Wcmix.$f.err\n";
    print JOB_FILE "#PBS -V\n";
    print JOB_FILE "#PBS -M shihcheng.guo\@gmail.com \n";
    print JOB_FILE "#PBS -m abe\n";
    print JOB_FILE "#PBS -A k4zhang-group\n";
    print JOB_FILE "cd $dir\n";
    my $f1=int(1386115*$f);
    my $f2=int(24850020*(1-$f));
    my $cmd1="perl ../randomSampleFromHaploInfo.pl $f1 ../6-T-Merge.hapinfo.txt > 6T.$f.txt";
    my $cmd2="perl ../randomSampleFromHaploInfo.pl $f2 ../WB.hapinfo.merge > WB.$f.txt";
    my $cmd3="cat WB.$f.txt 6T.$f.txt > Wcmix.$f.hapinfo";
    unlink "WB.$f.txt";
    unlink "6T.$f.txt";
    print JOB_FILE "$cmd1\n$cmd2\n$cmd3\n";
    system("qsub Wcmix.$f.hapinfo.job");
}
