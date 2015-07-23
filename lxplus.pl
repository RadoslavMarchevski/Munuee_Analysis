#!/usr/local/bin/perl
$edir=$ARGV[0];
$sdir=$ARGV[1];
$list=$ARGV[2];
$n1=$ARGV[3];
$n2=$ARGV[4];
$rdir=$ARGV[5];

chdir($edir);
unless(-d $sdir) {mkdir($sdir,0777);}

chdir($sdir);
unless(-e "compact") {symlink("$rdir/compact","compact");}
unless(-e "scompact.job") {symlink("$rdir/scompact.job","scompact.job");}
unless(-e "parameters.root") {symlink("$rdir/parameters.root","parameters.root");}

#unless(-d "userinc") {mkdir("userinc",0777);}
#system('cp ../userinc/histograms.txt userinc/histograms.txt');

open(IN,"$list");
open(OUT,">compact.list");

while(<IN>) {
  if($. >= $n1 && $. <= $n2) {print OUT;}
}
print OUT "*END*\n";
close(IN);
close(OUT);

$pwd=`pwd`;
chop($pwd);
#$cmd="cd $pwd; compact -l compact.list -so datasample.scmp";
#$cmd="scompact.job $pwd $sdir";
$cmd="lxplus.job $pwd ";

system("date");
system(' ssh lxplus "'.$cmd.'" &');
#system('bsub -q 1nh -o job.log -e job.log "'.$cmd.'"');
system("date");
