#!/usr/local/bin/perl

@pcs= (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
@pcsw=(1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1);
$npc=$#pcs+1;
$npcw=0;
foreach $pc (@pcs) {
    $npcw+=$pcsw[$pc-1]
}
#$nfi=1200;

$rdir=`pwd`;
chop($rdir);
$edir="$rdir";
$listfile="$rdir/dilution-pass2.castor.short";
$script2="$rdir/pp8nh.pl";
$script1="$rdir/pp1nd.pl";
$tag="sc";

$nfi=`wc -l $listfile`;
$nfi--;
#overrule nfi
#$nfi=12;

if($nfi<$npc) {$npc=$nfi;}
$nfij=int($nfi/$npcw);
$nfijm=int($npcw-($nfi-$npcw*$nfij));
$n1=1;

#$pwd=`pwd`;
#chop($pwd);

open(KUMAC,">$rdir/${tag}_merge.kumac");
print KUMAC "nt/hmerge $edir/${tag}_merge.hb _\n";

$npcw1=0;
for($i=0;$i<$npc;$i++) {
  unless($pcsw[$i]) {next;}
  if($i==$nfijm) {$nfij++;}
  $nfijw=int($nfij*$pcsw[$i]);
  $npcw1+=$pcsw[$i];  

  $n2=$n1+$nfijw-1;
  if($npcw1==$npcw) {$n2=$nfi;}
  if($n1>$n2) {next;}
  if($i<8) {
    $command="$script1 $edir ${tag}$pcs[$i] $listfile $n1 $n2 $rdir ";
  } else {
    $command="$script2 $edir ${tag}$pcs[$i] $listfile $n1 $n2 $rdir ";
  }


# execute
  print "$command\n";   
  system("$command");
#  print "$pc ".($n2-$n1+1)."\n";

  print KUMAC "      $edir/${tag}$pcs[$i]/super.hb _\n";

  $n1+=$nfijw;
}

print KUMAC "\n";
print KUMAC "h/fil 1 ${tag}_merge.hb\n";

close(KUMAC);

exit;



