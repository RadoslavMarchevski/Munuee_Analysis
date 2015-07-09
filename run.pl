#!/usr/bin/perl -w

@pcs=(1
      ,2,3,4,5
      ,5,6,7,8,9,10
      ,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40
      ,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64
      ,65,66,67,68,69,70
      ,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99
      ,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124
      ,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150);
      #,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175
      #,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200);
      #,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229
      #,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250
      #,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275
      #,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300
      #,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329
      #,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350
      #,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375
      #,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400);
@pcsw=(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,11,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1);


$npc=$#pcs+1;
$npcw=0;
foreach $pc (@pcs) {
    $npcw+=$pcsw[$pc-1]
}
#$nfi=1200;

$rdir=`pwd`;
chop($rdir);
$edir="$rdir";
#$script1="$rdir/pp2nd.pl";
$script1="$rdir/pp2nd.pl";
$script2="$rdir/lxplus.pl";
$script3="$rdir/pp8nh.pl";
$script4="$rdir/pp1nw.pl";

if($ARGV[0] eq "test") {
    $listfile="$rdir/testing_list.list";
    $tag="test";
} elsif($ARGV[0] eq "ss0123_NEW_v1") {
    $listfile="$rdir/SC_SS0123_pass5_tot.list";
    $tag="ss0123_NEW_v1";
} elsif($ARGV[0] eq "kminus_kmu3") {
    $listfile="~/cmc/MC_lists/kch_mu3_kminus_100kk_5.03.2013.list";
    $tag="kminus_kmu3";
} elsif($ARGV[0] eq "kplus_kmu3") {
    $listfile="~/cmc/MC_lists/kplus_kmu3_100kk_11.03.2013.list";
    $tag="kplus_kmu3";
} elsif($ARGV[0] eq "kplus_pipie") {
    $listfile="~/cmc/MC_lists/kch_pipie_10kk_kplus.list";
    $tag="kplus_pipie";
} elsif($ARGV[0] eq "kminus_k3pi") {
    $listfile="/afs/cern.ch/user/r/rmarchev/cmc/MC_lists/kch_3pic_kminus_1kkk.list";
    $tag="kminus_k3pi";
} elsif($ARGV[0] eq "kplus_k3pi") {
    $listfile="/afs/cern.ch/user/r/rmarchev/cmc/MC_lists/kch_k3pi_kplus_1kkk_9_12_2012.list";
    $tag="kplus_k3pi";
} elsif($ARGV[0] eq "kminus_3pin") {
    $listfile="/afs/cern.ch/user/r/rmarchev/cmc/MC_lists/kch_3pin_100kk_kminus_10_12_2012.list";
    $tag="kminus_3pin";
} elsif($ARGV[0] eq "kminus_piee") {
    $listfile="/afs/cern.ch/user/r/rmarchev/cmc/MC_lists/kch_piee_1kk_kminus_05.03.2013.list";
    $tag="kminus_piee";
} elsif($ARGV[0] eq "kplus_pipimu") {
    $listfile="~/cmc/MC_lists/2013-03-06-kch_kplus_pipimu_1k.list";
    $tag="kplus_pipimu";
} elsif($ARGV[0] eq "kminus_pipie") {
    $listfile="~/cmc/MC_lists/kch_pipie_10kk_minus.list";
    $tag="kminus_pipie";
} elsif($ARGV[0] eq "kplus_pipid") {
    $listfile="~/cmc/MC_lists/kplus_pipid_100kk_07.03.2013.list";
    $tag="kplus_pipid";
} elsif($ARGV[0] eq "kplus_munuee") {
    $listfile="~/cmc/MC_lists/kplusmunuee_30kk.list";
    $tag="kplus_munuee";
} elsif($ARGV[0] eq "kminus_pipid") {
    $listfile="~/cmc/MC_lists/kminus_pipid_100kk_9_03_2013.list";
    $tag="kminus_pipid";
} elsif($ARGV[0] eq "kminus_munuee") {
    $listfile="~/cmc/MC_lists/kmin_munuee_10kk.list";
    $tag="kminus_munuee";
} elsif($ARGV[0] eq "ss1") {
    $listfile="$rdir/SC_SS1_pass5_tot.list";
    $tag="ss1ana";
} elsif($ARGV[0] eq "ss2") {
    $listfile="$rdir/SC_SS2_pass5_tot.list";
    $tag="ss2ana";
} elsif($ARGV[0] eq "ss3") {
    $listfile="$rdir/SC_SS3_pass5_tot.list";
    $tag="ss3ana";
} elsif($ARGV[0] eq "ee") {
    $listfile="$rdir/SC_ee_SS0123_pass5_tot_sorted.list";
    $tag="ee";
} elsif($ARGV[0] eq "mumu") {
    $listfile="$rdir/2mu_SC_SS0123_pass5_tot.list";
    $tag="mumu";
} elsif($ARGV[0] eq "v39_resub") {
    $listfile="$rdir/repeat_v39.list";
    $tag="v39_resub";
}elsif($ARGV[0] eq "pass5_p") {
    $listfile="$rdir/na48p_ss0123.list";
    $tag="pass5_p";
}elsif($ARGV[0] eq "pass5_c_1st") {
    $listfile="$rdir/na48c_ss0123_0k_1k.list";
    $tag="pass5_c_1st";
} elsif($ARGV[0] eq "pass5_c_2nd") {
    $listfile="$rdir/na48c_ss0123_1k_2k.list";
    $tag="pass5_c_2nd";
} elsif($ARGV[0] eq "pass5_c_3rd") {
    $listfile="$rdir/na48c_ss0123_2k_3k.list";
    $tag="pass5_c_3rd";
} elsif($ARGV[0] eq "pass5_c_4th") {
    $listfile="$rdir/na48c_ss0123_3k_4k.list";
    $tag="pass5_c_4th";
} elsif($ARGV[0] eq "pass5_c_5th") {
    $listfile="$rdir/na48c_ss0123_4k_5k.list";
    $tag="pass5_c_5th";
} elsif($ARGV[0] eq "pass5_c_6th") {
    $listfile="$rdir/na48c_ss0123_5k_6k.list";
    $tag="pass5_c_6th";
}
else {
    print "Argument \"$ARGV[0]\" not recognized.\n";
    exit;
}

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

system("cat merge.cc > $rdir/${tag}_merge.cc");
open(ROOT,">>$rdir/${tag}_merge.cc");

print ROOT "void merge() { \n";
print ROOT "   Target = TFile::Open(\"${tag}_merged.root\",\"RECREATE\");\n";
print ROOT "   FileList = new TList();\n";

#open(CMPLIST,">$rdir/${tag}_list.list");
$npcw1=0;
for($i=0;$i<$npc;$i++) {
    unless($pcsw[$i]) {next;}
    if($i==$nfijm) {$nfij++;}
    $nfijw=int($nfij*$pcsw[$i]);
    $npcw1+=$pcsw[$i];

    $n2=$n1+$nfijw-1;
    if($npcw1==$npcw) {$n2=$nfi;}
    if($n1>$n2) {next;}
    if($i<300) {
	$command="$script3 $edir ${tag}$pcs[$i] $listfile $n1 $n2 $rdir ";
    }
    elsif($i<150){
	$command="$script2 $edir ${tag}$pcs[$i] $listfile $n1 $n2 $rdir ";
    }
    elsif($i<66){
	$command="$script1 $edir ${tag}$pcs[$i] $listfile $n1 $n2 $rdir ";
    }
    else{
	$command="$script4 $edir ${tag}$pcs[$i] $listfile $n1 $n2 $rdir ";
    }


# execute
    print "$command\n";
    system("$command");
#  print "$pc ".($n2-$n1+1)."\n";
#  system('sleep 5');
    print ROOT "   FileList->Add(TFile::Open(\"$edir/${tag}$pcs[$i]/histos.root\"));\n";
#  print CMPLIST "/shift/na48d010/data13/venelin/2003/data/${tag}$pcs[$i] \n";
    $n1+=$nfijw;
}

print ROOT "   MergeRootfile(Target,FileList);\n";
print ROOT "}\n";

close(ROOT);
#close(CMPLIST);
exit;
