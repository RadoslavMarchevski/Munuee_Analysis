#!/usr/bin/perl

if($#ARGV>-1) {
  $tag=$ARGV[0];
} else {
  die("no directory tag");
}
if($#ARGV>0) {
  $ftag=$ARGV[1];
} else {
  die("no file tag");
}

if($#ARGV>1) {
  $bckup=$ARGV[2];
} else {
  $bckup="backup";
}

$file="${tag}_merge.hb";
$ffile="${tag}_${ftag}.hb";
$bfile="${tag}_${bckup}.hb";

if(-e $file) {
  print " $file $bfile\n";
  system("mv $file $bfile");
}

if(-e $ffile) {
  die("file $ffile already exists");
}

open(SH,"|csh");
print SH <<EOF;
paw <<_eop_
0
  exec ${tag}_merge
exit
_eop_
EOF
close(SH);

system("mv $file $ffile");









