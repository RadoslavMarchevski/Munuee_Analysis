#!/usr/local/bin/perl

$t = $argv[0];

    mkdir("$t_run",0777);
    for($i=1;$i<=16;$i++) {
	system("mv $t$i $t_run");
    }




