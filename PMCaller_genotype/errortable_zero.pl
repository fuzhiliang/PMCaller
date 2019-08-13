#!/usr/bin/perl -w
use strict;
use warnings;

open I ,"$ARGV[0]" or die $!;
open O ,">$ARGV[1]" or die $!;
my $head=<I>;
print O $head ;
while (<I>) {
	chomp;
	my @tmp=split "\t",$_;
	my ($g1,$g2)=split "",$tmp[2];
	if (($g1 ne $g2) && ($tmp[3] ne $g1) && ($tmp[3] ne $g2)){
		$tmp[4]=0;
	}
	my $tmp=join "\t",@tmp;
	print O "$tmp\n";
}
close I;
close O;
