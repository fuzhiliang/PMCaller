#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
my $verbose ="v1.0";

my $i =$ARGV[0];
my $o =$ARGV[1];

open I ,"$i" or die $!;
open O ,">$o" or die $!;

while (<I>) {
	chomp ;
	my @tmp =split "\t",$_;
	my @gt=split "",$tmp[3];
	next if ($gt[0] ne $tmp[2] && $gt[1] ne $tmp [2]);
	print O "$_\n" ;
}
close I;
close O;
