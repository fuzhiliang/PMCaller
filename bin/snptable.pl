#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
my $verbose ="v1.0";

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
 #testv1/chr22.snp.gy2
my $in= $ARGV[0];
my $out= $ARGV[1];

my %hash_snp;
my %count_total;
my %hash_error;
my %count_error;
open I ,"$in" or die $!;
while (<I>) {
	chomp;
	my @tmp= split "\t", $_;
	my $length=length $tmp[4];
	my $count ="$tmp[2]"x$length;
	next if ($count eq $tmp[4]);
	#snp table
	$hash_snp{$tmp[2]}{$tmp[3]}++;
	$count_total{$tmp[2]}++;
	#error table
	my @base = split "", $tmp[4];
	foreach my $base (@base) {
		$hash_error{$tmp[3]}{$base}++;
		$count_error{$tmp[3]}++;
	}
}

close I;

open S ,">$out.snp.table" or die $!;  #snptable
print S "Ref\tGenotype\tFrequency\ttotal\t%F\n";
foreach my $ref (sort {$a cmp $b} keys %hash_snp) {
	foreach my $genotype ( keys %{$hash_snp{$ref}}) {
		my $p=$hash_snp{$ref}{$genotype}/$count_total{$ref};
		print S "$ref\t$genotype\t$hash_snp{$ref}{$genotype}\t$count_total{$ref}\t$p\n";
	}
}
close S;

open E ,">$out.error.table" or die $!;  #errortable
print E "Genotype\tBase\tFrequency\ttotal\t%F\n";
foreach my $genotype (sort {$a cmp $b} keys %hash_error) {
	foreach my $base ( keys %{$hash_error{$genotype}}) {
		my $p=$hash_error{$genotype}{$base}/$count_error{$genotype};
		print E "$genotype\t$base\t$hash_error{$genotype}{$base}\t$count_error{$genotype}\t$p\n";
	}
}
close E;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
