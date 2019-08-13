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
my $col=$ARGV[2];  #一长串base在第几列 
$col ||= 5;
$col = $col-1;
my %hash_snp;
my %count_total;
my %hash_error;
my %count_error;
open I ,"$in" or die $!;
my $length;
while (<I>) {
	chomp;
	my @tmp= split "\t", $_;
	$length=length $tmp[$col];
	my $count ="$tmp[2]"x$length;
	next if ($count eq $tmp[$col]);
	#snp table
	$hash_snp{$tmp[2]}{$tmp[3]}++;
	$count_total{$tmp[2]}++;
	#error table
	my @base = split "", $tmp[$col];
	foreach my $base (@base) {
		$hash_error{$tmp[3]}{$base}++;
		$count_error{$tmp[3]}++;
	}
}

close I;
my @base=("A","C","G","T");
my @genotype=("AA","AT","AC","AG","TT","CT","GT","CC","CG","GG");

open E ,">$out.$length.error.table" or die $!;  #errortable
open C ,">$out.$length.error.table.count" or die $!;  #errortable

my $genotype_list=join "\t",@genotype;
print E "Error_table\t$genotype_list\n";
print C "Error_table\t$genotype_list\n";

foreach my $base ( @base) {
		my $p_list;
		my $count_list;
	foreach my $genotype (@genotype) {
		my $p;
		my $count;
		if (exists $hash_error{$genotype}{$base}) {
			$p=$hash_error{$genotype}{$base}/$count_error{$genotype};
			$count=$hash_error{$genotype}{$base};
		}else{
			$p=0;
			$count=0;
		}
		$p_list .="\t$p"; 
		$count_list .="\t$count";
	}
	print E "$base$p_list\n";
	print C "$base$count_list\n";
}
close C;
close E;

open F ,">$out.$length.error.table.list" or die $!;  #errortable
print F "Genotype\tBase\tFrequency\ttotal\t%F\n";
foreach my $genotype (sort {$a cmp $b} keys %hash_error) {
	foreach my $base ( keys %{$hash_error{$genotype}}) {
		my $p=$hash_error{$genotype}{$base}/$count_error{$genotype};
		print F "$genotype\t$base\t$hash_error{$genotype}{$base}\t$count_error{$genotype}\t$p\n";
	}
}
close F;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";



