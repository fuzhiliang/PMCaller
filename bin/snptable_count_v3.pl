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
my @base=("A","C","G","T");
my @genotype=("AA","AT","AC","AG","TT","CT","GT","CC","CG","GG");
open S ,">$out.snp.table" or die $!;  #snptable
my $ATCG=join "\t",@base;
print S "Snp_table\t$ATCG\n";
foreach my $genotype ( @genotype) {
	my $p_list;
	foreach my $ref (@base) {
		my $p;
		if (exists $hash_snp{$ref}{$genotype}) {
$p=$hash_snp{$ref}{$genotype};
#$p=$hash_snp{$ref}{$genotype}/$count_total{$ref};
		}else{
			$p=0;
		}
		$p_list .="\t$p";
	}
	print S "$genotype$p_list\n";

}
close S;

open E ,">$out.error.table" or die $!;  #errortable
my $genotype_list=join "\t",@genotype;
print E "Error_table\t$genotype_list\n";
foreach my $base ( @base) {
		my $p_list;
	foreach my $genotype (@genotype) {
		my $p;
		if (exists $hash_error{$genotype}{$base}) {
			$p=$hash_error{$genotype}{$base};
#			$p=$hash_error{$genotype}{$base}/$count_error{$genotype};
		}else{
			$p=0;
		}
		$p_list .="\t$p"; 	
	}
	print E "$base$p_list\n";
}
close E;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";


# perl  /56T/Evan/PMCaller/script//bin/snptable_v3.pl /56T/Evan/PMCaller/mydata/all/DNA.C10T_1.ski /56T/Evan/PMCaller/mydata/all/DNA.C10T_1_freq