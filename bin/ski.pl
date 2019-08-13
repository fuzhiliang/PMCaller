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
my $sam= $ARGV[1];
my $out= $ARGV[2];
my $ski= $ARGV[3];
$ski ||= 10;


my %hash;
open I ,"$in" or die $!;
while (<I>) {
	chomp;
	my @tmp= split "\t", $_;
	my $genotype=$tmp[3];
	my $ref= $tmp[2];
	my $tmp_t= join "\t",($tmp[1],$tmp[2]);
	$hash{$tmp_t}="$ref\t$genotype";
}

close I;
open B ,"$sam" or die $!;
my $sam_base="";
my %number;
my %hash2;
while (<B>) {
	chomp;
	next if (/^@/);
	my @tmp= split "\t", $_;
	if ($tmp[5]=~m/([0-9]{2,3})M/){
		my $length=$1;
		my @seq_base=split "", $tmp[9];
		for (my $i=0 ;$i <$#seq_base ;$i++) {
			my $location = $tmp[3]+$i ;
			my $loc_t= "$tmp[2]\t$location";
			if (exists $hash{$loc_t}) {
				$number{$loc_t}++;
				next if ($number{$loc_t} >$ski);
				my $value= "$loc_t\t$hash{$loc_t}";
				if (exists $hash2{$value}) {
					$sam_base =$hash2{$value};
					$sam_base .=$seq_base[$i];
					$hash2{$value}=$sam_base;
				}else{
					$hash2{$value}=$seq_base[$i]; 
				}
			}
		
		}	
	}
}
open O ,">$out" or die $!;
foreach my $keys (keys %hash2) {
	print O "$keys\t$hash2{$keys}\n";
}
close I;
close B;
close O;


###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
