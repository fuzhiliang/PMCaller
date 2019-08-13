#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use List::Util;
my $verbose ="v1.0";
#toli@s3:/56T/Evan/PMCaller/mydata$ perl /56T/Evan/PMCaller/script/bin//vcf_split.pl chr22_samtools_mpileup.vcf testv1/chr22.snp.gy2 0.05 50

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


my $bam_vcf = $ARGV[0];
my $out =$ARGV[1];
my $cut_off =$ARGV[2];
my $cut_depth = $ARGV[3];

$cut_depth ||= 30;  #把这个值调低，dp输出，后面直接有这个中间文件
$cut_off ||=0.05;

open I ,"$bam_vcf" or die $!;
open O ,">$out" or die $!;
my %hash;
while (<I>) {
	chomp;
	next if (/^#/ || /^[[<]/ );	
	my @tmp= (split /\t/,$_);
	next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
	my @alt= split ",",$tmp[4];

#	next if (@alt != 2 || length($alt[0]) !=1);
	next if (length($alt[0]) !=1);
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];
######	C       T,G,<*> 0,103,6,255,206,6,255,235,12,6  530     377,151,2,0
	next if ($dp < $cut_depth); 	
	my @ad = split "," ,$ad;
	my $max_ad=&max(@ad);

	next if ($max_ad < $dp*(1-$cut_off) && $max_ad > $dp*(0.5+0.5*$cut_off));
	my $genotype;
	if ($ad[0] >= $dp*(1-$cut_off)) {
		$genotype="$tmp[3]$tmp[3]";
	}elsif ($ad[0] <= $dp*(0.5+0.5*$cut_off) && $ad[0] >= $dp*(0.5-0.5*$cut_off) ) {
		my @a=("$tmp[3]","$alt[0]");
		@a=sort @a;
		$genotype=join "",@a;
	}elsif ($ad[0] <= $dp*$cut_off && $ad[1] >= $dp*(1-$cut_off) ) {
		$genotype="$alt[0]$alt[0]";
#	}elsif (length($alt[0])==1 && length($alt[1]) ==1 ) {
#		$alt[0]=~s/(\w+)/\u$1/;
#		$alt[1]=~s/(\w+)/\u$1/;
#		my @b =("$alt[0]","$alt[1]");
#		@b=sort @b;
#		$genotype=join "",@b;
	}else{next;}
#	my $ads=join "\t",$ad[0];
	print O "$tmp[0]\t$tmp[1]\t$tmp[3]\t$genotype\t$dp\n";
}
close I;
close O;

#`touch "$out.check" ` ;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

###############################################################################
#Scripts usage and about.
# 程序的帮助文档，良好的描述是程序重用和共享的基础，也是程序升级和更新的前提
###############################################################################
sub max {
    my $currentMaxValue = shift @_;
    foreach ( @_ ) {
        if ( $_ > $currentMaxValue ) {
            $currentMaxValue = $_;
        }
    }
    return $currentMaxValue;
}

