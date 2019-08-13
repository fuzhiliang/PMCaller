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
###############################################################################
#Get the parameter and provide the usage.
###############################################################################
#my ($out,$bam_vcf,$cut_off,$cut_depth,);
#my %opts;
#GetOptions(
#    'verbose' => \$verbose,
#    'i=s' => \$bam_vcf,
#    'o=s' => \$out, 
#    'cut_depth=i' => \$cut_depth,
#    'cut_off=f' => \$cut_off,   ##0.05
#    'libs=i' => \@libs, ## 'libs=i@' => \$libs,
#    'define=s' => \%defines ## 'define=s%' => \$defines,
#) or die $!;
#&usage unless ( exists $bam_vcf );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
#print "Input file is $bam_vcf\nOutput file is $out\n";
#print "Database file is $opts{d}\n" if defined($opts{d});
#
#$opts{h}=1 unless defined($opts{h});
#$out ||= "./";


my $bam_vcf = $ARGV[0];
my $out =$ARGV[1];
my $cut_off =$ARGV[2];
my $cut_depth = $ARGV[3];

$cut_depth ||= 50;
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
#
#	next if (@alt != 2 || length($alt[0]) !=1);
	next if (length($alt[0]) !=1);
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];
	next if ($dp < $cut_depth||$dp > ($cut_depth+100)); 
	my (@ad) = (split ",",$ad);
	my $pro=100*$ad[0]/$dp;
	my $int = sprintf "%.f", $pro;
######	C       T,G,<*> 0,103,6,255,206,6,255,235,12,6  530     377,151,2,0
	#next if ($dp != $cut_depth);
	$hash{$tmp[3]}{$int}++;
#	$hash{$tmp[3]}{$ad[0]}++;
}
foreach my $key1 (keys %hash) {
	foreach my $key2 (keys %{$hash{$key1}}) {
		print O "$key1\t$key2\t$hash{$key1}{$key2}\n";
	}
}
close I;
close O;

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

sub usage {
    die(
        qq!
Usage:    perl $Bin/vcf_split.pl  -cut_depth 50 -cut_off 0.05 -i bam_vcf -o bam.high_depth.snp.txt ;
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -cut_depth filt low depth location [50] 
          -cut_off [0.05]	
          -h header line number, default 0
Author:   Evan Fu, liuyongxin_bio\@163.com, QQ:42789409
Version:  v1.0
Update:   2018/5/16
Notes:    
\n!
    )
}
