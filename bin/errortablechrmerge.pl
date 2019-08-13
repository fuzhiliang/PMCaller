#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
#use Getopt::Std;
use File::Basename;
use File::Glob;
use Data::Dumper;
use FindBin qw($Bin); 
use Cwd qw(abs_path getcwd);  
my $abs = abs_path(getcwd());  
my $verbose ="v1.0";
###############################################################################
#命令行参数据的定义和获取，记录程序初始时间，设置参数默认值
#Get the parameter and provide the usage.
###cd /56T/Evan/PMCaller/check_data && perl /56T/Evan/PMCaller/script/PMCaller_v1.pl -sam  chr22_head.P10T_1.sam -o test -H 40:50 -L 5:8:9:20
###############################################################################
my ($errorlist,$out,$sample);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'errorlist=s' => \$errorlist,
#	'sample=s' =>\$sample,
    'o=s' => \$out, 
) or die $!;
unless(defined $errorlist  ){&usage();exit 0;}

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $errorlist \nOutput file is $out\n";
#$out||='./';
#system "mkdir $out " unless (-d $out);

my @errorlist = glob ("$errorlist/*.errortable.list") ;
#print @errorlist ;
print "\n";
my %hash;
my %hash_total;
foreach my $errortable (@errorlist) {
	open E ,"$errortable" or die $!;
	while (<E>) {
		chomp;
		next if (/^H/);
		my @tmp=split "\t",$_;
		my $a="$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]";
		my $b="$tmp[0]\t$tmp[1]\t$tmp[2]";
		$hash{$a}+=$tmp[4];
		$hash_total{$b}+=$tmp[5];
	}
	close E;
}
print Dumper %hash;
open O , ">$out" or die $!;

foreach my $c (keys %hash) {
	print O "$c\t$hash{$c}\n";
#	print  "$c\t$hash{$c}\n";
}

sub usage {
    die(
        qq!
Usage:    perl $0 
Function: Template for Perl
Command:      
    'verbose' => \$verbose,
    'errorlist=s' => \$errorlist,
#	'sample=s' =>\$sample,
    'o=s' => \$out, 
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}