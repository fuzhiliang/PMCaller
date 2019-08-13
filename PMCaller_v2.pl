#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
#use Getopt::Std;
use File::Basename;
use File::Glob;
use FindBin qw($Bin); 
use Cwd qw(abs_path getcwd);  
my $abs = abs_path(getcwd());  
my $verbose ="v1.0";
###############################################################################
#命令行参数据的定义和获取，记录程序初始时间，设置参数默认值
#Get the parameter and provide the usage.
###############################################################################
my ($sam_path,$ski,$cut_off,$cut_depth,$vcf,$out,$genome);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'sam=s' => \$sam_path,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'reference=s' => \$genome,
    'cut_depth=i' => \$cut_depth,
    'ski=i' => \$ski,
    'cut_off=f' => \$cut_off,
) or die $!;
unless(defined $vcf && defined $sam_path ){&usage();exit 0;}
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $sam_path\nOutput file is $out\n";
my $samtools = "/usr/local/bin/samtools";
$genome ||= "/data/hg19/hg19_index/hg19.fa";
$cut_depth||= "50";
$ski ||= "10";
$cut_off ||= "0.05";
$out ||= "./";
`mkdir -p $out/shell` unless (-d "$out/shell");
my $vcf_name=basename($vcf);
open LOG ,">$out/shell/a.log";
print LOG "$start\n";
my @sam =glob "$sam_path/*.sam";
foreach my $sam (@sam){
	my $sam_name =basename($sam);
	$sam_name =~s/\.sam//;
	my $cmd="";
	open S ,">$out/shell/$sam_name.sh";
#	$cmd .="perl $Bin/bin/vcf_split.pl $vcf $out/$vcf_name.snp.gy $cut_off $cut_depth \n" ;
	$cmd .="perl $Bin/bin/ski.pl $out/$vcf_name.snp.gy $sam $out/$sam_name.skiout $ski \n"; 
	$cmd .="perl $Bin/bin/snptable.pl $out/$sam_name.skiout $out/$sam_name.list & \n";
	$cmd .="perl $Bin/bin/snptable_v1.pl $out/$sam_name.skiout $out/$sam_name & \n";
	print S $cmd ;
	close S;
}

unless (-f "$out/$vcf_name.snp.gy.check") {
` nohup perl $Bin/bin/vcf_split.pl $vcf $out/$vcf_name.snp.gy $cut_off $cut_depth  &`;
}
if (-f "$out/$vcf_name.snp.gy.check") {
	foreach my $sam (@sam){
	my $sam_name =basename($sam);
		$sam_name =~s/\.sam//;
		`nohup sh $out/shell/$sam_name.sh &`;
	}
}

#&runcmd($cmd);

###############################################################################
#Record the program running time!
# 输出程序运行时间
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
my $end=strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print LOG "End time id $end.\nThis compute totally consumed $duration_time s\.\n";
close LOG;
###############################################################################

sub usage {
    die(
        qq!
Usage:    perl /56T/Evan/PMCaller/script/PMCaller_v2.pl -sam /56T/Evan/PMCaller/mydata/DNA.C10T_1.sam.chr/ -vcf /56T/Evan/PMCaller/mydata/all/DNA.C10T_1.sorted.vcf -o /56T/Evan/PMCaller/mydata/all/C10T_out
Function: Template for Perl
Command:      
	-verbose	verbose
    -sam	sam_path
    -vcf	vcf_path
    -o		out 
    -reference		genome
    -cut_depth		cut_depth
    -ski			ski
    -cut_off		cut_off
Author:   Liu Yong-Xin, liuyongxin_bio\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}