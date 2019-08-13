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
###cd /56T/Evan/PMCaller/check_data && perl /56T/Evan/PMCaller/script/PMCaller_v1.pl -sam  chr22_head.P10T_1.sam -o test -H 40:50 -L 5:8:9:20
###############################################################################
my ($out_bam,$bam);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'o=s' => \$out_bam,
    'bam=s' => \$bam,
    ) or die $!;
unless(defined $out_bam || defined $bam ){&usage();exit 0;}

my $name=basename($bam);
$name=~s/(\.sort(ed)?)?\.bam$//;
my $head=`samtools view -H $bam`;
print $head ;
if ($head!~m/\@RG/) {
	print "GR\tID:$name\tLB:$name\tPL:illumina\tSM:$name\tPU:$name\n";
	`samtools addreplacerg -r "ID:$name"  -r "LB:$name"   -r "PL:illumina"     -r "SM:$name"   -r "PU:$name" -O BAM -o $out_bam.RG $bam `;
	$bam="$out_bam.RG";
}

#open I ,"samtools view -h --threads 30 $bam|";
open I ,"samtools view -h  $bam|";
open O ,">$out_bam.sam";
open F ,">$out_bam.filter";
print F $head ;
while (<I>) {
	chomp ;
	my $M=(split "\t",$_)[5] unless (/^@/); 
	if (defined $M && $M!~m/^([0-9]+)M$/){
		print F "$_\n";
	}else{
		print O "$_\n";
	}
}
close I ;
close O ;
my $genome;
if (-f "/data/data/hg19/hg19.fa"){
$genome="/data/data/hg19/hg19.fa" ;
}else {
$genome="/data/hg19/hg19_index/hg19.fa";
}

`samtools view -Sb $out_bam.sam  > $out_bam && rm $out_bam.sam`;
`samtools calmd  $out_bam  $genome   --output-fmt BAM  > $out_bam.MD.bam` ;
`rm $out_bam.RG  $out_bam`;
`ln -s $out_bam.MD.bam $out_bam ` ;

sub usage {
    die(
        qq!
Usage:    perl $0 -bam in.bam -o out.bam
		Function: Template for Perl 提取全局匹配的reads
Command:      
    'verbose' => \$verbose,
    'o=s' => \$out_bam,
    'bam=s' => \$bam,

Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}


