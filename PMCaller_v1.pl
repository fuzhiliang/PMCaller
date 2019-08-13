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
my ($sam,$bam,$L_list,$cut_off,$H_list,$vcf,$out,$genome,$threads);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'sam=s' => \$sam,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'reference=s' => \$genome,
    'H=s' => \$H_list,
    'L=s' => \$L_list,
    'cut_off=f' => \$cut_off,
	'threads=i' => \$threads,
) or die $!;
unless(defined $sam || defined $bam ){&usage();exit 0;}

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $sam\nOutput file is $out\n";

my $samtools = "/usr/local/bin/samtools";
$out ||= "./";
system "mkdir  -p $out" unless (-d $out);
$genome ||= "/data/hg19/hg19_index/hg19.fa";
my $genome_name=basename($genome);
unless (-f "$genome.fai") {
	system "mkdir  $out/genome/ && ln -s $genome $out/genome/";
	`$samtools faidx $out/genome/$genome_name `;
	$genome="$out/genome/$genome_name";
}
$H_list ||= ("20");
$L_list ||= ("10");
$cut_off ||= "0.05";
$threads ||=30;
`mkdir $out/shell` unless (-d "$out/shell");
my $sample;
if (-f $bam){
	$sample=basename($bam);
	$sample=~s/(\.sorted)?\.bam$//;
}elsif(-f $sam){
	$sample=basename($sam);
	$sample=~s/\.sam$//;
	print "$samtools view  --threads $threads -Sb $sam |$samtools sort --threads $threads > $out/$sample.sorted.bam \n";
	system "$samtools view  --threads $threads -Sb $sam |$samtools sort --threads $threads > $out/$sample.sorted.bam ";
	$bam = "$out/$sample.sorted.bam ";
}else{
	die "Input file error. please input format file of bam or sam !" ;
}

open LOG ,">$out/shell/$sample.log";
print LOG "$start\n";

my $cmd="";

unless (-f $vcf) {
	open S ,">$out/shell/$sample.mileup.sh";
	$cmd.="$samtools  mpileup -t DP,AD -uvf $genome $bam > $out/$sample.vcf \n ";
	$vcf=" $out/$sample.vcf";
	print S $cmd;
	close S;
	system "sh $out/shell/$sample.mileup.sh  ";
}

#samtools mpileup -t DP,AD -uvf /data/hg19/hg19_index/hg19.fa /56T/Evan/PMCaller/mydata/${sample}_1.sorted.bam > /56T/Evan/PMCaller/mydata/${sample}_1.sorted.vcf
#多个H 
$cmd="";
open G ,">$out/shell/$sample.error_table.sh";
$cmd .="perl $Bin/bin/vcf_split.pl $vcf $out/$sample.snp.gy $cut_off 20 && \n" unless (-s "$out/$sample.snp.gy");
#$cmd .="perl $Bin/bin/gy_sam_error.pl $out/$sample.snp.gy  $bam  $out/$sample.errortable.list $L_list $H_list  \n" ;
$cmd .="perl $Bin/bin/gy_sam_error_v1.pl $out/$sample.snp.gy  $bam  $out/$sample.errortable.list $L_list $H_list  \n" ;
print G $cmd;  
close G;
system "sh $out/shell/$sample.error_table.sh ";


#open E ">$out/shell/$sample.error_table.pic.sh";
#$cmd='';b
#$cmd.="perl $Bin/bin/error.pic_v2.pl ";
#perl /56T/Evan/PMCaller/script/bin/error.pic.pl  /56T/Evan/PMCaller/mydata/WES_error_table/${sample}/  /56T/Evan/PMCaller/mydata/WES_error_table/${sample}/pic.$sample.$h



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
Usage:    perl $0 
Function: Template for Perl
Command:      
    'verbose' => \$verbose,
    'sam=s' => \$sam,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'reference=s' => \$genome,
    'H=s' => \$H_list,
    'L=s' => \$L_list,
    'cut_off=f' => \$cut_off,
	'threads=i' => \$threads,
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}