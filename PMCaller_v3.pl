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
my ($sam,$bam,$L_list,$cut_off,$H_list,$vcf,$out,$genome,$threads,$cut_H,$H_min,$H_max,$quality);
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
    'cut_H=f' => \$cut_H,
	'Hmin=i' => \$H_min,
	'Hmax=i' => \$H_max,
	'threads=i' => \$threads,
	'quality=i' => \$quality,
) or die $!;
unless(defined $sam || defined $bam ){&usage();exit 0;}

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $sam\nOutput file is $out\n";
open LOG ,">$out/shell/a.log";
print LOG "$start\n";
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
$H_list ||= ("30:50:70:100");
$L_list ||= ("3:5:7:9:11:13:15:17:19:21");
$cut_off ||= "0.05";
$threads ||=30;
$cut_H=70;
$H_min=30;
$H_max=600;
`mkdir $out/shell` unless (-d "$out/shell");
my $sample;
if (-f $bam){
	$sample=basename($bam);
	$sample=~s/(\.sorted)?\.bam$//;
	unless ($bam=~/sort/) {
		system "$samtools sort  --threads $threads  $bam >  $out/$sample.sorted.bam";
		$bam="$out/$sample.sorted.bam";
	}
}elsif(-f $sam){
	$sample=basename($sam);
	$sample=~s/\.sam$//;
	print "$samtools view  --threads $threads -Sb $sam |$samtools sort --threads $threads > $out/$sample.sorted.bam \n";
	system "$samtools view  --threads $threads -Sb $sam |$samtools sort --threads $threads > $out/$sample.sorted.bam ";
	$bam = "$out/$sample.sorted.bam ";
}else{
	die "Input file error. please input format file of bam or sam !" ;
}

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
my $vcf_name=basename($vcf);
$vcf_name=~s/\.vcf$//;
my @L_list=(split ":",$L_list);
open G ,">$out/shell/$sample.error_table.sh";
$cmd .= "perl $Bin/bin/vcf_expansion_H60_120.pl   -vcf_list  $vcf -o  $out/snptable -Hmin $H_min -Hmax $H_max \n ";
foreach my $L (@L_list) {
$cmd .= "perl $Bin/bin/ski_Expansionvcf_bam_H70_L.pl -vcf $out/snptable/$vcf_name.exp.vcf  -o $out/snptable/H${cut_H}_L$L -cut_H $cut_H  -L $L &\n";
}
##error table 
$cmd .="perl $Bin/bin/vcf_split.pl $vcf $out/$sample.snp.gy $cut_off 20 \n" ; # unless (-s "$out/$sample.snp.gy"); #通过4个cutoff确定 的高可信度的位点
if (defined $quality) {
	$cmd .="perl $Bin/bin/gy_sam_error_quality_v1.pl -sam $bam -H_list $H_list -L_list $L_list -o $out -gy $out/$sample.snp.gy -quality $quality \n";
}else{
	$cmd .="perl $Bin/bin/gy_sam_error_v2.1.pl $out/$sample.snp.gy  $bam  $out/$sample.errortable.list $L_list $H_list  \n" ;
}
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
Usage:    perl $0 -bam  chr22.P10T_1.bam -o chr22 &
Function: Template for Perl  通过bam/sam 或bam和vcf 文件做snp table和error table
Command:      
	'verbose' => \$verbose,
	'sam=s' => \$sam,  
	'bam=s' => \$bam,			 #必须的
	'vcf=s' => \$vcf,			 #samtools mpileup 生成的结果，有就给
	'o=s' => \$out,  
	'reference=s' => \$genome,   # s1 服务器务必加上这个参数
	'H=s' => \$H_list,           #改少一点，会跑的快一些  ["30:50:70:100"] 
	'L=s' => \$L_list,          #削薄到多少x，“：“号分割   ["3:5:7:9:11:13:15:17:19:21"]
    'cut_off=f' => \$cut_off,   #error table 的cut off  [0.05]
	'cut_H=f' => \$cut_H,       #snp table  多少x以上的数据做削薄  [70]
	'Hmin=i' => \$H_min,        #取小一点，[30]
	'Hmax=i' => \$H_max,        #取大一点，[600] ,如果Hmax=Hmin 则是取单一深度的
	'threads=i' => \$threads,
 	'quality=i' => \$quality,   #分质量值，会生成两个error table   [20]
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}