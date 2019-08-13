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
my ($snptable,$errortablepre,$downsample,$matrix,$vcfmerge,$out,$genome,$tqual);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'vcfmerge=s' => \$vcfmerge,
    'o=s' => \$out, 
    'genome=s' => \$genome,
    'snptable=s' => \$snptable,
    'errortablepre=s' => \$errortablepre,
    'downsample=f' => \$downsample,
	'matrix' => \$matrix,
	'q=f'=> \$tqual, # the min threshold of base quality used for classify high quality reads
) or die $!;
unless(defined $vcfmerge){&usage();exit 0;}
unless(defined $snptable && defined $errortablepre) {&usage();exit 0;}
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Output file is $out\nsnptable is $snptable\nerrortablepre is $errortablepre\n";

`mkdir -p $out ` unless (-d $out);
#toli@s3:/56T/Evan/PMCaller/mydata/all$ head DNA.C10T_11_list.snp.table
#toli@s3:/56T/Evan/PMCaller/mydata/all$ head ../WES_error_table/DNA.C10T/DNA.C10T.30.eout.11.error.table.list
open N  ,"$snptable" or die $!;
my %SNP;
my %ERROR;

if($matrix){
	my %count;
	open S ,"$snptable" ;
	while (<S>) {
		chomp ;
		next if (/Snp/i);
		my @tmp=split "\t" , $_;
		for (my $i=1; $i<=4 ;$i++) {
			$count{$i}+=$tmp[$i];
		}
	}
	close S;
	while (<N>) {
		chomp;
		next if (/Snp/);
		my @tmp=split "\t", $_;
		my $genotype = $tmp[0];
		my @ref=("snptable","A","C","G","T");
		for (my $i=1;$i<=4 ;$i++) {
			my $fre=$tmp[$i]/$count{$i};
			$SNP{$ref[$i]}{$genotype}=$fre;
		}
	}
}else{
	while (<N>) {
		chomp;
		next if (/^Ref/);
		my ($ref,$genotype,$fre)=(split "\t",$_)[0,1,4];
		$SNP{$ref}{$genotype}=$fre;
	}
}
close N;
open EH,"$errortablepre.HighQual" or die $!;
while (<EH>) {
	chomp;
	next if (/^Genotype/);
	my ($genotype,$base,$fre)=(split "\t",$_)[0,1,4];
	$ERROR{'H'}{$base}{$genotype}=$fre;
}
close EH;
open EL,"$errortablepre.LowQual" or die $!;
while (<EL>) {
	chomp;
	next if (/^Genotype/);
	my ($genotype,$base,$fre)=(split "\t",$_)[0,1,4];
	$ERROR{'L'}{$base}{$genotype}=$fre;
}
close EL;

open V ,"$vcfmerge" or die $!;
open O ,">$out/genotype.HighQual_LowQual.qual$tqual" or die $!; 
#open C ,">$out/genotype.all" or die $!; 
open A ,">$out/check.list" ;
while (<V>) {
	chomp;
	next if (/TDepth/ );	
	my ($Chr,$Loc,$ref,$alt,$dp,$HQual_Dref,$HQual_DAlt,$LQual_Dref,$LQual_DAlt)= (split /\t/,$_);
	my $p_max=0;
	my $genotype;
	my $ref_UP=uc($ref);
	my $alt_UP=uc($alt);
	my @gt;
	if($ref_UP eq $alt_UP){
		$gt[0]="$ref_UP$ref_UP";
	}else{
		$gt[0]="$ref_UP$ref_UP";
		$gt[1]="$ref_UP$alt_UP";
		unless(exists($SNP{$alt_UP}{$gt[1]})){
			$gt[1]="$alt_UP$ref_UP";
		}
		$gt[2]="$alt_UP$alt_UP";
	}
	for (my $i=0;$i<=$#gt;$i++) {
		my $p=$SNP{$ref_UP}{$gt[$i]}*(($ERROR{'H'}{$ref_UP}{$gt[$i]}**$HQual_Dref)*($ERROR{'H'}{$alt_UP}{$gt[$i]}**$HQual_DAlt)*($ERROR{'L'}{$ref_UP}{$gt[$i]}**$LQual_Dref)*($ERROR{'L'}{$alt_UP}{$gt[$i]}**$LQual_DAlt));  #
		print  "$_\n$p=$SNP{$ref_UP}{$gt[$i]}*(($ERROR{'H'}{$ref_UP}{$gt[$i]}**$HQual_Dref)*($ERROR{'H'}{$alt_UP}{$gt[$i]}**$HQual_DAlt)*($ERROR{'L'}{$ref_UP}{$gt[$i]}**$LQual_Dref)*($ERROR{'L'}{$alt_UP}{$gt[$i]}**$LQual_DAlt)\n";
#		my $p=$SNP{$alt[0]}{$gt}*($ERROR{$alt[0]}{$gt}**$ad[0])*($ERROR{$alt[1]}{$gt}**$ad[1]);
		if ($p >= $p_max) {
			$genotype=$gt[$i];
			$p_max=$p;
		}
		print "dddddddd$genotype\n\n";

	}
	print "Ture genotype $genotype \n";
	print O "$Chr\t$Loc\t$ref\t$genotype\t$dp\t$HQual_Dref\t$HQual_DAlt\t$LQual_Dref\t$LQual_DAlt\t$alt_UP\n";
	
}
close V;
close O;
#close C;
close A;


###############################################################################
#Record the program running time!
# 输出程序运行时间
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
my $end=strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "End time id $end.\nThis compute totally consumed $duration_time s\.\n";
#close LOG;
###############################################################################

sub runcmd{
	my $name=shift @_;
	my $cmd=shift @_;
	`mkdir "$out/shell/"` unless (-d "$out/shell");
	open S ,">$out/shell/$name.sh" or die $!;
	print S "$cmd ";
	system "sh $out/shell/$name.sh ";
	close S;
}

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
Usage:    perl $0 -vcf /8T1/Evan/PMCaller/WES_filter/P#180/DNA.C10T_1/PMCaller_vs_gatk_snpwindows/down3/test/DNA.C10T_1.filter.sorted.DP0.06.MD.HighQual.LowQual.qual20.dubious.merge -o /8T1/Evan/PMCaller/WES_filter/P#180/DNA.C10T_1/PMCaller_vs_gatk_snpwindows/down3/test/ -snptable  /8T1/Evan/PMCaller/WES_filter/P#180/DNA.C10T_1/snptable/exp_ski/H60_ski40/snptable_H20_L3/DNA.C10T_1.filter.exp.snptable.list.snp.table_list  -errortablepre /8T1/Evan/PMCaller/WES_filter/P#180/DNA.C10T_1/PMCaller_vs_gatk_snpwindows/down3/test/DNA.C10T_1.filter.MD.errortable.list_L3 
Function: Template for Perl  
Command:      
    'verbose' => \$verbose,
    'vcf=s' => \$vcfmerge,
    'o=s' => \$out, 
    'genome=s' => \$genome,
    'snptable=s' => \$snptable,
    'errortablepre=s' => \$errortablepre,
    'downsample=f' => \$downsample,
	'matrix' => \$matrix,  如果snp tabel 是10行4列的格式就加上该参数
	'q=f' => \$tqual, # the min threshold of base quality used for classify high quality reads
Author:   Tania Wang,
Version:  v1.0
Update:   2018/07/16
Notes:    
\n!
    )
}
