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
###cd /56T/Evan/PMCaller/check_data && 
###############################################################################
my ($snptable,$errortablepre,$downsample,$matrix,$vcfexonint,$out,$genome,$tqual);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'vcfexonint=s' => \$vcfexonint,
    'o=s' => \$out, 
#    'genome=s' => \$genome,
    'snptable=s' => \$snptable,
    'errortablepre=s' => \$errortablepre,
#    'downsample=f' => \$downsample,
	'matrix' => \$matrix,
	'q=f'=> \$tqual, # the min threshold of base quality used for classify high quality reads
) or die $!;
unless(defined $vcfexonint){&usage();exit 0;}
unless(defined $snptable && defined $errortablepre) {&usage();exit 0;}
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Output file is $out\nsnptable is $snptable\nerrortablepre is $errortablepre\n";

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
open EH,"$errortablepre.EXON" or die $!;
print "$errortablepre.EXON \n";
while (<EH>) {
	chomp;
	next if (/^Genotype/);
	my ($genotype,$base,$fre)=(split "\t",$_)[0,1,4];
	$ERROR{'EXON'}{$base}{$genotype}=$fre;
}
close EH;
open EL,"$errortablepre.INT" or die $!;
while (<EL>) {
	chomp;
	next if (/^Genotype/);
	my ($genotype,$base,$fre)=(split "\t",$_)[0,1,4];
	$ERROR{'INT'}{$base}{$genotype}=$fre;
}
close EL;
#print Dumper %ERROR;
=c  准备工作  
perl /56T/Evan/PMCaller/script/bin/errortable_exon_int.pl /56T/Evan/PMCaller/NA12878/chr22/errortable/chr22.filter.gt  /56T/Evan/PMCaller/NA12878/chr22/errortable 3:5:9:15:19  chr22.filter  /data/hg19_anno/Region/refgene.nooverlap
awk '{if ($1==50 && $2==9){print $3"\t"$4"\t"$5"\t"$6"\t"$5/$6}}' /56T/Evan/PMCaller/NA12878/chr22//errortable/chr22.filter.EXON.errortable.list > /56T/Evan/PMCaller/NA12878/chr22//errortable/chr22.filter.errortable.list_9.EXON
awk '{if ($1==50 && $2==9){print $3"\t"$4"\t"$5"\t"$6"\t"$5/$6}}' /56T/Evan/PMCaller/NA12878/chr22//errortable/chr22.filter.INT.errortable.list > /56T/Evan/PMCaller/NA12878/chr22//errortable/chr22.filter.errortable.list_9.INT
perl  /56T/Evan/PMCaller/script/PMCaller_genotype/PMCaller_genotype_EXON_INT_preprocess.pl chr22.filter.vcf.m  # 输出 chr22.filter.vcf.m.EXONINT.modify
Chr     Loc     Ref     Alt     TDepth  DRef    DAlt    Location
chr22   16050159        C       T       5       4       1       INT
chr22   16050205        T       G       8       7       1       INT
=cut
open V ,"$vcfexonint" or die $!;
open O ,">$out/genotype.EXONINT" or die $!; 
#open C ,">$out/genotype.all" or die $!; 
open A ,">$out/check.list" ;
while (<V>) {
	chomp;
	next if (/TDepth/ );	
	my ($Chr,$Loc,$ref,$alt,$dp,$Dref,$DAlt,$type)= (split /\t/,$_);
	#print "$type\n" ; 
	my $p_max=0;
	my $genotype;
	$ref=uc($ref);
	$alt=uc($alt);
	my @gt;
	if($ref eq $alt){
		$gt[0]="$ref$ref";
	}else{
		$gt[0]="$ref$ref";
		$gt[1]="$ref$alt";
		unless(exists($SNP{$alt}{$gt[1]})){
			$gt[1]="$alt$ref";
		}
		$gt[2]="$alt$alt";
	}
	for (my $i=0;$i<=$#gt;$i++) {
		$ERROR{$type}{$ref}{$gt[$i]}=0 if (!exists $ERROR{$type}{$ref}{$gt[$i]});
		$ERROR{$type}{$alt}{$gt[$i]}=0 if (!exists $ERROR{$type}{$alt}{$gt[$i]});
		my $p=$SNP{$ref}{$gt[$i]}*($ERROR{$type}{$ref}{$gt[$i]}**$Dref)*($ERROR{$type}{$alt}{$gt[$i]}**$DAlt);  

		#print  "$p=$SNP{$ref}{$gt[$i]}*($ERROR{$type}{$ref}{$gt[$i]}**$Dref)*($ERROR{$type}{$alt}{$gt[$i]}**$DAlt)\n";
		if ($p >= $p_max) {
			$genotype=$gt[$i];
			$p_max=$p;
		}
	}
	print O "$Chr\t$Loc\t$ref\t$genotype\t$dp\t$Dref\t$DAlt\t$alt\n";

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
Usage:    perl $0 -vcfexonint  /56T/Evan/PMCaller/NA12878/chr22/down0.03/chr22.filter.vcf.m.EXONINT.modify -o /56T/Evan/PMCaller/NA12878/chr22/PMCallerdown9/    -snptable  /56T/Evan/PMCaller/NA12878/chr22/snptable/H40_L9/chr22.filter.exp.snp.table_list  -errortablepre /56T/Evan/PMCaller/NA12878/chr22/errortable/chr22.filter.errortable.list_9
Function: Template for Perl  
Command:      
    'verbose' => \$verbose,
    'vcfexonint=s' => \$vcfexonint,
    'o=s' => \$out, 
    'snptable=s' => \$snptable,
    'errortablepre=s' => \$errortablepre,
	'matrix' => \$matrix,  如果snp tabel 是10行4列的格式就加上该参数
Author:   Evan Fu,
Version:  v1.0
Update:   2018/08/08
Notes:    
\n!
    )
}
