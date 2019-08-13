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
###my $vcf=$ARGV[0];
#my $tag=$ARGV[1];

#perl /56T/Evan/PMCaller/script/PMCaller_genotype/PMCaller_genotype_High_Low_preprocess.pl  $root/$sample/PMCaller_vs_gatk_snpwindows/down$i/$sample.filter.sorted.DP$a.MD.HighQual.qual20.vcf
###############################################################################
my ($bam,$snptable,$errortable,$downsample,$matrix,$vcf,$out,$genome);
my %opts;
GetOptions(
    'verbose' => \$verbose,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'genome=s' => \$genome,
    'snptable=s' => \$snptable,
    'errortable=s' => \$errortable,
    'downsample=f' => \$downsample,
	'matrix' => \$matrix,
) or die $!;
unless(defined $vcf || defined $bam ){&usage();exit 0;}
unless(defined $snptable && defined $errortable) {&usage();exit 0;}
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $inputfile;
$inputfile=$bam if ($bam);
$inputfile=$vcf if ($vcf);
print "Input file is $inputfile\nOutput file is $out\nsnptable is $snptable\nerrortable is $errortable\n";
my $picard;
$picard="/soft/picard-tools-2.3.0/picard.jar" if (-e "/soft/picard-tools-2.3.0/picard.jar"); 
$picard="/soft/picard-tools-2.4.1/picard.jar" if (-e "/soft/picard-tools-2.4.1/picard.jar"); 

my $samtools ="/usr/local/bin/samtools";
$genome||="/data/hg19/hg19_index/hg19.fa" if (-e "/data/hg19/hg19_index/hg19.fa");
$genome||="/data/data/hg19/hg19.fa" if (-e "/data/data/hg19/hg19.fa");


my $threads||=30;
`mkdir -p $out ` unless (-d $out);
$out = abs_path($out);
$out||='./';
die "defined $downsample and !defined $bam" if (defined $downsample && !defined $bam );
warn "ignore bam file :$bam and downsample: $downsample because exists vcf file:$vcf .\n" if ((defined $downsample || defined $bam) && defined $vcf) ;
my $cmd;
if (defined $bam && !defined $vcf) {
	$bam=abs_path($bam);
	my $bam_name=basename$bam;
	$bam_name=~s/\.bam$//;
	if ($bam_name!~/sort/){
		$cmd .="$samtools sort --threads $threads $bam >$out/$bam_name.sort.bam \n";
		$bam="$out/$bam_name.sort.bam";
	}
	if (defined $downsample) {
		$cmd.="java -jar $picard DownsampleSam I=$bam O=$out/$bam_name.DP$downsample.bam P=$downsample \n";
		$bam="$out/$bam_name.DP$downsample.bam";
		#$cmd.="$samtools depth $bam >$bam.depth &\n";
	}
	$cmd .="$samtools  mpileup -t DP,AD -uvf $genome $bam > $out/$bam_name.DP$downsample.vcf \n";
	$vcf ="$out/$bam_name.DP$downsample.vcf";

	&runcmd("step1",$cmd);
}

######################################################
#精简vcf
unless ($vcf=~m/.modify$/) {
	open V,"$vcf" or die $!;
	open O,">$vcf.modify" or die $!;
	#print O "Chr\tLoc\tRef\tAlt\tTDepth_${tag}\tDRef_${tag}\tDAlt_${tag}\n"; 
	print O "Chr\tLoc\tRef\tAlt\tTDepth\tDRef\tDAlt\n";
	while (<V>) {
		chomp;
		next if (/^#/ || /^[[<]/ );	
		my @tmp= (split /\t/,$_);
		next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
		my $ref=uc($tmp[3]);
		my @alt= (split ",",$tmp[4])[0];
		next if (length($alt[0]) !=1);
		#pop @alt;
		my ($dp,$ad) = (split ":",$tmp[9])[1,2];	
		my @ad = split "," ,$ad;
		next if ($ad[0]*$ad[1]==0);
		$alt[0]=uc($alt[0]);
		print O "$tmp[0]\t$tmp[1]\t$tmp[3]\t$alt[0]\t$dp\t$ad[0]\t$ad[1]\n";
	}
	close V;
	close O;
	$vcf="$vcf.modify";
}
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
open E  ,"$errortable" or die $!;
while (<E>) {
	chomp;
	next if (/^Genotype/);
	my ($genotype,$base,$fre)=(split "\t",$_)[0,1,4];
	$ERROR{$base}{$genotype}=$fre;
}
close E;

open V ,"$vcf" or die $!;
open O ,">$out/genotype" or die $!; 
while (<V>) {
	chomp;
	next if (/TDepth/ );	
	my ($Chr,$Loc,$ref,$alt,$dp,$Dref,$DAlt)= (split /\t/,$_);
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
		my $p=$SNP{$alt_UP}{$gt[$i]}*(($ERROR{$ref_UP}{$gt[$i]}**$Dref)*($ERROR{$alt_UP}{$gt[$i]}**$DAlt));  #
		if ($p >= $p_max) {
			$genotype=$gt[$i];
			$p_max=$p;
		}
	}
	print O "$Chr\t$Loc\t$ref\t$genotype\t$dp\t$Dref\t$DAlt\t$alt_UP\n";
}
close V;
close O;

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
Usage:    perl $0 -vcf /56T/Evan/PMCaller/mydata/all/DNA.C10T_1.sorted.vcf -o  out  -snptable  /56T/Evan/PMCaller/mydata/all/DNA.C10T_11_list.snp.table -errortable  /56T/Evan/PMCaller/mydata/WES_error_table/DNA.C10T/DNA.C10T.30.eout.11.error.table.list 
		or    perl $0  -snptable  /56T/Evan/PMCaller/mydata/all/DNA.C10T_11_list.snp.table -errortable  /56T/Evan/PMCaller/mydata/WES_error_table/DNA.C10T/DNA.C10T.30.eout.11.error.table.list  -bam /56T/Tania/P#180/GATK_single/DNA.C10T_1.sorted.bam  -downsample 0.1 -o out3 
Function: Template for Perl  
Command:      
    'verbose' => \$verbose,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'genome=s' => \$genome,
    'snptable=s' => \$snptable,
    'errortable=s' => \$errortable,
    'downsample=f' => \$downsample,
	'matrix' => \$matrix,  如果snp tabel 是10行4列的格式就加上该参数
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}