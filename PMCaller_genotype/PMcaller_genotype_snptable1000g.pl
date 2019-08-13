#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
#use Getopt::Std;
use File::Basename;
use File::Glob;
use FindBin qw($Bin); 
use Data::Dumper;
use Cwd qw(abs_path getcwd);  
my $abs = abs_path(getcwd());  
my $verbose ="v1.0";
###############################################################################
#命令行参数据的定义和获取，记录程序初始时间，设置参数默认值
#Get the parameter and provide the usage.
###cd /56T/Evan/PMCaller/check_data && perl /56T/Evan/PMCaller/script/PMCaller_v1.pl -sam  chr22_head.P10T_1.sam -o test -H 40:50 -L 5:8:9:20
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
=c
#处理vcf
#perl  /56T/Evan/PMCaller/script/PMCaller_genotype/PMCaller_genotype_EXON_INT_preprocess.pl /56T/Evan/PMCaller/NA12878/averagedepth/chr10_exonint//down0.05/chr10.filter.vcf.m  /data/hg19_anno/Region/refgene.nooverlap
head /56T/Evan/PMCaller/NA12878/averagedepth/chr10_exonint//down0.05/chr10.filter.vcf.m.EXONINT.modify
Chr     Loc     Ref     Alt     TDepth  DRef    DAlt    Location
chr10   60410   C       T       19      18      1       INT
chr10   60563   T       A       10      9       1       INT
chr10   60565   A       T       11      10      1       INT

#准备snp table 和error table

#snp list
Chr	loc	AA	AC	AG	AT	CC	CG	CT	GG	GT	TT
chr22	16050159	0.1	0.2	0.3	0.04	0.01	0.04	0.01	0.1	0.1
chr22	16050195	0.1	0.2	0.3	0.04	0.01	0.04	0.01	0.1	0.1
chr22	16050205	0.1	0.2	0.3	0.04	0.01	0.04	0.01	0.1	0.1
chr22	16050256	0.1	0.2	0.3	0.04	0.01	0.04	0.01	0.1	0.1
chr22	16050290	0.1	0.2	0.3	0.04	0.01	0.04	0.01	0.1	0.1
chr22	16050299	0.1	0.2	0.3	0.04	0.01	0.04	0.01	0.1	0.1
chr22	16050396	0.1	0.2	0.3	0.04	0.01	0.04	0.01	0.1	0.1


=cut
my %SNP;
my %ERROR;
open N  ,"$snptable" or die $!;
my $head =<N>;
chomp $head;
my @head=(split "\t", $head);
shift @head;
shift @head;
while (<N>) {
	chomp ;
	next if (/AA/);
	my @tmp=split "\t", $_ ;
	my $chr=shift @tmp;
	my $loc=shift @tmp;
	for (my $i=0;$i<10 ;$i++) {
		$SNP{$chr}{$loc}{$head[$i]}=$tmp[$i];
	}
}

#print Dumper %SNP;
=c

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
=cut
open E  ,"$errortable" or die $!;
while (<E>) {
	chomp;
	next if (/^Genotype/);
	my ($genotype,$base,$fre)=(split "\t",$_)[0,1,4];
	$ERROR{$base}{$genotype}=$fre;
}
close E;

open O ,">$out/genotype_snp1000g" or die $!; 
open V ,"$vcf" or die $!;
while (<V>) {
	chomp ;
	next if (/Ref/);
	my ($chr,$loc,$ref,$alt,$depth,$dref,$dalt)=split "\t",$_ ;
	my $p_max=0;
	my $genotype;
	my @het=("$ref","$alt");
	@het=sort(@het);
	my $het=join "",@het;
	my @gt=("$ref$ref","$het","$alt$alt");
	print @gt;
	print "genotype\n\n\n";
	next if (!exists $SNP{$chr}{$loc});
	foreach my $gt (@gt) {
		$ERROR{$ref}{$gt}=0 if (!exists $ERROR{$ref}{$gt});
		$ERROR{$alt}{$gt}=0 if (!exists $ERROR{$alt}{$gt});
		my $p = $SNP{$chr}{$loc}{$gt}*($ERROR{$ref}{$gt}**$dref)*($ERROR{$alt}{$gt}**$dalt);
		print "$p = $SNP{$chr}{$loc}{$gt}*($ERROR{$ref}{$gt}**$dref)*($ERROR{$alt}{$gt}**$dalt) \n";
		if ($p >= $p_max) {
			$genotype=$gt;
			$p_max=$p;
		}
	}
	print O "$chr\t$loc\t$ref\t$genotype\t$depth\t$dref\t$dalt\n";
	die;
}

close O;
close V;

=c
open V ,"$vcf" or die $!;
open O ,">$out/genotype" or die $!; 
open C ,">$out/genotype.all" or die $!; 
open A ,">$out/check.list" ;
while (<V>) {
	chomp;
	next if (/^#/ || /^[[<]/ );	
	my @tmp= (split /\t/,$_);
	next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
	my $ref=$tmp[3];
	my @alt= (split ",",$tmp[4])[0];
	next if (length($alt[0]) !=1);
	#pop @alt;
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];	
	my @ad = split "," ,$ad;
	unshift (@alt,$ref);
#	print @alt;die;
	#ref=$ad[0]; ale[0]=$ad[1];
	my $p_max=0;
	my $genotype;
	for (my $i=0;$i<=$#alt ;$i++) {
		for (my $j=0;$j<=$#alt ;$j++) {
			my @gt=("$alt[$i]","$alt[$j]");
			@gt=sort@gt;
			my $gt=join "",@gt;
			$ERROR{$alt[0]}{$gt}=0 if (!exists $ERROR{$alt[0]}{$gt});
			$ERROR{$alt[1]}{$gt}=0 if (!exists $ERROR{$alt[1]}{$gt});
			my $p=$SNP{$alt[0]}{$gt}*($ERROR{$alt[0]}{$gt}**$ad[0])*($ERROR{$alt[1]}{$gt}**$ad[1]);  
			#print A "$p=\t$SNP{$alt[0]}{$gt}*($ERROR{$alt[0]}{$gt}**$ad[0])*($ERROR{$alt[1]}{$gt}**$ad[1])\n$alt[0]\t$alt[1]\t$gt\t$ad[0]\t$ad[1]\n$gt\t$p_max\t$p\t$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t@alt\t$tmp[9]\n";
			if ($p >= $p_max) {
				$genotype=$gt;
				$p_max=$p;
			}
			#print A "$genotype\t$gt\t$p_max\t$p\t$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t@alt\t$tmp[9]\n";
		}
	}
	print O "$tmp[0]\t$tmp[1]\t$tmp[3]\t$genotype\t$dp\t$ad[0]\t$ad[1]\n";
	#print C "$genotype\t$_\n";

}
close V;
close O;
close C;
close A;
=cut
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
Usage:    perl $0 -vcf /56T/Evan/PMCaller/NA12878/averagedepth/chr22_exonint/down0.05/chr22.filter.vcf.m.EXONINT.modify  -o /56T/Evan/PMCaller/NA12878/averagedepth/chr22_exonint/PMCallerdown5/  -snptable  1kg.list  -errortable /56T/Evan/PMCaller/NA12878/averagedepth/chr22_exonint/errortable/chr22.filter.errortable.list_5

Function: Template for Perl  
Command:      
    'verbose' => \$verbose,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'genome=s' => \$genome,
    'snptable=s' => \$snptable, #1KG频率
    'errortable=s' => \$errortable,
    'downsample=f' => \$downsample,
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}