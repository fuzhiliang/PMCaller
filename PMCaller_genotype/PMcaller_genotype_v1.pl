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
my ($bam,$snptable,$errortable_H,$errortable_L,$downsample,$matrix,$vcf,$vcf_H,$vcf_L,$out,$genome,$bam_H,$bam_L);
my %opts;
GetOptions(
    'verbose' => \$verbose,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'genome=s' => \$genome,
    'snptable=s' => \$snptable,
    'errortable_H=s' => \$errortable_H,
    'errortable_L=s' => \$errortable_L,
	'vcf_H=s' => \$vcf_H,
	'vcf_L=s' => \$vcf_L,
	'bam_H=s' => \$bam_H,
	'bam_L=s' => \$bam_L,
    'downsample=f' => \$downsample,
	'matrix' => \$matrix,
) or die $!;
unless(($bam ||$vcf) && ($bam_H ||$vcf_H) && ($bam_L ||$vcf_L)){&usage();exit 0;}
#unless(defined $snptable || defined $errortable_H && defined $errortable_L) {&usage();exit 0;}
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $vcf \nOutput file is $out\n";
my $picard;
$picard="/soft/picard-tools-2.3.0/picard.jar" if (-e "/soft/picard-tools-2.3.0/picard.jar"); 
$picard="/soft/picard-tools-2.4.1/picard.jar" if (-e "/soft/picard-tools-2.4.1/picard.jar"); 

my $samtools ="/usr/local/bin/samtools";
$genome||="/data/hg19/hg19_index/hg19.fa";
$out||='./';
my $threads||=30;
`mkdir $out ` unless (-d $out);
$out = abs_path($out);
##die "defined $downsample and !defined $" if (defined $downsample && !defined $bam );
warn "ignore bam file :$bam and downsample: $downsample because exists vcf file:$vcf .\n" if ((defined $downsample || defined $bam) && defined $vcf) ;

my $cmd;
if (defined $bam && !defined $vcf) {
	$bam=abs_path($bam);
	my $bam_name=basename$bam;
	$bam_name=~s/\.bam$//;
	if ($bam_name!~/sort(ed)?.bam$/){
		$cmd .="$samtools sort --threads $threads $bam >$out/$bam_name.sort.bam \n";
		$bam="$out/$bam_name.sort.bam";
	}
	if (defined $downsample) {
		#my $downsample_HL=$downsample/2;
		$cmd.="java -jar $picard DownsampleSam I=$bam O=$out/$bam_name.DP$downsample.bam P=$downsample \n";
		$bam="$out/$bam_name.DP$downsample.bam";
	}
	$cmd .="$samtools  mpileup -t DP,AD -uvf $genome $bam > $out/$bam_name.DP$downsample.vcf \n";
	$vcf ="$out/$bam_name.DP$downsample.vcf";

	&runcmd("step1.0",$cmd);
}

if (defined  $bam_H  && !defined $vcf_H) {
	$cmd="";
	$bam=abs_path($bam_H);
	my $bam_name=basename$bam;
	$bam_name=~s/\.[sb]am$//;
	if ($bam=~m/sam$/) {
		$cmd.="$samtools view -Sb $bam |$samtools sort --threads $threads  >$out/$bam_name.sort.bam \n";
		$bam="$out/$bam_name.sort.bam";
	}elsif ($bam_name!~/sort/){
		$cmd .="$samtools sort --threads $threads $bam >$out/$bam_name.sort.bam \n";
		$bam="$out/$bam_name.sort.bam";
	}
	if (defined $downsample) {
		my $downsample=$downsample*2;
		$cmd.="java -jar $picard DownsampleSam I=$bam O=$out/$bam_name.DP$downsample.bam P=$downsample \n";
		$bam="$out/$bam_name.DP$downsample.bam";
	}
	$cmd .="$samtools  mpileup -t DP,AD -uvf $genome $bam > $out/$bam_name.DP2x$downsample.vcf \n";
	$vcf_H ="$out/$bam_name.DP2x$downsample.vcf";

	&runcmd("step1.1",$cmd);
}

if (defined $bam_L && !defined $vcf_L) {
	$cmd="";
	$bam=abs_path($bam_L);
	my $bam_name=basename$bam;
	$bam_name=~s/\.[sb]am$//;
	if ($bam=~m/sam$/) {
		$cmd.="$samtools view -Sb $bam |$samtools sort --threads $threads >$out/$bam_name.sort.bam \n";
		$bam="$out/$bam_name.sort.bam";
	}elsif ($bam_name!~/sort/){
		$cmd .="$samtools sort --threads $threads $bam >$out/$bam_name.sort.bam \n";
		$bam="$out/$bam_name.sort.bam";
	}
	if (defined $downsample) {
		my $downsample=$downsample*2;
		$cmd.="java -jar $picard DownsampleSam I=$bam O=$out/$bam_name.DP$downsample.bam P=$downsample \n";
		$bam="$out/$bam_name.DP$downsample.bam";
		#$cmd.="$samtools depth $bam >$bam.depth &\n";
	}

	$cmd .="$samtools  mpileup -t DP,AD -uvf $genome $bam > $out/$bam_name.DP2x$downsample.vcf \n";
	$vcf_L ="$out/$bam_name.DP2x$downsample.vcf";

	&runcmd("step1.2",$cmd);
}



#把high quality 和low quality 的vcf存到hash
my %Qhigh;
my %Qlow;
open HIGH ,"$vcf_H" or die $!;
while (<HIGH>) {
	chomp;
	next if (/^#/);
	my @tmp= (split /\t/,$_);
	next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
	my $ref=$tmp[3];
	my $alt= (split ",",$tmp[4])[0];
	next if (length($alt) !=1);
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];
	my @ad = split "," ,$ad;
	my $loc="$tmp[0]\t$tmp[1]";
	#my $bb="$ref\t$alt\t$ad[0]\t$ad[1]";
	$Qhigh{$loc}{$ref}=$ad[0];
	$Qhigh{$loc}{$alt}=$ad[1];
}
close HIGH;
open LOW ,"$vcf_L" or die $!;
while (<LOW>) {
	chomp;
	next if (/^#/);
	my @tmp= (split /\t/,$_);
	next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
	my $ref=$tmp[3];
	my $alt= (split ",",$tmp[4])[0];
	next if (length($alt) !=1);
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];
	my @ad = split "," ,$ad;
	my $loc="$tmp[0]\t$tmp[1]";
	#my $bb="$ref\t$alt\t$ad[0]\t$ad[1]";
	$Qlow{$loc}{$ref}=$ad[0];
	$Qlow{$loc}{$alt}=$ad[1];
}
close LOW;

open N  ,"$snptable" or die $!;
my %SNP;


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

#处理High Quqlity的error table
my %ERROR_H;
open EH  ,"$errortable_H" or die $!;
while (<EH>) {
	chomp;
	next if (/^Genotype/);
	my ($genotype,$base,$fre)=(split "\t",$_)[0,1,4];
	$ERROR_H{$base}{$genotype}=$fre;
}
close EH;

#处理Low Quqlity的error table
my %ERROR_L;
open EL  ,"$errortable_L" or die $!;
while (<EL>) {
	chomp;
	next if (/^Genotype/);
	my ($genotype,$base,$fre)=(split "\t",$_)[0,1,4];
	$ERROR_L{$base}{$genotype}=$fre;
}
close EL;


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

	unshift (@alt,$ref);
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];
	my @ad = split "," ,$ad;
	my $loc="$tmp[0]\t$tmp[1]";
	my $p_max=0;
	my $genotype;
	for (my $i=0;$i<=$#alt ;$i++) {
		for (my $j=0;$j<=$#alt ;$j++) {
			my @gt=("$alt[$i]","$alt[$j]");
			@gt=sort@gt;
			my $gt=join "",@gt;

			$Qhigh{$loc}{$alt[0]}="0" if (!exists $Qhigh{$loc}{$alt[0]}) ;
			$Qhigh{$loc}{$alt[1]}="0"  if (!exists $Qhigh{$loc}{$alt[1]}) ;
			$Qlow{$loc}{$alt[0]}="0"  if (!exists $Qlow{$loc}{$alt[0]}) ;
			$Qlow{$loc}{$alt[1]}="0"  if (!exists $Qlow{$loc}{$alt[1]}) ;

			my $p=$SNP{$alt[0]}{$gt}*($ERROR_H{$alt[0]}{$gt}**$Qhigh{$loc}{$alt[0]})*($ERROR_H{$alt[1]}{$gt}**$Qhigh{$loc}{$alt[1]})*($ERROR_L{$alt[0]}{$gt}**$Qlow{$loc}{$alt[0]})*($ERROR_L{$alt[1]}{$gt}**$Qlow{$loc}{$alt[1]});  #

			#my $p=$SNP{$alt[0]}{$gt}*($ERROR{$alt[0]}{$gt}**$ad[0])*($ERROR{$alt[1]}{$gt}**$ad[1]);  #
			print A "$p=\t$SNP{$alt[0]}{$gt}*($ERROR_H{$alt[0]}{$gt}**$Qhigh{$loc}{$alt[0]})*($ERROR_H{$alt[1]}{$gt}**$Qhigh{$loc}{$alt[1]})*($ERROR_L{$alt[0]}{$gt}**$Qlow{$loc}{$alt[0]})*($ERROR_L{$alt[1]}{$gt}**$Qlow{$loc}{$alt[1]})\n$alt[0]\t$alt[1]\t$gt\t$ad[0]\t$ad[1]\n$gt\t$p_max\t$p\t$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t@alt\t$tmp[9]\n";
			if ($p >= $p_max) {
				$genotype=$gt;
				$p_max=$p;
			}
			#print A "$genotype\t$gt\t$p_max\t$p\t$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t@alt\t$tmp[9]\n";
		}
	}
	print O "$tmp[0]\t$tmp[1]\t$tmp[3]\t$genotype\t$dp\t$ad[0]\t$ad[1]\t$Qlow{$loc}{$alt[0]}\t$Qlow{$loc}{$alt[1]}\t$Qhigh{$loc}{$alt[0]}\t$Qhigh{$loc}{$alt[1]}\n";
	print C "$genotype\t$_\n";

#	my $max_ad=&max(@ad);
#	my $genotype;
#	if ($ad[0] >= $dp*(1-$cut_off)) {
#		$genotype="$tmp[3]$tmp[3]";
#	}elsif ($ad[0] <= $dp*(0.5+0.5*$cut_off) && $ad[0] >= $dp*(0.5-0.5*$cut_off) ) {
#		my @a=("$tmp[3]","$alt[0]");
#		@a=sort @a;
#		$genotype=join "",@a;
#	}elsif ($ad[0] <= $dp*$cut_off && $ad[1] >= $dp*(1-$cut_off) ) {
#		$genotype="$alt[0]$alt[0]";
#	}elsif (length($alt[1])==1 && length($alt[2]) ==1 ) {
#		$alt[1]=~s/(\w+)/\u$1/;
#		$alt[2]=~s/(\w+)/\u$1/;
#		my @b =("$alt[1]","$alt[2]");
#		@b=sort @b;
#		$genotype=join "",@b;
#	}else{next;}
#	print O "$tmp[0]\t$tmp[1]\t$tmp[3]\t$genotype\t$dp\n";
}
close V;
close O;
close C;
close A;
=c

cd /56T/Evan/PMCaller/GATK/DownSample
for sample in DNA.C10T DNA.C9T DNA.P10T DNA.P9T
do
for j in 0.20
do
echo "
#java -jar /soft/picard-tools-2.3.0/picard.jar DownsampleSam  I=/56T/Tania/P#180/GATK_single/${sample}_1.sorted.bam O=/56T/Evan/PMCaller/GATK/DownSample/$sample.DP${j}.bam P=${j}
#samtools depth /56T/Evan/PMCaller/GATK/DownSample/$sample.DP${j}.bam  >  /56T/Evan/PMCaller/GATK/DownSample/$sample.DP${j}.depth
samtools index  /56T/Evan/PMCaller/GATK/DownSample/$sample.DP${j}.bam 
java -Xmx25g -jar /soft/GATK/GenomeAnalysisTK.jar    -R  /data/hg19/hg19_index/hg19.fa   -T UnifiedGenotyper -I /56T/Evan/PMCaller/GATK/DownSample/$sample.DP${j}.bam -o /56T/Evan/PMCaller/GATK/DownSample/$sample.DP${j}.bam.gatk.vcf   -stand_call_conf 0.0 -stand_emit_conf 0.0 -dcov 200  -nt 8  -rf BadCigar 
" > /56T/Evan/PMCaller/GATK/DownSample/shell/$sample.DP.$j.sh
done
done
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
Usage:    perl $0 -vcf /56T/Evan/PMCaller/mydata/chr22_samtools_mpileup.vcf  -snptable  /56T/Evan/PMCaller/mydata/all/DNA.C10T_11_list.snp.table  -errortable_H DNA.C10T.30.eout.11.error.table.list -errortable_L DNA.C10T.30.eout.9.error.table.list  -vcf_H  /56T/Evan/PMCaller/mydata/chr22_samtools_mpileup.vcf -vcf_L /56T/Evan/PMCaller/mydata/chr22_samtools_mpileup.vcf  -o Qout 
 
Function: Template for Perl  分高质量和低质量reads call genotype 
Command:      
    'verbose' => \$verbose,
	'bam=s' => \$bam,  #不分b质量值的bam 
    'vcf=s' => \$vcf,   #不分质量值的vcf，存在则-bam失效
    'o=s' => \$out, 
    'genome=s' => \$genome,
    'snptable=s' => \$snptable,       #从高深度确定genotype,削薄后的snptable
    'errortable_H=s' => \$errortable_H,    
    'errortable_L=s' => \$errortable_L,   
	'vcf_H=s' => \$vcf_H,  #高质量reads 削薄后的vcf ，存在则-bam_H失效
	'vcf_L=s' => \$vcf_L,   #低质量reads 削薄后的vcf ，存在则-bam_L失效
	'bam_H=s' => \$bam_H,
	'bam_L=s' => \$bam_L,
    'downsample=f' => \$downsample,   #随机取bam中的多少reads   [0,1]
	'matrix' => \$matrix,    #snp table是列10行的时候加上该参数。
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}