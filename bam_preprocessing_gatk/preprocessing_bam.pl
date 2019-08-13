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
my ($bam,$out,$genome);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'bam=s' => \$bam,
    'o=s' => \$out, 
    'genome=s' => \$genome,
) or die $!;
unless(defined $bam){&usage();exit 0;}
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Output dir is $out\nbam is $bam \n";
`mkdir -p  $out ` unless (-d "$out");
$genome||="/data/hg19/hg19_index/hg19.fa";
my $picard="/soft/picard-tools-2.3.0/picard.jar";
$out ||="./"; 
my $GenomeAnalysisTK = "/soft/GATK/GenomeAnalysisTK.jar";
=c
for i in `ls /56T/P#180/DNA.*/*bam` 
do 
a=`echo $i|awk -F "/" '{print $4}'|awk -F "." '{print $1}'`
b=`echo $i|awk -F "/" '{print $5}'|awk -F "." '{print $1}'`
sample=$a"."$b
echo " 
#bwa mem -M -t 4 /data/hg19/hg19_index/hg19.fa /56T/Tania/P#180/FASTQ/$sample.trimmed.fastq.gz > /56T/Tania/P#180/GATK/${sample}.sam 
#samtools view $i -h |awk -F \"\\t\" '{print \$1\"\\t\"\$2\"\\tchr\"\$3\"\\t\"\$0}' |cut -f 1-3,7- > /56T/Tania/P#180/GATK/${sample}.sam 

#java -jar /soft/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups I=/56T/Tania/P#180/GATK/${sample}.sam O=/56T/Tania/P#180/GATK/${sample}.sorted.bam SO=coordinate RGID=${sample} RGLB=${sample} RGPL=illumina RGPU=${sample} RGSM=${sample} 
#java -jar /soft/picard-tools-2.3.0/picard.jar MarkDuplicates I=/56T/Tania/P#180/GATK/${sample}.sorted.bam O=/56T/Tania/P#180/GATK/${sample}.marked_duplicates.bam CREATE_INDEX=true M=/56T/Tania/P#180/GATK/${sample}.marked_dup_metrics.txt 
#java -jar /soft/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /data/hg19/hg19_index/hg19.fa -I /56T/Tania/P#180/GATK/${sample}.marked_duplicates.bam --known /soft/GATK/dbsnp/hg19/1000G_phase1.indels.hg19.vcf --known /soft/GATK/dbsnp/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o /56T/Tania/P#180/GATK/${sample}.mkdup.intervals -nt 8 
#java -jar /soft/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R /data/hg19/hg19_index/hg19.fa -I /56T/Tania/P#180/GATK/${sample}.marked_duplicates.bam -known /soft/GATK/dbsnp/hg19/1000G_phase1.indels.hg19.vcf -known /soft/GATK/dbsnp/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -targetIntervals /56T/Tania/P#180/GATK/${sample}.mkdup.intervals -o /56T/Tania/P#180/GATK/${sample}.mkdup.realigned.bam 
#java -jar /soft/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R /data/hg19/hg19_index/hg19.fa -I /56T/Tania/P#180/GATK/${sample}.mkdup.realigned.bam -knownSites /soft/GATK/dbsnp/hg19/All_20160601.vcf.gz -knownSites /soft/GATK/dbsnp/hg19/1000G_phase1.indels.hg19.vcf -knownSites /soft/GATK/dbsnp/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o /56T/Tania/P#180/GATK/${sample}.recal_data.table -nct 4 
#java -jar /soft/GATK/GenomeAnalysisTK.jar -T PrintReads -R /data/hg19/hg19_index/hg19.fa -I /56T/Tania/P#180/GATK/${sample}.mkdup.realigned.bam -BQSR /56T/Tania/P#180/GATK/${sample}.recal_data.table -o /56T/Tania/P#180/GATK/${sample}.recal_reads.bam -nct 4 

#java -jar /soft/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/hg19/hg19_index/hg19.fa -I /56T/Tania/P#180/GATK/${sample}.recal_reads.bam -o /56T/Tania/P#180/GATK/${sample}.vcf
#java -jar /soft/GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/hg19/hg19_index/hg19.fa -V /56T/Tania/P#180/GATK/${sample}.vcf -selectType SNP  -o /56T/Tania/P#180/GATK/${sample}.raw.snp.vcf 
java -jar /soft/GATK/GenomeAnalysisTK.jar -T VariantFiltration  -R /data/hg19/hg19_index/hg19.fa -V  /56T/Tania/P#180/GATK/${sample}.raw.snp.vcf --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"my_snp_filter\"  -o  /56T/Tania/P#180/GATK/${sample}.filter.snp.vcf 
grep -w -v my_snp_filter /56T/Tania/P#180/GATK/${sample}.filter.snp.vcf > /56T/Tania/P#180/GATK/${sample}.PASS.snp.vcf 
#java -jar /soft/GATK/GenomeAnalysisTK.jar -T SelectVariants -R /data/hg19/hg19_index/hg19.fa -V /56T/Tania/P#180/GATK/${sample}.vcf -selectType INDEL  -o /56T/Tania/P#180/GATK/${sample}.raw.InDel.vcf 
java -jar /soft/GATK/GenomeAnalysisTK.jar -T VariantFiltration  -R /data/hg19/hg19_index/hg19.fa -V  /56T/Tania/P#180/GATK/${sample}.raw.InDel.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"  --filterName \"my_indel_filter\"  -o  /56T/Tania/P#180/GATK/${sample}.filter.InDel.vcf 
grep -w -v my_indel_filter /56T/Tania/P#180/GATK/${sample}.filter.InDel.vcf > /56T/Tania/P#180/GATK/${sample}.PASS.InDel.vcf 
" > /56T/Tania/P#180/shell/$sample.GATK.sh 
done
=cut
my $sample=basename($bam);
$sample=~s/{(.sort(ed)?)?.bam$//;

my $cmd ="";
$cmd .= "java -jar $picard AddOrReplaceReadGroups I=$bam O=$out/$sample.sorted.bam  SO=coordinate RGID=${sample} RGLB=${sample} RGPL=illumina RGPU=${sample} RGSM=${sample} \n";
$cmd .= "java -jar $picard MarkDuplicates I=$out/$sample.sorted.bam  O=$out/${sample}.marked_duplicates.bam CREATE_INDEX=true M=$out/${sample}.marked_dup_metrics.txt \n";
$cmd .= "java -jar $GenomeAnalysisTK -T RealignerTargetCreator -R $genome -I $out/${sample}.marked_duplicates.bam --known /soft/GATK/dbsnp/hg19/1000G_phase1.indels.hg19.vcf --known /soft/GATK/dbsnp/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o $out/${sample}.mkdup.intervals -nt 8 \n";
$cmd .= "java -jar $GenomeAnalysisTK -T IndelRealigner -R $genome -I $out/${sample}.marked_duplicates.bam -known /soft/GATK/dbsnp/hg19/1000G_phase1.indels.hg19.vcf -known /soft/GATK/dbsnp/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -targetIntervals $out/${sample}.mkdup.intervals -o $out/${sample}.mkdup.realigned.bam \n"; 
$cmd .= "java -jar $GenomeAnalysisTK -T BaseRecalibrator -R $genome -I $out/${sample}.mkdup.realigned.bam -knownSites /soft/GATK/dbsnp/hg19/All_20160601.vcf.gz -knownSites /soft/GATK/dbsnp/hg19/1000G_phase1.indels.hg19.vcf -knownSites /soft/GATK/dbsnp/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o $out/${sample}.recal_data.table -nct 4 \n"; 
$cmd .= "java -jar $GenomeAnalysisTK -T PrintReads -R $genome -I $out/${sample}.mkdup.realigned.bam -BQSR $out/${sample}.recal_data.table -o $out/${sample}.recal_reads.bam -nct 4  \n";

&runcmd("bam_preprocessing",$cmd);
#perl /56T/Evan/PMCaller/script/PMCaller_exonint_v4.2.pl -o /56T/Evan/PMCaller/NA12878/averagedepth/chr22_reads1/ -L 3:5:9:15:19 -cut_H 300 -window 0 -bam /56T/Evan/PMCaller/NA12878/chr22/chr22.filter.sorted.bam -downsample 0.05 -quality 15 -vcf /56T/Evan/PMCaller/NA12878/chr22/chr22.filter.vcf -reads_cut 1  & 


sub runcmd{
	my $name=shift @_;
	my $cmd=shift @_;
	`mkdir "$out/shell/"` unless (-d "$out/shell");
	open S ,">$out/shell/$name.sh" or die $!;
	print S "$cmd ";
	print "Start analysis of $name... \n";
	system "sh $out/shell/$name.sh "  ; #unless ($printcmd) ;

	close S;
}

sub usage {
    die(
        qq!
Usage:    perl $0  -o   /56T/Evan/PMCaller/NA12878/probam/chr22_exonint/   -bam  /56T/Evan/PMCaller/NA12878/chr22_exonint/down0.05/chr22.filter.sorted.DP0.05.bam
Function: Template for Perl  处理bam文件，参考gatk，将bam文件根据已知的数据库进行预处理
Command:      
	verbose	str	verbose,
    'verbose' => \$verbose,
    'bam=s' => \$bam,
    'o=s' => \$out, 
    'genome=s' => \$genome,
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2018/8/14
Notes: 
\n!
    )
}