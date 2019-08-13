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
my ($L,$cut_peak2,$cut_peak3,$vcf,$out,$cut_H,$ski,$exp,$windows);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'vcf=s' => \$vcf,   #假的vcf，是扩增后vcf_expansion_H60_120.pl 的中间文件
    'o=s' => \$out, 
    'cut_peak2=s' => \$cut_peak2,
    'cut_peak3=s' => \$cut_peak3,
    'cut_H=s' => \$cut_H,  #多少x以上的数据可以用
    'L=s' => \$L,            #削薄到多少x
	'windows=s' => \$windows, 
	'exp' => \$exp,
	'ski' => \$ski,
) or die $!;
unless(defined $vcf  && $out){&usage();exit 0;}

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $vcf \nOutput file is $out\n";
$out||='./';
system "mkdir -p $out " unless (-d $out);

#my (@l)=(split ":" ,$L_list);

#my $maxski =&max(@l);
#warn "max ski is $L\n" if ($L >30);

#$h ||="50";
#gy:
#目标
#H	L	Genotype        Base    Frequency       total   %F
#30	1	AA      T       16070   60264   0.266660029204832
#		AA      C       20206   60264   0.335291384574539
#		AA      N       47      60264   0.000779901765564848
#		AA      A       5478    60264   0.090900039824771
#		AA      G       18463   60264   0.306368644630293
#		AC      G       15      408     0.0367647058823529
#		AC		A       216     408     0.529411764705882
#输入
#chrom   loction ref     alt     depth   ref_depth       alt_depth       exp_depth       exp_ref_depth   exp_alt_depth
#chrX    83664   G       C       80      79      1       80      79      1
#chrX    93703   G       A       80      79      1       80      79      1
#chrX    155178  G       A       80      59      21      80      59      21
#chrX    155379  A       C       80      79      1       80      79      1
#chrX    160306  G       C       80      77      2       80      78      2
#chrX    163281  C       T       80      78      2       80      78      2


print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $sample=basename($vcf);
$sample =~s/\.vcf$//;
my %hash;
my %hash2;
my %hash_snp;
my %count_total;
$cut_peak2||=0.25;
$cut_peak3||=0.73;
$cut_H ||=70;
$windows||=20;
$L||=$cut_H;

open I ,"$vcf" or die $!;
while (<I>) {
	chomp;
	next if (/^chrom/);
#	my ($chr,$loc,$ref,$alt,$depth,$ref_depth,$alt_depth)= (split "\t", $_)[0,1,2,3,7,8,9];
	my ($chr,$loc,$ref,$alt,$depth,$ref_depth,$alt_depth)= (split "\t", $_)[0,1,2,3,4,5,6];  #扩增前真实数据削薄
	$depth=$ref_depth+$alt_depth;
	my $CL="$chr\t$loc\t$ref\t$alt\t$depth\t$ref_depth\t$alt_depth";
	my $peak12||=$cut_peak2*$depth ;
	my $peak23||=$cut_peak3*$depth ;
	if (defined $ski) {  #只削薄
		next if ($depth > ($cut_H+$windows) || $depth < $cut_H);
		print "$depth\n" die;
	}elsif (defined $exp) {#只扩增
		next if ($depth > $cut_H || $depth < ($cut_H-$windows));
	}else{
		next if ($depth > ($cut_H+$windows) || $depth < ($cut_H-$windows) );
	}
	my $genotype;   #根据>H70的位点计算genotype
	if ($ref_depth<=$peak12){
	    $genotype="$alt$alt" ;
	}elsif($ref_depth<=$peak23){
		my @a=("$ref","$alt");
		@a=sort(@a);
		$genotype="$a[0]$a[1]";
	}else{
		$genotype="$ref$ref";
	}
	#削薄到$L x
	my %sns=();
	my %L_number;
	$L_number{$alt}=0;
	$L_number{$ref}=0;
	for(my $i = 0;$i <$L; $i ++) {
		my $no;
		#next if ($depth<$L);
		if ($depth<$L) {
			$no = int(rand($depth));
		}elsif($depth>=$L){
		#srand($depth+$i);
			do{	
				$no = rand($depth);
				$no = int($no)+1;   
			}while(exists $sns{"$no"});
			$sns{"$no"} = 1;
		}
		#print "$no\n";
		if ($no<=$ref_depth) {
			$L_number{$ref}++;
		}elsif($no<= ($ref_depth+$alt_depth)){
			$L_number{$alt}++;
			}
	}
	$hash{$genotype}{$CL}="$L\t$L_number{$ref}\t$L_number{$alt}";
	$hash2{$ref}{$L_number{$ref}}{$alt}++;          #非dubious位点存到hash了，但是后面没输出来
	$hash_snp{$ref}{$genotype}++ if ($L_number{$alt} != 0); 
	$count_total{$ref}++ if ($L_number{$alt} != 0); 
}

close I;
open N ,">$out/$sample.snptable.list" or die $!;
print N "chrom\tloction\tref\talt\tdepth\tref_depth\talt_depth\tski_depth\tski_ref_depth\tski_alt_depth\tgenotype\n";
foreach my $genotype (sort {$a cmp $b } keys %hash) {
	foreach my $CL (keys %{$hash{$genotype}}) {
		print N "$CL\t$hash{$genotype}{$CL}\t$genotype\n";
	}
}

my @ref=("A","C","G","T");

foreach my $ref (@ref) {
	my $fig=0;
	open N ,">$out/$sample.H$L.ref$ref.snptable.list" or die $!;
	my $head;
	for (my $i=1 ;$i<$L ;$i++) {
		my $alt_t="";
		$head="";
		foreach my $alt (@ref) {
			next if ($ref eq $alt);
			$head .="$alt\t";
			if (!exists $hash2{$ref}{$i}{$alt}) {
				$hash2{$ref}{$i}{$alt}=0;
			}
		$alt_t.="\t$hash2{$ref}{$i}{$alt}";
		}
	$head=~s/\t$//;
	print N "$head\n" if ($fig==0);
	$fig=1;
	print N "$i$alt_t\n";
	}
close N ;
}

my @genotype=("AA","AC","AG","AT","CC","CG","CT","GG","GT","TT");
open S ,">$out/$sample.snp.table" or die $!;  #snptable
open R ,">$out/$sample.snp.table_rate" or die $!;  #snptable
my $ATCG=join "\t",@ref;
print S "Snp_table\t$ATCG\n";
print R "Snp_table\t$ATCG\n";
foreach my $genotype ( @genotype) {
	my $p_list;
	my $rate_list;
	foreach my $ref (@ref) {
		my $p;
		my $rate;
		if (exists $hash_snp{$ref}{$genotype}) {
			$p=$hash_snp{$ref}{$genotype};
			$rate=$hash_snp{$ref}{$genotype}/$count_total{$ref};
		}else{
			$p=0;
			$rate=0;
		}
		$p_list .="\t$p";
		$rate_list .="\t$rate";
	}
	print S "$genotype$p_list\n";
	print R "$genotype$rate_list\n";
}
close S;

open L ,">$out/$sample.snp.table_list" or die $!;
print L "Ref\tGenotype\tFrequency\tTotal\t%F\n";
foreach my $ref (@ref) {
	foreach my $genotype (@genotype) {
		
		$hash_snp{$ref}{$genotype}=0 unless (exists $hash_snp{$ref}{$genotype});
		my $p=$hash_snp{$ref}{$genotype}/$count_total{$ref};
		$count_total{$ref}=0 unless (exists $count_total{$ref});
		print L "$ref\t$genotype\t$hash_snp{$ref}{$genotype}\t$count_total{$ref}\t$p\n";
	}
}

=c
对于SNP table  基于扩增后的”vcf“文件，削薄到指定L，根据cut_H 确定genotype，得到snptable list

for i in {3..25}
do 
nohup perl /56T/Evan/PMCaller/script/bin/ski_Expansionvcf_bam_H70_L.pl -vcf /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H50_160/H50_160.vcflist  -o /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H70up_L1..30/H70_L$i -L $i &
done

for i in {3..25}
do 
nohup perl /56T/Evan/PMCaller/script/bin/ski_Expansionvcf_bam_H70_L_ski.pl -vcf /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H50_160/H50_160.vcflist  -o /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H70up_L1..30/H70_L$i -L $i &
done
外显子数据

 perl /56T/Evan/PMCaller/script/bin/ski_Expansionvcf_bam_H70_L.pl -vcf /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H50_160/H50_160.vcflist -o /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H70up_L1..30/H70_L$i -L $i

=cut
###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
sub max {
    my $currentMaxValue = shift @_;
    foreach ( @_ ) {
        if ( $_ > $currentMaxValue ) {
            $currentMaxValue = $_;
        }
    }
    return $currentMaxValue;
}
sub min {
    my $currentMaxValue = shift @_;
    foreach ( @_ ) {
        if ( $_ < $currentMaxValue ) {
            $currentMaxValue = $_;
        }
    }
    return $currentMaxValue;
}
sub usage {
    die(
        qq!
Usage:    perl $0  
Function: Template for Perl 基于扩增得到的*exp.vcf 文件，当L==cut_H时，[cut_H,cut_H+windows]削薄到cut_H, [cut_H-windows,cut_H]扩增到$cut_H，得到做峰图数据 ，给了ski，只削薄，给了exp，只扩增
Command:      
    verbose      verbose
    'vcf=s' => \$vcf,   #假的vcf，是扩增后vcf_expansion_H60_120.pl 的中间文件
    'o=s' => \$out, 
    'cut_peak2=s' => \$cut_peak2, #第一个峰与二个峰x轴百分比  [0.25]
    'cut_peak3=s' => \$cut_peak3, #第二个峰与第三个峰x轴百分比  [0.75]
    'cut_H=s' => \$cut_H,         #假想的平均深度
    'L=s' => \$L,				 #不要给
	'windows=s' => \$windows,     #区间
	'exp' => \$exp,
	'ski' => \$ski,
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}