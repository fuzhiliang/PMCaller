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
my ($L_list,$peak23,$peak12,$vcf,$out,$threads);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'vcf=s' => \$vcf,   #假的vcf，是扩增后vcf_expansion_H60_120.pl 的中间文件
    'o=s' => \$out, 
    'peak12=s' => \$peak12,
    'peak23=s' => \$peak23,
    'L=s' => \$L_list,
	'threads=i' => \$threads,
) or die $!;
unless(defined $vcf  ){&usage();exit 0;}

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $vcf \nOutput file is $out\n";
$out||='./';
system "mkdir $out " unless (-d $out);

my (@l)=(split ":" ,$L_list);

my $maxski =&max(@l);
warn "max ski is $maxski\n" if ($maxski >30);

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
my $hh;

open I ,"$vcf" or die $!;
while (<I>) {
	chomp;
	next if (/^chrom/);
#	my ($chr,$loc,$ref,$alt,$depth,$ref_depth,$alt_depth)= (split "\t", $_)[0,1,2,3,7,8,9];
	my ($chr,$loc,$ref,$alt,$depth,$ref_depth,$alt_depth)= (split "\t", $_)[0,1,2,3,4,5,6];  #扩增前真实数据削薄
	$depth=$ref_depth+$alt_depth;
	my $CL="$chr\t$loc\t$ref\t$alt\t$depth\t$ref_depth\t$alt_depth";
	$peak12||=0.25*$depth ;
	$peak23||=0.75*$depth ;
	$hh=$depth;
	next if ($l[0] >$depth);
	$depth=$ref_depth+$alt_depth ;
	my $genotype;
	if ($ref_depth<=$peak12){
	    $genotype="$alt$alt" ;
	}elsif($ref_depth<=$peak23){
		my @a=("$ref","$alt");
		@a=sort(@a);
		$genotype="$a[0]$a[1]";
	}else{
		$genotype="$ref$ref";
	}
	#削薄
	my %sns=();
	my %L_number;
	$L_number{$alt}=0;
	$L_number{$ref}=0;
	for(my $i = 0;$i <$l[0]; $i ++) {
		my $no;
		my $range =$depth;
		do{
			$no = rand($depth);
			$no = int($no)+1;   
		}while(exists $sns{"$no"});
		$sns{"$no"} = 1;
		#print "$no\n";
		if ($no<=$ref_depth) {
			$L_number{$ref}++;
		}elsif($no<= ($ref_depth+$alt_depth)){
			$L_number{$alt}++;
			}
	}
	$hash{$genotype}{$CL}="$l[0]\t$L_number{$ref}\t$L_number{$alt}";
	$hash2{$ref}{$L_number{$ref}}{$alt}++;
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
	open N ,">$out/$sample.H$l[0].ref$ref.snptable.list" or die $!;
	my $head;
	for (my $i=1 ;$i<$l[0] ;$i++) {
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
my $ATCG=join "\t",@ref;
print S "Snp_table\t$ATCG\n";
foreach my $genotype ( @genotype) {
	my $p_list;
	foreach my $ref (@ref) {
		my $p;
		if (exists $hash_snp{$ref}{$genotype}) {
		$p=$hash_snp{$ref}{$genotype};
#$p=$hash_snp{$ref}{$genotype}/$count_total{$ref};
		}else{
			$p=0;
		}
		$p_list .="\t$p";
	}
	print S "$genotype$p_list\n";

}
close S;


=c
对于SNP table  基于扩增后的”vcf“文件，削薄到指定L，得到snptable list和做峰图数据
for i in 10  # {5..20}
do 
perl /56T/Evan/PMCaller/script/bin/ski_Expansionvcf_bam_H_L.pl -vcf /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H80_80/H80.vcflist -o b -peak12 13 -peak23 59 -L $i &
done

for i in {5..20} 
do
sh /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H80_80/shell/skiH80_L$i.sh 
done

for i in 5 10 15 20 30 50 
do 
echo " perl /56T/Evan/PMCaller/script/bin/ski_Expansionvcf_bam_H_L.pl -vcf /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H80_80/H80.vcflist -o /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H80_80/allskiH80_L$i -peak12 13 -peak23 59 -L $i & " >/56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H80_80/shell/skiH80_L$i.sh 
done
for i in 5 10 15 20 30 50 
do 
sh /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H80_80/shell/skiH80_L$i.sh 
done


for i in 50 60 70 80 90 100 110 120 130 
do 
nohup perl /56T/Evan/PMCaller/script/bin/ski_Expansionvcf_bam_H_L.pl -vcf /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H50_160/H50_160.vcflist  -o /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/v2/allskiH50_160_L$i -L $i &
done




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
Function: Template for Perl
Command:      
    'verbose' => \$verbose,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'peak12=s' => \$peak12,
    'peak23=s' => \$peak23,
    'L=s' => \$L_list,
	'threads=i' => \$threads,
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}