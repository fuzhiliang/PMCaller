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
###
###############################################################################
my ($Hmax,$Hmin,$vcf_list,$name,$out,$splitchr);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'vcf_list=s' => \$vcf_list,
    'o=s' => \$out, 
    'Hmin=i' => \$Hmin,
    'Hmax=i' => \$Hmax,
	'splitchr' => \$splitchr,
	'name=s' =>\$name,
) or die $!;
unless(defined $vcf_list || defined $out ){&usage();exit 0;}
my $start_time=time;
`mkdir $out ` unless (-d $out);
`mkdir "$out/shell" ` unless (-d "$out/shell");
$Hmin||=60;
$Hmax||=120;
my @vcf;
@vcf=(glob "$vcf_list/*.vcf")if (-d $vcf_list);
@vcf=("$vcf_list") if (-f $vcf_list);
my $hh=$Hmax;
foreach my $vcf (@vcf) {
	my $sample=basename($vcf);
	$sample=~s/\.vcf$//;
	$name||=$sample;
	open I ,"$vcf" or die $!;
	open O ,">$out/$sample.exp.vcf" or die $!;
	print O "chrom\tloction\tref\talt\tdepth\tref_depth\talt_depth\texp_depth\texp_ref_depth\texp_alt_depth\n";
	my %hash;
	while (<I>) {
		chomp;
		next if (/^#/ || /^[[<]/ );	
		my @tmp= (split /\t/,$_);
		next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
		my @alt= split ",",$tmp[4];
		next if (length($alt[0]) !=1);
		my ($dp,$ad) = (split ":",$tmp[9])[1,2];	
		my @ad = split "," ,$ad;
		next if ($dp<$Hmin || $dp > $Hmax);
		my $sum=$ad[0]+$ad[1];
		my $ref_ad_exp=$ad[0];
		my $alt_ad_exp=$ad[1];
		for(my $i=$sum;$i<$Hmax;$i++) {
			#my $numble=$Hmax-$sum;
			my $rand_integer = int(rand($sum));
			if($rand_integer <= $ad[0]){
				$ref_ad_exp++ ;
				}else{
				$alt_ad_exp++;
				}
		}
		my $col="$tmp[0]\t$tmp[1]\t$tmp[3]\t$alt[0]\t$dp\t$ad[0]\t$ad[1]\t$Hmax\t$ref_ad_exp\t$alt_ad_exp\n";
		print O "$col";
		$hash{$Hmax}{$tmp[3]}{$ref_ad_exp}{$alt[0]}++;

	}
	close O;
	my @ref=("A","C","G","T");
	#foreach my $hh (sort {$a <=>$b } keys %hash) {
	foreach my $ref (@ref) {
		my $fig=0;
		#open N ,">$out/$sample.H$hh.ref$ref.snptable.list" or die $!;
		my $head;
		for (my $i=0 ;$i<$hh ;$i++) {
			my $alt_t="";
			$head="";
			foreach my $alt (@ref) {
				next if ($ref eq $alt);
				$head .="$alt\t";
				if (!exists $hash{$hh}{$ref}{$i}{$alt}) {
					$hash{$hh}{$ref}{$i}{$alt}=0;
				}
			$alt_t.="\t$hash{$hh}{$ref}{$i}{$alt}";
			}
		$head=~s/\t$//;
		#print N "$head\n" if ($fig==0);
		$fig=1;
		#print N "$i$alt_t\n";
		}
	#close N ;
	}
}

=c
my @ref=("A","C","G","T");
if (defined $splitchr) {
	my $cmd;
#	foreach my $h (@h) {
		foreach my $ref (@ref) {
			$cmd .="perl $Bin/snpchr_merge.pl $out/\\*.H$hh.ref$ref.snptable.list $out/merge/ $name\n";
		}
#	}
	print $cmd;
	system "sh $cmd";
#for  h in  110 120 130 140 150 
#do
#for  b in A C G T 
#do
#a=SRR345592.sam.chr\\*.sam.H$h.ref$b.snptable.list
#echo "perl  /56T/Evan/PMCaller/script/bin/snpchr_merge.pl  /56T/Evan/PMCaller/ncbi_sra/li/4/BWA/SRR345592_chrom_snpout_H130_10pic/$a /56T/Evan/PMCaller/ncbi_sra/li/4/BWA/SRR345592_chrom_merge " >>merge.sh
#done
#done
}else{
#	foreach my $hh(@h) {
	foreach my $ref (@ref) {
		
			my $Rsh="";
			$Rsh .= "data <- read.table(\"$out/$name.H$hh.ref$ref.snptable.list\",sep='\\t',head=T) \n";
			$Rsh .= "pdf(\"$out/$name.H$hh.ref$ref.pdf\") \n";
			$Rsh .= "matplot(rownames(data),data, type =\"o\", cex=0.8, pch = 16, col = c(\"darkgreen\",\"gold\",\"darkblue\",\"magenta\"),main=\"SNP Table: sample $name H=$hh Ref=$ref \",xlab = \"# Reads supporting ref base\", ylab = \"# Position\")\n" ;
			$Rsh .= "legend(\"right\",bty=\"n\",legend=colnames(data),col = c\(\"darkgreen\",\"gold\",\"darkblue\",\"magenta\"),cex=0.8, pch =16) \n";
			$Rsh .= "dev.off()";

			open S ,">$out/shell/$name.H$hh.ref$ref.R";
			print S "$Rsh";
			system "nohup R CMD BATCH  $out/shell/$name.H$hh.ref$ref.R";
	}
#	}
}

=cut



#`touch "$out.check" ` ;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

###############################################################################
#Scripts usage and about.
# 程序的帮助文档，良好的描述是程序重用和共享的基础，也是程序升级和更新的前提
###############################################################################
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
Usage:    perl $0 
Function: Template for Perl
Command:      
    'verbose' => \$verbose,
    'vcf_list=s' => \$vcf_list,  vcf path or vcf format 
    'o=s' => \$out, 
    'Hmin=i' => \$Hmin,
    'Hmax=i' => \$Hmax,
	'splitchr' => \$splitchr,
	'name=s' =>\$name,
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/8
Notes:    
\n!
    )
}

=cut
for i in 50 60 70
do
for j in 110 120 130
do
echo "
perl /56T/Evan/PMCaller/script/bin/vcf_expansion_H60_120.pl -vcf_list  /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic  -o /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/H${i}_${j} --splitchr -name WGS_depth80x -Hmax $j -Hmin $i 
">/56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/shell/WGS_depth80x_H${i}_${j}.sh
done
done

=
for i in 50 60 70
do
for j in 110 120 130
do
nohup sh /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/merge/chrom_snpout_H130_10pic/expansion/shell/WGS_depth80x_H${i}_${j}.sh &
done
done

vcf不同染色体合并了
cd/56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/mergevcf &&  cat ../merge/chrom_snpout_H130_10pic/*.vcf > mergevcfgrep/wgs3.vcf &
echo "
nohup perl /56T/Evan/PMCaller/script/bin/vcf_expansion_H60_120.pl  -vcf_list /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/mergevcf/ -out  /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/mergevcf/expansionH30_300/ -Hmin 30 -Hmax 300 &
"> /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/mergevcf/shell/wgs3.expansionvcf_table.sh

for i in 30 40 50 60 70 80 90 100 110 120 130 
do 
echo " perl /56T/Evan/PMCaller/script/bin/ski_Expansionvcf_bam_H_L.pl -vcf /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/mergevcf/expansionH30_300/wgs3.exp.vcf  -o /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/mergevcf/expansionH30_300/wgs3Hx_L$i  -L $i &
>> /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/mergevcf/shell/wgs3.expansionvcf_table.sh
done


nohup sh /56T/Evan/PMCaller/ncbi_sra/li/3/BWA2/mergevcf/shell/wgs3.expansionvcf_table.sh & 



