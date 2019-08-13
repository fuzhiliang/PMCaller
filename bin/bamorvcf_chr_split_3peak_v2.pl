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
my ($sam,$bam,$wes,$vcf,$out,$genome,$windows,$stats,$H,$splitchr,$threads);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'sam=s' => \$sam,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
	'stats=s' => \$stats,
    'reference=s' => \$genome,
    'wes' => \$wes,
    'H=i' => \$H,
	'splitchr' => \$splitchr,
    'windows=i' => \$windows,
	'threads=i' => \$threads,
) or die $!;
unless(defined $sam || defined $bam ||defined $vcf){&usage();exit 0;}

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $sam $vcf $bam\nOutput file is $out\n";

$windows||=5 ;
$genome ||= "/data/hg19/hg19_index/hg19.fa";
my $samtools = "/usr/local/bin/samtools";
$out ||= "./";
$threads ||=30;
system "mkdir  -p $out" unless (-d $out);
`mkdir "$out/shell"` unless (-d "$out/shell");
my $name;
my $dir;
my $cmd;
if (-f $vcf && (defined $stats || defined $H)) {
	$name=basename($vcf);
	$dir=abs_path($vcf);
	$name=~s/\.vcf//;	
}elsif (-f $bam) {
	$name=basename($bam);
	$dir=abs_path($bam);
	$name=~s/\.sort(ed)?\.bam//;
	unless (-f $stats){
		$cmd .="$samtools stats --threads $threads $bam > $out/$name.mapping_info &\n"  unless ($H) ;
		$stats="$out/$name.mapping_info";
	}
	$cmd .="$samtools  mpileup -t DP,AD -uvf $genome $bam > $out/$name.vcf \n "; 
	
	$vcf ="$out/$name.vcf";
	open S ,">$out/shell/$name.stats_mpileup.sh";
	print S "$cmd";
	system "sh $out/shell/$name.stats_mpileup.sh ";
}elsif (-f $sam) {
	$name=basename($sam);
	$dir=abs_path($sam);
	$name=~s/\.sam//;
	$cmd .= "$samtools view --threads $threads -Sb $sam |$samtools sort --threads $threads >$out/$name.sorted.bam \n";
	$bam = "$out/$name.sorted.bam";
	unless (-f $stats){
		$cmd .="$samtools stats $bam > $out/$name.mapping_info &\n"  unless ($H) ;
		$stats="$out/$name.mapping_info";
	}
	$cmd .="$samtools  mpileup -t DP,AD -uvf $genome $bam > $out/$name.vcf \n "; 
	
	$vcf ="$out/$name.vcf";
	open S ,">$out/shell/$name.stats_mpileup.sh";
	print S "$cmd";
	system "sh $out/shell/$name.stats_mpileup.sh ";

}else{
die "aa \n";
}
$stats ||="$dir/$name.mapping_info";

if (!defined $H && defined $stats) {
	open S ,"$stats" or die $!;
	my $mapped;
	my $average_length;
	while (<S>) {
		chomp;
		$mapped=(split "\t",$_)[2] if ($_=~m/reads mapped/);
		$average_length=(split "\t" ,$_)[2] if ($_=~m/average length/);
		last;
	}
	close S;
	my $genome_size=3000000000;
	if (defined $wes) {
		$H=$mapped*$average_length*0.015/$genome_size ;
		$H=int($H);
	}else{
		$H=$mapped*$average_length/$genome_size ;
		$H=int($H);
	}
}
die "H$H must >=10 .\n" if ($H <=10);
my @h;
for (my $i=0; $i<5 ;$i++) {
	my $a=$H+($i-2)*$windows;
	push (@h ,$a);
}
#my @h=("$H-2*$windows","$H-$windows","$H","$H+$windows","$H+2*$windows" );
#chr1    10067   .       T       A,G,C   0,255,130,255,133,140,255,133,143,140:    176:  T173,  1,1,1
#chr1    10067   .       taa     ta      0,255,129:123:122,1
#chr1    10068   .       A       C,<*>   0,255,145,255,151,145:178:  A 176,  C2,0


open I ,"$vcf" or die $!;

my %hash;
while (<I>) {
	chomp;
	next if (/^#/ || /^[[<]/ );	
	my @tmp= (split /\t/,$_);
	next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
	my @alt= split ",",$tmp[4];
	next if (length($alt[0]) !=1);
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];
	my (@ad) = (split ",",$ad);
#	next if ($dp != $cut_depth);
	foreach my $hh (@h) {
		if ($dp == $hh) {
			$hash{$hh}{$tmp[3]}{$ad[0]}{$alt[0]}++;
		}
	}
}
my @ref=("A","C","G","T");
foreach my $hh (sort {$a <=>$b } keys %hash) {
	
	foreach my $ref (@ref) {
		my $fig=0;
		open O ,">$out/$name.H$hh.ref$ref.snptable.list" or die $!;
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
		print O "$head\n" if ($fig==0);
		$fig=1;
		print O "$i$alt_t\n";
		}
	}
	close O ;
}
close I;
if (defined $splitchr) {
	my $cmd;
	foreach my $h (@h) {
		foreach my $ref (@ref) {
			$cmd .="perl $Bin/snpchr_merge.pl $out/\\*.H$h.ref$ref.snptable.list $out/merge/ $name \n";
		}
	}
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
	foreach my $hh(@h) {
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
	}
}

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
Usage:    perl $Bin/vcf_split.pl -sam  /56T/Evan/PMCaller/EBI/SRS000116/BWA/ERR1044158.sam windows 10 -H 80 -o /56T/Evan/PMCaller/EBI/SRS000116/BWA/chrom_snpout_H80_10pic ;
Function: Template for Perl
Command:      'verbose' => \$verbose,
    'sam=s' => \$sam,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
	'stats=s' => \$stats,
    'reference=s' => \$genome,
    'wgs' => \$wgs,
    'H=i' => \$H,
	splitchr    is spint chrom
    'windows=i' => \$windows,
	--threads  [30]
Author:   Evan Fu, liuyongxin_bio\@163.com, QQ:42789409
Version:  v1.0
Update:   2018/5/16
Notes:    
\n!
    )
}
