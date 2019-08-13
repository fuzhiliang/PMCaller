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


my $i=$ARGV[0];
my $out=$ARGV[1];
my $name=$ARGV[2];
die "die$0\n" if (@ARGV<2);
$i=abs_path$i;
my $dir=`pwd`;
chomp ($dir);
$out="$dir/$out"; #获取绝对路径
print "outputpath:\t$out\n";
`mkdir $out` unless (-d $out);

my @chrlist=(glob "$i");
#print "  @chrlist\n";
my %hash;
my ($chr,$h,$ref,$head);
foreach my $names (@chrlist) {
	my $sample=basename($names);
	$sample=~m/(.*)\..*chr([0-9XY]{1,2})\..*H([0-9]{1,3})\.ref([ATCG])\..*/;
	$sample=~m/chr([0-9XY]{1,2})\..*H([0-9]{1,3})\.ref([ATCG])\..*/;

	$name||="a";
	$h=$2;
	$ref=$3;
#	print "$name\n$h\n$ref\n";
	open I ,"$names";
	$head=<I>;
	while (<I>) {
		chomp ;
        my @l=split "\t",$_;	
		for (my $j=1;$j<=3 ;$j++) {
			$hash{$l[0]}{$j}=0  if (!exists $hash{$l[0]}{$j}) ;
			my $suma=$l[$j]+$hash{$l[0]}{$j};
			$hash{$l[0]}{$j}=$suma;
		}
	}
}

open O ,">$out/$name.H$h.ref$ref.txt" or die $!;
print O "$head";
foreach my $l (sort {$a<=>$b} keys %hash) {
	my $col="";
	for (my $i=1;$i<=3 ;$i++) {
		$col .="\t$hash{$l}{$i}";
	}
	print O "$l$col\n";
}

close I ;
close O ;

#my @h =(110 , 120 , 130 , 140 , 150);
#my @ref=("A","C","G","T");
`mkdir "$out/shell"` unless (-d "$out/shell");
my $hh=$h;
#foreach my $hh(@h) {
#	foreach my $ref (@ref) {
		my $Rsh="";
		my $h1=$hh-2 ;
		$Rsh .= "data <- read.table(\"$out/$name.H$hh.ref$ref.txt\",sep='\\t',head=T) \n";
		$Rsh .= "data.m <- data[2:$h1,] \n";
		$Rsh .= "pdf(\"$out/$name.H$hh.ref$ref.pdf\") \n";
		$Rsh .= "matplot(rownames(data.m),data.m, type =\"o\", cex=0.8, pch = 16, col = c(\"darkgreen\",\"gold\",\"darkblue\",\"magenta\"),main=\"SNP Table: sample $name H=$hh REF=$ref\",xlab = \"# Reads supporting ref base\", ylab = \"# Position\")\n" ;
		$Rsh .= "legend(\"right\",bty=\"n\",legend=colnames(data),col = c\(\"darkgreen\",\"gold\",\"darkblue\",\"magenta\"),cex=0.8, pch =16) \n";
		$Rsh .= "dev.off()";
		open S ,">$out/shell/$name.H$hh.ref$ref.R";
		print S "$Rsh";
		close S ;
		system "cd $out/shell && nohup R CMD BATCH  $out/shell/$name.H$hh.ref$ref.R";
#	}
#}



#rm merge.sh
#for  h in  110 120 130 140 150 
#do
#for  b in A C G T 
#do
#a=SRR345592.sam.chr\\*.sam.H$h.ref$b.snptable.list
#echo "perl  /56T/Evan/PMCaller/script/bin/snpchr_merge.pl  /56T/Evan/PMCaller/ncbi_sra/li/4/BWA/SRR345592_chrom_snpout_H130_10pic/$a /56T/Evan/PMCaller/ncbi_sra/li/4/BWA/SRR345592_chrom_merge " >>merge.sh
#done
#done

#perl  /56T/Evan/PMCaller/script/bin/snpchr_merge.pl  SRR345592_chrom_snpout_H130_10pic/SRR345592.sam.chr\*.sam.H150.refT.snptable.list SRR345592_chrom_snpout_H130_10pic_merge