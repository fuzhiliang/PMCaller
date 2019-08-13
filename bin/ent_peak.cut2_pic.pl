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

my $indir=$ARGV[0];
my $out=$ARGV[1];


die "$0 indir outdir \n" if (@ARGV != 2 );
`mkdir $out ` unless (-d $out);
`mkdir "$out/shell" ` unless (-d "$out/shell");
$out=abs_path$out;
$indir=abs_path$indir;
print "$out\n$indir\n";
my @sample =(glob "$indir/*snptable.list");
my @H;

my $name;
foreach my $sample (@sample) {
	$sample=~m/.*\/(.*)\.H([0-9]{1,3})\..*/;
	$name=$1;
	my $h=$2;
	push (@H , $h);
}

#my @h= (30,35,40,45,50);
my @ref=("A","C","G","T");
foreach my $hh(@H) {
	foreach my $ref (@ref) {
		my $Rsh="";
		my $h1=$hh-2 ;
		$Rsh .= "data <- read.table(\"$indir/$name.H$hh.ref$ref.snptable.list\",sep='\\t',head=T) \n";
		$Rsh .= "data.m <- data[2:$h1,] \n";
		$Rsh .= "pdf(\"$out/$name.H$hh.ref$ref.pdf\") \n";
		$Rsh .= "matplot(rownames(data.m),data.m, type =\"o\", cex=0.8, pch = 16, col = c(\"darkgreen\",\"gold\",\"darkblue\",\"magenta\"),main=\"SNP Table: sample $name H=$hh REF=$ref\",xlab = \"# Reads supporting ref base\", ylab = \"# Position\")\n" ;
		$Rsh .= "legend(\"right\",bty=\"n\",legend=colnames(data),col = c\(\"darkgreen\",\"gold\",\"darkblue\",\"magenta\"),cex=0.8, pch =16) \n";
		$Rsh .= "dev.off()";
		open S ,">$out/shell/$name.H$hh.ref$ref.R";
		print S "$Rsh";
		close S ;
		system "cd $out/shell && nohup R CMD BATCH  $out/shell/$name.H$hh.ref$ref.R";
	}
}



