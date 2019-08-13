#! /usr/bin/perl -w
use strict;
use warnings;

my $i =$ARGV[0];
my $out = $ARGV[1];
open I ,"$i";
open O ,"$out";
my %hash;
while (<I>) {
	chomp;
	my ($h,$ref,$db,$total)=split "\t",$_;
	$hash{$h}{$db}{$ref}=$total;
}

my $out ||= "./snpout_filter/";
my $name="SRR345592";

my @ref=("A","C","G","T");
foreach my $hh (sort {$a <=>$b } keys %hash) {
	open O ,">$out/$name.H$hh.snptable.list" or die $!;
	print O "A\tC\tG\tT\n";
	for (my $i=1 ;$i<($hh-2) ;$i++) {
		my $ref_t="";
		foreach my $ref (@ref) {
			if (!exists $hash{$hh}{$i}{$ref}) {
				$hash{$hh}{$i}{$ref}=0;
			}
			$ref_t.="\t$hash{$hh}{$i}{$ref}";
		}
		print O "$i$ref_t\n";
	}
	close O ;
}
close I;
my @h= (30,35,40,45,50);
foreach my $hh(@h) {
	my $Rsh="";
	$Rsh .= "data <- read.table(\"$out/$name.H$hh.snptable.list\",sep='\\t',head=T) \n";
	$Rsh .= "pdf(\"$out/$name.H$hh.pdf\") \n";
	$Rsh .= "matplot(rownames(data),data, type =\"o\", cex=0.8, pch = 16, col = c(\"darkgreen\",\"gold\",\"darkblue\",\"magenta\"),main=\"SNP Table: sample $name H=$hh\",xlab = \"# Reads supporting ref base\", ylab = \"# Position\")\n" ;
	$Rsh .= "legend(\"right\",bty=\"n\",legend=colnames(data),col = c\(\"darkgreen\",\"gold\",\"darkblue\",\"magenta\"),cex=0.8, pch =16) \n";
	$Rsh .= "dev.off()";

	open S ,">$out/shell/$name.H$hh.R";
	print S "$Rsh";
	system "nohup R CMD BATCH  $out/shell/$name.H$hh.R";
}

