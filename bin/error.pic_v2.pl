#! usr/bin/perl -w
use strict;
use warnings;
use FileHandle;
my $in=$ARGV[0];  #/56T/Evan/PMCaller/mydata/WES_error_table/  目录  ${sample}.eout.*.error.table.list
#my $sample=$ARGV[1];
my $out=$ARGV[1];   #输出目录

die unless (-d $in);
system "mkdir -p $out" unless (-d $out);

my %pic;
#DNA.P9T.30.eout.21.error.table.list
my @list=glob ("$in/*.error.table.list");
#print @list;
open W ,">$out/low_Frequency.list";
print W "L\tGenotype\tBase\tFrequency\ttotal\t%F\n";
foreach my $list ( @list ){
#	print "$list\n";
	$list=~m/.*\.(\d{1,2})\..*/;
	my $ski=$1;
#	print "L=$ski\n";
	open I ,"$list" or dir $!;
	while (<I>) {
		chomp;
		next if (/^Genotype/);
		next if (/^genotype/);
		my ($genotype,$base,$F,$total,$frq)=(split "\t",$_)[0,1,2,3,4];
		if ($F<=10) {
			my ($a,$b)=split("",$genotype);
			if($a eq $base || $b eq $base){
				print "Error, L=$base,genotype=$genotype,base=$base,frequency=$frq.\n";
			}
			print W  "Error, L=$base,genotype=$genotype,base=$base,frequency=$frq.\n";
		}
		$pic{$genotype}{$ski}{$base}=$frq;
	}
}
open O ,">$out/error_pic.table" or die;
print O "Genotype\tL\tA\tC\tG\tT\n";
#my @bases =("A","C","G","T");
my @bases=("A","C","G","T");
#print  @bases;
foreach my $gt (sort{$a cmp $b} keys %pic) {
	foreach my $L ( keys %{$pic{$gt}}) {
		my $f_list="";
		foreach my $b (@bases) {

			if (exists  $pic{$gt}{$L}{$b}) {
				$f_list .="\t$pic{$gt}{$L}{$b}";
			}else{
				$f_list .="\t0";
			}
		}
		print O "$gt\t$L$f_list\n"; 
	}
}
close I; 
close O;

