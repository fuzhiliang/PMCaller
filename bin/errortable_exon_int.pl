#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
my $verbose ="v1.0";
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
 #testv1/chr22.snp.gy2
my $gt= $ARGV[0];   #碱基连在一起的，没有根据L分多列
my $out= $ARGV[1];
my $L=$ARGV[2];  #削薄后的L  用 ：号分割
my $sample=$ARGV[3];
my $gtf=$ARGV[4];   #/data/hg19_anno/Region/refgene.nooverlap 基因内区间 #chr1    11873   29370

#输入 chr22   24324590        TT      50      TTTTTTTTTTTTTTGTTTTTT    01111111111111110000
#输出 添加注释，是exon or int
#添加High low qual 生成两个table
die "perl $0 /56T/Evan/PMCaller/mydata/chr22/chr7.13.1.gt  /56T/Evan/PMCaller/mydata/chr22/chr22713/ 3:10 chr22\n" if (@ARGV<4);
#输出 exon  int 的两个error table
`mkdir $out` unless (-d $out);
my @ll= split ":",$L;
my $lmax=&max(@ll);
warn "max L is $lmax \n";

$gtf ||= "/data/hg19_anno/Region/refgene.nooverlap";
my %gtf;
open G ,"$gtf" or die $!;  #chr1    11873   29370
while (<G>) {
	chomp;
	my ($chr,$start,$end)=(split "\t",$_);
	$gtf{$chr}{$start}=$end;
}


my %hash_error;

my %count_error;
open I ,"$gt" or die $!;
open EI ,">$gt.exonint" or die $!;
while (<I>) {
	chomp;
	next if (/^H/);
	my $gene="INT";
	my ($chr,$loc)= (split "\t", $_)[0,1];
	my @tmp = split "\t", $_;
	if (exists $gtf{$chr} ) {
		foreach my $start (keys %{$gtf{$chr}}) {
			if ($loc >= $start && $loc <= $gtf{$chr}{$start}){
				$gene="EXON";		
				last;
			}
		}
	}
	print EI "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$gene\n";

	my $bb=(split "", $tmp[4])[0];
	foreach my $l (@ll) {
		next if (length($tmp[4]) < $l);
		
		my $baselist=substr($tmp[4],0,$l); #取前$l个base
		#my $quallist=substr($tmp[5],0,$l); #取前$l个质量值
		my $check_base="$bb"x$l ;
#		print "$baselist\t$check_base\n"; die ;
		next if ($baselist eq $check_base); #判断。
		my @base = split "" , $baselist;
#		my @qual = split "" , $quallist;
		#my $hqcount=$quallist=~ tr/0/0/;
		#my $lqcount=$quallist=~ tr/1/1/;
		for (my $i=0;$i<=$#base;$i++) {
		#	$hash_error{$tmp[3]}{$l}{$tmp[2]}{"High"}{$base[$i]}++  if ($qual[$i]==0);
		#	$hash_error{$tmp[3]}{$l}{$tmp[2]}{"Low"}{$base[$i]}++  if ($qual[$i]==1);
			$hash_error{$tmp[3]}{$l}{$tmp[2]}{$gene}{$base[$i]}++ ;
			$hash_error{$tmp[3]}{$l}{$tmp[2]}{$base[$i]}++;
			$count_error{$tmp[3]}{$l}{$tmp[2]}++;
		}
	}
}
	
close I;
my @base=("A","C","G","T");
my @genotype=("AA","AT","AC","AG","TT","CT","GT","CC","CG","GG");
open O ,">$out/$sample.errortable.list" or die $!;
open OE ,">$out/$sample.EXON.errortable.list" or die $!;
open OI ,">$out/$sample.INT.errortable.list" or die $!;

print O  "H\tL\tGenotype\tBase\tFrequency\ttotal\t%F\n";
foreach my $H (sort {$a cmp $b} keys %hash_error) {
	foreach my $L ( keys %{$hash_error{$H}}) {
		foreach my $gt (keys %{$hash_error{$H}{$L}}) {
			foreach my $base (@base) {
#				$hash_error{$H}{$L}{$gt}{"Low"}{$base}=0  if (!exists $hash_error{$H}{$L}{$gt}{"Low"}{$base}) ;
#				$hash_error{$H}{$L}{$gt}{"High"}{$base}=0  if (!exists $hash_error{$H}{$L}{$gt}{"High"}{$base}) ;
				$hash_error{$H}{$L}{$gt}{"EXON"}{$base}=0  if (!exists $hash_error{$H}{$L}{$gt}{"EXON"}{$base}) ;
				$hash_error{$H}{$L}{$gt}{"INT"}{$base}=0  if (!exists $hash_error{$H}{$L}{$gt}{"INT"}{$base}) ;
				$hash_error{$H}{$L}{$gt}{$base}=0  if (!exists $hash_error{$H}{$L}{$gt}{$base}) ;
				my $p_E= $hash_error{$H}{$L}{$gt}{"EXON"}{$base}/$count_error{$H}{$L}{$gt};
				my $p_I= $hash_error{$H}{$L}{$gt}{"INT"}{$base}/$count_error{$H}{$L}{$gt};
				my $p= $hash_error{$H}{$L}{$gt}{$base}/$count_error{$H}{$L}{$gt};
				
				print OE "$H\t$L\t$gt\t$base\t$hash_error{$H}{$L}{$gt}{'EXON'}{$base}\t$count_error{$H}{$L}{$gt}\t$p_E\n";
				print OI "$H\t$L\t$gt\t$base\t$hash_error{$H}{$L}{$gt}{'INT'}{$base}\t$count_error{$H}{$L}{$gt}\t$p_I\n";
				print O "$H\t$L\t$gt\t$base\t$hash_error{$H}{$L}{$gt}{$base}\t$count_error{$H}{$L}{$gt}\t$p\n";
			}
		}
	}
}
close OE;
close OI;
close O;


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
