#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

#输入：
#/8T4/Evan/PMCaller/WES_notQC_filter/PRJNA59817/SRR107087/BWA_notQC2/down5/SRR107087.filter.filter.sorted.DP0.1.HighQual.qual20.vcf.modify:chr1  117618500   T       NN      3       3       0
#/8T4/Evan/PMCaller/WES_notQC_filter/PRJNA59817/SRR107087/BWA_notQC2/down5/SRR107087.filter.filter.sorted.DP0.1.LowQual.qual20.vcf.modify:chr1   117618500   T       G       1       0       1
#输出:
#chr1   117618500   T       G	4       3       0	0	1

my $High_modify = $ARGV[0];
my $Low_modify = $ARGV[1];
my $out = $ARGV[2];

die "perl $0 /8T4/Evan/PMCaller/WES_notQC_filter/PRJNA59817/SRR107087/BWA_notQC2/down5/SRR107087.filter.filter.sorted.DP0.1.HighQual.qual20.vcf.modify SRR107087.filter.filter.sorted.DP0.1.LowQual.qual20.vcf.modify  SRR107087.filter.filter.sorted.DP0.1.HighQual.LowQual.qual20.dubious.merge " if (@ARGV<3);
my %High;
open H ,"$High_modify" or die $!;
while (<H>) {
	chomp ;
	my (@tmp)=split "\t",$_ ;
	my $loc = "$tmp[0]\t$tmp[1]\t$tmp[2]";
	#my $v_list = "$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]";
	$High{$loc}=$_;
}
open O ,">$out " or die $!;
open L ,"$Low_modify" or die $!;
while (<L>) {
	chomp ;
	my (@tmp)=split "\t",$_ ;
	my $loc = "$tmp[0]\t$tmp[1]\t$tmp[2]";
	if (exists $High{$loc}) {
		my $v_list=$High{$loc};
		my @v_list=split "\t",$v_list;
		next if ($tmp[5]+$v_list[5]==0);
		next if ($tmp[6]+$v_list[6]==0);
		next if ($tmp[3] eq "NN" && $v_list[3] eq "NN");
		my $alt=$tmp[3];
		$alt =  $v_list[3] unless ( $v_list[3] eq "NN");
		my $depth=$tmp[4] + $v_list[4];
		print O "$loc\t$alt\t$depth\t$v_list[5]\t$v_list[6]\t$tmp[5]\t$tmp[6]\n";
	}
}
close H;
close L;
close O;
