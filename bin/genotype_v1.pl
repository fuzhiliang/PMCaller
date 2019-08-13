#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
my $verbose ="v1.0";

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
#samtools mpileup -t DP,AD -uvf /data/hg19/hg19.fa chr22.P10T.ski.0.2.sort.bam > chr22.P10T.ski.0.2.vcf #准备文件
#toli@s3:/56T/Evan/PMCaller/mydata/depth_test$ ps fx |less
#nohup perl ../../script/bin/genotype_v1.pl chr22.P10T.ski.0.2.vcf /56T/Evan/PMCaller/mydata/testv1/table_list.snp.table /56T/Evan/PMCaller/mydata/testv1/table_list.error.table ./chr22_genotype_out


my $ski_bam_vcf = $ARGV[0];
my $snptable=$ARGV[1];
my $errortable=$ARGV[2];
my $out =$ARGV[3];

my %snptable;
my %errortable;
open S ,"$snptable" or die $!;
while (<S>) {
	chomp;
	next if (/^Ref/);
	my ($ref,$snp_genotype,$prior_probability)=(split "\t",$_)[0,1,4];
	$snptable{$ref}{$snp_genotype}=$prior_probability;
}

open E ,"$errortable" or die $!;
while (<E>) {
	chomp;
	next if (/^Genotype/);
	my ($error_genotype,$base,$contingent_probability)=(split "\t",$_)[0,1,4];
	$errortable{$base}{$error_genotype}=$contingent_probability;
}

open I ,"$ski_bam_vcf" or die $!;
open O ,">$out" or die $!;
my %hash;
while (<I>) {
	chomp;
	next if (/^#/ || /^[[<]/ );	
	my @tmp= (split /\t/,$_);
	next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
	$tmp[4]=~s/,<\*>//;
	my @alt= split ",",$tmp[4];

	next if (length($alt[0]) !=1);
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];
######	C       T,G,<*> 0,103,6,255,206,6,255,235,12,6  530:377,151,2,0
#	next if ($dp < $cut_depth); 	
	my @ad = split "," ,$ad;
	next if ($ad[1] == 0 );
	my @gg=("$tmp[3]","$alt[0]");
	@gg = sort @gg;
	my $gg= join "",@gg;
	my @genotype =("$tmp[3]$tmp[3]","$gg","$alt[0]$alt[0]");
#	print @genotype;
	my $genotype_p0=0;
	my $target_genotype;
	foreach my $genotype (@genotype) {
		my $snp_p= $snptable{$tmp[3]}{$genotype};
		my $e1=$errortable{$tmp[3]}{$genotype}**($ad[0]);
		my $e2=$errortable{$alt[0]}{$genotype}**($ad[1]);
		my $genotype_p = $snptable{$tmp[3]}{$genotype} * ($errortable{$tmp[3]}{$genotype}**($ad[0])) * ($errortable{$alt[0]}{$genotype} ** $ad[1]); 
#		print  "test\t$snp_p\t$e1\t$e2\t$genotype_p\n$snptable{$tmp[3]}{$genotype} * ($errortable{$tmp[3]}{$genotype}**($ad[0])) * ($errortable{$alt[0]}{$genotype}**$ad[1])\n";

		if ($genotype_p > $genotype_p0) {
			$target_genotype=$genotype;
			$genotype_p0=$genotype_p;
		}
	}

	print O "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t$dp\t$ad\t$target_genotype\n";

}
close I;
close O;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";



