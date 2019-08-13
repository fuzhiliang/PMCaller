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
my $in= $ARGV[0];   #碱基连在一起的，没有根据L分多列
my $out= $ARGV[1];
my $L=$ARGV[2];  #削薄后的L  用 ：号分割
my $sample=$ARGV[3];
#输入 chr22   24324590        TT      50      TTTTTTTTTTTTTTGTTTTTT
die "perl $0 /56T/Evan/PMCaller/mydata/chr22/chr7.13.1 /56T/Evan/PMCaller/mydata/chr22/chr22713/ 3:10 chr22\n" if (@ARGV<4);
#输出 H	L
`mkdir $out` unless (-d $out);
my @ll= split ":",$L;
my $lmax=&max(@ll);
warn "max L is $lmax \n";
my %hash_error;
my %count_error;
open I ,"$in" or die $!;
while (<I>) {
	chomp;
	next if (/^H/);
	my @tmp= split "\t", $_;
	my $bb=(split "", $tmp[4])[0];
	foreach my $l (@ll) {
		next if (length($tmp[4]) < $l);
		my $baselist=substr($tmp[4],0,$l);
		
		my $check_base="$bb"x$l ;
#		print "$baselist\t$check_base\n"; die ;
		next if ($baselist eq $check_base);
		my @base = split "", $baselist;
		foreach my $base (@base) {
			$hash_error{$tmp[3]}{$l}{$tmp[2]}{$base}++;
			$count_error{$tmp[3]}{$l}{$tmp[2]}++;
		}
	}
}
	
close I;
my @base=("A","C","G","T");
my @genotype=("AA","AT","AC","AG","TT","CT","GT","CC","CG","GG");
open O ,">$out/$sample.errortable.list" or die $!;
print O  "H\tL\tGenotype\tBase\tFrequency\ttotal\t%F\n";
foreach my $H (sort {$a cmp $b} keys %hash_error) {
	foreach my $L ( keys %{$hash_error{$H}}) {
		foreach my $gt (keys %{$hash_error{$H}{$L}}) {
			foreach my $base (@base) {
				$hash_error{$H}{$L}{$gt}{$base}=0  if (!exists $hash_error{$H}{$L}{$gt}{$base}) ;
				my $p= $hash_error{$H}{$L}{$gt}{$base}/$count_error{$H}{$L}{$gt};
				print O "$H\t$L\t$gt\t$base\t$hash_error{$H}{$L}{$gt}{$base}\t$count_error{$H}{$L}{$gt}\t$p\n";
			}
		}
	}
}

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
