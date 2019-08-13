#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
my $verbose ="v1.0";
my $start_time=time;

 #testv1/chr22.snp.gy2
my $in= $ARGV[0];	#genotype 
my $sam= $ARGV[1];	#bam or sam
my $out= $ARGV[2];  
my $l= $ARGV[3];
my $h=$ARGV[4];

die "perl $0 *snp.gy  /56T/Evan/PMCaller/mydata/DNA.C10T_1.sam|bam  /56T/Evan/PMCaller/mydata/all/errrortable_ski/error. 3:25  30:40:50:60:70" if (@ARGV <= 4) ;

my (@l)=(split ":" ,$l);
my (@h)=(split ":" ,$h);
my $maxski =&max(@l);
my $minski=&min(@l);
my $minH =&min(@h);
warn "H:$minH must >=30\n"if ($minH < 30);
warn "max ski is $maxski\n" if ($maxski >30);

#$h ||="50";
#gy:
#ƒø±Í
#H	L	Genotype        Base    Frequency       total   %F
#30	1	AA      T       16070   60264   0.266660029204832
#		AA      C       20206   60264   0.335291384574539
#		AA      N       47      60264   0.000779901765564848
#		AA      A       5478    60264   0.090900039824771
#		AA      G       18463   60264   0.306368644630293
#		AC      G       15      408     0.0367647058823529
#		AC		A       216     408     0.529411764705882
# ‰»Î
#chr1    10037   T       TT	200
#chr1    10052   C       CC 200
#chr1    10061   T       TT
#chr1    10079   T       TT
#chr1    11649   T       TT
#chr1    11688   T       TT
#chr1    11698   G       GG
#chr1    11726   T       TT
#chr1    11757   C       CC
#chr1    11763   C       CC


print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $sample=basename($sam);
$sample =~s/\.[bs]am$//;
my %hash;

open I ,"$in" or die $!;
while (<I>) {
	chomp;
	my @tmp= split "\t", $_;
	my $genotype=$tmp[3];
	my $ref= $tmp[2];

	my $tmp_t= join "\t",($tmp[0],$tmp[1]);
	foreach my $hh (@h) {
		if ($tmp[4]>=$hh) {
			$hash{$hh}{$tmp_t}="$ref\t$genotype";
		}
	}
}

close I;
if ($sam=~/sam$/) {
	open B ,"$sam" or die $!;
}else{
	open B ,"samtools view $sam |" or die $!;
}

my @base=("A","C","G","T","N");

#open B ,"$sam" or die $!;
print "open $sam file.\n";
my $sam_base="";
my %number;
#my %hash2;
my %error;
my %error_number;
while (<B>) {
	chomp;
	next if (/^@/);
	my @tmp= split "\t", $_;
	next if ($tmp[5]!~m/^([0-9]{2,3})M$/);
#	my $length=$1;
	my @seq_base=split "", $tmp[9];

	for (my $i=0 ;$i <$#seq_base ;$i++) {
		my $location = $tmp[3]+$i ;
		my $loc_t= "$tmp[2]\t$location";
		foreach my $hh (@h) {
			if (exists $hash{$hh}{$loc_t}) {
				my $ref_genotype= $hash{$hh}{$loc_t};
				my $gt=(split "\t",$ref_genotype)[1];
				$number{$hh}{$loc_t}++;
				my $num=$number{$hh}{$loc_t};
#				next if ($num < $minski);
#				next if ($num > $maxski);
#				$error{$hh}{$num}{$gt}{$seq_base[$i]}++;
#				$error_number{$hh}{$num}{$gt}++;

				foreach my $ll (@l) {
					if ($num <= $ll ) {
						$error{$hh}{$ll}{$gt}{$seq_base[$i]}++;
						$error_number{$hh}{$ll}{$gt}++;
					}
					foreach my $bbbb (@base) {
						if ($ll == $error{$hh}{$ll}{$gt}{$bbbb}) {
							$error{$hh}{$ll}{$gt}{$bbbb}-=$ll;
							$error_number{$hh}{$ll}{$gt}-=$ll;
						}
					}
				}

			}
		}
	}	
	
}
close B;
open O ,">$out" or die $!;
print O "H\tL\tgenotype\tbase\tfre\ttotal\n";
my @base=("A","C","G","T","N");
my @genotype=("AA","AT","AC","AG","TT","CT","GT","CC","CG","GG");
#my $ski_head=join "\t",(sort(@ski));
#my $head="Chr\tLocation\tref\tgenotype\t$ski_head\n" ;
#print O $head;
foreach my $H (keys %error) {
	foreach my $L (keys %{$error{$H}}) {
#		foreach my $gt (keys %{$error{$H}{$L}}) {
		foreach my $gt (@genotype ) {
			if (!exists $error_number{$H}{$L}{$gt}){
				$error_number{$H}{$L}{$gt}=0;
			}
#			foreach my $base (keys %{$error{$H}{$L}{$gt}}) {
			foreach my $base (@base) {
				if (exists $error{$H}{$L}{$gt}{$base}) {
					print O "$H\t$L\t$gt\t$base\t$error{$H}{$L}{$gt}{$base}\t$error_number{$H}{$L}{$gt}\n";
				}else {
					print O "$H\t$L\t$gt\t$base\t0\t$error_number{$H}{$L}{$gt}\n";
				}
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
sub min {
    my $currentMaxValue = shift @_;
    foreach ( @_ ) {
        if ( $_ < $currentMaxValue ) {
            $currentMaxValue = $_;
        }
    }
    return $currentMaxValue;
}