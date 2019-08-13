#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
my $verbose ="v1.0";
my $start_time=time;

 #testv1/chr22.snp.gy2
my $in= $ARGV[0];
my $sam= $ARGV[1];
my $out= $ARGV[2];
my $ski= $ARGV[3];
my $h=$ARGV[4];
#$ski ||= 10;
my @ski=split ":" ,$ski;
my $maxski =&max(@ski);
$h ||="50";
die "H:$h must >=30\n";
warn "max ski is $maxski\n" if ($maxski >20);
die "perl $0 /56T/Evan/PMCaller/mydata/all/DNA.C10T_1.snp.gy2  /56T/Evan/PMCaller/mydata/DNA.C10T_1.sam|bam  /56T/Evan/PMCaller/mydata/all/errrortable_ski/DNA.C10T.ski 5:7:9:11:13:15  40" if (@ARGV <= 4) ;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my %hash;
open I ,"$in" or die $!;
while (<I>) {
	chomp;
	my @tmp= split "\t", $_;
	my $genotype=$tmp[3];
	my $ref= $tmp[2];
	my $tmp_t= join "\t",($tmp[0],$tmp[1]);
	next if ($tmp[4] < $h );
	$hash{$tmp_t}="$ref\t$genotype";
}

close I;
if ($sam=~/sam$/) {
	open B ,"$sam" or die $!;
}else{
	open B ,"samtools view $sam |" or die $!;
}
#open B ,"$sam" or die $!;
my $sam_base="";
my %number;
my %hash2;
while (<B>) {
	chomp;
	next if (/^@/);
	my @tmp= split "\t", $_;
	if ($tmp[5]=~m/([0-9]{2,3})M/){
		my $length=$1;
		my @seq_base=split "", $tmp[9];
		for (my $i=0 ;$i <$#seq_base ;$i++) {
			my $location = $tmp[3]+$i ;
			my $loc_t= "$tmp[2]\t$location";
			if (exists $hash{$loc_t}) {
				$number{$loc_t}++;
#				next if ($number{$loc_t} >$ski);
				my $value= "$loc_t\t$hash{$loc_t}";
				if (exists $hash2{$value}) {
					$sam_base =$hash2{$value};
					$sam_base .=$seq_base[$i];
					$hash2{$value}=$sam_base;
				}else{
					$hash2{$value}=$seq_base[$i]; 
				}
			}
		
		}	
	}
}
open O ,">$out" or die $!;
my $ski_head=join "\t",(sort(@ski));
my $head="Chr\tLocation\tref\tgenotype\t$ski_head\n" ;
print O $head;
foreach my $keys (keys %hash2) {
	my $a="";
	foreach my $s (sort @ski) {
		my $seq =substr ($hash2{$keys}, 0, $s);
		$a .="\t$seq";
	}
	print O "$keys$a\n";
}

close B;
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