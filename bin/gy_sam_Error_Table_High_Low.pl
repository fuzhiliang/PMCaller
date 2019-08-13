#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin); 
my $verbose ="v1.0";
my $start_time=time;

 #testv1/chr22.snp.gy2
my $in= $ARGV[0];	#genotype 
my $sam= $ARGV[1];	#bam or sam
my $out= $ARGV[2];  
my $l= $ARGV[3];
my $h=$ARGV[4];
my $tqual=$ARGV[5];
$tqual ||=20;
die "perl $0 *snp.gy  /56T/Evan/PMCaller/mydata/DNA.C10T_1.sam|bam  /56T/Evan/PMCaller/mydata/all/errrortable_ski/error. 3:25  30:40:50:60:70  20" if (@ARGV <= 5) ;
# /56T/Evan/PMCaller/mydata/chr22 & perl  /56T/Evan/PMCaller/script/bin/gy_sam_Error_Table_High_Low.pl  chr22.P10T_1.snp.gy  chr22.P10T_1.sorted.bam splittmp 3:10:21 50
my (@l)=(split ":" ,$l);
my (@h)=(split ":" ,$h);
my $maxski =&max(@l);
my $minski=&min(@l);
my $minH =&min(@h);
warn "H:$minH must >=30\n"if ($minH < 30);
warn "max ski is $maxski\n" if ($maxski >30);

#$h ||="50";
#gy:
#输出
#H	L	Genotype        Base    Frequency       total   %F
#30	1	AA      T       16070   60264   0.266660029204832
#		AA      C       20206   60264   0.335291384574539
#		AA      N       47      60264   0.000779901765564848
#		AA      A       5478    60264   0.090900039824771
#		AA      G       18463   60264   0.306368644630293
#		AC      G       15      408     0.0367647058823529
#		AC		A       216     408     0.529411764705882
#输入
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
		if ($tmp[4]==$hh) {
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
#open O ,">$out" or die $! ;
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
	my $qual=$tmp[10];
	my $tag=IsLowQualityReads($qual,$tqual);
	for (my $i=0 ;$i <$#seq_base ;$i++) {
		my $location = $tmp[3]+$i ;
		my $loc_t= "$tmp[2]\t$location";
		my $baselist="";
		foreach my $hh (@h) {
			if (exists $hash{$hh}{$loc_t}) {
				my $ref_genotype= $hash{$hh}{$loc_t};
				my $gt=(split "\t",$ref_genotype)[1];
				$number{$hh}{$loc_t}++;
				my $num=$number{$hh}{$loc_t};
				next if ($num >$maxski) ;
				$error{$tag}{$hh}{$loc_t}{$gt} .=$seq_base[$i];
			}
		}
	}	
	
}
close B;

#my @genotype=("AA","AT","AC","AG","TT","CT","GT","CC","CG","GG");

open O ,">$out/$sample.HighQual.gt" or die $!;
print O "Chr\tLoc\tgenotype\tH\tBases\n";

foreach my $H (keys %{$error{0}}) {
	foreach my $loc_t (keys %{$error{0}{$H}}) {
		foreach my $gt (keys %{$error{0}{$H}{$loc_t}}){
			print O "$loc_t\t$gt\t$H\t$error{0}{$H}{$loc_t}{$gt}\n";
		}
	}
}

close O;

my $cmd="";
$cmd .= "perl $Bin/errortable_v1.pl $out/$sample.HighQual.gt $out $l $sample.HighQual  \n";
&runcmd("step1.0.HighQual",$cmd);

open O1 ,">$out/$sample.LowQual.gt" or die $!;
print O1 "Chr\tLoc\tgenotype\tH\tBases\n";

foreach my $H (keys %{$error{1}}) {
	foreach my $loc_t (keys %{$error{1}{$H}}) {
		foreach my $gt (keys %{$error{1}{$H}{$loc_t}}){
			print O1 "$loc_t\t$gt\t$H\t$error{1}{$H}{$loc_t}{$gt}\n";
		}
	}
}

close O1;

my $cmd1="";
$cmd1 .= "perl $Bin/errortable_v1.pl $out/$sample.LowQual.gt $out $l $sample.LowQual  \n";
&runcmd("step1.0.LowQual",$cmd1);
#perl /56T/Evan/PMCaller/script/bin/errortable_v1.pl chr7.13.1 chr22713/ 3:10 ch
#AAAAAAAAAAAAAAAAAAAAA
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
sub runcmd{
	my $name=shift @_;
	my $cmd=shift @_;
	`mkdir "$out/shell/"` unless (-d "$out/shell");
	open S ,">$out/shell/$name.sh" or die $!;
	print S "$cmd ";
	system "sh $out/shell/$name.sh ";
	close S;
}

sub IsLowQualityReads
{
        my ($qual,$threshold) = @_;
        #print "IsLowQualityReads\n";
        my $tag=0;
        my @chars=split("",$qual);
        my $len= scalar(@chars);
        #print "\$qual:$qual\n\$threshold:$threshold\n\$len:$len\n";
        for(0..$len-1){
                my $value=ord($chars[$_])-33;
                #print "$chars[$_]\t$value\n";
                $tag=1 if($value < $threshold );
        }
        return $tag;
}
