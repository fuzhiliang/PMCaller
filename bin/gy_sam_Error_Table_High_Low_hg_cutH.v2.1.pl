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
my $head_low_base=$ARGV[6]; #将每个reads中的前后n bp视为低质量base。
$tqual ||=20;
$head_low_base||=0;
die "perl $0 *snp.gy  /56T/Evan/PMCaller/mydata/DNA.C10T_1.bam  /56T/Evan/PMCaller/mydata/all/errrortable_ski/error  3:25  50  20  0 \n"  if (@ARGV <= 5) ;  # 将第一bp和最后一bp视为low quality
# /56T/Evan/PMCaller/mydata/chr22 & perl  /56T/Evan/PMCaller/script/bin/gy_sam_Error_Table_High_Low.pl  chr22.P10T_1.snp.gy  chr22.P10T_1.sorted.bam splittmp 3:10:21 50
my (@l)=(split ":" ,$l);
my (@h)=(split ":" ,$h);
my $maxski =&max(@l);
my $minski=&min(@l);
my $minH =&min(@h);
warn "H:$minH must >=30\n"if ($minH < 30);
warn "max ski is $maxski\n" if ($maxski >30);

`mkdir -p $out ` unless (-d $out);
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
$sample=~s/(\.sort(ed)?)?\.bam$//;
my %hash;

open I ,"$in" or die $!;
while (<I>) {
	chomp;
	my @tmp= split "\t", $_;
	my $genotype=$tmp[3];
	my $ref= $tmp[2];

	my $tmp_t= join "\t",($tmp[0],$tmp[1]);
	foreach my $hh (@h) {
#		if ($tmp[4]==$hh) {
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
#open O ,">$out" or die $! ;
#open B ,"$sam" or die $!;
print "open $sam file.\n";
my $sam_base="";
my %number;
#my %hash2;
my %error;
my %Qual;
while (<B>) {
	chomp;
	next if (/^@/);
	my @tmp= split "\t", $_;
	next if ($tmp[5]!~m/^([0-9]{2,3})M$/);
#	my $length=$1;
	my @seq_base=split "", $tmp[9];
	my $qual=$tmp[10];
	my $tag=IsLowQualityReads($qual,$tqual);
	for (my $i=0 ;$i <=$#seq_base ;$i++) {
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
				#$error{$tag}{$hh}{$loc_t}{$gt} .=$seq_base[$i];
				$error{$hh}{$loc_t}{$gt} .=$seq_base[$i];
#				print "$i\t$#seq_base\n";die;
				if ($i<$head_low_base || $i > ($#seq_base-$head_low_base)) {
					$Qual{$hh}{$loc_t}{$gt} .=1;	
				}else{
					$Qual{$hh}{$loc_t}{$gt} .=$tag;  #0=high  1=low  
				}
			}
		}
	}		
}
close B;

open O ,">$out/$sample.gt" or die $!;
print O "H\tL\tgenotype\tbase\tfre\ttotal\n";

foreach my $H (keys %error) {
	foreach my $loc_t (keys %{$error{$H}}) {
		foreach my $gt (keys %{$error{$H}{$loc_t}}){
			print O "$loc_t\t$gt\t$H\t$error{$H}{$loc_t}{$gt}\t$Qual{$H}{$loc_t}{$gt}\n";
		}
	}
}

close O;
my $cmd="";
#$cmd .= "perl $Bin/errortable_v1.pl $out/$sample.gt $out 3:11:25 $sample  \n";
$cmd .= "perl $Bin/errortable_highlow.pl $out/$sample.gt $out $l $sample  \n";  #gt 文件添加了一列质量值

&runcmd("step1.0",$cmd);


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

#修改gt文件，添加high low 标志