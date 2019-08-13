#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use POSIX qw(strftime);
#use Getopt::Std;
use File::Basename;
use File::Glob;
use FindBin qw($Bin); 
use Cwd qw(abs_path getcwd);  
my $abs = abs_path(getcwd());  
my $verbose ="v1.0";
###############################################################################
#命令行参数据的定义和获取，记录程序初始时间，设置参数默认值
#Get the parameter and provide the usage.
###cd /56T/Evan/PMCaller/check_data && perl /56T/Evan/PMCaller/script/PMCaller_v1.pl -sam  chr22_head.P10T_1.sam -o test -H 40:50 -L 5:8:9:20
###############################################################################
my ($snptable,$errortable,$out,$add_error,$matrix,$snp_rand,$error_rand);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'o=s' => \$out, 
    'add_error_rate=s' => \$add_error,
	'snp_rand=s' => \$snp_rand,
	'error_rand=s' => \$error_rand,
    'snptable=s' => \$snptable,
    'errortable=s' => \$errortable,
	'matrix' => \$matrix,
) or die $!;
#unless(defined $vcf || defined $bam ){&usage();exit 0;}
unless(defined $snptable || defined $errortable) {&usage();exit 0;}
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $errortable\nOutput file is $out\n";
#修改error_rand文件，只要纯合的四列

$add_error||=0.6;
$out||='./';

system "mkdir $out " unless (-d $out);

#snptalbe
if ($snp_rand && $snptable) {
	open A ,"$snp_rand" or die $!;
	my $head=<A>;
	chomp $head;
	my @head=split "\t",$head;
	#print @head ;
	my @group;
	while (<A>) {
		chomp;
		next if (/A/);
		push(@group,$_);
	}
	for (my $j=0;$j<=$#group ;$j++) {
		my @tmp=(split "\t",$group[$j]);
		
		my %snp;
		open E ,"$snptable" or die $!;
		while (<E>) {
			chomp;
			next if (/^Ref/i);
			my ($genotype,$base,$fre,$total)=(split "\t",$_)[0,1,2,3];
			$snp{$genotype}{$total}{$base}=$fre;
		}
		close E;


		for (my $i=0;$i<=$#tmp;$i++) {
			print "$i\n";
			foreach my $total (keys %{$snp{$head[$i]}}) {
				my $min=$snp{$head[$i]}{$total}{AA};
				my $max=$snp{$head[$i]}{$total}{AA};
				my $minbase="AA";				                  #
				my $maxbase="AA";
				foreach my $base ( keys %{$snp{$head[$i]}{$total}}) { #确定最小的，最大的base是什么，对应的值是什么
					next if ($snp{$head[$i]}{$total}{$base}==0);
	#				print "$min\t$snp{$head[$i]}{$total}{$base}\n";
					if ($snp{$head[$i]}{$total}{$base}<$min) {
						$min=$snp{$head[$i]}{$total}{$base};
						$minbase=$base;
					}
					if ($snp{$head[$i]}{$total}{$base}>$max) {
						$max=$snp{$head[$i]}{$total}{$base};
						$maxbase=$base;
					}
				}
				my $next_min=$max;  #
				my $next_minbase=$maxbase;  
				foreach my $base ( keys %{$snp{$head[$i]}{$total}}) { #确定第二小的base是什么，对应的值是什么
					next if ($snp{$head[$i]}{$total}{$base}==0);
	#				print "$min\t$snp{$head[$i]}{$total}{$base}\n";
					if ($snp{$head[$i]}{$total}{$base}>$min &&  $snp{$head[$i]}{$total}{$base}<=$next_min) {
						$next_min=$snp{$head[$i]}{$total}{$base};
						$next_minbase=$base;
					#	print "genotype\t$head[$i]\tref\t$base\n"
					}
				}
	#			print "minbase\t$minbase\n";
				if ($tmp[$i]==1){
					my $bb=int($min+$add_error*$min);
					my $bb2=int($next_min+$add_error*$next_min);
					my $cc=int($max-$add_error*$min-$add_error*$next_min);
					$snp{$head[$i]}{$total}{$minbase}=$bb;
					$snp{$head[$i]}{$total}{$next_minbase}=$bb2;
					$snp{$head[$i]}{$total}{$maxbase}=$cc;
				}else{
					my $bb=int($min-$add_error*$min);
					my $bb2=int($next_min-$add_error*$next_min);
					my $cc=int($max+$add_error*$min+$add_error*$next_min);
					$snp{$head[$i]}{$total}{$minbase}=$bb;
					$snp{$head[$i]}{$total}{$next_minbase}=$bb2;
					$snp{$head[$i]}{$total}{$maxbase}=$cc;
				}

			}
		}
		#print Dumper %snp; 
		open OE ,">$out/$j.snptable" or die $!;
		foreach my $genotype (sort {$a cmp $b} keys %snp) {
			foreach my $total (keys %{$snp{$genotype}}) {
				foreach my $base (sort {$snp{$genotype}{$total}{$a} cmp $snp{$genotype}{$total}{$b}} keys %{$snp{$genotype}{$total}}) {
					#print "$snp{$genotype}{$total}{$base}/$total\n";
					#my $c=$total;
					my $p=$snp{$genotype}{$total}{$base}/$total;
					print OE "$genotype\t$base\t$snp{$genotype}{$total}{$base}\t$total\t$p\n";
				}
			}
		}
		close OE;
	}
}






#errortable
if ($error_rand && $errortable) {

	open A ,"$error_rand" or die $!;
	my $head=<A>;
	chomp $head;
	my @head=split "\t",$head;
#	print @head ;
	my @group;
	while (<A>) {
		chomp;
		next if (/A/);
		push(@group,$_);
	}
#	print @group ;
	for (my $j=0;$j<=$#group ;$j++) {
		my @tmp=(split "\t",$group[$j]);
#		print @tmp ;
		my %error;
		open E ,"$errortable" or die $!;
		while (<E>) {
			chomp;
			next if (/^Genotype/i);
			my ($genotype,$base,$fre,$total)=(split "\t",$_)[0,1,2,3];
			$error{$genotype}{$total}{$base}=$fre;
		}
		close E;


		for (my $i=0;$i<=$#tmp;$i++) {
			print "$i\n";
			foreach my $total (keys %{$error{$head[$i]}}) {
				my $min=$error{$head[$i]}{$total}{A};
				
				my $max=$error{$head[$i]}{$total}{A};
				my $minbase="A";
				                   #
				my $maxbase="A";
				foreach my $base ( keys %{$error{$head[$i]}{$total}}) { #确定最小的，最大的base是什么，对应的值是什么
#					print "$min\t$error{$head[$i]}{$total}{$base}\n";
					if ($error{$head[$i]}{$total}{$base}<$min) {
						$min=$error{$head[$i]}{$total}{$base};
						$minbase=$base;
					}
					if ($error{$head[$i]}{$total}{$base}>$max) {
						$max=$error{$head[$i]}{$total}{$base};
						$maxbase=$base;
					}
				}
				my $next_min=$max;  #
				my $next_minbase=$maxbase; 
				foreach my $base ( keys %{$error{$head[$i]}{$total}}) { #确定第二小的base是什么，对应的值是什么
	#				print "base\t$base\tgenotype\t$head[$i]\n";
	#				print "$min\t$error{$head[$i]}{$total}{$base}\n";
					if ($error{$head[$i]}{$total}{$base}>$min  &&  $error{$head[$i]}{$total}{$base}<=$next_min) {
						$next_min=$error{$head[$i]}{$total}{$base};
						$next_minbase=$base;
			#			print "next_minbase\t$next_minbase\tgenotype\t$head[$i]\tP\t$next_min\n";
					}
				}
	#			print "minbase\t$minbase\n";
				if ($tmp[$i]==1){
					my $bb=int($min+$add_error*$min);
					my $bb2=int($next_min+$add_error*$next_min);
					my $cc=int($max-$add_error*$min-$add_error*$next_min);
					$error{$head[$i]}{$total}{$minbase}=$bb;
					$error{$head[$i]}{$total}{$next_minbase}=$bb2;
					$error{$head[$i]}{$total}{$maxbase}=$cc;
				}else{
					my $bb=int($min-$add_error*$min);
					my $bb2=int($next_min-$add_error*$next_min);
					my $cc=int($max+$add_error*$min+$add_error*$next_min);
					$error{$head[$i]}{$total}{$minbase}=$bb;
					$error{$head[$i]}{$total}{$next_minbase}=$bb2;
					$error{$head[$i]}{$total}{$maxbase}=$cc;
				}

			}
		}
		#print Dumper %error; 
		open OE ,">$out/$j.errortable" or die $!;
		foreach my $genotype (sort {$a cmp $b} keys %error) {
			foreach my $total (keys %{$error{$genotype}}) {
				foreach my $base (sort {$error{$genotype}{$total}{$a} cmp $error{$genotype}{$total}{$b}} keys %{$error{$genotype}{$total}}) {
					#print "$error{$genotype}{$total}{$base}/$total\n";
					#my $c=$total;
					my $p=$error{$genotype}{$total}{$base}/$total;
					print OE "$genotype\t$base\t$error{$genotype}{$total}{$base}\t$total\t$p\n";
				}
			}
		}
		close OE;
	}
}

my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
my $end=strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "End time id $end.\nThis compute totally consumed $duration_time s\.\n";




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

sub usage {
    die(
        qq!
Usage:   cd /56T/Evan/PMCaller/WES/P#180/DNA.C10T/PMCaller/add_error_0.4 &&  perl $0 -error_rand error_rand.txt -errortable L11_errortable   -snp_rand snp_rand.txt -snptable L11.snptable.list -o test2 
Function: Template for 引入误差，将snptable和errortable 按照特定方式（-snp_rand and -error_rand）引入误差，
Command:      
    'verbose' => \$verbose,
    'o=s' => \$out, 
    'add_error_rate=s' => \$add_error,  [0.4]
	'snp_rand=s' => \$snp_rand,         n行4列，必须有表头  
	'error_rand=s' => \$error_rand,     n行8列，必须有表头 
    'snptable=s' => \$snptable,         列表格式的snptable，必须有前四列
    'errortable=s' => \$errortable,     列表格式的errortable，必须有前四列
	'matrix' => \$matrix,
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2018/6/29
Notes:    
\n!
    )
}
