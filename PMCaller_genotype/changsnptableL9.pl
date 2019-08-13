#!/uer/bin/perl -w 

my $snptable=$ARGV[0];
my $type=$ARGV[1];

my $value=$ARGV[2];
my $o=$ARGV[3];
my $line=$ARGV[4];
die "perl $0 /56T/Evan/PMCaller/WGS/SRR6431668/BWA/chr1/snptable_0.9_L9/chr1.H40.snptable.list.snp.table_list  changsingsnptable/type.txt 1.5 changsingsnptable/L9_hetx1.5.snp 11\n" if (@ARGV!=5);
my %hash;
my @type;
open T ,"$type" or die $!;
while (<T>) {
	chomp;
	push (@type,$_);
	last if ($.==$line); 
	#my @b=split "",$gtture;
	#foreach my $b (@b) {
	#	if ($b eq $base);
	#	
	#}
}

close T ;
open I , "$snptable" or die $!;
open O , ">$o" or die $!;
while (<I>) {
	chomp;
	next if (/^Ref/i);
	my @tmp=split "\t" ,$_;
	$hash{$tmp[0]}{$tmp[1]}{$tmp[3]}=$tmp[2];
}

foreach my $ref (keys %hash) {
	foreach my $gt (keys %{$hash{$ref}}) {
		foreach my $total (keys %{$hash{$ref}{$gt}} ) {
		
			foreach my $type (@type) {

				my ($ref_type,$gtture,$gterror)=(split "\t",$type)[0,1,2];
				#$my @b=split "",$gtture;
				#my $alt;
#				foreach my $b (@b) {
#					$alt=$b if ($b ne $base);
#					#$ref=$b if ($b!=$base);
#				}
				if (($gt eq $gtture) && ($ref eq $ref_type)){
					#print "$gt\t$alt\n";
					my $cc=$hash{$ref}{$gt}{$total};
					my $dd=int($cc*$value);  #增大
					my $ee=$cc-$dd;
					$hash{$ref}{$gt}{$total}=$dd;
					$hash{$ref}{$gterror}{$total}+=$ee;
				}
			}

		}
	}
	
}

foreach my $ref (sort {$a cmp $b} keys %hash) {
	foreach my $gt ( keys %{$hash{$ref}}) {
		foreach my $total (keys %{$hash{$ref}{$gt}} ) {
			my $p=$hash{$ref}{$gt}{$total}/$total ;
			print O "$ref\t$gt\t$hash{$ref}{$gt}{$total}\t$total\t$p\n";
		}
	}
}
close I;
close O;

=cut
C       AA      168     9063    0.0185369083085071
C       AC      260     9063    0.0286880723822134
C       AG      0       9063    0
C       AT      0       9063    0
C       CC      6446    9063    0.711243517599029  降
C       CG      237     9063    0.0261502813637868
C       CT      1108    9063    0.122255323844202  升
C       GG      190     9063    0.0209643605870021
C       GT      0       9063    0
C       TT      654     9063    0.0721615359152599


C	CT	CC
A	AG	AA
T	CT	TT
C	AC	CC
G	AG	GG
A	AC	AA
T	GT	TT
G	GT	GG
C	CG	CC
T	CT	CC
G	CG	GG
