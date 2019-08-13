#!/uer/bin/perl -w 

my $errortable=$ARGV[0];
my $type=$ARGV[1];

my $value=$ARGV[2];
my $o=$ARGV[3];
my $line=$ARGV[4];
die "perl $0 /8T1/Evan/PMCaller/WGS_filter/SRR6431668/SRR6431668.filter.sort.errortable.list_L9 changsingeerrortable/type.txt 0.9 changsingeerrortable/L9_0.9.error 8\n" if (@ARGV!=5);
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
open I , "$errortable" or die $!;
open O , ">$o" or die $!;
while (<I>) {
	chomp;
	next if (/^E/i);
	my @tmp=split "\t" ,$_;
	$hash{$tmp[0]}{$tmp[1]}{$tmp[3]}=$tmp[2];
}

foreach my $gt (keys %hash) {
	foreach my $bb (keys %{$hash{$gt}}) {
		foreach my $total (keys %{$hash{$gt}{$bb}} ) {
		
			foreach my $type (@type) {

				my ($base,$gtture,$gterror)=(split "\t",$type)[0,1,2];
				my @b=split "",$gtture;
				my $alt;
				foreach my $b (@b) {
					$alt=$b if ($b ne $base);
					#$ref=$b if ($b!=$base);
				}
				if (($gt eq $gterror) && ($bb eq $alt)){
					#print "$gt\t$alt\n";
					my $cc=$hash{$gt}{$bb}{$total};
					my $dd=int($cc*$value);
					my $ee=$cc-$dd;
					$hash{$gt}{$bb}{$total}=$dd;
					$hash{$gt}{$base}{$total}+=$ee;
				}
			}

		}
	}
	
}

foreach my $gt (sort {$a cmp $b} keys %hash) {
	foreach my $bb ( keys %{$hash{$gt}}) {
		foreach my $total (keys %{$hash{$gt}{$bb}} ) {
			my $p=$hash{$gt}{$bb}{$total}/$total ;
			print O "$gt\t$bb\t$hash{$gt}{$bb}{$total}\t$total\t$p\n";
		}
	}
}
close I;
close O;

=cut
CC      T       12944   63813   0.202842680958425



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
