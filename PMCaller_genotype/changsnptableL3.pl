#!/uer/bin/perl -w 

my $errortable=$ARGV[0];
my $genotype=$ARGV[1];
my $base=$ARGV[2];
my $value=$ARGV[3];
my $o=$ARGV[4];

die "perl $0 /8T1/Evan/PMCaller/WGS_filter/SRR6431668/SRR6431668.filter.sort.errortable.list_L3 CC T 0.2 changsingeerrortable/L3CC0.2.error \n" if (@ARGV!=5);
open I , "$errortable" or die $!;
open O , ">$o" or die $!;
while (<I>) {
	chomp;
	next if (/^E/i);
	my @tmp=split "\t" ,$_;
	if ($tmp[0] eq $genotype){
		my $b=(split "",$tmp[0])[0];
		if ($tmp[1] eq $base) {
			$tmp[2] +=int($tmp[3]*$value);
		}
		if ($tmp[1] eq $b) {
			$tmp[2] -=int($tmp[3]*$value);
		}
	}
	$tmp[4] =$tmp[2]/$tmp[3];
	my $t=join "\t",@tmp;
	print O "$t\n";
}
close I;
close O;