#! /usr/bin/perl -w
my $vcf=$ARGV[0];
open V,"$vcf" or die $!;
open O,">$vcf.modify" or die $!;
#print O "Chr\tLoc\tRef\tAlt\tTDepth_${tag}\tDRef_${tag}\tDAlt_${tag}\n"; 
print O "Chr\tLoc\tRef\tAlt\tTDepth\tDRef\tDAlt\n";
while (<V>) {
	chomp;
	next if (/^#/ || /^[[<]/ );	
	my @tmp= (split /\t/,$_);
	next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
	my $ref=uc($tmp[3]);
	my @alt= (split ",",$tmp[4])[0];
	next if (length($alt[0]) !=1);
	#pop @alt;
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];	
	my @ad = split "," ,$ad;
	next if ($ad[0]*$ad[1]==0);
	$alt[0]=uc($alt[0]);
	print O "$tmp[0]\t$tmp[1]\t$tmp[3]\t$alt[0]\t$dp\t$ad[0]\t$ad[1]\n";
}
close V;
close O;