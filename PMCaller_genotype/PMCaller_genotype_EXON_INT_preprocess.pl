#! /usr/bin/perl -w 
my $vcf=$ARGV[0];
my $gtf=$ARGV[1];

$gtf ||= "/data/hg19_anno/Region/refgene.nooverlap";
my %gtf;
open G ,"$gtf" or die $!;  #chr1    11873   29370
while (<G>) {
	chomp;
	my ($chr,$start,$end)=(split "\t",$_);
	$gtf{$chr}{$start}=$end;
}

open V,"$vcf" or die $!;
open O,">$vcf.EXONINT.modify" or die $!;
print O "Chr\tLoc\tRef\tAlt\tTDepth\tDRef\tDAlt\tLocation\n"; 
#print O "Chr\tLoc\tRef\tAlt\tTDepth\tDRef\tDAlt\n";
while (<V>) {
	chomp;
	next if (/^#/ || /^[[<]/ );	
	my @tmp= (split /\t/,$_);
	next if ($tmp[3]!~/[ATCG]/i || length($tmp[3]) !=1);
	my $gene= &exonorint($tmp[0],$tmp[1]);
	my $ref=uc($tmp[3]);
	my @alt= (split ",",$tmp[4])[0];
	next if(length($alt[0]) !=1) ;
#	{
#		$alt[0]="NN";
#	}
	#pop @alt;
	my ($dp,$ad) = (split ":",$tmp[9])[1,2];	
	my @ad = split "," ,$ad;
	$alt[0]=uc($alt[0]);
	print O "$tmp[0]\t$tmp[1]\t$tmp[3]\t$alt[0]\t$dp\t$ad[0]\t$ad[1]\t$gene\n";
}
close V;
close O;


sub exonorint{ # ($chr,$loc)
	my $chr=shift @_;
	my $loc=shift @_;
	my $gene="INT";
	my @tmp = split "\t", $_;
	if (exists $gtf{$chr} ) {
		foreach my $start (keys %{$gtf{$chr}}) {
			if ($loc >= $start && $loc <= $gtf{$chr}{$start}){
				$gene="EXON";		
				last;
			}
		}
	}
	return $gene;
}