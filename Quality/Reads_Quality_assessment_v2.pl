#!/usr/bin/env perl
###################################
# Author: Tania
# Email: 515430717@qq.com
###################################
#!usr/bin/perl
use strict;
#use warnings;
use Getopt::Long;
use File::Basename qw(basename dirname);
sub usage{
        print <<USAGE;
usage:
                echo \" $0 -S \$sam  -o \$outdir
                -h help
                -s \$sam
                -o \$outdir
                -q \$tqual the min threshold of base quality used for classify high quality reads
USAGE
}
#=============================parameters===========================#
my ($help,$Sam,$outdir,$tqual);
GetOptions(
        "h"=>\$help,
        "s:s"=>\$Sam,
        "o:s"=>\$outdir,
        "q:n"=>\$tqual,
);
if($help){&usage();exit 0;}
unless(defined $Sam ){&usage();exit 0;}
unless(defined $outdir){$outdir="./";}
$tqual ||=15;
open SAM,$Sam or die $!;
my $basename=basename($Sam);
open OUT,">$outdir/$basename.HighQual.qual${tqual}";
open OUT2,">$outdir/$basename.LowQual.qual${tqual}";
open S1,">$outdir/$basename.HighQual.qual${tqual}.sam";
open S2,">$outdir/$basename.LowQual.qual${tqual}.sam";
open OUT3,">$outdir/$basename.Total.HL.num.qual${tqual}";
my (%num,%depth,%mismatch,$misatchstr,$MD);
while(<SAM>){
        chomp;
        if(/^@/){
                print S1 "$_\n";
                print S2 "$_\n";
        }
        my $line=$_;
        my ($chr,$start,$Mstr,$seq,$qual)=(split (/\s+/,$line))[2,3,5,9,10];
        my $tag=IsLowQualityReads($qual,$tqual);
        next if ($Mstr !~m/^([0-9]+)M$/);
        #print "\$tag:$tag\t$qual\n";
        if($tag == 0){
                print S1 "$_\n";
        }elsif($tag == 1){
                print S2 "$_\n";
        }
        $num{$tag}++;
        my @bases=split("",$seq);
        for(my $i=0;$i<@bases;$i++){
                my $Sloc=$i+1;
                $depth{$tag}{$Sloc}++;
        }
        my (@mis,%mispos,%del,%tmp);
        $line=~/(MD:Z:\S*)/;
        $MD=$1;
        #print "outter:$MD\n";
        $misatchstr=(split(":",$MD))[2];##10A5^AC6  49G40
        if($misatchstr =~ /[ATCGN]/i){
                #print "inner:$MD\n";
                #print "\$misatchstr:$misatchstr\n";
                @mis=$misatchstr=~/([0-9]+)([A-Z]|\^[A-Z]+)/ig; 
                my $pos=0;
                for(0..@mis-1){
                        #print "$_\t$mis[$_]\n";
                        if($mis[$_] =~/([0-9]+)/){
                                $pos+=$mis[$_];
                        #}elsif(length($mis[$_])>1){
                        #       $del{$pos}=$1;
                        #       $pos+=1;
                        }elsif(length($mis[$_]) == 1){
                                $pos++;
                                $mispos{$pos}=$mis[$_];
                        }
                }
                for my $loc (sort keys %mispos){
                        $mismatch{$tag}{$loc}++;
                        $tmp{$loc}++;
                }
                #for my $loc (sort keys %del){
                #       my $len=length($del{$loc});
                #       for(1..$len){
                #               my $index=$loc+$_;
                #               $mismatch{$index}++;
                #               $tmp{$index}++;
                #       }
                #}
                #for my $loc (sort keys %mismatch){
                #       print "misaa:$loc\t$mismatch{$loc}\n";
                #}
                #for my $loc (sort keys %tmp){
                #       print "tmpaa:$loc\t$tmp{$loc}\n";
                #}
        }
}
# HighQual
print OUT "ReadLoc\tMismatchNum\tDepth\tRate\n";
for my $loc (sort {$mismatch{0}{$a} <=> $mismatch{0}{$b}} keys %{$mismatch{0}}){
        my $rate=$mismatch{0}{$loc}/$depth{0}{$loc};
        print OUT "$loc\t$mismatch{0}{$loc}\t$depth{0}{$loc}\t$rate\n";
}
# LowQuals
print OUT2 "ReadLoc\tMismatchNum\tDepth\tRate\n";
for my $loc (sort{$mismatch{1}{$a} <=> $mismatch{1}{$b}} keys %{$mismatch{1}}){
        my $rate=$mismatch{1}{$loc}/$depth{1}{$loc};
        print OUT2 "$loc\t$mismatch{1}{$loc}\t$depth{1}{$loc}\t$rate\n";
}
print OUT3 "Low Quality reads num:$num{1}\tHigh Quality reads num:$num{0}\n";
close OUT;
close OUT2;
close OUT3;
open SH,">$outdir/$basename.qual${tqual}.Error.Profile.R";
print SH "filename='$outdir/$basename.HighQual.qual${tqual}';\n";
print SH "output='$outdir/$basename.HighQual.qual${tqual}.Error.Profile';\n";
print SH "data=read.table(filename,head=T)\n";
print SH "pdf(paste(output,'pdf',sep='.'))\n";
print SH "plot(data\$ReadLoc,data\$Rate*100,type='h',xlab='Read Position',ylab='Error rate (%)',main='Error Profile' )\n";
print SH "dev.off()\n";

print SH "filename='$outdir/$basename.LowQual.qual${tqual}';\n";
print SH "output='$outdir/$basename.LowQual.qual${tqual}.Error.Profile';\n";
print SH "data=read.table(filename,head=T)\n";
print SH "pdf(paste(output,'pdf',sep='.'))\n";
print SH "plot(data\$ReadLoc,data\$Rate*100,type='h',xlab='Position',ylab='Error rate (%)',main='Error Profile' )\n";
print SH "dev.off()\n";;
close SH;

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
