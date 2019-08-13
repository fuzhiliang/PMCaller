#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
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
###############################################################################
my ($sam,$bam,$L_list,$cut_off,$vcf,$out,$genome,$threads,$cut_H,$windows,$downsample,$quality,$reads_cut,$printcmd);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'sam=s' => \$sam,
	'bam=s' => \$bam,
    'vcf=s' => \$vcf,
    'o=s' => \$out, 
    'reference=s' => \$genome,
    'L=s' => \$L_list,
    'cut_off=f' => \$cut_off,
    'cut_H=f' => \$cut_H,
	'windows=i' => \$windows,
	'downsample=f' => \$downsample,
	'threads=i' => \$threads,
	'quality=i' => \$quality,
	'reads_cut=i'=> \$reads_cut,
	'printcmd'=> \$printcmd,
) or die $!;
unless(defined $sam || defined $bam ){&usage();exit 0;}

my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
my $start = strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is  $bam\nOutput file is $out\n";

my $samtools = "/usr/local/bin/samtools";
$out ||= "./";
system "mkdir  -p $out" unless (-d $out);

$genome ||= "/data/hg19/hg19_index/hg19.fa";
my $refgene ="/data/hg19_anno/Region/refgene.nooverlap";
my $genome_name=basename($genome);
unless (-f "$genome.fai") {
	system "mkdir  $out/genome/ && ln -s $genome $out/genome/";
	`$samtools faidx $out/genome/$genome_name `;
	$genome="$out/genome/$genome_name";
}
my $t=1;
$L_list ||= ("3:5:7:9:11:13:15:17:19:21");
$cut_off ||= "0.05";
$threads ||=30;
$cut_H||=40;
my $H_min||=30;
my $H_max||=600;
$windows||=0;
$reads_cut||=0;
`mkdir $out/shell` unless (-d "$out/shell");
open LOG ,">$out/shell/a.log";
print LOG "$start\n";

my $sample;
my $cmd="";
if (-f $bam){
	$sample=basename($bam);
	$sample=~s/(\.sort(ed)?)?\.bam$//;
	unless ($bam=~/sort/) {
		$cmd .= "$samtools sort  --threads $threads  $bam >  $out/$sample.sorted.bam \n";
		#system "$samtools sort  --threads $threads  $bam >  $out/$sample.sorted.bam";
		&runcmd("sort",$cmd);
		$bam="$out/$sample.sorted.bam";
	}
}elsif(-f $sam){
	$sample=basename($sam);
	$sample=~s/\.sam$//;
	$cmd .= "$samtools view  --threads $threads -Sb $sam |$samtools sort --threads $threads > $out/$sample.sorted.bam \n";
	#system "$samtools view  --threads $threads -Sb $sam |$samtools sort --threads $threads > $out/$sample.sorted.bam ";
	&runcmd("sam2bam",$cmd);
	$bam = "$out/$sample.sorted.bam ";
}else{
	die "Input file error. please input format file of bam or sam !" ;
}
$cmd ="";
unless (defined $vcf) {
	$cmd.="$samtools  mpileup -t DP,AD -uvf $genome $bam > $out/$sample.vcf \n ";
	$vcf=" $out/$sample.vcf";
	&runcmd("mileup",$cmd );
}


#samtools mpileup -t DP,AD -uvf /data/hg19/hg19_index/hg19.fa /56T/Evan/PMCaller/mydata/${sample}_1.sorted.bam > /56T/Evan/PMCaller/mydata/${sample}_1.sorted.vcf
#多个H 
$cmd="";
my $vcf_name=basename($vcf);
$vcf_name=~s/\.vcf$//;
my @L_list=(split ":",$L_list);
$cmd .= "perl $Bin/bin/vcf_expansion_H60_120.pl   -vcf_list  $vcf -o  $out/snptable -Hmin $H_min -Hmax $H_max \n ";
foreach my $L (@L_list) {
	#$cmd .= "perl $Bin/bin/ski_Expansionvcf_bam_H70_L.pl -vcf $out/snptable/$vcf_name.exp.vcf  -o $out/snptable/H${cut_H}_L$L -cut_H $cut_H  -L $L &\n";
	$cmd .= "perl $Bin/bin/ski_Expansionvcf_bam_H70_windows_x_L.pl -vcf $out/snptable/$vcf_name.exp.vcf   -o $out/snptable/H${cut_H}_L$L/ -cut_H $cut_H  -windows  $windows -L $L \n"; #  $out/snptable/H${cut_H}_L$L/$sample.snp.table_list
}
##error table 
$cmd .="perl $Bin/bin/vcf_split.pl $vcf $out/$sample.snp.gy $cut_off 20 \n" ; # unless (-s "$out/$sample.snp.gy"); #通过4个cutoff确定 的高可信度的位点
$cmd .= "awk '{if (\$5>=50){print \$0}}'  $out/$sample.snp.gy > $out/$sample.Hge50.snp.gy \n ";
#if (defined $quality) {
if (-f "$out/errortable/$sample.gt") {
	$cmd .="perl $Bin/bin/errortable_highlow.pl $out/errortable/$sample.gt $out/errortable $L_list $sample  \n";
}else{
	$cmd .="perl $Bin/bin/gy_sam_Error_Table_High_Low_v2.1.pl $out/$sample.snp.gy   $bam  $out/errortable/  $L_list $cut_H";
	if ($quality) {
		$cmd .= " $quality $reads_cut \n" ;  # 输出 $out/errortable/$sample.errortable.list
	}else{
		$cmd .= " 20 $reads_cut \n" ; 
	}
}
$cmd .="perl $Bin/bin/errortable_exon_int.pl $out/errortable/$sample.gt  $out/errortable $L_list  $sample  $refgene \n ";

foreach my $L (@L_list) {
	$cmd .="awk '{if (\$1==$cut_H && \$2==$L){print \$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$5/\$6}}' $out/errortable/$sample.errortable.list > $out/errortable/$sample.errortable.list_$L \n";
	$cmd .="awk '{if (\$1==$cut_H && \$2==$L){print \$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$5/\$6}}' $out/errortable/$sample.HighQual.errortable.list > $out/errortable/$sample.errortable.list_$L.HighQual \n";
	$cmd .="awk '{if (\$1==$cut_H && \$2==$L){print \$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$5/\$6}}' $out/errortable/$sample.LowQual.errortable.list > $out/errortable/$sample.errortable.list_$L.LowQual  \n";
	$cmd .="awk '{if (\$1==$cut_H && \$2==$L){print \$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$5/\$6}}' $out/errortable/$sample.EXON.errortable.list > $out/errortable/$sample.errortable.list_$L.EXON \n";
	$cmd .="awk '{if (\$1==$cut_H && \$2==$L){print \$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$5/\$6}}' $out/errortable/$sample.INT.errortable.list > $out/errortable/$sample.errortable.list_$L.INT  \n";
}
&runcmd("snp_and_errortable",$cmd );



my $picard;
if (-s  "/soft/picard-tools-2.4.1/picard.jar ") {
	$picard="/soft/picard-tools-2.4.1/picard.jar " ;
}else{
	$picard="/soft/picard-tools-2.3.0/picard.jar " ;
}
$downsample ||=0.1;
my $a=$downsample ;
#my $i = 5 ;
$cmd="";
$cmd .= "samtools index $bam \n";
$cmd .= "mkdir $out/down$a \n" unless (-d "$out/down$a"); 
$cmd .= "java -jar $picard  DownsampleSam I=$bam O=$out/down$a/$sample.sorted.DP$a.bam P=$a \n";
$cmd .= "samtools index $out/down$a/$sample.sorted.DP$a.bam \n";

$cmd .= "samtools  mpileup -t DP,AD -uvf $genome  $out/down$a/$sample.sorted.DP$a.bam > $out/down$a/$sample.vcf \n";
$cmd .= "grep \",<\\*>\"  $out/down$a/$sample.vcf >$out/down$a/$sample.vcf.m \n";
$cmd.= "perl  $Bin/PMCaller_genotype/PMCaller_genotype_EXON_INT_preprocess.pl $out/down$a/$sample.vcf.m  $refgene \n";   #输出 $out/down$a/$sample.vcf.m.EXONINT.modify

$cmd .= "java -Xmx25g -jar /soft/GATK/GenomeAnalysisTK.jar    -R  $genome   -T UnifiedGenotyper -I $out/down$a/$sample.sorted.DP$a.bam ";
$cmd .= " -o $out/down$a/$sample.DP$a.gatk.vcf   -stand_call_conf 0.0 -stand_emit_conf 0.0 -dcov 200  -nt 8  \n";
##java -Xmx25g -jar /soft/GATK/GenomeAnalysisTK.jar    -R  /data/hg19/hg19_index/hg19.fa   -T UnifiedGenotyper -I $root/PMCaller_vs_gatk_snpwindows/down$a/$sample.sorted.DP$a.bam -o $root/PMCaller_vs_gatk_snpwindows/down$a/$sample.sorted.DP$a.gatk.vcf   -stand_call_conf 0.0 -stand_emit_conf 0.0 -dcov 200  -nt 8 &
$cmd .= "cut -f 1-5,9,10  $out/down$a/$sample.DP$a.gatk.vcf|grep -v \"\#\"|sed 's/:/\\t/g'|cut -f 1-5,11,12| ";
$cmd .= "awk '{if(\$6==\"0/1\"){print \$0\"\\t\"\$4\$5}else if (\$6==\"1/1\"){print \$0\"\\t\"\$5\$5}else {print \$0\"\\t\"\$5}}'|";
$cmd .= "sed 's/,//g'|sed 's/CA/AC/'| sed 's/GA/AG/'| sed 's/TA/AT/'| sed 's/GC/CG/'|sed 's/TC/CT/' | sed 's/TG/GT/' > $out/down$a/$sample.DP$a.gatk.sed \n ";

&runcmd("downsample_and_gatk",$cmd);

foreach my $L (@L_list) {
	$cmd ="";
	$cmd .= "mkdir $out/PMCallerdown$L  \n" unless (-d "$out/PMCallerdown$L"); 
	$cmd .= "grep  \"DP=$L;\"  $out/down$a/$sample.vcf.m |awk '{print \$10\"\\t\"\$0}'|awk -F \"[:\\t]\" '{print \$3\"\\t\"\$0}' |cut -f 1,3-|awk -F \"[,\\t]\" '{if(\$2 != \"0\"){print \$0}}' |cut -f 2- | sed 's/:/\\t/g' | cut -f 1-5,13,14 | awk '{if(length(\$4)==1 ){print \$0}}' "; 
	$cmd .= " >$out/PMCallerdown$L/dubious_v2 \n";
	$cmd .= "cut -f 1,2,4,7  $out/PMCallerdown$L/dubious_v2 |sed 's/,/\t/g'|cut -f 1-5 | awk '{if (\$4==0 || \$5==0){}else{print \$0}}'> $out/PMCallerdown$L/dubious_v3 \n";
	$cmd .= "awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$3]=\$4}NR>FNR{if (\$1\"\\t\"\$2\"\\t\"\$3 in a){print \$0\"\\t\"a[\$1\"\\t\"\$2\"\\t\"\$3]}}' $out/$sample.Hge50.snp.gy  $out/PMCallerdown$L/dubious_v3 > $out/PMCallerdown$L/dubious_v4 \n";		
	&runcmd("dubious$L",$cmd);
}


if (defined $quality){
	$cmd="";
	$cmd .="perl $Bin/Reads_Quality_assessment_Split_High_Low_sam.pl -s $out/down$a/$sample.sorted.DP$a.bam  -q $quality  -o  $out/down$a/ \n";
	$cmd .="samtools view  $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.sam -b > $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.bam  && rm   $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.sam \n" ; 
	$cmd .="samtools  mpileup -t DP,AD -uvf $genome $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.bam  > $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.vcf \n";
	$cmd .="samtools view  $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.sam -b > $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.bam  && rm   $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.sam \n" ; 
	$cmd .="samtools  mpileup -t DP,AD -uvf $genome $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.bam  > $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.vcf \n";
	$cmd .="perl $Bin/PMCaller_genotype/PMCaller_genotype_High_Low_preprocess.pl  $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.vcf  HighQual \n";
	$cmd .="perl $Bin/PMCaller_genotype/PMCaller_genotype_High_Low_preprocess.pl  $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.vcf  LowQual \n";
	#$cmd .="awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$3]=\$0}NR>FNR{if(\$1\"\\t\"\$2\"\\t\"\$3 in a){print a[\$1\"\\t\"\$2\"\\t\"\$3]\"\\t\"\$0}}' $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.vcf.modify $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.vcf.modify|cut -f 1-7,11-|awk '{if(NR==1){printf \$0\"\\tTDepth\\tDRef_HighQual\\tDAlt_HighQual\\tDRef_LowQual\\tDAlt_LowQual\\n\"}else{printf \$0\"\\t\"\$5+\$9\"\\t\"\$6\"\\t\"\$7\"\\t\"\$10\"\\t\"\$11\"\\n\"}}'|cut -f 1-4,12- > $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.share.merge \n";
	$cmd .="perl /56T/Evan/PMCaller/script/PMCaller_genotype/merge.modify.pl   $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.vcf.modify $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.vcf.modify  $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.share.merge \n";
	$cmd .="awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$3]=\$0}NR>FNR{if(\$1\"\\t\"\$2\"\\t\"\$3 in a){}else{if (\$4 !=\"NN\"){print \$0\"\\t0\\t0\\t\"\$6\"\\t\"\$7}}}' $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.vcf $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.vcf.modify  |cut -f 1-5,8- > $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.LowQual.unique  \n ";
	$cmd .="awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$3]=\$0}NR>FNR{if(\$1\"\\t\"\$2\"\\t\"\$3 in a){}else{if (\$4 !=\"NN\"){print \$0\"\\t0\\t0\"}}}' $out/down$a/$sample.sorted.DP$a.LowQual.qual$quality.vcf.modify $out/down$a/$sample.sorted.DP$a.HighQual.qual$quality.vcf.modify > $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.HighQual.unique  \n" ;
	$cmd .="cat $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.share.merge  $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.LowQual.unique $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.HighQual.unique |awk '{if(NR==1){print \$0}else{if(\$6+\$8>=1 && \$7+\$9>=1){print \$0}}}' > $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.dubious.merge \n";
	&runcmd("quality",$cmd);
	foreach my $L (@L_list) {
		$cmd ="";
		$cmd .=" perl $Bin/PMCaller_genotype/PMCaller_genotype_High_Low.pl -vcfmerge  $out/down$a/$sample.sorted.DP$a.HighQual.LowQual.qual${quality}.dubious.merge  -o $out/PMCallerdown$L -snptable  $out/snptable/H${cut_H}_L$L/$sample.exp.snp.table_list  -errortablepre $out/errortable/$sample.errortable.list_$L  -q $quality \n";
		$cmd .=" awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$3]=\$4}NR>FNR{if (\$1\"\\t\"\$2\"\\t\"\$3 in a){print \$0\"\\t\"a[\$1\"\\t\"\$2\"\\t\"\$3]}else{print \$0\"\\tNN\"}}'  $out/PMCallerdown$L/genotype.HighQual_LowQual.qual${quality}   $out/PMCallerdown$L/dubious_v4  > $out/PMCallerdown$L/HighQual_LowQual.qual${quality}_pmcaller_vs_GT \n ";
		$cmd .=" awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$4]=\$8}NR>FNR{if (\$1\"\\t\"\$2\"\\t\"\$3 in a){print \$0\"\\t\"a[\$1\"\\t\"\$2\"\\t\"\$3]}else{print \$0\"\\t\"\$3\$3}}'  $out/down$a/$sample.DP$a.gatk.sed   $out/PMCallerdown$L/HighQual_LowQual.qual${quality}_pmcaller_vs_GT  >  $out/PMCallerdown$L/gatk_vs_HighQual_LowQual.qual${quality}_pmcaller_vs_GT  \n";
		$cmd .="mkdir $out/result \n" unless (-d "$out/result ");
		$cmd .="awk '{if (\$6==\$7){}else{print \$0}}' $out/PMCallerdown$L/gatk_vs_HighQual_LowQual.qual${quality}_pmcaller_vs_GT |wc -l |awk '{print \"PMCaller\\t\"\$0}'> $out/result/$L.HighLowErrorcount \n " ;
		$cmd .="awk '{if (\$6==\$8){}else{print \$0}}' $out/PMCallerdown$L/gatk_vs_HighQual_LowQual.qual${quality}_pmcaller_vs_GT |wc -l |awk '{print \"GATK\\t\"\$0}' >> $out/result/$L.HighLowErrorcount \n ";
		$cmd .="cat $out/PMCallerdown$L/gatk_vs_HighQual_LowQual.qual${quality}_pmcaller_vs_GT |wc -l |awk '{print \"Dubious\\t\"\$0}' >> $out/result/$L.HighLowErrorcount \n " ;
		&runcmd2("pmcaller_high_low.$L",$cmd);
	}
}

foreach my $L (@L_list) {
	$cmd="";
	$cmd .="perl $Bin/PMCaller_genotype/PMcaller_genotype.pl -vcf $out/down$a/$sample.vcf.m  -o $out/PMCallerdown$L/  -snptable  $out/snptable/H${cut_H}_L$L/$sample.exp.snp.table_list  -errortable $out/errortable/$sample.errortable.list_$L \n";
	$cmd .="awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$3]=\$4}NR>FNR{if (\$1\"\\t\"\$2\"\\t\"\$3 in a){print \$0\"\\t\"a[\$1\"\\t\"\$2\"\\t\"\$3]}else{print \$0\"\\tNN\"}}' $out/PMCallerdown$L/genotype  $out/PMCallerdown$L/dubious_v4  >  $out/PMCallerdown$L/pmcaller_vs_GT \n";
	$cmd .="awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$4]=\$8}NR>FNR{if (\$1\"\\t\"\$2\"\\t\"\$3 in a){print \$0\"\\t\"a[\$1\"\\t\"\$2\"\\t\"\$3]}else{print \$0\"\\t\"\$3\$3}}' $out/down$a/$sample.DP$a.gatk.sed   $out/PMCallerdown$L/pmcaller_vs_GT  >  $out/PMCallerdown$L/gatk_pmcaller_vs_GT \n";
	$cmd .="mkdir $out/result \n" unless (-d "$out/result ");
	$cmd .="awk '{if (\$6==\$7){}else{print \$0}}' $out/PMCallerdown$L/gatk_pmcaller_vs_GT |wc -l |awk '{print \"PMCaller\\t\"\$0}'>$out/result/$L.Errorcount \n " ;
	$cmd .="awk '{if (\$6==\$8){}else{print \$0}}' $out/PMCallerdown$L/gatk_pmcaller_vs_GT |wc -l |awk '{print \"GATK\\t\"\$0}'>> $out/result/$L.Errorcount \n " ;
	$cmd .="cat $out/PMCallerdown$L/gatk_pmcaller_vs_GT |wc -l |awk '{print \"Dubious\\t\"\$0}' >> $out/result/$L.Errorcount \n " ;
	&runcmd2("pmcaller$L",$cmd);
}

foreach my $L (@L_list) {
	$cmd ="";
	#$cmd.="perl $Bin/bin/errortable_exon_int.pl /56T/Evan/PMCaller/NA12878/chr22/errortable/chr22.filter.gt  /56T/Evan/PMCaller/NA12878/chr22/errortable $L_list  $sample  $refgene \n";
	$cmd.= "perl $Bin/PMCaller_genotype/PMCaller_genotype_Exon_Int.pl -vcfexonint $out/down$a/$sample.vcf.m.EXONINT.modify -o  $out/PMCallerdown$L/  -snptable $out/snptable/H${cut_H}_L$L/$sample.exp.snp.table_list -errortablepre $out/errortable/$sample.errortable.list_$L \n";

	$cmd .="awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$3]=\$4}NR>FNR{if (\$1\"\\t\"\$2\"\\t\"\$3 in a){print \$0\"\\t\"a[\$1\"\\t\"\$2\"\\t\"\$3]}else{print \$0\"\\tNN\"}}' $out/PMCallerdown$L/genotype.EXONINT  $out/PMCallerdown$L/dubious_v4  >  $out/PMCallerdown$L/EXONINT.pmcaller_vs_GT \n";
	$cmd .="awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$4]=\$8}NR>FNR{if (\$1\"\\t\"\$2\"\\t\"\$3 in a){print \$0\"\\t\"a[\$1\"\\t\"\$2\"\\t\"\$3]}else{print \$0\"\\t\"\$3\$3}}' $out/down$a/$sample.DP$a.gatk.sed   $out/PMCallerdown$L/EXONINT.pmcaller_vs_GT  >  $out/PMCallerdown$L/EXONINT.gatk_pmcaller_vs_GT \n";
	$cmd .="mkdir $out/result \n" unless (-d "$out/result ");
	$cmd .="awk '{if (\$6==\$7){}else{print \$0}}' $out/PMCallerdown$L/EXONINT.gatk_pmcaller_vs_GT |wc -l |awk '{print \"PMCaller\\t\"\$0}' >$out/result/$L.EXONINTErrorcount \n " ;
	$cmd .="awk '{if (\$6==\$8){}else{print \$0}}' $out/PMCallerdown$L/EXONINT.gatk_pmcaller_vs_GT |wc -l |awk '{print \"GATK\\t\"\$0}' >> $out/result/$L.EXONINTErrorcount \n " ;
	$cmd .="cat $out/PMCallerdown$L/EXONINT.gatk_pmcaller_vs_GT |wc -l |awk '{print \"Dubious\\t\"\$0}'>> $out/result/$L.EXONINTErrorcount \n " ;
	&runcmd2("pmcaller_EXONINT.$L",$cmd);
}

###############################################################################
#Record the program running time!
# 输出程序运行时间
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
my $end=strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print LOG "End time id $end.\nThis compute totally consumed $duration_time s\.\n";
close LOG;
###############################################################################

sub usage {
    die(
        qq!
Usage:    perl $0  -o   /56T/Evan/PMCaller/NA12878/chr22_exonint/  -L 3:5:9:15:19 -cut_H 40 -window 0 -bam  /56T/Evan/PMCaller/NA12878/chr22/chr22.filter.sorted.bam -downsample 0.05 -quality 20 -vcf /56T/Evan/PMCaller/NA12878/chr22/chr22.filter.vcf
Function: Template for Perl  通过bam或 bam和vcf文件call genotype，并跟GATK做比较  SNP 和Error table 都用平均深度
Command:      
	verbose	str	verbose,
	#sam		str sam,  
	-bam	str	bam,			 #必须的
	-vcf	str vcf,			 #samtools mpileup 生成的结果，有就给
	-o		str out,  
	-reference	str	genome,   # s1 服务器务必加上这个参数
	
	-L=s	int list	L_list,          #削薄到多少x，“：“号分割   ["3:5:7:9:11:13:15:17:19:21"]
    -cut_off	float	cut_off,   #error table 的cut off  [0.05]
	-cut_H	int	cut_H,       #snp table  削薄前的平均深度  [70]
	-downsample  float  downsample , #downsample的参数 [0.1]
	-threads	int	threads,
 	-quality	int quality,   #分质量值   [20]
	-reads_cut    \$reads_cut,
	-printcmd   #只输出shell，而不执行
Author:   Evan Fu, *\@163.com, QQ:42789409
Version:  v1.0
Update:   2018/8/8
Notes: 修改error table 取特定深度50 改为取接近平均深度   
\n!
    )
}

sub runcmd{
	my $name=shift @_;
	my $cmd=shift @_;
	`mkdir "$out/shell/"` unless (-d "$out/shell");
	open S ,">$out/shell/$name.sh" or die $!;
	print S "$cmd ";
	print "Start analysis of $name... \n";
	system "sh $out/shell/$name.sh "  unless ($printcmd) ;

	close S;
}

sub runcmd2{
	my $name=shift @_;
	my $cmd=shift @_;
	`mkdir "$out/shell/"` unless (-d "$out/shell");
	open S ,">$out/shell/$name.sh" or die $!;
	print S "$cmd ";
	print "Start analysis of $name... \n";
	system "nohup sh $out/shell/$name.sh &"  unless ($printcmd);

	close S;
}