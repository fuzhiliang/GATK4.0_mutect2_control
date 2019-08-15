#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
my $verbose ="v1.0";

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my ($tumor,$normal,$outdir,$hg38,$onlyprintcmd,$gatktype);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'tumor=s' => \$tumor,
    'normal=s' => \$normal,
    'outdir=s' => \$outdir, 
	'hg38' => \$hg38,
	'gatktype=i' => \$gatktype,
    'onlyprintcmd' => \$onlyprintcmd,
#    'libs=i' => \@libs, ## 'libs=i@' => \$libs,
#    'define=s' => \%defines, ## 'define=s%' => \$defines,
) or die $!;
#&usage unless ( exists $tumor && exists $outdir );
unless(defined $outdir && defined $tumor ){&usage();exit 0;}
my $SOFT="/home/fuzl/soft";
my $fastqc="$SOFT/FastQC/fastqc";

my ($genome,$INDEX,$dbSNP,$phasel_1KG,$Mills_and_1KG,$phasel_snp_1KG);
if (defined $hg38) { #dbSNP 库还有问题
	my $hg38_root="$SOFT/GATK/resources/bundle/hg38/";
	$genome="$hg38_root/Homo_sapiens_assembly38.fasta";
	$INDEX="$hg38_root/bwa_index/gatk_hg38";
	$dbSNP="$hg38_root/dbsnp_146.hg38.vcf.gz";
	$phasel_1KG="$hg38_root/1000G_phase1.snps.high_confidence.hg38.vcf";
	$Mills_and_1KG="$hg38_root/Mills_and_1000G_gold_standard.indels.hg38.vcf";
}else {
	my $hg19_root="$SOFT/GATK/resources/bundle/hg19";
	$genome="$hg19_root/ucsc.hg19.fasta";
	$INDEX="$hg19_root/bwa_index/gatk_hg19";
	$dbSNP="$hg19_root/dbsnp_138.hg19.vcf";	
	$phasel_snp_1KG="$hg19_root/1000G_phase1.snps.high_confidence.hg19.sites.vcf";
	$phasel_1KG="$hg19_root/1000G_phase1.indels.hg19.sites.vcf";
	$Mills_and_1KG="$hg19_root/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf";
}

my $picard="$SOFT/picard.jar";
my $bin_trim_galore ="$SOFT/TrimGalore-0.5.0/trim_galore";
#my $samtools="/data2/wangb/samtools";
my $samtools="$SOFT/samtools-1.9/samtools";
#my $GATK="$SOFT/gatk-4.0.8.1/gatk";
my $GATK="$SOFT/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar";
#my $bwa="/Share/home/tiangeng/.Work1/PERL/WGS/bwa-0.7.12/bwa";
my $bwa="$SOFT/bwa-0.7.17/bwa";

$outdir||='.';
&MKDIR("$outdir");
$outdir=abs_path("$outdir");
my $tumor_fq1= "${tumor}_1.fq.gz";
my $tumor_fq2= "${tumor}_2.fq.gz";
$tumor_fq1=abs_path($tumor_fq1);
$tumor_fq2=abs_path($tumor_fq2);

my $normal_fq1= "${normal}_1.fq.gz";
my $normal_fq2= "${normal}_2.fq.gz";
$normal_fq1=abs_path($normal_fq1);
$normal_fq2=abs_path($normal_fq2);


my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is: \ntumor:$tumor\nnormal:$normal\nOutput directory is $outdir\n";
print "Database file is $genome\n";

###############################################################################
my $cmd="";

my $genomename=basename($genome);
$genome=abs_path($genome);
my $genome_m=$genome;
$genome_m=~s/\.fa(sta)?$//;#绝对路径
my $genomename_m=$genomename;
$genomename_m=~s/\.fa(sta)?$//;#基本名

unless (-f "$genome.fai" && -f "$genome_m.dict" && -f "$INDEX.bwt")  {
#	system "mkdir $outdir/genome/ && ln -s $genome $outdir/genome/";
	&MKDIR("$outdir/genome/");
	`rm $outdir/genome/* && ln -s  $genome $outdir/genome/`;
	$genome="$outdir/genome/$genomename ";
	$cmd .="cd $outdir/genome/ \n";
	$cmd .="$samtools faidx $genome  \n" ;
	$cmd .="java -Xmx20G -jar $picard CreateSequenceDictionary R=$genome O=$genomename_m.dict \n" ;
	$cmd .="bwa index -a bwtsw -p $genomename_m $genome \n";
	$INDEX ="$outdir/genome/$genomename_m";
	&runcmd("Building database",$cmd);
}

#############################################
#QC
$cmd="";

my $tumor_name=basename($tumor_fq1);
$tumor_name=~s/_1\.f(ast)?q(\.gz)?$//;
my $normal_name=basename($normal_fq1);
$normal_name=~s/_1\.f(ast)?q(\.gz)?$//;

&MKDIR("$outdir/FASTQC");
$cmd .="$fastqc $tumor_fq1 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${tumor_name}_1_fastqc.zip \n";
$cmd .="$fastqc $tumor_fq2 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${tumor_name}_2_fastqc.zip \n ";

$cmd .="$fastqc $normal_fq1 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${normal_name}_1_fastqc.zip \n ";
$cmd .="$fastqc $normal_fq2 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${normal_name}_2_fastqc.zip \n ";

#$cmd .="mkdir $outdir/trim/ \n" unless (-d "$outdir/trim")
#$cmd .="$bin_trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 --paired -o $outdir/trim  $tumor $fq2 \n";

&MKDIR("$outdir/trim/");
$cmd .="java -jar $SOFT/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 12 -phred33 -trimlog $outdir/trim/logfile ";
$cmd .="$tumor_fq1 $tumor_fq2 $outdir/trim/${tumor_name}_1_paired.fq.gz  $outdir/trim/${tumor_name}_1_unpaired.fq.gz  $outdir/trim/${tumor_name}_2_paired.fq.gz $outdir/trim/${tumor_name}_2_unpaired.fq.gz ";
$cmd .=" HEADCROP:2 ";
#$cmd .=" CROP:147 "; # reads length 151 bp，切除3‘端2bp
$cmd .="ILLUMINACLIP:$SOFT/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \n";

$cmd .="java -jar $SOFT/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 12 -phred33 -trimlog $outdir/trim/logfile ";
$cmd .="$normal_fq1 $normal_fq2 $outdir/trim/${normal_name}_1_paired.fq.gz  $outdir/trim/${normal_name}_1_unpaired.fq.gz  $outdir/trim/${normal_name}_2_paired.fq.gz $outdir/trim/${normal_name}_2_unpaired.fq.gz ";
$cmd .=" HEADCROP:2 ";
#$cmd .=" CROP:147 "; # reads length 151 bp，切除3‘端2bp
$cmd .="ILLUMINACLIP:$SOFT/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \n";


$tumor_fq1="$outdir/trim/${tumor_name}_1_paired.fq.gz";
$tumor_fq2="$outdir/trim/${tumor_name}_2_paired.fq.gz";
$tumor_name=basename($tumor_fq1);
$tumor_name=~s/_1_paired.f(ast)?q(\.gz)?$//;


$normal_fq1="$outdir/trim/${normal_name}_1_paired.fq.gz";
$normal_fq2="$outdir/trim/${normal_name}_2_paired.fq.gz";
$normal_name=basename($normal_fq1);
$normal_name=~s/_1_paired.f(ast)?q(\.gz)?$//;

&MKDIR("$outdir/FASTQC_trim/");
#$cmd .="mkdir $outdir/FASTQC_trim \n" unless (-d "$outdir/FASTQC_trim");
$cmd .="$fastqc $tumor_fq1 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${tumor_name}_1_paired_fastqc.zip \n";
$cmd .="$fastqc $tumor_fq2 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${tumor_name}_2_paired_fastqc.zip \n ";

$cmd .="$fastqc $normal_fq1 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${normal_name}_1_paired_fastqc.zip \n";
$cmd .="$fastqc $normal_fq2 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${normal_name}_2_paired_fastqc.zip \n ";

&qsub("QC",$cmd,"10G");

my $vf||="40G";
#############################################################
#bwa
my $sample=$tumor_name;
&MKDIR("$outdir/bwa");

$cmd="";
#$cmd .="mkdir $outdir/bwa \n" unless (-d "$outdir/bwa");
$cmd .="$bwa mem -M -t 6 -R '\@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:Illumina' $INDEX  $tumor_fq1 $tumor_fq2 >$outdir/bwa/$sample.sam  \n";
&qsub("BWA",$cmd,$vf);

################################
&MKDIR("$outdir/GATK_pretreatment/");
$cmd="";
#java  -Xmx20G -jar $SOFT/picard.jar AddOrReplaceReadGroups I=/data2/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread_1.fq.sam  O=/data2/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread.sort.bam  SO=coordinate  RGLB="pe"  RGPU="HiSeq-2000" RGPL=illumina RGSM=lib0705-7-2_S19_100kread
$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $picard AddOrReplaceReadGroups I=$outdir/bwa/$sample.sam  O=$outdir/GATK_pretreatment/$sample.sort.bam  SO=coordinate  RGLB=$sample  RGPU=$sample RGPL=illumina RGSM=$sample  1>$outdir/GATK_pretreatment/log.sort 2>&1 \n";
$cmd .="$samtools flagstat $outdir/GATK_pretreatment/$sample.sort.bam > $outdir/GATK_pretreatment/${sample}.alignment.flagstat & \n ";
$cmd .="$samtools stats  $outdir/GATK_pretreatment/$sample.sort.bam > $outdir/GATK_pretreatment/${sample}.alignment.stat  & \n";
$cmd .="$samtools index $outdir/GATK_pretreatment/$sample.sort.bam \n";

#################################
$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $picard   MarkDuplicates  I=$outdir/GATK_pretreatment/$sample.sort.bam  O=$outdir/GATK_pretreatment/${sample}_marked.bam M=$outdir/GATK_pretreatment/${sample}.metrics 1>$outdir/GATK_pretreatment/log.mark  2>&1 \n";
my $bam="$outdir/GATK_pretreatment/${sample}_marked.bam";
$cmd .="$samtools index  $bam \n";

#################################
$cmd .="java  -Xmx$vf -Djava.io.tmpdir=./ -jar $GATK -T RealignerTargetCreator -R  $genome -I $bam -o $outdir/GATK_pretreatment/${sample}.intervals ";
$cmd .=" -known  $phasel_1KG " if ($phasel_1KG);
$cmd .=" -known  $Mills_and_1KG " if ($Mills_and_1KG);
$cmd .=" 1>$outdir/GATK_pretreatment/log.fix 2>&1 \n";

$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $GATK  -T IndelRealigner -R $genome -I $bam -o $outdir/GATK_pretreatment/${sample}.Realgn.bam  -targetIntervals $outdir/GATK_pretreatment/${sample}.intervals ";
$cmd .=" 1>$outdir/GATK_pretreatment/log.Realgn 2>&1 \n";
$bam="$outdir/GATK_pretreatment/${sample}.Realgn.bam";
#$cmd .="$samtools index  $bam \n";

##################################recal
my $recal=1;
if ($recal) {
	$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $GATK -T  BaseRecalibrator -R $genome  -I $bam -o $outdir/GATK_pretreatment/${sample}_recal_data.table  ";
	$cmd .="--knownSites   $dbSNP " ;
	$cmd .="--knownSites   $phasel_1KG " if ($phasel_1KG);
	$cmd .="--knownSites   $Mills_and_1KG " ;
	$cmd .=" 1>$outdir/GATK_pretreatment/log.recal 2>&1 \n";

	$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $GATK -T  PrintReads  -R $genome  -I $bam "; # 
	$cmd .="-BQSR $outdir/GATK_pretreatment/${sample}_recal_data.table "; 
	$cmd .="-o  $outdir/GATK_pretreatment/${sample}_recal.bam 1>$outdir/GATK_pretreatment/log.recal 2>&1 \n";

	$bam="$outdir/GATK_pretreatment/${sample}_recal.bam";
	$cmd .="$samtools index  $bam \n";
}
#################################
&qsub("GATK pretreatment",$cmd,$vf);


###########################


$sample=$normal_name;

&MKDIR("$outdir/bwa");

$cmd="";
#$cmd .="mkdir $outdir/bwa \n" unless (-d "$outdir/bwa");
$cmd .="$bwa mem -M -t 6 -R '\@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:Illumina' $INDEX  $normal_fq1 $normal_fq2 >$outdir/bwa/$sample.sam  \n";
&qsub("BWA normal",$cmd,$vf);

################################
&MKDIR("$outdir/GATK_pretreatment/");
$cmd="";
#java  -Xmx20G -jar $SOFT/picard.jar AddOrReplaceReadGroups I=/data2/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread_1.fq.sam  O=/data2/fuzl/project/508/demo/analysis/bwa/lib0705-7-2_S19_100kread.sort.bam  SO=coordinate  RGLB="pe"  RGPU="HiSeq-2000" RGPL=illumina RGSM=lib0705-7-2_S19_100kread
$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $picard AddOrReplaceReadGroups I=$outdir/bwa/$sample.sam  O=$outdir/GATK_pretreatment/$sample.sort.bam  SO=coordinate  RGLB=$sample  RGPU=$sample RGPL=illumina RGSM=$sample  1>$outdir/GATK_pretreatment/lognormal.sort 2>&1 \n";
$cmd .="$samtools flagstat $outdir/GATK_pretreatment/$sample.sort.bam > $outdir/GATK_pretreatment/${sample}.alignment.flagstat & \n ";
$cmd .="$samtools stats  $outdir/GATK_pretreatment/$sample.sort.bam > $outdir/GATK_pretreatment/${sample}.alignment.stat  & \n";
$cmd .="$samtools index $outdir/GATK_pretreatment/$sample.sort.bam \n";

#################################
$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $picard   MarkDuplicates  I=$outdir/GATK_pretreatment/$sample.sort.bam  O=$outdir/GATK_pretreatment/${sample}_marked.bam M=$outdir/GATK_pretreatment/${sample}.metrics 1>$outdir/GATK_pretreatment/log.mark  2>&1 \n";
my $normal_bam="$outdir/GATK_pretreatment/${sample}_marked.bam";
$cmd .="$samtools index  $normal_bam \n";

#################################
$cmd .="java  -Xmx$vf -Djava.io.tmpdir=./ -jar $GATK -T RealignerTargetCreator -R  $genome -I $bam -o $outdir/GATK_pretreatment/${sample}.intervals ";
$cmd .=" -known  $phasel_1KG " if ($phasel_1KG);
$cmd .=" -known  $Mills_and_1KG " if ($Mills_and_1KG);
$cmd .=" 1>$outdir/GATK_pretreatment/lognormal.fix 2>&1 \n";

$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $GATK  -T IndelRealigner -R $genome -I $bam -o $outdir/GATK_pretreatment/${sample}.Realgn.bam  -targetIntervals $outdir/GATK_pretreatment/${sample}.intervals ";
$cmd .=" 1>$outdir/GATK_pretreatment/lognormal.Realgn 2>&1 \n";
$normal_bam="$outdir/GATK_pretreatment/${sample}.Realgn.bam";
$cmd .="$samtools index  $normal_bam \n";

##################################recal
#my $recal=1;
if ($recal) {
	$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $GATK -T  BaseRecalibrator -R $genome  -I $bam -o $outdir/GATK_pretreatment/${sample}_recal_data.table  ";
	$cmd .="--knownSites   $dbSNP " ;
	$cmd .="--knownSites   $phasel_1KG " if ($phasel_1KG);
	$cmd .="--knownSites   $Mills_and_1KG " ;
	$cmd .=" 1>$outdir/GATK_pretreatment/lognormal.recal 2>&1 \n";

	$cmd .="java -Xmx$vf -Djava.io.tmpdir=./ -jar $GATK -T  PrintReads  -R $genome  -I $bam "; # 
	$cmd .="-BQSR $outdir/GATK_pretreatment/${sample}_recal_data.table "; 
	$cmd .="-o  $outdir/GATK_pretreatment/${sample}_recal.bam 1>$outdir/GATK_pretreatment/lognormal.recal 2>&1 \n";

	$normal_bam="$outdir/GATK_pretreatment/${sample}_recal.bam";
	$cmd .="$samtools index  $normal_bam \n";
}
#################################
&qsub("GATK pretreatment normal ",$cmd,$vf);





&MKDIR("$outdir/vcf/");
$gatktype||=2;
$sample="T_${tumor_name}_N_${normal_name}";
my $vcf ;
$cmd="";
if($gatktype == "1"){
	$cmd .="java -Xmx$vf -Djava.io.tmpdir=./  -jar $GATK    -T   HaplotypeCaller  -R $genome -I $bam ";
	#$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -o  $outdir/vcf/${sample}_gatk_HC_raw.vcf 1>$outdir/vcf/log.HC  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_gatk_HC_raw.vcf ";
}elsif($gatktype == "2"){
	$cmd .="java -Xmx$vf -Djava.io.tmpdir=./  -jar $GATK  -T  MuTect2  -R $genome -I:tumor $bam ";
	$cmd .=" -I:normal $normal_bam ";
	$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -o $outdir/vcf/${sample}_Mutect2_raw.vcf 1>$outdir/vcf/log.SM  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_Mutect2_raw.vcf ";
}else{
	$cmd .="java -Xmx$vf -Djava.io.tmpdir=./  -jar $GATK  -T UnifiedGenotyper  -R $genome -I $bam ";
	#$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -o $outdir/vcf/${sample}_gatk_UG.raw.vcf 1>$outdir/vcf/log.UG  2>&1 \n";
	$vcf="$outdir/vcf/${sample}_gatk_UG.raw.vcf";
}
################################
&qsub("Genotype calling",$cmd,$vf);

$cmd ="";
$cmd .="java -jar $GATK -T  SelectVariants -R $genome -V $vcf  -selectType SNP  -o $outdir/vcf/${sample}_gatk_raw.snp.vcf  \n";
$cmd .="java -jar $GATK -T  VariantFiltration -R $genome -V $outdir/vcf/${sample}_gatk_raw.snp.vcf   -filter  \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" -filterName \"my_snp_filter\"  ";
$cmd .=" -o $outdir/vcf/${sample}.filter.snp.vcf 1>$outdir/vcf/log.snp.filter  2>&1 \n";
$cmd .="grep -w -v my_snp_filter $outdir/vcf/${sample}.filter.snp.vcf > $outdir/vcf/${sample}.filter.PASS.snp.vcf \n";
$cmd .="java -jar $GATK -T  SelectVariants -R $genome -V $vcf  -selectType INDEL  -o $outdir/vcf/${sample}_gatk_raw.InDel.vcf \n";
$cmd .="java -jar $GATK -T  VariantFiltration -R $genome  -V $outdir/vcf/${sample}_gatk_raw.InDel.vcf    -filter \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"  -filterName \"my_indel_filter\"  ";
$cmd .="-o $outdir/vcf/${sample}.filter.InDel.vcf 1>$outdir/vcf/log.InDel.filter  2>&1 \n";
$cmd .="grep -w -v my_indel_filter $outdir/vcf/${sample}.filter.InDel.vcf > $outdir/vcf/${sample}.filter.PASS.InDel.vcf \n";

#grep -v "^#"  lib-FZ18-04229F_S8.filter.PASS.InDel.vcf |cut -f 1-5,10-|sed 's/:/\t/g' |sed 's/,/\t/'|awk '{print $8/($7+$8)"\t"$0}' > lib-FZ18-04229F_S8.filter.PASS.InDel.vcf_freq
#awk '{if ($7+$8 ==0){}else{print $8/($7+$8)"\t"$0}}

#$cmd .= "grep -v \"^#\"  $outdir/vcf/${sample}.filter.PASS.snp.vcf |cut -f 1-5,10-|sed \'s\/:\/\\t\/g\' |sed \'s\/,\/\\t/\'|awk \'\{print \$8\/\(\$7+\$8\)\"\\t\"\$0\}\' \> $outdir/vcf/${sample}.filter.PASS.snp.vcf_freq \n";
$cmd .= "grep -v \"^#\"  $outdir/vcf/${sample}.filter.PASS.snp.vcf |cut -f 1-5,10-|sed \'s\/:\/\\t\/g\' |sed \'s\/,\/\\t/\'|awk '\{if \(\$7\+\$8==0\)\{\}else\{print \$8\/\(\$7\+\$8\)\"\\t\"\$0\}\}\' \> $outdir/vcf/${sample}.filter.PASS.snp.vcf_freq \n";
$cmd .= "grep -v \"^#\"  $outdir/vcf/${sample}.filter.PASS.InDel.vcf |cut -f 1-5,10-|sed \'s\/:\/\\t\/g\' |sed \'s\/,\/\\t/\'|awk '\{if \(\$7\+\$8==0\)\{\}else\{print \$8\/\(\$7\+\$8\)\"\\t\"\$0\}\}\' \> $outdir/vcf/${sample}.filter.PASS.InDel.vcf_freq \n";
&qsub("Filter",$cmd,"10G");

#&MKDIR("$outdir/ANNOVAR/");
$cmd ="";
$cmd .= "perl  $SOFT/annovar/convert2annovar.pl -format  vcf4 $outdir/vcf/${sample}.filter.PASS.snp.vcf >$outdir/vcf/${sample}.filter.PASS.snp.vcf.avinput \n ";
$cmd .= "perl $SOFT/annovar/table_annovar.pl $outdir/vcf/${sample}.filter.PASS.snp.vcf.avinput   $SOFT/annovar/humandb/ -buildver hg19 -out  $outdir/vcf/${sample}_filter_PASS_snp_gatk -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . \n";

$cmd .= "perl  $SOFT/annovar/convert2annovar.pl -format  vcf4  $outdir/vcf/${sample}.filter.PASS.InDel.vcf  >  $outdir/vcf/${sample}.filter.PASS.InDel.vcf.avinput \n ";
$cmd .= "perl $SOFT/annovar/table_annovar.pl  $outdir/vcf/${sample}.filter.PASS.InDel.vcf.avinput  $SOFT/annovar/humandb/ -buildver hg19 -out  $outdir/vcf/${sample}_filter_PASS_InDel_gatk  -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . \n";

$cmd .="awk \'NR==FNR\{a[\$2\"\\t\"\$3]=\$1\"\\t\"\$8\"\\t\"\$9\}NR>FNR\{if \(FNR==1\)\{print \"Freq\\tRefDP\\tAltDP\\t\"\$0\}else\{if(\$1\"\\t\"\$2 in a ){print a[\$1\"\\t\"\$2]\"\\t\"\$0}}}' " ;
$cmd .=" $outdir/vcf/${sample}.filter.PASS.snp.vcf_freq  $outdir/vcf/${sample}_filter_PASS_snp_gatk.hg19_multianno.txt > $outdir/vcf/${sample}_filter_PASS_snp_gatk.hg19_multianno.txt.freq \n";

#$cmd .="awk \'NR==FNR{if\(length\(\$6\)\<length\(\$5\)\){a[\$2\"\\t\"\$3+1]=\$1\"\\t\"\$8\"\\t\"\$9}else{a[\$2\"\\t\"\$3]=\$1\"\\t\"\$8\"\\t\"\$9}}NR>FNR\{if \(FNR==1\)\{print \"Freq\\tRefDP\\tAltDP\\t\"\$0\}else\{if(\$1\"\\t\"\$2 in a ){print a[\$1\"\\t\"\$2]\"\\t\"\$0}}}' " ;
#$cmd .=" $outdir/vcf/${sample}.filter.PASS.InDel.vcf_freq   $outdir/vcf/${sample}_filter_PASS_InDel_gatk.hg19_multianno.txt > $outdir/vcf/${sample}_filter_PASS_InDel_gatk.hg19_multianno.txt.freq \n";
$cmd .="awk \'NR==FNR{if\(length\(\$6\)\<length\(\$5\)\){a[\$2\"\\t\"\$3+1]=\$1\"\\t\"\$8\"\\t\"\$9}else{a[\$2\"\\t\"\$3]=\$1\"\\t\"\$8\"\\t\"\$9}}NR>FNR\{if \(FNR==1\)\{print \"Freq\\tRefDP\\tAltDP\\t\"\$0\}else\{if(\$1\"\\t\"\$2 in a ){print a[\$1\"\\t\"\$2]\"\\t\"\$0}}}' " ;
$cmd .=" $outdir/vcf/${sample}.filter.PASS.InDel.vcf_freq   $outdir/vcf/${sample}_filter_PASS_InDel_gatk.hg19_multianno.txt > $outdir/vcf/${sample}_filter_PASS_InDel_gatk.hg19_multianno.txt.freq \n";

&qsub("ANNOVAR",$cmd,"10G");

###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

###############################################################################
sub usage {
    die(
        qq!
Usage:
	eg:
nohup perl $0 -tumor /home/fuzl/data/Gastric_demo/libccc1087163_S5_L007 -normal /home/fuzl/data/Gastric_demo/libcap1087168_S17_L008 -outdir /home/fuzl/data/Gastric_demo/analysis &
Function: Template for Perl FASTQC BWA GATK mutect2 pipeline .
Command:	-tumor str	tumor  sample name   [*_1.fq, *_2.fq]
			-normal str	normal sample name  
			-outdir	outdir
#			-hg38	refrence genome version. defalt [hg19]
			-gatktype gatk type, 0=UnifiedGenotyper, 1=HaplotypeCaller, 2=Mutect2  [2]
			-onlyprintcmd
Author:   Zhiliang Fu fuzl\@geneis.cn, QQ:594380908
Version:  v1.0
Update:   2018/8/8
Notes:    192.168.10.2 
\n!
    )
}

sub qsub {
	my $name=shift @_;
	my $cmd=shift @_;
	my $vf=shift @_;
	&MKDIR("$outdir/shell/");
	my $n=$name;
	$n=~s/\s+/_/g ;
	open CMD ,">$outdir/shell/$n.sh" or die $!;
	print CMD "$cmd";
	close CMD;
	if (-e "/dev/sg17"){ 
	# `qsub -cwd -l vf=4G -q all.q@sge_c17   $cmd `;
		system "cd $outdir/shell/ && perl /data2/fuzl/script/qsub.pl -l vf=$vf -q all.q\@sge_c17  -N  $n  -s 60  -b 100  $outdir/shell/$n.sh " unless (defined $onlyprintcmd);
	}else{
		system "sh $outdir/shell/$n.sh " unless (defined $onlyprintcmd) ;
	}
}


sub runcmd { # 
	my $name=shift @_;
	my $cmd=shift @_;
	&MKDIR("$outdir/shell/");
	print "Start $name analysis ... \n";
	my $n=$name;
	$n=~s/\s+/_/g ;
	open CMD ,">$outdir/shell/$n.sh" or die $!;
	print CMD "$cmd";
	close CMD;
	system "sh $outdir/shell/$n.sh ";
	print "End $name analysis !\n";
}

sub MKDIR{
	my $dir=shift @_;
	system "mkdir  -p $dir " unless (-d "$dir");
}

#perl //data2/fuzl/project/508/script/Fastqc_bwa_gatk3.8_508_qsub_v2.pl -tumor lib0705-7-2_S19_100kread_1.fq -fq2 lib0705-7-2_S19_100kread_2.fq -outdir /data2/fuzl/project/508/demo/analysis2