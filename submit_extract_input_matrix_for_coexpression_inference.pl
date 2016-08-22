##############################################################################################
#	File Name: submit_extract_input_matrix_for_coexpression_inference
#	> Author: Hua Yu
#	> Mail: huayu@genetics.ac.cn 
#	Created Time: 2014-08-20
##############################################################################################

#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($inputpath,$netfilepath,$runinfofile,$htseqpath,$cufflinkspath,$RAW,$FPKM,$suffix_htseq,$suffix_cufflinks,$queue,$help);
my $usage = <<USAGE;
Usage:
perl $0 --inputpath \$inputpath --netfilepath \$netfilepath --runinfofile \$runinfofile --htseqpath \$htseqpath --cufflinkspath \$cufflinkspath --RAW --FPKM --suffix_htseq \$suffix_htseq --suffix_cufflinks \$suffix_cufflinks --help
Options:
--inputpath: input file path
--netfilepath: network file path
--runinfofile: run information file
--htseqpath: relative htseq_count results path 
--cufflinkspath: relative cufflinks results path 
--RAW: whether or not extract raw expression value per gene
--FPKM: whether or not extract FPKM value per gene
--suffix_htseq: htseq result file suffix
--suffix_cufflinks: cufflinks result file suffix
--queue: PBS queue
--help: displaying help information
USAGE
;

GetOptions(
	"inputpath|i=s" => \$inputpath,
	"netfilepath=s" => \$netfilepath,
	"runinfofile=s" => \$runinfofile,
	"htseqpath=s" => \$htseqpath,
	"cufflinkspath=s" => \$cufflinkspath,
	"RAW!" => \$RAW,
	"FPKM!" => \$FPKM,
	"suffix_htseq=s" => \$suffix_htseq,
	"suffix_cufflinks=s" => \$suffix_cufflinks,
	"queue=s" => \$queue,
	"help!" => \$help,
) or die "$usage\n";

if ($help){
	print "$usage\n";
	exit(1);
}

if ($RAW){
	my $inputpath_htseq = "$inputpath/$htseqpath";
	open OUT,">extract_htseq_input_matrix_for_coexpression_inference.pbs" or die "couldn't not open the file extract_htseq_input_matrix_for_coexpression_inference.pbs:$!\n";
	print OUT <<PBS;
#PBS -N extract_htseq_input_matrix_for_coexpression_inference
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q $queue
#PBS -o extract_htseq_input_matrix_for_coexpression_inference.log

cd \$PBS_O_WORKDIR
date
Rscript --vanilla extract_input_matrix_for_coexpression_inference.R $inputpath_htseq $suffix_htseq $runinfofile
cp $inputpath_htseq/HTseqGene_Filtered_Raw_Count.txt $netfilepath
cp $inputpath_htseq/Median_htseq_DESeq_Normalized_Count.txt $netfilepath
cp $inputpath_htseq/VST_htseq_DESeq_Normalized_Count.txt $netfilepath
cp $inputpath_htseq/TMM_htseq_EdgeR_Normalized_Count.txt $netfilepath
cp $inputpath_htseq/UQ_htseq_EdgeR_Normalized_Count.txt $netfilepath
date
PBS
	close OUT;
	my $qsub = `qsub extract_htseq_input_matrix_for_coexpression_inference.pbs`;
	if($qsub ne ""){
		print "The task of extract_htseq_input_matrix_for_coexpression_inference.pbs is submitted in $qsub\n";
	}
}	

if($FPKM){
	my $inputpath_htseq = "$inputpath/$htseqpath";
	my $inputpath_cufflinks = "$inputpath/$cufflinkspath";
	open OUT,">extract_cufflinks_input_matrix_for_coexpression_inference.pbs" or die "couldn't not open the file extract_cufflinks_input_matrix_for_coexpression_inference.pbs:$!\n";
	print OUT <<PBS;
#PBS -N extract_cufflinks_input_matrix_for_coexpression_inference
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q $queue
#PBS -o extract_cufflinks_input_matrix_for_coexpression_inference.log

cd \$PBS_O_WORKDIR
date
Rscript --vanilla extract_input_matrix_for_coexpression_inference.R $inputpath_cufflinks $suffix_cufflinks $runinfofile $inputpath_htseq
cp $inputpath_cufflinks/FPKM_Cufflinks_Normalized_Count.txt $netfilepath
date
PBS
	close OUT;
	my $qsub = `qsub extract_cufflinks_input_matrix_for_coexpression_inference.pbs`;
	if($qsub ne ""){
		print "The task of extract_cufflinks_input_matrix_for_coexpression_inference.pbs is submitted in $qsub\n";
	}
}
