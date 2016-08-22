###########################################################################################
# This program submit the reads fixation program
# Author: Yu Hua 
# Created: 2014-7-16
###########################################################################################

#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my ($inputpath,$cleanpath,$pbspath,$queue,$threads,$help);
my $usage = <<USAGE;
Usage:
perl $0 --cleanpath \$cleanpath --pbspath \$pbspath --queue \$queue --threads \$threads --help \$help
Options: 
--cleanpath: clean file path
--pbspath: PBS script path
--queue: PBS queue
--threads: the number of threads
--help: displaying help information
USAGE
;

GetOptions(
	"inputpath|i=s" => \$inputpath,
	"cleanpath=s" => \$cleanpath,
	"pbspath=s" => \$pbspath,
	"queue=s" => \$queue,
	"threads=i" => \$threads,
	"help!" => \$help,
) or die "$usage\n";

if ($help){
	print "$usage\n";
}

my @samplereadsfiles =`find $cleanpath -name "*clean_PE_1.fastx.fastq.gz"`;
foreach my $samplereads1file (@samplereadsfiles){
	chomp $samplereads1file;
	$samplereads1file =~ /(.*)\/(.*)\/(.*)\.clean_PE_1\.fastx\.fastq\.gz/;
	my $study =$2;
	my $runid = $3;
	my $samplereads2file=$samplereads1file;
	$samplereads2file =~ s/clean_PE_1\.fastx\.fastq\.gz/clean_PE_2\.fastx\.fastq\.gz/;
	open(PBS,">$inputpath/$pbspath/$study/fix_fastx_clean_reads_$runid.pbs") or die "couldn't open the file fix_fastx_clean_reads_$runid.pbs:$!";
	print PBS <<EOF;
#PBS -N fix_fastx_clean_reads_$runid
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q $queue
#PBS -o $pbspath/$study/fix_fastx_clean_reads_$runid.log

cd \$PBS_O_WORKDIR
perl fixFASTXCleanReads.pl --inputpath $inputpath --cleanpath $cleanpath --samplereads1file $samplereads1file --samplereads2file $samplereads2file --threads $threads  
EOF
	close PBS;

	my $qsub = `qsub $pbspath/$study/fix_fastx_clean_reads_$runid.pbs`;
	if ($qsub ne ""){
		my $lc = localtime;
		print "$lc\tThe task of fix_fastx_clean_reads_$runid is submitted in $qsub\n"
	}
}


