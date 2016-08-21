###########################################################################################
# This program submit the reads fixation program
# Author: Yu Hua 
# Created: 2014-6-26
###########################################################################################

#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my ($inputdir,$cleanpath,$pbspath,$queue,$threads,$help);
my $usage = <<USAGE;
Usage:
perl $0 --inputdir \$inputdir --cleanpath \$cleanpath --pbspath \$pbspath --queue \$queue --threads \$threads --help \$help
Options:
--inputdir: input file directory
--cleanpath: clean file directory
--pbspath: pbs scripts directory
--queue: the computing queue
--threads: the number of threads
--help: display help information
USAGE
;

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"cleanpath=s" => \$cleanpath,
	"pbspath=s" => \$pbspath,
	"queue=s" => \$queue,
	"threads=i" => \$threads,
	"help!" => \$help,
) or die "$usage\n";

if ($help){
	print "$usage\n";
}

open(OUT,">$inputdir/$cleanpath/fastxCleanInfo.csv");
print OUT "Study,Run,PairedNum,PairedBase,UnPairedNum,UnPairedBase\n";
close OUT;

my @samples =`find $inputdir/$cleanpath -name "*clean_PE_1.fastx.fastq.gz"`;
foreach my $sample1 (@samples){
	chomp $sample1;
	$sample1 =~ /(.*)\/(.*)\/(.*)\.clean_PE_1\.fastx\.fastq\.gz/;
	my $study =$2;
	my $runid = $3;
	my $sample2=$sample1;
	$sample2 =~ s/clean_PE_1\.fastx\.fastq\.gz/clean_PE_2\.fastx\.fastq\.gz/;
	print "$sample2\n";
	open(PBS,">$inputdir/$pbspath/$study/fixFastxReadsPairs_$runid.pbs") or die "Couldn't open the file fixFastxReadsPairs_$runid.pbs:$!";
	print PBS <<EOF;
#PBS -N fixFastxReadsPairs_$runid
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q $queue
#PBS -o $inputdir/$pbspath/$study/fixFastxReadsPairs_$runid.log

cd \$PBS_O_WORKDIR
perl fixFASTXCleanReads.pl --inputdir $inputdir --cleanpath $cleanpath --sample1 $sample1 --sample2 $sample2 --threads $threads  
EOF
	close PBS;

	my $qsub = `qsub $inputdir/$pbspath/$study/fixFastxReadsPairs_$runid.pbs`;
	if ($qsub ne ""){
		my $lc = localtime;
		print "$lc\tThe task of fixFastxReadsPairs_$runid is submitted in $qsub\n"
	}
}


