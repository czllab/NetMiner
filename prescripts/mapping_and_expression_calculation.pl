######################################################################################################################
#	> File Name: mapping_and_expression_calculation.pl
#	> This program aligns the cleaned reads to reference genome and calculates the gene expression
#       > abundance by Tophat, htseq-count and cufflinks software
#	> Author: Hua Yu
#	> Mail: huayu@genetics.ac.cn 
#	> Created Time: 2014-08-02
######################################################################################################################

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputpath,$outputpath,$indexpath,$bowtiebase,$tophatpath,$picard,$htseq,$cufflinks,$htseqpath,$cufflinkspath,$nomultiplehits,$removeduplicate,$runinfofile,$gfffile,$coveragedir,$queue,$threads,$refgenome,$pbspath,$summarydir,$help);
my $usage = <<"USAGE"
Usage:
$0 --inputpath \$inputpath --outputpath \$outputpath --indexpath \indexpath --bowtiebase \$bowtiebase --tophatpath \$tophatpath --picardpath \$picardpath --htseqpath \$htseqpath --cufflinkspath \$cufflinkspath --nomultiplehits --removeduplicate --runinfofile \$runinfofile --gfffile \$gfffile --threads \$threads --refgenome \$refgenome --pbspath \$pbspath --help
Options:
--inputpath: input file path
--outputpath: output file path
--indexpath: bowtie index file path
--bowtiebase: bowtie index file basename
--tophatpath: tophat output file path
--picardpath: picard program path
--htseqpath: htseq-count output file path
--cufflinkspath: cufflinks output file path
--nomultiplehit: whether or not contain the multiple hits
--removeduplicate: whether or not remove the duplicate reads
--runinfofile: sequencing run information file
--gfffile: gff annotation file
--queue: computing queue
--threads: the number of threads
--refgenome: reference genome path
--pbspath: pbs script path
USAGE
;


GetOptions(
	"inputpath|i=s" => \$inputpath,
	"outputpath|o=s" => \$outputpath,
	"indexpath|index=s" => \$indexpath,
	"bowtiebase|base=s" => \$bowtiebase,
	"tophatpath=s" => \$tophatpath,
	"picardpath=s" => \$picardpath, 
	"htseqpath=s" => \$htseqpath,
	"cufflinkspath=s" => \$cufflinkspath,
	"nomultiplehits!" => \$nomultiplehits,
	"removeduplicate!" => \$removeduplicate,
	"runinfofile|runinfo=s" => \$runinfofile,
	"gfffile|gff=s" => \$gfffile,
	"coveragedir=s" => \$coveragedir,
	"queue=s" => \$queue,
	"threads=i" => \$threads,
	"refgenome=s" => \$refgenome,
	"pbspath=s" => \$pbspath,
	"help!" => \$help,
) or die "$usage";

if($help){
	print "$usage";
	exit(1);
}
my @samples = `find $inputpath -name "*E_1.fastx.fastq.gz"`;

my %InsertSize = &getMateInnerDist($runinfofile);

&createJob(\$picardpath,\@samples,\%InsertSize,\$threads,\$outputpath,\$indexpath,\$bowtiebase,\$tophatpath,\$htseqpath,\$cufflinkspath,\$pbspath,\$refgenome,\$gfffile,\$nomultiplehits,\$removeduplicate);

sub getMateInnerDist{
	my ($runinfofile) = @_;
	my (%nocols,%InsertSize);
	open RUN, "<$runinfofile" or die "couldn't open the file $runinfofile:$!\n";
	while(<RUN>){
		my $line = $_;
		$line =~ s/^\s+|\s+$//;
		if($line =~ /^[E|S|D]RR/){
			my @values = split /,/,$line;
			$InsertSize{$values[$nocols{"Run"}]} = $values[$nocols{"InsertSize"}];
		}else{
			my @fieldNames = split /,/,$line;
			for (my $i = 0; $i <= $#fieldNames; $i++){
				if ($fieldNames[$i] eq "Run"){
					$nocols{"Run"} = $i;
				}
				if ($fieldNames[$i] eq "InsertSize"){
					$nocols{"InsertSize"} = $i;
				}
			}
		}			
	}
	return %InsertSize;
}

sub submitJob{
	my ($samples,$outputpath,$pbspath,$marks) = @_;
	foreach my $sample (@$samples){
		my $taskNum =`qstat -uhyu | grep geneExpression | wc -l`; 
		while($taskNum > 100){
			print "The num of task remaining $taskNum\n";
			sleep 30;
			print `date`;
			$taskNum = `qstat -uhyu | grep geneExpression | wc -l`;
		}
		$sample =~ /.*\/(.*)\/(.*)\.clean.*fastq.gz/;
		my $study = $1;
		my $runid = $2;
		if($$marks{$runid} == 1){
			my $qsub = `qsub $$outputpath/$$pbspath/$study/geneExpressionAbundanceCalculation_$runid.pbs`;
			my $lc = localtime;
			if ($qsub ne ""){
				print "$lc\nThe task of geneExpressionAbundanceCalculation_$runid.pbs was submitted in $qsub\n";
			}
		}
	}
}

sub createJob{
	my ($picardpath,$samples,$InsertSize,$threads,$outputpath,$indexpath,$bowtiebase,$tophatpath,$htseqpath,$cufflinkspath,$pbspath,$refgenome,$gfffile,$nomultiplehits,$removeduplicate) = @_;
	my ($command,%marks);
	foreach my $sample (@$samples){
		$command="";
		chomp $sample;
		$sample =~ /.*\/(.*)\/(.*)\.clean.*fastq\.gz/;
		my $study = $1;
		my $runid = $2;
		my $sample_p1 = $sample;
		my $unsample_p1 = $sample;
		$sample_p1 =~ s/\.fastx\.fastq\.gz/\.paired_fastx\.fastq\.gz/;
		$unsample_p1 =~ s/\.fastx\.fastq\.gz/\.unpaired_fastx\.fastq\.gz/;
		my $sample_p2 = $sample_p1;
		my $unsample_p2 = $unsample_p1;
		$sample_p2 =~ s/E_1\.paired_fastx\.fastq\.gz/E_2\.paired_fastx\.fastq\.gz/;
		$unsample_p2 =~ s/E_1\.unpaired_fastx\.fastq\.gz/E_2\.unpaired_fastx\.fastq\.gz/;
		my $upsample_p1 = $sample;
		$upsample_p1 =~ s/.E\_1\.fastx\.fastq\.gz/UP\_1\.fastx\.fastq\.gz/;
		my $upsample_p2 = $upsample_p1;
		$upsample_p2 =~ s/UP\_1\.fastx\.fastq\.gz/UP\_2\.fastx\.fastq\.gz/;
		if ((-s $unsample_p1) == 63){
			unlink $unsample_p1;
		}
		if ((-s $unsample_p2) == 63){
			unlink $unsample_p2;
		}
		if ((-s $upsample_p1) == 63){
			unlink $upsample_p1;
		}
		if ((-s $upsample_p2) == 63){
			unlink $upsample_p2;
		}
		print "$upsample_p1\n$upsample_p2\n";
		if(!-e "$$outputpath/$$tophatpath/$study/$runid"){
			mkpath("$$outputpath/$$tophatpath/$study/$runid",0644);
			if($@){
				print "Make path $$outputpath/$$tophatpath/$study/$runid failed:\n$@";
				exit(1);
			}
		}
		if(!-e "$$outputpath/$$htseqpath/$study"){
			mkpath("$$outputpath/$$htseqpath/$study",0644);
			if($@){
				print "Make path $$outputpath/$$htseqpath/$study failed:\n";
				exit(1);
			}
		}
		if(!-e "$$outputpath/$$cufflinkspath/$study/$runid/"){
			mkpath("$$outputpath/$$cufflinkspath/$study/$runid/",0644);
			if($@){
				print "Make path $$outputpath/$$cufflinkspath/$study/$runid/ failed:\n";
				exit(1);
			}
		}
		
		if(!-e "$$outputpath/$$pbspath/$study/"){
			mkpath("$$outputpath/$$pbspath/$study/",0644);
			if($@){
				print "Make path $$outputpath/$$pbspath/$study/ failed:\n";
				exit(1);
			}
		}
		
		if(!-e "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam" || -z "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam"){
			if(-e $sample_p2){	
				if(-e $unsample_p1 && !-z $unsample_p1 && -e $unsample_p2 && !-z $unsample_p2 && -e $upsample_p1 && !-z $upsample_p1 && -e $upsample_p2 && !-z $upsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p1,$unsample_p2,$upsample_p1,$upsample_p2\n";
				}elsif(-e $unsample_p1 && !-z $unsample_p1 && -e $unsample_p2 && !-z $unsample_p2 && -e $upsample_p1 && !-z $upsample_p1){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid  --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p1,$unsample_p2,$upsample_p1\n";
				}elsif(-e $unsample_p1 && !-z $unsample_p1 && -e $unsample_p2 && !-z $unsample_p2 && -e $upsample_p2 && !-z $upsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p1,$unsample_p2,$upsample_p2\n";
				}elsif(-e $unsample_p1 && !-z $unsample_p1 && -e $upsample_p1 && !-z $upsample_p1 && -e $upsample_p2 && !-z $upsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p1,$upsample_p1,$upsample_p2\n";
				}elsif(-e $unsample_p2 && !-z $unsample_p2 && -e $upsample_p1 && !-z $upsample_p1 && -e $upsample_p2 && !-z $upsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p2,$upsample_p1,$upsample_p2\n";
				}elsif(-e $unsample_p1 && !-z $unsample_p1 && -e $unsample_p2 && !-z $unsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p1,$unsample_p2\n";
				}elsif(-e $unsample_p1 && !-z $unsample_p1 && -e $upsample_p1 && !-z $upsample_p1){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p1,$upsample_p1\n";
				}elsif(-e $unsample_p1 && !-z $unsample_p1 && -e $upsample_p2 && !-z $upsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p1,$upsample_p2\n";
				}elsif(-e $unsample_p2 && !-z $unsample_p2 && -e $upsample_p1 && !-z $upsample_p1){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p2,$upsample_p1\n";
				}elsif(-e $unsample_p2 && !-z $unsample_p2 && -e $upsample_p2 && !-z $upsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p2,$upsample_p2\n";
				}elsif(-e $upsample_p1 && !-z $upsample_p1 && -e $upsample_p2 && !-z $upsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$upsample_p1,$upsample_p2\n";
				}elsif(-e $unsample_p1 && !-z $unsample_p1){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p1\n";
				}elsif(-e $unsample_p2 && !-z $unsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$unsample_p2\n";
				}elsif(-e $upsample_p1 && !-z $upsample_p1){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$upsample_p1\n";
				}elsif(-e $upsample_p2 && !-z $upsample_p2){
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2,$upsample_p2\n";
				}else{
					$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --mate-inner-dist $$InsertSize{$runid} --rg-id $runid --rg-sample $runid --no-mixed --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample_p1 $sample_p2\n";
				}
			}else{
				$command .= "tophat2 -p $$threads --library-type fr-unstranded --solexa-quals --rg-id $runid --rg-sample $runid --no-convert-bam -o $$outputpath/$$tophatpath/$study/$runid $$indexpath/$$bowtiebase $sample\n";
			}
		}	
		if(!-e "$$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.sam" || -z "$$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.sam" || !-e "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam" || -z "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam"){
			if($$nomultiplehits){
				$command .= "grep -v -E -w 'NH:i:2|NH:i:3|NH:i:4|NH:i:5|NH:i:6|NH:i:7|NH:i:8|NH:i:9|NH:i:10|NH:i:11|NH:i:12|NH:i:13|NH:i:14|NH:i:15|NH:i:16|NH:i:17|NH:i:18|NH:i:19|NH:i:20' $$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam > $$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.sam\n";
			}
		}
		if(!-e "$$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.unique.sam" || -z "$$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.unique.sam" || !-e "$$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.sam" || -z "$$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.sam" || !-e "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam" || -z "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam"){
			if($$removeduplicate){
				$command .= "/public/software/jdk1.7.0_45/bin/java -Xmx15g -jar $$picard/MarkDuplicates.jar I=$$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.sam O=$$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.unique.sam METRICS_FILE=$$outputpath/$$tophatpath/$study/$runid/${runid}.metricsFile VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true ASSUME_SORTED=true\n";
			}
		}
		if (!-e "$$outputpath/$$htseqpath/$study/${runid}.htseq_count.txt" || -z "$$outputpath/$$htseqpath/$study/${runid}.htseq_count.txt" || !-e "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam" || -z "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam"){
			if($htseq){
				$command .= "htseq-count -s no -i ID -t gene -q $$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.unique.sam $$gfffile > $$outputpath/$$htseqpath/$study/${runid}.htseq_count.txt\n";
			}
		}
		if (!-e "$$outputpath/$$cufflinkspath/$study/$runid/transcripts.gtf" || !-e "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam" || -z "$$outputpath/$$tophatpath/$study/$runid/accepted_hits.sam"){
			if($cufflinks){
				$command .= "cufflinks -p $$threads --GTF $$gfffile -o $$outputpath/$$cufflinkspath/$study/$runid $$outputpath/$$tophatpath/$study/$runid/accepted_hits_NHi1.unique.sam\n";
			}
		}
		
		if($command ne ""){
			$marks{$runid} = 1;
			open OUT,">$$outputpath/$$pbspath/$study/mapping_and_expression_cal_$runid.pbs";
			print OUT <<"EOF";
#PBS -N mapping_and_expression_cal_$runid
#PBS -j oe
#PBS -l nodes=1:ppn=$$threads
#PBS -q $queue
#PBS -o $$outputpath/$$pbspath/$study/mapping_and_expression_cal_$runid.log

cd \$PBS_O_WORKDIR
date
$command
date
EOF
			close OUT; 
		}else{
			$marks{$runid} = 0;
		}
	}
	&submitJob($samples,$outputpath,$pbspath,\%marks);
}
		
