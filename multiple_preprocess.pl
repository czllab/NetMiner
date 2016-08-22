#########################################################################################################################################
#	> File Name: multiple_preprocess.pl
#	> This program includes three functions 
#	1) Converting SRA file to FASTQ file
#	2) Reporting reads quality of FASTQ file 
#	3) Trimming and filtering reads 
#	> Author: Hua Yu
#	> Mail: huayu@genetics.ac.cn 
#	Created Time: 2014-07-12
##########################################################################################################################################

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($runinfofile,$inputpath,$outputpath,$species,$pbspath,$rawfastqpath,$fastqcpath,$cleanpath,$trimmomatic,$queue,$threads,$help,%runHash);
my $QUAL = 10;
my $PERCENT = 50;
my $usage = <<USAGE;
Usage:
	perl $0 --runinfofile \$runinfofile --inputpath \$inputpath --outputpath \$outputpath --pbspath \$pbspath --rawfastqpath \$rawfastqpath --fastqcpath \$fastqcpath --cleanpath \$cleanpath --trimmomatic --queue \$queue --threads \$threads
Options:
	--runinfofile: NCBI run information file
	--inputpath: SRA file path 
	--outputpath: FASTQ file path
	--pbspath: PBS scripts path
	--rawfastqpath: fastq-dump result file path
	--fastqcpath: fastqc result file path
	--cleanpath: clean result file path
	--trimmomatic: whether or not use thrimmomatic software for reads trimming and filtering
	--queue: PBS queue
	--threads: the number of threads
	--help: display help information
USAGE
;

GetOptions(
	"runinfofile=s" => \$runinfofile,
	"inputpath=s" => \$inputpath,
	"outputpath=s" => \$outputpath,
	"pbspath=s" => \$pbspath,
	"rawfastqpath=s" => \$rawfastqpath,
	"fastqcpath=s" => \$fastqcpath,
	"cleanpath=s" => \$cleanpath,
	"trimmomatic!" => \$trimmomatic,
	"queue=s" => \$queue,
	"threads=i" => \$threads,
	"help!" => \$help,
) or die ($usage);

if ($help){
	print "$usage";
	exit(1);
}

my @runs=`find $inputpath -name "*.sra"`;
print "$_\n" foreach (@runs);

for my $run (@runs){
	chomp $run;
    my @fields=split/\//,$run;
    $runHash{$fields[-1]} = $fields[-2];  
}
my ($machine,$readLength) = &getMachineAndReadsLengthInfo(\$runinfofile);

&createJob(\%runHash,\$inputpath,\$outputpath,\$pbspath,\$rawfastqpath,\$fastqcpath,\$cleanpath,$machine,$readLength,\$trimmomatic,\$threads);

# creating preprocess jobs
sub createJob {
	my ($runHash,$inputpath,$outputpath,$pbspath,$rawfastqpath,$fastqcpath,$cleanpath,$machine,$readLength,$trimmomatic,$threads) = @_;
	my(%sra2fastqmarks,$in_p1,$in_p2,$type,$out_p1,$out_p2,$out_up1,$out_up2,$operators,$step);
	
	# Converting SRA file to FASTQ file
	my $jobname = "sra2fastq";
	foreach my $run (keys %$runHash){
		my $command = "";
		if(!-e "$$outputpath/$$rawfastqpath/$runHash{$run}"){
			mkpath("$$outputpath/$$rawfastqpath/$runHash{$run}",0644);
			if($@){
				print "Make path $$outputpath/$$rawfastqpath/$runHash{$run} failed:\n$@";
				exit(1);
			}
		}
		if(!-e "$$outputpath/$$fastqcpath/$runHash{$run}"){
			mkpath("$$outputpath/$$fastqcpath/$runHash{$run}",0644);
			if($@){
				print "Make path $$outputpath/$$fastqcpath/$runHash{$run} failed:\n$@";
				exit(1);
			}
		}
		if(!-e "$$outputpath/$$cleanpath/$runHash{$run}"){
			mkpath("$$outputpath/$$cleanpath/$runHash{$run}",0644);
			if($@){
				print "Make path $$outputpath/$$cleanpath/$runHash{$run} failed:\n$@";
				exit(1);
			}
		}
		if(!-e "$$outputpath/$$pbspath/$runHash{$run}"){
			mkpath("$$outputpath/$$pbspath/$runHash{$run}",0755);
			if($@){
				print "Make path $$outputpath/$$pbspath/$runHash{$run} failed:\n$@";
				exit(1);
			}
		}
		my $runid = $run;
		$runid =~ s/\.sra//;
		if (!-e "$$outputpath/$$rawfastqpath/$$runHash{$run}/${runid}_1.fastq.gz" ){
			$command .= "fastq-dump --split-files $$inputpath/$$runHash{$run}/$run -O $$outputpath/$$rawfastqpath/$$runHash{$run} --gzip\n";
			$command .= "fastqc -o $$outputpath/$$fastqcpath/$$runHash{$run} -t $$threads --extract -f fastq $$outputpath/$$rawfastqpath/$$runHash{$run}/$runid*.gz";
		}elsif(!-d "$$outputpath/$$fastqcpath/$$runHash{$run}/${runid}_1_fastqc"){
			$command = "fastqc -o $$outputpath/$$fastqcpath/$$runHash{$run} -t $$threads --extract -f fastq $$outputpath/$$rawfastqpath/$$runHash{$run}/$runid*.gz";
		}
		if ($command ne ""){
			$sra2fastqmarks{$runid} = 1;
			open OUT,">$$outputpath/$$pbspath/$runHash{$run}/$jobname\_$runid".".pbs";
			print OUT <<EOF;
#PBS -N $jobname\_$runid
#PBS -j oe
#PBS -l nodes=1:ppn=$$threads
#PBS -q $queue
#PBS -o $$outputpath/$$pbspath/$runHash{$run}/$jobname\_$runid.log

cd \$PBS_O_WORKDIR
date
$command
date
EOF
			close OUT;
		}else{
			$sra2fastqmarks{$runid} = 0;
		}
	}
	&submitJob($runHash,$outputpath,$pbspath,$jobname,\%sra2fastqmarks,$species,$queue);
	while (1){
		sleep (60);
		my $qstat = `qstat -uhyu | grep $jobname | grep $queue`;
		if ($qstat eq ""){
			last;
		}
	}
	
	# reads trimming and filtering
	my $jobname = "fastqfilter";
	my %marksfastqfilter;
	foreach my $run (keys %$runHash){
		my $runid = $run;
		$runid =~ s/\.sra//;
		print "$$outputpath/$$rawfastqpath/$$runHash{$run}/ -name $runid\_1.fastq.gz\n";
		$in_p1 = `find $$outputpath/$$rawfastqpath/$$runHash{$run} -name "$runid\_1.fastq.gz"`;
		chomp $in_p1;
		$in_p2 = `find $$outputpath/$$rawfastqpath/$$runHash{$run} -name "$runid\_2.fastq.gz"`;
		chomp $in_p2;
		if ($$trimmomatic){
			if(!-e "$$outputpath/$$cleanpath/$$runHash{$run}/$runid".".clean_PE_1.fastq.gz" && !-e "$$outputpath/$$cleanpath/$$runHash{$run}/$runid".".clean_SE_1.fastq.gz"){
				if(-e $in_p1){
					if(-e $in_p2){
						my $minLen = int (($$readLength{$runid}*7/20)+0.5);
						$type="PE";
						$out_p1="$$outputpath/$$cleanpath/$$runHash{$run}/$runid".".clean_PE_1.fastq.gz";
						$out_p2="$$outputpath/$$cleanpath/$$runHash{$run}/$runid".".clean_PE_2.fastq.gz";
						$out_up1="$$outputpath/$$cleanpath/$$runHash{$run}/$runid".".clean_UP_1.fastq.gz";
						$out_up2="$$outputpath/$$cleanpath/$$runHash{$run}/$runid".".clean_UP_2.fastq.gz";
						if ($$machine{$runid}=~/Illumina HiSeq/){
							$step="ILLUMINACLIP:/public/software/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:3 LEADING:2 TRAILING:2 SLIDINGWINDOW:3:2 MINLEN:$minLen TOPHRED33";
						}else{
							$step="ILLUMINACLIP:/public/software/Trimmomatic-0.32/adapters/TruSeq2-PE.fa:2:30:3 LEADING:2 TRAILING:2 SLIDINGWINDOW:3:2 MINLEN:$minLen TOPHRED33";
						}		
					}else{
						my $minLen = int (($$readLength{$runid}*7/10)+0.5);
						$type="SE"; 
						$out_p1="$$outputpath/$$cleanpath/$$runHash{$run}/$runid".".clean_SE_1.fastq.gz";
						if ($$machine{$runid}=~/Illumina HiSeq/){
							$step="ILLUMINACLIP:/public/software/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:3 LEADING:2 TRAILING:2 SLIDINGWINDOW:3:2 MINLEN:$minLen TOPHRED33";
						}else{
							$step="ILLUMINACLIP:/public/software/Trimmomatic-0.32/adapters/TruSeq2-SE.fa:2:30:3 LEADING:2 TRAILING:2 SLIDINGWINDOW:3:2 MINLEN:$minLen TOPHRED33";
						}
					}
				}else{
					print "couldn't find the fastq file of $$outputpath/$$rawfastqpath/$$runHash{$run}/${runid}_1.fastq.gz\n";
				}
				if(-e $in_p1 && $in_p2 eq ""){
					$operators = "$in_p1 $out_p1 $step";
				}elsif(-e $in_p1 && -e $in_p2){
					$operators = "$in_p1 $in_p2 $out_p1 $out_up1 $out_p2 $out_up2 $step"; 
				}else{
					$operators = "";
				}
			}else{
				$operators = "";
			}
			if($operators ne ""){
				open OUT,">$$outputpath/$$pbspath/$runHash{$run}/$jobname\_$runid".".pbs";
				print OUT <<EOF;
#PBS -N $jobname\_$runid
#PBS -j oe
#PBS -l nodes=1:ppn=$$threads
#PBS -q $queue
#PBS -o $$outputpath/$$pbspath/$runHash{$run}/$jobname\_$runid.log

cd \$PBS_O_WORKDIR

date   
/public/software/jdk1.7.0_45/bin/java -jar /public/software/Trimmomatic-0.32/trimmomatic-0.32.jar $type -threads $$threads $operators
for file in `find $$outputpath/$$cleanpath/$$runHash{$run} -name "${runid}*clean*.fastq.gz"`
	do 
		sample=\${file%.fastq.gz}
		sample=\${sample##\/*\/}
		gunzip -c \$file | fastq_quality_filter -v -q $QUAL -p $PERCENT -z -Q 33 -o $$outputpath/$$cleanpath/$$runHash{$run}/\$sample.fastx.fastq.gz
	done
fastqc -o $$outputpath/$$fastqcpath/$$runHash{$run} -t $$threads --extract -f fastq $$outputpath/$$cleanpath/$$runHash{$run}/$runid*.fastx.fastq.gz
date
EOF
				close OUT;
				$marksfastqfilter{$runid} = 1;
			}else{
				my @samples = `find $$outputpath/$$cleanpath/$$runHash{$run} -name "${runid}*clean*.fastx.fastq.gz"`;
				my $flagsample = $samples[0];
				chomp $flagsample;
				if(!-e $flagsample){
					open OUT,">$$outputpath/$$pbspath/$runHash{$run}/$jobname\_$runid".".pbs";
					print OUT <<EOF;
#PBS -N $jobname\_$runid
#PBS -j oe
#PBS -l nodes=1:ppn=$$threads
#PBS -q $queue
#PBS -o $$outputpath/$$pbspath/$runHash{$run}/$jobname\_$runid.log

cd \$PBS_O_WORKDIR
date
for file in `find $$outputpath/$$cleanpath/$$runHash{$run} -name "${runid}*clean*.fastq.gz"`
	do 
		sample=\${file%.fastq.gz}
		sample=\${sample##\/*\/}
		gunzip -c \$file | fastq_quality_filter -v -q $QUAL -p $PERCENT -z -Q 33 -o $$outputpath/$$cleanpath/$$runHash{$run}/\$sample.fastx.fastq.gz
	done
fastqc -o $$outputpath/$$fastqcpath/$$runHash{$run} -t $$threads --extract -f fastq $$outputpath/$$cleanpath/$$runHash{$run}/$runid*.fastx.fastq.gz
date
EOF
					close OUT;
					$marksfastqfilter{$runid} = 1;
				}else{
					$marksfastqfilter{$runid} = 0;
				}	
			}
		}
	}
	&submitJob($runHash,$outputpath,$pbspath,$jobname,\%marksfastqfilter,$species,$queue);
}

# getting the machines and reads lengths of RNA-seq samples
sub getMachineAndReadsLengthInfo{
	my ($runInfoFile) = @_;
	my (@fieldnames,%machine,%readLength,%nocol);
	open(RUNFILE,"<$$runInfoFile") or die "couldn't open the file $runinfofile: $!\n";
	while (<RUNFILE>){
		my $line = $_;
		$line =~ s/^\s+|\s+$//;
		if ($line =~ /^Run/){
			@fieldnames = split /,/, $line;
			for(my $i = 0; $i <= $#fieldnames; $i++){
				if($fieldnames[$i] eq "avgLength"){
					$nocol{"avgLength"} = $i;
				}
				if ($fieldnames[$i] eq "Model"){
					$nocol{"Model"} = $i;
				} 
			}
		}else{
			my @valueinfo = split /,/,$line;
			$machine{$valueinfo[0]} = $valueinfo[$nocol{"Model"}];
			$readLength{$valueinfo[0]} = $valueinfo[$nocol{"avgLength"}];
		}
	}
	return (\%machine,\%readLength);
}

# submit PBS jobs
sub submitJob{
	my ($runHash,$outputpath,$pbspath,$jobname,$marks,$species,$queue) = @_;
	foreach my $run (keys %$runHash){
		my $taskNum =`qstat -uhyu | grep $jobname | grep $queue | wc -l`; 
		while($taskNum > 100){
			print "The num of task remaining $taskNum\n";
			sleep 30;
			print `date`;
			$taskNum = `qstat -uhyu | grep $jobname | grep $queue | wc -l`;
		}
		my $qsub;
		my $runid = $run;
		$runid =~ s/\.sra//;
		print "$runid\n";
		if ($$marks{$runid} == 1){
			$qsub = `qsub $$outputpath/$$pbspath/$runHash{$run}/$jobname\_$runid.pbs`;
		}else{
			print "The task $jobname\_$runid is successfully finished\n";
		}
		my $lc = localtime;
		
		if ($qsub ne ""){
			print "$lc\nThe task of $jobname\_$runid was submitted in $qsub\n";
		}
	}
}
