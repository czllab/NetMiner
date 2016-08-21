#############################################################################################
# Program Name: fixFASTXReadsPairs.pl
# This program repair the cleaned reads of FASTX software
# Author: Yu Hua
# Created: 2014-6-26
#############################################################################################

#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my ($inputdir,$cleanpath,$sample1,$sample2,$threads,%header,$help);
my $usage = <<USAGE;
Usage: perl $0 --inputdir \$inputdir --cleanpath \$cleanpath --sample1 \$sample1 --sample2 \$sample2 --threads \$threads --help
Options:
--inputdir: inputdir file directory
--cleanpath: clean file path
--sample1: R1 paired reads file name
--sample2: R2 paired reads file name
--threads: the number of threads
--help: display help information
USAGE
;

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"cleanpath=s" => \$cleanpath,
	"sample1=s" => \$sample1,
	"sample2=s" => \$sample2,
	"threads=s" => \$threads,
	"help!" => \$help,
) or die "$usage\n";

if ($help){
	print "$usage\n";
}

$sample1 =~ /(.*)\/(.*)\/(.*).clean_PE_.\.fastx\.fastq\.gz/;
my $runid = $3;
my $study = $2;
my @samples = ($sample1,$sample2);
my $flag = 0;
foreach my $sample (@samples){
	my $pairedSample = $sample;
	my $unpairedSample = $sample;
	$pairedSample =~ s/\.fastx\.fastq\.gz/\.paired_fastx\.fastq/;
	$unpairedSample =~ s/\.fastx\.fastq\.gz/\.unpaired_fastx\.fastq/;
	my $sizeTotal = -s "$sample";
	my $sizePaired = -s "$pairedSample.gz";
	my $sizeUnpaired = -s "$unpairedSample.gz";
	if(!-e "$pairedSample.gz" || ($sizePaired + $sizeUnpaired) < $sizeTotal || $sizePaired < 10000){
		$flag = 1;
		last;
	}
}

if ($flag == 1){
	foreach my $sample (@samples){
		open(SAM,"gzip -dc $sample|") or die "Couldn't open the file $sample: $!\n";
		while(<SAM>){
			my $line = $_;
			chomp $line;
			if($line =~ /^\@.RR\d+(.*)\s+(.*)/){
				if($line =~/(.*)\s+length=\d+/){
					$header{$1} += 1;
				}else{
					$header{$line} += 1;
				}
			}
		}
		close SAM;
	}
}
# print "$_=>$header{$_}\n" foreach (keys %header);

my ($unpairednum,$unpairedbase,$pairedbase,$pairednum);
foreach my $sample (@samples){
	open(READSFILE,"gzip -dc $sample|");
	my $pairedSample = $sample;
	my $unpairedSample = $sample;
	$pairedSample =~ s/\.fastx\.fastq\.gz/\.paired_fastx\.fastq/;
	$unpairedSample =~ s/\.fastx\.fastq\.gz/\.unpaired_fastx\.fastq/;
	my $sizeTotal = -s "$sample";
	my $sizePaired = -s "$pairedSample.gz";
	my $sizeUnpaired = -s "$unpairedSample.gz";
	if(!-e "$pairedSample.gz" || ($sizePaired + $sizeUnpaired) < $sizeTotal || $sizePaired < 10000){
		if (-e "$pairedSample.gz"){
			unlink "$pairedSample.gz";
			unlink "$unpairedSample.gz";
		}
		open(SAM,">$pairedSample") or die "Couldn't open the file $pairedSample:$!\n";
		open(NOSAM,">$unpairedSample") or die "Couldn't open the file $unpairedSample:$!\n";
		while(defined(my $readInfo = <READSFILE>)){
			chomp $readInfo;
			my $sequence = <READSFILE>;
			chomp $sequence;
			my $readID = <READSFILE>;
			chomp $readID;
			my $readQual = <READSFILE>;
			chomp $readQual;
			if($readInfo =~ /^\@.RR\d+(.*)\s+(.*)/ && length($sequence) == length($readQual)){
				if($readInfo =~ /(.*)\s+length=\d+/){
					if($header{$1}==2){
						$pairednum += 1;
						$pairedbase += length($sequence);
						print SAM "$readInfo\n$sequence\n$readID\n$readQual\n";
					}else{
						$unpairednum += 1;
						$unpairedbase += length($sequence);
						print NOSAM "$readInfo\n$sequence\n$readID\n$readQual\n";
					}
				}else{
					if($header{$readInfo}==2){
						$pairednum += 1;
						$pairedbase += length($sequence);
						print SAM "$readInfo\n$sequence\n$readID\n$readQual\n";
					}else{
						$unpairednum += 1;
						$unpairedbase += length($sequence);
						print NOSAM "$readInfo\n$sequence\n$readID\n$readQual\n";
					}
				
				}
			}
		}
		close SAM;
		close NOSAM;
		close READSFILE;
		system("pigz -p $threads $pairedSample") == 0 or die "system calling failed\n";
		system("pigz -p $threads $unpairedSample") == 0 or die "system calling failed\n";
	}
}

open(OUT,">>$inputdir/$cleanpath/fastxCleanInfo.csv");
print OUT "$study,$runid,$pairednum,$pairedbase,$unpairednum,$unpairedbase\n";
close OUT;
