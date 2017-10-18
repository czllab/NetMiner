#############################################################################################
# Program Name: fixFASTXReadsPairs.pl
# This program repair paired-end cleaned reads produced by FASTX software
# Author: Yu Hua
# Created: 2014-7-16
#############################################################################################

#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my ($cleanpath,$samplereads1file,$samplereads2file,$threads,$help,%header);
my $usage = <<USAGE;
Usage: perl $0 --cleanpath \$cleanpath --samplereads1file \samplereads1file --samplereads2file \samplereads2file --threads \$threads --help
Options:
--cleanpath: clean file path
--samplereads1file: R1 reads file
--samplereads2file: R2 reads file
--threads: the number of threads
--help: displaying help information
USAGE
;

GetOptions(
	"cleanpath=s" => \$cleanpath,
	"samplereads1file=s" => \$samplereads1file,
	"samplereads2file=s" => \$samplereads2file,
	"threads=s" => \$threads,
	"help!" => \$help,
) or die "$usage\n";

if ($help){
	print "$usage\n";
}

$samplereads1file =~ /(.*)\/(.*)\/(.*).clean_PE_.\.fastx\.fastq\.gz/;
my $runid = $3;
my $study = $2;
my @samplereadsfiles = ($samplereads1file,$samplereads2file);
my $flag = 0;
foreach my $samplereadsfile (@samplereadsfiles){
	my $pairedsamplereads = $samplereadsfile;
	my $unpairedsamplereads = $samplereadsfile;
	$pairedsamplereads =~ s/\.fastx\.fastq\.gz/\.paired_fastx\.fastq.gz/;
	$unpairedsamplereads =~ s/\.fastx\.fastq\.gz/\.unpaired_fastx\.fastq.gz/;
	my $totalsize = -s "$samplereadsfile";
	my $pairedsize = -s "$pairedsamplereads";
	my $unpairedsize = -s "$unpairedsamplereads";
	if(!-e "$pairedsamplereads" || ($pairedsize + $unpairedsize) < $totalsize || $pairedsize < 10000){
		$flag = 1;
		last;
	}
}

if ($flag == 1){
	foreach my $samplereadsfile (@samplereadsfiles){
		open(SAM,"gzip -dc $samplereadsfile|") or die "couldn't open the file $samplereadsfile: $!\n";
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

my ($unpairednum,$unpairedbase,$pairedbase,$pairednum);
foreach my $samplereadsfile (@samplereadsfiles){
	my $pairedsamplereads = $samplereadsfile;
	my $unpairedsamplereads = $samplereadsfile;
	$pairedsamplereads =~ s/\.fastx\.fastq\.gz/\.paired_fastx\.fastq.gz/;
	$unpairedsamplereads =~ s/\.fastx\.fastq\.gz/\.unpaired_fastx\.fastq.gz/;
	my $totalsize = -s "$samplereadsfile";
	my $pairedsize = -s "$pairedsamplereads";
	my $unpairedsize = -s "$unpairedsamplereads";
	if(!-e "$pairedsamplereads" || ($pairedsize + $unpairedsize) < $totalsize || $pairedsize < 10000){
		if (-e "$pairedsamplereads"){
			unlink $pairedsamplereads;
			unlink $unpairedsamplereads;
		}
		open(PAIR,">$pairedsamplereads") or die "couldn't open the file $pairedsamplereads:$!\n";
		open(NONPAIR,">$unpairedsamplereads") or die "couldn't open the file $unpairedsamplereads:$!\n";
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
						print PAIR "$readInfo\n$sequence\n$readID\n$readQual\n";
					}else{
						print NONPAIR "$readInfo\n$sequence\n$readID\n$readQual\n";
					}
				}else{
					if($header{$readInfo}==2){
						print PAIR "$readInfo\n$sequence\n$readID\n$readQual\n";
					}else{
						print NONPAIR "$readInfo\n$sequence\n$readID\n$readQual\n";
					}
				}
			}
		}
		close PAIR;
		close NONPAIR;
		close READSFILE;
		system("pigz -p $threads $pairedsamplereads") == 0 or die "system calling pigz program failed\n";
		system("pigz -p $threads $unpairedsamplereads") == 0 or die "system calling pigz program failed\n";
	}
}
