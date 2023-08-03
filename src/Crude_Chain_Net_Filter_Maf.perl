#!/usr/bin/env perl

# Michael Hiller, 2013
# for the new genome alignment pipeline of psl filtering, chain patching, chain and net synteny filtering 
# This script performs the chain filtering, then nets the filtered chains and filters the nets again. 
# 2020: Can now be run on delta (our compute cluster) by running the netClass step (this requires a SQL database on a genome browser server) via ssh on genome (our genome browser server)
# The genome browser server must be ssh-able without password and must be specified in the variable $dbHost below

# FilterChains_Net_FilterNets modified to remove optional netClass and add conversion to maf

use strict;
use warnings;
use diagnostics;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use POSIX;
use Sys::Hostname;

$| = 1;		# == fflush(stdout)
my $verbose = 0;
my $inChain = "";				# filename + path
my $tTB = "";				# filename + path
my $qTB = "";				# filename + path
my $tSizes = "";				# filename + path
my $qSizes = "";
my $outPrefix = "";				# filename + path
my $Tassembly="";				# ref and query assembly
my $Qassembly="";
my $minScore="";				# can be an array, comma-separated
my $minTsize="";
my $minQsize="";
my $keepSynNetsWithScore = INT_MAX;      # keep nets classified as syn if the are above this score threshold
my $keepInvNetsWithScore = INT_MAX;      # keep nets classified as inv if the are above this score threshold

my $usage = "usage: $0 inChainFile referenceTwoBit queryTwoBit referenceSizes querySizes outPrefix minScores [-v|verbose -keepSynNetsWithScore int,-keepInvNetsWithScore int, -minSizeT comma-separated-string -minSizeQ comma-separated-string] \n
	inChain, referenceTwoBit, queryTwoBit, referenceSizes and querySizes should just be filenames including the path to the file
	inChain can be gzipped (DO NOT UNZIP)
	outPrefix should be the beginning of the output filenames including the  path to the file
	the to-be-produced files outFilteredChain and outFilteredNet must NOT have the .gz suffix ($0 will gzip them afterwards)
	minScores  can be a comma-separated list
	-keepSynNetsWithScore          keep nets classified as syn if the are above this score threshold
	-keepInvNetsWithScore          keep nets classified as inv if the are above this score threshold
	-minTsizes					   comma separated list of minimum target sizes
	-minQsizes					   comma separated list of minimum query sizes
\n";

# options
GetOptions ("v|verbose" => \$verbose, 
			"keepSynNetsWithScore=i" => \$keepSynNetsWithScore, 
			"keepInvNetsWithScore=i" => \$keepInvNetsWithScore,
			"minTsize=s" => \$minTsize,
			"minQsize=s" => \$minQsize) || die "$usage";	
#die "ERROR: you must be on genome to execute $0\n" if (hostname ne "genome");

if ($#ARGV < 6 ) {
	die "Not enough arguments\n$usage";
}
$inChain = $ARGV[0];
$tTB = $ARGV[1];				# filename + path
$qTB = $ARGV[2];
$tSizes = $ARGV[3];
$qSizes = $ARGV[4];
$outPrefix = $ARGV[5];
$minScore = $ARGV[6];

my $tmpDir = `set -o pipefail; mktemp -d`;
chomp($tmpDir);
print "set tmpDir to $tmpDir\n";

my $outFilteredChain = "${outPrefix}.filtered.chain";	# output filename + path
my $outFilteredNet = "${outPrefix}.filtered.net";	# output filename + path
my $outFilteredMaf = "${outPrefix}.filtered.maf";	# output filename + path

# split the minscores and sizes
my @minScores = split(/,/, $minScore);
my @minTsizes = split(/,/, $minTsize);
my @minQsizes = split(/,/, $minQsize);
die "ERROR: number of minScores differ from minTsizes\n" if ($minTsize and $#minScores != $#minTsizes);
die "ERROR: number of minScores differ from minQsizes\n" if ($minQsize and $#minScores != $#minQsizes);


my $call="";
# filter the chains
for (my $i=0; $i<=$#minScores; $i++) {
	$call="set -o pipefail; chainFilter $inChain -notQ=chrM -notT=chrM -minScore=$minScores[$i] ";
	if ($minTsize){
		$call.="-tMinSize=$minTsizes[$i] ";
	}
	if ($minQsize){
		$call.="-qMinSize=$minQsizes[$i] ";
	}
	$call.="> $tmpDir/filtered$i.chain";
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: chainFilter failed\ncall: $call\n";
}

# now cat all the files and chain
$call = "set -o pipefail; cat $tmpDir/filtered*.chain | chainSort stdin stdout | 
   chainPreNet stdin $tSizes $qSizes $outFilteredChain";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: chainPreNet command failed\ncall: $call\n";

# net and netSyntenic only
$call ="set -o pipefail; chainNet $outFilteredChain $tSizes $qSizes stdout /dev/null -minSpace=1 -rescore -linearGap=loose -tNibDir=$tTB -qNibDir=$qTB | netSyntenic stdin stdout > $tmpDir/filtered.net";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: chainNet failed. No $tmpDir/filtered.net is produced\ncall: $call\n";


# now filter the nets
# Use doScoreFilter if minSizeT and minSizeQ not set
my $filterMode = "";
if ($minTsize and $minQsize) {
	$filterMode = "-minScore $minScore -minSizeT $minTsize -minSizeQ $minQsize"
} else {
	$filterMode = "-doScoreFilter -minScore1 $minScore"
}

# only filter for syn/inv nets if the parameters are given above
my $syninvFilter = "";
if ($keepSynNetsWithScore < INT_MAX || $keepInvNetsWithScore < INT_MAX) {
	$syninvFilter = "-keepSynNetsWithScore $keepSynNetsWithScore -keepInvNetsWithScore $keepInvNetsWithScore";
}
$call ="set -o pipefail; NetFilterNonNested.perl $tmpDir/filtered.net $filterMode $syninvFilter > $outFilteredNet";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: NetFilterNonNested.perl failed\ncall $call\n";

$call ="set -o pipefail; netToAxt $outFilteredNet $outFilteredChain $tTB $qTB $tmpDir/filtered.axt";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: netToAxt failed\ncall: $call\n";

$call ="set -o pipefail; axtToMaf $tmpDir/filtered.axt $tSizes $qSizes $outFilteredMaf";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: axtToMaf failed\ncall: $call\n";

# gzip but remove any existing chain.gz and net.gz previous output files before
$call ="set -o pipefail; gzip -f $outFilteredChain $outFilteredNet $outFilteredMaf";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: call $call\n";

# cleanup
`set -o pipefail; rm -rf $tmpDir`;

print "DONE: result files: $outFilteredChain and $outFilteredNet and $outFilteredMaf\n";
#die "ERROR: $tmpDir\n";
