#!/usr/bin/perl
######################################################################
#	Name:		benchRun.pl
#	By:		John Johnston (johnj@msu.edu)
#	Created:	September 20, 2011
#	Last Mod:	09-21-2011
#	Purpose:	Automated execution of bioinformatics
#			benchmark tests
#	Output:		Benchmarking results in bioBenchResults.csv
#			Misc. runtime output in benchResults.txt
#			Application-specific output in output subdir
#			of each individual application.
#
#	Run Time:	Total number of CPU seconds spent in kernel
#			and user modes (user and sys combined).
#
######################################################################
#
#			GENERAL INSTRUCTIONS
#
#	First run the auto-build script "buildBench.sh".  Assuming
#	there are no errors, run this script using:
#	perl ./runBench.pl, or set the executable bit on the file and
#	run as:  ./runBench.pl
#
#	Informational output annouces each test, resulting run times,
#	and replicate.
#
#	CSV output produced in bioBenchResults.csv.
#
#	MODIFY "$replicates" variable below to control the number of
#	test iterations (default = 3).
#
#####################################################################
##
## SET TEST ITERATION VALUE BELOW!!
##
$replicates = 3;
##
## SET OMP_NUM_THREADS TO CONTROL CPU USAGE FOR VELVET!!
##
$ENV{'OMP_NUM_THREADS'}=1;
##
$runStart = `date`;
open (MYFILE, '>>bioBenchResults.csv');
$hostName = `hostname`;
@memoryInfo = `free -m`;
@totalMemory = split /\s+/, $memoryInfo[1];
@sysInfo = `cat /proc/cpuinfo`;
@vendorID = split /\s+/, $sysInfo[1];
@cpuFamily = split /\s+/, $sysInfo[2];
@model = split /\s+/, $sysInfo[3];
@modelName = split /\s+/, $sysInfo[4];
$compilerInfo = `gcc -v 2>&1`;
chomp($compilerInfo);
@compilerParse = split("version",$compilerInfo);
@compiler = split /\s+/, $compilerParse[1];
$modelLength = scalar(@modelName) - 1;
@sliceModelName = @modelName[3..$modelLength];
@speed = split /\s+/, $sysInfo[6];
@cache = split /\s+/, $sysInfo[7];
@cores = split /\s+/, $sysInfo[11];
chomp($hostName);
$hwInfo =  "$vendorID[2],$cpuFamily[3],$model[2],@sliceModelName,$speed[3]MHz,$cache[3]KB,$cores[3],$totalMemory[1]MB,gcc-$compiler[1]";
print "\nBeginning Bio-Benchmarking:\n";
print MYFILE "Run Start: $runStart\n\n";
print MYFILE "host,test,replicate,runTime,cpuVendor,cpuFamily,cpuModel,procSpeed,cacheSize,numCores,totalMem,compiler\n";

for ($i=1; $i <= $replicates; $i++) {

##
##	BEDTools v. 2.13.0 
##
print "\nRunning BEDTools, replicate # $i...\n\n";
system("sleep 5"); 
@timeResults = `(time -p ./BEDTools/bin/windowBed -l 1000 -r 1000 -a ./BEDTools/input/reads.gff -b ./BEDTools/input/start.gff | ./BEDTools/bin/fastaFromBed -fi ./BEDTools/input/yeast.fasta -bed stdin -fo stdout 1> ./BEDTools/output/results.txt 2>> benchResults.txt) 2>&1`;
chomp(@timeResults);
print "@timeResults\n";
@userTime = split /\s+/, $timeResults[1];
@sysTime = split /\s+/, $timeResults[2];
$totalTime = $userTime[1] + $sysTime[1];
print MYFILE "$hostName,bedtools,$i,$totalTime,$hwInfo\n";
##
##	CLUSTALW2 v. 2.1
##
print "\nRunning Clustalw2, replicate # $i...\n\n";
system("sleep 5");
@timeResults = `(time -p ./clustalw/src/clustalw2 -INFILE=./clustalw/input/tufa420.seq -OUTPUT=PHYLIP 1> ./clustalw/output/results.txt 2>> benchResults.txt) 2>&1`;
chomp(@timeResults);
print "@timeResults\n";
@userTime = split /\s+/, $timeResults[1];
@sysTime = split /\s+/, $timeResults[2];
$totalTime = $userTime[1] + $sysTime[1];
print MYFILE "$hostName,clustalw2,$i,$totalTime,$hwInfo\n";
##
##	BLAST v. 2.2.25
##
print "\nRunning BLAST, replicate # $i...\n\n";
system("sleep 5");
@timeResults = `(time -p ./BLAST/bin/blastall -p blastn -d ./BLAST/input/nt -i ./BLAST/input/batch2.fa  1> ./BLAST/output/results.txt 2>> benchResults.txt) 2>&1`;
chomp(@timeResults);
print "@timeResults\n";
@userTime = split /\s+/, $timeResults[1];
@sysTime = split /\s+/, $timeResults[2];
$totalTime = $userTime[1] + $sysTime[1];
print MYFILE "$hostName,blast,$i,$totalTime,$hwInfo\n";
##
##	HMMER v. 3.0
##
print "\nRunning HMMER, replicate # $i...\n\n";
system("sleep 5");
@timeResults = `(time -p ./HMMER/bin/hmmsearch ./HMMER/input/tufa420.hmm ./HMMER/input/uniprot_sprot.fasta  1> ./HMMER/output/results.txt 2>> benchResults.txt) 2>&1`;
chomp(@timeResults);
print "@timeResults\n";
@userTime = split /\s+/, $timeResults[1];
@sysTime = split /\s+/, $timeResults[2];
$totalTime = $userTime[1] + $sysTime[1];
print MYFILE "$hostName,hmmer,$i,$totalTime,$hwInfo\n";
##
##	MUMmer v. 3.0
##
print "\nRunning MUMmer, replicate # $i...\n\n";
system("sleep 5");
@timeResults = `(time -p ./MUMmer/mummer -b -c  ./MUMmer/input/hs_chrY.fa ./MUMmer/input/hs_chr17.fa  1> ./QuEST/output/results.txt 2>> benchResults.txt) 2>&1`;
chomp(@timeResults);
print "@timeResults\n";
@userTime = split /\s+/, $timeResults[1];
@sysTime = split /\s+/, $timeResults[2];
$totalTime = $userTime[1] + $sysTime[1];
print MYFILE "$hostName,mummer,$i,$totalTime,$hwInfo\n";
##
##	QuEST v. 2.4
##
print "\nRunning QuEST, replicate # $i...\n\n";
system("sleep 6");
@timeResults = `(time -p ./QuEST/generate_QuEST_parameters.pl -solexa_align_ChIP ./QuEST/input/GABP.align_25.hg18 -solexa_align_RX_noIP ./QuEST/input/Jurkat_RX_noIP.align_25.hg18 -gt ./QuEST/input/genome_table -ap ./QuEST/input/QuEST_analysis -ChIP_name GABP_Jurkat < ./QuEST/input/progCMD 1> ./QuEST/output/results.txt 2>> benchResults.txt) 2>&1`;
chomp(@timeResults);
print "@timeResults\n";
@userTime = split /\s+/, $timeResults[1];
@sysTime = split /\s+/, $timeResults[2];
$totalTime = $userTime[1] + $sysTime[1];
print MYFILE "$hostName,quest,$i,$totalTime,$hwInfo\n";
##
##	VELVET v. 1.1.05
##
print "\nRunning Velvet, replicate # $i...\n\n";
system("sleep 5");
@timeResults = `(time -p ((./velvet/velveth ./velvet/input/ASSEM 25 -shortPaired -fastq ./velvet/input/dmelRNA-reads.fastq 1> ./velvet/output/results.txt 2>> benchResults.txt) && (./velvet/velvetg ./velvet/input/ASSEM -exp_cov auto -cov_cutoff auto 1> ./velvet/output/results.txt 2>> benchResults.txt ))) 2>&1`;
chomp(@timeResults);
print "@timeResults\n";
@userTime = split /\s+/, $timeResults[1];
@sysTime = split /\s+/, $timeResults[2];
$totalTime = $userTime[1] + $sysTime[1];
print MYFILE "$hostName,velvet,$i,$totalTime,$hwInfo\n";

}  ## end replicate loop

print "DONE!\n";
print "Benchmarking results written to: bioBenchResults.csv\n";
print "Extraneous run output dumped to: benchResults.txt\n";
$runEnd = `date`;
print MYFILE "Run End: $runEnd\n";
close(MYFILE);
