# biobench2 README File #

Released: 2015-06-09
Last Modified: 2018-05-23

### Contents ###

+ A. OVERVIEW
+ B. INSTALLATION
+ C. RESULTS

***************************************

### OVERVIEW ###
Bioinformatics benchmarking package, based on the original BioBench 
developed by Albayraktaroglu et al, 2005 "BioBench: A benchmark suite of bioinformatics applications" 

Scripts developed and packages assembled by John Johnston
Michigan State University

Last Modified: 2018-05-23

***************************************

### INSTALLATION ###

After cloning the repo, perform the following steps:

1.  cd <biobench2 sourcedir>
2.  cd ../../QuEST/input
3.  Follow directions in README.1st
4.  cd ../../BLAST/input
5.  Follow directions in README.1st
12. cd ../../
13. sh ./buildBench.sh
14. perl ./runBench.pl

***************************************

### RESULTS ###

Your primary results file will be located in the root installation
directory, and is entitled "bioBenchResults.csv"

Program specific output can be found in the "output" subdirectory
of each application directory.

Results are presented in the format:

host,test,replicate,runTime,cpuVendor,cpuFamily,cpuModel,procSpeed,cacheSize,numCores,totalMem,compiler

These are grouped by replicate number and can be imported into
a spreadsheet program, or parsed using a script language.

