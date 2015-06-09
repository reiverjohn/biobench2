# biobench2
Bioinformatics benchmarking package, based on the original BioBench 
developed by Albayraktaroglu et al, 2005 "BioBench: A benchmark suite of bioinformatics applications" 

Scripts developed and packages assembled by John Johnston
Michigan State University

Last Modified: 2011

INSTALLATION

After cloning the repo, perform the following steps:

1.  cd <biobench2 sourcedir>
2.  cd HMMER/input
3.  tar xvf inputData.tgz
4.  cd ../../MUMmer/input
5.  tar xvf inputData.tgz
6.  cd ../../QuEST/input
7.  Follow directions in README.1st
8.  cd ../../BLAST/input
9.  Follow directions in README.1st
10. cd ../../velvet/input
11. Follow directions in README.1st
12  cd ../../
13. sh ./buildBench.sh
14. per ./runBench.pl

RESULTS

Your primary results file will be located in the root installation
directory, and is entitled "bioBenchResults.csv"

Program specific output can be found in the "output" subdirectory
of each application directory.

Results are presented in the format:

host,test,replicate,runTime,cpuVendor,cpuFamily,cpuModel,procSpeed,cacheSize,numCores,totalMem,compiler

These are grouped by replicate number and can be imported into
a spreadsheet program, or parsed using a script language.
