#!/bin/bash
######################################################################
#	Name:		buildBench.sh
#	By:		John Johnston (johnj@msu.edu)
#	Created:	September 20, 2011
#	Last Mod:	05/23/2018
#	Purpose:	Build applications for Biobench2 and prep
#			data sets.
######################################################################
basedir=`pwd`
echo
echo "*** Beginning Bio-Benchmarking Build ***"
echo
sleep 5 
echo "Building clustalw2..."
cd clustalw
./configure 
make
cd ..
echo
echo "Building QuEST..."
sleep 3
cd QuEST
./configure.pl
make clean
make
cd ..
echo
echo "Building BEDTools..."
sleep 3
cd BEDTools
make clean
make all
cd ..
echo
echo "BLAST is pre-built, skipping..."
sleep 3
echo
echo "Building PHYLIP..."
sleep 3
cd PHYLIP/src
make clean
make install
cd $basedir
echo
echo "Building HMMER..."
sleep 3
cd HMMER
INSTALLDIR=$basedir/HMMER
./configure --prefix=$INSTALLDIR
make clean
make
make check
make install
cd $basedir/HMMER/input
tar xvf inputData.tgz
cd $basedir
echo
echo "Building MUMmer..."
sleep 3
cd MUMmer
make clean
make check
make install
cd $basedir/MUMmer/input
tar xvf inputData.tgz
cd $basedir
echo
echo "Building Velvet..."
sleep 3
cd velvet
make clean
make
cd $basedir/velvet/input
cat dmelRNA-reads-* > dmelRNA-reads.tgz
tar xvf dmelRNA-reads.tgz
cd $basedir
echo
echo "DONE!"
echo
echo "NOTE: You will need to download the NCBI "nr" data set to BLAST/input for BLAST benchmarking!"
echo
echo "Now run ./runBench.pl to perform benchmarking"
