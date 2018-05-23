#!/bin/bash
#
## Purpose:  Build applications for Biobench2
## By:  John Johnston
## Last Mod: 2017-07-14
#
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
cd ../../
echo
echo "Building HMMER..."
sleep 3
cd HMMER
INSTALLDIR=.
./configure --prefix=$INSTALLDIR
make clean
make
make check
make install
cd ..
echo
echo "Building MUMmer..."
sleep3
cd MUMmer
make clean
make check
make install
cd ..
echo
echo "Building Velvet..."
sleep 3
cd velvet
make clean
make
cd input
cat dmelRNA-reads-* > dmelRNA-reads.tgz
tar xvf dmelRNA-reads.tgz
echo
echo "DONE!"
echo
echo
echo "Now Run ./runBench.pl to perform benchmarking"
