#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

use Cwd qw(realpath);
my $fullpath = realpath($0);
my @fullpath_fields = split(/\//,$fullpath);
my $exec_path = "";

for(my $i=0; $i<scalar(@fullpath_fields)-1; $i++){
    if($fullpath_fields[$i] ne ""){
        $exec_path = $exec_path."/".$fullpath_fields[$i];
    }
}


my $makefile_fname = $exec_path . "/Makefile";
my $misc_path = $exec_path . "/misc/";
my $misc_makefile_fname = $misc_path . "Makefile";
my $config_out_fname = $exec_path . "/config.mk";

my $uname_out = `uname -a`;

my @uname_out_fields = split(/ /, $uname_out);

my $OS = $uname_out_fields[0];
my $Darwin_ver = "";
if(scalar(@uname_out_fields) >= 3){
    my @Darwin_ver_fields = split(/\./, $uname_out_fields[2]);
    if(scalar(@Darwin_ver_fields) >= 3){
	$Darwin_ver = $Darwin_ver_fields[0];
    }
}

print "Testing OS. OS = $OS.\n";
print "Checking for g++ compiler."; $|++;
my $gpp_out = `which g++`;
chomp($gpp_out);

if( $gpp_out eq ""){
    print "Error: g++ compiler not found.\n";
    exit(1);
}
else{
    print " Found.\n";
}


print "Testing -pthread flag. "; $|++;


my $pthread_try = `g++ -pthread 2> out`;
my $pthread_out = `cat out | grep unrecognized | wc -l`;

chomp($pthread_out);

if($pthread_out == 0){
    print "-pthread tested ok. Including the flag.\n";
}
else{
    print "-pthread not recognized, omitting the flag\n";
}


if(-e $config_out_fname){
    unlink($config_out_fname);
}
open config_out_file, "> $config_out_fname" || die "$config_out_fname: $!\n";

my $CCFLAGS = "-Wall -ansi -pedantic -D DEBUG -O4 -funroll-loops"; # -pthread";



if( $OS eq "Darwin" ){
    if($Darwin_ver >= 9){
	$CCFLAGS = $CCFLAGS . " -m64";
	if($pthread_out == 0){
	    $CCFLAGS = $CCFLAGS . " -pthread";
	}
    }    
}
else{
    if($pthread_out == 0){
	$CCFLAGS = $CCFLAGS . " -pthread";
    }
}

print config_out_file "CCFLAGS = $CCFLAGS\n";

close config_out_file;
