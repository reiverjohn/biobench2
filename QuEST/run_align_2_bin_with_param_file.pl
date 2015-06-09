#!/usr/bin/perl

## This program is distributed under the free software
## licensing model. You have a right to use, modify and
## and redistribute this software in any way you like. 


## This program is implemented by Anton Valouev, Ph.D.
## as a part of QuEST analytical pipeline. If you use
## this software as a part of another software or for publication
## purposes, you must cite the publication anouncing
## QuEST release:

## the citation goes here

## This program is a wrapper that runs generate_peak_profile program
## on the entire genome one contig at a time

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    run_align_with_param_file.pl
	
    This program is a wrapper that runs the generate_peak_profile
    
    -----------------------
    mandatory parameters:
    -p <params_file>               a file containing parameters calibrate_peak_shift

    -----------------------
    optional parameters:
    -e <exec_path>                 path to the calibrate_peak_shift executable
    -h                             to display this help

};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory argmuments
my $params_fname = "";
## /mandatory arguments

## setting optional arguments
use Cwd qw(realpath);
my $fullpath = realpath($0);
my @fullpath_fields = split(/\//,$fullpath);
my $exec_path = "";

for(my $i=0; $i<scalar(@fullpath_fields)-1; $i++){
    if($fullpath_fields[$i] ne ""){
	$exec_path = $exec_path."/".$fullpath_fields[$i];
    }
}

$exec_path = $exec_path . "/align_2_bin";

if( ! -e $exec_path ){
    print "Error in align_2_bin batch script:\n";
    print "Failed to locate executable $exec_path.\n";
    print "Aborting.\n";
    exit(0);
}

## /setting optional arguments

## reading arguments
while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;    
    if ( $this_arg eq '-h') {die "$usage\n";}
    elsif ( $this_arg eq '-p') {$params_fname = shift @ARGV;}
    elsif ( $this_arg eq '-e') {$exec_path = shift @ARGV;}
    else{ print "Warning: unknown parameter: $this_arg\n"; }
}

my $die_pars_bad = "false";       ## die if parameters are bad
my $bad_par = "";                 ## store bad parameter here

print "\n=====================\n\n";
print "mandatory parameters: \n\n";
print "params_file:      $params_fname\n";
print "exec_path:        $exec_path\n";
if( $params_fname eq ''){ $die_pars_bad = "true"; $bad_par." "."$params_fname"; }
if( $exec_path eq ''){ $die_pars_bad = "true"; $bad_par." "."$exec_path"; }
print "\n=====================\n\n";

if( $die_pars_bad eq "true"){
    print "$usage";
    print "Bad parameters: $bad_par\n";
    exit(0);
}
## /reading arguments

## top-level script
open params_file, "< $params_fname" || die "$params_fname: $\n";
my @params = <params_file>;
close params_file;

my $QuEST_align_file =     "** missing **";
my $genome_table_fname =   "** missing **";
my $output_path =    "** missing **";
my $collapse_reads_flag = "false";

for(my $i=0; $i<scalar(@params); $i++){
    my $cur_param = $params[$i];
    chomp($cur_param);
    
    my @cur_par_fields = split(/ /, $cur_param);
    if(scalar(@cur_par_fields >= 2)){
	my $cur_par_name = $cur_par_fields[0];
	if ($cur_par_name eq "align_file"){
	    $QuEST_align_file = $cur_par_fields[2];	
	}
	elsif( $cur_par_name eq "genome_table" ){
	    $genome_table_fname= $cur_par_fields[2];
	}
	elsif($cur_par_name eq "output_path"){
	    $output_path = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "collapse_reads"){
	    $collapse_reads_flag = $cur_par_fields[2];
	}
	else{
	    if($params[$i] ne ""){
		print "Warning: unrecognized parameter: $cur_par_name";
	    }
	}
    }
}

print "Read the following parameters: \n\n";

print "align_file:                        $QuEST_align_file\n";
print "genome_table:                      $genome_table_fname\n";
print "output_path:                       $output_path\n";
print "\n";

my $optional_param_string = "";

if($collapse_reads_flag eq "true"){
    $optional_param_string = $optional_param_string . " collapse_reads=true";
}

my $system_command = $exec_path . " align_file=$QuEST_align_file genome_table=$genome_table_fname output_path=$output_path" . $optional_param_string;
print "system_command: $system_command\n";
my $error_code = system("$system_command\n");
if($error_code != 0){
    print "run_align_2_bin_with_param_file.pl: Caught an error message (code $error_code), passing error message to the top script.\n";
    exit(3);
}


