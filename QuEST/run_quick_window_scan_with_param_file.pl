#!/usr/bin/perl

## This program is distributed under the free software
## licensing model. You have a right to use, modify and
## and redistribute this software in any way you like. 


## This program is implemented by Anton Valouev, Ph.D.
## as a part of QuEST analytical pipeline. If you use
## this software as a part of another software or for publication
## purposes, you must cite the publication anouncing
## QuEST release:

## Valouev A, Johnson DS, Sundquist A, Medina C, Anton E, Batzoglou S,
## Myers RM, Sidow A (2008). Genome-wide analysis of transcription 
## factor binding sites based on ChIP-Seq data.
## Nat Methods 5, 829-834 (2008)


## This program is a wrapper that runs quick_window_scan
## on the entire genome one contig at a time

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    run_quick_window_scan_with_param_file.pl
	
    This program is a wrapper that runs the quick_window_scan
    once contig at a time
    
    -----------------------
    mandatory parameters:
    -p <params_file>               a file containing parameters for quick_window_scan

    -----------------------
    optional parameters:
    -e <exec_path>                 path to the quick window scan executable
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

$exec_path = $exec_path."/quick_window_scan";
## /setting optional arguments

if(-e $exec_path){
    ## executable exists do nothing
}
else{
    print "Failed to locate executable $exec_path.\n";
    print "Aborting quick_window_scan batch script.\n";
    exit(0);
}

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

my $bin_align_ChIP_path = "";
my $bin_align_RX_noIP_path = "";

my $genome_table_fname = "";
my $output_file = "";
my $calc_window = "";
my $ChIP_tags_threshold = "";
my $ChIP_to_background_tag_ratio = "";
my $enrichment_ratio = "";

for(my $i=0; $i<scalar(@params); $i++){
    my $cur_param = $params[$i];
    chomp($cur_param);
    my @cur_par_fields = split(/ /, $cur_param);
    my $cur_par_name = $cur_par_fields[0];
    if( $cur_par_name eq "bin_align_ChIP_path" ){
	$bin_align_ChIP_path = $cur_par_fields[2];	
    }
    elsif( $cur_par_name eq "bin_align_RX_noIP_path" ){
	$bin_align_RX_noIP_path = $cur_par_fields[2];	
    }
    elsif( $cur_par_name eq "genome_table" ){
	$genome_table_fname = $cur_par_fields[2];
    }
    elsif($cur_par_name eq "output_file"){
	$output_file  = $cur_par_fields[2];
    }
    elsif($cur_par_name eq "calc_window"){
	$calc_window = $cur_par_fields[2];
    }
    elsif($cur_par_name eq "ChIP_tags_threshold"){
	$ChIP_tags_threshold = $cur_par_fields[2];
    }
    elsif($cur_par_name eq "ChIP_to_background_tag_ratio"){
	$ChIP_to_background_tag_ratio = $cur_par_fields[2];
    }
    else{
	if($params[$i] ne ""){
	    print "Warning: unrecognized parameter: $cur_par_name\n";
	    exit(1);
	}
    }
}

print "Read the following parameters: \n\n";
print "bin_align_ChIP_path:   $bin_align_ChIP_path\n";
print "bin_align_RX_noIP_path:$bin_align_RX_noIP_path\n";

print "genome_table_file:            $genome_table_fname\n";
print "output_file:                  $output_file\n";
print "calc_window:                  $calc_window\n";
print "ChIP_tags_threshold:          $ChIP_tags_threshold\n";
print "ChIP_to_background_tag_ratio: $ChIP_to_background_tag_ratio\n";
print "\n";

my $optional_param_string = "";
if($calc_window ne ""){
    $optional_param_string = "calc_window=$calc_window";
}
if($ChIP_tags_threshold ne ""){
    $optional_param_string = $optional_param_string . " ChIP_tags_threshold=$ChIP_tags_threshold ChIP_to_background_tag_ratio=$ChIP_to_background_tag_ratio";
}

if( ! -e $genome_table_fname){
    print "Error in Quick Window Scan batch script:\n";
    print "Failed to locate genome table $genome_table_fname.\n";
    print "Aborting.\n";
    exit(0);
}

open genome_table_file, "< $genome_table_fname" || die "$genome_table_fname: $\n";
my @genome_table = <genome_table_file>;
close genome_table_file;

my @contig_names;
my @contig_sizes;
my $contig_counter = 0;

for(my $i=0; $i<scalar(@genome_table); $i++){
    my $cur_entry = $genome_table[$i];
    chomp($cur_entry);

    my @cur_entry_fields = split(/ /, $cur_entry);
    if( scalar(@cur_entry_fields) == 2 ){
	$contig_names[$contig_counter] = $cur_entry_fields[0];
	$contig_sizes[$contig_counter] = $cur_entry_fields[1];
	$contig_counter++;
    }
}

for(my $i=0; $i<scalar(@contig_names); $i++){
    my $cur_contig_name = $contig_names[$i];
    my $cur_contig_size = $contig_sizes[$i];
    my $cur_bin_align_ChIP_fname = $bin_align_ChIP_path . "/" . $cur_contig_name . ".align.bin";
    my $cur_bin_align_RX_noIP_fname;
    if($bin_align_RX_noIP_path eq "NA"){
	$cur_bin_align_RX_noIP_fname  = "NA";
    }
    else{
	$cur_bin_align_RX_noIP_fname =  $bin_align_RX_noIP_path . "/" . $cur_contig_name . ".align.bin";
    }

    my $system_command = "$exec_path contig_id=$cur_contig_name contig_size=$cur_contig_size ChIP_bin_align_file=$cur_bin_align_ChIP_fname RX_noIP_bin_align_file=$cur_bin_align_RX_noIP_fname output_file=$output_file.$cur_contig_name $optional_param_string";
    print "\n======================================\n\n";
    print "system_command: $system_command\n\n";
    
    my $error_code = system("$system_command");
    if( $error_code != 0 ){
	print "run_quick_window_scan_with_param_file.pl: Caught an error message (code $error_code). Passing the error message to the top script.\n";
	exit(3);
    }
}

my $output_file_list = "";
for(my $i=0; $i<scalar(@contig_names); $i++){
    my $cur_contig = $contig_names[$i];
    if($i>0){
	$output_file_list = $output_file_list . " ";
    }
    $output_file_list = $output_file_list . "$output_file.$cur_contig"
}

unlink("$output_file");
unlink("$output_file.sorted");
system("cat $output_file_list > $output_file");
system("sort -nr -k3 $output_file > $output_file.sorted");
#system("sort -nr +2 $output_file > $output_file.sorted");
system("rm $output_file_list");

## /top-level script
