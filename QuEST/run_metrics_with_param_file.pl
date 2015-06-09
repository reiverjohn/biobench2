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
    run_metrics_with_param_file.pl
	
    This program is a wrapper that runs metrics
    
    -----------------------
    mandatory parameters:
    -p <params_file>               a file containing parameters for this script

    -----------------------
    optional parameters:
    -e <exec_path>                 path to the metrics executable
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

$exec_path = $exec_path."/metrics";

if( ! -e $exec_path ){
    print "Error in generate profile batch script:\n";
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

my $missing = "** missing **";

my $genome_table_fname=     "** missing **";

my $ChIP_bin_align_path =   "** missing **";
my $background_bin_align_path = "** missing **";

my $QuEST_calls_fname     =     "** missing **";
my $kde_bandwidth     =     "** missing **";
my $region_width      =     "** missing **";
my $output_fname      =     "** missing **";
my $dist_threshold    =     "** missing **";

my $ChIP_tags = "** missing **";
my $background_tags = "** missing **";

my $ChIP_basal_level = "** missing **";
my $background_basal_level = "** missing **";

my $ChIP_basal_tag_rate = "** missing **";
my $background_basal_tag_rate = "** missing **";

my $report_fname = $missing;

for(my $i=0; $i<scalar(@params); $i++){
    my $cur_param = $params[$i];
    chomp($cur_param);
    
    my @cur_par_fields = split(/ /, $cur_param);
    if(scalar(@cur_par_fields >= 2)){
	my $cur_par_name = $cur_par_fields[0];
	if ($cur_par_name eq "genome_table"){
	    $genome_table_fname = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "ChIP_tags"){
	    $ChIP_tags = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "background_tags"){
	    $background_tags = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "ChIP_basal_level"){
	    $ChIP_basal_level = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "background_basal_level"){
	    $background_basal_level = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "ChIP_basal_tag_rate"){
	    $ChIP_basal_tag_rate = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "background_basal_tag_rate"){
	    $background_basal_tag_rate = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "ChIP_bin_align_path"){
	    $ChIP_bin_align_path = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "background_bin_align_path"){
	    $background_bin_align_path = $cur_par_fields[2];	
	}
	elsif( $cur_par_name eq "QuEST_calls_file" ){
	    $QuEST_calls_fname= $cur_par_fields[2];
	}
	elsif($cur_par_name eq "kde_bandwidth"){
	    $kde_bandwidth = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "region_width"){
	    $region_width = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "dist_threshold"){
	    $dist_threshold = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "output_file"){
	    $output_fname = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "report_file"){
	    $report_fname = $cur_par_fields[2];
	}
	else{
	    if($params[$i] ne ""){
		print "Warning: unrecognized parameter: $cur_par_name";
	    }
	}
    }
}

print "exec:                              $exec_path\n\n";

print "Read the following parameters: \n\n";

print "genome_table:                      $genome_table_fname\n";
print "ChIP_bin_align_path:               $ChIP_bin_align_path\n";
print "background_bin_align_path:         $background_bin_align_path\n";

print "ChIP_tags:                         $ChIP_tags\n";
print "background_tags:                   $background_tags\n";

print "QuEST_calls_file:                  $QuEST_calls_fname\n";
print "output_file:                       $output_fname\n";
print "report_file:                       $report_fname\n";

print "kde_bandwidth:                     $kde_bandwidth\n";     
print "region_width:                      $region_width\n";     
print "dist_threshold:                    $dist_threshold\n";     


print "\n";

if(-e $report_fname){
    unlink($report_fname);
}
if(-e $output_fname){
    unlink($output_fname);
}



## reading genome table

my $system_command = $exec_path . " genome_table=$genome_table_fname ChIP_bin_align_path=$ChIP_bin_align_path background_bin_align_path=$background_bin_align_path QuEST_calls_file=$QuEST_calls_fname output_file=$output_fname region_width=$region_width kde_bandwidth=$kde_bandwidth dist_threshold=$dist_threshold ChIP_basal_level=$ChIP_basal_level background_basal_level=$background_basal_level ChIP_basal_tag_rate=$ChIP_basal_tag_rate background_basal_tag_rate=$background_basal_tag_rate ChIP_tags=$ChIP_tags background_tags=$background_tags report_file=$report_fname";
    print "system command: $system_command\n";
    my $error_code = system("$system_command\n");
    if($error_code != 0){
	print "run_strand_correlation_with_param_file.pl: Caught an error message (code $error_code), passing error message to the top script.\n";
	exit(3);
    }
#    #print "error_code: $error_code\n";
#    #exit(0);#
#
#}


