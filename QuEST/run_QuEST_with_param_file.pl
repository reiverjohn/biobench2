#!/usr/bin/perl

## This program is distributed under the free software
## licensing model. You have a right to use, modify and
## and redistribute this software in any way you like. 
## However, any redistribution of this software should
## indicate that this code was derived from the QuEST
## analytical software and a citation:

## Valouev A, Johnson DS, Sundquist A, Medina C, Anton E, Batzoglou S,
## Myers RM, Sidow A (2008). Genome-wide analysis of transcription
## factor binding sites based on ChIP-Seq data.
## Nat Methods 5, 829-834 (2008)

## Any publication that uses results from this code's
## deriviative should cite this paper.

## This program is implemented by Anton Valouev, Ph.D.
## as a part of QuEST analytical pipeline. If you use
## this software as a part of another software or for publication
## purposes, you must cite the publication anouncing
## QuEST release:

## Valouev A, Johnson DS, Sundquist A, Medina C, Anton E, Batzoglou S,
## Myers RM, Sidow A (2008). Genome-wide analysis of transcription
## factor binding sites based on ChIP-Seq data.
## Nat Methods 5, 829-834 (2008)

## This program is a master scritp for the QuEST analysis pipeline

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    run_QuEST_with_param_file.pl
	
    This program is a wrapper that runs QuEST analytical pipeline
    
    -----------------------
    mandatory parameters:
    
    -ap <analysis_path>            QuEST analysis path directory
    
    or

    -p <params_file>               a file containing QuEST batch parameters    

    -----------------------
    optional parameters:

    -h                             for help

};

my $help_msg = qq{
    run_QuEST_with_param_file.pl
	
    This program is a wrapper that runs QuEST analytical pipeline
    
    -----------------------
    mandatory parameters:

    -ap <analysis_path>            QuEST analysis path directory
    
    or

    -p <params_file>               a file containing QuEST batch parameters    

    -----------------------
    optional parameters:

    -h

    (Use directives below to override execution directives in QuEST batch parameter file)

    -quick_window_scan                 to run quick window scan
    -calibrate_peak_shift              to estimate peak shift

    -generate_profile_ChIP             to generate ChIP profile
    -generate_profile_background       to generate background profile
    -generate_profile_pseudo_ChIP      to generate RX_noIP profile    

    -report_bad_control_regions        to report regions with enrichment in background data
    
    -peak_call_ChIP                    to call peaks in ChIP data
    -peak_call_pseudo_CHIP             to call peaks in pseudo_ChIP_data 

    -metrics                           to generate metrics for the peaks and regions
    -filter_ChIP_calls                 to filter ChIP_regions based on the metrics       

    -bed_tracks_ChIP                   to generate bed files for the ChIP data
    -bed_tracks_RX_noIP                to generate bed files for the RX_noIP data

    -wig_tracks_ChIP                   to generate wig tracks for the ChIP data
    -wig_tracks_background             to generate wig tracks for the backgrond data

    -bed_graph_tracks_ChIP             to generate bed graph tracks for the ChIP data

    -calls_tracks                      to generate bed tracks for the QuEST calls
    
    -h                                 to display this help

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

my $QuEST_info_fname = $exec_path . "/QuEST.info";
open QuEST_info_file, "< $QuEST_info_fname" || die "$QuEST_info_fname: $!\n";
while(<QuEST_info_file>){
    chomp;
    print "$_\n";
}
close QuEST_info_file;

my $override_QuEST_directives = "false";

#my $override_align_2_bin_ChIP = "false";
#my $override_align_2_bin_background = "false";
#my $override_align_2_bin_RX_noIP = "false";
#my $override_align_2_bin_pseudo_ChIP = "false";

my $override_quick_window_scan = "false";
my $override_calibrate_peak_shift = "false";
my $override_generate_profile_ChIP = "false";
my $override_generate_profile_background = "false";
my $override_generate_profile_pseudo_ChIP = "false";

my $override_report_bad_control_regions = "false";

my $override_peak_call_ChIP = "false";
my $override_peak_call_pseudo_ChIP = "false";

#my $override_q_values = "false";
my $override_metrics_ChIP = "false";
my $override_filter_ChIP_calls = "false";

my $override_wig_tracks_ChIP = "false";
my $override_wig_tracks_RX_noIP = "false";
my $override_wig_tracks_background = "false";
my $override_calls_tracks = "false";
my $override_bed_tracks_ChIP = "false";
my $override_bed_tracks_RX_noIP = "false";

my $override_bed_graph_tracks_ChIP = "false";


## /setting optional arguments

## reading arguments
my $analysis_path = "";
while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;    
    if ( $this_arg eq '-h') {print "$help_msg\n"; exit(0)}
    elsif ( $this_arg eq '-p') {$params_fname = shift @ARGV;}
    elsif ( $this_arg eq '-ap'){
	$analysis_path = shift @ARGV ;
	
	if(length($analysis_path) == 0){
	    print "bad analysis path $analysis_path\n";
	    exit(4);
	}

	while(substr($analysis_path,length($analysis_path)-1,1) eq "/"){
	    $analysis_path  = substr($analysis_path, 0, length($analysis_path)-1);
	}
	
	#if(substr($analysis_path,0,1) eq "/"){
	    # absolute path
	$params_fname = $analysis_path . "/parameters/QuEST.batch.pars";
	#}
	#else
    }
#    elsif ( $this_arg eq '-align_2_bin_ChIP') {$override_align_2_bin_ChIP = "true";}
#    elsif ( $this_arg eq '-align_2_bin_background') {$override_align_2_bin_background = "true";}
#    elsif ( $this_arg eq '-align_2_bin_pseudo_ChIP') {$override_align_2_bin_pseudo_ChIP = "true";}
#    elsif ( $this_arg eq '-align_2_bin_RX_noIP') {$override_align_2_bin_RX_noIP = "true";}

    elsif ( $this_arg eq '-quick_window_scan') {$override_quick_window_scan = "true";}
    elsif ( $this_arg eq '-calibrate_peak_shift') {$override_calibrate_peak_shift = "true";}
    elsif ( $this_arg eq '-generate_profile_ChIP') {$override_generate_profile_ChIP = "true";}
    elsif ( $this_arg eq '-generate_profile_background') {$override_generate_profile_background = "true";}
    elsif ( $this_arg eq '-generate_profile_pseudo_ChIP') {$override_generate_profile_pseudo_ChIP = "true";}
    elsif ( $this_arg eq '-report_bad_control_regions') {$override_report_bad_control_regions = "true";}
    elsif ( $this_arg eq '-peak_call_ChIP') {$override_peak_call_ChIP = "true";}
    elsif ( $this_arg eq '-peak_call_pseudo_ChIP') {$override_peak_call_pseudo_ChIP = "true";}
    elsif ( $this_arg eq '-metrics') {$override_metrics_ChIP = "true";}
    elsif ( $this_arg eq '-bed_tracks_ChIP') {$override_bed_tracks_ChIP = "true";}
    elsif ( $this_arg eq '-bed_tracks_RX_noIP') {$override_bed_tracks_RX_noIP = "true";}
    elsif ( $this_arg eq '-wig_tracks_ChIP') {$override_wig_tracks_ChIP = "true";}
    elsif ( $this_arg eq '-wig_tracks_RX_noIP') {$override_wig_tracks_RX_noIP = "true";}
    elsif ( $this_arg eq '-wig_tracks_background') {$override_wig_tracks_background = "true";}
    elsif ( $this_arg eq '-bed_graph_tracks_ChIP') {$override_bed_graph_tracks_ChIP = "true";}
    elsif ( $this_arg eq '-filter_ChIP_calls') {$override_filter_ChIP_calls = "true";}
    elsif ( $this_arg eq '-calls_tracks') {$override_calls_tracks = "true";}
#    elsif ( $this_arg eq '-q-values') { $override_q_values = "true";}
    else{ print "Warning: unknown parameter: $this_arg\n"; exit(4); }
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

if($override_quick_window_scan eq "true" ||
   $override_calibrate_peak_shift eq "true" ||
   $override_generate_profile_ChIP eq "true" ||
   $override_generate_profile_background eq "true" ||
   $override_generate_profile_pseudo_ChIP eq "true" ||
   $override_report_bad_control_regions eq "true" || 
   $override_peak_call_ChIP eq "true" ||
   $override_peak_call_pseudo_ChIP eq "true" ||
   $override_metrics_ChIP eq "true" ||
   $override_bed_tracks_ChIP eq "true" ||
   $override_bed_tracks_RX_noIP eq "true" ||
   $override_wig_tracks_ChIP eq "true" ||
   $override_wig_tracks_RX_noIP eq "true" ||   
   $override_wig_tracks_background eq "true" ||
   $override_filter_ChIP_calls eq "true" ||
   $override_calls_tracks eq "true" ||
   $override_bed_graph_tracks_ChIP eq "true"){ 
    $override_QuEST_directives = "true";
}

## top-level script
if( ! -e $params_fname){
    print "Failed to open $params_fname. Aborting.\n";
    exit(1);
}
open params_file, "< $params_fname" || die "$params_fname: $!\n";
my @params = <params_file>;
close params_file;

my $analysis_with_FDR = "";

my $quick_window_scan_param_fname = "";
my $calibrate_peak_shift_param_fname = "";

my $ChIP_generate_profile_param_fname = "";
my $background_generate_profile_param_fname = "";
my $pseudo_ChIP_generate_profile_param_fname = "";

my $report_bad_control_regions_param_fname = "";

my $ChIP_peak_caller_param_fname = "";
my $pseudo_ChIP_peak_caller_param_fname = "";

my $metrics_ChIP_param_fname = "";
my $filter_ChIP_calls_param_fname = "";

my $bed_graph_tracks_ChIP_param_fname = "";

my $ChIP_reads_bed_param_fname = "";
my $ChIP_normalized_wig_param_fname = "";
my $ChIP_unnormalized_wig_param_fname = "";

my $RX_noIP_reads_bed_param_fname = "";
my $background_normalized_wig_param_fname = "";
my $background_unnormalized_wig_param_fname = "";

my $ChIP_calls_tracks_param_fname = "";

my $QuEST_parametrization_report_fname = "";
my $metrics_report_fname = "";
my $filter_report_fname = "";

my $quick_window_scan_ChIP = "false";
my $calibrate_peak_shift_ChIP = "false";
my $generate_profile_ChIP = "false";
my $generate_profile_background = "false";
my $generate_profile_pseudo_ChIP = "false";
my $report_bad_control_regions = "true";
my $peak_caller_ChIP_vs_background = "false";
my $peak_caller_pseudo_ChIP_vs_background = "false";
my $metrics_ChIP = "false";
my $filter_ChIP_calls = "false";

my $bed_tracks_ChIP = "false";
my $bed_tracks_RX_noIP = "false";
my $wig_tracks_ChIP = "false";
my $wig_tracks_RX_noIP = "false";
my $wig_tracks_background = "false";
my $calls_tracks = "false";
my $bed_graph_tracks_ChIP = "false";

my $log_fname = "";
my $output_fname = "";

my $undef = "undef";

my $scores_dir = $undef;


for(my $i=0; $i<scalar(@params); $i++){
    my $cur_param = $params[$i];
    chomp($cur_param);
    
    if( length($cur_param) > 0 ){
	if( substr($cur_param,0,1) ne "#" ){
	    my @cur_par_fields = split(/ /, $cur_param);
	    
	    if(scalar(@cur_par_fields >= 2)){
		my $cur_par_name = $cur_par_fields[0];       
		if ($cur_par_name eq "analysis_with_FDR"){
		    $analysis_with_FDR = $cur_par_fields[2];	
		}
		elsif( $cur_par_name eq "quick_window_scan:" ){
		    if( $cur_par_fields[1] eq "ChIP" ){
			$quick_window_scan_ChIP = "true";
		    }
		}
		elsif( $cur_par_name eq "quick_window_scan:" ){
		    if( $cur_par_fields[1] eq "ChIP" ){
			$quick_window_scan_ChIP = "true";
		    }
		}
		elsif( $cur_par_name eq "calibrate_peak_shift:" ){
		    if( $cur_par_fields[1] eq "ChIP" ){
			$calibrate_peak_shift_ChIP = "true";
		    }
		}
		elsif( $cur_par_name eq "generate_profile:"){
		    if( $cur_par_fields[1] eq "ChIP" ){
			$generate_profile_ChIP = "true";
		    }
		    elsif( $cur_par_fields[1] eq "background" ){
			$generate_profile_background = "true";
		    }
		    elsif( $cur_par_fields[1] eq "pseudo_ChIP" ){
			$generate_profile_pseudo_ChIP = "true";
		    }
		}
		elsif( $cur_par_name eq "report_bad_control_regions:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$report_bad_control_regions = "true";
		    }
		}
		elsif( $cur_par_name  eq "peak_caller:" ){
		    if( scalar(@cur_par_fields) == 4 ){
			if( $cur_par_fields[1] eq "ChIP" &&
			    $cur_par_fields[2] eq "vs" &&
			    $cur_par_fields[3] eq "background"){
			    $peak_caller_ChIP_vs_background = "true";
			}
			elsif( $cur_par_fields[1] eq "pseudo_ChIP" &&
			    $cur_par_fields[2] eq "vs" &&
			    $cur_par_fields[3] eq "background"){
			    $peak_caller_pseudo_ChIP_vs_background = "true";
			}
		    }
		}
		elsif( $cur_par_name eq "metrics_ChIP:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$metrics_ChIP = "true";
		    }
		}
		elsif( $cur_par_name eq "bed_graph_tracks_ChIP:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$bed_graph_tracks_ChIP = "true";
		    }
		}
		elsif( $cur_par_name eq "bed_tracks_ChIP:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$bed_tracks_ChIP = "true";
		    }
		}
		elsif( $cur_par_name eq "bed_tracks_RX_noIP:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$bed_tracks_RX_noIP = "true";
		    }
		}
		elsif( $cur_par_name eq "wig_tracks_ChIP:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$wig_tracks_ChIP = "true";
		    }
		}
		elsif( $cur_par_name eq "wig_tracks_RX_noIP:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$wig_tracks_RX_noIP = "true";
		    }
		}
		elsif( $cur_par_name eq "wig_tracks_background:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$wig_tracks_background = "true";
		    }
		}
		elsif( $cur_par_name eq "filter_ChIP_calls:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$filter_ChIP_calls = "true";
		    }
		}
		elsif( $cur_par_name eq "calls_tracks:" ){
		    if( $cur_par_fields[1] eq "true" ){
			$calls_tracks = "true";
		    }
		}		
		
#		elsif ($cur_par_name eq "align_2_bin_ChIP_param_file"){
#		    $align_2_bin_ChIP_param_fname = $cur_par_fields[2];	
#		}
#		elsif ($cur_par_name eq "align_2_bin_background_param_file"){
#		    $align_2_bin_background_param_fname = $cur_par_fields[2];	
#		}
#		elsif ($cur_par_name eq "align_2_bin_RX_noIP_param_file"){
#		    $align_2_bin_RX_noIP_param_fname = $cur_par_fields[2];	
#		}
#		elsif ($cur_par_name eq "align_2_bin_pseudo_ChIP_param_file"){
#		    $align_2_bin_pseudo_ChIP_param_fname = $cur_par_fields[2];	
#		}

		elsif( $cur_par_name eq "scores_dir"){
		    $scores_dir = $cur_par_fields[2];
		}
		elsif ($cur_par_name eq "quick_window_scan_param_file"){
		    $quick_window_scan_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "calibrate_peak_shift_param_file"){
		    $calibrate_peak_shift_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "ChIP_generate_profile_param_file"){
		    $ChIP_generate_profile_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "background_generate_profile_param_file"){
		    $background_generate_profile_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "pseudo_ChIP_generate_profile_param_file"){
		    $pseudo_ChIP_generate_profile_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "report_bad_control_regions_param_file"){
		    $report_bad_control_regions_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "ChIP_peak_caller_param_file"){
		    $ChIP_peak_caller_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "pseudo_ChIP_peak_caller_param_file"){
		    $pseudo_ChIP_peak_caller_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "metrics_ChIP_param_file"){
		    $metrics_ChIP_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "filter_ChIP_calls_param_file"){
		    $filter_ChIP_calls_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "bed_graph_tracks_ChIP_param_file"){
		    $bed_graph_tracks_ChIP_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "ChIP_reads_bed_param_file"){
		    $ChIP_reads_bed_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "RX_noIP_reads_bed_param_file"){
		    $RX_noIP_reads_bed_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "ChIP_normalized_wig_param_file"){
		    $ChIP_normalized_wig_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "ChIP_unnormalized_wig_param_file"){
		    $ChIP_unnormalized_wig_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "background_normalized_wig_param_file"){
		    $background_normalized_wig_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "background_unnormalized_wig_param_file"){
		    $background_unnormalized_wig_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "filtered_calls_tracks_param_file"){
		    $ChIP_calls_tracks_param_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "log_file"){
		    $log_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "QuEST_parametrization_report_file"){
		    $QuEST_parametrization_report_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "metrics_report_file"){
		    $metrics_report_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "filter_report_file"){
		    $metrics_report_fname = $cur_par_fields[2];	
		}
		elsif ($cur_par_name eq "output_file"){
		    $output_fname = $cur_par_fields[2];	
		}
		else{
		    if($params[$i] ne ""){
			print "Warning: unrecognized parameter: $cur_par_name in $params_fname\n"; exit(4);
		    }
		}
	    }
	}
    }
}

if($override_QuEST_directives eq "true"){
    if($override_quick_window_scan eq "true" ){
	$quick_window_scan_ChIP = "true";
    }
    else{ $quick_window_scan_ChIP = "false"; }

    if($override_calibrate_peak_shift eq "true" ){
	$calibrate_peak_shift_ChIP = "true";
    }
    else{ $calibrate_peak_shift_ChIP = "false"; }

    if($override_generate_profile_ChIP eq "true" ){
	$generate_profile_ChIP = "true";
    }
    else{ $generate_profile_ChIP = "false"; }

    if($override_generate_profile_background eq "true" ){
	$generate_profile_background = "true";
    }
    else{ $generate_profile_background = "false"; }

    if( $override_generate_profile_pseudo_ChIP eq "true" ){
	$generate_profile_pseudo_ChIP = "true";
    }
    else{ $generate_profile_pseudo_ChIP = "false"; }

    if( $override_report_bad_control_regions eq "true" ){
	$report_bad_control_regions = "true";
    }
    else{ $report_bad_control_regions = "false"}
    

    if( $override_peak_call_ChIP eq "true" ){
	$peak_caller_ChIP_vs_background = "true";
	$metrics_ChIP = "true";
    }
    else{ $peak_caller_ChIP_vs_background = "false"; }
    
    if( $override_peak_call_pseudo_ChIP eq "true" ){
	$peak_caller_pseudo_ChIP_vs_background = "true";
	$analysis_with_FDR = "true";
    }
    else{
	$peak_caller_pseudo_ChIP_vs_background = "false";
    }

    if( $override_metrics_ChIP eq "true" ){
	$metrics_ChIP = "true";
    }
    else{
	$metrics_ChIP = "false";
    }

    if( $override_calls_tracks eq "true" ){
	$calls_tracks = "true";
    }
    else{
	$calls_tracks = "false";
    }

    if( $override_filter_ChIP_calls eq "true" ){
	$filter_ChIP_calls = "true";
    }
    else{
	$filter_ChIP_calls = "false";
    }

    if( $override_bed_graph_tracks_ChIP eq "true" ){
	$bed_graph_tracks_ChIP = "true";
    }
    else{
	$bed_graph_tracks_ChIP = "false";
    }

    if( $override_bed_tracks_ChIP eq "true" ){
	$bed_tracks_ChIP = "true";
    }
    else{
	$bed_tracks_ChIP = "false";
    }

    if( $override_bed_tracks_RX_noIP eq "true" ){
	$bed_tracks_RX_noIP = "true";
    }
    else{
	$bed_tracks_RX_noIP = "false";
    }

    if( $override_wig_tracks_ChIP eq "true" ){
	$wig_tracks_ChIP = "true";
    }
    else{
	$wig_tracks_ChIP = "false";
    }
    
    if( $override_wig_tracks_RX_noIP eq "true" ){
	$wig_tracks_RX_noIP = "true";
    }
    else{
	$wig_tracks_RX_noIP = "false";
    }
    if( $override_wig_tracks_background eq "true" ){
	$wig_tracks_background = "true";
    }
    else{
	$wig_tracks_background = "false";
    }
}

print "override_QuEST_directives:               $override_QuEST_directives\n";
print "\n";
print "\n";
print "quick_window_scan_ChIP:                  $quick_window_scan_ChIP\n";
print "calibrate_peak_shift_ChIP:               $calibrate_peak_shift_ChIP\n";
print "generate_profile_ChIP:                   $generate_profile_ChIP\n";
print "generate_profile_pseudo_ChIP:            $generate_profile_pseudo_ChIP\n";
print "generate_profile_background:             $generate_profile_background\n";
print "peak_caller_ChIP_vs_background:          $peak_caller_ChIP_vs_background\n";
print "peak_caller_pseudo_ChIP_vs_background:   $peak_caller_pseudo_ChIP_vs_background\n";
print "report_bad_control_regions:              $report_bad_control_regions\n";
print "metrics_ChIP:                            $metrics_ChIP\n";
print "filter_ChIP_calls:                       $filter_ChIP_calls\n";
print "bed_tracks_ChIP:                         $bed_tracks_ChIP\n";
print "bed_tracks_RX_noIP:                      $bed_tracks_RX_noIP\n";
print "bed_graph_tracks_ChIP:                   $bed_graph_tracks_ChIP\n";
print "wig_tracks_ChIP:                         $wig_tracks_ChIP\n";
print "wig_tracks_RX_noIP:                      $wig_tracks_RX_noIP\n";
print "wig_tracks_background:                   $wig_tracks_background\n";
print "calls_tracks:                            $calls_tracks\n";

#my $kbhit = <STDIN>;

#print "\n";
#print "Read the following parameters: \n\n";
#
#print "analysis_with_FDR:                        $analysis_with_FDR\n";
#print "\n";
#print "quick_window_scan_param_file:               $quick_window_scan_param_fname\n";
#print "calibrate_peak_shift_param_file:          $calibrate_peak_shift_param_fname\n";
#print "\n";
#print "ChIP_generate_profile_param_file:         $ChIP_generate_profile_param_fname\n";
#print "background_generate_profile_param_file:   $background_generate_profile_param_fname\n";
#if($analysis_with_FDR eq "yes"){
#    print "pseudo_ChIP_generate_profile_param_file:  $pseudo_ChIP_generate_profile_param_fname\n";
#}
#
#print "\n";
#print "ChIP_peak_caller_param_file:              $ChIP_peak_caller_param_fname\n";
#if($analysis_with_FDR eq "yes"){
#    print "pseudo_ChIP_peak_caller_param_file:       $pseudo_ChIP_peak_caller_param_fname\n";
#}
#
#print "calls_tracks_param_file:         	    $ChIP_calls_tracks_param_fname\n";
#
#print "\n";
#print "log_file:                                 $log_fname\n";
#print "output_file:                              $output_fname\n";

## QuEST execution

## some necessary checks/setups
if(-e $log_fname){
    unlink($log_fname);   
}
if(-e $output_fname){
    unlink($output_fname);
}
## some necessary checks/setups
#my $align_2_bin_batcher = $exec_path . "/run_align_2_bin_with_param_file.pl";
my $quick_window_scan_batcher = $exec_path . "/run_quick_window_scan_with_param_file.pl";
my $calibrate_peak_shift_batcher = $exec_path . "/run_calibrate_peak_shift_with_param_file.pl";
my $generate_profile_batcher = $exec_path . "/run_generate_profile_with_param_file.pl";
my $report_bad_control_regions_batcher = $exec_path . "/run_report_bad_control_regions_with_param_file.pl";

my $peak_caller_batcher = $exec_path . "/run_peak_caller_with_param_file.pl";
my $metrics_batcher = $exec_path . "/run_metrics_with_param_file.pl";
#my $qv_calculation_batcher = $exec_path . "/run_qv_calculation_with_param_file.pl";
my $filter_regions_batcher = $exec_path . "/filter_regions.pl";
my $profile_2_wig_batcher = $exec_path . "/run_profile_2_wig_with_param_file.pl";
my $regions_2_tracks_batcher = $exec_path . "/regions_2_tracks.pl";
my $align_2_bed_batcher = $exec_path . "/QuEST_align_2_bed.pl";
my $bed_graph_tracks_batcher = $exec_path . "/run_bin_align_2_bedgraph_with_param_file.pl";

if($bed_tracks_ChIP eq "true"){
    my $cur_system_command = "$align_2_bed_batcher -p $ChIP_reads_bed_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command &");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running ChIP_align_2_bed. Aborting.\n";
	exit(4);
    }
}

if($bed_tracks_RX_noIP eq "true"){
    my $cur_system_command = "$align_2_bed_batcher -p $RX_noIP_reads_bed_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command &");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running RX_noIP_align_2_bed. Aborting.\n";
	exit(4);
    }
}

if($bed_graph_tracks_ChIP eq "true"){
    my $cur_system_command = "$bed_graph_tracks_batcher -p $bed_graph_tracks_ChIP_param_fname\n";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running bed_graph_tracks_batcher. Aborting.\n";
	exit(4);
    }
}

if($generate_profile_ChIP eq "true" || $generate_profile_background eq "true" || $generate_profile_pseudo_ChIP eq "true"){
    my $scores_path = $analysis_path . "/scores";
    my $ChIP_scores_path = $scores_path . "/ChIP/";
    my $background_scores_path = $scores_path . "/background";
    my $pseudo_ChIP_scores_path = $scores_path . "/pseudo_ChIP";

    if( -d $scores_path){}else{ mkdir $scores_path or die "Failed to create the scores path $scores_path.\n"}
    if( -d $ChIP_scores_path){}else{ mkdir $ChIP_scores_path or die "Failed to create ChIP scores path $ChIP_scores_path.\n"}
    if( -d $pseudo_ChIP_scores_path){}else{ mkdir $pseudo_ChIP_scores_path or die "Failed to create pseudo ChIP scores path $pseudo_ChIP_scores_path.\n"}
    if( -d $background_scores_path){}else{ mkdir $background_scores_path or die "Failed to create background scores path $background_scores_path.\n"}
}

if( $quick_window_scan_ChIP eq "true"){
    my $cur_system_command = "$quick_window_scan_batcher -p $quick_window_scan_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running quick_window_scan. Aborting.\n";
	exit(4);
    }
}
if( $calibrate_peak_shift_ChIP eq "true" ){
    my $cur_system_command = "$calibrate_peak_shift_batcher -p $calibrate_peak_shift_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running calibrate_peak_shift. Aborting.\n";
	exit(4);
    }
}
if( $generate_profile_ChIP eq "true" ){
    my $cur_system_command = "$generate_profile_batcher -p $ChIP_generate_profile_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running generate_ChIP_profile. Aborting.\n";
	exit(4);
    }
}
if( $generate_profile_background eq "true" ){
    my $cur_system_command = "$generate_profile_batcher -p $background_generate_profile_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running generate_background_profile. Aborting.\n";
	exit(4);
    }
}
if( ($analysis_with_FDR eq "true") && $generate_profile_pseudo_ChIP eq "true"){
    my $cur_system_command = "$generate_profile_batcher -p $pseudo_ChIP_generate_profile_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running generate_pseudo_ChIP_profile. Aborting.\n";
	exit(4);
    }
}

if( $report_bad_control_regions eq "true"){
    my $cur_system_command = "$report_bad_control_regions_batcher -p $report_bad_control_regions_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running report_bad_control_regions batch script. Aborting.\n";
	exit(4);
    }
}

if( $peak_caller_ChIP_vs_background eq "true" ){
    my $cur_system_command = "$peak_caller_batcher -p $ChIP_peak_caller_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running peak_caller_ChIP_vs_background. Aborting.\n";
	exit(4);
    }
}

if( ($analysis_with_FDR eq "true") && $peak_caller_pseudo_ChIP_vs_background eq "true" ){
    my $cur_system_command = "$peak_caller_batcher -p $pseudo_ChIP_peak_caller_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running peak_caller_pseudo_ChIP_vs_background. Aborting.\n";
	exit(4);
    }
}

if( $metrics_ChIP eq "true"){
    my $cur_system_command  = "$metrics_batcher -p $metrics_ChIP_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running ChIP metrics Module. Aborting.\n";
	exit(4);
    }
}

if( $filter_ChIP_calls eq "true"){
    my $cur_system_command ="$filter_regions_batcher -p $filter_ChIP_calls_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
   
   if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running filter ChIP calls Module. Aborting.\n";
	exit(4);
    }
}

open log_file, "> $log_fname" || die "$log_fname: $\n";
open output_file, "> $output_fname" || die "$output_fname: $\n";

# creating QuEST summaries

my $ChIP_peak_calls_fname = "";
my $accepted_ChIP_peak_calls_fname = "";
my $rejected_ChIP_peak_calls_fname = "";

my $pseudo_ChIP_peak_calls_fname = "";

if(-e $ChIP_peak_caller_param_fname){
    open ChIP_peak_caller_param_file, "< $ChIP_peak_caller_param_fname";
    
    while( <ChIP_peak_caller_param_file> ){
	chomp;
	my @cur_fields = split(/ /, $_);
	if(scalar(@cur_fields) > 0){
	    if($cur_fields[0] eq "output_file"){
		$ChIP_peak_calls_fname = $cur_fields[2];
	    }
	}
    }    
    close ChIP_peak_caller_param_file;
}
else{
    print "Error in run_QuEST_with_param_file.pl: Failed to locate ChIP_peak_caller_param_fname: $ChIP_peak_caller_param_fname";
    print log_file"Error in run_QuEST_with_param_file.pl: Failed to locate ChIP_peak_caller_param_fname: $ChIP_peak_caller_param_fname";
}

if(-e $filter_ChIP_calls_param_fname){
    open filter_ChIP_calls_param_file, "< $filter_ChIP_calls_param_fname";
    
    while( <filter_ChIP_calls_param_file> ){
	chomp;
	my @cur_fields = split(/ /, $_);
	if(scalar(@cur_fields) > 0){
	    if($cur_fields[0] eq "output_file_prefix"){
		$accepted_ChIP_peak_calls_fname = $cur_fields[2] . ".accepted";
		$rejected_ChIP_peak_calls_fname = $cur_fields[2] . ".rejected";
	    }
	}
    }    
    close filter_ChIP_calls_param_file;
}
else{
    print "Error in run_QuEST_with_param_file.pl: Failed to locate $filter_ChIP_calls_param_fname\n";
    print log_file"Error in run_QuEST_with_param_file.pl: Failed to locate $filter_ChIP_calls_param_fname\n";
}

if( $analysis_with_FDR eq "true"){
    if(-e $pseudo_ChIP_peak_caller_param_fname){
	open pseudo_ChIP_peak_caller_param_file, "< $pseudo_ChIP_peak_caller_param_fname";
	
	while( <pseudo_ChIP_peak_caller_param_file> ){
	    chomp;
	    my @cur_fields = split(/ /, $_);
	    if(scalar(@cur_fields) > 0){
		if($cur_fields[0] eq "output_file"){
		    $pseudo_ChIP_peak_calls_fname = $cur_fields[2];
		}
	    }
	}    
	close pseudo_ChIP_peak_caller_param_file;
    }
    else{
	print "Error in run_QuEST_with_param_file.pl: Failed to locate pseudo_ChIP_peak_caller_param_fname: $pseudo_ChIP_peak_caller_param_fname";
	print log_file "Error in run_QuEST_with_param_file.pl: Failed to locate pseudo_ChIP_peak_caller_param_fname: $pseudo_ChIP_peak_caller_param_fname";
    }
}

# counting ChIP peaks
my $ChIP_peaks = -1;
my $ChIP_regions = -1;

#if( $peak_caller_ChIP_vs_background eq "true" ){
if(-e $ChIP_peak_calls_fname){
    $ChIP_peaks = 0;
    $ChIP_regions = 0;
    
    open ChIP_peak_calls_file, "< $ChIP_peak_calls_fname";
    
    while( <ChIP_peak_calls_file> ){
	my $cur_entry = $_;
	chomp($cur_entry);
	if(length($cur_entry) > 0){
	    if(substr($cur_entry, 0, 1) eq "R"){
		$ChIP_regions++;
	    }
	    elsif(substr($cur_entry, 0, 1) eq "P"){
		$ChIP_peaks++;
	    }
	}
    }
    close ChIP_peak_calls_file;
}
else{
    print "Failed to locate ChIP_peak_calls_file: $ChIP_peak_calls_fname\n";
    print log_file "Failed to locate ChIP_peak_calls_file: $ChIP_peak_calls_fname\n";
}

my $ChIP_peaks_accepted = -1;
my $ChIP_regions_accepted = -1;

if(-e $accepted_ChIP_peak_calls_fname){
    $ChIP_peaks_accepted = 0;
    $ChIP_regions_accepted = 0;
    open accepted_ChIP_peak_calls_file, "< $accepted_ChIP_peak_calls_fname";
    
    while( <accepted_ChIP_peak_calls_file> ){
	my $cur_entry = $_;
	chomp($cur_entry);
	if(length($cur_entry) > 0){
	    if(substr($cur_entry, 0, 1) eq "R"){
		$ChIP_regions_accepted++;
	    }
	    elsif(substr($cur_entry, 0, 1) eq "P"){
		$ChIP_peaks_accepted++;
	    }
	}
    }
    close accepted_ChIP_peak_calls_file;
}
else{
    print "Failed to locate ChIP_peak_calls_file: $accepted_ChIP_peak_calls_fname\n";
    print log_file "Failed to locate ChIP_peak_calls_file: $accepted_ChIP_peak_calls_fname\n";
}

my $ChIP_peaks_rejected = -1;
my $ChIP_regions_rejected = -1;

if(-e $rejected_ChIP_peak_calls_fname){
    $ChIP_peaks_rejected = 0;
    $ChIP_regions_rejected = 0;
    open rejected_ChIP_peak_calls_file, "< $rejected_ChIP_peak_calls_fname";
    
    while( <rejected_ChIP_peak_calls_file> ){
	my $cur_entry = $_;
	chomp($cur_entry);
	if(length($cur_entry) > 0){
	    if(substr($cur_entry, 0, 1) eq "R"){
		$ChIP_regions_rejected++;
	    }
	    elsif(substr($cur_entry, 0, 1) eq "P"){
		$ChIP_peaks_rejected++;
	    }
	}
    }
    close rejected_ChIP_peak_calls_file;
}
else{
    print "Failed to locate ChIP_peak_calls_file: $rejected_ChIP_peak_calls_fname\n";
    print log_file "Failed to locate ChIP_peak_calls_file: $rejected_ChIP_peak_calls_fname\n";
}

#}
#print "Found: $ChIP_regions regions and $ChIP_peaks\n";

# /counting ChIP peaks

# counting pseudo_ChIP peaks
my $pseudo_ChIP_peaks = -1;
my $pseudo_ChIP_regions = -1;
#if( $peak_caller_pseudo_ChIP_vs_background eq "true" ){
if( $analysis_with_FDR eq "true"){
    $pseudo_ChIP_peaks = 0;
    $pseudo_ChIP_regions = 0;
    if(-e $pseudo_ChIP_peak_calls_fname){
	open pseudo_ChIP_peak_calls_file, "< $pseudo_ChIP_peak_calls_fname";
	
	while( <pseudo_ChIP_peak_calls_file> ){
	    my $cur_call = $_;
	    chomp($cur_call);;
	    if(length($cur_call) > 0){
		if(substr($cur_call, 0 ,1 ) eq "R"){
		    $pseudo_ChIP_regions++;
		}
		elsif(substr($cur_call, 0, 1) eq "P"){
		    $pseudo_ChIP_peaks++;
		}
	    }
	}
	close pseudo_ChIP_peak_calls_file;
    }
    else{
	print "Failed to locate pseudo_ChIP_peak_calls_file: $pseudo_ChIP_peak_calls_fname\n";
	print log_file "Failed to locate pseudo_ChIP_peak_calls_file: $pseudo_ChIP_peak_calls_fname\n";
    }
}
#}
# /counting pseudo_ChIP peaks

# /creating QuEST summaries

print output_file "\#\#   please cite: \n";
print output_file "\#\#   Valouev A, Johnson DS, Sundquist A, Medina C, Anton E, Batzoglou S,\n";
print output_file "\#\#   Myers RM, Sidow A\n";
print output_file "\#\#   Genome-wide analysis of transcription factor binding sites based\n";
print output_file "\#\#   on ChIP-Seq data.\n";
print output_file "\#\#   Nat Methods. 2008 Sep; 5:(9):829-35\n";
print output_file "\n";

if($ChIP_peaks >= 0){
    print output_file "ChIP peaks: $ChIP_peaks\n";
}
else{
    print output_file "ChIP peaks: not available\n";
}

if($ChIP_peaks_accepted >= 0){
    print output_file "ChIP peaks accepted: $ChIP_peaks_accepted\n";
}
else{
    print output_file "ChIP peaks accepted: not available\n";
}

if($ChIP_peaks_rejected >= 0){
    print output_file "ChIP peaks rejected: $ChIP_peaks_rejected\n";
}
else{
    print output_file "ChIP peaks rejected: not available\n";
}

print output_file "\n";

if($ChIP_regions >= 0){
    print output_file "ChIP regions: $ChIP_regions\n";
}
else{
    print output_file "ChIP regions: not available\n";
}

if($ChIP_regions_accepted >= 0){
    print output_file "ChIP regions accepted: $ChIP_regions_accepted\n";
}
else{
    print output_file "ChIP regions accepted: not available\n";
}

if($ChIP_regions_rejected >= 0){
    print output_file "ChIP regions rejected: $ChIP_regions_rejected\n";
}
else{
    print output_file "ChIP regions rejected: not available\n";
}

if( $analysis_with_FDR eq "true"){
    print output_file "\n";
    if($pseudo_ChIP_peaks >= 0){
	print output_file "pseudo_ChIP peaks: $pseudo_ChIP_peaks\n";
	print output_file "pseudo_ChIP regions: $pseudo_ChIP_regions\n";

	print output_file "\n";

	if($ChIP_peaks != 0){
	    printf output_file "ChIP peak FDR estimate: %.2f percent\n", 100 * $pseudo_ChIP_peaks / $ChIP_peaks;
	}
	else{
	    print output_file "ChIP peak FDR estimate: NA\n";
	}
    }
    else{
	print "ChIP peak FDR estimate: not available\n";
    }
    if($pseudo_ChIP_regions >= 0){
	if($ChIP_regions != 0){
	    printf output_file "ChIP regions FDR estimate: %.2f percent\n", 100 * $pseudo_ChIP_regions / $ChIP_regions;
	}
	else{
	    print output_file "ChIP regions FDR estimate: NA\n";
	}
    }
    else{
	print "ChIP regions FDR estimate: not available\n";
    }
}



close output_file;
close log_file;

## output tracks

if( $calls_tracks eq "true"){
    my $cur_system_command = "$regions_2_tracks_batcher -p $ChIP_calls_tracks_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when calls_2_tracks batcher. Aborting.\n";
	exit(4);
    }
}

## ChIP
if($wig_tracks_ChIP eq "true"){
    my $ChIP_normalized_wig_error_code = system("$profile_2_wig_batcher -p $ChIP_normalized_wig_param_fname");
    if($ChIP_normalized_wig_error_code != "0"){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $ChIP_normalized_wig_error_code) when generate normalized ChIP wig profile. Aborting.\n";
	exit(4);

    }
    my $ChIP_wig_error_code = system("$profile_2_wig_batcher -p $ChIP_unnormalized_wig_param_fname");
    if($ChIP_wig_error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $ChIP_wig_error_code) when generate normalized ChIP wig profile. Aborting.\n";
	exit(4);
	
    }
}
## /ChIP

## background
if($wig_tracks_background eq "true"){
    my $cur_system_command = "$profile_2_wig_batcher -p $background_normalized_wig_param_fname";
    print "system_command: $cur_system_command\n";
    my $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running wig_gracks_normalized_background_batcher. Aborting.\n";
	exit(4);
    }    

    $cur_system_command = "$profile_2_wig_batcher -p $background_unnormalized_wig_param_fname";
    print "system_command: $cur_system_command\n";
    $error_code = system("$cur_system_command");
    if($error_code != 0){
	print "run_QuEST_with_param_file.pl: Caught an error message (code $error_code) when running wig_tracks_unnormalized_background_batcher. Aborting.\n";
	exit(4);
    }    
}
## /background

## /output tracks

if(-d $scores_dir){
    use File::Path;
    rmtree($scores_dir);
}

## /QuEST execution

