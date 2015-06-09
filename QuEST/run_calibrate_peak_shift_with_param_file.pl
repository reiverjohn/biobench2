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

## This program is a wrapper that runs calibrate_peak_shift
## on the entire genome one contig at a time

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    run_calibrate_peak_shift_with_param_file.pl
	
    This program is a wrapper that runs the calibrate_peak_shift
    
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

$exec_path = $exec_path."/calibrate_peak_shift";
## /setting optional arguments

if( ! -e $exec_path ){
    print "Error in calibrate peak shift master script.\n";
    print "Couldn't locate executable $exec_path. Aborting.\n";
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

my $silent_flag = "off";
my $genome_table_fname =                      "** missing **";
my $regions_fname =                     "** missing **";
my $ChIP_bin_align_path =                "** missing **";
my $output_fname =                        "** missing **";

my $kde_bandwidth = "";                     # optional
my $region_width = "";                     # optional
my $top_regions_passed = "";                # optional
my $dist_threshold = "";
my $correlation_lower_threshold = "";
my $peak_shift_lower_threshold = "";
my $peak_shift_upper_threshold = "";
my $os_upper_threshold = "";

my $estimation_method = "";

my $alternative_ps = "undef";


my $verbose = "yes";

for(my $i=0; $i<scalar(@params); $i++){
    my $cur_param = $params[$i];
    chomp($cur_param);
    
    my @cur_par_fields = split(/ /, $cur_param);
    if(scalar(@cur_par_fields >= 2)){
	my $cur_par_name = $cur_par_fields[0];	
	if ($cur_par_name eq "regions_file"){
	    $regions_fname = $cur_par_fields[2];	
	}
	elsif($cur_par_name eq "genome_table"){
	    $genome_table_fname = $cur_par_fields[2];	
	}
	elsif($cur_par_name eq "bin_align_path"){
	    $ChIP_bin_align_path = $cur_par_fields[2];	
	}
	elsif($cur_par_name eq "output_file"){
	    $output_fname = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "kde_bandwidth"){
	    $kde_bandwidth = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "region_width"){
	    $region_width = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "top_regions_passed"){
	    $top_regions_passed = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "dist_threshold"){
	    $dist_threshold = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "correlation_lower_threshold"){
	    $correlation_lower_threshold = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "peak_shift_lower_threshold"){
	    $peak_shift_lower_threshold = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "peak_shift_upper_threshold"){
	    $peak_shift_upper_threshold = $cur_par_fields[2];
	}
	elsif($cur_par_name eq "os_upper_threshold"){
	    $os_upper_threshold = $cur_par_fields[2];
	}
	elsif ($cur_par_name eq "estimation_method"){
	    $estimation_method = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "alternative_ps"){
	    $alternative_ps = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "verbose"){
	    $verbose = $cur_par_fields[2];	
	}
	else{
	    if($params[$i] ne ""){
		print "Warning: unrecognized parameter: $cur_par_name";
	    }
	}
    }
}

print "Read the following parameters: \n\n";

print "genome_table:                      $genome_table_fname\n";
print "regions_file:                      $regions_fname\n";
print "bin_align_path:                    $ChIP_bin_align_path\n";
print "output_file:                       $output_fname\n";
print "kde_bandwidth:                     $kde_bandwidth\n";
print "region_width:                      $region_width\n";
print "top_regions_passed:                $top_regions_passed\n";
print "dist_threshold:                    $dist_threshold\n";
print "correlation_lower_threshold:       $correlation_lower_threshold\n";
print "peak_shift_lower_threshold:        $peak_shift_lower_threshold\n";
print "peak_shift_upper_threshold:        $peak_shift_upper_threshold\n";
print "os_upper_threshold:                $os_upper_threshold\n";
print "estimation_method:                 $estimation_method\n";
print "verbose:                           $verbose\n";
print "\n";


my $output_fname_tmp = $output_fname . ".with_stats";

my $mandatory_param_string = "genome_table=$genome_table_fname bin_align_path=$ChIP_bin_align_path regions_file=$regions_fname output_file=$output_fname_tmp kde_bandwidth=$kde_bandwidth region_width=$region_width top_regions_passed=$top_regions_passed dist_threshold=$dist_threshold correlation_lower_threshold=$correlation_lower_threshold peak_shift_lower_threshold=$peak_shift_lower_threshold peak_shift_upper_threshold=$peak_shift_upper_threshold os_upper_threshold=$os_upper_threshold estimation_method=$estimation_method";

my $system_command = $exec_path . " " . $mandatory_param_string; # . $optional_param_string;

print "system_command: $system_command\n";
my $error_code = system("$system_command\n");

if($error_code != 0){
    print "run_calibrate_peak_shift_with_param_file.pl: Caught an error message (code $error_code), passing the error message to the top script.\n";
    exit(2);
}

if(!-e $output_fname_tmp){
    print "Error in run_calibrate_peak_shift_with_param_file.pl: Failed to locate $output_fname_tmp. Aborting.\n";
    exit(1);
}

open output_file_tmp, "< $output_fname_tmp";

#my $regions = -1;
my $used_regions = -1;
my $peak_shift_estimate = -1;
my $estimation_method_detected = "";

while( <output_file_tmp> ){
    chomp;
    if(length($_) > 0){
	if( substr($_,0,1) ne "#" ){
	    my @cur_entry_fields = split(/ /, $_);
	    if( scalar(@cur_entry_fields) >= 2){
		#if($cur_entry_fields[0] eq "regions:"){
		#    $regions = $cur_entry_fields[1];
		#}
		if($cur_entry_fields[0] eq "used_regions:"){
		    $used_regions = $cur_entry_fields[1];
		}
		elsif($cur_entry_fields[0] eq "peak_shift_estimate:"){
		    $peak_shift_estimate = $cur_entry_fields[1];
		}
		elsif($cur_entry_fields[0] eq "estimation_method:"){
		    $estimation_method_detected = $cur_entry_fields[1];
		}
	    }
	}
    }
}

close output_file_tmp;

my $recommended_used_regions = 10;
my $peak_shift_recommended_min = 20;
my $peak_shift_recommended_value = 50;
my $peak_shift_recommended_max = 130;


if( $verbose eq "yes" ){
    my $suggest_override = "false";
#    if( $used_regions <= $recommended_used_regions ){
#	print "\nQuEST has determined that the fraction of used regions $used_regions for peak\n";
#	print "estimate is low. It does not mean that the estimate is bad, but\n";
#	print "the peak shift estimate may not be very accurate.\n";
#	$suggest_override = "true";
#    }


    if( $peak_shift_estimate eq "NA"){
	print "\nQuEST failed to identify regions to determine peak shift estimate.\n";
	$suggest_override = "true";
    }
    elsif( $peak_shift_estimate < $peak_shift_recommended_min ){
	printf "\nQuEST has determined that the peak shift estimate ( %.0f )\n", $peak_shift_estimate;
	print "is too low. This may be normal in some cases, but also may be indicative\n";
	print "of read mapping artifacts into repetetive regions. You might consider\n";
	print "revising the list of candidate regions to exclude the ones falling into the\n";
	print "regions annotated as repeats or try using median-based estimate for the peak shift\n";
	print "to achieve a more robust estimate.\n";
	$suggest_override = "true";
    }
    elsif( $used_regions <= $recommended_used_regions ){
	print "\nQuEST has determined that the number of used regions ( $used_regions ) for peak\n";
	print "estimate is low. It does not mean that the estimate is bad, but\n";
	print "the peak shift estimate may not be very accurate.\n";
	$suggest_override = "true";
    }
	
    if( $suggest_override eq "true" ){
	if($alternative_ps eq "undef"){
	    my $cur_answer = "";
	    while($cur_answer ne "y" && $cur_answer ne "n" ){
		print "Do you want to override the value of peak shift? (y/n): ";
		$cur_answer = <STDIN>;	    
		chomp($cur_answer);
	    }
	    if( $cur_answer eq "y"){
		$estimation_method_detected = "user_input";
		print "Please enter the value of peak shift (recommended $peak_shift_recommended_value): ";
		$cur_answer = <STDIN>;
		chomp($cur_answer);
		$peak_shift_estimate = $cur_answer;
	    }
	    else{
		$peak_shift_estimate = $peak_shift_recommended_value;
	    }
	}
	else{
	    $peak_shift_estimate = $alternative_ps;
	}
    }

    
    if( $peak_shift_estimate > $peak_shift_recommended_max ){
	printf "\nQuEST has determined that the peak shift estimate ( %.0f )\n", $peak_shift_estimate;
	print "is too high. This may be normal in some cases especially \n";
	print "if the insert size in the library is large, however \n";
	print "typically values around $peak_shift_recommended_value are expected.\n";
	print "regions annotated as repeats or try using median-based estimate for the peak shift.\n";
	my $cur_answer = "";
	while($cur_answer ne "y" && $cur_answer ne "n" ){
	    print "Do you want to override the value of peak shift? (y/n): ";
	    $cur_answer = <STDIN>;
	    chomp($cur_answer);
	}
	if( $cur_answer eq "y"){
	    $estimation_method_detected = "user_input";
	    my $number_given = "false";

	    while($number_given eq "false"){
		print "Please enter the value of peak shift (recommended $peak_shift_recommended_value): ";
		$cur_answer = <STDIN>;
		chomp($cur_answer);
		if($cur_answer =~ /^\d+$/ ){
		    $number_given = "true";
		}
		else{
		    print "Please provide a positive integer (no spaces at all)\n";
		}
	    }
	    $peak_shift_estimate = $cur_answer;
	}	
    }
}
else{
    if($peak_shift_estimate eq "NA"){
	$peak_shift_estimate = $peak_shift_recommended_value;
	$estimation_method_detected = "failed_to_find_regions_so_used_suggested_peak_shift_value";
    }
}

if(-e $output_fname){
    unlink($output_fname);
}

#if(-e $output_fname_tmp){
#    unlink($output_fname_tmp);
#}

open output_file, "> $output_fname";
print output_file "used_regions: $used_regions\n";
#print output_file "used: $used_regions\n";
print output_file "peak_shift_estimate: $peak_shift_estimate\n";
print output_file "estimation_method: $estimation_method_detected\n";
close output_file;
