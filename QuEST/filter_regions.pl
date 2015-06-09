#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    filter_regions.pl

    -----------------------------
    mandatory parameters:
    
    -i regions_file
    -o output_file_prefix            
    
    or 

    -p parameter_file

    -----------------------------
    optional parameters:
    
    -h (help)

    -ps_upper  <peak_shift upper thershold>
    -ps_lower  <peak_shift lower threshold>

    -cor_upper <correlation upper threshold>
    -cor_lower <correlation lower threshold>

    -log10_qv  <log10 qv threshold>

};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

my $undef = "undef";

## mandatory arguments

my $regions_fname = "";

my $bad_control_regions_fname = $undef;
my $using_bad_control_regions = "false";

my $output_fname_prefix = "";
my $accept_fname = "";
my $reject_fname = "";

## optional arguments
my $correlation_upper_threshold = $undef;
my $correlation_lower_threshold = $undef;

my $peak_shift_lower_threshold = $undef;
my $peak_shift_upper_threshold = $undef;

my $log10_qv_threshold = $undef;

my $report_peaks = "all"; 

## can be "all" if all peaks within the region to be reported
## can be "best" where it will output a single peak for each enriched region

my $filter_out_empty_regions = "true"; ## if yes, will filter out regions without any peaks


my $param_fname = $undef;
my $report_fname = $undef;


## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-i') {$regions_fname = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_fname_prefix = shift @ARGV;}

    elsif ( $this_arg eq '-ps_upper') {$peak_shift_upper_threshold = shift @ARGV;}
    elsif ( $this_arg eq '-ps_lower') {$peak_shift_lower_threshold = shift @ARGV;}

    elsif ( $this_arg eq '-cor_upper') {$correlation_upper_threshold = shift @ARGV;}
    elsif ( $this_arg eq '-cor_lower') {$correlation_lower_threshold = shift @ARGV;}
    
    elsif ( $this_arg eq '-log10_qv') {$log10_qv_threshold = shift @ARGV;}

    elsif ( $this_arg eq '-p') {$param_fname = shift @ARGV;}

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}



## print parameters

if($param_fname ne "undef"){
    open param_file, "< $param_fname" || die "$param_fname: $!\n";
    while( <param_file> ){
	chomp;
	my $cur_par_entry = $_;
	my @cur_par_entry_fields = split(/ /,$cur_par_entry);
	if(scalar(@cur_par_entry_fields) >= 3){
	    my $cur_par_name = $cur_par_entry_fields[0];
	    my $cur_par_value = $cur_par_entry_fields[2];

	    if($cur_par_name eq "input_file"){
		$regions_fname = $cur_par_value;
	    }
	    elsif($cur_par_name eq "output_file_prefix"){
		$output_fname_prefix = $cur_par_value;
	    }
	    elsif($cur_par_name eq "log10_qv"){
		$log10_qv_threshold = $cur_par_value;
	    }
	    elsif($cur_par_name eq "ps_upper_threshold"){
		$peak_shift_upper_threshold = $cur_par_value;
	    }
	    elsif($cur_par_name eq "ps_lower_threshold"){
		$peak_shift_lower_threshold = $cur_par_value;
	    }
	    elsif($cur_par_name eq "cor_upper_threshold"){
		$correlation_upper_threshold = $cur_par_value;
	    }
	    elsif($cur_par_name eq "cor_lower_threshold"){
		$correlation_lower_threshold = $cur_par_value;
	    }
	    elsif($cur_par_name eq "report_peaks"){
		$report_peaks = $cur_par_value;
	    }
	    elsif($cur_par_name eq "filter_out_empty_regions"){
		$filter_out_empty_regions = $cur_par_value;
	    }
	    elsif($cur_par_name eq "bad_control_regions_file"){
		$bad_control_regions_fname = $cur_par_value;
	    }
	    elsif($cur_par_name eq "report_file"){
		$report_fname = $cur_par_value;
	    }
	    else{
		print "Bad parameter $cur_par_name\n";
		exit;
	    }
	}
    }
    close param_file;
}

if ( $regions_fname eq ""){
    die "you should specify the regions file\n";
}
if( $output_fname_prefix eq ""){
    die "you should specify output file prefix\n";
}

$reject_fname = $output_fname_prefix . ".rejected";
$accept_fname = $output_fname_prefix . ".accepted";

if( -e $bad_control_regions_fname){
    $using_bad_control_regions = "true";
}
else{
    $using_bad_control_regions = "false";
}

print "\n-----------------\n\n";
print "regions file                  :    $regions_fname\n";
print "bad_control_regions_file      :    $regions_fname\n";
print "using_bad_control_regions     :    $using_bad_control_regions\n";
print "output file prefix            :    $output_fname_prefix\n";
print "accept file                   :    $accept_fname\n";
print "reject file                   :    $reject_fname\n";
print "report file                   :    $report_fname\n";
print "\n";

print "ps_upper          :    $peak_shift_upper_threshold\n";
print "ps_lower          :    $peak_shift_lower_threshold\n";
print "\n";

print "cor_upper         :    $correlation_upper_threshold\n";
print "cor_lower         :    $correlation_lower_threshold\n";
print "\n";

print "log10_qv          :    $log10_qv_threshold\n";
print "\n";

print "report_peaks      :    $report_peaks\n";
print "\n";

print "filter_out_empty_regions      :    $filter_out_empty_regions\n";
print "\n";

print "\n-----------------\n\n";

if(-e $report_fname){
    unlink($report_fname);
}
open report_file, "> $report_fname" || die "$report_fname: $!\n";

my @bad_control_regions;
my $bad_control_regions_counter = 0;

if($using_bad_control_regions eq "true"){
    open bad_control_regions_file, "< $bad_control_regions_fname" || die "$bad_control_regions_fname: $!\n";
    while(<bad_control_regions_file>){
	chomp;
	my $cur_region = $_;
	if(substr($cur_region,0,1) eq "R"){
	    my @cur_region_fields = split(/ /, $cur_region);
	    if(scalar(@cur_region_fields) == 7){
		my @cur_region_range_fields = split(/-/, $cur_region_fields[2]);
		if(scalar(@cur_region_range_fields) == 2){
		    my $new_bad_control_region = {
			entry => $cur_region,
			chrom => $cur_region_fields[1],
			start => $cur_region_range_fields[0],
			end => $cur_region_range_fields[1],
		    };
		    $bad_control_regions[$bad_control_regions_counter] = $new_bad_control_region;
		    $bad_control_regions_counter++;
		}
	    }
	}
    }
    close bad_control_regions_file;
}

printf "Read %d bad control regions.\n", scalar(@bad_control_regions);


if(-e $accept_fname){ unlink("$accept_fname");}
open accept_file, "> $accept_fname" || die "$accept_fname: $!\n";


if(-e $reject_fname){ unlink("$reject_fname");}
open reject_file, "> $reject_fname" || die "$reject_fname: $!\n";

open regions_file, "< $regions_fname" || die "$regions_fname: $!\n";

my $regions_passed = 0;
my $regions_failed = 0;

my $regions_overlapping_with_bad_control_regions = 0;

my $peaks_failed_within_regions_passed = 0;
my $peaks_passed_within_regions_passed = 0;

while(<regions_file>){
    chomp;
    my $cur_entry = $_;

    if( length($cur_entry) > 0 ){
	if(substr($cur_entry, 0, 1) ne "#"){
	    
	    if(substr($cur_entry,0,1) eq "R"){		
		my $cur_region = $cur_entry;		
		my @cur_region_fields = split(/ /, $cur_region);
		if(scalar(@cur_region_fields) >= 17){
		    my $region_failed_reason = $undef;
		    my $cur_region_chrom = $cur_region_fields[1];
		    my @cur_region_range_fields = split(/-/, $cur_region_fields[2]);

		    my $cur_region_start = $cur_region_range_fields[0];
		    my $cur_region_end = $cur_region_range_fields[1];

		    my $cur_region_ps = $cur_region_fields[18];
		    my $cur_region_cor = $cur_region_fields[20];
		    my $cur_region_log10_qv = - $cur_region_fields[22];

		    my @local_peaks;
		    my $local_peaks_counter = 0;
		    my $end_of_region = "false";
#		    my @local_peaks_passed;

		    while($end_of_region eq "false"){
			my $cur_peak_entry = <regions_file>;
			chomp($cur_peak_entry);
			
			if( length($cur_peak_entry) > 0){
			    if( substr($cur_peak_entry,0,1) eq "#"){
			    }
			    elsif( substr($cur_peak_entry,0,1) eq "P"){
				my $cur_peak = {
				    entry => $cur_peak_entry,
				    passed => "undef",
				    score => "undef",
				    ps => "undef",
				    cor => "undef",
				    log10_qv => "undef",
				};
				$local_peaks[$local_peaks_counter] = $cur_peak;
				$local_peaks_counter++;
			    }
			    else{
				$end_of_region = "true";
			    }
			}
			else{
			    $end_of_region = "true";
			}
		    }
		    
		    my $region_log10_qv_passed = "undef";
		    
		    my $region_ps_upper_passed = "undef";
		    my $region_ps_lower_passed = "undef";
		    
		    my $region_cor_upper_passed = "undef";
		    my $region_cor_lower_passed = "undef";
		    
		    if($log10_qv_threshold ne $undef){
			if($cur_region_log10_qv < $log10_qv_threshold){
			    $region_log10_qv_passed = "true";
			}
			else{
			    $region_log10_qv_passed = "false";
			}
		    }
		    else{
			$region_log10_qv_passed = "true"; # pass everything if no threshold specified
		    }
		    
		    
		    if($peak_shift_lower_threshold ne $undef){
			if($cur_region_ps >= $peak_shift_lower_threshold){
			    $region_ps_lower_passed = "true";
			}
			else{
			    $region_ps_lower_passed = "false";
			}
		    }
		    else{ $region_ps_lower_passed = "true"; }  # pass everything if no threshold specified
		    
		    if($peak_shift_upper_threshold ne $undef){
			if($cur_region_ps <= $peak_shift_upper_threshold){
			    $region_ps_upper_passed = "true";
			}
			else{
			    $region_ps_upper_passed = "false";
			}
		    }
		    else{ $region_ps_upper_passed = "true"; }  # pass everything if no threshold specified
		    
		    
		    if($correlation_lower_threshold ne $undef){
			if($cur_region_cor eq "nan"){
			    $region_cor_lower_passed = "false";
			}
			elsif($cur_region_cor > $correlation_lower_threshold){
			    $region_cor_lower_passed = "true";
			}
			else{
			    $region_cor_lower_passed = "false";
			}
		    }
		    else{ $region_cor_lower_passed = "true"; }  # pass everything if no threshold specified
		    
		    if($correlation_upper_threshold ne $undef){
			if($cur_region_cor eq "nan"){
			    $region_cor_lower_passed = "false";
			}
			elsif($cur_region_cor < $correlation_upper_threshold){
			    $region_cor_upper_passed = "true";
			}
			else{
			    $region_cor_upper_passed = "false";
			}
		    }
		    else{ $region_cor_upper_passed = "true"; }  # pass everything if no threshold specified
		    
		    my $region_passed = "true";
		    if($region_log10_qv_passed ne "true" || 
		       $region_ps_upper_passed ne "true" || $region_ps_lower_passed ne "true" ||
		       $region_cor_upper_passed ne "true" || $region_cor_lower_passed ne "true"){
			$region_passed = "false";
		    }
		    
		    if($region_log10_qv_passed eq "undef" ||
		       $region_ps_upper_passed eq "undef" || $region_ps_lower_passed eq "undef" ||
		       $region_cor_upper_passed eq "undef" || $region_cor_lower_passed eq "undef"){
			
			$region_passed = "undef";
		    }
		    
		    if($region_passed eq "undef"){
			print "region_passed: undef\n\n";
			print "$cur_region\n\n";

			print "region_log10_qv_passed: $region_log10_qv_passed\n";
			print "region_ps_upper_passed: $region_ps_upper_passed\n";
			print "region_ps_lower_passed: $region_ps_lower_passed\n";
			print "region_cor_upper_passed: $region_cor_upper_passed\n";
			print "region_cor_lower_passed: $region_cor_lower_passed\n";

			print "\n";
			print "cur_region_log10_qv: $cur_region_log10_qv\n";
			print "cur_region_ps: $cur_region_ps\n";
			print "\n";
			exit(4);
		    }

		    
		    for(my $i=0; $i<scalar(@local_peaks); $i++){
			my $cur_peak_entry = $local_peaks[$i]->{entry};
			my @cur_peak_entries = split(/ /, $cur_peak_entry);
			if( scalar(@cur_peak_entries) >= 17){
			    my $cur_peak_ps = $cur_peak_entries[12];
			    my $cur_peak_cor = $cur_peak_entries[14];
			    my $cur_peak_log10_qv = - $cur_peak_entries[16];
			    my $cur_peak_score = $cur_peak_entries[4];

			    $local_peaks[$i]->{ps} = $cur_peak_ps;
			    $local_peaks[$i]->{cor} = $cur_peak_cor;
			    $local_peaks[$i]->{log10_qv} = $cur_peak_log10_qv;
			    $local_peaks[$i]->{score} = $cur_peak_score;
			    
			    my $peak_log10_qv_passed = $undef;
			    
			    my $peak_ps_upper_passed = "undef";
			    my $peak_ps_lower_passed = "undef";
			    
			    my $peak_cor_upper_passed = "undef";
			    my $peak_cor_lower_passed = "undef";
			    
			    if($log10_qv_threshold ne $undef){
				if($cur_peak_log10_qv eq "inf"){
				}
				elsif($cur_peak_log10_qv < $log10_qv_threshold){
				    $peak_log10_qv_passed = "true";
				}
				else{
				    $peak_log10_qv_passed = "false";
				}
			    }
			    else{ $peak_log10_qv_passed = "true" } # pass everything if no threshold specified
			    
			    
			    if($peak_shift_lower_threshold ne $undef){
				if($cur_peak_ps > $peak_shift_lower_threshold){
				    $peak_ps_lower_passed = "true";
				}
				else{
				    $peak_ps_lower_passed = "false";
				}
			    }
			    else{ $peak_ps_lower_passed = "true"; } # pass everything if no threshold specified
			    
			    if($peak_shift_upper_threshold ne $undef){
				if($cur_peak_ps < $peak_shift_upper_threshold){
				    $peak_ps_upper_passed = "true";
				}
				else{
				    $peak_ps_upper_passed = "false";
				}
			    }
			    else{ $peak_ps_upper_passed = "true"; } # pass everything if no threshold specified
			    
			    if($correlation_lower_threshold ne $undef){
				if($cur_peak_cor eq "nan"){
				    $peak_cor_lower_passed = "false";
				}
				elsif($cur_peak_cor > $correlation_lower_threshold){
				    $peak_cor_lower_passed = "true";
				}
				else{
				    $peak_cor_lower_passed = "false";
				}
			    }
			    else{ $peak_cor_lower_passed = "true"; } # pass everything if no threshold specified
			    
			    if($correlation_upper_threshold ne $undef){
				if($cur_peak_cor eq "nan"){
				    $peak_cor_lower_passed = "false";
				}
				elsif($cur_peak_cor < $correlation_upper_threshold){
				    $peak_cor_upper_passed = "true";
				}
				else{
				    $peak_cor_upper_passed = "false";
				}
			    }
			    else{ $peak_cor_upper_passed = "true"; } # pass everything if no threshold specified
			    
			    my $peak_passed = "true";
			    if($peak_log10_qv_passed ne "true" ||
			       $peak_ps_upper_passed ne "true" || $peak_ps_lower_passed ne "true" ||
			       $peak_cor_upper_passed ne "true" || $peak_cor_lower_passed ne "true"){
				$peak_passed = "false";
			    }
			    
			    if($peak_log10_qv_passed eq "undef" || 
			       $peak_ps_upper_passed eq "undef" || $peak_ps_lower_passed eq "undef" ||
			       $peak_cor_upper_passed eq "undef" || $peak_cor_lower_passed eq "undef"){
				$peak_passed = "undef";
			    }
			    
			    if($peak_passed eq "undef"){
				print "peak_passed: undef\n";
				print "$cur_peak_entry\n";
				exit(4);
			    }
			    
			    if($peak_passed eq "true"){
				#$peaks_passed_within_regions_passed++;
				$local_peaks[$i]->{passed} = "true";
			    }
			    else{
				$peaks_failed_within_regions_passed++;
				$local_peaks[$i]->{passed} = "false";
				# fail this peak
			    }
			}
		    }
			

		    #if at least one peak passes, pass the region, but reject failed peaks

		    if($report_peaks eq "best"){
			my $best_peak_ind = -1;
			my $best_peak_score = -1;
			for(my $i=0; $i<scalar(@local_peaks); $i++){
			    if($local_peaks[$i]->{score} > $best_peak_score){
				$best_peak_score = $local_peaks[$i]->{score};
				$best_peak_ind = $i;
			    }
			    #$local_peaks[$i]->{passed} = "false";
			}
			
			for(my $i=0; $i<scalar(@local_peaks); $i++){
			    if($i != $best_peak_ind){
				$local_peaks[$i]->{passed} = "false";
			    }
			}
		    }

		    if($filter_out_empty_regions eq "true"){
			$region_passed = "false";
		    }
		    
		    for(my $i=0; $i<scalar(@local_peaks); $i++){
			if($local_peaks[$i]->{passed} eq "true"){
			    $region_passed = "true";
			}
		    }

		    
		    
		    for(my $i=0; $i<scalar(@local_peaks); $i++){
			if($local_peaks[$i]->{passed} eq "true"){
			    $region_passed = "true";
			}
		    }

		    ## comparing against contaminating regions
		    for(my $i=0; $i<scalar(@bad_control_regions); $i++){
			if($cur_region_chrom eq $bad_control_regions[$i]->{chrom}){
			    # if there are any common bases, mark the region as tainted
			    my $bad_region = $bad_control_regions[$i]->{entry};
			    my $bad_region_start = $bad_control_regions[$i]->{start};
			    my $bad_region_end = $bad_control_regions[$i]->{end};
			    my $overlap_exists = "false";

			    if($cur_region_start <= $bad_region_start){
				if($cur_region_end >= $bad_region_start){
				    $overlap_exists = "true"
				}
			    }
			    else{
				if($cur_region_start <= $bad_region_end){
				    $overlap_exists = "true";
				}
			    }
			    if($overlap_exists eq "true"){
				$region_passed = "false";
				if($region_failed_reason eq $undef){
				    $regions_overlapping_with_bad_control_regions++;
				}
				$region_failed_reason = "overlaps_with_bad_control_regions\n";
				print "Removed: ";
				print "$cur_region\n";
				print "Bad region: $bad_region\n";
				print "\n";
				#exit(0);
			    }
			}
		    }
		    ## /comparing against contaminating regions

		    if($region_passed eq "true"){
			my $cur_failed_peaks = 0;
			print accept_file "$cur_region\n";
			$regions_passed++;
			for(my$i=0; $i<scalar(@local_peaks); $i++){
			    if($local_peaks[$i]->{passed} eq "true"){
				printf accept_file "%s\n", $local_peaks[$i]->{entry};
				$peaks_passed_within_regions_passed++;
			    }
			    else{
				printf reject_file "%s\n", $local_peaks[$i]->{entry};
				$peaks_failed_within_regions_passed++;
				$cur_failed_peaks++;
			    }
			}
			if($cur_failed_peaks > 0){
			    print reject_file "\n";
			}
			print accept_file "\n";
		    }
		    else{
			$regions_failed++;
			print reject_file "$cur_region\n";
			for(my $i=0; $i<scalar(@local_peaks); $i++){
			    printf reject_file "%s\n", $local_peaks[$i]->{entry};
			}
			print reject_file "\n";
		    }
		}		
	    }
	    else{
		print "skipping entry: $cur_entry\n";
		print reject_file "$cur_entry\n";
		#fail this
	    }
	}
    }
}

print report_file "regions passed : $regions_passed\n";
print report_file "regions failed : $regions_failed\n";
print report_file "\n";

print report_file "peaks passed within passed regions : $peaks_passed_within_regions_passed\n";
print report_file "peaks failed within passed regions : $peaks_failed_within_regions_passed\n";

print report_file "\n";
print report_file "regions failed because they overlap with enriched control regions : $regions_overlapping_with_bad_control_regions\n";

print "regions passed : $regions_passed\n";
print "regions failed : $regions_failed\n";
print "\n";

print "peaks passed within passed regions : $peaks_passed_within_regions_passed\n";
print "peaks failed within passed regions : $peaks_failed_within_regions_passed\n";
print "\n";

close regions_file;
close accept_file;
close reject_file;
close report_file;
