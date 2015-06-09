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

## This program is a wrapper that runs peak_caller program
## on the entire genome one contig at a time

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    run_report_bad_control_regions_with_param_file.pl
	
    This program is a wrapper for the reprort_bad_control_regions executable
    
    -----------------------
    mandatory parameters:
    -p <params_file>               a file containing parameters calibrate_peak_shift

    -----------------------
    optional parameters:
    -e <exec_path>                 path to the calibrate_peak_shift executable
    -h                             to display this help

};

my $citation = qq{\#\#   please cite: 
\#\#   Valouev A, Johnson DS, Sundquist A, Medina C, Anton E, Batzoglou S,
\#\#   Myers RM, Sidow A. 
\#\#   Genome-wide analysis of transcription factor binding sites based on
\#\#   ChIP-Seq data.
\#\#   Nat Methods 5, 829-834 (2008)
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

$exec_path = $exec_path."/report_bad_control_regions";
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

my $background_score_path       = $missing;
my $output_fname                = $missing;
my $region_seeding_threshold    = $missing;
my $region_extension_threshold  = $missing;
my $score_normalizer            = $missing;

my $genome_table_fname =            "** missing **";


for(my $i=0; $i<scalar(@params); $i++){
    my $cur_param = $params[$i];
    chomp($cur_param);
    
    my @cur_par_fields = split(/ /, $cur_param);
    if(scalar(@cur_par_fields >= 2)){
	my $cur_par_name = $cur_par_fields[0];
	if ($cur_par_name eq "background_score_path"){
	    $background_score_path = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "output_file"){
	    $output_fname = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "region_seeding_threshold"){
	    $region_seeding_threshold = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "region_extension_threshold"){
	    $region_extension_threshold = $cur_par_fields[2];	
	}
	elsif($cur_par_name eq "genome_table"){
	    $genome_table_fname = $cur_par_fields[2];
	}
	elsif ($cur_par_name eq "score_normalizer"){
	    $score_normalizer = $cur_par_fields[2];	
	}
	else{
	    if($params[$i] ne ""){
		print "Warning: unrecognized parameter: $cur_par_name\n";
		exit(1);
	    }
	}
    }
}

print "Read the following parameters: \n\n";

print "background_score_path:             $background_score_path\n";
print "output_file:                       $output_fname\n";
print "region_seeding_threshold:          $region_seeding_threshold\n";
print "region_extension_threshold:        $region_extension_threshold\n";
print "genome_table_file:                 $genome_table_fname\n";
print "score_normalizer:                  $score_normalizer\n";

print "\n";

#my $optional_param_string = "";

if($background_score_path eq "NA"){
    if(-e $output_fname){
	unlink($output_fname);
    }
    system("touch $output_fname");
}

my @files_to_trash;
my $files_to_trash_counter = 0;

my @region_files;
my $region_files_counter = 0;


## read genome table

if( ! -e $genome_table_fname){
    print "Error in run peak caller script: couldn't find genome table $genome_table_fname.\n";
    print "Aborting.\n";
    exit(0);
}

open genome_table_file, "< $genome_table_fname" || die "$genome_table_fname: $\n";
my @genome_table = <genome_table_file>;

my @contigs;
my $contig_counter = 0;

for(my $i=0; $i<scalar(@genome_table); $i++){
    my $cur_genome_table_entry = $genome_table[$i];
    chomp($cur_genome_table_entry);

    my @cur_genome_table_entry_fields = split(/ /, $cur_genome_table_entry);
    if( scalar(@cur_genome_table_entry_fields) == 2){
	$contigs[$contig_counter] = $cur_genome_table_entry_fields[0];
	$contig_counter++;
    }
}

## /read genome table

for(my $i=0; $i<scalar(@contigs); $i++){
    my $cur_contig = $contigs[$i];
    my $system_command;
    my $cur_contig_regions_file = $output_fname . "." . $cur_contig . ".tmp";
    $system_command = $exec_path . 
	" score_file=$background_score_path/$cur_contig.score" . 
	" output_file=$cur_contig_regions_file" .
	" region_seeding_threshold=$region_seeding_threshold" .
	" region_extension_threshold=$region_extension_threshold" .
	" score_normalizer=$score_normalizer" .
	" contig_id=$cur_contig";

    print "\nsystem_command: $system_command\n";
    
    $files_to_trash[$files_to_trash_counter] = $cur_contig_regions_file;
    $files_to_trash_counter++;

    $region_files[$region_files_counter] = $cur_contig_regions_file;
    $region_files_counter++;
    my $error_code = 0;
    $error_code = system("$system_command");
    if($error_code != 0){
	print "run_report_bad_control_regions_with_param_file.pl: Caught an error message (code $error_code), passing the error message to the top script.\n";
	exit(2);
    }
}

# merging
print "\n";
print "Merging artifactual regions into a single file\n";

my $cat_string = "";
#
#
#my $output_regions_fname = $output_fname_prefix; # . ".regions.no_metrics";

my $merged_regions_fname = $output_fname . ".tmp1";
my $sorted_regions_fname = $output_fname;
#my $extended_tmp2_fname = $output_fname_prefix . ".tmp2";
#my $extended_tmp3_fname = $output_fname_prefix . ".tmp3";

$files_to_trash[$files_to_trash_counter] = $merged_regions_fname;
$files_to_trash_counter++;

if(-e $merged_regions_fname){
    unlink($merged_regions_fname);
}
system("touch $merged_regions_fname");

#if(-e $extended_tmp1_fname){
#    unlink($extended_tmp1_fname);
#}
#system("touch $extended_tmp1_fname");
#
#if(-e $extended_tmp2_fname){
#    unlink($extended_tmp2_fname);
#}
#system("touch $extended_tmp2_fname");
#
#if(-e $extended_tmp3_fname){
#    unlink($extended_tmp3_fname);
#}
#system("touch $extended_tmp3_fname");

for(my $i=0; $i<scalar(@region_files); $i++){
    if($i!=0){
	$cat_string = $cat_string . " ";
    }
    $cat_string = $cat_string . "$region_files[$i]";
}

system("cat $cat_string > $merged_regions_fname");
# /merging
#
# sorting regions
print "Sorting artifactual regions\n";

open sorted_regions_file, "> $sorted_regions_fname" or die "can't open $sorted_regions_fname: $!\n";
open merged_regions_file, "< $merged_regions_fname" or die "can't open $merged_regions_fname: $!\n";

my @regions;
my $region_counter = 0;

while(<merged_regions_file>){
    my $cur_region = $_;
    if(length($cur_region) > 0){
	if(substr($cur_region, 0, 1) eq "R"){
	    chomp($cur_region);
	    my $cur_region_output_string = $cur_region;
	    my @cur_region_entries = split(/ /, $cur_region);
	    my $cur_region_max_score = $cur_region_entries[6];

#	    my $end_of_region_output_found = "false";
#	    while($end_of_region_output_found eq "false"){
#		my $cur_peak = <extended_tmp1_file>;
#		if( length($cur_peak) > 0 ){
#		    if(substr($cur_peak, 0, 1) eq "P"){
#			$cur_region_output_string = $cur_region_output_string . $cur_peak;
#		    }
#		    else{
#			$end_of_region_output_found = "true";
#		    }
#		}
#		else{
#		    $end_of_region_output_found = "true";
#		}
#	    }
#
	    my $cur_region = {
		entry => $cur_region_output_string,
		score => $cur_region_max_score,
	    };

	    $regions[$region_counter] = $cur_region;
	    $region_counter++;		     
	}
    }
}

print "\nA total of $region_counter artifactual region(s) were detected from control data\n\n";

close merged_regions_file;

my @regions_sorted = reverse sort{ $a->{score} <=> $b->{score} } @regions;

for (my $i=0; $i<scalar(@regions_sorted); $i++){
    printf sorted_regions_file "%s\n", $regions_sorted[$i]->{entry};
}

close sorted_regions_file;

# /sorting regions

# adding count tags to the regions and peaks

#open extended_tmp2_file, "< $extended_tmp2_fname" or die "can't open $extended_tmp2_fname: $!\n";
#open extended_tmp3_file, "> $extended_tmp3_fname" or die "can't open $extended_tmp3_fname: $!\n";
#
#$region_counter = 1;
#while(<extended_tmp2_file>){
#    my $cur_region = $_;
#    if(length($cur_region) > 0){
#	if(substr($cur_region, 0, 1) eq "R"){
#
#	    chomp($cur_region);
#	    my @cur_region_entries = split(/ /, $cur_region);
#	    $cur_region_entries[0] = $cur_region_entries[0] . "-$region_counter";
#
#	    for(my $i=0; $i<scalar(@cur_region_entries); $i++){
#		if($i>0){ print extended_tmp3_file " "; }
#		print extended_tmp3_file "$cur_region_entries[$i]";
#	    }
#	    print extended_tmp3_file "\n";
#    
#	    my $peak_counter = 1;
#	    my $end_of_region_output_found = "false";
#	    while($end_of_region_output_found eq "false"){
#		my $cur_peak = <extended_tmp2_file>;
#		chomp($cur_peak);
#		if( length($cur_peak) > 0 ){
#		    if(substr($cur_peak, 0, 1) eq "P"){
#			my @cur_peak_entries = split(/ /, $cur_peak);
#			$cur_peak_entries[0] = $cur_peak_entries[0] . "-$region_counter-$peak_counter";
#			$peak_counter++;
#			for(my $i=0; $i<scalar(@cur_peak_entries); $i++){
#			    if($i>0){ print extended_tmp3_file " "; }
#			    print extended_tmp3_file "$cur_peak_entries[$i]";
#			}
#			print extended_tmp3_file "\n";
#		    }
#		    else{
#			$end_of_region_output_found = "true";
#		    }
#		}
#		else{
#		    $end_of_region_output_found = "true";
#		}
#	    }
#	    $region_counter++;		
#	    print extended_tmp3_file "\n";
#	}
#    }
#}
#
# /adding count tags to the regions and peaks

#close extended_tmp2_file;
#close extended_tmp3_file;

# /sorting regions
# output
#
#open output_regions_file, "> $output_regions_fname" or die "can't open $output_regions_fname: $!\n";
#
#open extended_tmp3_file, "< $extended_tmp3_fname" or die "can't open $extended_tmp3_fname: $!\n";
#
#my @peaks;
#my $peak_count = 0;
#
#while(<extended_tmp3_file>){
#    chomp;
#
#    my $cur_entry = $_;
#    if(length($cur_entry) > 0){
#	if(substr($cur_entry, 0, 1) eq "R" || substr($cur_entry, 0, 1) eq "P"){
#	    my @cur_entry_fields = split(/ /, $cur_entry);
#	    if(scalar(@cur_entry_fields) >= 7){
#		my $cur_ChIP_score = $cur_entry_fields[4];
#		my $cur_control_score = $cur_entry_fields[6];
#		
#		my $cur_ef = "inf";
#		if($cur_control_score != 0 && $ChIP_reads > 0 && $background_reads > 0){
#		    $cur_ef = ($cur_ChIP_score/$ChIP_reads) / ($cur_control_score/$background_reads);
#
#		    printf output_regions_file "%s ef: %.3f\n", $cur_entry, $cur_ef;
#		}
#		else{
#		    $cur_ef = $cur_ChIP_score / $ChIP_basal_level;
#		}
#		printf output_regions_file "%s ef: %.3f\n", $cur_entry, $cur_ef;
#	    }
#
#	    if(substr($cur_entry, 0, 1) eq "P"){		
#		my @cur_entry_fields = split(/ /, $cur_entry);
#		if(scalar(@cur_entry_fields) >= 7){
#		    my $cur_ChIP_score = $cur_entry_fields[4];
#		    my $cur_control_score = $cur_entry_fields[6];
#		    
#		    my $cur_peak = {
#			entry => $cur_entry,
#			ChIP_score => $cur_ChIP_score,
#			control_score => $cur_control_score,
#		    };
#		    $peaks[$peak_count] = $cur_peak;
#		    $peak_count++;		
#		}
#		else{
#		    print "Error: couldn't understand the peak entry:\n$cur_entry\n";
#		    print "Skipping\n\n";
#		}
#	    }
#	}
#	else{ print output_regions_file "$cur_entry\n"; }
#    }
#    else{ print output_regions_file "$cur_entry\n"; }
#}
#
#close output_regions_file;#
#
#system("rm $extended_tmp1_fname");
#system("rm $extended_tmp2_fname");
#system("rm $extended_tmp3_fname");

print "Removing temporarary files.\n";

for(my $i=0; $i<scalar(@files_to_trash); $i++){
    unlink($files_to_trash[$i]);
}

print "done\n";
