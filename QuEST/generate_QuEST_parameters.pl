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

## This program generates parameters for individual QuEST modules

use strict;
use warnings;
use diagnostics;

use Cwd qw(realpath);

my $usage = qq{
    generate_QuEST_parameters.pl

    --------------------------------
    Citation:

    Valouev A, Johnson DS, Sundquist A, Medina C, Anton E, Batzoglou S,
    Myers RM, Sidow A (2008). Genome-wide analysis of transcription factor binding sites based on ChIP-Seq data.
    Nat Methods 5, 829-834 (2008)

    -------------------------------

    QuEST parametrization script

    -----------------------
    mandatory arguments:
    -----------------------

    for the mapped ChIP data use one of the following:

    -QuEST_align_ChIP <QuEST_align_file>
    -solexa_align_ChIP <solexa_ChIP_align_file>      
    -eland_align_ChIP <eland_file>
    -eland_extended_ChIP <eland_extended_file>
    -bowtie_align_ChIP <bowtie_file>
    -maq_align_ChIP <maq_file>
    -sam_align_ChIP <sam_file>

    - - - - - - - - - - - - 

    for the mapped background data (RX_noIP, input, total chromatin, control) use one of the following:

    -QuEST_align_RX_noIP <QuEST_align_RX_noIP_file>
    -solexa_align_RX_noIP  <solexa_RX_noIP_align_file>   
    -eland_align_RX_noIP <eland_file>
    -eland_extended_RX_noIP <eland_extended_file>
    -bowtie_align_RX_noIP <bowtie_file>
    -maq_align_RX_noIP <maq_file>
    -sam_align_RX_noIP <sam_file>

    - - - - - - - - - - - 

    for the reference genome use one of the following:

    -rp <reference_path>              a path to the reference genome fasta files
                                      formated in UCSC Genome Browser database fashion
                                      (typical naming convention is chr1.fa chr2.fa ...
				      chrM.fa)

    -gt <genome_table>                a genome table containing contig names followed by
                                      their sizes				      

    - - - - - - - - - - - 				      

    for the analysis:				      
				      
    -ap <analysis_path>               specify a directory where all QuEST analysis results
                                      will be reported

    -----------------------
    optional arguments:
    -----------------------

    -silent                          assume all default parameters and do not alert
                                     to override any values (good for automation)
    -advanced                        configure all of QuEST parameters (-silent flag
                                     should be off)
    -ChIP_name <TF name>             will appear in the names of browser tracks and some files

};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

use Cwd qw(realpath);
my $QuEST_fullpath = realpath($0);
my @QuEST_fullpath_fields = split(/\//,$QuEST_fullpath);
my $QuEST_exec_path = "";

for(my $i=0; $i<scalar(@QuEST_fullpath_fields)-1; $i++){
    if($QuEST_fullpath_fields[$i] ne ""){
        $QuEST_exec_path = $QuEST_exec_path."/".$QuEST_fullpath_fields[$i];
    }
}

my $QuEST_info_fname = $QuEST_exec_path . "/QuEST.info";
open QuEST_info_file, "< $QuEST_info_fname" || die "$QuEST_info_fname: $!\n";
while(<QuEST_info_file>){
    chomp;
    print "$_\n";
}
close QuEST_info_file;

my $generate_profile_exec = $QuEST_exec_path . "/generate_profile";
my $peak_caller_exec = $QuEST_exec_path . "/peak_caller";
my $quick_window_scan_exec = $QuEST_exec_path . "/quick_window_scan";
my $calibrate_peak_shift_exec = $QuEST_exec_path . "/calibrate_peak_shift";
my $generate_profile_views_exec = $QuEST_exec_path . "/generate_profile_views";
my $metrics_exec = $QuEST_exec_path . "/metrics";
my $profile_2_wig_exec = $QuEST_exec_path . "/profile_2_wig";
my $align_2_QuEST_exec = $QuEST_exec_path . "/align_2_QuEST";
my $QuEST_configure_exec = $QuEST_exec_path . "/configure.pl";
my $align_2_bin_exec = $QuEST_exec_path . "/collapse_reads";
my $report_bad_control_regions_exec = $QuEST_exec_path . "/report_bad_control_regions";

my $run_QuEST_script = $QuEST_exec_path . "/run_QuEST_with_param_file.pl";

my $missing_execs = "";

if(! -e $generate_profile_exec){ $missing_execs = $missing_execs . " " . $generate_profile_exec; }
if(! -e $peak_caller_exec){ $missing_execs = $missing_execs . " " . $peak_caller_exec; }
if(! -e $quick_window_scan_exec){ $missing_execs = $missing_execs . " " . $quick_window_scan_exec; }
if(! -e $align_2_bin_exec){ $missing_execs = $missing_execs . " " . $align_2_bin_exec; }
if(! -e $calibrate_peak_shift_exec){ $missing_execs = $missing_execs . " " . $calibrate_peak_shift_exec; }
if(! -e $generate_profile_views_exec){ $missing_execs = $missing_execs . " " . $generate_profile_views_exec; }
if(! -e $metrics_exec){ $missing_execs = $missing_execs . " " . $metrics_exec; }
if(! -e $profile_2_wig_exec){ $missing_execs = $missing_execs . " " . $profile_2_wig_exec; }
if(! -e $align_2_QuEST_exec){ $missing_execs = $missing_execs . " " . $align_2_QuEST_exec; }
if(! -e $report_bad_control_regions_exec){ $missing_execs = $missing_execs . " " . $report_bad_control_regions_exec; }

if( $missing_execs ne ""){
    print "The following QuEST executables are missing: \n";
    print "$missing_execs\n";
    print "It is likely that you did not compile QuEST.\n";
    print "Let's perform QuEST compilation.";

    system("$QuEST_configure_exec");
    system("cd $QuEST_exec_path; make");
}

$missing_execs = "";

if(! -e $generate_profile_exec){ $missing_execs = $missing_execs . " " . $generate_profile_exec; }
if(! -e $peak_caller_exec){ $missing_execs = $missing_execs . " " . $peak_caller_exec; }
if(! -e $quick_window_scan_exec){ $missing_execs = $missing_execs . " " . $quick_window_scan_exec; }
if(! -e $align_2_bin_exec){ $missing_execs = $missing_execs . " " . $align_2_bin_exec; }
if(! -e $calibrate_peak_shift_exec){ $missing_execs = $missing_execs . " " . $calibrate_peak_shift_exec; }
if(! -e $generate_profile_views_exec){ $missing_execs = $missing_execs . " " . $generate_profile_views_exec; }
if(! -e $metrics_exec){ $missing_execs = $missing_execs . " " . $metrics_exec; }
if(! -e $profile_2_wig_exec){ $missing_execs = $missing_execs . " " . $profile_2_wig_exec; }
if(! -e $align_2_QuEST_exec){ $missing_execs = $missing_execs . " " . $align_2_QuEST_exec; }
if(! -e $report_bad_control_regions_exec){ $missing_execs = $missing_execs . " " . $report_bad_control_regions_exec; }

if( $missing_execs ne ""){
    print "The following QuEST executables are still missing: \n";
    print "$missing_execs\n";
    print "QuEST has likely failed to compile.\n";
    print "Make sure that g++ and make are installed on your system\n";
    exit(4);
}


my $undef = "undef";

## mandatory argmuments
my $missing = "** missing **";

my $ChIP_align_fname =        $missing;
my $ChIP_file_type =          "";

my $RX_noIP_align_fname =     $missing;
my $RX_noIP_file_type =       "";

my $reference_path =                 $missing;
my $genome_table_source_fname =      $missing;

my $analysis_path =                  $missing;

## /mandatory arguments

## optional arguments
my $ChIP_name = "ChIP";
## /optional arguments

## flag options

#my $solexa_align_ChIP_flag =          "off";
#my $QuEST_align_ChIP_flag =           "off";
#
#my $solexa_align_RX_noIP_flag =       "off";
#my $QuEST_align_RX_noIP_flag =        "off";

my $rp_flag =                         "off";
my $gt_flag =                         "off";
my $ap_flag =                         "off";

my $silent_flag =                     "off";
my $advanced_flag =                   "off";

my $alternative_ps = "undef";

## reading arguments
while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;    
    if ( $this_arg eq '-h') {print "$usage\n"; exit(4)}
    elsif ( $this_arg eq '-solexa_align_ChIP') {$ChIP_align_fname = shift @ARGV; $ChIP_file_type = "solexa_align";}
    elsif ( $this_arg eq '-QuEST_align_ChIP') {$ChIP_align_fname = shift @ARGV; $ChIP_file_type = "QuEST_align";}
    elsif ( $this_arg eq '-eland_align_ChIP') {$ChIP_align_fname = shift @ARGV; $ChIP_file_type = "eland_align";}
    elsif ( $this_arg eq '-eland_extended_ChIP') {$ChIP_align_fname = shift @ARGV; $ChIP_file_type = "eland_extended_align";}
    elsif ( $this_arg eq '-bowtie_align_ChIP') {$ChIP_align_fname = shift @ARGV; $ChIP_file_type = "bowtie_align";}
    elsif ( $this_arg eq '-maq_align_ChIP') {$ChIP_align_fname = shift @ARGV; $ChIP_file_type = "maq_align";}
    elsif ( $this_arg eq '-sam_align_ChIP') {$ChIP_align_fname = shift @ARGV; $ChIP_file_type = "sam_align";}

    elsif ( $this_arg eq '-solexa_align_RX_noIP') {$RX_noIP_align_fname = shift @ARGV; $RX_noIP_file_type = "solexa_align";}
    elsif ( $this_arg eq '-QuEST_align_RX_noIP') {$RX_noIP_align_fname = shift @ARGV; $RX_noIP_file_type = "QuEST_align";}
    elsif ( $this_arg eq '-eland_align_RX_noIP') {$RX_noIP_align_fname = shift @ARGV; $RX_noIP_file_type = "eland_align";}
    elsif ( $this_arg eq '-eland_extended_RX_noIP') {$RX_noIP_align_fname = shift @ARGV; $RX_noIP_file_type = "eland_extended_align";}
    elsif ( $this_arg eq '-bowtie_align_RX_noIP') {$RX_noIP_align_fname = shift @ARGV; $RX_noIP_file_type = "bowtie_align";}
    elsif ( $this_arg eq '-maq_align_RX_noIP') {$RX_noIP_align_fname = shift @ARGV; $RX_noIP_file_type = "maq_align";}
    elsif ( $this_arg eq '-sam_align_RX_noIP') {$RX_noIP_align_fname = shift @ARGV; $RX_noIP_file_type = "sam_align";}

    elsif ( $this_arg eq '-rp') {$reference_path = shift @ARGV; $rp_flag = "on";}
    elsif ( $this_arg eq '-gt') { $genome_table_source_fname = shift @ARGV; $gt_flag = "on";}
    elsif ( $this_arg eq '-ap') {$analysis_path = shift @ARGV; $ap_flag = "on";}
    elsif ( $this_arg eq '-silent') {$silent_flag = "on";}
    elsif ( $this_arg eq '-advanced') {$advanced_flag = "on";}
    elsif ( $this_arg eq '-ChIP_name') {$ChIP_name = shift @ARGV;}
    elsif ( $this_arg eq '-alternative_ps') {$alternative_ps = shift @ARGV;}
    else{ print "Warning: unknown parameter: $this_arg\n"; exit(4); }
}

my $die_pars_bad = "false";       ## die if parameters are bad
my $bad_par = "";                 ## store bad parameter here

my $cur_path = "";

my $fullpath = realpath(".");

my @fullpath_fields = split(/\//,$fullpath);
for(my $i=0; $i<scalar(@fullpath_fields); $i++){
    if($fullpath_fields[$i] ne ""){
        $cur_path = $cur_path."/".$fullpath_fields[$i];
    }
}

if( substr($analysis_path,0,2) eq "./" ){
    $analysis_path = $cur_path . "/" . substr($analysis_path,2,length($analysis_path)-2);
    while(substr($analysis_path,length($analysis_path)-1) eq "/"){
	$analysis_path = substr($analysis_path,0,length($analysis_path)-1);
    }
}

print "\n=====================\n\n";
print "parameters: \n\n";
print "ChIP_name:                     $ChIP_name\n";
print "ChIP_align_file:               $ChIP_align_fname\n";
print "ChIP_file_type:                $ChIP_file_type\n";
print "\n";
print "RX_noIP_align_file:            $RX_noIP_align_fname\n";
print "RX_noIP_file_type:             $RX_noIP_file_type\n";
print "\n";

print "reference path:                $reference_path\n";
print "genome_table:                  $genome_table_source_fname\n";

print "\n";

print "analysis_path:                 $analysis_path\n";

print "\n";

print "silent:                        $silent_flag\n";
print "advanced:                      $advanced_flag\n";


if( $analysis_path eq $missing ){ $die_pars_bad = "true"; $bad_par = $bad_par." "."analysis_path"; }

print "\n=====================\n\n";

if( $die_pars_bad eq "true"){
    print "$usage";
    print "Bad parameters: $bad_par\n";
    exit(0);
}
## /reading arguments

if( $reference_path eq $missing and $genome_table_source_fname eq $missing){
    print "You should provide either the path to the reference genome or \n";
    print "the genome table containing contig names and their sizes. Aborting.\n";
    exit(0);
}

if( $ChIP_file_type eq ""){
    print "You must supply ChIP align file. Aborting.\n";
    exit(0);
}

if( $RX_noIP_align_fname eq $missing){
    print "************************************************************\n";
    print "**    Warning: you have not provided any control data     **\n";
    print "**       Analysis results without control data are        **\n";
    print "**     in some cases very unreliable. If possible, use    **\n";
    print "**    use control data to obtain more accurate results.   **\n";
    print "**                                                        **\n";
    print "**               Press any button to continue.            **\n";
    print "************************************************************\n";

    my $kbhit = <STDIN>;
}


my $ChIP_tracks_base_priority = int((1+rand(9))*10);

while ($ChIP_tracks_base_priority <= 2){
    $ChIP_tracks_base_priority = int((1+rand(9))*10);
}

my $ChIP_bedgraph_track_priority = $ChIP_tracks_base_priority + 0.1;
my $ChIP_reads_track_priority = $ChIP_tracks_base_priority + 0.5;
my $ChIP_normalized_wig_track_priority = $ChIP_tracks_base_priority + 1;
my $ChIP_unnormalized_wig_track_priority = $ChIP_tracks_base_priority + 2;
my $ChIP_calls_track_priority = $ChIP_tracks_base_priority + 3;


my $RX_noIP_tracks_base_priority =  (1+rand(9))*10;

while($RX_noIP_tracks_base_priority > $ChIP_tracks_base_priority){
    $RX_noIP_tracks_base_priority =  (1+rand(9))*10;
}
my $RX_noIP_reads_track_priority = $RX_noIP_tracks_base_priority + 0;
my $RX_noIP_normalized_wig_track_priority = $RX_noIP_tracks_base_priority + 1;
my $RX_noIP_unnormalized_wig_track_priority = $RX_noIP_tracks_base_priority + 2;
my $RX_noIP_calls_track_priority = $RX_noIP_tracks_base_priority + 3;


## /setting up analysis directory and subdrirectories
print "creating necessary directories.\n\n";
print "- - - - - - - - - - - - - - - -\n\n";

my $logs_path = $analysis_path . "/logs";
my $module_outputs_path = $analysis_path . "/module_outputs";
my $parameters_path = $analysis_path . "/parameters";
my $calls_path = $analysis_path . "/calls";
my $data_path = $analysis_path . "/data";
my $scores_path = $analysis_path . "/scores";
my $tracks_path = $analysis_path . "/tracks";

my $bin_align_path = $analysis_path . "/bin_align";
my $bin_align_ChIP_path = $bin_align_path . "/ChIP";
my $bin_align_background_path = $bin_align_path . "/background";
my $bin_align_RX_noIP_path = $bin_align_path . "/RX_noIP";
my $bin_align_pseudo_ChIP_path = $bin_align_path . "/pseudo_ChIP";

my $bed_graph_tracks_path = $tracks_path . "/bed_graph";
my $bed_graph_tracks_by_chr_path = $bed_graph_tracks_path . "/by_chr";
my $bed_graph_tracks_by_chr_ChIP_path = $bed_graph_tracks_by_chr_path . "/ChIP";

my $data_bed_tracks_path = $tracks_path . "/data_bed_files";
my $data_bed_tracks_by_chr_path = $data_bed_tracks_path . "/by_chr";
my $ChIP_bed_tracks_by_chr_path = $data_bed_tracks_by_chr_path . "/ChIP";
my $RX_noIP_bed_tracks_by_chr_path = $data_bed_tracks_by_chr_path . "/RX_noIP";

my $wig_tracks_path = $tracks_path . "/wig_profiles";
my $wig_tracks_by_chr_path = $wig_tracks_path . "/by_chr";
my $wig_tracks_ChIP_normalized_by_chr_path = $wig_tracks_by_chr_path . "/ChIP_normalized";
my $wig_tracks_ChIP_unnormalized_by_chr_path = $wig_tracks_by_chr_path . "/ChIP_unnormalized";
my $wig_tracks_background_normalized_by_chr_path = $wig_tracks_by_chr_path . "/background_normalized";
my $wig_tracks_background_unnormalized_by_chr_path = $wig_tracks_by_chr_path . "/background_unnormalized";

my $ChIP_scores_path = $scores_path . "/ChIP";
my $pseudo_ChIP_scores_path = $scores_path . "/pseudo_ChIP";
my $background_scores_path = $scores_path . "/background";

my $cur_log_fname = $logs_path . "/generate_QuEST_pars.log";
open cur_log_file, "< $cur_log_fname" || die "$cur_log_fname: $!\n";


print "analysis directory: $analysis_path:";
if(-d $analysis_path){
    print " already exists\n";

    # test write
    my $dummy_fname = $analysis_path . "/dummy_file.tmp";
    system("touch $dummy_fname");
    if( -e $dummy_fname){
    }
    else{
	print "\n";
	print "Unable to create files in the analysis directory. Check your permissions.\n";
	exit(1);
    }
    unlink($dummy_fname);
    # /test write
}
else{
    mkdir $analysis_path or die "Failed to create the analysis directory\n";
    print "   created\n";
}

print "logs directory: $logs_path: ";
if(-d $logs_path){
    print "already exists\n";
}
else{
    print "created\n";
    mkdir($logs_path);
}

print "module outputs directory: $module_outputs_path: ";
if(-d $module_outputs_path){
    print "already exists\n";
}
else{
    print "created\n";
    mkdir($module_outputs_path);
}

print "parameters directory: $parameters_path: ";
if(-d $parameters_path){
    print "already exists\n";
}
else{
    print "created\n";
    mkdir($parameters_path);
}


print "calls directory: $calls_path: ";
if(-d $calls_path){
    print "already exists\n";
}
else{
    print "created\n";
    mkdir($calls_path);
}

print "data directory: $data_path: ";
if(-d $data_path){
    print "already exists\n";
}
else{
    print "created\n";
    mkdir($data_path);
}

print "scores directory: $scores_path: ";

if(-d $scores_path){
    print "already exists\n";
}
else{
    print "created\n";
    mkdir($scores_path);
}


print "tracks directory: $tracks_path: ";
if(-d $tracks_path){
    print "already exists\n";
}
else{
    mkdir($tracks_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "bin_align directory: $bin_align_path: ";
if(-d $bin_align_path){
    print "already exists\n";
}
else{
    mkdir($bin_align_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "bin_align ChIP directory: $bin_align_ChIP_path: ";
if(-d $bin_align_ChIP_path){
    print "already exists\n";
}
else{
    mkdir($bin_align_ChIP_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "bin_align background directory: $bin_align_background_path: ";
if(-d $bin_align_background_path){
    print "already exists\n";
}
else{
    mkdir($bin_align_background_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "bin_align RX_noIP directory: $bin_align_RX_noIP_path: ";
if(-d $bin_align_RX_noIP_path){
    print "already exists\n";
}
else{
    mkdir($bin_align_RX_noIP_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "bin_align pseudo_ChIP directory: $bin_align_pseudo_ChIP_path: ";
if(-d $bin_align_pseudo_ChIP_path){
    print "already exists\n";
}
else{
    mkdir($bin_align_pseudo_ChIP_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "bed graph tracks directory: $bed_graph_tracks_path: ";
if(-d $bed_graph_tracks_path){
    print "already exists\n";
}
else{
    mkdir($bed_graph_tracks_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "bed graph tracks by chromosome directory: $bed_graph_tracks_by_chr_path: ";
if(-d $bed_graph_tracks_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($bed_graph_tracks_by_chr_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "bed graph tracks by chromosome ChIP directory: $bed_graph_tracks_by_chr_ChIP_path: ";
if(-d $bed_graph_tracks_by_chr_ChIP_path){
    print "already exists\n";
}
else{
    mkdir($bed_graph_tracks_by_chr_ChIP_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "data bed tracks directory: $data_bed_tracks_path: ";
if(-d $data_bed_tracks_path){
    print "already exists\n";
}
else{
    mkdir($data_bed_tracks_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "data bed tracks by chromosome directory: $data_bed_tracks_by_chr_path: ";
if(-d $data_bed_tracks_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($data_bed_tracks_by_chr_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "ChIP bed tracks by chromosome directory: $ChIP_bed_tracks_by_chr_path: ";
if(-d $ChIP_bed_tracks_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($ChIP_bed_tracks_by_chr_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "RX_noIP bed tracks by chromosome directory: $RX_noIP_bed_tracks_by_chr_path: ";
if(-d $RX_noIP_bed_tracks_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($RX_noIP_bed_tracks_by_chr_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "wig tracks directory: $wig_tracks_path: ";
if(-d $wig_tracks_path){
    print "already exists\n";
}
else{
    mkdir($wig_tracks_path) or die "Can't create directory: $!\n";
    print "created\n";
}

print "wig tracks by chromosome directory: $wig_tracks_by_chr_path: ";
if(-d $wig_tracks_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($wig_tracks_by_chr_path) or die "Can't create directory $wig_tracks_by_chr_path: $!\n";
    print "created\n";
}

print "wig tracks by chromosome directory for normalized ChIP profiles: $wig_tracks_ChIP_normalized_by_chr_path: ";
if(-d $wig_tracks_ChIP_normalized_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($wig_tracks_ChIP_normalized_by_chr_path) or die "Can't create directory $wig_tracks_ChIP_normalized_by_chr_path: $!\n";
    print "created\n";
}

print "wig tracks by chromosome directory for unnormalized ChIP profiles: $wig_tracks_ChIP_unnormalized_by_chr_path: ";
if(-d $wig_tracks_ChIP_unnormalized_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($wig_tracks_ChIP_unnormalized_by_chr_path) or die "Can't create directory $wig_tracks_ChIP_unnormalized_by_chr_path: $!\n";
    print "created\n";
}

print "wig tracks by chromosome directory for normalized backround profiles: $wig_tracks_background_normalized_by_chr_path: ";
if(-d $wig_tracks_background_normalized_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($wig_tracks_background_normalized_by_chr_path) or die "Can't create directory $wig_tracks_background_normalized_by_chr_path: $!\n";
    print "created\n";
}

print "wig tracks by chromosome directory for unnormalized background profiles: $wig_tracks_background_unnormalized_by_chr_path: ";
if(-d $wig_tracks_background_unnormalized_by_chr_path){
    print "already exists\n";
}
else{
    mkdir($wig_tracks_background_unnormalized_by_chr_path) or die "Can't create directory $wig_tracks_background_unnormalized_by_chr_path: $!\n";
    print "created\n";
}


print "ChIP scores directory: $ChIP_scores_path: ";
if(-d $ChIP_scores_path){
    print "already exists\n";
}
else{
    mkdir($ChIP_scores_path)  or die "Can't create directory: $!\n";
    print "created\n";
}

print "background scores directory: $background_scores_path: ";
if(-d $background_scores_path){
    print "already exists\n";
}
else{
    mkdir($background_scores_path)  or die "Can't create directory: $!\n";
    print "created\n";
}

print "pseudo_ChIP scores directory: $pseudo_ChIP_scores_path: ";
if(-d $pseudo_ChIP_scores_path){
    print "already exists\n";
}
else{
    mkdir($pseudo_ChIP_scores_path)  or die "Can't create directory: $!\n";
    print "created\n";
}



print "\n";
print "- - - - - - - - - - - - - - - -\n\n";

## /setting up analysis directory and subdrirectories

## reading reference files

my @contig_names;
my @contig_sizes;
my $contig_counts = 0;

if( $reference_path ne "** missing **" ){
    opendir(REF_DIR, $reference_path) or die "can't open $reference_path: $!";
    
    my @fasta_fnames;
    my $fasta_fnames_counter = 0;
    while (defined(my $cur_fname = readdir(REF_DIR))){
	if ( substr($cur_fname,length($cur_fname)-3,3) eq ".fa" ){
	    $fasta_fnames[$fasta_fnames_counter] = $cur_fname; 
	    $fasta_fnames_counter++;
	}
    }
    close(REF_DIR);
    print "\n";
    print "The following fasta files were found in your reference path:\n\n";
    for(my $i=0; $i<scalar(@fasta_fnames); $i++){
	
	print "$fasta_fnames[$i]\n";
    }
    print "\nEvaluating the size of the genome. This may take several minutes.\n";
    
    print "\nFound the following contigs: \n\n";
    
    
    for(my $i=0; $i<scalar(@fasta_fnames); $i++){
	my $cur_contig_fname = $reference_path . "/$fasta_fnames[$i]";
	my $cur_contig_name = "";
	my $cur_contig_size = 0;
	
	if(-e $cur_contig_fname){
	    open cur_contig_file, "< $cur_contig_fname";
	    while(<cur_contig_file>){
		chomp;
		if( substr($_,0,1) eq ">" ){
		    if($cur_contig_name eq ""){
			# first contig, don't save anything
		    }
		    else{			
			$contig_names[$contig_counts] = $cur_contig_name;
			$contig_sizes[$contig_counts] = $cur_contig_size;
			$contig_counts++;
			
			print "$cur_contig_name: $cur_contig_size bps\n";
		    }
		    $cur_contig_name = substr($_,1,length($_)-1);
		}
		elsif( substr($_,0,1) eq "#" ){
		    # skip, comment
		}
		else{
		    $cur_contig_size = $cur_contig_size + length($_);
		    #increment the size
		}
	    }
	    if($cur_contig_name ne ""){
		$contig_names[$contig_counts] = $cur_contig_name;
		$contig_sizes[$contig_counts] = $cur_contig_size;
		$contig_counts++;
		
		print "$cur_contig_name: $cur_contig_size bps\n";
	    }
	    close cur_contig_file;
	}    
    }
}
else{
    print "Found the following contigs in the genome table: \n\n";
    if( ! -e $genome_table_source_fname ){
	print "Failed to locate genome table: $genome_table_source_fname\n";
	print "Aborting.\n";
	exit(4)
    }
    
    open genome_table_source_file, "< $genome_table_source_fname" || die "$genome_table_source_fname: $!\n";

    my @gt = <genome_table_source_file>;

    for(my $j=0; $j < scalar(@gt); $j++){
	my $cur_gt_entry = $gt[$j];
	chomp($cur_gt_entry);
	if(substr($cur_gt_entry,0,1) ne "#"){
	    my @cur_gt_entry_fields = split(/ /, $cur_gt_entry);
	    if( scalar(@cur_gt_entry_fields) == 2){
		$contig_names[$contig_counts] = $cur_gt_entry_fields[0];
		$contig_sizes[$contig_counts] = $cur_gt_entry_fields[1];
		$contig_counts++;

		print "$cur_gt_entry_fields[0] $cur_gt_entry_fields[1]\n";
	    }
	}
    }

    close genome_table_source_file;
}

if (scalar(@contig_names) == 0){
    print "There were not contigs found. Please correct reference path or genome table.\n";
    exit(1);
}

my $duplicate_contigs = "false";
for(my $i=0; $i<scalar(@contig_names); $i++){
    my $cur_contig_name = $contig_names[$i];
    for(my $j=0; $j<$i; $j++){
	if($cur_contig_name eq $contig_names[$j]){
	    $duplicate_contigs = "true";
	}
    }
}

if($duplicate_contigs eq "true"){
    print "It appears that there are duplicate contigs in the genome that you provided.";
    print "Please revise the genome or genome table\n";
    exit(1);
}
else{
    print "\nYour genome looks fine and does not appear to contain any\n";
    print "duplicated chromosome entries.\n";
}


my $contigs_string = "";

my $genome_size = 0;

for(my $i=0; $i<scalar(@contig_names); $i++){
    if($i!=0){ $contigs_string = $contigs_string . ","; }
    $contigs_string = $contigs_string . "$contig_names[$i]";
    $genome_size = $genome_size + $contig_sizes[$i];
}

printf "\nEstimated size of the genome is %.1f Kb ( %.2f Gbps ).\n", $genome_size/1000, $genome_size/1000000000;

my $contig_table_fname = $parameters_path . "/genome_table";

if( -e $contig_table_fname ){
    unlink($contig_table_fname);
}

open contig_table_file, "> $contig_table_fname";

for(my $i=0; $i<scalar(@contig_sizes); $i++){
    print contig_table_file "$contig_names[$i] $contig_sizes[$i]\n";
}

close contig_table_file;

my $ChIP_uncollapsed_align_fname = $data_path . "/ChIP_uncollapsed.QuEST";

my $answer = "";

if( ! -e $ChIP_align_fname ){
    print "ChIP_align_file: $ChIP_align_fname does not exsist. Aborting.\n";
    exit(1);
}

my $ChIP_tags = 0;

print "\nConverting ChIP reads file into the QuEST format.\n";
print "If the counter below stays at zero, the program can't\n";
print "match the ChIP alignments to the reference genome that you provided.\n";
print "In this case check that the sequences you have provided match the formats\n";
print "accepted by QuEST (see README.txt) and that names of contigs in the alignment\n";
print "files are *exactly* the same as in your genome or genome table.\n\n";


my $cur_system_command;
my $ChIP_uncollapsed_read_num_fname = $data_path . "/ChIP_read_num_uncollapsed.txt";
my $ChIP_align_2_QuEST_conversion_report_fname = $module_outputs_path . "/align_2_QuEST_ChIP.report.txt";

if(-e $ChIP_uncollapsed_read_num_fname){
    unlink($ChIP_uncollapsed_read_num_fname);
}

if($ChIP_file_type eq "solexa_align"){
    $cur_system_command = "$align_2_QuEST_exec align_file=$ChIP_align_fname output_file=$ChIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=solexa_align report_file=$ChIP_uncollapsed_read_num_fname error_report_file=$ChIP_align_2_QuEST_conversion_report_fname";
}
elsif($ChIP_file_type eq "QuEST_align"){
        $cur_system_command = "$align_2_QuEST_exec align_file=$ChIP_align_fname output_file=$ChIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=QuEST report_file=$ChIP_uncollapsed_read_num_fname error_report_file=$ChIP_align_2_QuEST_conversion_report_fname";
}
elsif($ChIP_file_type eq "eland_align"){
        $cur_system_command = "$align_2_QuEST_exec align_file=$ChIP_align_fname output_file=$ChIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=eland report_file=$ChIP_uncollapsed_read_num_fname error_report_file=$ChIP_align_2_QuEST_conversion_report_fname";
}
elsif($ChIP_file_type eq "eland_extended_align"){
        $cur_system_command = "$align_2_QuEST_exec align_file=$ChIP_align_fname output_file=$ChIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=eland_extended report_file=$ChIP_uncollapsed_read_num_fname error_report_file=$ChIP_align_2_QuEST_conversion_report_fname";
}
elsif($ChIP_file_type eq "bowtie_align"){
        $cur_system_command = "$align_2_QuEST_exec align_file=$ChIP_align_fname output_file=$ChIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=bowtie report_file=$ChIP_uncollapsed_read_num_fname error_report_file=$ChIP_align_2_QuEST_conversion_report_fname";
}
elsif($ChIP_file_type eq "maq_align"){
        $cur_system_command = "$align_2_QuEST_exec align_file=$ChIP_align_fname output_file=$ChIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=maq report_file=$ChIP_uncollapsed_read_num_fname error_report_file=$ChIP_align_2_QuEST_conversion_report_fname";
}
elsif($ChIP_file_type eq "sam_align"){
        $cur_system_command = "$align_2_QuEST_exec align_file=$ChIP_align_fname output_file=$ChIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=sam report_file=$ChIP_uncollapsed_read_num_fname error_report_file=$ChIP_align_2_QuEST_conversion_report_fname";
}
else{
    print "Failed to understand the aligment type\n";
    exit;
}

my $cur_system_command_out = system($cur_system_command);

if($cur_system_command_out != 0){
    print "\n\n";
    print "Error during the execution of $cur_system_command. Error code: $cur_system_command_out.\n\n";
    exit(1);
}

open ChIP_uncollapsed_read_num_file, "< $ChIP_uncollapsed_read_num_fname" || die "Failed to open $ChIP_uncollapsed_read_num_fname\n";
my @ChIP_uncollapsed_read_num_file_contents = <ChIP_uncollapsed_read_num_file>;
close(ChIP_uncollapsed_read_num_file);

my $ChIP_uncollapsed_tags = $ChIP_uncollapsed_read_num_file_contents[0];
chomp($ChIP_uncollapsed_tags);

my $recommended_stack_p_value_threshold = 0.001;
my $recommended_collapse_window = 100;
my $recommended_percent_positions_hit_threshold = 30;
my $recommended_stack_count_threshold = 2;
my $new_stack_size = 1;

my $recommended_collapse_tags = "true";

my $stack_p_value_threshold = $recommended_stack_p_value_threshold;
my $collapse_window = $recommended_collapse_window;
my $percent_positions_hit_threshold = $recommended_percent_positions_hit_threshold;
my $stack_count_threshold = $recommended_stack_count_threshold;

my $collapse_tags = $recommended_collapse_tags;


if($advanced_flag eq "on"){

    print "QuEST recommends collapsing unusually high stacks of sequence reads because\n";
    print "these typically result from either repeated sampling of the same template DNA\n";
    print "or sequencing repeat regions.\n";

    print "\n";

    print "Collapsing stacks usually makes it easier to detect correct regions and reduce\n";
    print "the number of false positive peaks and regions.\n";

    print "\n";

    my $cur_answer = "";
    while( $cur_answer ne "y" && $cur_answer ne "n" ){
	print "Would you like QuEST to collapse stacks? (y/n, recommended y): ";
	
	if($silent_flag eq "on"){
	    $cur_answer = "y";
	    print "$cur_answer\n";
	}
	else{
	    $cur_answer = <STDIN>;
	}
	chomp($cur_answer);
    }

    if($cur_answer eq "y"){
	$collapse_tags = "true";
    }
    else{
	$collapse_tags = "false";
    }

    if($collapse_tags eq "true"){
	print "\n";
	print "When collapsing tags, a window is used to calculate adaptive\n";
	print "tag frequency. We recommend using the default $recommended_collapse_window bp window.\n";
	
	my $cur_answer1 = "";
	while($cur_answer1 ne "y" && $cur_answer1 ne "n"){
	    print "Would you like to override the size of collapse window? (y/n, recommended n): ";
	    if($silent_flag eq "on"){
		$cur_answer1 = "n";
		print "$cur_answer1\n";
	    }
	    else{
		$cur_answer1 = <STDIN>;
	    }
	    chomp($cur_answer1);
	}

	if($cur_answer1 eq "y"){
	    my $cur_answer2 = "";
	    
	    print "Please specify the width of collapse window (recommended $recommended_collapse_window bps): ";
	    $cur_answer2 = <STDIN>;
	    
	    chomp($cur_answer2);
	    $collapse_window = $cur_answer2;
	}
	
	print "\n";
	print "Stacks are collapsed to a single instance of a read if stack p-value is smaller than stack_p_value_threshold.\n";
	print "QuEST uses stack_p_value_threshold = $recommended_stack_p_value_threshold for collapsing stacks.\n";
	
	$cur_answer1 = "";
	
	while($cur_answer1 ne "y" && $cur_answer1 ne "n"){
	    print "Would you like to override this p-value? (y/n, recommended n): ";
	    if($silent_flag eq "on"){
		$cur_answer1 = "n";
		print "$cur_answer1\n";
	    }
	    else{
		$cur_answer1 = <STDIN>;
		chomp($cur_answer1);	    
	    }
	}
	
	if($cur_answer1 eq "y"){
	    my $cur_answer2 = "";
	    print "\n";
	    print "Please specify stack_p_value_threshold (recommended $stack_p_value_threshold): ";
	    $cur_answer2 = <STDIN>;
	    
	    chomp($cur_answer2);
	    $stack_p_value_threshold = $cur_answer2;
	}

	print "\n";
	print "Stacks are collapsed when tags are sufficiently sparse within the collapse window\n";
	print "since true binding sites have both high stacking of tags and high frequency of sampled positions.\n";
	print "QuEST recommends collapsing tags when $recommended_percent_positions_hit_threshold or less percent";
	print " of positions in the collapse window\n";
	print "are sampled by read starts.\n";
	
	$cur_answer1 = "";

	while($cur_answer1 ne "y" && $cur_answer1 ne "n"){
	    print "Would you like to change the percent of positions hit threshold? (y/n, recommended n): ";
	    if($silent_flag eq "on"){
		$cur_answer1 = "n";
		print "$cur_answer1\n";
	    }
	    else{
		$cur_answer1 = <STDIN>;
		chomp($cur_answer1);
	    }
	}
	
	if($cur_answer1 eq "y"){
	    my $cur_answer2 = "";
	    print "\n";
	    print "Please specify percent of positions hit threshold: (0-100, recommended $recommended_percent_positions_hit_threshold): ";
	    $cur_answer2 = <STDIN>;
	    chomp($cur_answer2);
	    $percent_positions_hit_threshold = $cur_answer2;
	}

	print "\n";
	print "Stacks are collapsed if the minimum number of tags in the stack is met.\n";
	print "QuEST uses $recommended_stack_count_threshold for the minimum number of tags in the stack.\n";
	
	$cur_answer1 = "";
	while($cur_answer1 ne "y" && $cur_answer1 ne "n"){
	    print "Would you like to override the count threhsold of a stack? (y/n, recommended n): ";
	    if($silent_flag eq "on"){
		$cur_answer1 = "n";
		print "$cur_answer1\n";
	    }
	    else{
		$cur_answer1 = <STDIN>;
		chomp($cur_answer1);
	    }
	}

	if($cur_answer1 eq "y"){
	    print "Please specify the minimum stack height (recommended $recommended_stack_count_threshold): ";
	    my $cur_answer2 = <STDIN>;
	    chomp($cur_answer2);
	    $stack_count_threshold = $cur_answer2;
	}

	print "\n";
	print "Stacks can be replaced with either one instance of a tag or removed\n";
	print "altogether. We recommend collapsing to a single instance.\n";
	print "\n";

	$cur_answer1 = "";
	while($cur_answer1 ne "y" && $cur_answer1 ne "n"){
	    print "Do you agree? (y/n, recommended y): ";
	    $cur_answer1 = <STDIN>;
	    chomp($cur_answer1);
	}
	
	if($cur_answer eq "n"){
	    $new_stack_size = 0;
	}
    }
}


if($collapse_tags eq "true"){
    print "Collapsing ChIP data...\n";
}
else{
    print "Converting ChIP data to binary format...\n";
}

my $ChIP_collapsed_align_fname = $data_path . "/ChIP.QuEST";
my $ChIP_collapsed_read_num_fname = $data_path . "/ChIP_read_num.txt";

$cur_system_command = "$align_2_bin_exec align_file=$ChIP_uncollapsed_align_fname output_path=$bin_align_ChIP_path genome_table=$contig_table_fname QuEST_collapsed_file=$ChIP_collapsed_align_fname collapse_reads=$collapse_tags collapse_window=$collapse_window count_threshold=$stack_count_threshold stack_p_value_threshold=$stack_p_value_threshold percent_positions_hit_threshold=$percent_positions_hit_threshold report_file=$ChIP_collapsed_read_num_fname new_stack_size=$new_stack_size";

$cur_system_command_out = system($cur_system_command);

if($cur_system_command_out != 0){
    print "\n\n";
    print "Error during the execution of command: \n";
    print "$cur_system_command\n";
    exit(1);
}

if($collapse_tags eq "false"){
    $cur_system_command = "cp $ChIP_uncollapsed_align_fname $ChIP_collapsed_align_fname";
    system($cur_system_command);
}

open ChIP_read_num_file, "< $ChIP_collapsed_read_num_fname" || die "Failed to open $ChIP_collapsed_read_num_fname\n";
my @ChIP_read_num_file_contents = <ChIP_read_num_file>;
close(ChIP_read_num_file);

$ChIP_tags = $ChIP_read_num_file_contents[0];

chomp($ChIP_tags);
print "After collapsing, there are $ChIP_tags ChIP read alignments\n";

if($ChIP_tags <= 0){
    print "No ChIP reads were found in your data. Please check your alignment files.\n";
    exit(4);
}

## creating RX_noIP file

my $RX_noIP_uncollapsed_align_fname = $data_path . "/RX_noIP_uncollapsed.QuEST";

my $RX_noIP_collapsed_align_fname = "";
my $RX_noIP_collapsed_read_num_fname = "";

unlink($RX_noIP_uncollapsed_align_fname);

my $RX_noIP_tags = 0;
my $RX_noIP_uncollapsed_tags = 0;

if($RX_noIP_align_fname ne $missing){
    
    print "\nConverting RX_noIP reads file into the QuEST format.\n";
    print "If the counter below stays at zero, the program can't\n";
    print "match the RX_noIP alignments to the reference genome that you provided.\n";
    print "In this case check that the sequences you have provided match the formats\n";
    print "accepted by QuEST (see README.txt) and that names of contigs in the alignment\n";
    print "files are *exactly* the same as in your genome or genome table.\n\n";
    
    
    my $RX_noIP_uncollapsed_read_num_fname = $data_path . "/RX_noIP_uncollapsed_read_num.txt";
    if(-e $RX_noIP_uncollapsed_read_num_fname){
	unlink($RX_noIP_uncollapsed_read_num_fname);
    }
    
    if($RX_noIP_file_type eq "solexa_align"){
	$cur_system_command = "$align_2_QuEST_exec align_file=$RX_noIP_align_fname output_file=$RX_noIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=solexa_align report_file=$RX_noIP_uncollapsed_read_num_fname";
    }
    elsif($RX_noIP_file_type eq "QuEST_align"){
	$cur_system_command = "$align_2_QuEST_exec align_file=$RX_noIP_align_fname output_file=$RX_noIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=QuEST report_file=$RX_noIP_uncollapsed_read_num_fname";
    }
    elsif($RX_noIP_file_type eq "eland_align"){
	$cur_system_command = "$align_2_QuEST_exec align_file=$RX_noIP_align_fname output_file=$RX_noIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=eland report_file=$RX_noIP_uncollapsed_read_num_fname";
    }
    elsif($RX_noIP_file_type eq "eland_extended_align"){
	$cur_system_command = "$align_2_QuEST_exec align_file=$RX_noIP_align_fname output_file=$RX_noIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=eland_extended report_file=$RX_noIP_uncollapsed_read_num_fname";
    }
    elsif($RX_noIP_file_type eq "bowtie_align"){
	$cur_system_command = "$align_2_QuEST_exec align_file=$RX_noIP_align_fname output_file=$RX_noIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=bowtie report_file=$RX_noIP_uncollapsed_read_num_fname";
    }
    elsif($RX_noIP_file_type eq "maq_align"){
	$cur_system_command = "$align_2_QuEST_exec align_file=$RX_noIP_align_fname output_file=$RX_noIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=maq report_file=$RX_noIP_uncollapsed_read_num_fname";
    }
    elsif($RX_noIP_file_type eq "sam_align"){
	$cur_system_command = "$align_2_QuEST_exec align_file=$RX_noIP_align_fname output_file=$RX_noIP_uncollapsed_align_fname genome_table=$contig_table_fname align_type=sam report_file=$RX_noIP_uncollapsed_read_num_fname";
    }
    else{
	print "Failed to understand the aligment type.\n";
	exit;
    }
    
    $cur_system_command_out = system($cur_system_command);
    if($cur_system_command_out != 0){
	print "Error during the execution of $cur_system_command. Error code: $cur_system_command_out.\n";
	exit(1);
    }

    open RX_noIP_uncollapsed_read_num_file, "< $RX_noIP_uncollapsed_read_num_fname" || die "Failed to open $RX_noIP_uncollapsed_read_num_fname\n";
    my @RX_noIP_uncollapsed_read_num_file_contents = <RX_noIP_uncollapsed_read_num_file>;
    close(RX_noIP_uncollapsed_read_num_file);
    
    $RX_noIP_uncollapsed_tags = $RX_noIP_uncollapsed_read_num_file_contents[0];
    chomp($RX_noIP_uncollapsed_tags);
    
    print "Collapsing RX_noIP (control) data...\n";
    
    $RX_noIP_collapsed_align_fname = $data_path . "/RX_noIP.QuEST";
    $RX_noIP_collapsed_read_num_fname = $data_path . "/RX_noIP_read_num.txt";
    
    # just generate a new QuEST align file for now, convert to binary in the FDR step
    $cur_system_command = "$align_2_bin_exec align_file=$RX_noIP_uncollapsed_align_fname output_path=not_applicable genome_table=$contig_table_fname QuEST_collapsed_file=$RX_noIP_collapsed_align_fname collapse_reads=$collapse_tags collapse_window=$collapse_window count_threshold=$stack_count_threshold stack_p_value_threshold=$stack_p_value_threshold percent_positions_hit_threshold=$percent_positions_hit_threshold report_file=$RX_noIP_collapsed_read_num_fname new_stack_size=$new_stack_size";
    
    $cur_system_command_out = system($cur_system_command);
    
    if($cur_system_command_out != 0){
	print "\n\n";
	print "Error during the execution of command: \n";
	print "$cur_system_command\n";
	exit(1);
    }
    
    open RX_noIP_read_num_file, "< $RX_noIP_collapsed_read_num_fname" || die "Failed to open $RX_noIP_collapsed_read_num_fname\n";
    my @RX_noIP_read_num_file_contents = <RX_noIP_read_num_file>;
    close(RX_noIP_read_num_file);
    
    $RX_noIP_tags = $RX_noIP_read_num_file_contents[0];
    chomp($RX_noIP_tags);

}
else{
#    $RX_noIP_collapsed_align_fname = $data_path . "/RX_noIP.QuEST";
#    $RX_noIP_collapsed_read_num_fname = $data_path . "/RX_noIP_read_num.txt";
    $RX_noIP_tags = 0;
    system("touch $RX_noIP_uncollapsed_align_fname");
    $RX_noIP_collapsed_align_fname = $data_path . "/RX_noIP.QuEST";
    system("touch $RX_noIP_collapsed_align_fname");
#    $cur_system_command = "$align_2_bin_exec align_file=$RX_noIP_uncollapsed_align_fname output_path=not_applicable genome_table=$contig_table_fname QuEST_collapsed_file=$RX_noIP_collapsed_align_fname collapse_reads=$collapse_tags collapse_window=$collapse_window count_threshold=$stack_count_threshold stack_p_value_threshold=$stack_p_value_threshold percent_positions_hit_threshold=$percent_positions_hit_threshold report_file=$RX_noIP_collapsed_read_num_fname new_stack_size=$new_stack_size";
#    
#    $cur_system_command_out = system($cur_system_command);
#    
#    if($cur_system_command_out != 0){
#	print "\n\n";
#	print "Error during the execution of command: \n";
#	print "$cur_system_command\n";
#	exit(1);
#    }
}
print "After collapsing, there were $RX_noIP_tags RX_noIP read alignments\n";


## creating RX_noIP file

## /creating RX_noIP file 

## make assessment to the number of tags

my $ChIP_recommended_tag_count = int( ( $genome_size / 3080436051 ) * 5000000);
my $ChIP_excess_fold = $ChIP_tags / $ChIP_recommended_tag_count;

my $RX_noIP_recommended_tag_count = ( $genome_size / 3080436051 ) * 7000000;
my $RX_noIP_excess_fold = $RX_noIP_tags / $RX_noIP_recommended_tag_count;
my $RX_noIP_excess_fold_with_FDR = ($RX_noIP_tags - $ChIP_tags ) / $RX_noIP_recommended_tag_count;

print "\n";
printf "You seem to have at least %.2f fold the recommended amount of ChIP reads.\n", $ChIP_excess_fold;
printf "You seem to have at least %.2f fold the recommended amount of RX_noIP reads without FDR estimate.\n", $RX_noIP_excess_fold;
if($ChIP_tags > $RX_noIP_tags){
    print "You do not have enough RX_noIP for FDR estimation. Will procede without FDR.\n";
}
else{
    printf "You seem to have at least %.2f fold the recommended amount of RX_noIP reads with FDR estimate.\n", $RX_noIP_excess_fold_with_FDR;
}
print "\n";

if($ChIP_excess_fold < 1){
    print "\nWarning: you seem to not have sufficiently many ChIP sequencing reads.\n";
    print "We recommend at least $ChIP_recommended_tag_count aligned reads,\n";
    print "however you only provided $ChIP_tags. QuEST will still work, but your results\n";
    print "may not be very accurate. \n";
    my $answer = "";
    while( $answer ne "y" && $answer ne "n" ){
	print "Do you want to continue? (y/n): ";
	if($silent_flag eq "on"){
	    $answer = "y";
	    print "$answer\n";
	}
	else{
	    $answer = <STDIN>;
	}
	chomp($answer);
	if($answer eq "y"){
	    # do nothing
	}
	elsif($answer eq "n"){
	    print "You have aborted QuEST configuration.\n";
	    exit(0);
	}
    }
}
else{
    printf "You seem to have at least %.2f fold the recommended amount of ChIP sequence reads.\n", $ChIP_excess_fold;
}

if($RX_noIP_excess_fold < 1){

    if($RX_noIP_tags > 0){
	print "\nWarning: you we seem to not have sufficiently many RX_noIP sequencing reads.\n";
	printf "We recommend at least %.0f aligned reads,\n", $RX_noIP_recommended_tag_count;
	print "however you only provided $RX_noIP_tags. QuEST will still work, but your results\n";    
	print "may not be very accurate. \n";
	
	my $answer = "";
	while( $answer ne "y" && $answer ne "n" ){
	    print "Do you want to continue? (y/n): ";
	    if($silent_flag eq "on"){
		$answer = "y";
		print "$answer\n";
	    }
	    else{
		$answer = <STDIN>;
	    }	
	    chomp($answer);
	    if($answer eq "y"){
		# do nothing
	    }
	    elsif($answer eq "n"){
		print "You have aborted QuEST configuration.\n";
		exit(0);
	    }
	}
    }
    else{
	print "\n";
	print "Because you have not provided any control data, QuEST will use uniform control\n";
	print "distribution for peak calling.\n";
	print "\n";
    }
}
else{
    if( $RX_noIP_excess_fold_with_FDR >=1 ){
	print "You have enough RX_noIP reads to run the analysis with FDR estimation from control data.\n";
    }
    else{
	
	print "You may not have enought RX_noIP reads to run the analysis with FDR.\n";
	#print "However, QuEST will still allocate the pseudo_ChIP data as you requested.\n";
    }
}

my $analysis_with_FDR = "not defined";
my $background_tags = -1;
my $pseudo_ChIP_tags = -1;


if( $RX_noIP_align_fname eq $missing){ ## sanity check for positive tag count in the control data
    
}
else{
    if($RX_noIP_tags == 0){
	print "No background tags were found in your RX_noIP file.\n";
	print "If you don't have any control data, omit RX_noIP parameter when running this script.\n";
	exit(0);
    }
}

if($RX_noIP_tags <= $ChIP_tags){
    print "Since the number of Rx_noIP tags is less than ChIP tags, there is no\n";
    print "possibility to provide an FDR estimate. \n";
    print "QuEST will be configured without FDR estimation.\n";
    
    $analysis_with_FDR = "false";
    $background_tags = $RX_noIP_tags;
    
}
else{    
    
    $answer = "";
    while( $answer ne "y" && $answer ne "n" ){	
	print "\n";
	print "Would you like to run QuEST with FDR analysis? (y/n): ";
	if($silent_flag eq "on"){
	    $answer = "y";
	    print "$answer\n";
	}
	else{
	    $answer = <STDIN>;
	    chomp($answer);
	}
    }
    
    if($answer eq "n"){
	$analysis_with_FDR = "false";
	$background_tags = $RX_noIP_tags;
    }
    elsif($answer eq "y"){
	$analysis_with_FDR = "true";
	$pseudo_ChIP_tags = $ChIP_tags;
	$background_tags = $RX_noIP_tags - $pseudo_ChIP_tags;
    }
}

## /make assessment to the number of tags

## allocating RX_noIP_tags
my $background_RX_noIP_fname = "";
my $pseudo_ChIP_RX_noIP_fname = "";

if($analysis_with_FDR eq "true"){
    print "Randomizing and allocating pseudo_ChIP and background data\n";
    
    $background_RX_noIP_fname = $data_path . "/background_RX_noIP.align.txt";
    $pseudo_ChIP_RX_noIP_fname = $data_path . "/pseudo_ChIP_RX_noIP.align.txt";
    
    my @RX_noIP_allocation_mask;
    
    for(my $i=0; $i<$RX_noIP_tags; $i++){
	if($i<$pseudo_ChIP_tags){
	    $RX_noIP_allocation_mask[$i] = "p";
	}
	else{ $RX_noIP_allocation_mask[$i] = "b";}
    }
    
    # permuting
    for(my $i=0; $i<$RX_noIP_tags; $i++){
	if($i%10000 == 0){
	    printf "\r%.2f M / %.2f M", $i/1000000, $RX_noIP_tags/1000000;
	    $|++;
	}
	my $random_number=  int(rand($RX_noIP_tags));
	my $tmp_var = $RX_noIP_allocation_mask[$i];
	$RX_noIP_allocation_mask[$i] = $RX_noIP_allocation_mask[$random_number];
	$RX_noIP_allocation_mask[$random_number] = $tmp_var;
    }
    # /permuting

    open background_RX_noIP_file, "> $background_RX_noIP_fname" || 
	die "Failed to open $background_RX_noIP_fname\n";
    open pseudo_ChIP_RX_noIP_file, "> $pseudo_ChIP_RX_noIP_fname" || 
	die "Failed to open $pseudo_ChIP_RX_noIP_fname\n";
    open RX_noIP_align_file, "< $RX_noIP_collapsed_align_fname" || 
	die "Failed to open $RX_noIP_collapsed_align_fname\n";
 

    #print "opening $RX_noIP_collapsed_align_fname...\n";

    my $cur_line_counter = 0;
    while(<RX_noIP_align_file>){
	#print "$RX_noIP_allocation_mask[$cur_line_counter]\n";
	chomp;
	if($RX_noIP_allocation_mask[$cur_line_counter] eq "p"){
	    print pseudo_ChIP_RX_noIP_file "$_\n";	    
	}
	elsif($RX_noIP_allocation_mask[$cur_line_counter] eq "b"){
	    print background_RX_noIP_file "$_\n";	    
	}
	$cur_line_counter++;
    }
    close background_RX_noIP_file;
    close pseudo_ChIP_RX_noIP_file;
   # close updated_RX_noIP_align_file;

    # create binary files now
    
    $cur_system_command = "$align_2_bin_exec align_file=$pseudo_ChIP_RX_noIP_fname output_path=$bin_align_pseudo_ChIP_path genome_table=$contig_table_fname QuEST_collapsed_file=not_applicable collapse_reads=false report_file=not_applicable";

    $cur_system_command_out = system($cur_system_command);
    
    if($cur_system_command_out != 0){
	print "\n\n";
	print "Error during the execution of command: \n";
	print "$cur_system_command\n";
	exit(1);
    }

        
    $cur_system_command = "$align_2_bin_exec align_file=$background_RX_noIP_fname output_path=$bin_align_background_path genome_table=$contig_table_fname QuEST_collapsed_file=not_applicable collapse_reads=false report_file=not_applicable";
    
    $cur_system_command_out = system($cur_system_command);
    
    if($cur_system_command_out != 0){
	print "\n\n";
	print "Error during the execution of command: \n";
	print "$cur_system_command\n";
	exit(1);
    }
    
#    $cur_system_command = "$align_2_bin_exec align_file=$RX_noIP_collapsed_align_fname output_path=$bin_align_RX_noIP_path genome_table=$contig_table_fname QuEST_collapsed_file=not_applicable collapse_reads=false report_file=not_applicable";
#    
#    $cur_system_command_out = system($cur_system_command);
#    
#    if($cur_system_command_out != 0){
#	print "\n\n";
#	print "Error dudring the execution of command: \n";
#	print "$cur_system_command\n";
#	exit(1);
#    }

}
elsif($analysis_with_FDR eq "false"){ # && $RX_noIP_tags > 0){ 
#    $background_tags = $RX_noIP_tags;
#    $background_RX_noIP_fname = $RX_noIP_uncollapsed_align_fname;
    $background_RX_noIP_fname = $RX_noIP_collapsed_align_fname;

    # create binary files now
    
    $cur_system_command = "$align_2_bin_exec align_file=$background_RX_noIP_fname output_path=$bin_align_RX_noIP_path genome_table=$contig_table_fname QuEST_collapsed_file=not_applicable collapse_reads=false report_file=not_applicable";

    $cur_system_command_out = system($cur_system_command);
    
    if($cur_system_command_out != 0){
	print "\n\n";
	print "Error during the execution of command: \n";
	print "$cur_system_command\n";
	exit(1);
    }
 
    my @bin_files = <$bin_align_RX_noIP_path/*>;
    #foreach $file in (@bin_files){
    for(my $i=0; $i<scalar(@bin_files); $i++){
	my $file=$bin_files[$i];
	$cur_system_command = "cp $file $bin_align_background_path";
	print "$cur_system_command\n";
	my $cur_system_command_out = system($cur_system_command);
	if($cur_system_command_out != 0){
	    print "\n\n";
	    print "Error during the execution of command: \n";
	    print "$cur_system_command\n";
	    exit(1);
	}
    }
    #$cur_system_command = "cp $bin_align_RX_noIP_path/* $bin_align_background_path";
    #$cur_system_command_out = system($cur_system_command);
    #if($cur_system_command_out != 0){
#	print "\n\n";
#	print "Error during the execution of command: \n";
#	print "$cur_system_command\n";
#	exit(1);
#    }
}

## calculating cutoffs/thresholds

print "\nWe will now determine necessary QuEST thresholds.\n\n";

my $screen_msg = "";

my $positive_region_size;
my $recommended_positive_region_size = 300;

my $KDE_bandwidth = -1;
my $recommended_KDE_bandwidth = 30;

my $ps_lower_threshold = $undef;
my $ps_upper_threshold = $undef;

my $cor_lower_threshold = $undef;

my $filter_out_empty_regions = "false";
my $report_peaks  = "all";

$screen_msg = qq{
    ----------- QuEST parameter configuration -------------
};

print "$screen_msg";

$screen_msg =  qq{
What kind of ChIP experiment is this?

1. Transcription factor with defined motif and narrow (punctate) binding site
   resulting in regions of enrichment 100-300 bp wide.
   (bandwidth = 30 bp, region_size = 300 bp)

2. PolII-like factor resulting in regions of enrichment 300-1000 bp with
   some narrow binding sites and some wide sites possibly occupying
   the entire gene length.
   (bandwidth = 60 bp, region_size = 600 bp)

3. Histone-type ChIP resulting in wide regions of enrichment
   1Kb and up possibly occupying multiple genes.
   (bandwidth = 100 bp, region_size = 1000 bp)

4. Neither, I would like to configure individual parametrs myself.

};

print "$screen_msg";

my $selection = "";
while( $selection ne "1" && $selection ne "2" &&
       $selection ne "3" && $selection ne "4"){
    print "Please specify your selection [1-4]: ";
    $selection = <STDIN>;
    chomp($selection);

}

if($selection eq "1"){
    $positive_region_size = 300;
    $KDE_bandwidth = 30;
    $ps_lower_threshold = 25;
    $filter_out_empty_regions = "true";
    $report_peaks = "best";
#    $ps_upper_threshold = 125;
#    $cor_lower_threshold = 0.5;
}
elsif($selection eq "2"){
    $positive_region_size = 600;
    $KDE_bandwidth = 60;
}
elsif($selection eq "3"){
    $positive_region_size = 1000;
    $KDE_bandwidth = 100;
}
elsif($selection eq "4"){
    print "\n";
    print "KDE bandwidth defines smoothing parameter for the ChIP density\n";
    print "estimation. Typical values are in the range of 30-100 bps.\n";
    print "\n";

    $answer = "";
    print "Please specify the value of KDE bandwidth (recommended 30-100): ";
    $answer = <STDIN>;
    chomp($answer);

    $KDE_bandwidth = $answer;
    
    print "\n";
    print "ChIP regions are used to calculate enrichment statistics, and\n";
    print "estimate peak shift parameter that is then applied for the entire\n";
    print "data. \n";
    print "\n";
    
    print "A good way to evaluate the size of the region is by loading the\n";
    print "ChIP-Seq data into the genome browser and inspecting likely targets\n";
    print "to see how big the enriched regions are.\n";
    print "\n";
   
    print "Please specify the size of the region (recommended 300-1000 ): ";
    $answer = "";
    $answer = <STDIN>;
    
    chomp($answer);
    
    $positive_region_size = $answer;
}


#my ($ChIP_threshold, $background_threshold, $er, $rr, $dip_fraction);

#my $imaginable_background_tags = 10000000;

#my $recommended_ChIP_threshold = 0.3 * ( 3080436051 / $genome_size ) * ( $ChIP_tags / 7000000 ); #old, convservative threshold

my $recommended_mappable_genome_fraction = 0.75;
my $mappable_genome_fraction = -1;

if($advanced_flag eq "on"){
    print "Please enter the mappable genome fraction (0-1, recommeded $recommended_mappable_genome_fraction): ";
    my $answer1 = <STDIN>;
    chomp($answer1);
    if($answer1 >=0 && $answer1 <=1 ){
	$mappable_genome_fraction = $answer1;
    }
    else{
	print "Error: mappble genome fraction should be between 0 and 1, but you provided $answer1\n";
	exit(1);
    }
}
else{
    $mappable_genome_fraction = $recommended_mappable_genome_fraction;
}

my $ChIP_basal_level = $ChIP_tags / ($genome_size * $mappable_genome_fraction);
my $background_basal_level = -1;
if($background_tags == 0){
    $background_basal_level = "NA"; 
    #$imaginable_background_tags / ($genome_size * $mappable_genome_fraction);
}
else{
    $background_basal_level = $background_tags / ($genome_size * $mappable_genome_fraction);    
}



my $pseudo_ChIP_basal_level =  "NA";
if($analysis_with_FDR eq "true"){
    my $pseudo_ChIP_basal_level = $pseudo_ChIP_tags / ($genome_size * $mappable_genome_fraction);
}

my $background_used = "false";
if($background_tags > 0){
    $background_used = "true";
}
else{
    $background_used = "false";
}

# stringent parameters
my $stringent_ChIP_enrichment = 50;
my $stringent_ChIP_extension_enrichment = 3;

my $stringent_ChIP_to_background_enrichment;
if($background_used eq "true"){
    $stringent_ChIP_to_background_enrichment = 3.0;
}
else{
    $stringent_ChIP_to_background_enrichment = "NA";
}
#my $stringent_ChIP_to_background_rescue_ratio = 3.0;

my $recommended_ChIP_enrichment = 30;
my $recommended_ChIP_extension_enrichment = 3.0;
my $recommended_ChIP_to_background_enrichment;
if($background_used eq "true"){
   $recommended_ChIP_to_background_enrichment = 3.0;
}
else{
    $recommended_ChIP_to_background_enrichment = "NA";
}
#y $recommended_ChIP_to_background_rescue_ratio = 3.0;

my $relaxed_ChIP_enrichment = 10.0;
my $relaxed_ChIP_extension_enrichment = 3.0;
my $relaxed_ChIP_to_background_enrichment;
if($background_used eq "true"){
    $relaxed_ChIP_to_background_enrichment = 2.5;
}
else{
    $relaxed_ChIP_to_background_enrichment = "NA";
}

my $recommended_dip_fraction = 0.1;

my $ChIP_threshold = "";
my $ChIP_extension_threshold = "";
my $ChIP_to_background_ratio = "NA";

#my $background_threshold = "NA";
my $dip_fraction = "";
#my $rr = "";

$dip_fraction = $recommended_dip_fraction;

print "\n";
print "Let\'s now determine QuEST peak calling parameters.\n";
print "(Parameters below are optimized for human and mouse data\n";
print "when default stack collapsing is performed).\n";
print "\n";
print "Choose one of the following options:\n";
print "\n";
print "1. Stringent peak calling parameters.\n";
print "   ChIP enrichment                 = $stringent_ChIP_enrichment\n";
print "   ChIP to background enrichment   = $stringent_ChIP_to_background_enrichment\n";
print "   ChIP extension enrichment       = $stringent_ChIP_extension_enrichment\n";

print "\n";    
print "2. Recommended peak calling parameters.\n";
print "   ChIP enrichment                 = $recommended_ChIP_enrichment\n";
print "   ChIP to background enrichment   = $recommended_ChIP_to_background_enrichment\n";
print "   ChIP extension enrichment       = $recommended_ChIP_extension_enrichment\n";

print "\n";
print "3. Relaxed peak calling parameters.\n";
print "   ChIP enrichment                 = $relaxed_ChIP_enrichment\n";
print "   ChIP to background enrichment   = $relaxed_ChIP_to_background_enrichment\n";
print "   ChIP extension enrichment       = $relaxed_ChIP_extension_enrichment\n";

print "\n";
print "4. Neither, I want to specify peak calling parameters myself.\n";
print "\n";


$selection = "";
while( $selection ne "1" && $selection ne "2" &&
    $selection ne "3" && $selection ne "4"){
    print "Please specify your selection [1-4]: ";
    $selection = <STDIN>;
    chomp($selection);
}

my $ChIP_seeding_enrichment = -1;
my $ChIP_extension_enrichment = -1;
my $ChIP_to_background_enrichment = "NA";

if($selection eq "1" || $selection eq "2" || $selection eq "3"){
    

    if($selection eq "1"){
	$ChIP_seeding_enrichment = $stringent_ChIP_enrichment;
	$ChIP_extension_enrichment  = $stringent_ChIP_extension_enrichment;
	if($background_tags > 0){
	    $ChIP_to_background_enrichment = $stringent_ChIP_to_background_enrichment;
	}
    }
    if($selection eq "2"){
	$ChIP_seeding_enrichment = $recommended_ChIP_enrichment;
	$ChIP_extension_enrichment  = $recommended_ChIP_extension_enrichment;
	if($background_tags > 0){
	    $ChIP_to_background_enrichment = $recommended_ChIP_to_background_enrichment;
	}
    }
    if($selection eq "3"){
	$ChIP_seeding_enrichment = $relaxed_ChIP_enrichment;
	$ChIP_extension_enrichment  = $relaxed_ChIP_extension_enrichment;
	if($background_tags > 0){
	    $ChIP_to_background_enrichment = $relaxed_ChIP_to_background_enrichment;
	}
    }
}
elsif($selection eq "4"){
    $answer = "";
    print "Please specifty ChIP enrichment (1 or greater, recommended $recommended_ChIP_enrichment): ";
    $answer = <STDIN>;
    chomp($answer);
    
    my $user_ChIP_enrichment = $answer;
    
    #if($answer > 1){
    #}
    #else{
#	print "The value cannot be less than 1, you have specified $answer\n. Exiting.\n";
#	exit(1);
#    }
    
    my $user_ChIP_to_background_enrichment = "NA";
    if($background_tags > 0){
	print "Please specifty ChIP to background enrichment (1 or greater, recommended $recommended_ChIP_to_background_enrichment): ";
	
	$answer = <STDIN>;
	chomp($answer);
	
	$user_ChIP_to_background_enrichment = $answer;
	
	if($answer >= 1){	    
	}
	else{
	    print "The value cannot be less than 1, you have specified $answer\n. Exiting.\n";
	    exit(1);
	}
    }
    
    print "Please specify the ChIP extension enrichment (1 or greater, recommended $recommended_ChIP_extension_enrichment): ";
    $answer = <STDIN>;
    chomp($answer);
    
    my $user_ChIP_extension_enrichment = $answer;
#    if($answer >= 1){
#    }
#    else{
#	print "The value cannot be less than 1, you have specified $answer\n. Exiting.\n";
#	exit(1);
#    }

    
    
    $ChIP_seeding_enrichment = $user_ChIP_enrichment;
    $ChIP_extension_enrichment = $user_ChIP_extension_enrichment;
    if($background_tags > 0){
	$ChIP_to_background_enrichment = $user_ChIP_to_background_enrichment;
    }
}
else{
    print "Wrong selection option\n";
    exit(0);
}

$ChIP_threshold = $ChIP_basal_level * $ChIP_seeding_enrichment;
$ChIP_extension_threshold = $ChIP_basal_level * $ChIP_extension_enrichment;
if($background_tags > 0){
    $ChIP_to_background_ratio = ($ChIP_tags / $background_tags) * $ChIP_to_background_enrichment;
}

$screen_msg = qq{
The following parameters were specified for running QuEST:

positive_region_size               :  $positive_region_size
KDE_bandwidth                      :  $KDE_bandwidth
ChIP_seeding_fold_enrchment        :  $ChIP_seeding_enrichment
ChIP_extension_fold_enrichment     :  $ChIP_extension_enrichment
ChIP_to_background_fold_enrichment :  $ChIP_to_background_enrichment
};
print "$screen_msg";
## /calculating Peak Caller thresholds
    
## calculating Quick Window Scan thresholds

#print "\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#print "~~~~~~~~~~~~~~ Quick Window Scan Module ~~~~~~~~~~~~~~~~~~~~~\n";


my $recommended_quick_window_scan_calc_window = $positive_region_size;


my $recommended_ChIP_tag_enrichment = 20.0;
my $recommended_ChIP_to_background_tag_enrichment = 3.0;

my $ChIP_basal_tag_rate = $ChIP_tags / ($genome_size * $mappable_genome_fraction);
my $background_basal_tag_rate = "NA";
if($background_tags > 0){
    $background_basal_tag_rate = $background_tags / ($genome_size * $mappable_genome_fraction);
}



my $ChIP_tag_threshold;
my $ChIP_to_background_tag_ratio = "NA";
my $quick_window_scan_calc_window;

if($advanced_flag eq "off"){
    $quick_window_scan_calc_window = $recommended_quick_window_scan_calc_window;
    my $basal_ChIP_tag_count = $quick_window_scan_calc_window * $ChIP_basal_tag_rate;    
    $ChIP_tag_threshold = $basal_ChIP_tag_count * $recommended_ChIP_tag_enrichment;

    if($background_tags > 0){
	$ChIP_to_background_tag_ratio = ($ChIP_tags / $background_tags) * $recommended_ChIP_to_background_tag_enrichment;
    }
    else{
	$ChIP_to_background_tag_ratio = "NA";
    }
}
elsif($advanced_flag eq "on"){
    
    print "\n";
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    print "~~~~~~~~~~~~~~ Quick Window Scan Module ~~~~~~~~~~~~~~~~~~~~~\n";
    
    
    print "\n";
    print "Quick window scan is the module that identifies a training\n";
    print "set of potential peaks for estimating peak shift parameter.\n";
    print "To do so, it needs the size of the window to count the tags in\n";
    print "and a count threshold to identify positive regions where that\n";
    print "threshold is exceeded.\n\n";

#    print "Based on the size of the genome and the number of ChIP tags, we\n";
#    print "recommend using $recommended_quick_window_scan_calc_window for the window size and ";
#    printf "%.0f for the the ", $recommended_quick_window_scan_count_threshold;
#    print "tag count threshold.\n";
    
    $answer = "";

    print "Please specify the size of the quick scan window (recommended $recommended_quick_window_scan_calc_window): ";
    $answer = <STDIN>;
    chomp($answer);

    if($answer > 0){
	$quick_window_scan_calc_window = $answer;	
    }
    else{
	print "You have to provide a positive number, but you gave $answer.\nExiting.\n";
	exit(1);
    }
    
    my $basal_ChIP_tag_count = $quick_window_scan_calc_window * $ChIP_basal_tag_rate;
    my $recommended_ChIP_tag_threshold = $basal_ChIP_tag_count * $recommended_ChIP_tag_enrichment;

    print "Please specify the ChIP tags enrichment\n";
    printf "(1 or above, recommended %.2f corresponding to %.0f ChIP tags): ", $recommended_ChIP_tag_enrichment, $recommended_ChIP_tag_threshold;
    $answer =  <STDIN>;
    chomp($answer);
    
    if($answer > 1){
	my $ChIP_tag_enrichment = $answer;

	$ChIP_tag_threshold = $basal_ChIP_tag_count * $ChIP_tag_enrichment;

	if($background_tags > 0){
	    $ChIP_to_background_tag_ratio = ($ChIP_tags / $background_tags) * $recommended_ChIP_to_background_tag_enrichment;
	}
	else{
	    $ChIP_to_background_tag_ratio = "NA";
	}
    }
   
    print "\n";
    print "Specified parameters:\n\n";
    print "quick_window_scan_calc_window:         $quick_window_scan_calc_window\n";
    print "ChIP_tag_threshold:                    $ChIP_tag_threshold\n";
    print "ChIP_to_background_tag_ratio:          $ChIP_to_background_tag_ratio\n";
    
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
    
## /calculating Quick Window Scan thresholds
    
}
## calculating Peak Shift Calibrator thresholds

#print "\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#print "~~~~~~~~~~~~~ Peak Shift Calibrator Module ~~~~~~~~~~~~~~~~~~\n";


my $recommended_PSC_region_width = $positive_region_size;
my $recommended_top_regions_passed = 200;
my $recommended_correlation_lower_threshold = 0.7;
my $recommended_peak_shift_lower_threshold = 25;
my $recommended_peak_shift_upper_threshold = 150;
my $recommended_dist_threshold = 2*($recommended_peak_shift_upper_threshold + 50);
my $recommended_os_upper_threshold = 10;

my $PSC_region_width = $recommended_PSC_region_width;
my $top_regions_passed = $recommended_top_regions_passed;
my $correlation_lower_threshold = $recommended_correlation_lower_threshold;
my $peak_shift_lower_threshold = $recommended_peak_shift_lower_threshold;
my $peak_shift_upper_threshold = $recommended_peak_shift_upper_threshold;
my $dist_threshold = $recommended_dist_threshold;
my $os_upper_threshold = $recommended_os_upper_threshold;
my $peak_shift_estimation_method = "mode";
#my $peak_shift_estimation_method = "median";

if($advanced_flag eq "off"){
    # don't do anything, all the defaults are already set
}
elsif($advanced_flag eq "on"){
    
    print "\n";
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    print "~~~~~~~~~~~~~ Peak Shift Calibrator Module ~~~~~~~~~~~~~~~~~~\n";
    
    
    print "\nCalibrate Peak Shift Module estimates the value of peak shift based on the\n";
    print "the list of the regions provided by the Quick Window Scan Module.\n\n";
    

    $answer = "";

    
    print "\nYou can limit the total number of regions on which to base your estimate.\n";
    print "We recommend limiting this number to $recommended_top_regions_passed.\n";

    $answer = "";
    while( $answer ne "y" && $answer ne "n" ){
	print "Do you agree? (y/n): ";
	
	if($silent_flag eq "on"){
	    $answer = "y";
	    print "$answer\n";
	}
	else{
	    $answer = <STDIN>;
	}
	chomp($answer);
	
	if( $answer eq "y" ){
	    $top_regions_passed = $recommended_top_regions_passed
	    }
	elsif( $answer eq "n" ){
	    print "Please specify the number of regions to use for peak_shift estimation (recommended $recommended_top_regions_passed): ";
	    my $cur_answer = <STDIN>;
	    chomp($cur_answer);
	    if( $cur_answer >= 0 ){
		$top_regions_passed = $recommended_top_regions_passed;
	    }
	    else{ print "You should specify a non-negative number. Aborting.\n"; exit(0); }    
	}
    }    

    print "\n";
    print "Regions can be filtered based on tag count correlation between the strands.\n";
    print "We recommend using $recommended_correlation_lower_threshold for the correlation threshold.\n";
    
    $answer = "";
    while( $answer ne "y" && $answer ne "n" ){
	print "Do you agree? (y/n): ";
	if($silent_flag eq "on"){
	    $answer = "y";
	    print "$answer\n";
	}
	else{
	    $answer = <STDIN>;
	}
	chomp($answer);
	
	if( $answer eq "y" ){
	    $correlation_lower_threshold = $recommended_correlation_lower_threshold;
	}
	elsif( $answer eq "n" ){
	    print "Please specify the correlation threshold: (recommended $recommended_correlation_lower_threshold): ";
	    $answer = <STDIN>;
	    chomp($answer);
	    if( $answer >= 0 ){
		$correlation_lower_threshold = $answer;
	    }
	    else{ print "You should specify a non-negative number. Aborting.\n"; exit(0); }    
	}
    }

    print "\n";
    print "Specified parameters for the Peak Shift Calibration Module:\n\n";

    print "top_regions_passed:                     $top_regions_passed\n";
    print "correlation_lower_threshold:            $correlation_lower_threshold\n";
    print "peak_shift_lower_threshold:             $peak_shift_lower_threshold\n";
    print "peak_shift_upper_threshold:             $peak_shift_upper_threshold\n";
    print "dist_threshold:                         $dist_threshold\n";

    print "peak_shift_estimation_method:           $peak_shift_estimation_method\n";
	
}
## /calculating Peak Shift Calibrator thresholds

my $strand_correlation_region_width = $positive_region_size;

## system info

my $processors = -1;
my $recommended_threads = -1;
my $threads = -1;

if(-e "/proc/cpuinfo"){
    $processors = `cat /proc/cpuinfo | grep processor | wc -l`;
    chomp($processors);
}
elsif(-e "/usr/sbin/sysctl"){
    $processors = `/usr/sbin/sysctl -n hw.ncpu`;
}
else{
    print "I failed to determine the number of processors on your platform.\n";
    print "Please specify it manually: ";
    $answer = <STDIN>;
    chomp($answer);
    if($answer <= 0){
	print "You need to specify a positive number. Aborting.\n";
	exit(0);
    }
    $processors = $answer;
}

if($processors == 0){
    print "I failed to determine the number of processors on your platform.\n";
    print "Please specify it manually: ";
    $answer = <STDIN>;
    chomp($answer);
    if($answer <= 0){
	print "You need to specify a positive number. Aborting.\n";
	exit(0);
    }
    $processors = $answer;
}
#print "\nNumber of processors: $processors\n";

if($processors == 1){
    $recommended_threads = 1;
}
else{ $recommended_threads = $processors - 1; }

$threads = $recommended_threads;

if($advanced_flag eq "off"){
}
elsif($advanced_flag eq "on"){
    
    print "\nThe Profiler Module supports thread-level parallelism. We suggest utilizing \n";
    print "$recommended_threads for the number of threads to distribute the computation over.\n";
    
    $answer = "";
    while( $answer ne "y" && $answer ne "n" ){
	print "Do you agree? (y/n): ";
	if($silent_flag eq "on"){
	    $answer = "y";
	    print "$answer\n";
	}
	else{
	    $answer = <STDIN>;
	}
	
	chomp($answer);
	
	if($answer eq "y"){
	    $threads = $recommended_threads;
	}
	elsif($answer eq "n"){
	    print "Please specify desired number of threads: ";
	    my $cur_answer2 = <STDIN>;
	    chomp($cur_answer2);
	    $threads = $cur_answer2;
	}
    }
    print "\n";
}
## /system info

## /calculating cutoffs/thresholds


## writing Quick Window Scan parameters

my $quick_window_scan_par_fname = $parameters_path . "/quick_window_scan.batch.pars";
my $quick_window_scan_output_fname = $module_outputs_path . "/quick_window_scan.out";

unlink($quick_window_scan_par_fname);

open quick_window_scan_par_file, "> $quick_window_scan_par_fname" 
    or die "Failed to open $quick_window_scan_par_fname\n";
print quick_window_scan_par_file "genome_table = $contig_table_fname\n";
print quick_window_scan_par_file "bin_align_ChIP_path = $bin_align_ChIP_path\n";
if($background_tags > 0){
    if($analysis_with_FDR eq "false"){
	print quick_window_scan_par_file "bin_align_RX_noIP_path = $bin_align_RX_noIP_path\n";
    }
    elsif($analysis_with_FDR eq "true"){
	print quick_window_scan_par_file "bin_align_RX_noIP_path = $bin_align_background_path\n";
    }
}
else{
    print quick_window_scan_par_file "bin_align_RX_noIP_path = NA\n";
}
#print quick_window_scan_par_file "ChIP_align_file = $updated_ChIP_align_fname\n";
#print quick_window_scan_par_file "RX_noIP_align_file = $updated_RX_noIP_align_fname\n";
print quick_window_scan_par_file "output_file = $quick_window_scan_output_fname\n";
print quick_window_scan_par_file "calc_window = $quick_window_scan_calc_window\n";
printf quick_window_scan_par_file "ChIP_tags_threshold = %.0f\n", $ChIP_tag_threshold;
if($background_tags > 0){
    printf quick_window_scan_par_file "ChIP_to_background_tag_ratio = %f\n", $ChIP_to_background_tag_ratio;
}
else{
    print quick_window_scan_par_file "ChIP_to_background_tag_ratio = NA\n";
}
# = $ChIP_tags_threshold\n";

#printf quick_window_scan_par_file "count_threshold = %d\n", $quick_window_scan_count_threshold;
#printf quick_window_scan_par_file "enrichment_ratio = %.3f\n", $quick_window_scan_enrichment_ratio;

close quick_window_scan_par_file;

#print "Saved Quick Window Scan batch parameter file.\n";

## /writing Quick Window Scan parameters

## writing Calibrate Peak Shift parameters

my $calibrate_peak_shift_par_fname = $parameters_path . "/calibrate_peak_shift.batch.pars";
my $calibrate_peak_shift_output_fname = $module_outputs_path . "/calibrate_peak_shift.out";

open calibrate_peak_shift_par_file, "> $calibrate_peak_shift_par_fname"
    or die "Failed to open $calibrate_peak_shift_par_fname\n";

if($silent_flag eq "on"){
    print calibrate_peak_shift_par_file "verbose = off\n";
}
else{
    print calibrate_peak_shift_par_file "verbose = yes\n";
}
#print calibrate_peak_shift_par_file "quick_scan_file = $quick_window_scan_output_fname.sorted\n";
#print calibrate_peak_shift_par_file "ChIP_align_file = $updated_ChIP_align_fname\n";
print calibrate_peak_shift_par_file "genome_table = $contig_table_fname\n";
print calibrate_peak_shift_par_file "bin_align_path = $bin_align_ChIP_path\n";
print calibrate_peak_shift_par_file "regions_file = $quick_window_scan_output_fname.sorted\n";
print calibrate_peak_shift_par_file "output_file = $calibrate_peak_shift_output_fname\n";
print calibrate_peak_shift_par_file "kde_bandwidth = $KDE_bandwidth\n";
print calibrate_peak_shift_par_file "region_width = $PSC_region_width\n";
print calibrate_peak_shift_par_file "top_regions_passed = $top_regions_passed\n";
print calibrate_peak_shift_par_file "dist_threshold = $dist_threshold\n";
print calibrate_peak_shift_par_file "correlation_lower_threshold = $correlation_lower_threshold\n";
print calibrate_peak_shift_par_file "peak_shift_lower_threshold = $peak_shift_lower_threshold\n";
print calibrate_peak_shift_par_file "peak_shift_upper_threshold = $peak_shift_upper_threshold\n";
print calibrate_peak_shift_par_file "os_upper_threshold = $os_upper_threshold\n";
print calibrate_peak_shift_par_file "estimation_method = $peak_shift_estimation_method\n";

if($alternative_ps ne "undef"){
    print calibrate_peak_shift_par_file "alternative_ps = $alternative_ps\n";
}
close calibrate_peak_shift_par_file;

#print "Saved Calibrate Peak Shift batch parameter file.\n";

## /writing Calibrate Peak Shift parameters

## writing Generate Peak Profile parameters

#my $Profiler_calc_window = $KDE_bandwidth * 10;

#my $scores_path = $analysis_path . "/scores";
#
#if(-d $scores_path){
#    print "Warning: score directory $scores_path already exists\n";
#}
#else{
#    mkdir($scores_path);
#}

#my $ChIP_scores_path = $scores_path . "/ChIP/";
#my $pseudo_ChIP_scores_path = $scores_path . "/pseudo_ChIP/";
#my $background_scores_path = $scores_path . "/background/";

if(-d $ChIP_scores_path){ 
#    print "Warning: ChIP scores directory $ChIP_scores_path already exists\n"; 
#    print "When the Profiler Module will run the files will be overwritten\n";
}
else{ mkdir($ChIP_scores_path); }


if( $analysis_with_FDR eq "true"){
    if(-d $pseudo_ChIP_scores_path){ 
#	print "Warning: pseudo_ChIP scores directory $pseudo_ChIP_scores_path already exists\n"; 
#	print "When the Profiler Module will run the files will be overwritten\n";
    }
    else{ mkdir($pseudo_ChIP_scores_path); }
}

if(-d $background_scores_path){ 
#    print "Warning: backgrounds scores directory $background_scores_path already exists\n"; 
#    print "When the Profiler Module will run the files will be overwritten\n";
}
else{ mkdir($background_scores_path); }

my $ChIP_generate_profile_par_fname = $parameters_path . "/generate_profile.ChIP.batch.pars";
#my $ChIP_generate_profile_output_fname = $parameters_path . "/generate_profile.ChIP.out";

my $pseudo_ChIP_generate_profile_par_fname = $parameters_path . "/generate_profile.pseudo_ChIP.batch.pars";
#my $pseudo_ChIP_generate_profile_output_fname = $parameters_path . "/generate_profile.pseudo_ChIP.out";

my $background_generate_profile_par_fname = $parameters_path . "/generate_profile.background.batch.pars";
#my $background_generate_profile_output_fname = $analysis_path . "/generate_profile.background.out";

## writing ChIP Profiler param file

unlink($ChIP_generate_profile_par_fname);

open ChIP_generate_profile_ChIP_par_file, "> $ChIP_generate_profile_par_fname"
    or die "Failed to open $ChIP_generate_profile_par_fname";

print ChIP_generate_profile_ChIP_par_file "genome_table = $contig_table_fname\n";
print ChIP_generate_profile_ChIP_par_file "bin_align_path = $bin_align_ChIP_path\n";
print ChIP_generate_profile_ChIP_par_file "output_score_path = $ChIP_scores_path\n";
print ChIP_generate_profile_ChIP_par_file "peak_shift_source_file = $calibrate_peak_shift_output_fname\n";
print ChIP_generate_profile_ChIP_par_file "kde_bandwidth = $KDE_bandwidth\n";
print ChIP_generate_profile_ChIP_par_file "threads = $threads\n";

close ChIP_generate_profile_ChIP_par_file;

#print "Saved Profiler batch parameter file for the ChIP data\n";

## /writing ChIP Profiler param file

## writing background Profiler param file

unlink($background_generate_profile_par_fname);
if($background_tags > 0){
    
    open background_generate_profile_ChIP_par_file, "> $background_generate_profile_par_fname"
	or die "Failed to open $background_generate_profile_par_fname";

    print background_generate_profile_ChIP_par_file "genome_table = $contig_table_fname\n";
    print background_generate_profile_ChIP_par_file "bin_align_path = $bin_align_background_path\n";
    print background_generate_profile_ChIP_par_file "output_score_path = $background_scores_path\n";
    print background_generate_profile_ChIP_par_file "peak_shift_source_file = $calibrate_peak_shift_output_fname\n";
    print background_generate_profile_ChIP_par_file "kde_bandwidth = $KDE_bandwidth\n";
    print background_generate_profile_ChIP_par_file "threads = $threads\n";
    
    close background_generate_profile_ChIP_par_file;
}



## /writing background Profiler param file

## writing pseudo_ChIP Profiler param file

if($analysis_with_FDR eq "true"){

    unlink($pseudo_ChIP_generate_profile_par_fname);

    open pseudo_ChIP_generate_profile_ChIP_par_file, "> $pseudo_ChIP_generate_profile_par_fname"
	or die "Failed to open $pseudo_ChIP_generate_profile_par_fname";
 
    print pseudo_ChIP_generate_profile_ChIP_par_file "genome_table = $contig_table_fname\n";
    print pseudo_ChIP_generate_profile_ChIP_par_file "bin_align_path = $bin_align_pseudo_ChIP_path\n";
    print pseudo_ChIP_generate_profile_ChIP_par_file "output_score_path = $pseudo_ChIP_scores_path\n";
    print pseudo_ChIP_generate_profile_ChIP_par_file "peak_shift_source_file = $calibrate_peak_shift_output_fname\n";
    print pseudo_ChIP_generate_profile_ChIP_par_file "kde_bandwidth = $KDE_bandwidth\n";
    print pseudo_ChIP_generate_profile_ChIP_par_file "threads = $threads\n";

    close pseudo_ChIP_generate_profile_ChIP_par_file;

#    print "Saved Profiler batch parameter file for the pseudo_ChIP data\n";
}

## /writing pseudo_ChIP Profiler param file

## writing report_bad_control_regions param file

my $report_bad_control_regions_par_fname = $parameters_path . "/report_bad_control_regions.batch.pars";
my $bad_control_regions_fname = $module_outputs_path . "/bad_control_regions.txt";

if(-e $report_bad_control_regions_par_fname){
    unlink($report_bad_control_regions_par_fname);
}

my $report_bad_control_regions = "undefined";
if($background_tags > 0){
    $report_bad_control_regions = "true";
}
else{
    $report_bad_control_regions = "false";
}

open report_bad_control_regions_par_file, ">$report_bad_control_regions_par_fname" or
    die "Failed to open $report_bad_control_regions_par_fname\n";

if($report_bad_control_regions eq "true"){
    
    print report_bad_control_regions_par_file "genome_table = $contig_table_fname\n";
    print report_bad_control_regions_par_file "background_score_path = $background_scores_path\n";
    print report_bad_control_regions_par_file "region_seeding_threshold = 50\n";
    print report_bad_control_regions_par_file "region_extension_threshold = 5\n";
    print report_bad_control_regions_par_file "score_normalizer = $background_basal_level\n";
    print report_bad_control_regions_par_file "output_file = $bad_control_regions_fname\n"
}
else{
    print report_bad_control_regions_par_file "skip\n";
}
close report_bad_control_regions_par_file;


## /writing report_bad_control_regions param file


## writing ChIP Peak Caller param file

my $ChIP_peak_caller_par_fname = $parameters_path . "/peak_caller.ChIP.batch.pars";
my $ChIP_peak_calls_fname = $calls_path . "/peak_caller.ChIP.out";

my $pseudo_ChIP_peak_caller_par_fname = $parameters_path . "/peak_caller.pseudo_ChIP.batch.pars";
my $pseudo_ChIP_peak_calls_fname = $calls_path . "/peak_caller.pseudo_ChIP.out";

unlink($ChIP_peak_caller_par_fname);

open ChIP_peak_caller_par_file, "> $ChIP_peak_caller_par_fname"
    or die "Failed to open $ChIP_peak_caller_par_fname\n";

print ChIP_peak_caller_par_file "genome_table = $contig_table_fname\n";
print ChIP_peak_caller_par_file "chip_profile_path = $ChIP_scores_path\n";
if($background_tags > 0){
    print ChIP_peak_caller_par_file "background_profile_path = $background_scores_path\n";
}
else{
    print ChIP_peak_caller_par_file "background_profile_path = NA\n";
}
print ChIP_peak_caller_par_file "output_file = $ChIP_peak_calls_fname\n";
print ChIP_peak_caller_par_file "ChIP_threshold = $ChIP_threshold\n";
print ChIP_peak_caller_par_file "ChIP_extension_threshold = $ChIP_extension_threshold\n";

if($background_tags > 0){
    print ChIP_peak_caller_par_file "ChIP_to_background_ratio = $ChIP_to_background_ratio\n";
}
else{
    print ChIP_peak_caller_par_file "ChIP_to_background_ratio = NA\n";
    
}
#print ChIP_peak_caller_par_file "background_threshold = $background_threshold\n";
#print ChIP_peak_caller_par_file "rescue_ratio = $rr\n";
print ChIP_peak_caller_par_file "local_maximum_radius = 10\n";
print ChIP_peak_caller_par_file "dip_fraction = $dip_fraction\n";
print ChIP_peak_caller_par_file "ChIP_reads = $ChIP_tags\n";
print ChIP_peak_caller_par_file "background_reads = $background_tags\n";
print ChIP_peak_caller_par_file "ChIP_basal_level = $ChIP_basal_level\n";

close ChIP_peak_caller_par_file;

 #print "Saved Peak Caller batch parameter file for the ChIP data\n";
## /writing ChIP Peak Caller param file

## writing pseudo_ChIP Peak Caller param file
if($analysis_with_FDR eq "true"){
    unlink($pseudo_ChIP_peak_caller_par_fname);
    
    open pseudo_ChIP_peak_caller_par_file, "> $pseudo_ChIP_peak_caller_par_fname"
	or die "Failed to open $pseudo_ChIP_peak_caller_par_fname\n";
 
    print pseudo_ChIP_peak_caller_par_file "genome_table = $contig_table_fname\n";
    print pseudo_ChIP_peak_caller_par_file "chip_profile_path = $pseudo_ChIP_scores_path\n";
    print pseudo_ChIP_peak_caller_par_file "background_profile_path = $background_scores_path\n";
    print pseudo_ChIP_peak_caller_par_file "output_file = $pseudo_ChIP_peak_calls_fname\n";
    print pseudo_ChIP_peak_caller_par_file "ChIP_threshold = $ChIP_threshold\n";
    print pseudo_ChIP_peak_caller_par_file "ChIP_extension_threshold = $ChIP_extension_threshold\n";
    print pseudo_ChIP_peak_caller_par_file "ChIP_to_background_ratio = $ChIP_to_background_ratio\n";
    #print pseudo_ChIP_peak_caller_par_file "background_threshold = $background_threshold\n";
    #print pseudo_ChIP_peak_caller_par_file "rescue_ratio = $rr\n";
    print pseudo_ChIP_peak_caller_par_file "local_maximum_radius = 10\n";
    print pseudo_ChIP_peak_caller_par_file "dip_fraction = $dip_fraction\n";
    
    #print pseudo_ChIP_peak_caller_par_file "ChIP_reads = $pseudo_ChIP_tags\n";
    #print pseudo_ChIP_peak_caller_par_file "background_reads = $background_tags\n";
    
    print pseudo_ChIP_peak_caller_par_file "ChIP_reads = $pseudo_ChIP_tags\n";
    print pseudo_ChIP_peak_caller_par_file "background_reads = $background_tags\n";
#    print pseudo_ChIP_peak_caller_par_file "ChIP_basal_level = $pseudo_ChIP_basal_level\n";
    print pseudo_ChIP_peak_caller_par_file "ChIP_basal_level = $ChIP_basal_level\n";
    
    
    close pseudo_ChIP_peak_caller_par_file;

    print "Saved Peak Caller batch parameter file for the pseudo_ChIP data.\n";
}

## /writing pseudo_ChIP Peak Caller param file

my $metrics_ChIP_par_fname = $parameters_path . "/metrics.ChIP.batch.pars";
my $metrics_ChIP_report_fname = $module_outputs_path . "/metrics.ChIP.report.txt";
my $regions_fname = $ChIP_peak_calls_fname;# . ".regions.no_metrics";
my $metrics_ChIP_out_fname = $ChIP_peak_calls_fname . ".with_metrics";


unlink($metrics_ChIP_par_fname);

open metrics_ChIP_par_file, "> $metrics_ChIP_par_fname" or die
    "Failed to open $metrics_ChIP_par_fname\n";

print metrics_ChIP_par_file "genome_table = $contig_table_fname\n";
print metrics_ChIP_par_file "ChIP_bin_align_path = $bin_align_ChIP_path\n";
print metrics_ChIP_par_file "background_bin_align_path = $bin_align_background_path\n";
print metrics_ChIP_par_file "QuEST_calls_file = $regions_fname\n";
print metrics_ChIP_par_file "output_file = $metrics_ChIP_out_fname\n";
print metrics_ChIP_par_file "region_width = $strand_correlation_region_width\n";
print metrics_ChIP_par_file "kde_bandwidth = $KDE_bandwidth\n";
print metrics_ChIP_par_file "dist_threshold = $dist_threshold\n";
print metrics_ChIP_par_file "ChIP_tags = $ChIP_tags\n";
print metrics_ChIP_par_file "background_tags = $background_tags\n";
print metrics_ChIP_par_file "ChIP_basal_level = $ChIP_basal_level\n";
print metrics_ChIP_par_file "background_basal_level = $background_basal_level\n";
print metrics_ChIP_par_file "ChIP_basal_tag_rate = $ChIP_basal_tag_rate\n";
print metrics_ChIP_par_file "background_basal_tag_rate = $background_basal_tag_rate\n";

print metrics_ChIP_par_file "report_file = $metrics_ChIP_report_fname\n";

close metrics_ChIP_par_file;

#print "Saved Strand Correlation batch parameter file for the ChIP data.\n";

my $filter_ChIP_regions_par_fname = $parameters_path . "/peak_filter.ChIP.batch.pars";
my $filter_ChIP_calls_report_fname = $module_outputs_path . "/peak_filter.ChIP.report.txt";
my $ChIP_peaks_filtered_fname_prefix = $ChIP_peak_calls_fname;

unlink($filter_ChIP_regions_par_fname);

open filter_ChIP_regions_par_file, "> $filter_ChIP_regions_par_fname" or die
    "Failed to open $filter_ChIP_regions_par_fname\n";


print filter_ChIP_regions_par_file "input_file = $metrics_ChIP_out_fname\n";
print filter_ChIP_regions_par_file "output_file_prefix = $ChIP_peaks_filtered_fname_prefix\n";
#print filter_ChIP_regions_par_file "log10_qv_threshold = $filter_ChIP_regions_log10_qv_threshold\n";

if($ps_lower_threshold ne $undef){
    print filter_ChIP_regions_par_file "ps_lower_threshold = $ps_lower_threshold\n";
}
if($ps_upper_threshold ne $undef){
    print filter_ChIP_regions_par_file "ps_upper_threshold = $ps_upper_threshold\n";
}
if($cor_lower_threshold ne $undef){
    print filter_ChIP_regions_par_file "cor_lower_threshold = $cor_lower_threshold\n";
}
print filter_ChIP_regions_par_file "filter_out_empty_regions = $filter_out_empty_regions\n";
print filter_ChIP_regions_par_file "report_peaks = $report_peaks\n";
print filter_ChIP_regions_par_file "bad_control_regions_file = $bad_control_regions_fname\n";
print filter_ChIP_regions_par_file "report_file = $filter_ChIP_calls_report_fname\n";

close filter_ChIP_regions_par_file;

## writing output calls as a bed file parameter file

my $ChIP_filtered_bed_par_fname = $parameters_path . "/calls_tracks.ChIP.filtered.batch.pars";
my $ChIP_filtered_calls_fname = $ChIP_peaks_filtered_fname_prefix . ".accepted";
my $ChIP_filtered_calls_bed_fname = $tracks_path . "/ChIP_calls.filtered.bed";
my $ChIP_filtered_calls_track_name = $ChIP_name . "_filtered";
open ChIP_filtered_bed_par_file, "> $ChIP_filtered_bed_par_fname" or die
    "Failed to open $ChIP_filtered_bed_par_fname\n";

print ChIP_filtered_bed_par_file "regions_file = $ChIP_filtered_calls_fname\n";
print ChIP_filtered_bed_par_file "output_file = $ChIP_filtered_calls_bed_fname\n";
print ChIP_filtered_bed_par_file "track_name = $ChIP_filtered_calls_track_name\n";
print ChIP_filtered_bed_par_file "track_priority = $ChIP_calls_track_priority\n";

close ChIP_filtered_bed_par_file;

## /writing output calls as a bed file parameter file

## writing bed graph files for ChIP

my $bed_graph_tracks_ChIP_par_fname = $parameters_path . "/bedgraph_ChIP.batch.pars";
my $bed_graph_tracks_ChIP_track_name = $ChIP_name . "_tag_start_counts";
my $bed_graph_ChIP_output_fname = $bed_graph_tracks_path . "/" . $ChIP_name . ".bedGraph";
my $bed_graph_tracks_by_chr_ChIP_prefix = $bed_graph_tracks_by_chr_ChIP_path . "/" . $ChIP_name;

open bed_graph_tracks_by_chr_ChIP_file, "> $bed_graph_tracks_ChIP_par_fname"
    or die "Failed to open $bed_graph_tracks_ChIP_par_fname\n";


print bed_graph_tracks_by_chr_ChIP_file "genome_table = $contig_table_fname\n";
print bed_graph_tracks_by_chr_ChIP_file "track_name = $bed_graph_tracks_ChIP_track_name\n";
print bed_graph_tracks_by_chr_ChIP_file "output_file = $bed_graph_ChIP_output_fname\n";
print bed_graph_tracks_by_chr_ChIP_file "track_priority = $ChIP_bedgraph_track_priority\n";
print bed_graph_tracks_by_chr_ChIP_file "count_threshold = 1\n";
print bed_graph_tracks_by_chr_ChIP_file "bin_align_prefix = $bin_align_ChIP_path/\n";
print bed_graph_tracks_by_chr_ChIP_file "by_chr_prefix = $bed_graph_tracks_by_chr_ChIP_prefix\n";

## /writing bed graph files for ChIP

## writing track generation param file
## ChIP
my $ChIP_reads_bed_par_fname = $parameters_path . "/ChIP.bed.batch.pars";
my $ChIP_reads_bed_fname = $data_bed_tracks_path . "/$ChIP_name.bed";
my $ChIP_reads_bed_by_chr_prefix = $ChIP_bed_tracks_by_chr_path . "/$ChIP_name";
#my $ChIP_reads_bed_track_priority = $ChIP_bed_track_priority;

open ChIP_reads_bed_par_file, "> $ChIP_reads_bed_par_fname"
	or die "Failed to open $ChIP_reads_bed_par_fname\n";

print ChIP_reads_bed_par_file "genome_table = $contig_table_fname\n";
print ChIP_reads_bed_par_file "input_file = $ChIP_collapsed_align_fname\n";
print ChIP_reads_bed_par_file "output_file = $ChIP_reads_bed_fname\n";
print ChIP_reads_bed_par_file "track_name = $ChIP_name";
print ChIP_reads_bed_par_file "_tags\n";
print ChIP_reads_bed_par_file "output_by_chr_prefix = $ChIP_reads_bed_by_chr_prefix\n";
print ChIP_reads_bed_par_file "track_priority = $ChIP_reads_track_priority\n";
print ChIP_reads_bed_par_file "gz = yes\n";

close ChIP_reads_bed_par_file;

my $ChIP_normalized_wig_par_fname = $parameters_path . "/ChIP_normalized.wig.batch.pars";
my $ChIP_normalized_wig_fname = $wig_tracks_path . "/" . $ChIP_name . "_normalized.profile.wig";

my $ChIP_wig_normalizer = $ChIP_basal_level;
my $ChIP_wig_threshold =  10; 

open ChIP_normalized_wig_par_file, "> $ChIP_normalized_wig_par_fname"
	or die "Failed to open $ChIP_normalized_wig_par_fname\n";

print ChIP_normalized_wig_par_file "genome_table = $contig_table_fname\n";
print ChIP_normalized_wig_par_file "profile_path = $ChIP_scores_path\n";
print ChIP_normalized_wig_par_file "output_file = $ChIP_normalized_wig_fname\n";
print ChIP_normalized_wig_par_file "output_by_chr_path = $wig_tracks_ChIP_normalized_by_chr_path\n";
print ChIP_normalized_wig_par_file "profile_threshold = $ChIP_wig_threshold\n";
print ChIP_normalized_wig_par_file "normalizer = $ChIP_wig_normalizer\n";
print ChIP_normalized_wig_par_file "track_name = $ChIP_name";
print ChIP_normalized_wig_par_file "_normalized\n";
print ChIP_normalized_wig_par_file "track_priority = $ChIP_normalized_wig_track_priority\n";
print ChIP_normalized_wig_par_file "gz = yes\n";
			  
close ChIP_normalized_wig_par_file;

my $ChIP_unnormalized_wig_par_fname = $parameters_path . "/ChIP_unnormalized.wig.batch.pars";
my $ChIP_unnormalized_wig_fname = $wig_tracks_path . "/" . $ChIP_name . "_unnormalized.profile.wig";
my $ChIP_unnormalized_wig_threshold = 5.0 * ( 1.0/sqrt(2.0 * 3.14159 * $KDE_bandwidth * $KDE_bandwidth) );

open ChIP_unnormalized_wig_par_file, "> $ChIP_unnormalized_wig_par_fname"
	or die "Failed to open $ChIP_unnormalized_wig_par_fname\n";

print ChIP_unnormalized_wig_par_file "genome_table = $contig_table_fname\n";
print ChIP_unnormalized_wig_par_file "profile_path = $ChIP_scores_path\n";
print ChIP_unnormalized_wig_par_file "output_file = $ChIP_unnormalized_wig_fname\n";
print ChIP_unnormalized_wig_par_file "output_by_chr_path = $wig_tracks_ChIP_unnormalized_by_chr_path\n";
print ChIP_unnormalized_wig_par_file "profile_threshold = $ChIP_unnormalized_wig_threshold\n";
print ChIP_unnormalized_wig_par_file "normalizer = 1\n";
print ChIP_unnormalized_wig_par_file "track_name = $ChIP_name";
print ChIP_unnormalized_wig_par_file "_unnormalized\n";
print ChIP_unnormalized_wig_par_file "track_priority = $ChIP_unnormalized_wig_track_priority\n";
print ChIP_unnormalized_wig_par_file "gz = yes\n";
			  
close ChIP_unnormalized_wig_par_file;
## /ChIP

## background 
my $RX_noIP_reads_bed_par_fname = $parameters_path . "/RX_noIP.bed.batch.pars";
my $RX_noIP_reads_bed_fname = $data_bed_tracks_path . "/RX_noIP.bed";
my $RX_noIP_reads_bed_by_chr_prefix = $RX_noIP_bed_tracks_by_chr_path . "/RX_noIP";

open RX_noIP_reads_bed_par_file, "> $RX_noIP_reads_bed_par_fname"
	or die "Failed to open $RX_noIP_reads_bed_par_fname\n";

print RX_noIP_reads_bed_par_file "genome_table = $contig_table_fname\n";
print RX_noIP_reads_bed_par_file "input_file = $RX_noIP_collapsed_align_fname\n";
print RX_noIP_reads_bed_par_file "output_file = $RX_noIP_reads_bed_fname\n";
print RX_noIP_reads_bed_par_file "track_name = RX_noIP";
print RX_noIP_reads_bed_par_file "_tags\n";
print RX_noIP_reads_bed_par_file "output_by_chr_prefix = $RX_noIP_reads_bed_by_chr_prefix\n";
print RX_noIP_reads_bed_par_file "track_priority = $RX_noIP_reads_track_priority\n";
print RX_noIP_reads_bed_par_file "gz = yes\n";

close RX_noIP_reads_bed_par_file;

my $background_normalized_wig_par_fname = "NA";
my $background_normalized_wig_fname = "NA";
my $background_unnormalized_wig_par_fname = "NA";
my $background_unnormalized_wig_fname = "NA";

if($background_tags > 0){
    $background_normalized_wig_par_fname = $parameters_path . "/background_normalized.wig.batch.pars";
    $background_normalized_wig_fname = $wig_tracks_path . "/background_normalized.profile.wig";

    my $background_wig_normalizer = $background_basal_level;
    my $background_wig_threshold = 10;

    open background_normalized_wig_par_file, "> $background_normalized_wig_par_fname"
	or die "Failed to open $background_normalized_wig_par_fname\n";
    
    print background_normalized_wig_par_file "genome_table = $contig_table_fname\n";
    print background_normalized_wig_par_file "profile_path = $background_scores_path\n";
    print background_normalized_wig_par_file "output_file = $background_normalized_wig_fname\n";
    print background_normalized_wig_par_file "output_by_chr_path = $wig_tracks_background_normalized_by_chr_path\n";
    print background_normalized_wig_par_file "profile_threshold = $background_wig_threshold\n";
    print background_normalized_wig_par_file "normalizer = $background_wig_normalizer\n";
    print background_normalized_wig_par_file "track_name = $ChIP_name";
    print background_normalized_wig_par_file "_background_normalized\n";
    print background_normalized_wig_par_file "track_priority = $RX_noIP_normalized_wig_track_priority\n"; 
    print background_normalized_wig_par_file "gz = yes\n";
    
    close background_normalized_wig_par_file;
    
    $background_unnormalized_wig_par_fname = $parameters_path . "/background_unnormalized.wig.batch.pars";
    $background_unnormalized_wig_fname = $wig_tracks_path . "/background_unnormalized.profile.wig";

    my $background_unnormalized_wig_threshold = 5.0 * ( 1.0/sqrt(2.0 * 3.14159 * $KDE_bandwidth * $KDE_bandwidth) );
    my $background_unnormalized_wig_track_priority = 22;
    
    open background_unnormalized_wig_par_file, "> $background_unnormalized_wig_par_fname"
	or die "Failed to open $background_unnormalized_wig_par_fname\n";
    
    print background_unnormalized_wig_par_file "genome_table = $contig_table_fname\n";
    print background_unnormalized_wig_par_file "profile_path = $background_scores_path\n";
    print background_unnormalized_wig_par_file "output_file = $background_unnormalized_wig_fname\n";
    print background_unnormalized_wig_par_file "output_by_chr_path = $wig_tracks_background_unnormalized_by_chr_path\n";
    print background_unnormalized_wig_par_file "profile_threshold = $background_unnormalized_wig_threshold\n";
    print background_unnormalized_wig_par_file "normalizer = 1\n";
    print background_unnormalized_wig_par_file "track_name = $ChIP_name";
    print background_unnormalized_wig_par_file "_background_unnormalized\n";
    print background_unnormalized_wig_par_file "track_priority = $RX_noIP_unnormalized_wig_track_priority\n";
    print background_unnormalized_wig_par_file "gz = yes\n";
    
    close background_unnormalized_wig_par_file;
}
## /background
    
## /writing track generation param file

## writing a superscript param file

my $QuEST_param_fname = $parameters_path . "/QuEST.batch.pars";
my $QuEST_log_fname = $logs_path . "/QuEST.log";
my $QuEST_output_fname = $module_outputs_path . "/QuEST.out";
my $QuEST_parametrization_report_fname = $module_outputs_path . "/QuEST_parametrization.report.txt";

if(-e $QuEST_parametrization_report_fname){
    unlink($QuEST_parametrization_report_fname);
}

open QuEST_parametrization_report_file, "> $QuEST_parametrization_report_fname" or die "Failed to open $QuEST_parametrization_report_fname\n";

print QuEST_parametrization_report_file "-----------------------------\n";
print QuEST_parametrization_report_file "Genome statistics\n";
print QuEST_parametrization_report_file "\n";
print QuEST_parametrization_report_file "Genome size: $genome_size\n";
print QuEST_parametrization_report_file "Mappable genome fraction: $mappable_genome_fraction\n";
printf QuEST_parametrization_report_file "chromosomes: %d\n", scalar(@contig_names);
print QuEST_parametrization_report_file "\n";


print QuEST_parametrization_report_file "-----------------------------\n";
print QuEST_parametrization_report_file "Input data statistics\n";
print QuEST_parametrization_report_file "\n";

if($RX_noIP_align_fname eq $missing){
    print QuEST_parametrization_report_file "control data: missing\n";
}

print QuEST_parametrization_report_file "ChIP reads before collapsing: $ChIP_uncollapsed_tags\n";
print QuEST_parametrization_report_file "RX_noIP reads before collapsing: $RX_noIP_uncollapsed_tags\n";
print QuEST_parametrization_report_file "\n";

print QuEST_parametrization_report_file "Effective ChIP reads: $ChIP_tags\n";
print QuEST_parametrization_report_file "Effective RX_noIP reads: $RX_noIP_tags\n";
print QuEST_parametrization_report_file "Effective background reads: $background_tags\n";
print QuEST_parametrization_report_file "\n";

printf QuEST_parametrization_report_file "Eliminated %d ", 100 * ($ChIP_uncollapsed_tags - $ChIP_tags) / $ChIP_uncollapsed_tags;
print QuEST_parametrization_report_file " % of ChIP data due to stacking\n";
if($RX_noIP_tags > 0){
    printf QuEST_parametrization_report_file "Eliminated %d ", 100 * ($RX_noIP_uncollapsed_tags - $RX_noIP_tags ) / $RX_noIP_uncollapsed_tags;
    print QuEST_parametrization_report_file " % of RX_noIP data due to stacking\n";
}
printf QuEST_parametrization_report_file "\n";


print QuEST_parametrization_report_file "Analysis with FDR: $analysis_with_FDR\n";
if($analysis_with_FDR eq "true"){
    print QuEST_parametrization_report_file "Effective pseudo ChIP reads: $pseudo_ChIP_tags\n";
}
print QuEST_parametrization_report_file "\n";

print QuEST_parametrization_report_file "-----------------------------\n";
print QuEST_parametrization_report_file "Peak calling parameters\n";
print QuEST_parametrization_report_file "\n";

print QuEST_parametrization_report_file "KDE bandwidth: $KDE_bandwidth\n";
print QuEST_parametrization_report_file "\n";

print QuEST_parametrization_report_file "ChIP seeding fold enrichment: $ChIP_seeding_enrichment\n";
print QuEST_parametrization_report_file "ChIP CDP seeding threshold: $ChIP_threshold\n";
print QuEST_parametrization_report_file "\n";

print QuEST_parametrization_report_file "ChIP extension fold enrichment: $ChIP_extension_enrichment\n";
print QuEST_parametrization_report_file "ChIP CDP extension threshold: $ChIP_extension_threshold\n";
print QuEST_parametrization_report_file "\n";

if($background_tags > 0){
    print QuEST_parametrization_report_file "ChIP-to-background fold enrichment: $ChIP_to_background_enrichment\n";
    print QuEST_parametrization_report_file "ChIP-to-background CDP ratio: $ChIP_to_background_ratio\n";
    print QuEST_parametrization_report_file "\n";
}



$ChIP_threshold = $ChIP_basal_level * $ChIP_seeding_enrichment;
$ChIP_extension_threshold = $ChIP_basal_level * $ChIP_extension_enrichment;
if($background_tags > 0){
    $ChIP_to_background_ratio = ($ChIP_tags / $background_tags) * $ChIP_to_background_enrichment;
}

close QuEST_parametrization_report_file;


if( -e $QuEST_param_fname ){ unlink($QuEST_param_fname); }

open QuEST_param_file, "> $QuEST_param_fname" 
    or die "Failed to open $QuEST_param_fname\n";

print QuEST_param_file "analysis_with_FDR = $analysis_with_FDR\n";

#print QuEST_param_file "align_2_bin_ChIP_param_file = $bin_align_ChIP_par_fname\n";
#print QuEST_param_file "align_2_bin_background_param_file = $bin_align_background_par_fname\n";
#print QuEST_param_file "align_2_bin_RX_noIP_param_file = $bin_align_RX_noIP_par_fname\n";

print QuEST_param_file "quick_window_scan_param_file = $quick_window_scan_par_fname\n";
print QuEST_param_file "calibrate_peak_shift_param_file = $calibrate_peak_shift_par_fname\n";
print QuEST_param_file "ChIP_generate_profile_param_file = $ChIP_generate_profile_par_fname\n";
if($background_tags > 0){
    print QuEST_param_file "background_generate_profile_param_file = $background_generate_profile_par_fname\n";
}

print QuEST_param_file "report_bad_control_regions_param_file = $report_bad_control_regions_par_fname\n";

print QuEST_param_file "ChIP_peak_caller_param_file = $ChIP_peak_caller_par_fname\n";

if($analysis_with_FDR eq "true"){
    print QuEST_param_file "pseudo_ChIP_generate_profile_param_file = $pseudo_ChIP_generate_profile_par_fname\n";
    print QuEST_param_file "pseudo_ChIP_peak_caller_param_file = $pseudo_ChIP_peak_caller_par_fname\n";
#    print QuEST_param_file "align_2_bin_pseudo_ChIP_param_file = $bin_align_pseudo_ChIP_par_fname\n";
}

print QuEST_param_file "bed_graph_tracks_ChIP_param_file = $bed_graph_tracks_ChIP_par_fname\n";
print QuEST_param_file "metrics_ChIP_param_file = $metrics_ChIP_par_fname\n";
print QuEST_param_file "filter_ChIP_calls_param_file = $filter_ChIP_regions_par_fname\n";
print QuEST_param_file "ChIP_reads_bed_param_file = $ChIP_reads_bed_par_fname\n";
print QuEST_param_file "ChIP_normalized_wig_param_file = $ChIP_normalized_wig_par_fname\n";
print QuEST_param_file "ChIP_unnormalized_wig_param_file = $ChIP_unnormalized_wig_par_fname\n";
print QuEST_param_file "RX_noIP_reads_bed_param_file = $RX_noIP_reads_bed_par_fname\n";
print QuEST_param_file "background_normalized_wig_param_file = $background_normalized_wig_par_fname\n";
print QuEST_param_file "background_unnormalized_wig_param_file = $background_unnormalized_wig_par_fname\n";
print QuEST_param_file "filtered_calls_tracks_param_file = $ChIP_filtered_bed_par_fname\n";
print QuEST_param_file "scores_dir = $scores_path\n";


print QuEST_param_file "log_file = $QuEST_log_fname\n";
print QuEST_param_file "output_file = $QuEST_output_fname\n";

print QuEST_param_file "QuEST_parametrization_report_file = $QuEST_parametrization_report_fname\n";
print QuEST_param_file "metrics_report_file = $metrics_ChIP_report_fname\n";
print QuEST_param_file "filter_report_file = $filter_ChIP_calls_report_fname\n";

print QuEST_param_file "\n";
print QuEST_param_file "# ** Pipeline steps, use \# to comment out **\n";
print QuEST_param_file "\n";

#print QuEST_param_file "align_2_bin: ChIP\n";
#print QuEST_param_file "align_2_bin: background\n";
#
#print QuEST_param_file "align_2_bin: RX_noIP\n";
#if($analysis_with_FDR eq "true"){
#    print QuEST_param_file "align_2_bin: pseudo_ChIP\n";
#}

print QuEST_param_file "quick_window_scan: ChIP\n";
print QuEST_param_file "calibrate_peak_shift: ChIP\n";
print QuEST_param_file "generate_profile: ChIP\n";

if($background_tags > 0){
    print QuEST_param_file "generate_profile: background\n";
}
print QuEST_param_file "report_bad_control_regions: $report_bad_control_regions\n";

if($analysis_with_FDR eq "true"){
    print QuEST_param_file "generate_profile: pseudo_ChIP\n";
}
print QuEST_param_file "peak_caller: ChIP vs background\n";
if($analysis_with_FDR eq "true" || $analysis_with_FDR eq "simulation"){
    print QuEST_param_file "peak_caller: pseudo_ChIP vs background\n";
}

print QuEST_param_file "metrics_ChIP: true\n";

#if($background_tags > 0){
#    print QuEST_param_file "q-values: false\n";
#    print QuEST_param_file "q-values: true\n";
#}
#else{
#    print QuEST_param_file "q-values: false\n";
#}

print QuEST_param_file "bed_graph_tracks_ChIP: true\n";

print QuEST_param_file "bed_tracks_ChIP: true\n";
if($RX_noIP_tags > 0){
    print QuEST_param_file "bed_tracks_RX_noIP: true\n";
}
else{
    print QuEST_param_file "bed_tracks_RX_noIP: false\n";
}

print QuEST_param_file "wig_tracks_ChIP: true\n";
if($background_tags > 0){
    print QuEST_param_file "wig_tracks_background: true\n";
}
#print QuEST_param_file "wig_tracks_RX_noIP: true\n";
print QuEST_param_file "filter_ChIP_calls: true\n";
print QuEST_param_file "calls_tracks: true\n";
    

close QuEST_param_file;

## /writing a superscript param file

## /writing parameter files
close cur_log_file;

print "\n";

$answer = "";

while($answer ne "y" && $answer ne "n"){
    print "Would you like to run QuEST analysis now? (y/n): ";
    $answer = <STDIN>;
    chomp($answer);
    
}

if($answer eq "y"){
    my $run_QuEST_system_command = "$run_QuEST_script -ap $analysis_path";
    system($run_QuEST_system_command);
}
