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
    run_profile_2_wig_with_param_file.pl
	
    This program is a wrapper that runs the peak_caller
    
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

$exec_path = $exec_path."/profile_2_wig";
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

my $profile_path =           "** missing **";
my $output_fname =           "** missing **";
my $output_by_chr_path =     "** missing **";
my $profile_threshold =      "** missing **";
my $genome_table_fname =     "** missing **";
my $normalizer =             "** missing **";
my $track_name =             "** missing **";
my $track_priority =         "** missing **";

my $gz_flag = "no";

my $track_color = "";

my $red_value = int(rand(256));
my $green_value = int(rand(256));
my $blue_value = int(rand(256));

$track_color = $red_value . "," . $green_value . "," . $blue_value;

for(my $i=0; $i<scalar(@params); $i++){
    my $cur_param = $params[$i];
    chomp($cur_param);
    
    my @cur_par_fields = split(/ /, $cur_param);
    if(scalar(@cur_par_fields >= 2)){
	my $cur_par_name = $cur_par_fields[0];
	if ($cur_par_name eq "profile_path"){
	    $profile_path = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "output_file"){
	    $output_fname = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "output_by_chr_path"){
	    $output_by_chr_path = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "output_file"){
	    $output_fname = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "profile_threshold"){
	    $profile_threshold = $cur_par_fields[2];	
	}
	elsif($cur_par_name eq "genome_table"){
	    $genome_table_fname = $cur_par_fields[2];
	}
	elsif ($cur_par_name eq "normalizer"){
	    $normalizer = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "track_name"){
	    $track_name = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "track_priority"){
	    $track_priority = $cur_par_fields[2];	
	}
	elsif ($cur_par_name eq "gz"){
	    $gz_flag = $cur_par_fields[2];	
	}
	else{
	    if($params[$i] ne ""){
		print "Warning: unrecognized parameter: $cur_par_name";
		exit(4);
	    }
	}
    }
}

my $output_fname_gz = $output_fname . ".gz";

print "Read the following parameters: \n\n";

print "profile_path:                      $profile_path\n";
print "output_file:                       $output_fname\n";
print "output_by_chr_path:                $output_by_chr_path\n";
print "profile_threshold:                 $profile_threshold\n";
print "genome_table_file:                 $genome_table_fname\n";
print "normalizer:                        $normalizer\n";
print "track_name:                        $track_name\n";
print "track_priority:                    $track_priority\n";
print "gz:                                $gz_flag\n";

print "\n";

my @files_to_trash;
my $files_to_trash_counter = 0;

my $optional_param_string = "";
#if($local_maximum_radius ne ""){
#    $optional_param_string = $optional_param_string . " local_maximum_radius=$local_maximum_radius";
#}
#if($dip_fraction ne ""){
#    $optional_param_string = $optional_param_string . " dip_fraction=$dip_fraction";
#}
#

## read genome table

if( ! -e $genome_table_fname){
    print "Error in run peak caller script: couldn't find genome table $genome_table_fname.\n";
    print "Aborting.\n";
    exit(0);
}

open genome_table_file, "< $genome_table_fname" || die "$genome_table_fname: $\n";
my @genome_table = <genome_table_file>;
close genome_table_file;

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

if(-e $output_fname){
    unlink($output_fname);
    system("touch $output_fname");
}

open output_file, "> $output_fname" || die "$output_fname: $\n";
print output_file "track type=wiggle_0 name=\"$track_name\" description=\"$track_name\" visibility=full color=$track_color priority=$track_priority maxHeightPixels=50:50:11\n";
close output_file;


for(my $i=0; $i<scalar(@contigs); $i++){
    my $cur_contig = $contigs[$i];
    my $cur_output_fname = "$profile_path/$cur_contig.wig.tmp";
    my $system_command = 
	$exec_path . 
	" score_file=$profile_path/$cur_contig.score" . 
	" output_file=$cur_output_fname" .
	" profile_threshold=$profile_threshold" .
	" contig_id=$cur_contig" .
	" normalizer=$normalizer" . 
	" track_name=$track_name" .
	" track_color=$track_color";
#	$optional_param_string;
    print "\nsystem_command: $system_command\n";
    
    $files_to_trash[$files_to_trash_counter] = $cur_output_fname;
    $files_to_trash_counter++;

    my $error_code = 0;
    $error_code = system("$system_command");
    if($error_code != 0){
	print "run_profile_2_wig_with_param_file.pl: Caught an error message (code $error_code), passing the error message to the top script.\n";
	exit(2);
    }

    my $cur_chr_wig_fname = $output_by_chr_path . "/$cur_contig.wig";
    my $cur_chr_wig_fname_gz = $cur_chr_wig_fname . ".gz";
    my $cur_chr_track_name = $track_name . "_$cur_contig";
    if(-e $cur_chr_wig_fname){
	unlink($cur_chr_wig_fname);
    }
    if(-e $cur_chr_wig_fname_gz){
	if($gz_flag eq "yes"){
	    unlink($cur_chr_wig_fname_gz);
	}
    }

    
    
    open cur_chr_output_file, "> $cur_chr_wig_fname" || die "$cur_chr_wig_fname: $\n";
    print cur_chr_output_file "track type=wiggle_0 name=\"$cur_chr_track_name\" description=\"$cur_chr_track_name\" visibility=full color=$track_color priority=$track_priority maxHeightPixels=50:50:11\n";
    close cur_chr_output_file;
    system("cat $cur_output_fname >> $cur_chr_wig_fname");
    if($gz_flag eq "yes"){
	print "\ngzip will run in the background. \n";
	print "Please allow enough time to finish before opening \n";
	print "$cur_chr_wig_fname_gz\n";
	system("gzip $cur_chr_wig_fname &");
    }
    system("cat $cur_output_fname >> $output_fname");
}

if(-e $output_fname_gz){
    unlink($output_fname_gz);
}

if($gz_flag eq "yes"){

    print "\ngzip will run in the background. \n";
    print "Please allow enough time to finish before opening \n";
    print "$output_fname.gz\n";
    print "\n";
    $|++;
    system("gzip $output_fname &");
    print "done\n";
}


for(my $i=0; $i<scalar(@files_to_trash); $i++){                                                            
    system("rm $files_to_trash[$i]");                                                                      
}
