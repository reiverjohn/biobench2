#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    bottleneck_metrick.pl

    -----------------------------
    mandatory parameters:
    
    -r regions_file
    -o output_file            
    -n track_name

    or 

    -p param_file

    -----------------------------
    optional parameters:
    -h (help)
    -tp track_priority;

};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments

my $undef = "undef";

my $regions_fname = "";
my $output_fname = "";
my $track_name = "";
my $track_priority = 10;

my $regions_color = "0,191,255";
my $peaks_color = "178,34,34";

my $param_fname = $undef;
## optional arguments

## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-r') {$regions_fname = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_fname = shift @ARGV;}
    elsif ( $this_arg eq '-n') {$track_name = shift @ARGV;}
    elsif ( $this_arg eq '-tp') {$track_priority = shift @ARGV;}
    elsif ( $this_arg eq '-p') {$param_fname = shift @ARGV;}
    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

print "param file: $param_fname\n";

if($param_fname ne $undef){
    open param_file, "< $param_fname" || die "$param_fname: $!\n";
    while( <param_file> ){
        chomp;
        my $cur_par_entry = $_;
        my @cur_par_entry_fields = split(/ /,$cur_par_entry);
        if(scalar(@cur_par_entry_fields) >= 3){
            my $cur_par_name = $cur_par_entry_fields[0];
            my $cur_par_value = $cur_par_entry_fields[2];
	    
            if($cur_par_name eq "regions_file"){
                $regions_fname = $cur_par_value;
            }
            if($cur_par_name eq "output_file"){
                $output_fname = $cur_par_value;
            }
	    if($cur_par_name eq "track_name"){
                $track_name = $cur_par_value;
            }
	    if($cur_par_name eq "track_priority"){
                $track_priority = $cur_par_value;
            }
	}
    }
}

if ( $regions_fname eq ""){
    die "you should specify peak file\n";
}
if( $output_fname eq ""){
    die "you should specify output file\n";
}

## print parameters

print "\n-----------------\n\n";
print "regions file: $regions_fname\n";
print "output file: $output_fname\n";
print "track name:  $track_name\n";
print "\n-----------------\n\n";


my $regions_name = $track_name . "_regions";

if(-e $output_fname){ unlink("$output_fname");}
open output_file, "> $output_fname" || die "$output_fname: $!\n";
print output_file "track name=$track_name description=$regions_name itemRgb=\"On\" priority=$track_priority  visibility=1\n";

open regions_file, "< $regions_fname" || die "$regions_fname: $!\n";

while(<regions_file>){
    chomp;
    my $cur_entry = $_;

    if( length($cur_entry) > 0 ){
	if(substr($cur_entry, 0, 1) ne "#"){
	    
	    my @cur_entry_fields = split(/ /, $cur_entry);
	    if(scalar(@cur_entry_fields) >= 9){
		
		if(substr($cur_entry,0,1) eq "R"){
		    my $cur_region_chr = $cur_entry_fields[1];
		    my $cur_region_range = $cur_entry_fields[2];
		    my @cur_region_range_fields = split(/-/, $cur_region_range);
		    if(scalar(@cur_region_range_fields) == 2){
			my $cur_region_start_coord = $cur_region_range_fields[0] + 1;
			my $cur_region_end_coord = $cur_region_range_fields[1] + 1;
			my $cur_region_id = $cur_entry_fields[0];
			my $cur_region_score = $cur_entry_fields[4];
			#print "$cur_region_chr $cur_region_start_coord $cur_region_end_coord\n";
			print output_file "$cur_region_chr $cur_region_start_coord $cur_region_end_coord $cur_region_id $cur_region_score + $cur_region_start_coord $cur_region_end_coord $regions_color\n";
		    }
		}		
		elsif(substr($cur_entry,0,1) eq "P"){
		    my $cur_peak_id = $cur_entry_fields[0];
		    
		    my $cur_peak_chr = $cur_entry_fields[1];
		    my $cur_peak_coord = $cur_entry_fields[2]+1;
		    my $cur_peak_start_coord = $cur_peak_coord-2;
		    my $cur_peak_end_coord = $cur_peak_coord+2;;
		    my $cur_peak_score = $cur_entry_fields[4];
		    
#		    print output_file "$cur_peak_chr $cur_peak_start_coord $cur_peak_end_coord $cur_peak_id $cur_peak_score + $cur_peak_start_coord $cur_peak_end_coord $peaks_color\n";
		}
	    }
        }
    }
}
    
#print "read $peak_counter individual peaks\n";

close regions_file;

my $peaks_track_name = $track_name . "_peaks";
my $peaks_track_priority=$track_priority+1;

print output_file "track name=$peaks_track_name description=$peaks_track_name itemRgb=\"On\" priority=$peaks_track_priority visibility=1\n";


open regions_file, "< $regions_fname" || die "$regions_fname: $!\n";


while(<regions_file>){
    chomp;
    my $cur_entry = $_;

    if( length($cur_entry) > 0 ){
	if(substr($cur_entry, 0, 1) ne "#"){
	    
	    my @cur_entry_fields = split(/ /, $cur_entry);
	    if(scalar(@cur_entry_fields) >= 9){
		
		if(substr($cur_entry,0,1) eq "R"){
		    my $cur_region_chr = $cur_entry_fields[1];
		    my $cur_region_range = $cur_entry_fields[2];
		    my @cur_region_range_fields = split(/-/, $cur_region_range);
		    if(scalar(@cur_region_range_fields) == 2){
			my $cur_region_start_coord = $cur_region_range_fields[0] + 1;
			my $cur_region_end_coord = $cur_region_range_fields[1] + 1;
			my $cur_region_id = $cur_entry_fields[0];
			my $cur_region_score = $cur_entry_fields[4];
			#print "$cur_region_chr $cur_region_start_coord $cur_region_end_coord\n";
			#print output_file "$cur_region_chr $cur_region_start_coord $cur_region_end_coord $cur_region_id $cur_region_score + $cur_region_start_coord $cur_region_end_coord $regions_color\n";
		    }
		}		
		elsif(substr($cur_entry,0,1) eq "P"){
		    my $cur_peak_id = $cur_entry_fields[0];
		    
		    my $cur_peak_chr = $cur_entry_fields[1];
		    my $cur_peak_coord = $cur_entry_fields[2]+1;
		    my $cur_peak_start_coord = $cur_peak_coord-2;
		    my $cur_peak_end_coord = $cur_peak_coord+2;;
		    my $cur_peak_score = $cur_entry_fields[4];
		    
		    print output_file "$cur_peak_chr $cur_peak_start_coord $cur_peak_end_coord $cur_peak_id $cur_peak_score + $cur_peak_start_coord $cur_peak_end_coord $peaks_color\n";
		}
	    }
        }
    }
}
    
#print "read $peak_counter individual peaks\n";

close regions_file;


close output_file;

#close output_file;
#close reject_file;
