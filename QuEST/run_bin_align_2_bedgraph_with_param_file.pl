#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    bottleneck_metrick.pl

    -----------------------------
    mandatory parameters:
    
    -b  bin_align_prefix
    -o  output_file            
    -n  track_name
    -gt genome_table
    -ct count_threshold
    -cp by_chr_output_prefix

    or 

    -p  param_file

    -----------------------------
    optional parameters:
    -h (help)
    -tp track_priority;    

};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

use Cwd qw(realpath);
my $fullpath = realpath($0);
my @fullpath_fields = split(/\//,$fullpath);
my $exec_path = "";

for(my $i=0; $i<scalar(@fullpath_fields)-1; $i++){
    if($fullpath_fields[$i] ne ""){
        $exec_path = $exec_path."/".$fullpath_fields[$i];
    }
}

$exec_path = $exec_path."/bin_align_2_bedgraph";

if( ! -e $exec_path ){
    print "Error in generate profile batch script:\n";
    print "Failed to locate executable $exec_path.\n";
    print "Aborting.\n";
    exit(0);
}

## mandatory arguments

my $undef = "undef";

my $bin_align_prefix = "";
my $output_fname = "";
my $genome_table_fname = "";
my $track_name = "";
my $track_priority = 10;
my $count_threshold = 1;
my $param_fname = $undef;
my $by_chr_prefix = "";
## optional arguments

## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-b') {$bin_align_prefix = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_fname = shift @ARGV;}
    elsif ( $this_arg eq '-n') {$track_name = shift @ARGV;}
    elsif ( $this_arg eq '-tp') {$track_priority = shift @ARGV;}
    elsif ( $this_arg eq '-gt') {$genome_table_fname = shift @ARGV;}
    elsif ( $this_arg eq '-ct') {$count_threshold = shift @ARGV;}
    elsif ( $this_arg eq '-cp') {$by_chr_prefix = shift @ARGV;}
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
	    
	    if($cur_par_name eq "bin_align_prefix"){
                $bin_align_prefix = $cur_par_value;
            }
	    if($cur_par_name eq "by_chr_prefix"){
                $by_chr_prefix = $cur_par_value;
            }
	    if($cur_par_name eq "genome_table"){
                $genome_table_fname = $cur_par_value;
            }
	    if($cur_par_name eq "count_threshold"){
                $count_threshold = $cur_par_value;
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


if( $output_fname eq ""){
    die "you should specify output file\n";
}

## print parameters

print "\n-----------------\n\n";
print "bin_align_prefix: $bin_align_prefix\n";
print "genome_table: $genome_table_fname\n";
print "output file: $output_fname\n";
print "track name:  $track_name\n";
print "count_threshold:  $count_threshold\n";
print "by_chr_prefix: $by_chr_prefix\n";
print "\n-----------------\n\n";

my @by_chr_files = <$by_chr_prefix*>;
#print @by_chr_files;
#for(my $i=0; $i<scalar(@by_chr_files); $i++){
    #unlink($by_chr_files[$i]);
#}
foreach my $file (@by_chr_files){
    #print $file . "\n";
    unlink($file);
}


if(-e $output_fname){
    unlink($output_fname);
}
if(-e $output_fname . ".gz"){
    unlink($output_fname . ".gz");
}

#my $system_command = "$exec_path bin_align_file_prefix=$bin_align_prefix output_file=$output_fname track_name=$track_name genome_table=$genome_table_fname track_priority=$track_priority count_threshold=$count_threshold by_chr_path=$by_chr_path; rm $by_chr_path*gz; rm $output_fname*gz; gzip $by_chr_path*; gzip $output_fname";
my $system_command = "$exec_path bin_align_file_prefix=$bin_align_prefix output_file=$output_fname track_name=$track_name genome_table=$genome_table_fname track_priority=$track_priority count_threshold=$count_threshold by_chr_prefix=$by_chr_prefix; gzip $by_chr_prefix*; gzip $output_fname";

print "> $system_command\n";
my $error_code = system("$system_command\n");
if($error_code != 0){
    print "run_bin_align_2_bedgraph_with_param_file.pl: Caught an error message (code $error_code), passing error message to the top script.\n";
    exit(3);
}

