#!/usr/bin/perl

## this program takes align file from the solexa output
## and outputs the bed track readable into the 
## ucsc genome browser
use strict;
use warnings;
use diagnostics;

#use lib '/home/shogo/tools/ucsc_tracks';

my @color_alias;
my @color_value;

sub init_colors{
    push(@color_alias,"DarkSlateGray"); push(@color_value, "047,079,079");
    push(@color_alias,"DarkSlateGrey"); push(@color_value, "047,079,079");
    push(@color_alias,"DimGray"); push(@color_value, "105,105,105");
    push(@color_alias,"DimGrey"); push(@color_value, "105,105,105");
    push(@color_alias,"SlateGray"); push(@color_value, "112,128,144");
    push(@color_alias,"SlateGrey"); push(@color_value, "112,128,144");
    push(@color_alias,"LightSlateGray"); push(@color_value, "119,136,153");
    push(@color_alias,"LightSlateGrey"); push(@color_value, "119,136,153");
    push(@color_alias,"MidnightBlue"); push(@color_value, "025,025,112");
    push(@color_alias,"navy"); push(@color_value, "000,000,128");
    push(@color_alias,"NavyBlue"); push(@color_value, "000,000,128");
    push(@color_alias,"CornflowerBlue"); push(@color_value, "100,149,237");
    push(@color_alias,"DarkSlateBlue"); push(@color_value, "072,061,139");
    push(@color_alias,"SlateBlue"); push(@color_value, "106,090,205");
    push(@color_alias,"MediumSlateBlue"); push(@color_value, "123,104,238");
    push(@color_alias,"LightSlateBlue"); push(@color_value, "132,112,255");
    push(@color_alias,"MediumBlue"); push(@color_value, "000,000,205");
    push(@color_alias,"RoyalBlue"); push(@color_value, "065,105,225");
    push(@color_alias,"blue"); push(@color_value, "000,000,255");
    push(@color_alias,"DodgerBlue"); push(@color_value, "030,144,255");
    push(@color_alias,"DeepSkyBlue"); push(@color_value, "000,191,255");
    push(@color_alias,"SkyBlue"); push(@color_value, "135,206,235");
    push(@color_alias,"LightSkyBlue"); push(@color_value, "135,206,250");
    push(@color_alias,"SteelBlue"); push(@color_value, "070,130,180");
    push(@color_alias,"LightSteelBlue"); push(@color_value, "176,196,222");
    push(@color_alias,"LightBlue"); push(@color_value, "173,216,230");
    push(@color_alias,"PowderBlue"); push(@color_value, "176,224,230");
    push(@color_alias,"PaleTurquoise"); push(@color_value, "175,238,238");
    push(@color_alias,"DarkTurquoise"); push(@color_value, "000,206,209");
    push(@color_alias,"MediumTurquoise"); push(@color_value, "072,209,204");
    push(@color_alias,"turquoise"); push(@color_value, "064,224,208");
    push(@color_alias,"cyan"); push(@color_value, "000,255,255");
    push(@color_alias,"CadetBlue"); push(@color_value, "095,158,160");
    push(@color_alias,"MediumAquamarine"); push(@color_value, "102,205,170");
    push(@color_alias,"aquamarine"); push(@color_value, "127,255,212");
    push(@color_alias,"DarkGreen"); push(@color_value, "000,100,000");
    push(@color_alias,"DarkOliveGreen"); push(@color_value, "085,107,047");
    push(@color_alias,"DarkSeaGreen"); push(@color_value, "143,188,143");
    push(@color_alias,"SeaGreen"); push(@color_value, "046,139,087");
    push(@color_alias,"MediumSeaGreen"); push(@color_value, "060,179,113");
    push(@color_alias,"LightSeaGreen"); push(@color_value, "032,178,170");
    push(@color_alias,"PaleGreen"); push(@color_value, "152,251,152");
    push(@color_alias,"SpringGreen"); push(@color_value, "000,255,127");
    push(@color_alias,"LawnGreen"); push(@color_value, "124,252,000");
    push(@color_alias,"green"); push(@color_value, "000,255,000");
    push(@color_alias,"chartreuse"); push(@color_value, "127,255,000");
    push(@color_alias,"MediumSpringGreen"); push(@color_value, "000,250,154");
    push(@color_alias,"LimeGreen"); push(@color_value, "050,205,050");
    push(@color_alias,"YellowGreen"); push(@color_value, "154,205,050");
    push(@color_alias,"ForestGreen"); push(@color_value, "034,139,034");
    push(@color_alias,"OliveDrab"); push(@color_value, "107,142,035");
    push(@color_alias,"DarkKhaki"); push(@color_value, "189,183,107");
    push(@color_alias,"PaleGoldenrod"); push(@color_value, "238,232,170");
    push(@color_alias,"gold"); push(@color_value, "255,215,000");
    push(@color_alias,"LightGoldenrod"); push(@color_value, "238,221,130");
    push(@color_alias,"goldenrod"); push(@color_value, "218,165,032");
    push(@color_alias,"DarkGoldenrod"); push(@color_value, "184,134,011");
    push(@color_alias,"RosyBrown"); push(@color_value, "188,143,143");
    push(@color_alias,"IndianRed"); push(@color_value, "205,092,092");
    push(@color_alias,"SaddleBrown"); push(@color_value, "139,069,019");
    push(@color_alias,"sienna"); push(@color_value, "160,082,045");
    push(@color_alias,"peru"); push(@color_value, "205,133,063");
    push(@color_alias,"burlywood"); push(@color_value, "222,184,135");
    push(@color_alias,"wheat"); push(@color_value, "245,222,179");
    push(@color_alias,"SandyBrown"); push(@color_value, "244,164,096");
    push(@color_alias,"tan"); push(@color_value, "210,180,140");
    push(@color_alias,"chocolate"); push(@color_value, "210,105,030");
    push(@color_alias,"firebrick"); push(@color_value, "178,034,034");
    push(@color_alias,"brown"); push(@color_value, "165,042,042");
    push(@color_alias,"DarkSalmon"); push(@color_value, "233,150,122");
    push(@color_alias,"salmon"); push(@color_value, "250,128,114");
    push(@color_alias,"LightSalmon"); push(@color_value, "255,160,122");
    push(@color_alias,"orange"); push(@color_value, "255,165,000");
    push(@color_alias,"DarkOrange"); push(@color_value, "255,140,000");
    push(@color_alias,"coral"); push(@color_value, "255,127,080");
    push(@color_alias,"LightCoral"); push(@color_value, "240,128,128");
    push(@color_alias,"tomato"); push(@color_value, "255,099,071");
    push(@color_alias,"OrangeRed"); push(@color_value, "255,069,000");
    push(@color_alias,"red"); push(@color_value, "255,000,000");
    push(@color_alias,"HotPink"); push(@color_value, "255,105,180");
    push(@color_alias,"DeepPink"); push(@color_value, "255,020,147");
    push(@color_alias,"pink"); push(@color_value, "255,192,203");
    push(@color_alias,"LightPink"); push(@color_value, "255,182,193");
    push(@color_alias,"PaleVioletRed"); push(@color_value, "219,112,147");
    push(@color_alias,"maroon"); push(@color_value, "176,048,096");
    push(@color_alias,"MediumVioletRed"); push(@color_value, "199,021,133");
    push(@color_alias,"VioletRed"); push(@color_value, "208,032,144");
    push(@color_alias,"magenta"); push(@color_value, "255,000,255");
    push(@color_alias,"violet"); push(@color_value, "238,130,238");
    push(@color_alias,"plum"); push(@color_value, "221,160,221");
    push(@color_alias,"orchid"); push(@color_value, "218,112,214");
    push(@color_alias,"MediumOrchid"); push(@color_value, "186,085,211");
    push(@color_alias,"DarkOrchid"); push(@color_value, "153,050,204");
    push(@color_alias,"DarkViolet"); push(@color_value, "148,000,211");
    push(@color_alias,"BlueViolet"); push(@color_value, "138,043,226");
    push(@color_alias,"MediumPurple"); push(@color_value, "147,112,219");
    push(@color_alias,"thistle"); push(@color_value, "216,191,216");
};

sub get_color{
    my $cur_color_alias = shift;
    for(my $i=0; $i<scalar(@color_alias); $i++){
        if($cur_color_alias eq $color_alias[$i]){
            return $color_value[$i];
        }
    }
    return "not found";
};

sub get_random_color{
    use POSIX qw(ceil floor);
    my $cur_rand_int = floor(rand() * (scalar(@color_value) ) );
    return $color_value[$cur_rand_int];
};


my $usage = qq{
    QuEST_align_2_bed.pl

    -----------------------------
    mandatory parameters:
    
    -i  input_QuEST_align_file     QuEST file containing read alignments
    -o  output_prefix              output bed file  
    -gt genome_table                      

    or

    -p param_file                 parameter file
    
    -----------------------------
    optional parameters:
   
    -s  feature_shift             how much to the right shift all features
    -fs feature_size              size of the feature in the browser
    -gz gzip_enabled              specify this flag to get gz file

    -fc forward_color             specify "x,y,z" corresponding to + features
    -rc reverse_color             specify "x,y,z" corresponding to - features

    -t  track_name                 the string that will appear on the ucsc track
    
    -c  upper_count                limit the number of reads to this number

};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments

my $input_fname = "";
my $output_fname = "";
my $param_fname = "";
my $genome_table_fname = "";
my $output_by_chr_fname_prefix = "";

my $no_tag_limit = "true";

## optional arguments

my $feature_shift = 0; #-1;
my $feature_size = 5;
my $gz_flag = "yes";

init_colors();

my $forward_color;
my $reverse_color;

my $min_dist = 200;
my $cur_dist = 0;

while($cur_dist < $min_dist){
    $forward_color = get_random_color();
    $reverse_color = get_random_color();
    
    my @fc_fields = split(/,/,$forward_color);
    my @rc_fields = split(/,/,$reverse_color);
   
    my $fc1 = $fc_fields[0];
    my $fc2 = $fc_fields[1];
    my $fc3 = $fc_fields[2];

    my $rc1 = $rc_fields[0];
    my $rc2 = $rc_fields[1];
    my $rc3 = $rc_fields[2];

    $cur_dist = sqrt( ($fc1-$rc1)*($fc1-$rc1) + ($fc2-$rc2)*($fc2-$rc2) + ($fc3-$rc3)*($fc3-$rc3) );
}

#print "fc: $forward_color rc: $reverse_color\n";
while($reverse_color eq $forward_color){ $reverse_color = main::get_random_color(); }

my $track_name = "QuEST_align_track";
my $track_priority = "0";
my $tag_limit;



## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    #printf "$this_arg\n";
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-i') {$input_fname = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_by_chr_fname_prefix = shift @ARGV;}
    elsif ( $this_arg eq '-gt') {$genome_table_fname = shift @ARGV;}
    elsif ( $this_arg eq '-p') {$param_fname = shift @ARGV;}

    elsif ( $this_arg eq '-s') {$feature_shift = shift @ARGV;}
    elsif ( $this_arg eq '-fs') {$feature_size = shift @ARGV;}
    elsif ( $this_arg eq '-gz') {$gz_flag = "yes";}

    elsif ( $this_arg eq '-fc') {$forward_color = shift @ARGV;}
    elsif ( $this_arg eq '-rc') {$reverse_color = shift @ARGV;}

    elsif ( $this_arg eq '-t') {$track_name = shift @ARGV;}

    elsif ( $this_arg eq '-c') {$tag_limit = shift @ARGV; $no_tag_limit = "false";}

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n"; exit(1); }
}


## print parameters

if($param_fname ne ""){
    open param_file, "< $param_fname" || die "$param_fname: $!\n";
    
    while(<param_file>){
	chomp;
	my $cur_line = $_;
	if(length($cur_line) > 0){
	    
	    if(substr($cur_line, 0, 1) ne "#"){
		my @cur_entry_fields = split(/ /, $cur_line);
		if(scalar(@cur_entry_fields) >= 3 ){
		    
		    my $cur_param_name = $cur_entry_fields[0];
		    my $cur_eq_symbol = $cur_entry_fields[1];
		    my $cur_param_value = $cur_entry_fields[2];

		    if($cur_eq_symbol eq "="){
			if($cur_param_name eq "genome_table"){
			    $genome_table_fname = $cur_param_value;			    
			}
			elsif($cur_param_name eq "input_file"){
			    $input_fname = $cur_param_value;
			}
			elsif($cur_param_name eq "output_file"){
			    $output_fname = $cur_param_value;
			}
			elsif($cur_param_name eq "track_name"){
			    $track_name = $cur_param_value;
			}
			elsif($cur_param_name eq "output_by_chr_prefix"){
			    $output_by_chr_fname_prefix = $cur_param_value;			    
			}
			elsif($cur_param_name eq "track_priority"){
			    $track_priority = $cur_param_value;			    
			}
			elsif($cur_param_name eq "gz"){
			    $gz_flag = $cur_param_value;
			}
			else{
			    print "Error in QuEST_align_2_bed.pl:\n";
			    print "Failed to understand parameter $cur_line in $param_fname\n";
			    exit(4);
			}
		    }
		}
	    }
	}
    }
    
    close param_file;
}
else{
    $output_fname = $output_by_chr_fname_prefix . ".bed";
}

if ( $input_fname eq '' || $output_fname eq ''){
    die "Either input or output file wasn't specified in QuEST_align_2_bed.pl\n";
}

print "\n-----------------\n\n";
print "input file: $input_fname\n";
print "output file: $output_fname\n";
print "param file: $param_fname\n";
print "genome_table: $genome_table_fname\n";
print "output_by_chr_prefix: $output_by_chr_fname_prefix\n";

print "feature shift: $feature_shift\n";
print "feature size: $feature_size\n";
print "gz: $gz_flag\n";

print "forward color: $forward_color\n";
print "reverse color: $reverse_color\n";
printf "color_dist: %.2f\n", $cur_dist;

print "track name: $track_name\n";
print "track priority: $track_priority\n";

if($no_tag_limit eq "false"){ print "tag_limit: $tag_limit\n";}
else{ print "tag_limit: NA\n"; }

my @chr_names;
my @chr_sizes;
my $chr_counter = 0;

open genome_table_file, "< $genome_table_fname" || die "$genome_table_fname: $!\n";
while(<genome_table_file>){
    chomp;
    my $cur_line = $_;
    if(length($cur_line) > 0){
	if(substr($cur_line,0,1) ne "#"){
	    my @cur_entry_fields = split(/ /,$cur_line);
	    if(scalar(@cur_entry_fields) >= 2){
		my $cur_chr = $cur_entry_fields[0];
		my $cur_chr_size = $cur_entry_fields[1];
		if($cur_chr_size =~ /^[+-]?\d+$/){
		    # ok, it's a number
		}
		else{
		    printf "Error in the QuEST_align_2_bed.pl: Bad chr size\n";
		    printf "Expected an integer, but found \"$cur_chr_size\"\n";
		    printf "In the line $cur_line\n";
		    exit(1);
		}
		$chr_names[$chr_counter] = $cur_chr;
		$chr_sizes[$chr_counter] = $cur_chr_size;
		$chr_counter++;		
		
		#print "$cur_chr\n";
	    }
	}
    }
}

close genome_table_file;



print "\n-----------------\n\n";

open input_file, "< $input_fname" || die "$input_fname: $!\n";
my $output_fname_gz = $output_fname . ".gz";

if(-e $output_fname){
    unlink("$output_fname");
}
if(-e $output_fname_gz){
    unlink($output_fname_gz);
}
open output_file, "> $output_fname" || die "$output_fname: $!\n";

#print output_file "browser position chr1: 5000000-5200000\n";
#print output_file "browser hide all\n";
print output_file "track name=\"$track_name\" description=\"$track_name\" priority=$track_priority visibility=4 itemRgb=\"On\"\n";

my @reads_by_chr_for;
my @reads_by_chr_for_sizes;

my @reads_by_chr_rev;
my @reads_by_chr_rev_sizes;

for(my $i=0; $i<scalar(@chr_names); $i++){
    my @dummy_array_for;
    my @dummy_array_rev;

    $reads_by_chr_for_sizes[$i] = 0;
    $reads_by_chr_for[$i] = \@dummy_array_for;

    $reads_by_chr_rev_sizes[$i] = 0;
    $reads_by_chr_rev[$i] = \@dummy_array_rev;
}

my $mapped_features = 0;
my $unmapped_features = 0;
my $reads_found = 0;

#print "analyzing...\n";
while(<input_file>){
    chomp;
    # print "$_\n";
    my @cur_fields = split(/ /, $_);

#    if($reads_found % 1000 == 0){
#	printf "\r%.3f M reads found, %d unmapped", $reads_found/1000000, $unmapped_features;
#	$|++;
#    }

    if(scalar(@cur_fields) == 3){
	#my $cur_read = "solexa_read";
	my $cur_chr = $cur_fields[0];
	my $cur_coord = $cur_fields[1];
	my $cur_orient = $cur_fields[2];    		

	if($cur_coord =~ /^[+-]?\d+$/){
	    ## ok, it's a number
	}
	else{
	    printf "Error in QuEST_align_2_bed.pl: \n";
	    printf "Expected a number, but found $cur_coord for the alignment \"$_\"\n";
	    exit(1);
	}
	
	my $cur_chr_ind = -1;
	for(my $i=0; $i<scalar(@chr_names); $i++){
	    if($cur_chr eq $chr_names[$i]){
		$cur_chr_ind = $i;
	    }
	}
	
	if($cur_chr_ind >= 0){
	    if($cur_orient eq "+"){
		$reads_by_chr_for[$cur_chr_ind][$reads_by_chr_for_sizes[$cur_chr_ind]] = $cur_coord;
		$reads_by_chr_for_sizes[$cur_chr_ind]++;
		$mapped_features++;
	    }
	    elsif($cur_orient eq "-"){
		$reads_by_chr_rev[$cur_chr_ind][$reads_by_chr_rev_sizes[$cur_chr_ind]] = $cur_coord;
		$reads_by_chr_rev_sizes[$cur_chr_ind]++;
		$mapped_features++;
	    }
	    else{
		printf "Error in QuEST_align_2_bed.pl\n";
		printf "Expected +/-, but found $cur_orient for the aligment line \"$_\"\n";
		exit(1);
		$unmapped_features++;
	    }
	}

	$reads_found++;
    }
}
close input_file;

print "\n\n";

my $skipped_reads = 0;

for(my $i=0; $i<scalar(@chr_names); $i++){
    my $cur_chr = $chr_names[$i];
    my $cur_chr_size = $chr_sizes[$i];
    my $cur_chr_output_fname = $output_by_chr_fname_prefix . ".$cur_chr.bed";
    my $cur_chr_output_fname_gz = $output_by_chr_fname_prefix . ".$cur_chr.bed.gz";
    if(-e $cur_chr_output_fname){
	unlink($cur_chr_output_fname);
    }    
    if(-e $cur_chr_output_fname_gz){
	unlink($cur_chr_output_fname_gz);	
    }

    my $cur_track_name = $track_name . "_$cur_chr";
    
    open cur_chr_output_file, "> $cur_chr_output_fname" || die "$cur_chr_output_fname: $!\n";
    print cur_chr_output_file "track name=\"$cur_track_name\" description=\"$cur_track_name\" priority=$track_priority visibility=4 itemRgb=\"On\"\n";

    for(my $j=0; $j<$reads_by_chr_for_sizes[$i]; $j++){
	my $cur_coord = $reads_by_chr_for[$i][$j];
	my $cur_read = "QuEST_align";
	my $start_coord = $cur_coord + 1 + $feature_shift;
	my $end_coord = $cur_coord + 1 + $feature_shift + $feature_size - 1;

	if($end_coord <= $start_coord){ die "wrong order start->end: $end_coord < $start_coord\n";}
	if($start_coord <= $end_coord && $start_coord >= 0 && $end_coord < $cur_chr_size){
	    printf cur_chr_output_file "%s %d %d %s 0 + %d %d %s\n", $cur_chr, 
	    $start_coord, $end_coord, $cur_read, $start_coord, $end_coord, $forward_color;
	    printf output_file "%s %d %d %s 0 + %d %d %s\n", $cur_chr, 
	    $start_coord, $end_coord, $cur_read, $start_coord, $end_coord, $forward_color;

	}
	else{ $skipped_reads++; }
    }
    
    for(my $j=0; $j<$reads_by_chr_rev_sizes[$i]; $j++){
	my $cur_coord = $reads_by_chr_rev[$i][$j];
	my $cur_read = "QuEST_align";
	
	my $start_coord = $cur_coord + 1 - $feature_size + $feature_shift;
	my $end_coord = $cur_coord + 1 + $feature_shift;

	if($end_coord <= $start_coord){ die "wrong order start->end: $end_coord < $start_coord\n";}
	if($start_coord <= $end_coord && $start_coord >= 0 && $end_coord < $cur_chr_size){
	    printf cur_chr_output_file "%s\t%d\t%d\t%s\t0\t-\t%d\t%d\t%s\n", $cur_chr,
	    $start_coord, $end_coord, $cur_read, $start_coord, $end_coord, $reverse_color;	    
	    printf output_file "%s\t%d\t%d\t%s\t0\t-\t%d\t%d\t%s\n", $cur_chr,
	    $start_coord, $end_coord, $cur_read, $start_coord, $end_coord, $reverse_color;	    
	}
	else{ $skipped_reads++; }
    }

    close cur_chr_output_file;
    system("gzip $cur_chr_output_fname &");
}


close output_file;
close input_file;


if ($gz_flag eq "yes"){ 
    unlink("$output_fname.gz");
    system("gzip $output_fname &");
}

