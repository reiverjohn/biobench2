QuEST 2.4

Table of Contents

I. License and Terms of Use
II. Disclosures
III. QuEST Reference and Home Page [renumber subsequent sections]
IV. Changes
V. Software summary
VI. Quick start guide
VII. File formats
VIII. QuEST overview
IX. Explanation of QuEST flags
X. QuEST parametrization
XI. FAQ

===============================================================
I. Licence and Terms of Use
===============================================================

This code is written by Anton Valouev and is copyrighted under
the Lesser GPL:

Copyright (C) 2009 Anton Valouev

 This program is free software; you can redistribute it and/or modify 
 it under the terms of the GNU Lesser General Public License as       
 published by the Free Software Foundation; version 2.1 or later.    
 This program is distributed in the hope that it will be useful,     
 but WITHOUT ANY WARRANTY; without even the implied warranty of      
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU Lesser General Public License for more details. 
 You should have received a copy of the GNU Lesser General Public 
 License along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 MA 02111-1307, USA.
 
The author may be contacted via email at: valouev@stanford.edu and
 valouev@gmail.com

==================================================================
II. Disclosures
==================================================================

current version: 2.4
current development status: stable release

#################################################################
#                                                               #
#   If you discover any bugs or have any suggestions,           #
#   please contact Anton Valouev (valouev@stanford.edu).        #
#                                                               #
#################################################################

==================================================================
III. QuEST Reference and Home Page

==================================================================
QuEST is a statistical package for analysis of ChIP-Seq data.

If you use QuEST or any derivatives of its results, or any part of
QuEST for any purpose, please cite the following in your 
publications, patents and code:

 Valouev A, Johnson DS, Sundquist A, Medina C, Anton E, Batzoglou S,
 Myers RM, Sidow A (2008). Genome-wide analysis of transcription
 factor binding sites based on ChIP-Seq data.
 Nat Methods 5, 829-834 (2008)

QuEST Home Page:
http://mendel.stanford.edu/downloads/quest/

==============================================================
IV. Changes
==============================================================
QuEST 2.4
   1. Incorporated region filter. If background data is provided, QuEST
      will scan the background data to detect artifactual regions where
      where enrichment fold of the control data exceeds 100 fold and
      will extend until 5 fold enrichment boundary. Any ChIP regions
      overlapping artifactual control regions will be filtered out
   2. Added eland_extended format for the input data
   3. Added bowtie format for the input data (bowtie alignments generated
      with -z flag are not supported. -z flag will cause bowtie to report
      alignments in non-contiguous order. Currently stick to bowtie
      without -z flag).
   4. Added maq alignment code. Output of maq mapview command can now
      be directly feeded into QuEST.
   5. Added report files for importing alignments (check align_2_QuEST
      report file under module outputs for the import error log)
   6. Added report file for the metrics module (check metrics.ChIP.report.txt
      under the module outputs directory). You can now find out how
      much ChIP/control data false into the called regions. This is a good
      metric for the ChIP quality
   7. Added report file for the peak filter module (check the peak_filter.ChIP.report.txt
      under the module outputs directory)
   8. Added report file for the QuEST parametrization (check the QuEST_parametrization.report.txt
      for the complete set of QuEST parametrization parameters)
      
QuEST 2.3
   1. QuEST ChIP scores (after normalization) now represent fold enrichment
      and are thus more interpretable
   2. QuEST configuration script now operates with fold enrichments rather
      than KDE score cutoffs. Internally QuEST still uses KDE scores, but
      they are now not visible to the user directly and parametrization
      occurs using normalized ChIP scores (fold enrichment)
   3. Eliminated rescue fold, stringency relative to background (input) can
      be set using ChIP_to_background fold enrichment parameter during parametrization
   4. Added a filter to the "TF with motif and a narrow (punctate) binding site"
      option that only reports a single peak per region and removes regions with
      no peaks in the "peak_caller.ChIP.out.accepted" file. All initial regions
      and peaks before filtering can be found in the "peak_caller.ChIP.out.with_metrics"
      file
   5. Added tag collapsing that removes stacks from the input data. 
   6. Removed oversampling metric throughout as it wasn't found to be very useful
   7. Substituted Bonferroni correction with Benjamini-Hochberg correction for p-values
   8. Added -log10_pv field to the output file to provide p-values

QuEST 2.2
   1. added a new data visualization output through the bed graph UCSC genome browser display
   2. improved region boundary detection by applying "ChIP extension threshold". Now
      regions are seeded using "ChIP threshold" and boundaries are defined
      by the "ChIP extension threshold".
      

QuEST 2.0
   Parts of QuEST were rewritten to significantly speed up the execution.
      
   QuEST now incorporates a number of improvements and additional
   features listed below. QuEST now produces bed and wig files 
   that can be viewed on the UCSC genome browser. 

   1. Global peak shift estimation.
      Global peak shift estimation is now completely re-written. In
      the initial pass candidate regions are identified based on the 
      tag count threshold and tag count enrichment threshold compared to control
      experiment. For each candidated region local peak shift estimate
      is given by the distance for which correlation between opposit
      strand profiles is maximized. The global peak shift estimate
      is given by the mode of the peak shift density estimate
      obtained using KDE on local peak shift estimates from candidate
      regions. This approach provides substantial improvement over 
      previously implemented method.

   2. Local peak shift estimate
      Local peak shift estimate is a new metric that can be used 
      to filter out regions where peak shift falls outside the 
      expected range (typically 30-100 bp). Local peak shift estimate
      is given by the shift between opposite strand density profiles
      where correlation between profile values within a window is maximized.

   3. Strand correlation values
      Strand correlation values reported by QuEST is a new metric
      that can be used to filter out regions and peaks. Strand correlation
      measures similarity between density profiles of +/- strands
      after correcting for the local peak shift. Regions with low values
      of correlation (0.5 or less) are likely to be artifacts and not
      real regions of enrichment.

   4. Tag oversampling
      Tag oversampling is a new metric that can be used to filter out regions
      resulting from mapping artifacts and repetetive regions. Oversampling
      measures the ratio of observed mean of tag instances at positions
      with mapped tags within the region compared to expectation of tag instances
      under uniform distribution. High values of oversampling usually indicate
      outliers such as loci with only a few positions where reads start but
      with many reads starting at these few positions. Outliers usually appear
      to have oversampling values of 10 or more. Good ChIP-enrichment regions
      usually have oversampling values of 1-5. Be aware that bottlenecked libraries
      (i.e. where the same template is typically sequenced multiple times) will
      result in high values of oversampling even in the correct regions.

   5. Q-values for regions
      Q-values for regions represent Poisson p-values based on the control
      experiment corrected for the number of tests (genome size). 

   6. Q-values for peaks
      Q-values for peaks represent approximated FDR for the peak intensity.
      Control data profile is first analyzed to build the peak num vs profile
      threshold curve which is than approximated by the exponential curve
      in the form c*exp(-score/mu) for large values of score. This curve is
      then applied with the scale to the ChIP data to approximate FDN for
      any score value and the q-value is then calculated by deviding FDN
      by the total number of peaks identified.


===============================================================
V. Software Summary
===============================================================

Target OS: Unix-type systems

QuEST was tested on Linux ES with kernel 2.6.9, and Mac OS 10.5.
Software requirements: perl, make, g++ compiler
Hardware requirements: 1Gb of RAM, ~ 20 Gb of disc space for 
temporary files. 

===============================================================

What you will need to run QuEST: 

1. ChIP data

   Reads from the ChIP-Seq experiment. 
   
2. RX_noIP data
  
   Reads from the RX_noIP experiment. This is also known
   as total chromatin or input.

The following alignment formats are currently supported in QuEST (see below for the format descriptions): 
	     - Solexa align
	     - Eland align
	     - Eland extended align
	     - Bowtie align (except for generated with -z option)
	     - MAQ align (generated by running MAQ mapview command)	
	     - QuEST align

3. Reference genome or genome table (see below for the format 
   description)

   For the reference genome, there should be one contig per .fa file.
   Alternatively, you can provide a genome table containing chromosome
   names and their sizes (see below for the detailed format 
   description). Supported formats: fasta, genome table. Note that 
   you need either the analysis path or the genome table, but not 
   both, or else the script will not run.

==============================================================
VI. Quick Start Guide
==============================================================

To run QuEST, the following commands should execute on your system:
make, g++, perl. 

1. Configure QuEST to run with your data. 
   This will generate all the necessary configuration files and 
   thresholds:

   command: generate_QuEST_parameters.pl

   flags:   
	    -solexa_align_ChIP <file containing ChIP reads in
			       solexa align format (see below)>

	    -solexa_align_RX_noIP <file containing RX_noIP 
				  (aka input/total chromatin) 
				  reads in solexa align format
				  (see below)>

	    -rp <path to the reference genome containing .fa files>
	    -ap <analysis path where the QuEST results and 
		all files will be located> 

            -silent (this will assume all the default parameters 
		    and will not prompt the user)

   example:

   $/home/agent_smith/QuEST/generate_QuEST_parameters.pl -solexa_align_ChIP /Volumes/Projects/GABP/GABP.align_25.hg18.txt -solexa_align_RX_noIP /Volumes/Projects/background/Jurkat_RX_noIP.align_25.hg18.txt -rp /Volumes/Projects/Genomes/Human/hg18/ -ap /Volumes/Projects/GABP/analysis/

   all flags and help:

   to list all flags execute "generate_QuEST_parameters.pl"

2. Run QuEST:

   command: run_QuEST_with_param_file.pl

   flags: 
          -ap <path to the analysis directory>

   example:

   $/home/agent_smith/QuEST/run_QuEST_with_param_file.pl -ap /Volumes/Projects/GABP/analysis/

3. QuEST will output the following files of interest:

QuEST summary: QuEST.out

QuEST calls:

  QuEST calls for ChIP peaks and regions: analysis_directory/calls/peak_caller.ChIP.out
  QuEST calls with additional metrics: analysis_direcotry/calls/peak_caller.ChIP.out.with_metrics
  QuEST calls that passed the quality filter: analysis_directory/calls/peak_caller.ChIP.out.accepted
  QuEST calls that did not pass the quality filter: analysis_directory/calls/peac_caller.ChIP.out.rejected

  pseudo_ChIP peak calls: peak_caller.pseudo_ChIP.out	

UCSC genome browser tracks (these are good for visualisation of data and analysis):

  QuEST calls bed track: analysis_director/tracks/ChIP_calls.filtered.bed

  ChIP data bed tracks: analysis_directory/tracks/data_bed_files/ChIP.bed.gz
  (these files can sometimes be too large to load into UCSC genome browser,
  but you can use by chromosome files instead)

  ChIP data bed tracks by chormosome: analysis_directory/tracks/data_bed_files/by_chr/ChIP

  RX_noIP data bed tracks: analysis_directory/tracks/data_bed_files/RX_noIP.bed.gz	
  (these files can sometimes be too large to load into UCSC genome browser,	
  but you can use by chromosome files instead)

  RX_noIP data bed tracks by chormosome: analysis_directory/tracks/data_bed_files/by_chr/RX_noIP

  normalized ChIP enrichment profiles: analysis_directory/tracks/wig_profiles/ChIP_normalized.profile.wig.gz
  normalized ChIP enrichment profiles by chromosome: analysis_directory/tracks/wig_profiles/by_chr/ChIP_normalized/

  unnormalized ChIP enriment profiles: analysis_directory/tracks/wig_profiles/ChIP_unnormalized.profile.wig.gz
  unnormalized ChIP enrichment profiles by chromosome: analysis_directory/tracks/wig_profiles/by_chr/ChIP_unnormalized/	

  normalized background enrichment profiles: analysis_directory/tracks/wig_profiles/background_normalized.profile.wig.gz
  normalized background enrichment profiles by chromosome: analysis_directory/tracks/wig_profiles/by_chr/background_normalized/

  unnormalized background enriment profiles: analysis_directory/tracks/wig_profiles/background_unnormalized.profile.wig.gz
  unnormalized background enrichment profiles by chromosome: analysis_directory/tracks/wig_profiles/by_chr/background_unnormalized/	

Report files:
   module_outputs/align_2_QuEST.ChIP.report.txt
   module_outpus/align_2_QuEST.RX_noIP.report.txt
   these two files contain information about problems that QuEST encountered when importing
   alignments

   module_outputs/metrics.ChIP.report.txt
   contains information about what percentage of ChIP and background fell within the called regions

   module_outputs/QuEST_parametrization.report.txt
   contains information about parametrization choices provided for running QuEST

   module_outputs/peak_filter.ChIP.report.txt
   contains information about result of application of QuEST filters on regions and peaks

==============================================================
VII. File Formats
==============================================================
The following input file formats are supported in QuEST:

1.  Read aligment formats
1.a Solexa/Illumina read alignment files (standard Solexa 2007 align file).
    The fields are space delimited (X means QuEST does not care about
    what is in that particular field as long as it is present). Each
    line corresponds to the alignment of one read and the coordinate
    gives the first base in the reference genome to which the read
    aligns. F/R specifies whether the alignment is to the sense or 
    antisense strand of the reference genome. Coordinates begin with 
    position 1.

Definition:
<DNA sequence> X X <chromosome:coordinate> <F/R>

Example:
GACTACACCTCTAACCACTACCATC 1250 1 chr4:72788542 F

1.b Solexa/Illumina eland/vmf alignment files (extended eland alignment
    files are not yet supported). Each line corresponds to a 
    single alignment, and aligment coordinate corresponds to the first
    aligned base relative to the reference genome.
    Fields are tab delimited (X means QuEST does not care about
    the value of a particular field as long as it is non-empty).
    F/R specifies whether the alignment is to the sense or antisense
    strand of the reference genome. Coordinates begin with position 1.

Definition:
X	<DNA sequence>	X	X	X	X	<chrom_name.fa>	<coordinate>	<F/R>

Example:
>HWI-EAS105_4:4:1:1664:1790     GAAAAGAACTGAAGGATCTCAAAAAT      U1      0       1       0       chr1.fa 37883447        F       DD 1C


1.c QuEST align file
    Sequencing read QuEST align files. Fields are space delimited; 
    each line corresponds to alignment of one read. Coordinate gives 
    position of the 5 prime base coordinate of the read relative to 
    the reference genome. +/- specifies whether the alignment is to 
    the sense or antisense strand of the reference genome. 
    Coordinates begin with position 0.

Definition:
<chromosome> <coordinate> <+/->

1.d Eland extended align file.
    will add description

Example:
>HWUSI-EAS270R:5:1:0:372#0/1    NGTCTCACTGGCAAATCCATTCCCTCACTGACATGA    1:1:1   chr4.fa:135597163RA13G11G3TGG1TG,chr8.fa:113733968FG35

1.e Bowtie align file (except for ones generated using -z option)
    QuEST will select the best alignment

Example: 
18213   +       chrX    19506600        CTGGTGTTAATGTTCAAGCTGTCTCTGGCTGGAGCT    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    21


1.f MAQ align file generated using MAQ mapview command
    will add description

Example:
HWUSI-EAS270R:5:59:47:243#0/1   chr10   3000454 -       0       0       0       0       0       0       0       9       18      36      AAGCTGTTTTCCTGCTACTCTCTGAGCCTGCTACCC    _\^`^`WP]\^`]^`\V`Z`\`\`Z```XX``T```

2. Reference genome files
2.a QuEST genome table file. Specifies sizes of the chromosomes
    to which the aligments were made. The entries are organized one
    chromosome per line, fields are space separated. Note that 
    chromosome names of the aligments and genome table file should 
    match or else QuEST will not use those reads.

Definition:
<chromosome> <size in bp>

Example:
chr18 76117153

2.b Reference genome fasta files. QuEST will accept standard fasta
    files with the reference sequence such as those that can be 
    obtained from UCSC genome browser. Chromosome names in fasta 
    files following '>' symbol must match chromosome names in the 
    alignment files or else QuEST will dismiss those reads. Space 
    is not allowed in the chromosome names. Chromosomes can be 
    organized one per file or all in one file, but make sure no two 
    same chromosomes are present in .fa files in the path provided to
    QuEST or else the size of the genome may be overestimated leading
    to improper choice of QuEST parameters. Using fasta files QuEST
    will produce a genome table file (see above) in the analysis path 
    directory which can be used for this genome in future analysis 
    instead of the fasta files. You can check the genome table file 
    for the correct names of the chromosome and presence of any 
    duplicate chromosome entries and edit it to your liking. 

The following output file formats are adopted in QuEST:

1. QuEST ChIP calls. QuEST calls provide information about regions
   of significant enrichment compared to control experiment. (these
   may be slightly outdated, please check your peak call file)

   calls/peak_caller.ChIP.out
   This file contains initial unfiltered results of the ChIP-Seq caller.
   Entries represent significant regions and peaks and are sorted
   according to the enrichment intensity of the strongest peak
   within the region. Regions are separated by the empty line.
   Each region entry is followed by one or more peak entries 
   corresponding to peak within that region.

   Definition:

   R-<region_id> <chromosome> <region start>-<region end> ChIP: <ChIP enrichment intensity of the strongest peak> control: <control enrichment intensity at the position of the strongest ChIP peak> max_pos: <coordinate of the highest peak in the region> ef: <normalized enrichment fold at the position of the strongest peak>
   P-<region_id>-<peak_id> <chromosome> <peak coordinate> ChIP: <ChIP enrichment intensity> control: <control enrichment intensity> region: <region start>-<region end> ef: <normalized_enrichment fold>

   Example:

   R-1 chr8 123640927-123641979 ChIP: 0.58 control: 0.0027 max_pos: 123641348 ef: 262.007
   P-1-1 chr8 123641348 ChIP: 0.58 control: 0.0027 region: 123640927-123641979 ef: 262.007

   R-2 chr8 87322905-87324677 ChIP: 0.533 control: 3.7e-05 max_pos: 87323637 ef: 17570.083
   P-2-1 chr8 87323192 ChIP: 0.431 control: 0 region: 87322905-87324677 ef: inf
   P-2-2 chr8 87323637 ChIP: 0.533 control: 3.7e-05 region: 87322905-87324677 ef: 17570.083
   P-2-3 chr8 87323991 ChIP: 0.474 control: 0.00356 region: 87322905-87324677 ef: 162.397
   P-2-4 chr8 87324464 ChIP: 0.292 control: 0.0069 region: 87322905-87324677 ef: 51.616

   In this example two regions are reported one containing a
   single local peak and another containing 4 local peaks.

   calls/peak_caller.ChIP.out.with_metrics

   This file contains the same entries as the peak_caller.ChIP.out file
   with additional metrics characterizing regions and peaks.

   Definition:

   R-<region_id> <chromosome> <region start>-<region end> ChIP: <ChIP enrichment intensity of the strongest peak> control: <control enrichment intensity at the position of the strongest ChIP peak> max_pos: <coordinate of the highest peak in the region> ef: <normalized enrichment fold at the position of the strongest peak> ChIP_tags: <ChIP reads starting within the region> background_tags: <background reads starting within the region> tag_ef: <normalized tag enrichment fold within the region> ps: <local peak shift estimate> cor: <local strand correlation value> os: <local oversampling value> z_score: <z-value for deviating from the Poisson expectation> n_lg_qv: <negative log_10 of Poisson p-value corrected for the number of tests>
   P-<region_id>-<peak_id> <chromosome> <peak coordinate> ChIP: <ChIP enrichment intensity> control: <control enrichment intensity> region: <region start>-<region end> ef: <normalized_enrichment fold> ChIP_tags: <ChIP reads starting within the region> background_tags: <background reads starting within the region> tag_ef: <normalized tag enrichment fold within the region> ps: <local peak shift estimate> cor: <local strand correlation value> os: <local oversampling value> z_score: NA n_lg_qv: <negative log_10 of approximated FDR>

   Example:

   R-1 chr8 123640927-123641979 ChIP: 0.58 control: 0.0027 max_pos: 123641348 ef: 262.007 ChIP_tags: 660 background_tags: 16 tag_ef: 50.3 ps: 34 cor: 0.824 os: 1.54 z_score: 178.603 n_lg_qv: 40.731
   P-1-1 chr8 123641348 ChIP: 0.58 control: 0.0027 region: 123640927-123641979 ef: 262.007 ChIP_tags: 267 background_tags: 6 tag_ef: 54.3 ps: 34 cor: 0.832 os: 1.54 z_score: NA n_lg_qv: 2.118

   R-2 chr8 87322905-87324677 ChIP: 0.533 control: 3.7e-05 max_pos: 87323637 ef: 17570.083 ChIP_tags: 717 background_tags: 6 tag_ef: 146 ps: 34 cor: 0.754 os: 1.56 z_score: 321.054 n_lg_qv: 40.731
   P-2-1 chr8 87323192 ChIP: 0.431 control: 0 region: 87322905-87324677 ef: inf ChIP_tags: 185 background_tags: 0 tag_ef: inf ps: 45 cor: 0.905 os: 1.68 z_score: NA n_lg_qv: 2.109
   P-2-2 chr8 87323637 ChIP: 0.533 control: 3.7e-05 region: 87322905-87324677 ef: 17570.083 ChIP_tags: 252 background_tags: 0 tag_ef: inf ps: 34 cor: 0.745 os: 1.49 z_score: NA n_lg_qv: 2.115
   P-2-3 chr8 87323991 ChIP: 0.474 control: 0.00356 region: 87322905-87324677 ef: 162.397 ChIP_tags: 221 background_tags: 1 tag_ef: 270 ps: -8 cor: 0.665 os: 1.42 z_score: NA n_lg_qv: 2.111
   P-2-4 chr8 87324464 ChIP: 0.292 control: 0.0069 region: 87322905-87324677 ef: 51.616 ChIP_tags: 127 background_tags: 2 tag_ef: 77.5 ps: 61 cor: 0.768 os: 1.57 z_score: NA n_lg_qv: 2.101

   calls/peak_caller.ChIP.out.accepted, calls/peak_caller.ChIP.out.rejected
   
   These files contain entries in the same format as 
   calls/peak_caller.ChIP.out.with_metrics with additional filters applied
   based on the peak shift and oversampling values. Default filter values
   are os <= 5 and ps>= 0. .accpeted file will contain entries passing the 
   filter and .rejected file will contain regions and peaks for which 
   the filter didn't pass either for a region or for any one peak within the
   region.
   

2. QuEST summary file module_outputs/QuEST.out contains information about
   the number of identified regions and peaks and FDR if applicable.

   Example:

   ChIP peaks: 2510
   ChIP peaks accepted: 1762
   ChIP peaks rejected: 748

   ChIP regions: 2319
   ChIP regions accepted: 1767
   ChIP regions rejected: 552

3. UCSC Genome Browser track files

   tracks/ChIP.filtered.bed
   
   Track file displaying regions and individual peak positions.
   
   tracks/data_bed_files
   
   directory containing bed files for both ChIP and RX_noIP data
   for the entire genome as well as for individual chromosomes
   
   tracks/wig_profiles

   directory containing wig files for both normalized and 
   unnormalized CDP profiles for ChIP and background for the entire
   genome and by chromosome.

==============================================================
VIII. QuEST software overview
==============================================================

QuEST performs analysis of ChIP and control sequencing data to
identify positions in the genome that are enriched in the ChIP data
compared to the control data. If sufficient number of the control
reads is provided, QuEST also performs a false discovery rate (FDR)
for the choice of QuEST parameters.

The QuEST workflow consists of several modules that are implemented as 
stand-alone executables consolidated into the QuEST pipeline using 
wrapper and master scripts.

QuEST starts with the generate_QuEST_parameters.pl script which 
converts the data to usable formats, creates parameter files and 
defines necessary thresholds. Specifically, the sequencing data are
converted into QuEST alignment formats, the genome is formatted into 
a genome table, and the control data is split into background and 
pseudo ChIP data according to user’s needs. Parameter files are
created for individual QuEST modules that contain all the necessary 
input/output definitions and thresholds for the modules to be executed.

The QuEST master script is called run_QuEST_with_param_file.pl. It 
performs the QuEST analysis by executing necessary modules according
to the user specifications.

Modules are implemented by wrappers that run executables by 
specifying parameters from the parameter files. Quick Window Scan
and Calibrate Peak Shift Modules implement the parameter estimation
part of QuEST.  The modules CDP and Peak Caller implement the peak 
finding part of QuEST.  The FDR estimation is implemented by 
run_QuEST_with_param_file.pl based on the output of peak calling for 
ChIP and pseudo ChIP data.

Here is a brief description of how QuEST works:

1. Peak shift estimation. 
To estimate peak shift, QuEST initially selects candidate regions
to perform this estimation. Initial candidates are identified using
Quick Window Scan module, that is a sliding window that finds positions
where a window contains a lot of tags and is enriched compared to control.
Regions are then ranked by the tag count and the top 200 hundred after
some light filtering is used for peak shift estimation. For each region
the peak shift estimate is given as the distance for which forward and
reverse tag density profile show maximum of correlation after being shifted
towards each other by that distance. Global peak shift estimate is then
given by the mode of local peak shift estimates across 200 regions. This mode
is calculated after applying kernel density estimation to the 200 sample
points (bw=10).

2. CDP calculation (controlled by 
   parameters/generate_profile.ChIP.batch.pars, 
   parameters/generate_profile.background.batch.pars)

Internally, QuEST calculates enrichment profiles (also referred to as CDP scores)
for both ChIP and control data that quantify enrichment at every bp in the genome. 
CDP profiles H(i) represent kernel density estimates of tag distribution and are 
given by the formula

H_plus(i) = Sum(i-3*h < j < i+3*h){ K( (j-i)/h ) * C_plus(j) / h }
H_minus(i) = Sum(i-3*h < j < i+3*h){ K( (j-i)/h ) * C_minus(j) / h }

H(i) = H_plus(i-lambda) + H_minus(i+lambda)

where

 K(x) is a Gaussian kernel K(x) = exp( -x*x/2) / (2*pi)^{0.5} ,
 C_plus(j) and C_minus(j) represent numbers of 5' read starts mapped
 to the position j on + or - of the reference genome respectively, 
 h is a kernel density bandwidth,
 lambda is a peak shift.

3. Density normalization

Although internally QuEST operates with raw CDP values H(i), peak calling stringency
is controlled by fold enrichment thresholds which are relevant to the normalized CDP values
H_norm(i). Internally, QuEST converts fold enrichment thresholds to CDP thresholds
by scaling them according to the number of tags in the epxeriment, the size of the genome 
and its mappable fraction. This is equivalent to applying fold enrichment thresholds to 
the normalized CDP values given by:

H_norm(i) = H(i) * (G*M/T)

where 
      G - the size of the genome
      M - the mappable fraction of the genome
      T - the number of effectively used mapped reads in the experiment

4. Peak calling (controlled by parameters/peak_caller.ChIP.batch.pars)

Peak caller first finds "seeds" and then extends them in both directions
to obtain regions of enrichment.

a "seed" at the position i must satisfy the following set of conditions:

     a. H_ChIP_norm(i) > ChIP_fold_enrichment_threshold
     b. H_ChIP_norm(i) > ChIP_to_background_fold_enrichment_threshold * H_RX_noIP_norm(i) 

Condition (b) is not used in the absense of control data.

For each seed i, the region of enrichment is defined as the interval [l,r] that contains 
seed i such that for every j from [l,r] the following enrichment condition is sastisfied:

    H_ChIP_norm(j) > ChIP_fold_enrichment_extension_threshold

All overlapping regions are then merged and peaks k within regions are identified
as positions for which the following conditions are satisfied:

   a. position k is a local maximum of normalized density H_ChIP_norm(j)
   b. positions k is "enriched", i.e. H_ChIP_norm(k) > ChIP_fold_enrichment_threshold
   c. Maximum of a region satisfying (a,b) is a peak
   d. Between the current peak and the next higher peak within the same region 
   there exists a dip in the values of enrichment profile

5. Metrics (controlled by parameters/peak_filter.ChIP.batch.pars)
Metrics are calculated for each region and peak and include local peak shift
estimate, correlation of strands at the global peak shift estimate, Q-values,
fold enrichment, tag fold enrichment, and rank by Q-value.

6. Peak and region filter (controlled by parameters/peak_filter.ChIP.batch.pars)
Peak and region filter are applied to eliminate potentially artifactual peaks
and regions. Currently, peak and region filter is applied when using default
QuEST option of motif-specific TF (with punctate peak). All subordinated
peaks within the region will be eliminated (because these tend to be small
"spilling" peaks from a nearby big peak), this option is controlled by report_peaks = best.
All peaks with a "wrong" value of peak shift will be removed (controlled by 
ps_lower_threshold = 25). All regions containing no peaks after peak filter is
applied will be removed (controlled by filter_out_empty_regions = true).


Filter applies only in some cases (such as punctate peak option) and includes
filter on local peak shift () to remove repeat mismapping effects, 
best peak per region, removing regions that don't have any peaks after filtering

7. P-values and Q-values
P-values for regions are calculated based on the number of reads observed within
a region of enrichment under Poisson expected distribution. Poisson lambda is
estimated per bp adaptively using control data in a narrow (~300 bp), medium (~3Kb)
and large (~10 Kb) region centered at the mid point of a given region. The largest
of the 3 lambdas and uniform model lambda is selected to provide a conservative p-value
estimate.

P-values for peaks are calculated based on a normal z-score of a kernel density estimation
statistic. The mean of expected KDE value is given by expected tags/bp and variance 
is given by expected tags/bp  times the sum of square of a kernel function. If significant
deviation from the expected tags/bp is observed in the control data at the peak
position, control data is instead used to estimate expected tags/bp.

Q-values are calculated using Benjamini-Hochberg correction of p-values

8. Stack collapsing. Mapping and PCR artifacts are often present in ChIP-seq samples
and can result in false peaks. To deal with this issue, QuEST implements collapsing
of sparse read stacks. Read stack is 2 or more reads that have the same 5' read start
on the same strand (don't confuse this with staggered reads). By default, QuEST will 
identify all stacks and then calculate a p-value of each stack using binomial distribution of reads 
obtained from a 100 bp window centered around a stack. If the p-value is less than 0.001,
and the stack occurs in a 100 bp window with less than 30% of positions having read starts,
the stack is reduced to a single instance of a read. More control over the parameters
of stack collapsing can be obtained using -advanced flag during QuEST parametrization.
Data bed files produced by QuEST will contain data after the stack collapsing is applied.

9. Visualizing enrichment profiles (controlled by parameters/ChIP_normalized.wig.batch.pars)

Normalized enrichment profiles representing normalized CDP values H_norm(i) can
be visualized through the UCSC genome browser custom tracks. By default, QuEST will only
write wig values for positions with fold enrichment values exceeding 10. In case of weak
signal this may result in poor visualization of results. In this case it is recommended
to lower the display cutoff (controlled by profile_threshold parameter) and re-run QuEST
portion that regenerates wig profiles. This is done changing the default threshold
to a desired value, e.g 3 (profile_threshold = 3), and then running run_QuEST_with_param_file.pl
with flags -generate_profile_ChIP -wig_tracks_ChIP.

==============================================================
IX. Explanation of QuEST flags
==============================================================

generate_QuEST_parameters.pl

ChIP data:
     -QuEST_align_ChIP (QuEST alignment format)
     -solexa_align_ChIP (the very first solexa alignment format)
     -eland_align_ChIP (eland output format)
     -bowtie_align_ChIP (bowtie native format)
     -maq_align_ChIP (maq native format)
     -sam_align_ChIP (sam format, by bowtie or maq)	

RX_noIP data, aka input, total chromatin, background, control:
     -QuEST_align_RX_noIP (QuEST alignment format)
     -solexa_align_RX_noIP (the very first solexa alignment format)
     -eland_align_RX_noIP (eland output format)
     -bowtie_align_RX_noIP (bowtie native format)
     -maq_align_RX_noIP (maq native format)
     -maq_align_RX_noIP (sam format, by bowtie or maq)	

==============================================================
X. QuEST parametrization
==============================================================

1. bandwidth
   Bandwidth parameter affects the amount of smoothing applied
   to the tag map data during calculation of enrichment profiles
   (CDP). Larger bandwidth results in more smoothing, more significant
   p-values at the cost of peak precision and reduced enrichment
   values. When working with punctate peaks, bandwidth of 30 is 
   recommended. Histone IPs are more complicated and usually require 
   more smoothing. How much? That depends on the
   type of histone mod, cell type, antibody, etc. Values between
   100 and 1000 bp usually work pretty well.

==============================================================
XI. FAQ
==============================================================

1. Q: When trying to load bed/wig files into the genome browser
      the browser just hangs there and never displays the tracks. What
      should I do?

   A: UCSC genome browser currently has a set limit on the number of
      features that can be displayed in the custom track (typically about
      10 M for the bed file). Try using individual chromosome file
      to overcome this limitation.

2. Q: I've only got 300 peaks after running QuEST with recommended (low)
      stringency, although I was expecting more than that. Can I trust my results?

   A: It depends, and one needs to evaluate a combination
      of factors here. Are the right motifs associated with those sites?
      Do the genes that are bound make sense in the context of a protein function?
      Do the sites validate with qPCR? (this has to be done on ChIP DNA, not the library)
      If this is a punctate-site type protein, do the sites have the right
      tag distribution (i.e. + and - tags cluster on opposite sides of the site)?
      Typically low number of sites is indicative of low quality of antibody or
      that ChIP didn't work. If you see consistenly low enrichment values at the
      binding sites that's indicative of suboptimally performing antibody.

