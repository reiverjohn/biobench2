#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <list>
#include <math.h>
#include <algorithm>

#include <assert.h>
#include <time.h>
#include <pthread.h>

#include "string_utils.h"
#include "seq_contig.h"
#include "params.h"

using namespace std;

// this program takes eland, solexa or QuEST align files,
// compares them to the genome table, and stores in 
// QuEST format

struct contig{
  string name;
  int size;
};

class hit{
public:
  bool matched; //true if alignment to anything is reported false otherwise
  string chrom;
  char strand;
  int coord;   // 5' coordinate of the read
  int mismatches;
  bool gt_match;
  
  bool offending_string;

  //int ind; //index of this read (to locate alignments for the same reads scattered across the alignment file)
  string alignment_id;
  

  hit(){
    chrom = "na";
    strand = 'n';
    coord = -1;
    mismatches = -1;
    gt_match = false;
    offending_string = false;
  }

  hit& operator= (const hit& h){
    chrom = h.chrom;
    strand = h.strand;
    coord = h.coord;
    gt_match = h.gt_match;
    offending_string = h.offending_string;
    //ind = h.ind;
    alignment_id = h.alignment_id;

    return *this;
  }

  void import_sam_alignment(string& alignment_string, ofstream& error_str, bool using_error_stream, vector<contig>& contigs){
    char tab_symbol = char(9);

    if(alignment_string.length() > 0){
      if(alignment_string[0] != '@'){
	vector<string> sam_alignment_fields = split(alignment_string, tab_symbol);
	if(sam_alignment_fields.size() >= 9){
	  unsigned short flag = (unsigned short) atoi(sam_alignment_fields[1].c_str());
	  //cout<<endl;
	  //cout<<alignment_string<<endl;
	  //cout<<"flag: "<<flag<<endl;
	  //cout<<(flag & 0x0004)<<endl;
	  //exit(0);
	  if( (flag & 0x0004) == 0x0004){
	    //cout<<"here"<<endl;
	    //exit(0);
	    //unmapped
	    matched = false;
	  }
	  else{ //mapped
	    matched = true;
	    if( (flag & 0x0010) == 0x0010){
	      strand = '-';
	    }
	    else{
	      strand = '+';
	    }
	    chrom = sam_alignment_fields[2];
	    if(strand == '+'){
	      coord = atoi(sam_alignment_fields[3].c_str());	      
	    }
	    else{
	      coord = atoi(sam_alignment_fields[3].c_str()) + sam_alignment_fields[9].length() - 1;
	    }
	    
	    gt_match = false;
	    for(unsigned int i=0; i<contigs.size(); i++){
	      if(contigs[i].name == chrom){
		if(contigs[i].size > coord){
		  gt_match = true;
		}
		else{
		  if(using_error_stream){
		    error_str<<"Alignment: "<<alignment_string<<endl;
		    error_str<<"Coordinate "<<coord<<" out of range: [0,"<<contigs[i].size<<")"<<endl;
		    error_str<<endl;
		  }
		}
	      }
	    }

	    if(gt_match == false){
	      //cout<<"Alginment: "<<alignment_string<<endl;
	      if(using_error_stream){
		error_str<<"Alignment: "<<alignment_string<<endl;
		error_str<<"Failed to find contig "<<chrom<<" in the genome table"<<endl;
		
	      }
	    }
	    else{
	      //cout<<chrom<<" "<<coord<<" "<<strand<<endl;
	      //exit(0);
	    }

	  }
	}
	else{
	  offending_string = true;	 
	  if(using_error_stream){
	    error_str<<"Alignment: "<<alignment_string<<endl;
	    error_str<<"Expected 12 or more fields but found "<<sam_alignment_fields.size()<<endl;
	    error_str<<endl;
	  }
	}
      }
    }
  }

  void import_bowtie_alignment(string& alignment_string, ofstream& error_str, bool using_error_stream){
    char tab_symbol = char(9);
    char comma_char = ',';

    vector<string> alignment_string_fields = split(alignment_string, tab_symbol);
    if(alignment_string_fields.size() >=7){
      //ind = atoi(alignment_string_fields[0].c_str());
      alignment_id = alignment_string_fields[0];

      if(alignment_string_fields[1] == "+" || alignment_string_fields[1] == "-"){
	chrom = alignment_string_fields[2];
	int _coord = atoi(alignment_string_fields[3].c_str());
	string read = alignment_string_fields[4];
	if(alignment_string_fields[1] == "+"){
	  strand = '+';
	  coord = _coord;
	}
	else{
	  strand = '-';
	  coord = _coord + read.length() -1;
	}

	if(alignment_string_fields.size() >= 8){
	  string mismatch_string = alignment_string_fields[7];
	  vector<string> mismatch_string_fields = split(mismatch_string, comma_char);
	  mismatches = mismatch_string_fields.size();
	}
	else{
	  mismatches = 0;
	}
	offending_string = false;
      }
      else{
	offending_string = true;
	if(using_error_stream){
	  error_str<<"Alignment: "<<alignment_string<<endl;
	  error_str<<"Expected +/- in the 2nd field, but found: "<<alignment_string_fields[1]<<endl;
	  error_str<<endl;
	  return;
	}
      }
      //cout<<"alignment ok"<<endl;
    }
    else{
      offending_string = true;
      if(using_error_stream){
	error_str<<"Alignment: "<<alignment_string<<endl;
	error_str<<"expected at least 8 tab-delimited fields, but found "<<alignment_string_fields.size()<<endl;
	error_str<<endl;
	return;
      }
    }
  }

  void import_eland_extended_alignment(string& read, string& chrom_fname, string& hit_string, vector<contig>& genome_table){     //for eland_extended
    char dot_symbol = '.';
    
    int read_length = read.length();
    gt_match = false;
    vector<string> chrom_filename_fields = split(chrom_fname, dot_symbol);
    //string chrom;
    //cout<<"Parsing the hit string"<<endl;
    if(chrom_filename_fields.size() <= 2){
      string parsed_chrom = chrom_filename_fields[0];
      
      bool chrom_matched = false;
      int matched_chrom_ind = -1;
      for(unsigned int i=0; i<genome_table.size() && !chrom_matched; i++){
	if(genome_table[i].name == parsed_chrom){
	  chrom_matched = true;
	  chrom = parsed_chrom;
	  matched_chrom_ind = i;
	}
      }
      if(chrom_matched){
	gt_match = true;	
      }
      else{
	chrom = parsed_chrom;
	coord = -1;
      }
      
      int F_or_R_ind = -1;      
      char _strand = '*';
      for(unsigned int i=0; i<hit_string.size(); i++){
	if(hit_string[i] == 'F' || hit_string[i] == 'R'){
	  F_or_R_ind = i;
	  _strand = hit_string[i];
	}
      }
      
      if(F_or_R_ind >= 0){
	
	string coord_string = hit_string.substr(0, F_or_R_ind);
	string edit_string = hit_string.substr(F_or_R_ind+1, hit_string.length()-F_or_R_ind-1);
	
	int _coord = atoi(coord_string.c_str());
	if(_strand == 'F'){
	  coord = _coord - 1;	  
	  strand = '+';
	}
	else{
	  coord = _coord - 1 + (read_length-1);
	  strand = '-';
	}

	if(chrom_matched){
	  if(coord >= genome_table[matched_chrom_ind].size){
	    gt_match = false;
	    cout<<"Warning: read 5 prime end out of boundaries. Calculated coord: "<<coord<<" contig size: "<<genome_table[matched_chrom_ind].size<<endl;
	  }
	}
	//	cout<<"Hit: "<<chrom<<" "<<strand<<" "<<coord<<" edit string: "<<edit_string<<endl;
	
	//now parse the edit string
	mismatches= 0;
	for(unsigned int i=0; i<edit_string.size(); i++){
	  char this_edit_string_char = edit_string[i];
	  if(this_edit_string_char == 'A' || this_edit_string_char == 'C' ||
	     this_edit_string_char == 'G' || this_edit_string_char == 'T'){
	    mismatches++;
	  }
	}
	
	//cout<<"Hit: "<<chrom<<" "<<strand<<" "<<coord<<" edit string: "<<edit_string<<" mm: "<<mismatches<<endl;
      }
    }
  }
  ~hit(){};
};

int main(int argc, char* argv[]){
  
  params pars(argc, argv);
  pars.require("align_file","align_file",STRING_TYPE);
  pars.require("output_file", "output_file", STRING_TYPE);
  pars.require("genome_table", "genome_table", STRING_TYPE);
  pars.require("align_type", "eland/eland_extended/QuEST/solexa_align/bowtie/maq/sam", STRING_TYPE);
  pars.optional("report_file","report_file","none",STRING_TYPE);
  pars.optional("error_report_file","error_report_file","none",STRING_TYPE);

  if(!pars.enforce()){
    exit(1);
  }

  //cout<<endl;
  //pars.list_all_params();

  string align_fname = pars.get_string_value("align_file");
  string output_fname = pars.get_string_value("output_file");  
  string genome_table_fname = pars.get_string_value("genome_table");
  string report_fname = pars.get_string_value("report_file");
  string error_report_fname = pars.get_string_value("error_report_file");

  string align_type = pars.get_string_value("align_type");
 
  if(align_type != "eland" && align_type != "QuEST" && align_type != "solexa_align" && align_type != "eland_extended" && align_type != "bowtie" && align_type != "maq" && align_type != "sam"){
    cerr<<"Bad align type. Expected eland/eland_extended/QuEST/solexa_align/bowtie/maq/sam but found "<<align_type<<". Aborting."<<endl;
    exit(4);
  }


  cout<<"align_file:         "<<align_fname<<endl;
  cout<<"output_file:        "<<output_fname<<endl;
  cout<<"genome_table:       "<<genome_table_fname<<endl;
  cout<<"alignment_type:     "<<align_type<<endl;
  cout<<"report_file:        "<<report_fname<<endl;
  cout<<"error_report_file:  "<<error_report_fname<<endl;
  cout<<endl;

  ifstream genome_table_str(genome_table_fname.c_str());
  if(!genome_table_str.good()){
    cerr<<"bad file name: "<<genome_table_fname<<endl;
    exit(1);
  }

  bool using_report_file;
  if(report_fname == "none"){
    using_report_file = false;
  }
  else{
    using_report_file = true;
  }

  bool using_error_report_file;
  if(error_report_fname == "none"){
    using_error_report_file = false;
  }
  else{
    using_error_report_file = true;
  }

  ofstream report_str;
  if(using_report_file){
    report_str.open(report_fname.c_str());
    
    if(!report_str.good()){
      cerr<<"bad file name: "<<report_fname<<endl;
      exit(1);
    }
  }
  ofstream error_report_str;
  if(using_error_report_file){
    error_report_str.open(error_report_fname.c_str());
    
    if(!error_report_str.good()){
      cerr<<"bad file name: "<<error_report_fname<<endl;
      exit(1);
    }
  }


  //vector<string> contig_names;
  //vector<int> contig_sizes;

  vector<contig> contigs;

  //char dot_symbol = '.';
  char comma_symbol = ',';
  char gap_symbol = ' ';
  char tab_symbol = char(9);
  char column_symbol = ':';

  while(genome_table_str.good()){
    string dummy_string;
    getline(genome_table_str, dummy_string);
    if(genome_table_str.good()){
      if(dummy_string.length() > 0){
	if(dummy_string[0] != '#'){
	  vector<string> cur_contig_fields = split(dummy_string, gap_symbol);
	  if(cur_contig_fields.size() >= 2){
	    string cur_contig_name = cur_contig_fields[0];
	    int cur_contig_size  = atoi(cur_contig_fields[1].c_str());
	    
	    //assert(cur_contig_size >= 0);
	    if(cur_contig_size <= 0){
	      cerr<<"Bad contig size "<<cur_contig_size<<" in "<<dummy_string<<". Skipping."<<endl;	      
	    }
	    else{
	      //check for gaps and tabs in the contig name
	      bool symbols_replaced_in_the_contig_name = false;
	      for(unsigned int i=0; i<cur_contig_name.length(); i++){
		if(cur_contig_name[i] == gap_symbol || cur_contig_name[i] == tab_symbol){
		  cur_contig_name[i] = '_';
		  symbols_replaced_in_the_contig_name = true;
		}
	      }

	      if(symbols_replaced_in_the_contig_name){
		cerr<<"some gap and tab symbols were replaced by _ in the contig name $cur_contig_fields[0]"<<endl;
	      }
	      
	      bool contig_name_already_exists = false;
	      for(unsigned int i=0; i<contigs.size(); i++){
		if(contigs[i].name == cur_contig_name){
		  contig_name_already_exists = true;
		}
	      }

	      if(!contig_name_already_exists){
		contig new_contig;
		new_contig.name = cur_contig_name;
		new_contig.size = cur_contig_size;
		contigs.push_back(new_contig);
	      }
	    }
	  }
	}
      }
    }
  }

  genome_table_str.close();  

  //checking for duplicates

  int line_counter = 0;
  int reads_matched = 0;
  int reads_not_matched = 0;
  int reads_not_in_genome_table = 0;
  int offending_strings = 0;

  remove(output_fname.c_str());
  ofstream out_str(output_fname.c_str());
  if(!out_str.good()){
    cerr<<"bad file name: "<<output_fname<<endl;
    exit(1);
  }

  ifstream align_str(align_fname.c_str());
  if(!align_str.good()){
    cerr<<"Failed to read the file "<<align_fname<<endl;
    exit(1);
  }

  string dummy_string;
  
  stringstream dummy_stream;

  if(align_type == "QuEST"){
    while(align_str.good()){
      getline(align_str, dummy_string);
      if(line_counter % 10000 == 0){
	printf ("\rreads read: %.2f M, matched : %.2f M", ((double)line_counter/1000000.0), ((double)reads_matched/1000000.0) );
	cout.flush();				      
      }
      line_counter++;
      
      if(align_str.good()){
	if(dummy_string.length() > 0){
	  if( dummy_string[0] != '#'){    
	    
	    vector<string> cur_hit_entry_fields = split(dummy_string, gap_symbol);
	    if(cur_hit_entry_fields.size() >= 3){

	      string cur_hit_contig_name = cur_hit_entry_fields[0];
	      int cur_hit_5p_coord = atoi(cur_hit_entry_fields[1].c_str());
	      string cur_hit_orient = cur_hit_entry_fields[2];
	      
	      if(cur_hit_orient != "+" && cur_hit_orient != "-"){
		cerr<<"Warning: expected +/- for orientation of read but found "<<cur_hit_orient;
		cerr<<". Skipping."<<endl;
	      }
	      else{	      

		bool reference_contig_found = false;
		for(unsigned int i=0; i<contigs.size(); i++){
		  if(cur_hit_contig_name == contigs[i].name){
		    reference_contig_found = true;
		    int contig_size = contigs[i].size;
		    
		    if(cur_hit_5p_coord < 0 || cur_hit_5p_coord >= contig_size){
		      cerr<<"Warning read coordinate "<<cur_hit_5p_coord<<" is out of boundaries [ 0,";
		      cerr<<contig_size<<" ]. Skipping."<<endl;		
		    }
		    else{
		      dummy_stream<<cur_hit_contig_name<<" "<<cur_hit_5p_coord<<" "<<cur_hit_orient<<endl;
		      reads_matched++;
		    }
		  }
		}
		if(!reference_contig_found){
		  if(using_error_report_file){
		    error_report_str<<"line: "<<line_counter<<" Align_string "<<dummy_string<<endl;
		    error_report_str<<"Failed to find "<<cur_hit_contig_name<<" in the genome table."<<endl;
		    error_report_str<<endl;
		  }
		}
	      }
	    }
	  }
	}   
      }
    }
    out_str<<dummy_stream.str();
  }

  if(align_type == "solexa_align"){
    while(align_str.good()){
      getline(align_str, dummy_string);
      if(line_counter % 10000 == 0){
	printf ("\raligments: %.2f M, matched : %.2f M", ((double)line_counter/1000000.0), ((double)reads_matched/1000000.0 ) );
	cout.flush();				      
      }
      line_counter++;
      
      if(align_str.good()){
	if(dummy_string.length() > 0){
	  if( dummy_string[0] != '#'){    
	    
	    vector<string> cur_hit_entry_fields = split(dummy_string, gap_symbol);
	    if(cur_hit_entry_fields.size() >= 5){
	      string cur_hit_pos = cur_hit_entry_fields[3];
	      vector<string> cur_hit_pos_fields = split(cur_hit_pos, column_symbol);
	      if(cur_hit_pos_fields.size() == 2){
		
		string cur_read = cur_hit_entry_fields[0];
		string cur_hit_contig_name = cur_hit_pos_fields[0];
		int cur_hit_5p_coord = atoi(cur_hit_pos_fields[1].c_str());
		string cur_hit_orient = cur_hit_entry_fields[4];
		
		if(cur_hit_orient != "F" && cur_hit_orient != "R"){
		  cerr<<endl;
		  cerr<<"Warning: expected F/R for orientation of read but found "<<cur_hit_orient;
		  cerr<<". Skipping."<<endl;
		}
		else{	      
		  if(cur_hit_orient == "F"){ 
		    cur_hit_orient = "+"; 
		    cur_hit_5p_coord = cur_hit_5p_coord - 1;
		  }
		  else{
		    if(cur_hit_orient == "R"){
		      cur_hit_orient = "-";
		      cur_hit_5p_coord = cur_hit_5p_coord + cur_read.length() - 2;
		    }
		  }
		  
		  bool reference_contig_found = false;
		  for(unsigned int i=0; i<contigs.size(); i++){
		    if(cur_hit_contig_name == contigs[i].name){
		      reference_contig_found = true;
		      int contig_size = contigs[i].size;
		      
		      if(cur_hit_5p_coord < 0 || cur_hit_5p_coord >= contig_size){
			cerr<<endl;
			cerr<<"for read: "<<dummy_string<<endl;
			cerr<<"Warning read coordinate "<<cur_hit_5p_coord<<" is out of boundaries [ 0,";
			cerr<<contig_size<<" ]. Skipping."<<endl;		
		      }
		      else{
			dummy_stream<<cur_hit_contig_name<<" "<<cur_hit_5p_coord<<" "<<cur_hit_orient<<endl;
			reads_matched++;
		      }
		    }
		  }
		  if(!reference_contig_found){
		    cerr<<endl;
		    cerr<<"Failed to identify reference contig in the genome table for the alignment: "<<endl;
		    cerr<<dummy_string<<endl;		    
		  }
		}
	      }
	      else{
		cerr<<"Failed to understand: "<<dummy_string<<endl;
		cerr<<"Expected "<<cur_hit_pos<<" to have 2 fields separated by column, but found "<<cur_hit_pos_fields.size()<<" fields."<<endl;
		cerr<<"Skipping the alignment."<<endl;
	      }
	    }
	  }
	}
      }
    }
    out_str<<dummy_stream.str();
  }

  if(align_type == "eland"){
    while(align_str.good()){
      getline(align_str, dummy_string);
      if(line_counter % 10000 == 0){
	printf ("\raligments: %.2f M, matched : %.2f M, not in gt: %.3f K", ((double)line_counter/1000000.0), ((double)reads_matched/1000000.0 ), ((double)reads_not_in_genome_table/1000.0));
	//printf ("\raligments: %.2f M, matched : %.2f M", ((double)line_counter/1000000.0), ((double)reads_matched/1000000.0 ) );
	cout.flush();				      
      }
      line_counter++;
      
      if(align_str.good()){
	if(dummy_string.length() > 0){
	  if( dummy_string[0] != '#'){    
	    vector<string> cur_hit_entry_fields = split(dummy_string, tab_symbol);
	    if(cur_hit_entry_fields.size() >= 9){
	      int cur_read_field_num = -1;
	      int cur_hit_contig_fname_field_num = -1;
	      int cur_unique_flag_field_num = -1;
	      int cur_hit_5p_coord_field_num = -1;
	      int cur_hit_orient_field_num = -1;

	      bool alignment_fields_understood = false;
	      if(cur_hit_entry_fields[2].substr(0,1) == "U"){ 
		cur_read_field_num = 1;
		cur_hit_contig_fname_field_num = 6;
		cur_unique_flag_field_num = 2;
		cur_hit_5p_coord_field_num = 7;
		cur_hit_orient_field_num = 8;
		alignment_fields_understood = true;
	      }
	      else{
		if(cur_hit_entry_fields[3].substr(0,1) == "U" && cur_hit_entry_fields.size() >= 10){ 
		  cur_read_field_num = 1;
		  cur_hit_contig_fname_field_num = 7;
		  cur_unique_flag_field_num = 3;
		  cur_hit_5p_coord_field_num = 8;
		  cur_hit_orient_field_num = 9;
		  alignment_fields_understood = true;
		}
	      }
	      
	      if(!alignment_fields_understood){	      
		cerr<<"Failed to identify the field identity in the alignment: "<<endl;
		cerr<<dummy_string<<endl;
		cerr<<"Expecting Uniqueness flag in either field 3 or 4 but didn't find any."<<endl;
	      }
	      else{
		
		string cur_read = cur_hit_entry_fields[cur_read_field_num];
		string cur_hit_contig_fname = cur_hit_entry_fields[cur_hit_contig_fname_field_num];
		string cur_unique_flag = cur_hit_entry_fields[cur_unique_flag_field_num];
		int cur_hit_5p_coord = atoi(cur_hit_entry_fields[cur_hit_5p_coord_field_num].c_str());
		string cur_hit_orient = cur_hit_entry_fields[cur_hit_orient_field_num];
		
		if(cur_unique_flag.substr(0,1) == "U"){		

		  if(cur_hit_contig_fname.size() >= 3){
		    if(cur_hit_contig_fname.substr(cur_hit_contig_fname.size() - 3) == ".fa"){
		      
		      string cur_hit_contig_name = cur_hit_contig_fname.substr(0, cur_hit_contig_fname.length()-3);
		      if(cur_hit_orient != "F" && cur_hit_orient != "R"){
			cerr<<endl;
			cerr<<"Warning: expected +/- for orientation of read but found "<<cur_hit_orient;
			cerr<<". Skipping."<<endl;
		      }
		      else{	      
			if(cur_hit_orient == "F"){ 
			  cur_hit_orient = "+"; 
			  cur_hit_5p_coord = cur_hit_5p_coord - 1;
			}
			else{
			  if(cur_hit_orient == "R"){
			    cur_hit_orient = "-";
			    cur_hit_5p_coord = cur_hit_5p_coord + cur_read.length() - 2;
			  }
			}
			
			bool reference_contig_found = false;
			for(unsigned int i=0; i<contigs.size(); i++){
			  if(cur_hit_contig_name == contigs[i].name){
			    reference_contig_found = true;
			    int contig_size = contigs[i].size;
			    
			    if(cur_hit_5p_coord < 0 || cur_hit_5p_coord >= contig_size){
			      cerr<<endl;
			      cerr<<"for read: "<<dummy_string<<endl;
			      cerr<<"Warning read coordinate "<<cur_hit_5p_coord<<" is out of boundaries [ 0,";
			      cerr<<contig_size<<" ]. Skipping."<<endl;					      
			    }
			    else{
			      dummy_stream<<cur_hit_contig_name<<" "<<cur_hit_5p_coord<<" "<<cur_hit_orient<<endl;
			      reads_matched++;
			    }
			  }
			}
			if(!reference_contig_found){
			  reads_not_in_genome_table++;
			  if(using_error_report_file){
			    error_report_str<<"Line: "<<line_counter<<endl;
			    error_report_str<<"String: "<<dummy_string<<endl;
			    error_report_str<<"Failed to identify reference contig in the genome table for the alignment: "<<endl;
			    error_report_str<<endl;
			  }
			}
		      }
		    }
		    else{
		      reads_not_in_genome_table++;
		      if(using_error_report_file){
			error_report_str<<"Line: "<<line_counter<<endl;
			error_report_str<<"String: "<<dummy_string<<endl;
			error_report_str<<"Failed to understand contig file name "<<cur_hit_contig_fname<<" in this alignment. "<<endl;
			error_report_str<<"Expecting the file name to end with \".fa\" but found otherwise. Skipping"<<endl;
			error_report_str<<endl;
		      }
		    }
		  }
		  else{
		    cerr<<endl;
		    cerr<<"Failed to understand contig file name "<<cur_hit_contig_fname<<" in the alignment:"<<endl;
		    cerr<<dummy_string<<endl;
		    cerr<<"Expecting 3 or more symbols in the name, but found "<<cur_hit_contig_fname.length()<<endl;
		    cerr<<"Skipping the alignment."<<endl;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    out_str<<dummy_stream.str();
  }

  if(align_type == "eland_extended"){
    while(align_str.good()){
      getline(align_str, dummy_string);
      if(line_counter % 10000 == 0){
	printf ("\raligments: %.2f M, matched : %.2f M, not in gt: %.3f M", ((double)line_counter/1000000.0), ((double)reads_matched/1000000.0 ), ((double)reads_not_in_genome_table/1000000.0));
	cout.flush();				      
      }
      line_counter++;
      
      if(align_str.good()){
	if(dummy_string.length() > 0){
	  if( dummy_string[0] != '#'){    
	    vector<string> cur_hit_entry_fields = split(dummy_string, tab_symbol);
	    //cout<<endl<<"Alignment: "<<dummy_string<<endl;
	    if(cur_hit_entry_fields.size() >= 4 && cur_hit_entry_fields[2] != "NM"){
	      string cur_read = cur_hit_entry_fields[1];
	      //int match_length = cur_hit_entry_fields[1].length();
	      vector<string> match_string_fields = split(cur_hit_entry_fields[2],column_symbol);
	      if(match_string_fields.size() != 3){
		//cerr<<"alignment: "<<dummy_string<<endl;
		//cerr<<"match_string: "<<cur_hit_entry_fields[2]<<endl;
		//assert(match_string_fields.size() == 3);
	      }
	      else{
		int mm0 = atoi(match_string_fields[0].c_str());
		int mm1 = atoi(match_string_fields[1].c_str());
		int mm2 = atoi(match_string_fields[2].c_str());
		
		if( mm0 == 1 ||
		    ( mm0 == 0 && mm1 == 1) ||
		    ( mm0 == 0 && mm1 == 0 && mm2 == 1)){
		  
		  vector<hit> this_read_hits;
		  
		  vector<string> chrom_hits = split(cur_hit_entry_fields[3], column_symbol);
		  if(chrom_hits.size() > 0){
		    //cerr<<endl<<cur_hit_entry_fields[3]<<endl;
		    //for(unsigned int i=0; i<chrom_hits.size(); i++){
		      //cerr<<i<<" -> "<<chrom_hits[i]<<endl;
		    //}
		    string cur_chrom_fname = chrom_hits[0];
		    
		    
		    for(unsigned int i=1; i<chrom_hits.size(); i++){
		      //string cur_chrom_fname = chrom_hits[i];
		      vector<string> cur_chrom_hits = split(chrom_hits[i], comma_symbol);
		      for(unsigned int j=0; j<cur_chrom_hits.size(); j++){
			string cur_hit_string = "not defined";
			if(i<chrom_hits.size()-1){
			  if(j<cur_chrom_hits.size()-1){
			    cur_hit_string = cur_chrom_hits[j];
			    //cout<<"cur_chr: "<<cur_chrom_fname<<" , hit_string: "<<cur_hit_string<<endl;
			    //string d_string;
			    //cin>>d_string;
			  }
			  else{
			    cur_chrom_fname = cur_chrom_hits[j];
			  }
			}
			else{
			  cur_hit_string = cur_chrom_hits[j];
			  //cout<<"cur_chr:"<<cur_chrom_fname<<" , hit_string: "<<cur_hit_string<<endl;
			  //string d_string;
			  //cin>>d_string;
			}
			
			if(cur_hit_string != "not defined"){
			  hit new_hit;
			  new_hit.import_eland_extended_alignment(cur_read, cur_chrom_fname, cur_chrom_hits[j], contigs);
			  this_read_hits.push_back(new_hit);
			}
		      }
		    }
		  }
		  
		  if(this_read_hits.size() > 0){
		    int min_hit_ind = 0;
		    int min_hit_mm = this_read_hits[0].mismatches;
		    bool min_hit_unique = true;
		    for(int i=1; i<(int)this_read_hits.size(); i++){
		      if(this_read_hits[i].mismatches == min_hit_mm){
			min_hit_unique = false;
		      }
		      if(this_read_hits[i].mismatches < min_hit_mm){
			min_hit_mm = this_read_hits[i].mismatches;
			min_hit_ind = i;
			min_hit_unique = true;
		      }
		    }
		    if(min_hit_unique){
		      assert(min_hit_ind < (int)this_read_hits.size());
		      if(this_read_hits[min_hit_ind].gt_match){
			
			dummy_stream<<this_read_hits[min_hit_ind].chrom<<" ";
			dummy_stream<<this_read_hits[min_hit_ind].coord<<" ";
			dummy_stream<<this_read_hits[min_hit_ind].strand<<endl;
			reads_matched++;
		      }
		      else{
			reads_not_in_genome_table++;
			if(using_error_report_file){
			  error_report_str<<"line: "<<line_counter<<" Align_string: "<<dummy_string<<endl;
			  error_report_str<<"Failed to find "<<this_read_hits[min_hit_ind].chrom<<" in the genome table or coordinate out of boundaries."<<endl;
			  error_report_str<<endl;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    out_str<<dummy_stream.str();
  }

  if(align_type == "bowtie"){
    hit cur_best_hit;
    cur_best_hit.alignment_id = "not_defined";
    bool unique_best_hit_exists = false; 
    //vector<int> hit_inds;
    vector<string> alignment_ids;

    while(align_str.good()){
      getline(align_str, dummy_string);
      if(line_counter % 10000 == 0){
	printf ("\raligments: %.2f M, matched : %.2f M, not in gt: %.3f M, offending: %.3f M", ((double)line_counter/1000000.0), ((double)reads_matched/1000000.0 ), ((double)reads_not_in_genome_table/1000000.0), ((double) offending_strings/1000000.0)) ;
	cout.flush();				      
      }
      //line_counter++;
      
      if(align_str.good()){
	line_counter++;
	if(dummy_string.length() > 0){
	  if( dummy_string[0] != '#'){    
	    vector<string> cur_alignment_fields = split(dummy_string, tab_symbol);
	    if(cur_alignment_fields.size() >= 7){
	      //cout<<"importing"<<endl;
	      hit new_hit;
	      new_hit.import_bowtie_alignment(dummy_string, error_report_str, using_error_report_file);
	      
	      //cout<<new_hit.ind<<" best: "<<cur_best_hit.ind<<endl;
	      

	      if(new_hit.offending_string){
		offending_strings++;
	      }
	      else{
		//if(cur_best_hit.ind == new_hit.ind){
		if(cur_best_hit.alignment_id == new_hit.alignment_id){
		  //cout<<"same"<<endl;
		  // same read, choose a better alignment
		  if(cur_best_hit.mismatches > new_hit.mismatches){
		    // this is a better alignemnt, so keep
		    cur_best_hit = new_hit;
		    unique_best_hit_exists = true;
		  }
		  else{
		    if(cur_best_hit.mismatches == new_hit.mismatches){
		      unique_best_hit_exists = false;
		    }
		  }
		}
		else{		  
		  //cout<<"new"<<endl;
		  // this is a new read, time to save the best alignment for the previous read
		  if(unique_best_hit_exists){
		    //cout<<"saved"<<endl;
		    reads_matched++;
		    dummy_stream<<cur_best_hit.chrom<<" ";
		    dummy_stream<<cur_best_hit.coord<<" ";
		    dummy_stream<<cur_best_hit.strand<<endl;

		    //hit_inds.push_back(cur_best_hit.ind);
		    alignment_ids.push_back(cur_best_hit.alignment_id);
		  }
		  //else{
		  //  cout<<"not saved"<<endl;
		  //}
		  cur_best_hit = new_hit;
		  unique_best_hit_exists = true;
		}
	      }
	    }
	    else{
	      offending_strings++;
	      if(using_error_report_file){
		error_report_str<<"Line: "<<line_counter<<" string: "<<dummy_string<<endl;
		error_report_str<<"Too few fields"<<endl;
		error_report_str<<endl;
	      }
	      //cout<<"too few fields: "<<endl;
	      //cout<<dummy_string<<endl;
	      //cout<<"fields: "<<cur_alignment_fields.size()<<endl;
	    }
	  }
	}
      }
      else{
	// save the last alignment
	if(unique_best_hit_exists){
	  reads_matched++;
	  dummy_stream<<cur_best_hit.chrom<<" ";
	  dummy_stream<<cur_best_hit.coord<<" ";
	  dummy_stream<<cur_best_hit.strand<<endl;
	  //hit_inds.push_back(cur_best_hit.ind);
	  alignment_ids.push_back(cur_best_hit.alignment_id);
	}
      }
    }
    out_str<<dummy_stream.str();
    cout<<endl<<endl;

    cout<<"Lines read: "<<line_counter<<endl;
    cout<<"Alignments matched: "<<reads_matched<<endl;
    cout<<"Offending lines: "<<offending_strings<<endl;
    cout<<endl;

    cout<<"Checking for duplicate hits..."<<endl;
    cout<<"[ This screens for contiguity of bowtie alignment order ]"<<endl;

    sort(alignment_ids.begin(), alignment_ids.end());


    if(alignment_ids.size() > 0){

      string prev_alignment_id = alignment_ids[0];
      bool same_ind_found = false;
      for(unsigned int i=1; i<alignment_ids.size(); i++){
	if(prev_alignment_id == alignment_ids[i]){
	  same_ind_found = true;
	  if(using_error_report_file){
	    error_report_str<<"Found multiple instances of"<<prev_alignment_id<<" that were out of order"<<endl;
	  }
	}
	else{
	  //  prev_ind = hit_inds[i];
	  prev_alignment_id = alignment_ids[i];
	}
      }
      if(same_ind_found){
	cerr<<"It seems that the bowtie alignment file was out of order."<<endl;
	cerr<<"This usually is the case when -z flag is used when running Bowtie."<<endl;
	cerr<<"Currently out-of-order alignments are not supported by QuEST."<<endl;
	exit(1);
      }
      else{
	cerr<<"Your bowtie alignment file seems to be ok."<<endl;
      }

    }
    
    else{
      cout<<"Seems like there are no alignments to check. Are you sure your bowtie file is ok?\n";
    }
  }

  if(align_type == "sam"){
    while(align_str.good()){
      getline(align_str, dummy_string);
      if(line_counter % 10000 == 0){
	printf ("\raligments: %.2f M, matched : %.2f M, no match: %.2f M, not in gt: %.3f M, offending: %.3f M", ((double)line_counter/1000000.0), ((double)reads_matched/1000000.0 ), ((double)reads_not_matched/1000000.0), ((double)reads_not_in_genome_table/1000000.0), ((double) offending_strings/1000000.0)) ;
	cout.flush();				      
      }
      
      if(align_str.good()){
	line_counter++;
	if(dummy_string.length() > 0){
	  if( dummy_string[0] != '@' && dummy_string[0] != '#'){    

	    hit new_hit;
	    new_hit.import_sam_alignment(dummy_string, error_report_str, using_error_report_file, contigs);
	    
	    if(new_hit.offending_string){
	      offending_strings++;
	    }
	    else{
	      if(new_hit.matched){
		//reads_matched++;
		if(new_hit.gt_match){
		  reads_matched++;
		  dummy_stream<<new_hit.chrom<<" ";
		  dummy_stream<<new_hit.coord<<" ";
		  dummy_stream<<new_hit.strand<<endl;
		}
		else{
		  reads_not_in_genome_table++;
		}
	      }
	      else{
		reads_not_matched++;
	      }
	    }
	  }
	}
      }
    }
    
    out_str<<dummy_stream.str();
    cout<<endl<<endl;

    cout<<"Lines read: "<<line_counter<<endl;
    cout<<"Alignments matched: "<<reads_matched<<endl;
    cout<<"Offending lines: "<<offending_strings<<endl;
    cout<<endl;
  }
	    

  if(align_type == "maq"){
    while(align_str.good()){
      getline(align_str, dummy_string);
      if(line_counter % 10000 == 0){
	printf ("\raligments: %.2f M, matched : %.2f M, not in gt: %.3f M, offending: %.3f M", ((double)line_counter/1000000.0), ((double)reads_matched/1000000.0 ), ((double)reads_not_in_genome_table/1000000.0), ((double) offending_strings/1000000.0)) ;
	cout.flush();				      
      }
      
      if(align_str.good()){
	line_counter++;
	if(dummy_string.length() > 0){
	  if( dummy_string[0] != '#'){    
	    vector<string> cur_alignment_fields = split(dummy_string, tab_symbol);
	    if(cur_alignment_fields.size() >= 15){
	      string chrom = cur_alignment_fields[1];	      
	      int _coord = atoi(cur_alignment_fields[2].c_str());
	      int contig_ind = -1;

	      bool chrom_matched = false;
	      for(unsigned int i=0; i<contigs.size() && !chrom_matched; i++){
		if(contigs[i].name == chrom){
		  contig_ind = (int) i;
		  chrom_matched = true;
		}
	      }
	      //	      if(contig_ind >= 0 && contig_ind < contigs.size()){
	      if(chrom_matched){
		string strand = cur_alignment_fields[3];
		int coord = 0;
		if(strand == "+" || strand == "-"){
		  if(strand == "+"){
		    coord = _coord - 1;		    
		  }
		  else{
		    coord = (_coord-1) + cur_alignment_fields[14].length() - 1;
		  }
		  if(coord >= 0 && coord < contigs[contig_ind].size){
		    dummy_stream<<chrom<<" ";
		    dummy_stream<<coord<<" ";
		    dummy_stream<<strand<<endl;
		    reads_matched++;
		  }
		  else{
		    offending_strings++;
		    if(using_error_report_file){
		      error_report_str<<"line: "<<line_counter<<endl;
		      error_report_str<<dummy_string<<endl;
		      error_report_str<<"Coordinate "<<cur_alignment_fields[2]<<" out of contig boundaries [ 0, "<<contigs[contig_ind].size<<" ]"<<endl;
		      error_report_str<<endl;
		    }
		  }
		}
		else{
		  offending_strings++;
		  if(using_error_report_file){		      
		    error_report_str<<"line: "<<line_counter<<endl;
		    error_report_str<<dummy_string<<endl;
		    error_report_str<<"Expected +/- in the 4th field, but found: "<<strand<<endl;
		    error_report_str<<endl;
		  }
		}
	      }
	      else{
		reads_not_in_genome_table++;
		if(using_error_report_file){
		  error_report_str<<"line: "<<line_counter<<endl;
		  error_report_str<<dummy_string<<endl;
		  error_report_str<<"Failed to find "<<chrom<<" among the genome_table contigs"<<endl;
		  error_report_str<<endl;
		}
	      }
	    }
	    else{
	      offending_strings++;
	      if(using_error_report_file){
		error_report_str<<"line: "<<line_counter<<endl;
		error_report_str<<dummy_string<<endl;
		error_report_str<<"Rejected because too few fields"<<endl;
		error_report_str<<endl;
	      }
	    }
	  }
	}
      }
    }
    out_str<<dummy_stream.str();
  }
  
  align_str.close();  
  out_str.close();
  cout<<endl;

  cout<<endl;

  if(using_report_file){
    report_str<<reads_matched<<endl;
    report_str.close();
  }

  if(using_error_report_file){
    error_report_str<<reads_matched<<endl;
    error_report_str.close();
  }

  return 0;

}
