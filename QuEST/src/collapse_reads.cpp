#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>

#include <stdlib.h>

#include "string_utils.h"
#include "seq_contig.h"
#include "params.h"

using namespace std;

int int_max(int i, int j){
  if(i>j) return i;
  return j;
}

int int_min(int i, int j){
  if(i<j) return i;
  return j;
}


double Poisson_p_value(int k, double lambda, double precision ){
  //double precision = 0.00000001;
  
  double delta = 1;
  
  double cur_sum = 0;
  double prev_sum = 0;
  double cur_i = k;

  while(delta > precision){
    double log_cur_p = -lambda + cur_i*log(lambda) ;
    for(int j=1; j<=cur_i; j++){
      log_cur_p = log_cur_p - log((double)j);
    }
    
    double cur_p = exp(log_cur_p);
    cur_sum += cur_p;
    delta = cur_p;//cur_sum - prev_sum;
    
    prev_sum = cur_sum;
    cur_i++;
  }
  return cur_sum;
}

double small_binomial_p_value(int k, double p, int n){
  assert(k==2 || k==3);
  if(k==2){
    double res = 1.0 - pow(1-p, n) - ((double)n)*p*pow(1-p,n);
    return res;
  }
  if(k==3){
    double res = 1.0 - pow(1-p, n) - ((double)n)*p*pow(1-p,n-1) - ((double)n)*((double)n-1)*p*p*pow(1-p,n-2);
    return res;
  }
  assert(false);
  return -1;
}

int main(int argc, char* argv[]){
  
  params pars(argc, argv);
  pars.require("align_file","align_file",STRING_TYPE);
  pars.require("output_path", "output_path", STRING_TYPE);
  pars.require("genome_table", "genome_table", STRING_TYPE);
  pars.require("QuEST_collapsed_file","QuEST_collapsed_file",STRING_TYPE);
  pars.require("report_file","report_file",STRING_TYPE);
  
  pars.optional("new_stack_size","new_stack_size","1",INT_TYPE);
  pars.optional("collapse_reads","collapse_reads","true",STRING_TYPE);  
  pars.optional("collapse_window","collapse_window","100",INT_TYPE);
  pars.optional("count_threshold","count_threhsold","3",INT_TYPE);
  pars.optional("stack_p_value_threshold","stack_p_value_threshold","0.00001",DOUBLE_TYPE);
  pars.optional("percent_positions_hit_threshold","percent_positions_hit_threshold","30",DOUBLE_TYPE);

  if(!pars.enforce()){
    exit(1);
  }

  cout<<endl;
  pars.list_all_params();  

  string align_fname = pars.get_string_value("align_file");
  string output_path = pars.get_string_value("output_path");
  string genome_table_fname = pars.get_string_value("genome_table");

  string QuEST_collapsed_fname = pars.get_string_value("QuEST_collapsed_file");
  string report_fname = pars.get_string_value("report_file");
  
  string collapse_reads = pars.get_string_value("collapse_reads");  
  int collapse_window = pars.get_int_value("collapse_window");
  
  int count_threshold = pars.get_int_value("count_threshold");
  double stack_p_value_threshold = pars.get_double_value("stack_p_value_threshold");
  double percent_positions_hit_threshold = pars.get_double_value("percent_positions_hit_threshold");

  double p_v_precision = stack_p_value_threshold/10;
  
  int new_stack_size = pars.get_int_value("new_stack_size");
  
  /*
  if(report_fname != "not_applicable"){
    remove(report_fname.c_str());
    ofstream report_ofstr(report_fname.c_str());
    
    if(!report_ofstr.good()){
      cerr<<"Bad file name "<<report_fname<<endl;
      exit(1);
    }
  }
  */

  ofstream QuEST_collapsed_ofstr;
  if(QuEST_collapsed_fname != "not_applicable"){
    QuEST_collapsed_ofstr.open(QuEST_collapsed_fname.c_str());
    if(!QuEST_collapsed_ofstr.good()){
      cerr<<"Bad file name: "<<QuEST_collapsed_fname<<endl;
      exit(1);
    }
  }

  /*
  cout<<"align_file: "<<align_fname<<endl;
  cout<<"output_path: "<<output_path<<endl;
  cout<<"genome_table: "<<genome_table_fname<<endl;
  cout<<"collapse_reads: "<<collapse_reads<<endl;
  cout<<endl;
  */

  ifstream genome_table_str(genome_table_fname.c_str());
  if(!genome_table_str.good()){
    cerr<<"bad file name: "<<genome_table_fname<<endl;
    exit(1);
  }

  vector<string> contigs;
  vector<int> contig_sizes;

  char gap_symbol = ' ';
  while(genome_table_str.good()){
    string dummy_string;
    getline(genome_table_str, dummy_string);
    if(genome_table_str.good()){
      if(dummy_string.length() > 0){
	if(dummy_string[0] != '#'){
	  vector<string> cur_contig_fields = split(dummy_string, gap_symbol);
	  if(cur_contig_fields.size() >= 2){
	    string cur_contig_id = cur_contig_fields[0];
	    int cur_contig_size  = atoi(cur_contig_fields[1].c_str());
	    assert(cur_contig_size >= 0);

	    contigs.push_back(cur_contig_id);
	    contig_sizes.push_back(cur_contig_size);
	  }
	}
      }
    }
  }

  genome_table_str.close();
  
  vector< vector <int> > pos_hits;
  vector< vector <int> > neg_hits; 

  vector<int> dummy_int_vec;
  
  for(unsigned int i=0; i<contigs.size(); i++){
    pos_hits.push_back(dummy_int_vec);
    neg_hits.push_back(dummy_int_vec);
  }

  genome_table_str.close();

  
  int line_counter = 0;

  char sep1 = ' ';

  if(align_fname != "NA"){
    ifstream align_str(align_fname.c_str());
    if(!align_str.good()){
      cerr<<"Failed to read the file "<<align_fname<<endl;
      exit(1);
    }
    
    
    string dummy_string;
    
    while(align_str.good()){
      getline(align_str, dummy_string);
      if(line_counter % 10000 == 0){
	printf ("\rread %.2f M reads ", ((double)line_counter/1000000.0) );
	cout.flush();				      
      }
      line_counter++;
      
      if(align_str.good()){
	if(dummy_string.length() > 0){
	  if( dummy_string[0] != '#'){    
	    
	    vector<string> cur_hit_entry_fields = split(dummy_string, sep1);
	    if(cur_hit_entry_fields.size() >= 3){
	      string cur_hit_contig_name = cur_hit_entry_fields[0];
	      int cur_hit_5p_coord = atoi(cur_hit_entry_fields[1].c_str());
	      string cur_hit_orient = cur_hit_entry_fields[2];
	      
	      if(cur_hit_orient != "+" && cur_hit_orient != "-"){
		cout<<"Warning: expected +/- for orientation of read but found "<<cur_hit_orient;
		cout<<". Skipping."<<endl;
	      }
	      else{	      
		
		for(unsigned int i=0; i<contigs.size(); i++){
		  if(cur_hit_contig_name == contigs[i]){
		    int contig_size = contig_sizes[i];
		    
		    if(cur_hit_5p_coord < 0 || cur_hit_5p_coord >= contig_size){
		      cout<<"Warning read coordinate "<<cur_hit_5p_coord<<" is out of boundaries [ 0,";
		      cout<<contig_size<<" ]. Skipping."<<endl;		
		    }
		    else{
		      int start_coord = cur_hit_5p_coord;
		      
		      if(cur_hit_orient == "+"){
			pos_hits[i].push_back(start_coord);
		      }
		      else{
			neg_hits[i].push_back(start_coord);
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
    
    align_str.close();  
  }
  cout<<endl<<endl;

  cout<<"sorting hits"<<endl;
  for(unsigned int i=0; i<contigs.size(); i++){
    sort(pos_hits[i].begin(), pos_hits[i].end());
    sort(neg_hits[i].begin(), neg_hits[i].end());
  }
  
  //cout<<"saving the binary align files..."<<endl;

  int collapsed_reads = 0;
  
  for(unsigned int i=0; i<contigs.size(); i++){

    ofstream cur_out_str;
    if(output_path != "not_applicable"){
      string cur_output_fname = output_path + "/" + contigs[i] + ".align.bin";
      remove(cur_output_fname.c_str());
      cur_out_str.open(cur_output_fname.c_str());
      if(!cur_out_str.good()){
	cerr<<"bad file name: "<<cur_output_fname<<endl;
	exit(1);
      }
    }


    string cur_contig = contigs[i];
    int cur_entries_pos = (int) pos_hits[i].size();
    int cur_entries_neg = (int) neg_hits[i].size();
    
    
    cout<<endl;
    cout<<"-------------------------"<<endl;
    cout<<"contig: "<<cur_contig<<endl;
    cout<<endl;
    cout<<"+ reads: "<<cur_entries_pos<<endl;
    cout<<"- reads: "<<cur_entries_neg<<endl;

    
    
    int cur_contig_name_size = cur_contig.size();
    //  cout<<endl<<"contig_name_size: "<<contig_name_size<<endl;
    
    if(collapse_reads == "false"){
      
      if(output_path != "not_applicable"){
	cur_out_str.write((char*) (&cur_contig_name_size), sizeof(cur_contig_name_size));
	cur_out_str.write((char*) &(cur_contig[0]), cur_contig_name_size*sizeof(char));
	
	cur_out_str.write((char*)(&cur_entries_pos), sizeof(cur_entries_pos));
	cur_out_str.write((char*)(&cur_entries_neg), sizeof(cur_entries_neg));
	
	cur_out_str.write((char*)(&pos_hits[i][0]), sizeof(int)*cur_entries_pos);
	cur_out_str.write((char*)(&neg_hits[i][0]), sizeof(int)*cur_entries_neg);
      }
      collapsed_reads += cur_entries_pos + cur_entries_neg;
    } 

    else{
      vector<short> pos_hit_counts(contig_sizes[i]);
      vector<short> neg_hit_counts(contig_sizes[i]);

      vector<int> pos_hits_collapsed;
      vector<int> neg_hits_collapsed;

      int stacks_collapsed = 0;
      int reads_in_collapsed_stacks = 0;

      //cout<<"Mapping counts"<<endl;
      cout<<endl;
      
      for(int j=0; j<contig_sizes[i]; j++){
	pos_hit_counts[j] = 0;
	neg_hit_counts[j] = 0;
      }

      for(unsigned int j=0; j<pos_hits[i].size(); j++){
	pos_hit_counts[pos_hits[i][j]] ++;
      }

      for(unsigned int j=0; j<neg_hits[i].size(); j++){
	neg_hit_counts[neg_hits[i][j]] ++;
      }

      for(int j=0; j<(int)pos_hit_counts.size(); j++){
	if(pos_hit_counts[j] >= count_threshold){
	  int reads_in_the_collapse_window =0;
	  int positions_hit = 0;
	  for(int k=int_max(0,j-collapse_window/2); k<=int_min(j+collapse_window/2,contig_sizes[i]-1); k++){
	    if(pos_hit_counts[k] > 0){
	      reads_in_the_collapse_window += pos_hit_counts[k];
	      positions_hit++;
	    }    
	  }
	  	  
	  double percent_positions_hit = 100 * ((double)positions_hit)/((double)(collapse_window+1));
	  int cur_stack = (int) pos_hit_counts[j];

	  //double p_v_precision = stack_p_value_threshold/10;
	  double cur_stack_p_value;
	  if(cur_stack <=3 && cur_stack>=2){
	    cur_stack_p_value = small_binomial_p_value(cur_stack, 1.0/((double)collapse_window), reads_in_the_collapse_window);
	  }
	  else{
	    cur_stack_p_value = 
	      Poisson_p_value(cur_stack, ((double)reads_in_the_collapse_window) / ((double)collapse_window),p_v_precision);
	  }

	  if(cur_stack_p_value <= stack_p_value_threshold && percent_positions_hit <= percent_positions_hit_threshold){
	    
	    if(QuEST_collapsed_fname != "not_applicable"){
	      QuEST_collapsed_ofstr<<contigs[i]<<" "<<j<<" +"<<endl;
	    }
	    for(int s=0; s<new_stack_size; s++) pos_hits_collapsed.push_back(j);
	    stacks_collapsed++;
	    reads_in_collapsed_stacks += cur_stack;
	  }
	  else{
	    for(int k=0; k<cur_stack; k++){
	      pos_hits_collapsed.push_back(j);
	      QuEST_collapsed_ofstr<<contigs[i]<<" "<<j<<" +"<<endl;
	    }
	  }
	}
	else{
	  if(pos_hit_counts[j] > 0){
	    for(int k=0; k<pos_hit_counts[j]; k++){
	      pos_hits_collapsed.push_back(j);
	      if(QuEST_collapsed_fname != "not_applicable"){
		QuEST_collapsed_ofstr<<contigs[i]<<" "<<j<<" +"<<endl;
	      }
	    }
	  }
	}
      }
      
      for(int j=0; j<(int)neg_hit_counts.size(); j++){
	if(neg_hit_counts[j] >= count_threshold){
	  int reads_in_the_collapse_window =0;
	  int positions_hit = 0;
	  for(int k=int_max(0,j-collapse_window/2); k<=int_min(j+collapse_window/2,contig_sizes[i]-1); k++){
	    if(neg_hit_counts[k] > 0){
	      reads_in_the_collapse_window += neg_hit_counts[k];
	      positions_hit++;
	    }    
	  }
	  	  
	  double percent_positions_hit = 100*((double)positions_hit)/((double)(collapse_window+1));
	  int cur_stack = (int) neg_hit_counts[j];

	  double cur_stack_p_value;
	  if(cur_stack <=3 && cur_stack>=2){
	    cur_stack_p_value = small_binomial_p_value(cur_stack, 1.0/((double)collapse_window), reads_in_the_collapse_window);
	  }
	  else{
	    cur_stack_p_value = 
	      Poisson_p_value(cur_stack, ((double)reads_in_the_collapse_window) / ((double)collapse_window),p_v_precision);
	  }
	  //   Poisson_p_value(cur_stack, ((double)reads_in_the_collapse_window) / ((double)collapse_window), p_v_precision);
	  if(cur_stack_p_value <= stack_p_value_threshold && percent_positions_hit <= percent_positions_hit_threshold){
	    if(QuEST_collapsed_fname != "not_applicable"){
	      QuEST_collapsed_ofstr<<contigs[i]<<" "<<j<<" -"<<endl;
	    }
	    for(int s=0; s<new_stack_size; s++) neg_hits_collapsed.push_back(j);
	    stacks_collapsed++;
	    reads_in_collapsed_stacks += cur_stack;
	  }
	  else{
	    for(int k=0; k<cur_stack; k++){
	      neg_hits_collapsed.push_back(j);
	      if(QuEST_collapsed_fname != "not_applicable"){
		QuEST_collapsed_ofstr<<contigs[i]<<" "<<j<<" -"<<endl;
	      }
	    }
	  }
	}
	else{
	  if(neg_hit_counts[j] > 0){
	    for(int k=0; k<neg_hit_counts[j]; k++){
	      neg_hits_collapsed.push_back(j);
	      if(QuEST_collapsed_fname != "not_applicable"){
		QuEST_collapsed_ofstr<<contigs[i]<<" "<<j<<" -"<<endl;
	      }
	    }
	  }
	}	
      }

      //exit(0);
      
      /*
      vector<int> pos_hits_collapsed;
      vector<int> neg_hits_collapsed;
      
      
      int last_coord = -1;
      for(unsigned int j=0; j<pos_hits[i].size(); j++){
	if(j>0){
	  if(last_coord == pos_hits[i][j]){ //do nothing
	  }
	  if(last_coord < pos_hits[i][j]){
	    pos_hits_collapsed.push_back(last_coord);
	  }
	}
	last_coord = pos_hits[i][j];
      }


      last_coord = -1;
      for(unsigned int j=0; j<neg_hits[i].size(); j++){
	if(j>0){
	  if(last_coord == neg_hits[i][j]){ //do nothing
	  }
	  if(last_coord < neg_hits[i][j]){
	    neg_hits_collapsed.push_back(last_coord);
	  }
	}
	last_coord = neg_hits[i][j];
      }
      */


      int cur_entries_pos_collapsed = pos_hits_collapsed.size();
      int cur_entries_neg_collapsed = neg_hits_collapsed.size();

      collapsed_reads += (cur_entries_pos_collapsed + cur_entries_neg_collapsed);
      
      cout<<"+ reads after collapsing: "<<cur_entries_pos_collapsed<<endl;
      cout<<"- reads after collapsing: "<<cur_entries_neg_collapsed<<endl;
      cout<<endl;
      
      cout<<"stacks collapsed: "<<stacks_collapsed<<endl;
      cout<<"reads in collapsed stacks: "<<reads_in_collapsed_stacks<<endl;

      
      if(output_path != "not_applicable"){
	cur_out_str.write((char*) (&cur_contig_name_size), sizeof(cur_contig_name_size));
	cur_out_str.write((char*) &(cur_contig[0]), cur_contig_name_size*sizeof(char));
	
	cur_out_str.write((char*)(&cur_entries_pos_collapsed), sizeof(cur_entries_pos_collapsed));
	cur_out_str.write((char*)(&cur_entries_neg_collapsed), sizeof(cur_entries_neg_collapsed));
	
	cur_out_str.write((char*)(&pos_hits_collapsed[0]), sizeof(int)*cur_entries_pos_collapsed);
	cur_out_str.write((char*)(&neg_hits_collapsed[0]), sizeof(int)*cur_entries_neg_collapsed);
      }

    }

    if(output_path != "not_applicable"){
      cur_out_str.close();
    }
  }  

  if(report_fname != "not_applicable"){
    remove(report_fname.c_str());
    ofstream report_ofstr(report_fname.c_str());
    
    if(!report_ofstr.good()){
      cerr<<"Bad file name "<<report_fname<<endl;
      exit(1);
    }
    
    report_ofstr<<collapsed_reads<<endl;
    report_ofstr.close();
  }

  QuEST_collapsed_ofstr.close();
  
  cout<<"done!"<<endl;
  cout<<endl;

  return 0;

}
