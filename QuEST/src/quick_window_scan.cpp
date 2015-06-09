#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <list>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "string_utils.h"
#include "seq_contig.h"
#include "params.h"

using namespace std;

#define BLACK_BACKGROUND

int el(const int* counts, unsigned int i, unsigned int j){
  return counts[i*3+j];
}

double weight(int x){
  
  if(-50<= x && x<=50) return x/15000.0+1.0/300.0;
  if(50<x && x<=250) return -1*x/(30000) + 1.0/120.0;
  return 0;
  
}

double weight_for(int x){
  return weight(-x);
}
double weight_rev(int x){
  return weight(x);
}

double normal_weight(int x, double sigma){
  double pi=3.14;
  return exp(-x*x/(2*sigma*sigma))/(sqrt(2*pi)*sigma);
}

int main(int argc, char* argv[]){

  params pars(argc, argv);
  pars.require("ChIP_bin_align_file","align_file",STRING_TYPE);
  pars.require("RX_noIP_bin_align_file","align_file",STRING_TYPE);
  pars.require("contig_id","contig_id",STRING_TYPE);
  pars.require("contig_size","contig_size",INT_TYPE);
  pars.require("output_file","output_summary_file",STRING_TYPE);
 
  pars.optional("calc_window","calc_window","300",INT_TYPE);
  pars.optional("ChIP_tags_threshold", "ChIP_tags_threshold", "600", INT_TYPE);
  pars.optional("ChIP_to_background_tag_ratio", "ChIP_to_background_tag_ratio", "NA", DOUBLE_TYPE);
  //pars.optional("background_tags_threshold", "background_tags_threshold", "NA", INT_TYPE);
  //  pars.optional("count_threshold","count_threshold","600",INT_TYPE);
  //  pars.require("enrichment_ratio","enrichment_ratio",DOUBLE_TYPE);

  if(!pars.enforce()){
    pars.list_all_params();
    exit(1);
  }

  pars.list_all_params();
  string ChIP_bin_align_fname = pars.get_string_value("ChIP_bin_align_file");
  string RX_noIP_bin_align_fname = pars.get_string_value("RX_noIP_bin_align_file");

  string contig_id = pars.get_string_value("contig_id");
  int contig_size = pars.get_int_value("contig_size");
  string output_fname = pars.get_string_value("output_file");

  int calc_window = pars.get_int_value("calc_window");
  int ChIP_tags_threshold = pars.get_int_value("ChIP_tags_threshold");

  double ChIP_to_background_tag_ratio;
  string ChIP_to_background_tag_ratio_string = pars.get_string_value("ChIP_to_background_tag_ratio");
  if(ChIP_to_background_tag_ratio_string == "NA"){
    ChIP_to_background_tag_ratio = pars.get_double_value("ChIP_to_background_tag_ratio");
  }
  else{
    ChIP_to_background_tag_ratio = -1;
  }

  //int count_threshold = pars.get_int_value("count_threshold");
  //double enrichment_ratio = pars.get_double_value("enrichment_ratio");

  cout<<endl;
  cout<<"Quick Window Scan will use the following parameters:"<<endl;
  cout<<endl;
  cout<<"- - - - - - - - - - - - "<<endl;
  
  cout<<endl;

  cout<<"ChIP_bin_align_file:            "<<ChIP_bin_align_fname<<endl;
  cout<<"RX_noIP_bin_align_file:         "<<RX_noIP_bin_align_fname<<endl;
  cout<<"output_file:                "<<output_fname<<endl;
  cout<<endl;
  

  cout<<"contig_id:                 "<<contig_id<<endl;
  cout<<"contig_size:               "<<contig_size<<endl;
  cout<<endl;

  cout<<"calc_window:               "<<calc_window<<endl;
  cout<<"ChIP_tags_threshold:       "<<ChIP_tags_threshold<<endl;
  
  //cout<<"enrichment_ratio:          "<<enrichment_ratio<<endl;
  cout<<endl;
  
  cout<<"- - - - - - - - - - - - "<<endl;

  if(contig_size <= 0){
    cout<<"Bad size of the contig: "<<contig_size<<endl;
    cout<<"Aborting."<<endl;
    exit(0);
  }

  ifstream ChIP_bin_align_str(ChIP_bin_align_fname.c_str());
  if(!ChIP_bin_align_str.good()){
    cerr<<"bad file name: "<<ChIP_bin_align_fname<<endl;
    exit(1);
  }
  

  // reading ChIP data
  int ChIP_entries_pos;
  int ChIP_entries_neg;

  int ChIP_contig_name_size;

  ChIP_bin_align_str.read((char*)(&ChIP_contig_name_size), sizeof(ChIP_contig_name_size));

  string ChIP_cur_contig(ChIP_contig_name_size, ' ');
  ChIP_bin_align_str.read((char*)(&ChIP_cur_contig[0]), ChIP_contig_name_size*sizeof(char));

  ChIP_bin_align_str.read((char*)(&ChIP_entries_pos), sizeof(ChIP_entries_pos));
  ChIP_bin_align_str.read((char*)(&ChIP_entries_neg), sizeof(ChIP_entries_neg));

  vector<int> ChIP_pos_hits(ChIP_entries_pos);
  vector<int> ChIP_neg_hits(ChIP_entries_neg);

  ChIP_bin_align_str.read((char*)(&ChIP_pos_hits[0]), sizeof(int)*ChIP_entries_pos);
  ChIP_bin_align_str.read((char*)(&ChIP_neg_hits[0]), sizeof(int)*ChIP_entries_neg);

  ChIP_bin_align_str.close();

  cout<<"ChIP:    "<<ChIP_cur_contig<<" +: "<<ChIP_entries_pos<<" -: "<<ChIP_entries_neg<<endl;

  //  /reading ChIP data

  // reading RX_noIP

  
  bool background_data_used = false;

  if(RX_noIP_bin_align_fname == "NA"){
    background_data_used = false;
  }
  else{
    background_data_used = true;
  }

  ifstream RX_noIP_bin_align_str;
  if(background_data_used){
    RX_noIP_bin_align_str.open(RX_noIP_bin_align_fname.c_str());
    if(!RX_noIP_bin_align_str.good()){
      cerr<<"bad file name: "<<RX_noIP_bin_align_fname<<endl;
      exit(1);
    }
  }
  
  int RX_noIP_entries_pos = 0;
  int RX_noIP_entries_neg = 0;
  
  string RX_noIP_cur_contig = "";
  
  if(background_data_used){
    //cout<<"reading control data..."<<endl;
    int RX_noIP_contig_name_size;
    
    RX_noIP_bin_align_str.read((char*)(&RX_noIP_contig_name_size), sizeof(RX_noIP_contig_name_size));
    //cout<<"Contig name size: "<<RX_noIP_contig_name_size<<endl;
    
    string tmp_RX_noIP_cur_contig(RX_noIP_contig_name_size, ' ');
    
    RX_noIP_bin_align_str.read((char*)(&tmp_RX_noIP_cur_contig[0]), RX_noIP_contig_name_size*sizeof(char));
    RX_noIP_cur_contig = tmp_RX_noIP_cur_contig;
    //cout<<"Read contig name: "<<RX_noIP_cur_contig<<endl;
    
    RX_noIP_bin_align_str.read((char*)(&RX_noIP_entries_pos), sizeof(RX_noIP_entries_pos));
    RX_noIP_bin_align_str.read((char*)(&RX_noIP_entries_neg), sizeof(RX_noIP_entries_neg));

    RX_noIP_cur_contig = tmp_RX_noIP_cur_contig;
    //cout<<"loaded RX_noIP binary entries"<<endl;
  }    
  vector<int> RX_noIP_pos_hits(RX_noIP_entries_pos);
  vector<int> RX_noIP_neg_hits(RX_noIP_entries_neg);
  
  if(background_data_used){
    RX_noIP_bin_align_str.read((char*)(&RX_noIP_pos_hits[0]), sizeof(int)*RX_noIP_entries_pos);
    RX_noIP_bin_align_str.read((char*)(&RX_noIP_neg_hits[0]), sizeof(int)*RX_noIP_entries_neg);
    
    RX_noIP_bin_align_str.close();
    
    
    cout<<"RX_noIP: "<<RX_noIP_cur_contig<<" +: "<<RX_noIP_entries_pos<<" -: "<<RX_noIP_entries_neg<<endl;
    //
    //assert(ChIP_cur_contig == RX_noIP_cur_contig);
  }
  
    // /reading RX_noIP data
  

  remove(output_fname.c_str());
  ofstream out_str(output_fname.c_str());
  if(!out_str.good()){
    cerr<<"bad file name: "<<output_fname<<endl;
    exit(1);
  }

  vector<int> ChIP_read_counts(contig_size);
  vector<int> RX_noIP_read_counts(contig_size);

  for(int i=0; i<contig_size; i++){
    ChIP_read_counts[i] = 0;
    RX_noIP_read_counts[i] = 0;
  }


  for(unsigned int i=0; i<ChIP_pos_hits.size(); i++){
    int cur_read_5p_coord = ChIP_pos_hits[i];
    assert(cur_read_5p_coord < contig_size);
    ChIP_read_counts[cur_read_5p_coord]++;
  }

  for(unsigned int i=0; i<ChIP_neg_hits.size(); i++){
    int cur_read_5p_coord = ChIP_neg_hits[i];
    assert(cur_read_5p_coord < contig_size);
    ChIP_read_counts[cur_read_5p_coord]++;
  }

  for(unsigned int i=0; i<RX_noIP_pos_hits.size(); i++){
    int cur_read_5p_coord = RX_noIP_pos_hits[i];
    assert(cur_read_5p_coord < contig_size);
    RX_noIP_read_counts[cur_read_5p_coord]++;
  }

  for(unsigned int i=0; i<RX_noIP_neg_hits.size(); i++){
    int cur_read_5p_coord = RX_noIP_neg_hits[i];
    assert(cur_read_5p_coord < contig_size);
    RX_noIP_read_counts[cur_read_5p_coord]++;
  }
  

  cout<<"outputting positive regions"<<endl;

  int cur_ChIP_read_count = 0;
  int cur_RX_noIP_read_count = 0;
  int cur_status = -1; //uninit
  int cur_max = 0;
  double cur_max_ratio = -1;
  unsigned int cur_max_ind = 0;

  
  
  ///
  ///
  double cur_ChIP_to_background_tag_ratio = -1;

  for(int i=0; i<calc_window; i++){
    cur_ChIP_read_count += ChIP_read_counts[i];
    cur_RX_noIP_read_count += RX_noIP_read_counts[i];

    if(cur_RX_noIP_read_count > 0){
      cur_ChIP_to_background_tag_ratio = ((double) cur_ChIP_read_count ) / ((double) cur_RX_noIP_read_count);
    }
    else{ cur_ChIP_to_background_tag_ratio = -1; }
  }

  if(cur_ChIP_read_count >= ChIP_tags_threshold){

    //if(cur_ChIP_to_background_tag_ratio >= ChIP_to_background_tag_ratio){// || cur_ef  == -1)){ //in the peak
    if(cur_ChIP_read_count > cur_max){
      if(background_data_used){
	if(cur_ChIP_to_background_tag_ratio > ChIP_to_background_tag_ratio || cur_ChIP_to_background_tag_ratio == -1){
	  cur_status = 1;
	  cur_max = cur_ChIP_read_count;
	  cur_max_ind = calc_window/2 - 1;
	}
	else{
	  cur_status = 0;
	}
      }
      else{
	cur_status = 1;
	cur_max = cur_ChIP_read_count;
	//cur_max_ratio = cur_ChIP_to_background_tag_ratio;
	cur_max_ind = calc_window/2-1;    
      }
    }
    else{
      cur_status = 1;
      cur_max = cur_ChIP_read_count;
      cur_max_ratio = cur_ChIP_to_background_tag_ratio;
      cur_max_ind = calc_window/2-1;    
    }
  }
  else{ cur_status = 0; } //waiting to get over the threshold
  
  int regions_found = 0;
  
  for(int i=calc_window/2; i<(int) contig_size-(int)calc_window; i++){
    
    cur_ChIP_read_count += ChIP_read_counts[i+calc_window/2] - ChIP_read_counts[i-calc_window/2];
    
    cur_RX_noIP_read_count += RX_noIP_read_counts[i+calc_window/2] - RX_noIP_read_counts[i-calc_window/2];
    
    if(cur_RX_noIP_read_count > 0){      
      cur_ChIP_to_background_tag_ratio = ((double) cur_ChIP_read_count ) / ((double) cur_RX_noIP_read_count);
    }
    else{ cur_ChIP_to_background_tag_ratio = -1; }
    
    switch (cur_status){
    case 0: {
      if (cur_ChIP_read_count >= ChIP_tags_threshold){
	if(background_data_used){ 
	  //(cur_ef >= enrichment_ratio || cur_ef == -1)){
	  if(cur_ChIP_to_background_tag_ratio >= ChIP_to_background_tag_ratio || cur_ChIP_to_background_tag_ratio == -1){
	    cur_status = 1;
	    cur_max = cur_ChIP_read_count;
	    cur_max_ind = i;
	  }
	  else{
	    cur_status = 0;
	  }
	}
	else{
	  cur_status = 1; 
	  cur_max = cur_ChIP_read_count;
	  cur_max_ind = i;	
	  //cur_max_ef = cur_ef;
	}
      }
      break;
    }
    case 1:{
      if(cur_ChIP_read_count < ChIP_tags_threshold){
	//if(cur_max_ef >= 0){
	out_str<<contig_id<<" "<<cur_max_ind<<" "<<cur_max<<" "<<endl;//cur_max_ef<<endl;
	//}
	//else{
	// out_str<<contig_id<<" "<<cur_max_ind<<" "<<cur_max<<" inf"<<endl;
	  //}
	cur_status = 0;	
	cur_max = 0;
	//cur_ef = 0;
	regions_found++;
      }
      else{
	if(cur_ChIP_read_count > cur_max){
	  if(background_data_used){
	    if(cur_ChIP_to_background_tag_ratio >= ChIP_to_background_tag_ratio){
	      cur_max = cur_ChIP_read_count;
	      cur_max_ind = i;
	    }
	    else{
	      out_str<<contig_id<<" "<<cur_max_ind<<" "<<cur_max<<" "<<endl;//cur_max_ef<<endl;
	      
	      cur_status = 0;
	      cur_max = 0;
	      regions_found++;
	    }
	  }
	  else{
	    cur_max = cur_ChIP_read_count;
	    //cur_max_ef = cur_ef;
	    cur_max_ind = i;
	  }
	}
      }
      break;
    }
    }
  }

  cout<<endl<<"identified "<<regions_found<<" candidated regions"<<endl<<endl;
  
  out_str.close();
  
  return 0;

}
