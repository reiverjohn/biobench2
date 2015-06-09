#include <new>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "params.h"

//using namespace std;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::setprecision;


using std::ios;


void new_handler(){
  cout<<"The new_handler is called from peak_caller.cpp"<<endl;
  cout<<"This is usually a sign of unsuccessful memory allocation"<<endl;
  cout<<"Exiting!\n";
  
  exit(1);
}


int main(int argc, char* argv[]){
  
  params pars(argc, argv);
  pars.require("score_file1","score_file1",STRING_TYPE);
  pars.require("score_file2","score_file2",STRING_TYPE);
  pars.require("output_file","output_file",STRING_TYPE);
  pars.require("ChIP_threshold","ChIP_threshold",DOUBLE_TYPE);
  pars.require("ChIP_extension_threshold", "ChIP_extension_threshold", DOUBLE_TYPE);
  pars.require("ChIP_to_background_ratio", "ChIP_to_background_ratio", DOUBLE_TYPE);
  pars.optional("peak_caller_method","maximum/local_maximum","local_maximum",STRING_TYPE);
  pars.optional("local_maximum_radius","local_maximum_radius","10",INT_TYPE);
  pars.require("contig_id","contig_id",STRING_TYPE);
  //pars.require("background_threshold","background_threshold",DOUBLE_TYPE);
  pars.optional("dip_fraction","dip_fraction","0.1",DOUBLE_TYPE);
  //pars.require("rescue_ratio","rescue_ratio",DOUBLE_TYPE);

  if(!pars.enforce()) exit(1);  

  string scores_fname1 = pars.get_string_value("score_file1");
  string scores_fname2 = pars.get_string_value("score_file2");
  string output_fname = pars.get_string_value("output_file");  

  string contig_id = pars.get_string_value("contig_id");

  //double score_threshold = pars.get_double_value("chip_threshold");
  double ChIP_threshold = pars.get_double_value("ChIP_threshold");
  double ChIP_extension_threshold = pars.get_double_value("ChIP_extension_threshold");
  //double background_threshold = pars.get_double_value("background_threshold");
  double ChIP_to_background_ratio = pars.get_double_value("ChIP_to_background_ratio");
  double dip_fraction = pars.get_double_value("dip_fraction");
  //double above_threshold_mult = pars.get_double_value("rescue_ratio");

  //string peak_caller_method = pars.get_string_value("peak_caller_method");
  int local_maximum_radius = pars.get_int_value("local_maximum_radius");

  bool background_used;
  if(scores_fname2 == "NA"){
    background_used = false;
  }
  else{
    background_used = true;
  }


  /*
  if(peak_caller_method == "maximum" || peak_caller_method == "local_maximum"){}
  else{
    cerr<<"bad peak_caller_method: "<<peak_caller_method<<endl;
    exit(1);
  }
  */

  cout<<endl<<"contig: "<<contig_id<<endl;
  cout<<"score_file1: "<<scores_fname1<<endl;
  cout<<"score_file2: "<<scores_fname2<<endl;
  cout<<"ChIP_threshold: "<<ChIP_threshold<<endl;
  cout<<"ChIP_extension_threshold: "<<ChIP_extension_threshold<<endl;
  cout<<"ChIP_to_background_ratio: "<<ChIP_to_background_ratio<<endl;
  //  cout<<"background_threshold: "<<background_threshold<<endl;
  cout<<endl;

  ifstream scores_str1(scores_fname1.c_str());
  if(!scores_str1.good()){
    cerr<<"bad file: "<<scores_fname1<<endl;
    exit(1);
  }
  long begin1, end1, file_size1;

  begin1 = scores_str1.tellg();
  scores_str1.seekg(0, ios::end);
  end1 = scores_str1.tellg();
  scores_str1.seekg(0, ios::beg);
  file_size1 = end1 - begin1;

  long begin2, end2, file_size2 = 0; //for checking file sizes
  ifstream scores_str2;
  if(background_used){
    scores_str2.open(scores_fname2.c_str());
    if(!scores_str2.good()){
      cerr<<"bad file: "<<scores_fname2<<endl;
      exit(1);
    }
    
    begin2 = scores_str2.tellg();
    scores_str2.seekg(0, ios::end);
    end2 = scores_str2.tellg();
    scores_str2.seekg(0, ios::beg);
    file_size2 = end2 - begin2;
    
    if(file_size1 != file_size2){
      cout<<"Error:\n"<<endl;
      cout<<"File sizes for score ChIP and backgound files do not match."<<endl;
      cout<<"This usually happens when not enough disk space is available."<<endl;
      cout<<endl;
      cout<<"File1: "<<scores_fname1<<endl;
      cout<<"Size: "<<file_size1<<endl;
      cout<<endl;
      cout<<"File2: "<<scores_fname2<<endl;
      cout<<"Size: "<<file_size2<<endl;
      cout<<endl;
      cout<<"Exiting."<<endl;
      exit(1);
    }
  }

  cout<<"loading the scores"<<endl;

  unsigned int _size1, _size2;
  int entries1, entries2;

  scores_str1.read((char*) (&_size1), sizeof(_size1));
  if(background_used){
    scores_str2.read((char*) (&_size2), sizeof(_size2));
  }

  cout<<endl;
  cout<<"_size1: "<<file_size1<<endl;
  if(background_used){
    cout<<"_size2: "<<file_size2<<endl;
  }
  cout<<endl;


  cout<<"File1 entries: "<<_size1<<endl;
  if(background_used){
    cout<<"File2 entries: "<<_size2<<endl;
  }
   cout<<endl;


  //  assert(_size1 == _size2);
  
  if(background_used){
    if(_size1 != _size2){
      cout<<"Error: \n"<<endl;
      cout<<"Number of entries do not match between the ChIP and background score files."<<endl;
      cout<<endl;
      cout<<"File1: "<<scores_fname1<<endl;
      cout<<"_size1: "<<_size1<<endl;
      cout<<endl;
      cout<<"File2: "<<scores_fname2<<endl;
      cout<<"_size2: "<<_size2<<endl;
      cout<<endl;
      cout<<"Exiting."<<endl;
      exit(1);
    }
  }

  scores_str1.read((char*) (&entries1), sizeof(entries1));
  if(background_used){
    scores_str2.read((char*) (&entries2), sizeof(entries2));
  }

  long expected_file1_size = 2*sizeof(unsigned int) + _size1*sizeof(float);
  if(file_size1 != expected_file1_size){
    cout<<"Bad file size for the file: "<<scores_fname1<<endl;
    cout<<"Expected:    "<<expected_file1_size<<endl;
    cout<<"Actual size: "<<file_size1<<endl;
    cout<<"Exiting"<<endl;
    exit(1);
  }

  if(background_used){
    long expected_file2_size = 2*sizeof(unsigned int) + _size2*sizeof(float);
    if(file_size2 != expected_file2_size){
      cout<<"Bad file size for the file: "<<scores_fname2<<endl;
      cout<<"Expected:    "<<expected_file2_size<<endl;
      cout<<"Actual size: "<<file_size2<<endl;
      cout<<"Exiting"<<endl;
      exit(1);
    }
  }

  // cout<<"the scores file contains "<<_size1<<" entries"<<endl;

  cout<<"score file 1: "<<entries1<<" tags."<<endl;
  if(background_used){
    cout<<"score file 2: "<<entries2<<" tags."<<endl;
  }

  cout<<endl;

  cout<<"allocating memory for the ChIP scores."<<endl;
  //float* scores_ChIP = new float[_size1];
  vector<float> scores_ChIP(_size1);
  //std::set_new_handler(&new_handler); //make sure allocation was ok

  cout<<"allocating memory for the background scores."<<endl;
  vector<float> scores_background(_size2);
  //float* scores_background = new float[_size2];
  //std::set_new_handler(&new_handler);

  cout<<"allocating memory for the mask."<<endl;
  vector<bool> mask(_size1);
  //bool* mask = new bool[_size1];
  //std::set_new_handler(&new_handler);

  cout<<endl;
  cout<<"loading the ChIP score values..."<<endl;
  scores_str1.read((char*) &(scores_ChIP[0]), sizeof(scores_ChIP[0])*_size1);
  cout<<"file1 loaded."<<endl;
  scores_str2.read((char*) &(scores_background[0]), sizeof(scores_background[0])*_size2);
  cout<<"file2 loaded."<<endl;
  cout<<endl;

  scores_str1.close();
  scores_str2.close();

  unsigned int _size = _size1;

  remove(output_fname.c_str());
  ofstream out_str(output_fname.c_str());
  if(!out_str.good()){
    cerr<<"bad file name: "<<output_fname<<endl;
    exit(1);
  }
  
  cout<<"detecting seeds..."<<endl;
  //detect initial seeds
  for(unsigned int i=0; i<_size; i++){
    if(scores_ChIP[i] >= ChIP_threshold){
      mask[i] = true;
    }
    else{
      mask[i] = false;
    }
  }

  cout<<"extending seeds..."<<endl;
  //extend the seeds
  for(unsigned int i=1; i<_size-1; i++){
    if(mask[i-1] == false && mask[i] == true){
      //can be extended to the left
      bool extend = true;
      int j = i-1;
      while(extend){
	if(j==0){
	  extend = false;
	}
	if(scores_ChIP[j] >= ChIP_extension_threshold){
	  mask[j] = true;
	}
	else{
	  extend = false;
	}
	j--;
      }
    }
    if(mask[i] == true && mask[i+1] == false){
      //can be extended to the right 
      bool extend = true;
      unsigned int j = i+1;
      while(extend){
	if(j>=_size-1){
	  extend = false;
	}
	if(scores_ChIP[j] >= ChIP_extension_threshold){
	  mask[j] = true;
	}
	else{
	  extend = false;
	}
	j++;
      }
    }
  }
  

  cout<<"finding regions..."<<endl;

  vector<unsigned int> region_starts;
  vector<unsigned int> region_ends;

  bool prev_mask = false;
  
  for(unsigned int i=0; i<_size-1; i++){
    if(mask[i] == true && prev_mask == false){      

      //left boundary of a region      
      //find right boundary of a region
      bool right_boundary_found = false;
      int region_end = -1;
      unsigned int j=i+1;

      while(!right_boundary_found){
	if(j>=_size-1){
	  right_boundary_found = true;
	  region_end = _size;
	}
	else{
	  if(mask[j] == false){
	    right_boundary_found = true;
	    region_end = j;
	  }
	}
	j++;
      }
      if(region_end > 0 && region_end < (int)_size){
	region_starts.push_back(i);
	region_ends.push_back(region_end);
      }
    }

    prev_mask = mask[i];
  }

  cout<<"removing unsatisfactory peaks..."<<endl;

  //find peaks
  int peaks_passed = 0;
  for(unsigned int i=0; i<region_starts.size(); i++){
    unsigned int region_begin = region_starts.at(i);
    unsigned int region_end = region_ends.at(i);

    vector<int> peak_candidates; //coordinates of peak candidates
    for(unsigned int j=region_begin; j<region_end; j++){
      
      bool local_maximum = true;

      local_maximum_radius = 1;
      for(unsigned int k=j-(unsigned int)local_maximum_radius; k<=j+(unsigned int)local_maximum_radius; k++){

	if(k>0 && k<_size){// && k!=j){
	  if(scores_ChIP.at(k) > scores_ChIP.at(j)){
	    local_maximum = false;
	  }
	}
      }


      if(local_maximum){

	if(scores_ChIP.at(j) >= ChIP_threshold){
	  if(background_used){
	    double cur_ChIP_to_background_ratio = -1;
	    if(scores_background.at(j) > 0){
	      cur_ChIP_to_background_ratio = scores_ChIP.at(j) / scores_background.at(j);
	    }
	    if(cur_ChIP_to_background_ratio == -1 || 
	       cur_ChIP_to_background_ratio >= ChIP_to_background_ratio){
	      peak_candidates.push_back(j);
	    }
	  }
	  else{
	    peak_candidates.push_back(j);
	  }
	}
      }
    }

    if(peak_candidates.size() > 0){
      vector<int> peak_candidates_mask;
      for(unsigned int j=0; j<peak_candidates.size(); j++){
	peak_candidates_mask.push_back(0);
      }
      
      //first find maximum
      
      int max_ind = -1;
      float max_score = -1;

      for(unsigned int j=0; j<peak_candidates.size(); j++){
	if(scores_ChIP.at(peak_candidates.at(j)) > max_score){
	  max_score = scores_ChIP.at(peak_candidates.at(j)); 
	  max_ind = j;
	}
      }

      assert(max_ind >= 0);
      peak_candidates_mask.at(max_ind) = 1; //pass the top peak
      float region_max_score = max_score;
      int region_max_coord = peak_candidates.at(max_ind);

      bool continue_with_candidates = true;
      
      
      //keep iterating through the list in the decreasing order of peaks
      while(continue_with_candidates){
	continue_with_candidates = false;
	int new_max_ind = -1; 
	float new_max_score = -1;
	
	for(unsigned int j=0; j<peak_candidates.size(); j++){
	  if(peak_candidates_mask.at(j) == 0){ //not been analyzed
	    if(scores_ChIP.at(peak_candidates.at(j)) > new_max_score){
	      new_max_ind = j;
	      new_max_score = scores_ChIP.at(peak_candidates.at(j));
	    }
	  }
	}
	if(new_max_ind >= 0){
	  continue_with_candidates = true;

	  int left_tmp_coord, right_tmp_coord;
	  double min_score = (1-dip_fraction)*scores_ChIP.at(peak_candidates.at(new_max_ind));
	  if(new_max_ind > max_ind){
	    left_tmp_coord = peak_candidates.at(max_ind);
	    right_tmp_coord = peak_candidates.at(new_max_ind);

	  }
	  else{
	    left_tmp_coord = peak_candidates.at(new_max_ind);
	    right_tmp_coord = peak_candidates.at(max_ind);
	  }

	  bool min_threshold_satisfied = false;
	  for(int k=left_tmp_coord; k<=right_tmp_coord; k++){
	    if(scores_ChIP.at(k) <= min_score){
	      min_threshold_satisfied = true;
	    }
	  }
	  if(min_threshold_satisfied){
	    peak_candidates_mask.at(new_max_ind) = 1;
	    max_ind = new_max_ind;
	  }
	  else{
	    peak_candidates_mask.at(new_max_ind) = -1;
	  }
	}
      }

      

      //cout<<"outputting the calls..."<<endl;

      out_str<<"R ";
      out_str<<contig_id<<" "<<region_begin<<"-"<<region_end-1;
      out_str<<" ChIP: "<<setprecision(3)<<region_max_score;
      if(background_used){
	out_str<<" control: "<<setprecision(3)<<scores_background.at(region_max_coord);
      }
      else{
	out_str<<" control: 0";
      }
      out_str<<" max_pos: "<<region_max_coord;
      out_str<<endl;

      for(unsigned int j=0; j<peak_candidates.size(); j++){
	if(peak_candidates_mask.at(j) == 1){
	  peaks_passed++;
	  out_str<<"P ";
	  out_str<<contig_id<<" "<<peak_candidates.at(j);
	  out_str<<" ChIP: "<<setprecision(3)<<scores_ChIP.at(peak_candidates.at(j));
	  if(background_used){
	    out_str<<" control: "<<setprecision(3)<<scores_background.at(peak_candidates.at(j));
	  }
	  else{
	    out_str<<" control: 0";
	  }
	  out_str<<" region: "<<region_begin<<"-"<<region_end;
	  out_str<<endl;
	}
      }
      out_str<<endl;
    }
    
  }

  cout<<endl;
  cout<<"found: "<<region_starts.size()<<" regions and "<<peaks_passed<<" peaks"<<endl;

  out_str.close();

  return 0;
}
