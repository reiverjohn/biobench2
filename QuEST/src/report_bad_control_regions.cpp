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

struct region{
  int region_begin;
  int region_end;  
  double mean_fe;
  double max_fe;
};

int main(int argc, char* argv[]){
  
  params pars(argc, argv);
  pars.require("score_file","score_file",STRING_TYPE);
  pars.require("output_file","output_file",STRING_TYPE);

  pars.require("score_normalizer", "score_normalizer", DOUBLE_TYPE);
  pars.require("region_seeding_threshold", "region_seeding_threshold", DOUBLE_TYPE);
  pars.require("region_extension_threshold", "region_extension_threshold", DOUBLE_TYPE);

  pars.require("contig_id","contig_id",STRING_TYPE);

  if(!pars.enforce()) exit(1);  

  string scores_fname = pars.get_string_value("score_file");
  string output_fname = pars.get_string_value("output_file");  

  string contig_id = pars.get_string_value("contig_id");

  //double score_threshold = pars.get_double_value("chip_threshold");
  double score_normalizer = pars.get_double_value("score_normalizer");  
  double region_seeding_threshold = pars.get_double_value("region_seeding_threshold"); //applied after normalization
  double region_extension_threshold = pars.get_double_value("region_extension_threshold"); //applied after normalization

  cout<<endl<<"contig: "<<contig_id<<endl;
  cout<<"score_file: "<<scores_fname<<endl;
  cout<<"region_seeding_threshold: "<<region_seeding_threshold<<endl;
  cout<<"region_extension_threshold: "<<region_extension_threshold<<endl;
  cout<<endl;

  ifstream scores_str(scores_fname.c_str());
  if(!scores_str.good()){
    cerr<<"bad file: "<<scores_fname<<endl;
    exit(1);
  }

  long begin, end, file_size;

  begin = scores_str.tellg();
  scores_str.seekg(0, std::ios::end);
  end = scores_str.tellg();
  scores_str.seekg(0, std::ios::beg);
  file_size = end - begin;

  //cout<<"File size: "<<file_size<<endl;

  cout<<"loading control data score file"<<endl;

  unsigned int _size;
  int entries; 

  scores_str.read((char*) (&_size), sizeof(_size));

  cout<<"Score file is expected to contain "<<_size<<" entries."<<endl;
  cout<<endl;

  scores_str.read((char*) (&entries), sizeof(entries));
  long expected_file_size = 2*sizeof(unsigned int) + _size*sizeof(float);

  if(file_size != expected_file_size){
    cout<<"Bad file size for the file: "<<scores_fname<<endl;
    cout<<"Expected:    "<<expected_file_size<<endl;
    cout<<"Actual size: "<<file_size<<endl;
    cout<<"Exiting"<<endl;
    exit(1);
  }

  cout<<"score file: "<<entries<<" tags."<<endl;
  cout<<endl;

  cout<<"allocating memory for the scores."<<endl;
  vector<float> scores(_size);

  cout<<"allocating memory for the mask."<<endl;
  vector<bool> mask(_size);

  cout<<endl;
  cout<<"loading the score values..."<<endl;
  scores_str.read((char*) &(scores[0]), sizeof(scores[0])*_size);
  cout<<"Score file loaded."<<endl;
  cout<<endl;

  scores_str.close();

  remove(output_fname.c_str());
  ofstream out_str(output_fname.c_str());
  if(!out_str.good()){
    cerr<<"bad file name: "<<output_fname<<endl;
    exit(1);
  }
  
  cout<<"detecting seeds..."<<endl;
  //detect initial seeds
  for(unsigned int i=0; i<_size; i++){
    if(scores[i] >= region_seeding_threshold * score_normalizer){
      mask[i] = true;
    }
    else{
      mask[i] = false;
    }
  }

  cout<<"extending seeds..."<<endl;
  //  cout<<"extenging to the left..."<<endl;
  for(unsigned int i=1; i<_size; i++){
    if(mask[i-1] == true && scores[i] >= region_extension_threshold * score_normalizer){
      mask[i] = true; //extending to the right
    }
  }

  //cout<<"extending to the right..."<<endl;
  for(unsigned int i=_size-1; i>=1; i--){
    if(mask[i] == true && scores[i-1] >= region_extension_threshold * score_normalizer){
      mask[i-1] = true; //extending to the left
    }
  }
  
  cout<<"finding regions..."<<endl;

  vector<region> artifactual_regions;


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
	region new_region;
	new_region.region_begin = i;
	new_region.region_end = region_end;
	artifactual_regions.push_back(new_region);
      }
    }
    prev_mask = mask[i];
  }

  cout<<"collecting stats on artifactual regions..."<<endl;
  for(unsigned int i=0; i<artifactual_regions.size(); i++){
    double score_mean = 0;
    double score_max = 0;
    int region_size = artifactual_regions[i].region_end - artifactual_regions[i].region_begin;
    assert(region_size > 0);
    for(int j=artifactual_regions[i].region_begin; j<artifactual_regions[i].region_end; j++){
      score_mean += scores[j] / (double(region_size));
      if(scores[j] > score_max) score_max = scores[j];
    }
    artifactual_regions[i].mean_fe = score_mean / score_normalizer;
    artifactual_regions[i].max_fe = score_max / score_normalizer;
  }

  cout<<endl;
  cout<<"Outputting regions... "<<endl;
  for(unsigned int i=0; i<artifactual_regions.size(); i++){
    
    out_str<<"R ";
    out_str<<contig_id<<" "<<artifactual_regions[i].region_begin<<"-"<<artifactual_regions[i].region_end;
    out_str<<" mean_fe: "<<setprecision(3)<<artifactual_regions[i].mean_fe;
    out_str<<" max_fe: "<<setprecision(3)<<artifactual_regions[i].max_fe;
    out_str<<endl;
    
  }
  out_str.close();

  cout<<endl;
  cout<<"-------------------------------------"<<endl;
  cout<<"|  found: "<<artifactual_regions.size()<<" artifactual region(s)   |"<<endl;
  cout<<"-------------------------------------"<<endl;

  return 0;
}
