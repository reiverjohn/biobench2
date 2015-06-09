#include <iostream>
#include <iomanip>
#include <iomanip>

#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <math.h>
#include <assert.h>
#include <time.h>

using namespace std;

#define BLACK_BACKGROUND

#include "params.h"

#define read_size 25

int main(int argc, char* argv[]){
  
  srand( time(NULL) );
  int r_color = rand() % 256;
  int g_color = rand() % 256;
  int b_color = rand() % 256;

  ostringstream dummy_stream;
  dummy_stream<<r_color<<","<<g_color<<","<<b_color;
  
  string color_string = dummy_stream.str();
  
  params pars(argc, argv);
  pars.require("score_file","score_file",STRING_TYPE);
  pars.require("output_file","output_file",STRING_TYPE);
  
  pars.require("profile_threshold","profile_threshold", DOUBLE_TYPE);
  pars.require("normalizer", "normalizer", DOUBLE_TYPE);
  pars.require("contig_id","contig_id",STRING_TYPE);
  pars.optional("track_name", "track_name", "ChIP", STRING_TYPE);
  pars.optional("track_color", "track_color", color_string.c_str(), STRING_TYPE);

  if(!pars.enforce()) exit(1);

  string scores_fname = pars.get_string_value("score_file");
  string output_fname = pars.get_string_value("output_file");
  double profile_threshold = pars.get_double_value("profile_threshold");
  string contig_id = pars.get_string_value("contig_id");
  string track_name = pars.get_string_value("track_name");
  string track_color = pars.get_string_value("track_color");
  
  double normalizer = pars.get_double_value("normalizer");

  cout<<endl<<"contig: "<<contig_id<<endl;
  cout<<"score_file: "<<scores_fname<<endl;
  cout<<"profile_threshold: "<<profile_threshold<<endl;
  cout<<"normalizer: "<<normalizer<<endl;  
  cout<<"track_name: "<<track_name<<endl;
  cout<<"track_color: "<<track_color<<endl;

  remove(output_fname.c_str());
  ofstream out_str(output_fname.c_str());
  if(!out_str.good()){
    cerr<<"bad file name: "<<output_fname<<endl;
    exit(1);
  }

  //  out_str<<"browser positions "<<contig_id<<"1-10"<<endl;

  //out_str<<"track type=wiggle_0 name=\""<<track_name<<"\" description=\""<<track_name<<"\" ";
  //out_str<<"visibility=full color="<<track_color<<" priority=10 maxHeightPixels=50:50:11"<<endl;
  out_str<<"variableStep chrom="<<contig_id<<" span=1"<<endl;

  ifstream scores_str(scores_fname.c_str());
  if(!scores_str.good()){
    cerr<<"bad file: "<<scores_fname<<endl;
    exit(1);
  }

  cout<<"loading the scores"<<endl;

  unsigned int _size;
  int entries;

  scores_str.read((char*) (&_size), sizeof(_size));

  cout<<"_size: "<<_size<<endl;

  scores_str.read((char*) (&entries), sizeof(entries));

  cout<<"the scores file contains "<<_size<<" positions"<<endl;

  cout<<"score file : "<<entries<<" tags."<<endl;

  vector<float> scores(_size);

  cout<<"loading the score values...";
  cout.flush();
  scores_str.read((char*) &(scores[0]), sizeof(scores[0])*_size);
  cout<<"loaded."<<endl;

  scores_str.close();

  cout<<"outputting the wig file..."<<endl;

  int half_perc = _size / 200;
  
  stringstream tmp_stream;

  for(unsigned int i=0; i<_size; i++){
    if(i%half_perc == 0){
      int cur_perc = (int) (100.0 * ((double) i) / ((double) _size));
      cout<<"\rfinished "<<setprecision(1)<<cur_perc<<" %          ";
      cout.flush();
    }
    //    if(scores[i] >= profile_threshold){
    double cur_value = scores[i] / normalizer;
    if(cur_value > profile_threshold){
      tmp_stream<<i+1<<" "<<setprecision(3)<<cur_value<<endl;
    }
  }

  out_str<<tmp_stream.str();

  cout<<endl;
  
  out_str.close();
  return 0;
}
