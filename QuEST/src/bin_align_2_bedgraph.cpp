#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>

//#include <stdio.h>
//#include <conio.h>

#include <assert.h>

//using namespace std;

#include "params.h"
#include "string_utils.h"

using std::ifstream;
using std::ofstream;

using std::istringstream;
using std::cerr;
using std::endl;
using std::cout;

class genome_table{
public:
  vector<string> chr;
  vector<int> chr_size;

  unsigned int size();
  bool load(ifstream& istr);
  int chr_ind(const string& inp_chr);
};

unsigned int genome_table::size(){
  assert(chr.size() == chr_size.size());
  return chr.size();
	 
}

bool genome_table::load(ifstream& istr){
  if(!istr.good()) return false;
  //char sep = ' ';
  while(istr.good()){
    string cur_line;
    getline(istr, cur_line);
    if(istr.good()){
      //cout<<"read: "<<cur_line<<endl;
      istringstream cur_line_str(cur_line, istringstream::in);

      string cur_chr;
      int cur_chr_size;

      assert(cur_line_str.good());
      if(cur_line.length() > 0){
	if( cur_line.substr( 0, 1) != "#"){
	  cur_line_str>>cur_chr;
	  assert(cur_line_str.good());
	  cur_line_str>>cur_chr_size;
	  
	  chr.push_back(cur_chr);
	  chr_size.push_back(cur_chr_size);
	}
      }
    }
  }
  
  //cout<<endl;
  //cout<<"Read "<<size()<<" chromosomes"<<endl;
  
  return true;
}

int genome_table::chr_ind(const string& inp_chr){
  int res = -1;
  for(unsigned int i=0; i<chr.size(); i++){
    if(chr[i] == inp_chr) return i;
  }
  return res;
}

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("bin_align_file_prefix","bin_align_file_prefix",STRING_TYPE);
  pars.optional("bin_align_file_suffix","bin_align_file_suffix",".align.bin",STRING_TYPE);

  pars.require("output_file","output_file", STRING_TYPE);  
  pars.require("by_chr_prefix","by_chr_prefix",STRING_TYPE);

  pars.require("genome_table","genome_table", STRING_TYPE);
  pars.optional("track_name","track_name","bar_graph_track",STRING_TYPE);
  pars.optional("track_color","track_color","10,210,10",STRING_TYPE);
  pars.optional("track_priority","track_priority","0",STRING_TYPE);
  pars.optional("count_threshold","count_threshold","1",INT_TYPE);

  if(!pars.enforce()){
    exit(1);
  }


  string gt_fname = pars.get_string_value("genome_table");
  string output_fname = pars.get_string_value("output_file");
  string by_chr_prefix = pars.get_string_value("by_chr_prefix");

  string output_tmp_for_fname = output_fname + ".tmp.for";
  string output_tmp_rev_fname = output_fname + ".tmp.rev";

  string bin_align_file_prefix = pars.get_string_value("bin_align_file_prefix");
  string bin_align_file_suffix = pars.get_string_value("bin_align_file_suffix");
  
  string track_name = pars.get_string_value("track_name");
  string track_color = pars.get_string_value("track_color");
  string track_priority = pars.get_string_value("track_priority");

  int count_threshold = pars.get_int_value("count_threshold");

  ifstream gt_str(gt_fname.c_str());
  if(!gt_str.good()){
    cerr<<"Bad genome table: "<<gt_fname<<endl;
    exit(1);
  }
  
  genome_table gt;
  bool gt_loaded = gt.load(gt_str);
  
  if(!gt_loaded){
    cerr<<"failed to load genome table: bad file"<<endl;
    exit(1);
  }
  
  gt_str.close();
  
  
  double genome_size = 0;
  for(unsigned int i=0; i<gt.size(); i++){
    genome_size += (double) gt.chr_size[i];
  }

  remove(output_fname.c_str());
  remove(output_tmp_for_fname.c_str());
  remove(output_tmp_rev_fname.c_str());

  ofstream output_tmp_for_str(output_tmp_for_fname.c_str());
  if(!output_tmp_for_str.good()){
    cerr<<"Bad file: "<<output_tmp_for_fname.c_str();
    exit(1);
  }

  ofstream output_tmp_rev_str(output_tmp_rev_fname.c_str());
  if(!output_tmp_rev_str.good()){
    cerr<<"Bad file: "<<output_tmp_rev_fname.c_str();
    exit(1);
  }

  
  ofstream output_str(output_fname.c_str());
  if(!output_str.good()){
    cerr<<"Bad file: "<<output_fname.c_str();
    exit(1);
  }

  // writing a header

  string track_name_for = track_name + " + strand";
  string track_name_rev = track_name + " - strand";
  
  string track_color_for = "10,210,10";
  string track_color_rev = "10,10,210";
  

  output_tmp_for_str<<"track type=bedGraph name=\""<<track_name_for<<"\" description=\""<<track_name_for<<"\"";
  output_tmp_for_str<<" visibility=full color="<<track_color_for<<" priority="<<track_priority<<" maxHeightPixels=50:50:11";
  output_tmp_for_str<<" graphType=bar";
  output_tmp_for_str<<endl;

  output_tmp_rev_str<<"track type=bedGraph name=\""<<track_name_rev<<"\" description=\""<<track_name_rev<<"\"";
  output_tmp_rev_str<<" visibility=full color="<<track_color_rev<<" priority="<<track_priority<<" maxHeightPixels=50:50:11";
  output_tmp_rev_str<<" graphType=bar";
  output_tmp_rev_str<<endl;

  
  // /writing a header
  
  
  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chr = gt.chr[k];
    int cur_chr_size = gt.chr_size[k];

    string cur_chr_bedGraph_fname = by_chr_prefix + "." + gt.chr[k] + ".bedGraph";
    remove(cur_chr_bedGraph_fname.c_str());

    ofstream cur_chr_bedGraph_ofstr(cur_chr_bedGraph_fname.c_str());
    if(!cur_chr_bedGraph_ofstr.good()){
      cerr<<"Bad file name: "<<cur_chr_bedGraph_fname<<endl;
      exit(4);
    }
    
    string cur_chr_bedGraph_track_name_for = track_name + " " + gt.chr[k] + " + strand";
    string cur_chr_bedGraph_track_name_rev = track_name + " " + gt.chr[k] + " - strand";
    /*
    cur_chr_bedGraph_ofstr<<"track type=bedGraph name=\""<<cur_chr_bedGraph_track_name_for<<"\" description=\""<<track_name_for<<"\"";
    cur_chr_bedGraph_ofstr<<" visibility=full color="<<track_color_for<<" priority="<<track_priority<<" maxHeightPixels=50:50:11";
    cur_chr_bedGraph_ofstr<<" graphType=bar";
    cur_chr_bedGraph_ofstr<<endl;
    
    cur_chr_bedGraph_ofstr<<"track type=bedGraph name=\""<<cur_chr_bedGraph_track_name_rev<<"\" description=\""<<track_name_rev<<"\"";
    cur_chr_bedGraph_ofstr<<" visibility=full color="<<track_color_rev<<" priority="<<track_priority<<" maxHeightPixels=50:50:11";
    cur_chr_bedGraph_ofstr<<" graphType=bar";
    cur_chr_bedGraph_ofstr<<endl;
    */
    string bin_align_fname = bin_align_file_prefix + gt.chr[k] + bin_align_file_suffix;    

    ifstream bin_align_str(bin_align_fname.c_str());    

    if(!bin_align_str.good()){
      cerr<<"Bad file name: "<<bin_align_fname<<endl;
      exit(1);
    }
    
    int entries_pos;
    int entries_neg;
    
    int contig_name_size;
    
    bin_align_str.read((char*)(&contig_name_size), sizeof(contig_name_size));

    string cur_contig(contig_name_size, ' ');
    bin_align_str.read((char*)(&cur_contig[0]), contig_name_size*sizeof(char));

    bin_align_str.read((char*)(&entries_pos), sizeof(entries_pos));
    bin_align_str.read((char*)(&entries_neg), sizeof(entries_neg));

    vector<int> pos_hits(entries_pos);
    vector<int> neg_hits(entries_neg);

    bin_align_str.read((char*)(&pos_hits[0]), sizeof(int)*entries_pos);
    bin_align_str.read((char*)(&neg_hits[0]), sizeof(int)*entries_neg);

    bin_align_str.close();

    //cout<<cur_contig<<" +: "<<entries_pos<<" -: "<<entries_neg<<endl;

    assert(cur_contig == gt.chr[k]);
    
    vector<unsigned short> counts_for(cur_chr_size);
    vector<unsigned short> counts_rev(cur_chr_size);

    for(unsigned int i=0; i<counts_for.size(); i++){
      counts_for[i] = 0;
      counts_rev[i] = 0;
    }

    for(unsigned int i=0; i<pos_hits.size(); i++){
      /*
      if(pos_hits[i] < 0 || pos_hits[i] > cur_chr_size){
	cerr<<"Bad coordinate: "<<pos_hits[i]<<endl;
	cerr<<"Anticipated ragne for "<<gt.chr[k]<<" is [ 0, "<<cur_chr_size<<" )"<<endl;
	exit(1);
      }
      */
      counts_for[pos_hits[i]]++;      
    }
    for(unsigned int i=0; i<neg_hits.size(); i++){
      counts_rev[neg_hits[i]]++;
    }

    cur_chr_bedGraph_ofstr<<"track type=bedGraph name=\""<<cur_chr_bedGraph_track_name_for<<"\" description=\""<<track_name_for<<"\"";
    cur_chr_bedGraph_ofstr<<" visibility=full color="<<track_color_for<<" priority="<<track_priority<<" maxHeightPixels=50:50:11";
    cur_chr_bedGraph_ofstr<<" graphType=bar";
    cur_chr_bedGraph_ofstr<<endl;

    for(unsigned int i=0; i<counts_for.size(); i++){
      if(counts_for[i] >= count_threshold){
	output_tmp_for_str<<cur_chr<<" "<<i+1<<" "<<i+2<<" "<<counts_for[i]<<endl;
	cur_chr_bedGraph_ofstr<<cur_chr<<" "<<i+1<<" "<<i+2<<" "<<counts_for[i]<<endl;
      }      
    }

    cur_chr_bedGraph_ofstr<<endl;

    cur_chr_bedGraph_ofstr<<"track type=bedGraph name=\""<<cur_chr_bedGraph_track_name_rev<<"\" description=\""<<track_name_rev<<"\"";
    cur_chr_bedGraph_ofstr<<" visibility=full color="<<track_color_rev<<" priority="<<track_priority<<" maxHeightPixels=50:50:11";
    cur_chr_bedGraph_ofstr<<" graphType=bar";
    cur_chr_bedGraph_ofstr<<endl;

    for(unsigned int i=0; i<counts_rev.size(); i++){
      if(counts_rev[i] >= count_threshold){
	output_tmp_rev_str<<cur_chr<<" "<<i+1<<" "<<i+2<<" "<<counts_rev[i]<<endl;
	cur_chr_bedGraph_ofstr<<cur_chr<<" "<<i+1<<" "<<i+2<<" "<<counts_rev[i]<<endl;
      }      
    }
    cur_chr_bedGraph_ofstr.close();
  }
  
  output_tmp_for_str.close();
  output_tmp_rev_str.close();
  
  ifstream tmp_for_ifstr(output_tmp_for_fname.c_str());
  if(!tmp_for_ifstr.good()){
    cerr<<"Bad file: "<<output_tmp_for_fname.c_str();
    exit(1);
  }
  
  ifstream tmp_rev_ifstr(output_tmp_rev_fname.c_str());
  if(!tmp_rev_ifstr.good()){
    cerr<<"Bad file: "<<output_tmp_rev_fname.c_str();
    exit(1);
  }

  while(tmp_for_ifstr.good()){
    string dummy_line;
    getline(tmp_for_ifstr, dummy_line);
    if(tmp_for_ifstr.good()){
      output_str<<dummy_line<<endl;
    }
  }

  while(tmp_rev_ifstr.good()){
    string dummy_line;
    getline(tmp_rev_ifstr, dummy_line);
    if(tmp_rev_ifstr.good()){
      output_str<<dummy_line<<endl;
    }
  }

  cout<<"fusing bargraph files..."<<endl;
  
  tmp_for_ifstr.close();
  tmp_rev_ifstr.close();
  
  output_str.close();

  remove(output_tmp_for_fname.c_str());
  remove(output_tmp_rev_fname.c_str());
  
  return 0;
}
