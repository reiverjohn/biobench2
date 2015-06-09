#include <iostream>
#include <fstream>
#include <string>

#ifndef _ASSERT_H
#include "assert.h"
#endif

using std::string;
using std::cerr;
using std::endl;

#include "seq_contig.h"
#include "sequence_utility_functions.h"

#ifndef FAILED_TO_LOAD_CONTIG_FROM_FILE
#define FAILED_TO_LOAD_CONTIG_FROM_FILE -1
#endif

#ifndef SUCCEEDED_TO_LOAD_CONTIG_FROM_FILE
#define SUCCEEDED_TO_LOAD_CONTIG_FROM_FILE 1
#endif

namespace Error{
  struct bad_ifstream_in_the_function{
    const char* thrown_from;
    bad_ifstream_in_the_function(const char* _thrown_from){ 
      thrown_from = _thrown_from;
    }
  };
}

seq_contig::seq_contig(unsigned int _size){
   seq.reserve(_size);
}
void seq_contig::read_contig(std::ifstream& o){
  string dummy_str;
  getline(o, dummy_str);
  assert(fa_header(dummy_str));

  seq = "";
  header = dummy_str;
  
  bool continue_reading = true;
  while(continue_reading){
    getline(o, dummy_str);
    if(o.good() && !fa_header(dummy_str)){
      assert(valid_seq_string(dummy_str));
      seq += dummy_str;
    }
    else{ continue_reading = false; }
  }
}

void seq_contig::read_contig(const char* _fname){
  std::ifstream _seq_fname(_fname);
  if(!_seq_fname.good()){
    cerr<<"Bad filename: "<<_fname<<endl;
  }
  read_contig(_seq_fname);
}

int seq_contig::load_contig_from_file(const string& contig_id, std::ifstream& o){
  unsigned int contig_counter = 0;
  if(!o.good()) throw Error::bad_ifstream_in_the_function("int seq_contig::load_contig_from_file");

  std::ios::pos_type save_pos = o.tellg(); //save position to go back later
  o.seekg(0, std::ios::beg);
  
  std::ios::pos_type contig_start_pos;
  bool continue_search = true;
  while(continue_search){
    string dummy_string;
    std::ios::pos_type current_pos = o.tellg();
    getline(o, dummy_string);
    if(o.good()){
      if(dummy_string[0] == '>'){ //fasta header
	//std::cerr<<"header found: "<<dummy_string<<"\n";
	if(dummy_string.substr(1, dummy_string.length()-1) == contig_id){
	  continue_search = false;
	  contig_start_pos = current_pos;
	  contig_ind = contig_counter;
	}
	contig_counter++;
      }
    }
    else{
      continue_search = false;
      o.clear();
      o.seekg(save_pos);
      return FAILED_TO_LOAD_CONTIG_FROM_FILE;
    }
  }
  o.seekg(contig_start_pos);
  read_contig(o);
  o.seekg(save_pos);
  return SUCCEEDED_TO_LOAD_CONTIG_FROM_FILE;
}

seq_contig::seq_contig(std::ifstream& o){
  read_contig(o);
}

seq_contig::seq_contig(const char* _fname){
  read_contig(_fname);
}


void seq_contig::ag_recode(){
  string nseq = seq.substr(0,seq.size()-1);
  for(unsigned int i=0; i<seq.size()-1; i++){
    char cur_code = 'n';
    if(valid_nucl(seq[i]) && valid_nucl(seq[i]))
      cur_code = ag_encode(seq[i],seq[i+1]);
    nseq[i] = cur_code;
  }
  seq = nseq;
}

void seq_contig::save(std::ofstream& ostr){
  unsigned int str_length = 50;
  ostr<<header<<endl;
  for(unsigned int i=0; i<(seq.size()-seq.size()%str_length)/str_length; i++){
    int cur_start = i*str_length;
    string cur_str = seq.substr(cur_start, str_length);
    ostr<<cur_str<<endl;
  }
  if(seq.size()%str_length != 0){
    string cur_str = seq.substr(seq.size()-seq.size()%str_length, seq.size()%str_length);
    ostr<<cur_str<<endl;
  }
}

void seq_contig::save(std::ofstream& ostr, int region_start, int region_end){
  assert(region_start >= 0);
  assert(region_start < region_end);
  assert(region_end <= (int)(seq.length()));
  int region_size = region_end - region_start;
  unsigned int str_length = 50;
  ostr<<header<<endl;
  for(int i=0; i<(int)((region_size-region_size%str_length)/str_length); i++){
    int cur_start = region_start+i*str_length;
    string cur_str = seq.substr(cur_start, str_length);
    ostr<<cur_str<<endl;
  }
  if(seq.size()%str_length != 0){
    string cur_str = seq.substr(region_start + region_size - region_size%str_length,
				region_size%str_length);
    ostr<<cur_str<<endl;
  }
}

int seq_contig::size(){
  return seq.length();
}

void seq_contig::validate(){
  for(unsigned int i=0; i<seq.size(); i++){
    assert(valid_letter(seq[i]));
  }
}

char& seq_contig::operator[](unsigned int pos){
#ifdef DEBUG
  assert(pos < seq.size());
#endif
  return seq[pos];
}

void seq_contig::capitalize(){
  for(unsigned int i=0; i<seq.length(); i++){
    if(seq[i] == 'a'){
      seq[i] = 'A';
    }
    if(seq[i] == 'c'){
      seq[i] = 'C';
    }
    if(seq[i] == 'g'){
      seq[i] = 'G';
    }
    if(seq[i] == 't'){
      seq[i] = 'T';
    }
  }
}

