#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <assert.h>

using std::vector;
using std::string;
using std::cerr;
using std::endl;

//#include "utils.h"
/*
std::vector<string> split(const string& inp_string, const string& sep){
  std::vector<string> res;
  if(inp_string.length() == 0)  return res;
  if(inp_string.length() < sep.length()){
    res.push_back(inp_string);
    return res;
  }
  unsigned int last_i = 0;
  for(unsigned int i=0; i<inp_string.length()-sep.length(); i++){
    //std::cout<<"comp: "<<inp_string.substr(i,sep.length())<<" to "<<sep<<std::endl;
    if(inp_string.substr(i,sep.length()) == sep){
      if(i>last_i+1){ 
	res.push_back(inp_string.substr(last_i, i-last_i));
      }
	last_i = i+1;
	//}
    }
  }
  if(last_i < inp_string.length()-1) 
    res.push_back(inp_string.substr(last_i, inp_string.length()-last_i));

  if(res.size() == 0 && inp_string[0] != sep[0]) res.push_back(inp_string);
  return res;
}
*/
std::vector<string> split(const string& inp_string, char sep){
  std::vector<string> res;
  if(inp_string.length() == 0)  return res;

  unsigned int last_i = 0;
  for(unsigned int i=0; i<inp_string.length()/*-1*/; i++){
    //if(inp_string.substr(i,sep.length()) == sep){
    if(inp_string[i] == sep){
      if(i>last_i){ 
	res.push_back(inp_string.substr(last_i, i-last_i));
      }
      last_i = i+1;
      //}
    }
  }
  if(last_i < inp_string.length()) 
    res.push_back(inp_string.substr(last_i, inp_string.length()-last_i));

  if(res.size() == 0 && inp_string[0] != sep) res.push_back(inp_string);
  return res;
}

unsigned int get_string_count(std::ifstream& ifs){
  using std::ios;

  assert(ifs.good());
  std::ios::pos_type cur_pos = ifs.tellg();
  ifs.seekg(0, ios::beg);

  string dummy_string;
  unsigned int string_counter = 0;
  while(ifs.good()){
    getline(ifs,dummy_string);
    if(ifs.good()) string_counter++;    
  }

  ifs.clear();
  ifs.seekg(cur_pos);
  return string_counter;
}

void chomp(string& str){
  //erases bad characters
  bool continue_chomp = true;

  while(continue_chomp){
    switch(str[str.length()-1]){
    case char(13):{ str.erase(str.length()-1,1); break;}
    default:{ continue_chomp = false; break; }
    }
  }  
}


std::ios::pos_type perc_file_length(std::ifstream& ifs){
  //assert(ifs.good());
  if(!ifs.good()){ cerr<<"bad file in perc_file_length"<<endl; std::terminate();}
  std::ios::pos_type cur_pos = ifs.tellg();
  ifs.seekg(0, std::ios::end);
  std::ios::pos_type _length = ifs.tellg();
  ifs.seekg(cur_pos);
  return (_length - _length%100) / 100;
}
