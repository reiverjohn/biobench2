#include <cstdlib>

#include <string>
#include <iostream>
#include <vector>

#include "assert.h"


using std::string;
using std::vector;



#include "params.h"

namespace Error{
  struct bad_argc{
    int argc;
    bad_argc(int _argc){ argc = _argc; }
  };
  struct bad_param{
    const char* value;
    bad_param(const char* _value){ value = _value; }
  };
  struct bad_param_type{
    param_type bp_type;
    bad_param_type(param_type _bp_type){ bp_type = _bp_type;}
  };
  struct param_missing{
    const char* par_missing;
    param_missing(const char* _par_missing){ par_missing = _par_missing; }
  };
}

params::params(int _argc, char* _argv[]){
  if(_argc < 1) throw Error::bad_argc(_argc);

  par_names.resize(_argc);
  par_values.resize(_argc);

  par_names[0] = "executable";
  par_values[0] = _argv[0];
  
  for(int i=1; i<_argc; i++){
    string cur_par_name="";
    string cur_par_value="";
    try{ 
      parse(_argv[i], cur_par_name, cur_par_value); 
      par_names[i] = cur_par_name;
      par_values[i] = cur_par_value;
    }
    catch(Error::bad_param bp){
      std::cerr<<"Warning:\ncouldn't parse the parameter: "<<bp.value<<"\nskipping\n";
      par_names[i] = "";
      par_values[i] = "";
    }    
  }
  //std::cout<<"this can print\n";
}

void params::parse(const char* this_par, string& cur_par_name, string& cur_par_value){
  string dummy_string = this_par;
  unsigned int eq_pos = 0;
  bool eq_pos_found = false;

  for(unsigned int i=0; i<dummy_string.length() && !eq_pos_found; i++){
    if(dummy_string[i] == '='){
      eq_pos_found = true;
      eq_pos = i;
    }
  }
  
  if(!eq_pos_found || eq_pos == dummy_string.length()-1){
    throw Error::bad_param(this_par);
  }

  cur_par_name = dummy_string.substr(0, eq_pos);
  cur_par_value = dummy_string.substr(eq_pos+1, dummy_string.length()-eq_pos-1);

}

void params::require(const string& par_name, const string& par_legend, param_type par_type){
  required_params.push_back(par_name);
  required_params_legend.push_back(par_legend);
  required_params_type.push_back(par_type);
}

void params::require(const char* this_par_name, const char* this_par_legend, 
		     param_type this_par_type){
  string this_par_name_str = this_par_name;
  string this_par_legend_str = this_par_legend;
  
  this->require(this_par_name_str, this_par_legend_str, this_par_type);
}

void params::optional(const char* this_par_name, const char* this_par_legend, 
		      const char* default_value, param_type par_type){
  string _this_par_name = this_par_name;
  string _this_par_legend = this_par_legend;
  string _default_value = default_value;
  
  optional_params.push_back(_this_par_name);
  optional_params_legend.push_back(_this_par_legend);
  optional_params_type.push_back(par_type);
  optional_params_default_value.push_back(default_value);

  bool par_found = false;
  unsigned int par_index = 0;
  for(unsigned int i=0; i<par_names.size(); i++){
    if(par_names.at(i) == _this_par_name){
      par_found = true;
      par_index = i;
    }
  }
  if(!par_found){
    par_names.push_back(_this_par_name);
    par_values.push_back(default_value);  
  }
}

string stringify_type(param_type par_type){
  string dummy_string;
  switch(par_type){
  case CHAR_TYPE:
    dummy_string = "CHAR_TYPE";    
    break;
    
  case STRING_TYPE:
    dummy_string = "STRING_TYPE"; 
    break;

  case INT_TYPE:
    dummy_string = "INT_TYPE";    
    break;

  case DOUBLE_TYPE:
    dummy_string = "DOUBLE_TYPE";    
    break;
  default: 
    std::cerr<<"bad par_type: "<<par_type<<"\n";
    throw Error::bad_param_type(par_type);
  }
  return dummy_string;
}

bool good_type(const string& cur_par, param_type req_type){
  switch(req_type){
  case CHAR_TYPE:
    if(cur_par.length() == 1) return true; break;
  case STRING_TYPE:
    return true; break;
  case INT_TYPE:
    return true; break;
  case DOUBLE_TYPE:
    return true; break;
  default:
    std::cerr<<"bad param type: "<<req_type<<"\n";
    throw Error::bad_param_type(req_type);
    break;
  }
  return false;
}

bool params::enforce(){
  bool all_enforced = true;
  for(unsigned int i=0; i<required_params.size(); i++){
    string cur_req_param = required_params.at(i);
    string cur_par_legend = required_params_legend.at(i);
    param_type cur_req_type = required_params_type.at(i);

    //std::cout<<"enforcing "<<cur_req_param<<"\n";

    bool param_found = false;
    unsigned int param_pos = 0;

    for(unsigned j=0; j<par_names.size() && !param_found; j++){
      if(par_names.at(j) == cur_req_param){
	param_found = true;
	param_pos = j;
      }
    }
    if(!param_found){
      all_enforced = false;
      std::cout<<"    Required:\t"<<cur_req_param<<" = <"<<cur_par_legend<<">\t"<<stringify_type(cur_req_type)<<" [ parameter missing ]\n";
    }
    else{ //type_check now
      string cur_par_value = par_values.at(param_pos);
      if(!good_type(cur_par_value, cur_req_type)){
	all_enforced = false;
	std::cout<<"    Required:\t"<<cur_req_param<<" = <"<<cur_par_legend<<">\t"<<stringify_type(cur_req_type)<<" [ unexpected type ]\n";
      }
    }
  }
  assert(optional_params.size() == optional_params_legend.size());
  assert(optional_params.size() == optional_params_default_value.size());
  assert(optional_params.size() == optional_params_type.size());

  /*
  for(unsigned int i=0; i<optional_params.size(); i++){
    string cur_param = optional_params.at(i);
    string cur_legend = optional_params_legend.at(i);
    string cur_default = optional_params_default_value.at(i);
    param_type cur_type = optional_params_type.at(i);
    
    bool param_found = false;
    unsigned int param_pos = 0;
    
    for(unsigned j=0; j<par_names.size() && !param_found; j++){
      if(par_names.at(j) == cur_param){
	param_found = true;
	param_pos = j;
      }
    }
    
    //std::cout<<"param_found: "<<param_found<<"\n";

    if(!param_found){
      std::cerr<<"    Optional:\t"<<cur_param<<" = <"<<cur_legend<<">\t"<<stringify_type(cur_type)<<" [ default: "<<cur_default<<" ]\n";
    }
    else{ //type_check now
      string cur_par_value = par_values.at(param_pos);
      if(!good_type(cur_par_value, cur_type)){
	all_enforced = false;
	std::cerr<<"    Optional:\t"<<cur_param<<" = <"<<cur_legend<<">\t"<<stringify_type(cur_type)<<" [ unexpected type ]\n";
      }
    }
  }
  */
  
  return all_enforced;
}

void params::list_all_params(){
  assert(par_names.size() == par_values.size());
  std::cout<<"Specified parameters:\n";
  for(unsigned int i=0; i<par_names.size(); i++){
    std::cout<<par_names.at(i)<<"\t"<<par_values.at(i)<<"\n";
  }
}

int params::get_int_value(const char* par_name){
  string dummy_string = par_name;
  for(unsigned int i=0; i<par_names.size(); i++){
    if(par_name == par_names.at(i)) return atoi(par_values.at(i).c_str());
  }
#ifdef DEBUG
  std::cerr<<"bad par_name: "<<par_name<<"\n";
#endif
  throw Error::param_missing(par_name);
}

double params::get_double_value(const char* par_name){
  string dummy_string = par_name;
  for(unsigned int i=0; i<par_names.size(); i++){
    if(par_name == par_names.at(i)) return atof(par_values.at(i).c_str());
  }
#ifdef DEBUG
  std::cerr<<"bad par_name: "<<par_name<<"\n";
#endif
  throw Error::param_missing(par_name);
}

string params::get_string_value(const char* par_name){
  string dummy_string = par_name;
  for(unsigned int i=0; i<par_names.size(); i++){
    if(par_name == par_names.at(i)) return par_values.at(i);
  }
#ifdef DEBUG
  std::cerr<<"bad par_name: "<<par_name<<"\n";
#endif
  throw Error::param_missing(par_name);
}
char params::get_char_value(const char* par_name){
  string dummy_string = par_name;
  for(unsigned int i=0; i<par_names.size(); i++){
    if(par_name == par_names.at(i)){
      string cur_par_value = par_values.at(i);
      return cur_par_value.at(0);
    }
  }
  throw Error::param_missing(par_name);
}
