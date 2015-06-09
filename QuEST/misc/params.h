#ifndef PARAMS_H
#define PARAMS_H

#include <vector>
#include <string>
#include <cstdlib>

#ifndef param_type
#define param_type unsigned int
#endif

#ifndef PARAM_TYPES
#define PARAM_TYPES

#define CHAR_TYPE 0
#define STRING_TYPE 10
#define INT_TYPE 20
#define DOUBLE_TYPE 30

#endif

using std::vector;
using std::string;

class params{
 private:
  //unsigned int argc;
  //vector<string> argv;

  vector<string> required_params;
  vector<string> required_params_legend;
  vector<param_type> required_params_type;

  vector<string> optional_params;
  vector<string> optional_params_legend;
  vector<param_type> optional_params_type;
  vector<string> optional_params_default_value;

  vector<string> par_names;
  vector<string> par_values;

  void parse(const char* this_par, string& cur_par_name, string& cur_par_value);
 public:
  params(int _argc, char* _argv[]);
  void require(const char* this_par_name, const char* this_par_legend, param_type this_par_type); //specify what sort of parameters are required
  void require(const string& par_name, const string& par_legend, param_type par_type); //same
  void optional(const char* this_par_name, const char* this_par_legend, const char* default_value, param_type par_type);

  void list_all_params();
  bool enforce(); //enforces the required parameters

  string get_string_value(const char* par_name);
  char get_char_value(const char* par_name);
  int get_int_value(const char* par_name);
  double get_double_value(const char* par_name);
};

#endif
