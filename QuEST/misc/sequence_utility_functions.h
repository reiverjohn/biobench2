#ifndef SEQUENCE_UTILITY_FUNCTIONS_H_
#define SEQUENCE_UTILITY_FUNCTIONS_H_

using std::string;

bool fa_header(string _str);
string ag_recode(const string& str);
bool valid_seq_string(string _str);
bool valid_nucl(char c);
bool valid_col(char c);
bool valid_letter(char c);
bool seq_missing(char c);
bool col_missing(char c);

char ag_encode(char c1, char c2);
string ag_encode(const string& str, int start, int length);
string ag_encode(const string& str);
short ag_encode_char(char c);

unsigned short encode_nucl(char c);
unsigned short encode_nucl_complement(char c);

char decode_nucl(unsigned int code);

unsigned int seq_space_tuple_index(const string& str, unsigned int starts_at_pos, unsigned int tuple_size);
unsigned int seq_space_tuple_index_rev_compl(const string& str, unsigned int starts_at_pos, unsigned int tuple_size);
unsigned int color_space_tuple_index(const string& str, int starts_at_pos, int tuple_size);
int color_space_tuple_index_rev(const string& str, int starts_at_pos, int tuple_size);
char capitalize(char c);
string capitalize(const string& _seq);
string rev_compl(const string& str);
char complement(char c);
string invert(const string& str);

namespace Error{
  struct Bad_ag_encode{
    char first_char, second_char;
    Bad_ag_encode(char _first, char _second){
      first_char = _first; second_char = _second;
    }
  };
  struct Bad_nucleotide{
    char c;
    Bad_nucleotide(char _c){ c=_c;}
  };
  struct Bad_fasta_letter{
    char c;
    Bad_fasta_letter(char _c){c=_c;}
  };
  struct Failed_to_compute_seq_tuple_index{
    string bad_tuple;
    Failed_to_compute_seq_tuple_index(const string& _str){
      bad_tuple = _str;
    }
  };
}

#endif
