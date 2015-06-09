#include <string>
#include <iostream>

#include <assert.h>

using std::string;
using std::cerr;
using std::endl;

#include "sequence_utility_functions.h"

static unsigned int pow4[] =
  {1,4,16,64,256,1024,4096,16384,65536,262144,1048576,
   4194304,16777216,67108864,268435456,1073741824};
/*
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
*/
char complement(char c){
  switch(c){
  case 'a':{ return 't'; break;}
  case 't':{ return 'a'; break;}
  case 'g':{ return 'c'; break;}
  case 'c':{ return 'g'; break;}
  case 'A':{ return 'T'; break;}
  case 'T':{ return 'A'; break;}
  case 'G':{ return 'C'; break;}
  case 'C':{ return 'G'; break;}
  case 'n':{ return 'n'; break;}
  case 'N':{ return 'N'; break;}
  case 'M':{ return 'M'; break;}
  case 'm':{ return 'm'; break;}
  case 'R':{ return 'R'; break;}
  case 'r':{ return 'r'; break;}
  default:{ cerr<<"bad char in complement: "<<c<<endl;
  assert(false); return 'Q'; break;}
  }
  //assert(valid_letter(c));                                                     
  //char nc = capitalize(c);
  /*
  if(nc == 'N') return 'N';
  if(nc == 'T') return 'A';
  if(nc == 'A') return 'T';
  if(nc == 'C') return 'G';
  if(nc == 'G') return 'C';
  assert(false);
  */
}


bool small(char c){
  if(c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'n') return true;
  return false;
}

string capitalize(const string& _seq){
  string res = _seq;
  for(unsigned int i=0; i<_seq.length(); i++){
    switch(_seq[i]){
    case 'a':{res[i] = 'A'; break;}
    case 'A':{break;}
    case 't':{res[i] = 'T'; break;}
    case 'T':{break;}
    case 'g':{res[i] = 'G'; break;}
    case 'G':{break;}
    case 'c':{res[i] = 'C'; break;}
    case 'C':{break;}      
    case 'n':{res[i] = 'N'; break;}
    case 'N':{break;}
    default:{ throw Error::Bad_nucleotide(_seq[i]); break;}
    }
  }
  return res;
}
char capitalize(char c){
  if(c == 'n') return 'N';
  if(c == 'a') return 'A';
  if(c == 't') return 'T';
  if(c == 'g') return 'G';
  if(c == 'c') return 'C';

  if(c == 'N') return 'N';
  if(c == 'A') return 'A';
  if(c == 'T') return 'T';
  if(c == 'G') return 'G';
  if(c == 'C') return 'C';

  //cerr<<c;
  assert(false);
}

string rev_compl(const string& _str){
  //string dummy_str = _str;                                                                                                      
  //reverse(dummy_str.begin(), dummy_str.end());                                                                                  
  string res;
  for(int i=_str.size()-1; i>=0; i--){
    char cur_char = _str[i];
    //cout<<cur_char;                                                                                                             
    //assert(false);                                                                                                              
    res.push_back(complement(cur_char));
  }
  return res;
}

bool valid_nucl(char c){
  if(c == 'a' || c == 'A' || c == 't' || c == 'T' ||
     c == 'g' || c == 'G' || c == 'c' || c == 'C')
    return true;
  else return false;
}
bool valid_col(char c){
  if(c == '0' || c == '1' || c == '2' || c == '3')
    return true;
  else return false;
}
bool valid_letter(char c){
  if(c == 'a' || c == 'A' || c == 't' || c == 'T' ||
     c == 'g' || c == 'G' || c == 'c' || c == 'C' ||
     c == 'n' || c == 'N' || c == 'm' || c == 'M' ||
     c == 'r' || c == 'R')
    return true;
  else return false;
}

bool seq_missing(char c){
  //if(c== 'n' || c== 'N' || c=='m' || c=='M'|| c) return true;  
  if(!valid_nucl(c)) return true;
  return false;
}

bool col_missing(char c){
  //if(c== 'n' || c== 'N' || c=='m' || c=='M'|| c) return true;  
  if(! (valid_col(c))) return true;
  return false;
}

bool fa_header(string _str){
  if(_str[0] == '>') return true;
  return false;
}
string ag_recode(const string& str){
  string nseq =str.substr(0,str.size()-1);
  for(unsigned int i=0; i<str.size()-1; i++){
    char cur_code = 'n';
    if(valid_nucl(str[i]) && valid_nucl(str[i]))
      cur_code = ag_encode(str[i],str[i+1]);
    nseq[i] = cur_code;
  }
  return nseq;
}
bool valid_seq_string(string _str){
  for(unsigned int i=0; i<_str.size(); i++){
    if(!valid_letter(_str[i])){

      cerr<<"warning: bad letter ["<<_str[i]<<"] code: "<<(int) _str[i]<<" at position "<<i<<endl;
      cerr<<" in the expected sequence string "<<_str<<endl;
      throw(Error::Bad_fasta_letter(_str[i]));
      return false;
    }
  }
  return true;
}

char ag_encode(char c1, char c2){
  if(c1 == 'N' || c1 == 'n' || c2 == 'N' || c2 == 'n') return 'N';

  if(((c1 == 'A' || c1 == 'a') && (c2 == 'A' || c2 == 'a')) ||
     ((c1 == 'C' || c1 == 'c') && (c2 == 'C' || c2 == 'c')) ||
     ((c1 == 'G' || c1 == 'g') && (c2 == 'G' || c2 == 'g')) ||
     ((c1 == 'T' || c1 == 't') && (c2 == 'T' || c2 == 't'))) return '0';

  if(((c1 == 'A' || c1 == 'a') && (c2 == 'C' || c2 == 'c')) ||
     ((c1 == 'C' || c1 == 'c') && (c2 == 'A' || c2 == 'a')) ||
     ((c1 == 'G' || c1 == 'g') && (c2 == 'T' || c2 == 't')) ||
     ((c1 == 'T' || c1 == 't') && (c2 == 'G' || c2 == 'g'))) return '1';

  if(((c1 == 'A' || c1 == 'a') && (c2 == 'G' || c2 == 'g')) ||
     ((c1 == 'C' || c1 == 'c') && (c2 == 'T' || c2 == 't')) ||
     ((c1 == 'G' || c1 == 'g') && (c2 == 'A' || c2 == 'a')) ||
     ((c1 == 'T' || c1 == 't') && (c2 == 'C' || c2 == 'c'))) return '2';

  if(((c1 == 'A' || c1 == 'a') && (c2 == 'T' || c2 == 't')) ||
     ((c1 == 'C' || c1 == 'c') && (c2 == 'G' || c2 == 'g')) ||
     ((c1 == 'G' || c1 == 'g') && (c2 == 'C' || c2 == 'c')) ||
     ((c1 == 'T' || c1 == 't') && (c2 == 'A' || c2 == 'a'))) return '3';

  throw Error::Bad_ag_encode(c1,c2);

  //cerr<<c1<<c2<<endl;
  //assert(false);
  //return '?';
}string ag_encode(const string& str, int start, int length){
  assert(start >= 0 && (unsigned int) (start+length) <= str.length());
  string res;
  for(int i=start; i<start+length-1; i++){
    try{ res += ag_encode(str[i], str[i+1]); }
    catch( Error::Bad_ag_encode e ){
      cerr<<"failed to encode ( "<<e.first_char<<","<<e.second_char<<" ) ";
      cerr<<" at pos: "<<i<<endl;
      //skip();
      //assert(false);
    }
  }
  return res;
}

string ag_encode(const string& str){
  return ag_encode(str, 0, str.length());
}

unsigned short encode_nucl(char c){
  switch(c){
  case 'A':{ return 0; break;}
  case 'a':{ return 0; break;}
  case 'C':{ return 1; break;}
  case 'c':{ return 1; break;}
  case 'G':{ return 2; break;}
  case 'g':{ return 2; break;}
  case 'T':{ return 3; break;}
  case 't':{ return 3; break;}
  case 'N':{ return 0; break;}
  case '.':{ return 0; break;}
  case 'n':{ return 0; break;}
  default:{ throw Error::Bad_nucleotide(c);}
  }
}
unsigned short encode_nucl_complement(char c){
  switch(c){
  case 'A':{ return 3; break;}
  case 'a':{ return 3; break;}
  case 'C':{ return 2; break;}
  case 'c':{ return 2; break;}
  case 'G':{ return 1; break;}
  case 'g':{ return 1; break;}
  case 'T':{ return 0; break;}
  case 't':{ return 0; break;}
  case 'N':{ return 3; break;}
  case 'n':{ return 3; break;}
  case '.':{ return 3; break;}
  default:{ throw Error::Bad_nucleotide(c);}
  }
}

char decode_nucl(unsigned int code){
  assert(code < 4);
  switch(code){
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'T';
  }
  return -1;
}

unsigned int seq_space_tuple_index(const string& str, unsigned int starts_at_pos, unsigned int tuple_size){
  //  cerr<<"starts: "<<starts_at_pos<<" t_s: "<<tuple_size<<" l: "<<str.length()<<endl;
  //a-0, c-1, g-2, t-3
  assert(starts_at_pos >= 0);
  assert(starts_at_pos + tuple_size <= (str.length()));

  unsigned int res = 0;

  for(unsigned int i = 0; i<tuple_size; i++){
    try{
      res += encode_nucl(str[starts_at_pos+i])*pow4[i];
    }
    catch(Error::Bad_nucleotide e){
      cerr<<"seq_space_tuple_index: Error [ "<<e.c<<" ] in encode char at pos: "<<starts_at_pos<<" ";
      for(unsigned int j=0; j<tuple_size; j++){
	cerr<<str[starts_at_pos + j];
        //cerr<<endl;
	//throw Error::Failed_to_compute_seq_tuple_index(str.substr(starts_at_pos, tuple_size));
      }
      cerr<<endl;
      throw Error::Failed_to_compute_seq_tuple_index(str.substr(starts_at_pos, tuple_size));
    }
  }

  assert(res < pow4[tuple_size]);
  return res;
}

unsigned int seq_space_tuple_index_rev_compl(const string& str, unsigned int starts_at_pos, unsigned int tuple_size){
  //  cerr<<"starts: "<<starts_at_pos<<" t_s: "<<tuple_size<<" l: "<<str.length()<<endl;
  assert(starts_at_pos >= 0);
  assert(starts_at_pos < (str.length()));
  assert(starts_at_pos - tuple_size + 1 >= 0);

  unsigned int res = 0;

  for(unsigned int i = 0; i<tuple_size; i++){
    try{
      //setbits(res, str[starts_at_pos-i],i);
      res += /*ag_encode_char*/encode_nucl_complement(str[starts_at_pos-i])*pow4[i];//int_pow(4,i);            
    }
    catch(Error::Bad_nucleotide e){
      cerr<<"seq_space_tuple_index_rev_compl: Error [ "<<e.c<<" ] in encode char at pos: "<<starts_at_pos<<" ";
      for(unsigned int j=0; j<tuple_size; j++){
	cerr<<str[starts_at_pos - j];
      }
      cerr<<endl;
    }
  }

  assert(res < pow4[tuple_size]);
  return res;
} 

unsigned int color_space_tuple_index(const string& str, int starts_at_pos, int tuple_size){
  //  cerr<<"starts: "<<starts_at_pos<<" t_s: "<<tuple_size<<" l: "<<str.length()<<endl;
  assert(starts_at_pos >= 0);
  assert(starts_at_pos + tuple_size <= (int)(str.length()));

  unsigned int res = 0;

  for(int i = 0; i<tuple_size; i++){
    try{
      res += ag_encode_char(str[starts_at_pos+i])*pow4[i];//int_pow(4,i);                                         
    }
    catch(int a){
      if(a == 20){
        cerr<<"Error in encode char at pos: "<<starts_at_pos<<" ";
        for(int j=0; j<tuple_size; j++){
          cerr<<str[starts_at_pos + j];
        }
        cerr<<endl;
	throw 21;
      }
    }
  }

  assert(res < pow4[tuple_size]);
  return res;
}

int color_space_tuple_index_rev(const string& str, int starts_at_pos, int tuple_size){
  //  cerr<<"starts: "<<starts_at_pos<<" t_s: "<<tuple_size<<" l: "<<str.length()<<endl;
  assert(starts_at_pos >= 0);
  assert(starts_at_pos < (int)(str.length()));
  assert(starts_at_pos - tuple_size + 1 >= 0);

  int res = 0;

  for(int i = 0; i<tuple_size; i++){
    try{
      res += ag_encode_char(str[starts_at_pos-i])*pow4[i];//int_pow(4,i); 
    }
    catch(Error::Bad_ag_encode e){
      cerr<<"Failed to color encode: "<<e.first_char<<" "<<e.second_char<<endl;
      assert(false);
    }
  }
  return res;
}

short ag_encode_char(char c){
  switch(c){
  case '0':{ return 0; break;}
  case '1':{ return 1; break;}
  case '2':{ return 2; break;}
  case '3':{ return 3; break;}
  default:{throw 20; return -1; break;assert(false); break;} //you shouldn't be here
  }
}


string invert(const string& str){
  string res;
  for(int i=str.length()-1; i>=0; i--){
    res += str[i];
  }
  return res;
}
