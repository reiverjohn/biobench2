#ifndef SEQ_CONTIG_H
#define SEQ_CONTIG_H

#define FAILED_TO_LOAD_CONTIG_FROM_FILE -1
#define SUCCEEDED_TO_LOAD_CONTIG_FROM_FILE 1


class seq_contig{
public:
  std::string header;
  std::string seq;
  unsigned int contig_ind; //the index of contig in the file it's read from

  char& operator[](unsigned int pos);

  void read_contig(std::ifstream& o);
  void read_contig(const char* fname);
  
  int size();
  void validate();
  seq_contig(){}
  seq_contig(unsigned int _size);//{ seq.reserve(_size);}
  seq_contig(std::string _header, std::string _seq){ header = _header; seq = _seq; validate();}
  seq_contig(std::ifstream& o);
  seq_contig(const char* fname);
  
  void ag_recode();
  void save(std::ofstream& ostr);
  void save(std::ofstream& ostr, int region_start, int region_end);
  int load_contig_from_file(const std::string& contig_id, std::ifstream& o);
  int get_contig_ind();//const std::string& contig_id, std::ifstream& o);
  void capitalize();
};

#endif
