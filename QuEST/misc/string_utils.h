#ifndef UTILS_H
#define UTILS_H

std::vector<std::string> split(const std::string& inp_string, const std::string& sep);
std::vector<std::string> split(const std::string& inp_string, char sep);

void chomp(std::string& str);
std::ios::pos_type perc_file_length(std::ifstream& ifs);
unsigned int get_string_count(std::ifstream& ifs);
#endif
