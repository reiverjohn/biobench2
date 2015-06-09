#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <fstream>
//#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
//#include <list>
#include <math.h>
#include <assert.h>
//#include <time.h>

#include "string_utils.h"
//#include "seq_contig.h"
#include "params.h"

using namespace std;

struct mode{
  double value;
  double pos;
};

double kernel(double x, double bw){
  double c = 0.39894228;  
  return (c/bw)*exp(-x*x/(2.0*bw*bw));
}

int __int_min(int x, int y){
  if(x<y) return x;
  return y;
}

int __int_max(int x, int y){
  if(x>y) return x;
  return y;
}

double sqr(double x){
  return x*x;
}

double median(const vector<double>& vec){
  vector<double> tmp_vec = vec;
  sort(tmp_vec.begin(), tmp_vec.end());
  int cur_size = vec.size();
  assert(cur_size > 0);

  double res = 0;
  
  if(cur_size > 0){
    if(cur_size % 2 == 1){
      res = tmp_vec[(cur_size-1)/2+1];
    }
    
    else{
      res = 0.5 * (tmp_vec[cur_size/2-1] + tmp_vec[cur_size/2]);
    }
  }

  return res;
}

vector<mode> modes(const vector<double>& vec){
  double bw = 10;
  vector<mode> res;
  if(vec.size() == 0){ return res;}
  if(vec.size() == 1){
    mode new_mode;
    new_mode.value = kernel(0,bw);
    new_mode.pos = vec[0];
    res.push_back(new_mode);
    return res;
  }

  vector<double> tmp_vec = vec;
  sort(tmp_vec.begin(), tmp_vec.end());
  double min_value = tmp_vec[0];
  double max_value = tmp_vec[tmp_vec.size()-1];

  if(min_value == max_value){
    mode new_mode;
    new_mode.value = vec.size() * kernel(0,bw);
    new_mode.pos = vec[0];
    res.push_back(new_mode);
    return res;
  }

  int min_value_int = (int) floor(min_value);
  int max_value_int = (int) ceil(max_value);

  vector<double> density_estimate(max_value_int - min_value_int + 1);
  for(unsigned int i=0; i<density_estimate.size(); i++){
    density_estimate[i] = 0;
  }
  
  for(unsigned int i=0; i<tmp_vec.size(); i++){
    double cur_pos = tmp_vec[i];
    for(int j=min_value_int; j<=max_value_int; j++){
      density_estimate[j-min_value_int] += kernel((double)(cur_pos-j), bw);
    }
  }

  for(unsigned int i=0; i<density_estimate.size(); i++){
    bool this_is_mode = false;
    double cur_density_value = density_estimate[i];
    if(i!=0 && i !=density_estimate.size()-1){
      double next_density_value = density_estimate[i+1];
      double prev_density_value = density_estimate[i-1];
      if(cur_density_value > next_density_value && cur_density_value > prev_density_value){
	this_is_mode = true;
      }
    }
    else{
      if(i==0){
	double next_density_value = density_estimate[i+1];
	if(cur_density_value > next_density_value){
	  this_is_mode = true;
	}
      }
      else{
	if(i == density_estimate.size() -1){
	  double prev_density_value = density_estimate[i-1];
	  if(cur_density_value > prev_density_value){
	    this_is_mode = true;
	  }
	}
      }
    }

    if(this_is_mode){
      mode new_mode;
      new_mode.value = cur_density_value;
      new_mode.pos = i+min_value;
      res.push_back(new_mode);
    }
  }

  return res;
}

double mean(const vector<double>& vec){
  int vec_size = vec.size();
  assert(vec_size > 0);
  double res = 0;

  for(unsigned int i=0; i<vec.size(); i++){
    res += vec[i] / ((double) vec_size);    
  }
  return res;
}

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("genome_table","genome_table",STRING_TYPE);
  pars.require( "bin_align_path", "bin_align_file", STRING_TYPE);

  pars.require("regions_file","regions_file", STRING_TYPE);
  pars.require("output_file","output_file",STRING_TYPE);

  pars.optional("kde_bandwidth","kde_bandwidth","30",INT_TYPE);
  pars.optional("region_width","region_width", "300", INT_TYPE);
  pars.optional("dist_threshold","dist_threshold","300", INT_TYPE);
  pars.optional("correlation_lower_threshold","correlation_threshold","0.8",DOUBLE_TYPE);
  pars.optional("peak_shift_lower_threshold","peak_shift_lower_threshold","0",DOUBLE_TYPE);
  pars.optional("peak_shift_upper_threshold","peak_shift_upper_threshold","5000",DOUBLE_TYPE);
  pars.optional("os_upper_threshold","os_upper_threshold","10",DOUBLE_TYPE);
  pars.optional("top_regions_passed","top_regions_passed","200",INT_TYPE);
  pars.optional("estimation_method","estimation_method","median",STRING_TYPE);

  if(!pars.enforce()){
    pars.list_all_params();
    exit(1);
  }

  pars.list_all_params();
  string genome_table_fname = pars.get_string_value("genome_table");
  string ChIP_bin_align_path = pars.get_string_value("bin_align_path");

  string regions_fname = pars.get_string_value("regions_file");
  string output_fname = pars.get_string_value("output_file");

  int bandwidth =  pars.get_int_value("kde_bandwidth");
  int region_width = pars.get_int_value("region_width");
  int dist_threshold = pars.get_int_value("dist_threshold");

  double correlation_lower_threshold = pars.get_double_value("correlation_lower_threshold");
  double peak_shift_lower_threshold = pars.get_double_value("peak_shift_lower_threshold");
  double peak_shift_upper_threshold = pars.get_double_value("peak_shift_upper_threshold");
  double os_upper_threshold = pars.get_double_value("os_upper_threshold");
  
  int top_regions_passed = pars.get_int_value("top_regions_passed");

  string estimation_method = pars.get_string_value("estimation_method");
  
  cout<<"genome_table: "<<genome_table_fname;
  cout<<"ChIP_bin_align_path : "<<ChIP_bin_align_path<<endl;
  
  cout<<"regions file : "<<regions_fname<<endl;
  cout<<"output file : "<<output_fname<<endl;

  cout<<"bandwidth: "<<bandwidth<<endl;
  cout<<"region_width: "<<region_width<<endl;
  cout<<"dist_threshold: "<<dist_threshold<<endl;

  ifstream genome_table_str(genome_table_fname.c_str());
  if(!genome_table_str.good()){
    cerr<<"bad file name: "<<genome_table_fname<<endl;
    exit(1);
  }

  vector<string> contigs;
  vector<int> contig_sizes;

  char gap_symbol = ' ';
  while(genome_table_str.good()){
    string dummy_string;
    getline(genome_table_str, dummy_string);
    if(genome_table_str.good()){
      if(dummy_string.length() > 0){
        if(dummy_string[0] != '#'){
          vector<string> cur_contig_fields = split(dummy_string, gap_symbol);
          if(cur_contig_fields.size() >= 2){
            string cur_contig_id = cur_contig_fields[0];
            int cur_contig_size  = atoi(cur_contig_fields[1].c_str());
            assert(cur_contig_size >= 0);

            contigs.push_back(cur_contig_id);
            contig_sizes.push_back(cur_contig_size);
          }
        }
      }
    }
  }

  genome_table_str.close();

  vector< vector<int> > ChIP_read_coordinates_for;
  vector< vector<int> > ChIP_read_coordinates_rev;

  for(unsigned int i=0; i<contigs.size(); i++){
    string cur_contig = contigs[i];

    string bin_align_fname = ChIP_bin_align_path + "/" + cur_contig + ".align.bin";

    ifstream bin_align_str(bin_align_fname.c_str());
    if(!bin_align_str.good()){
      cerr<<"bad file name: "<<bin_align_fname<<endl;
      exit(1);
    }

    int entries_pos;
    int entries_neg;

    int contig_name_size;

    bin_align_str.read((char*)(&contig_name_size), sizeof(contig_name_size));

    string cur_ChIP_contig(contig_name_size, ' ');
    bin_align_str.read((char*)(&cur_ChIP_contig[0]), contig_name_size*sizeof(char));

    assert(cur_contig == cur_ChIP_contig);

    bin_align_str.read((char*)(&entries_pos), sizeof(entries_pos));
    bin_align_str.read((char*)(&entries_neg), sizeof(entries_neg));

    vector<int> pos_hits(entries_pos);
    vector<int> neg_hits(entries_neg);

    bin_align_str.read((char*)(&pos_hits[0]), sizeof(int)*entries_pos);
    bin_align_str.read((char*)(&neg_hits[0]), sizeof(int)*entries_neg);

    bin_align_str.close();

    cout<<cur_contig<<" +: "<<entries_pos<<" -: "<<entries_neg<<endl;

    ChIP_read_coordinates_for.push_back( pos_hits);
    ChIP_read_coordinates_rev.push_back( neg_hits);

  }
  /*
  ifstream ChIP_align_str(ChIP_align_fname.c_str());
  if(!ChIP_align_str.good()){
    cerr<<"bad file name: "<<ChIP_align_fname<<endl;
    exit(1);
  }
  */
  /*
  ifstream RX_noIP_align_str(RX_noIP_align_fname.c_str());
  if(!RX_noIP_align_str.good()){
    cerr<<"bad file name: "<<RX_noIP_align_fname<<endl;
    exit(1);
  }
  */

  ifstream regions_str(regions_fname.c_str());
  if(!regions_str.good()){
    cerr<<"bad file name: "<<regions_fname<<endl;
    exit(1);
  }

  remove(output_fname.c_str());
  ofstream out_str(output_fname.c_str());
  if(!out_str.good()){
    cerr<<"bad file name: "<<output_fname<<endl;
    exit(1);
  }

  cout<<"\nreading alignments...\n"<<endl;
  
  /*
  //vector<string> contig_names;
  vector< vector<int> > ChIP_read_coordinates_for;
  vector< vector<int> > ChIP_read_coordinates_rev;

  int alignment_counter = 0;
  char sep = ' ';

  while(ChIP_align_str.good()){
    
    if(alignment_counter % 10000 == 0){
      printf("\rread %.2f M alignments", ((double)alignment_counter)/(1000000.0));
      cout.flush();
    }

    string cur_alignment;
    getline(ChIP_align_str, cur_alignment);

    vector<string> cur_alignment_fields = split(cur_alignment, sep);

    //    cout<<cur_alignment_fields[0]<<endl;

    if(cur_alignment_fields.size() == 3){
      
      string cur_contig = cur_alignment_fields[0];
      int cur_coord = atoi(cur_alignment_fields[1].c_str());      
      string cur_orient = cur_alignment_fields[2];

      //      cout<<cur_contig<<endl;
      alignment_counter++;

      bool contig_found = false;
      int contig_ind = -1;
      for(unsigned int i=0; i<contig_names.size(); i++){
	if(contig_names[i] == cur_contig){
	  contig_found = true;
	  contig_ind = i;
	}
      }

      if(!contig_found){
	  contig_names.push_back(cur_contig);
	  vector<int> dummy_vec;
	  ChIP_read_coordinates_for.push_back(dummy_vec);
	  ChIP_read_coordinates_rev.push_back(dummy_vec);
	  contig_ind = ChIP_read_coordinates_for.size() - 1;
      }

      if(cur_orient == "+"){
	ChIP_read_coordinates_for[contig_ind].push_back(cur_coord);
      }
      else if(cur_orient == "-"){
	ChIP_read_coordinates_rev[contig_ind].push_back(cur_coord);
      }
    }
  }

  assert(ChIP_read_coordinates_for.size() == ChIP_read_coordinates_rev.size());
  
  cout<<endl;
  ChIP_align_str.close();

  cout<<endl;

  cout<<"sorting coordinates...";
  cout.flush();

  for(unsigned int i=0; i<ChIP_read_coordinates_for.size(); i++){
    sort(ChIP_read_coordinates_for[i].begin(), ChIP_read_coordinates_for[i].end());
    sort(ChIP_read_coordinates_rev[i].begin(), ChIP_read_coordinates_rev[i].end());
  }

  cout<<"done!"<<endl;
  
  cout<<"found the following "<<contig_names.size()<<" contigs:"<<endl<<endl;;
  for(unsigned int i=0; i<contig_names.size(); i++){
    cout<<contig_names[i]<<" +: "<<ChIP_read_coordinates_for[i].size()<<" -: ";
    cout<<ChIP_read_coordinates_rev[i].size()<<endl;
  }

  cout<<endl;
  cout<<"processing regions"<<endl;
  */

  char space_char = ' ';
  //char dash_char = '-';

  int passed_regions_counter = 0;
  vector<double> peak_shifts;

  while(regions_str.good() && passed_regions_counter < top_regions_passed){
    string cur_locus;
    getline(regions_str, cur_locus);
    vector<string> cur_locus_fields = split(cur_locus, space_char);

    if(regions_str.good()){
      if(cur_locus_fields.size() > 0){
	string cur_locus_id = cur_locus_fields[0];
	string cur_contig;
	int start_coord, end_coord;
	//if(cur_locus_id.substr(0, 1) == "R" || cur_locus_id.substr(0, 1) == "P"){
	//region encountered, contig in the field 1, coordinate range in the field 2
	
	assert(cur_locus_fields.size() >= 3);
	cur_contig = cur_locus_fields[0];
	
	int cur_locus_pos = atoi(cur_locus_fields[1].c_str());
	
	//if(cur_locus_id.substr(0, 1) == "R"){
	// cur_locus_pos = atoi(cur_locus_fields[8].c_str());
	//}
	//else{
	// cur_locus_pos = atoi(cur_locus_fields[2].c_str());
	//}
	
	start_coord = cur_locus_pos - bandwidth*3 - region_width - dist_threshold;
	end_coord = cur_locus_pos + bandwidth*3 + region_width + dist_threshold;	 	  

	int start_coord_for_os = cur_locus_pos - region_width/2;
	int end_coord_for_os = cur_locus_pos + region_width/2;
	
	int region_size = end_coord - start_coord; //both ends are inclusive, so keep track of <=       
	int os_region_size = end_coord_for_os - start_coord_for_os + 1;
	
	int center_coord = region_size / 2 + 1;
	
	vector<int> start_counts_for(region_size+1);
	vector<int> start_counts_rev(region_size+1);

	vector<int> start_counts_os_for(os_region_size);
	vector<int> start_counts_os_rev(os_region_size);

	vector<double> profile_for(region_size+1);
	vector<double> profile_rev(region_size+1);
	

	for(unsigned int i=0; i<start_counts_for.size(); i++){
	  start_counts_for[i] = 0;
	  start_counts_rev[i] = 0;
	  
	  profile_for[i] = 0.0;
	  profile_rev[i] = 0.0;
	}

	for(unsigned int i=0; i<start_counts_os_for.size(); i++){
	  start_counts_os_for[i] = 0;
	  start_counts_os_rev[i] = 0;
	}
	
	bool contig_found = false;
	int contig_ind = -1;
	
	for(unsigned int i=0; i<contigs.size(); i++){
	  if(contigs[i] == cur_contig){
	    contig_found = true;
	    contig_ind = i;
	  }
	}
	
	if( !(contig_ind < (int)ChIP_read_coordinates_for.size() && contig_ind >=0) ){
	  cout<<"failed to locate the following region: "<<endl;
	  cout<<cur_contig<<" "<<cur_locus_pos<<endl;
	  exit(1);
	}
	//assert(contig_ind < (int)ChIP_read_coordinates_for.size() && contig_ind >=0);
	
	if(contig_found){
	  for(unsigned int i=0; i<ChIP_read_coordinates_for[contig_ind].size(); i++){
	    int cur_read_coord = ChIP_read_coordinates_for[contig_ind][i];
	    
	    if(cur_read_coord >= start_coord && cur_read_coord <= end_coord){
	      start_counts_for[cur_read_coord - start_coord]++;
	    }
	    if(cur_read_coord >= start_coord_for_os && cur_read_coord <= end_coord_for_os){
	      start_counts_os_for[cur_read_coord - start_coord_for_os]++;
	    }
	  }
	  
	  for(unsigned int i=0; i<ChIP_read_coordinates_rev[contig_ind].size(); i++){
	    int cur_read_coord = ChIP_read_coordinates_rev[contig_ind][i];
	    
	    if(cur_read_coord >= start_coord && cur_read_coord <= end_coord){
	      start_counts_rev[cur_read_coord - start_coord]++;
	    }
	    if(cur_read_coord >= start_coord_for_os && cur_read_coord <= end_coord_for_os){
	      start_counts_rev[cur_read_coord - start_coord_for_os]++;
	    }
	  }
	}
	else{
	  cerr<<"failed to find contig"<<endl;
	  //print the failed entry into the file
	}
	
	//calculating the profile;
	for(int i=0; i<(int)start_counts_for.size(); i++){
	  if(start_counts_for[i] > 0){
	    for(int j=__int_max(0, i-5*bandwidth); 
		j<__int_min(start_counts_for.size()-1, i+5*bandwidth); j++){
	      //for(int j=0; j<(int)start_counts_for.size(); j++){
	      profile_for[j] += start_counts_for[i]*kernel((double)(j-i), bandwidth);
	    }
	  }
	  if(start_counts_rev[i] > 0){
	    for(int j=__int_max(0, i-5*bandwidth); 
		j<__int_min(start_counts_for.size()-1, i+5*bandwidth); j++){
	      profile_rev[j] += start_counts_rev[i]*kernel((double)(j-i), bandwidth);
	    }
	  }
	}
	
	//calculating the correlation
	
	vector<double> rho(dist_threshold+1);
	
	for(int i=-dist_threshold/2; i<=dist_threshold/2; i++){
	  int start_for = center_coord-region_width/2-i;
	  int end_for = center_coord+region_width/2-i+1;
	  
	  int start_rev = center_coord-region_width/2+i;
	  int end_rev = center_coord+region_width/2+i+1;
	  
	  double mu_x=0, var_x=0; //forward strand
	  double mu_y=0, var_y=0; //reverse strand
	  
	  for(int j=start_for; j<end_for; j++){
	    mu_x += profile_for[j] / ((double) region_width + 1);
	  }
	  
	  for(int j=start_rev; j<end_rev; j++){
	    mu_y += profile_rev[j] / ((double) region_width + 1);
	  }
	  
	  for(int j=start_for; j<end_for; j++){
	    var_x += sqr(profile_for[j]-mu_x) / ((double) region_width + 1);
	  }
	  
	  for(int j=start_rev; j<end_rev; j++){
	    var_y += sqr(profile_rev[j]-mu_y) / ((double) region_width + 1);
	  }
	  
	  double s_x = sqrt(var_x);
	  double s_y = sqrt(var_y);
	  
	  
	  double cur_rho = 0;
	  
	  for(int j=0; j<=region_width; j++){
	    cur_rho += (profile_for[start_for + j] - mu_x) * 
	      (profile_rev[start_rev+j] - mu_y) / (s_x*s_y*((double)region_width+1.0));
	  }	
	  
	  rho[dist_threshold/2 + i] = cur_rho;
	  
	  //out_str<<i<<"\t"<<mu_x<<"\t"<<s_x<<"\t"<<mu_y<<"\t"<<s_y<<"\t"<<cur_rho<<endl;
	  
	  //cout<<"ps: "<<i<<" mu_x: "<<mu_x<<" s_x: "<<s_x;
	  //cout<<" mu_y: "<<mu_y<<" s_y: "<<s_y<<" rho: "<<cur_rho<<endl;
	}	  
	
	//exit(0);
	
	double rho_max = 0;
	int rho_max_ind = 0;
	
	for(int i=0; i<(int) rho.size(); i++){
	  if(i==0){
	    rho_max = rho[i];
	    rho_max_ind = i;
	  }
	  else if(rho_max < rho[i]){
	    rho_max = rho[i];
	    rho_max_ind = i;
	  }
	}
	
	int peak_shift = (rho_max_ind - dist_threshold/2);
	
	/*
	  for(unsigned int i=0; i<profile_for.size(); i++){
	  out_str<<i<<"\t"<<setprecision(3)<<profile_for[i]<<"\t"<<profile_rev[i]<<endl;
	  }
	  
	  exit(0);
	*/
	// begin calculate oversampling
	double cur_os = -1.0; //instead of +inf

	int reads_in_the_region = 0;
	int pos_with_reads = 0;
	int pos = 2*os_region_size;

	assert(pos > 0);

	for(unsigned int i=0; i<start_counts_os_for.size(); i++){
	  if(start_counts_os_for[i] > 0){
	    pos_with_reads++;
	    reads_in_the_region += start_counts_os_for[i];
	  }
	}
	for(unsigned int i=0; i<start_counts_os_rev.size(); i++){
	  if(start_counts_os_rev[i] > 0){
	    pos_with_reads++;
	    reads_in_the_region += start_counts_os_rev[i];
	  }
	}
	
	if(pos_with_reads > 0 && pos > 0){
	  double lambda = ((double) reads_in_the_region) / ((double)pos);
	  double nu = lambda / ( 1.0 - exp( reads_in_the_region * log( 1.0-1.0/((double)pos) ) ) );
	  double mu = ((double)reads_in_the_region) / ((double)pos_with_reads);
	  
	  cur_os = mu / nu;

	  bool region_passed = false;
	  if(cur_os <= os_upper_threshold &&
	     rho_max >= correlation_lower_threshold &&
	     peak_shift <= peak_shift_upper_threshold &&
	     peak_shift >= peak_shift_lower_threshold){
	    peak_shifts.push_back(peak_shift);
	    passed_regions_counter++;
	    region_passed = true;
	  }
	  
	  
	  // end calculate oversampling
	  
	  //if(cur_os >= 0){
	  if(region_passed){
	    cout<<passed_regions_counter<<": "<<cur_locus<<" ps: "<<peak_shift<<" cor: "<<setprecision(3)<<rho_max<<" os: "<<cur_os<<endl;
	    out_str<<cur_locus<<" ps: "<<peak_shift<<" cor: "<<setprecision(3)<<rho_max<<" os: "<<cur_os<<endl;
	  }
	  //else{
	  //  cout<<cur_locus<<" ps: "<<peak_shift<<" cor: "<<setprecision(3)<<rho_max<<" os: inf"<<endl;
	  //  out_str<<cur_locus<<" ps: "<<peak_shift<<" cor: "<<setprecision(3)<<rho_max<<" os: inf"<<endl;
	  //}
	}
	
	//exit(0);
      }
      else{
	//out_str<<cur_locus<<endl;
      }
    }
    else{
      //out_str<<cur_locus<<endl;
    }
  }

  double ps_mean = -1;
  if (peak_shifts.size() > 0){
    ps_mean = mean(peak_shifts);
  }
  double ps_median = -1;
  if(peak_shifts.size() > 0){
    ps_median = median(peak_shifts);
  }
  vector<mode> ps_modes = modes(peak_shifts);

  double largest_mode_pos = -1;
  double largest_mode_value = -1;
  bool largest_mode_found = false;
  
  
  if(ps_modes.size() > 0){
    largest_mode_pos = ps_modes[0].pos;
    largest_mode_value = ps_modes[0].value;
    largest_mode_found = true;
  }
  for(unsigned int i=1; i<ps_modes.size(); i++){
    double cur_mode_pos = ps_modes[i].pos;
    double cur_mode_value = ps_modes[i].value;

    if(cur_mode_value > largest_mode_value){
      largest_mode_value = cur_mode_value;
      largest_mode_pos = cur_mode_pos;
    }
  }


  out_str<<"used_regions: "<<peak_shifts.size()<<endl;
  if(peak_shifts.size() > 0){
    for(unsigned int i=0; i<ps_modes.size(); i++){
      out_str<<"mode: "<<ps_modes[i].pos<<" value: "<<ps_modes[i].value<<endl;
      cout<<"mode: "<<ps_modes[i].pos<<" value: "<<ps_modes[i].value<<endl;
    }
    out_str<<endl;

    out_str<<"peak_shift mean: "<<ps_mean<<endl;
    out_str<<"peak_shift median: "<<ps_median<<endl;
    if(largest_mode_found){
      out_str<<"mode: "<<largest_mode_pos<<endl;
    }
    else{
      out_str<<"mode: NA"<<endl;
    }

    if(estimation_method == "mean"){
      out_str<<"peak_shift_estimate: "<<ps_mean<<endl;    
    }
    if(estimation_method == "median"){
      out_str<<"peak_shift_estimate: "<<ps_median<<endl;    
    }
    if(estimation_method == "mode"){
      if(largest_mode_found){
	out_str<<"peak_shift_estimate: "<<largest_mode_pos<<endl;    
      }
      else{
	out_str<<"peak_shift_estimate: NA"<<endl;    
      }
    }
  }
  else{
    out_str<<"peak_shift mode: NA"<<endl;
    out_str<<"peak_shift mean: NA"<<endl;
    out_str<<"peak_shift median: NA"<<endl;
    out_str<<"peak_shift_estimate: NA"<<endl;    
  }
  out_str<<"estimation_method: "<<estimation_method<<endl;


  out_str.close();

  if(peak_shifts.size() > 0){    
    cout<<"peak shift mean: "<<ps_mean<<endl;
    cout<<"peak shift median: "<<ps_median<<endl;
    if(largest_mode_found){
      cout<<"peak_shift_estimate: "<<largest_mode_pos<<endl;    
    }
    else{
      cout<<"peak_shift_estimate: NA"<<endl;    
    }
  }
  else{
    cout<<"peak shift mean: NA"<<endl;
    cout<<"peak shift median: NA"<<endl;
  }
  
  return 0;

}
