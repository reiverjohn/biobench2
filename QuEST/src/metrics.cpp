#include <iostream>
#include <iomanip>
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

double _sqr(double x){
  return x*x;
}

double kernel_sum_square(double sigma){
  double res = 0;
  for(int i= - 5*((int)sigma); i<5*((int)sigma); i++){
    res += _sqr(1.0/sqrt(2.0*3.1415*sigma*sigma) * exp(-((double)i*i)/(2*sigma*sigma)));
  }
  return res;
}

double Poisson_p_value(int k, double lambda){
  //double precision = 0.00000001;                                                                                               

  //double delta = 1;

  //double cur_sum = 0;
  //double prev_sum = 0;
  //double cur_i = k;

  //  while(delta > precision){
  
  double res = 1.0;
  for(int i=0; i<k; i++){
    double cur_i = i;
    
    double log_cur_p = -lambda + cur_i*log(lambda) ;
    for(int j=1; j<=cur_i; j++){
      log_cur_p = log_cur_p - log((double)j);
    }

    double cur_p = exp(log_cur_p);
    //cur_sum += cur_p;
    //delta = cur_p;//cur_sum - prev_sum;                                                                                          

    //prev_sum = cur_sum;
    //cur_i++;
    //cerr<<"i: "<<i<<" cur_p: "<<cur_p<<" res: "<<res<<endl;
    res -= cur_p;
  }
  return res;
}

double n(double x){

  // !!! this is the possible performance bottleneck

  double A = 1.0/sqrt(2.0 * 3.1415);
  return A * exp(-x*x*0.5);

}

double gaussian_tail_weight(double x){
  double a1 = 0.4361836;
  double a2 = -0.1201676;
  double a3 = 0.9372980;

  
  double k = 1.0/(1.0 + (0.33267 * x));
  
  if (x >= 0.0){
      return n(x)* (a1*k + (a2*k*k) + (a3*k*k*k));
  }
  assert(false);
  return -1;
}

double log_gaussian_tail_weight(double x){  
  double a1 = 0.4361836;

  if(x>=3){
    //return lognx + log(a1) + log(k);
    return (-0.5*log(2.0*3.1415) - x*x*0.5 ) + log(a1) - log(1.0+0.33267*x);
  }
  if(x>=0 && x<3){
    return log(gaussian_tail_weight(x));
  }

  assert(false);
  return -1;
}

double log_Poisson_pv_approximation(int k, double lambda){
  double tmp_var = ((double)k) + 0.5 - lambda;
  if(k<3*sqrt(lambda) || k<5 || tmp_var <= 0){
    return log(Poisson_p_value(k, lambda));
  }

  return log_gaussian_tail_weight(((double)k + 0.5 - lambda) / sqrt(lambda) );
}

double Poisson_pv_approximation(int k, double lambda){
  return gaussian_tail_weight(((double)k + 0.5 - lambda) / sqrt(lambda) );
}


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

double double_max(double x, double y){
  if (x > y) return x;
  return y;
}

double double_min(double x, double y){
  if(x<y) return x;
  return y;
}

double sqr(double x){
  return x*x;
}

struct QuEST_region{  
  string entry_string;
  
  string QuEST_id;

  string chrom;

  int start_coord;
  int end_coord;

  int rank;

  int ps;
  double ps_cor;

  int ChIP_tags;
  int background_tags; 
  double basal_tags;

  double ChIP_value;
  double background_value;

  string max_pos;
  string ef;

  int region_size_used;

  double log_pv;
  double log_qv;
  

  double tag_ef;
};

struct QuEST_peak{
  string entry_string;

  string QuEST_id;
  
  string chrom;
  
  int coord;
  
  string region;
  string ef;
  
  int rank;
  
  double ChIP_value;
  double background_value;

  double log_pv;
  double log_qv;

  int ps;
  double ps_cor;
};

struct QuEST_call{
  //consists of QuEST region and peaks contained within
  QuEST_region region;
  vector<QuEST_peak> peaks;
};

int main(int argc, char* argv[]){

  params pars(argc, argv);

  //  pars.optional("mappable_genome_fraction","mappable_genome_fraction","0.7", DOUBLE_TYPE);
  pars.require("genome_table","genome_table",STRING_TYPE);
  
  pars.require("ChIP_bin_align_path", "ChIP_bin_align_path", STRING_TYPE);
  pars.require("background_bin_align_path", "background_bin_align_path", STRING_TYPE);

  pars.require("QuEST_calls_file","QuEST_calls_file", STRING_TYPE);
  pars.require("output_file","output_file",STRING_TYPE);
  pars.require("report_file","report_file",STRING_TYPE);
  
  pars.require("ChIP_tags", "ChIP_tags", INT_TYPE);
  pars.require("background_tags", "background_tags", INT_TYPE);

  pars.require("ChIP_basal_level", "ChIP_basal_level", DOUBLE_TYPE);
  pars.require("ChIP_basal_tag_rate", "ChIP_basal_tag_rate", DOUBLE_TYPE);

  pars.optional("bandwidth","bandwidth","30",INT_TYPE);
  pars.optional("region_width","region_width", "300", INT_TYPE);
  pars.optional("dist_threshold","dist_threshold","300", INT_TYPE); //for calculating the peak shift
  

  if(!pars.enforce()){
    pars.list_all_params();
    exit(1);
  }

  pars.list_all_params();
  
  string genome_table_fname = pars.get_string_value("genome_table");
  //  double mappable_genome_fraction = pars.get_double_value("mappable_genome_fraction");
  
  string ChIP_bin_align_path = pars.get_string_value("ChIP_bin_align_path");
  string background_bin_align_path = pars.get_string_value("background_bin_align_path");

  string QuEST_calls_fname = pars.get_string_value("QuEST_calls_file");
  string output_fname = pars.get_string_value("output_file");
  string report_fname = pars.get_string_value("report_file");

  int bandwidth =  pars.get_int_value("bandwidth");
  int region_width = pars.get_int_value("region_width");
  int dist_threshold = pars.get_int_value("dist_threshold");

  int ChIP_tags = pars.get_int_value("ChIP_tags");
  int background_tags = pars.get_int_value("background_tags");

  double ChIP_basal_level = pars.get_double_value("ChIP_basal_level");
  double ChIP_basal_tag_rate = pars.get_double_value("ChIP_basal_tag_rate");

  string correction_type = "Benjamini-Hochberg";

  double background_basal_level = pars.get_double_value("background_basal_level");

  bool background_used = false;
  if(background_tags > 0){
    background_used = true;
  }
  else{
    background_used = false;
  }

  int ChIP_tags_within_regions = 0;
  int background_tags_within_regions = 0;
  
  //  double background_basal_tag_rate = pars.get_double_value("background_basal_tag_rate");

  // p-value parameters
  
  // the maximum of background estimate rates 
  // is calculated based on windows below
  int w1 = region_width;
  int w2 = 10*region_width;
  int w3 = 10000;   

  cout<<"genome_table: "<<genome_table_fname<<endl;
  //cout<<"mappable genome fraction: "<<mappable_genome_fraction<<endl;
  cout<<"ChIP_tags: "<<ChIP_tags<<endl;
  cout<<"background_tags: "<<background_tags<<endl;

  cout<<"ChIP_bin_align_path : "<<ChIP_bin_align_path<<endl;
  cout<<"background_bin_align_path : "<<background_bin_align_path<<endl;  

  cout<<"QuEST calls file : "<<QuEST_calls_fname<<endl;
  cout<<"output file : "<<output_fname<<endl;

  cout<<"bandwidth: "<<bandwidth<<endl;
  cout<<"region_width: "<<region_width<<endl;
  cout<<"dist_threshold: "<<dist_threshold<<endl;

  vector<string> contig_names;
  vector<int> contig_sizes;

  ifstream genome_table_str(genome_table_fname.c_str());
  if(!genome_table_str.good()){
    cerr<<"bad file name: "<<genome_table_fname<<endl;
    exit(1);
  }

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

            contig_names.push_back(cur_contig_id);
            contig_sizes.push_back(cur_contig_size);
          }
        }
      }
    }
  }

  genome_table_str.close();


  double genome_size = 0;
  for(int i=0; i<(int)contig_sizes.size(); i++){
    genome_size += ((double)contig_sizes[i]);
  }


  //double mappable_genome_fraction = 0.7;

  vector< vector<int> > ChIP_coordinates_for;
  vector< vector<int> > ChIP_coordinates_rev;

  for(unsigned int i=0; i<contig_names.size(); i++){
    string cur_contig = contig_names[i];

    string ChIP_bin_align_fname = ChIP_bin_align_path + "/" + cur_contig + ".align.bin";

    ifstream ChIP_bin_align_str(ChIP_bin_align_fname.c_str());
    if(!ChIP_bin_align_str.good()){
      cerr<<"bad file name: "<<ChIP_bin_align_fname<<endl;
      exit(1);
    }

    int entries_pos;
    int entries_neg;

    int contig_name_size;

    ChIP_bin_align_str.read((char*)(&contig_name_size), sizeof(contig_name_size));

    string cur_ChIP_contig(contig_name_size, ' ');
    ChIP_bin_align_str.read((char*)(&cur_ChIP_contig[0]), contig_name_size*sizeof(char));

    assert(cur_contig == cur_ChIP_contig);

    ChIP_bin_align_str.read((char*)(&entries_pos), sizeof(entries_pos));
    ChIP_bin_align_str.read((char*)(&entries_neg), sizeof(entries_neg));

    vector<int> pos_hits(entries_pos);
    vector<int> neg_hits(entries_neg);

    ChIP_bin_align_str.read((char*)(&pos_hits[0]), sizeof(int)*entries_pos);
    ChIP_bin_align_str.read((char*)(&neg_hits[0]), sizeof(int)*entries_neg);

    ChIP_bin_align_str.close();

    cout<<cur_contig<<" +: "<<entries_pos<<" -: "<<entries_neg<<endl;

    ChIP_coordinates_for.push_back( pos_hits);
    ChIP_coordinates_rev.push_back( neg_hits);

  }

  vector< vector<int> > background_coordinates_for;
  vector< vector<int> > background_coordinates_rev;
  
  for(unsigned int i=0; i<contig_names.size(); i++){
    string cur_contig = contig_names[i];
    
    string background_bin_align_fname = background_bin_align_path + "/" + cur_contig + ".align.bin";
    
    ifstream background_bin_align_str(background_bin_align_fname.c_str());
    if(!background_bin_align_str.good()){
      cerr<<"bad file name: "<<background_bin_align_fname<<endl;
      exit(1);
    }

    int entries_pos;
    int entries_neg;

    int contig_name_size;

    background_bin_align_str.read((char*)(&contig_name_size), sizeof(contig_name_size));

    string this_contig(contig_name_size, ' ');
    background_bin_align_str.read((char*)(&this_contig[0]), contig_name_size*sizeof(char));

    assert(cur_contig == this_contig);

    background_bin_align_str.read((char*)(&entries_pos), sizeof(entries_pos));
    background_bin_align_str.read((char*)(&entries_neg), sizeof(entries_neg));

    vector<int> pos_hits(entries_pos);
    vector<int> neg_hits(entries_neg);

    background_bin_align_str.read((char*)(&pos_hits[0]), sizeof(int)*entries_pos);
    background_bin_align_str.read((char*)(&neg_hits[0]), sizeof(int)*entries_neg);

    background_bin_align_str.close();

    cout<<cur_contig<<" +: "<<entries_pos<<" -: "<<entries_neg<<endl;

    background_coordinates_for.push_back( pos_hits);
    background_coordinates_rev.push_back( neg_hits);
  }
  
  double ChIP_2_background_tag_ratio = ((double) ChIP_tags) / ((double) background_tags);
  
  ifstream QuEST_calls_str(QuEST_calls_fname.c_str());
  if(!QuEST_calls_str.good()){
    cerr<<"bad file name: "<<QuEST_calls_fname<<endl;
    exit(1);
  }

  //reading QuEST calls

  char space_char = ' ';
  char dash_char = '-';

  vector<QuEST_call> ChIP_calls;

  while(QuEST_calls_str.good()){
    string dummy_string;
    getline(QuEST_calls_str, dummy_string);
    if(QuEST_calls_str.good()){
      if(dummy_string.length() > 0){
	string first_char = dummy_string.substr(0,1);

	if(first_char != "#"){
	  if(first_char == "R"){	    
	    // reading a region
	    //QuEST_region cur_region;
	    QuEST_call cur_ChIP_call;

	    cur_ChIP_call.region.entry_string = dummy_string;
	    vector<string> cur_region_fields = split(dummy_string, space_char);
	    assert(cur_region_fields.size()  >= 3);
	    cur_ChIP_call.region.QuEST_id = cur_region_fields[0];
	    cur_ChIP_call.region.chrom = cur_region_fields[1];
	    
	    vector<string> cur_region_range = split(cur_region_fields[2], dash_char);
	    assert(cur_region_range.size() == 2);
	    
	    cur_ChIP_call.region.start_coord = atoi(cur_region_range[0].c_str());
	    cur_ChIP_call.region.end_coord = atoi(cur_region_range[1].c_str());	    	    	    

	    cur_ChIP_call.region.ChIP_value = atof(cur_region_fields[4].c_str());
	    cur_ChIP_call.region.background_value = atof(cur_region_fields[6].c_str());
	    
	    cur_ChIP_call.region.max_pos = cur_region_fields[8];
	    cur_ChIP_call.region.ef = cur_region_fields[10];
	    
	    ChIP_calls.push_back(cur_ChIP_call);
	  }
	  if(first_char == "P"){
	    assert(ChIP_calls.size() > 0);
	    // reading a peak
	    
	    QuEST_peak cur_peak;
	    cur_peak.entry_string = dummy_string;

	    vector<string> cur_peak_fields = split(dummy_string, space_char);
	    
	    assert(cur_peak_fields.size() >= 7);

	    cur_peak.QuEST_id = cur_peak_fields[0];
	    cur_peak.chrom = cur_peak_fields[1];	    	    
	    cur_peak.coord = atoi(cur_peak_fields[2].c_str());
	    
	    cur_peak.ChIP_value = atof(cur_peak_fields[4].c_str());
	    cur_peak.background_value = atof(cur_peak_fields[6].c_str());

	    cur_peak.region = cur_peak_fields[8];
	    cur_peak.ef = cur_peak_fields[10];

	    assert(cur_peak.chrom == ChIP_calls[ChIP_calls.size()-1].region.chrom);
	    assert(cur_peak.coord >= ChIP_calls[ChIP_calls.size()-1].region.start_coord);
	    assert(cur_peak.coord <= ChIP_calls[ChIP_calls.size()-1].region.end_coord);

	    ChIP_calls[ChIP_calls.size() - 1].peaks.push_back(cur_peak);
	    
	  }
	}
      }
    }
  }

  
  int peaks_read = 0;
  int regions_read = 0;
  for(unsigned int i=0; i<ChIP_calls.size(); i++){
    regions_read++;
    peaks_read += ChIP_calls[i].peaks.size();
  }
  
  cout<<endl;
  cout<<"regions read: "<<regions_read<<endl;
  cout<<"peaks read: "<<peaks_read<<endl;
  cout<<endl;

  remove(output_fname.c_str());
  ofstream out_str(output_fname.c_str());
  if(!out_str.good()){
    cerr<<"bad file name: "<<output_fname<<endl;
    exit(1);
  }

  remove(report_fname.c_str());
  ofstream report_str(report_fname.c_str());
  if(!report_str.good()){
    cerr<<"bad file name: "<<report_fname<<endl;
    exit(1);
  }

  cout<<"calculating peak p-values"<<endl;
  for(unsigned int k=0; k<ChIP_calls.size(); k++){
    for(unsigned int l=0; l<ChIP_calls[k].peaks.size(); l++){
      double cur_ChIP_value = ChIP_calls[k].peaks[l].ChIP_value;

      double ChIP_mu = ChIP_basal_level;
      double ChIP_sigma_sq = ChIP_mu * kernel_sum_square(bandwidth);
      double ChIP_z = (cur_ChIP_value - ChIP_mu ) / sqrt(ChIP_sigma_sq);
      double cur_log_pv = log_gaussian_tail_weight(ChIP_z);

      if(background_tags > 0){
	double cur_background_value = ChIP_calls[k].peaks[l].background_value;
	double background_mu = background_basal_level;
	double background_sigma_sq = background_mu * kernel_sum_square(bandwidth);
	double background_z = (cur_background_value - background_mu) / sqrt(background_sigma_sq);
	double new_log_pv = -1;
	if(background_z >= 3.0){
	  ChIP_mu = cur_background_value * ((double)ChIP_tags) / ((double)background_tags);
	  ChIP_sigma_sq = ChIP_mu * kernel_sum_square(bandwidth);
	  double new_ChIP_z = (cur_ChIP_value - ChIP_mu ) / sqrt(ChIP_sigma_sq);
	  if(ChIP_z > 0){
	    new_log_pv = log_gaussian_tail_weight(new_ChIP_z);
	  }
	  else{
	    new_log_pv = 0;
	  }
	  //cur_log_pv = new_log_pv;
	}

	//cout<<"ChIP: "<<cur_ChIP_value<<" ChIP_z: "<<ChIP_z<<" background_z: "<<background_z<<" cur_log_pv: "<<cur_log_pv<<" new_log_pv: "<<new_log_pv<<" kss: "<<kernel_sum_square(30)<<endl;
	if(background_z >= 3.0){
	  cur_log_pv = new_log_pv;
	}
      }
      else{      
	//cout<<"ChIP: "<<cur_ChIP_value<<" ChIP_z: "<<ChIP_z<<" cur_log_pv: "<<cur_log_pv<<" kss: "<<kernel_sum_square(30)<<endl;
      }
      // char c;
      //     cin>>c;

            /*
      double lambda_basal = ((double)ChIP_tags) / (mappable_genome_fraction * genome_size);
      double ChIP_value_normalizer = (1/sqrt(2*3.1415*bandwidth*bandwidth));

      double normalized_ChIP_value = cur_ChIP_value / ChIP_value_normalizer;
      double normalized_lambda_basal  = lambda_basal / ChIP_value_normalizer;

      double cur_log_pv = log_Poisson_pv_approximation((int)(normalized_ChIP_value), double_max(normalized_lambda_basal, 1.0));
      */

      //double cur_log_qv = cur_log_pv + log((double)ChIP_calls[k].peaks[l].rank);

      //cout<<"ChIP_p_value: "<<cur_ChIP_value<<" basal: "<<lambda_basal<<" normalized: "<<normalized_ChIP_value<<" log_pv: "<<cur_log_pv<<" log_qv: "<<cur_log_qv<<endl;
      ChIP_calls[k].peaks[l].log_pv = cur_log_pv;            
    }
  }

  cout<<endl<<"calculationg region p-values"<<endl;
  
  for(unsigned int k=0; k<ChIP_calls.size(); k++){
    //double lamda0, lambda1, lambda2, lambda3;

    printf("\r%d of %d", k+1, (int)ChIP_calls.size());
    cout.flush();
    
    int contig_ind = -1;
    for(unsigned int m=0; m<contig_names.size(); m++){
      if(contig_names[m] == ChIP_calls[k].region.chrom){
	contig_ind = m;
      }
    }
    
    assert(contig_ind >= 0);
    
    int cur_ChIP_reads = 0;
    
    int region_start_coord = ChIP_calls[k].region.start_coord;
    int region_end_coord = ChIP_calls[k].region.end_coord;

    int region_center = (region_end_coord + region_start_coord)/2; 

    if(region_end_coord - region_start_coord < region_width){
      region_start_coord = region_center - region_width /2;
      region_end_coord = region_center + region_width / 2;
    }


    double max_background_lambda = -1.0;

    for(unsigned int i=0; i<ChIP_coordinates_for[contig_ind].size(); i++){
      int cur_read_coord = ChIP_coordinates_for[contig_ind][i];
      
      if(cur_read_coord >= region_start_coord && cur_read_coord <= region_end_coord){
	cur_ChIP_reads++;
      }
    }
    
    for(unsigned int i=0; i<ChIP_coordinates_rev[contig_ind].size(); i++){
      int cur_read_coord = ChIP_coordinates_rev[contig_ind][i];
      
      if(cur_read_coord >= region_start_coord && cur_read_coord <= region_end_coord){
	cur_ChIP_reads++;
      }
    }

    ChIP_calls[k].region.ChIP_tags = cur_ChIP_reads;

    //    double ChIP_lambda_basal = ((double)ChIP_tags)  * ((double)(region_end_coord - region_start_coord+1)) / (genome_size*mappable_genome_fraction);
    double ChIP_lambda_basal = ((double)(region_end_coord - region_start_coord+1)) * ChIP_basal_tag_rate;

    for(int l=0; l<=3; l++){
      int tmp_region_start = -100;
      int tmp_region_end = -1000;
      if(l==0){ // the region itself
	tmp_region_start = region_start_coord;
	tmp_region_end = region_end_coord;
      }
      if(l==1){ //the region size given by w1
	tmp_region_start = region_center - w1/2;
	tmp_region_end = region_center + w1/2;
      }
      if(l==2){ //the region size given by w2;
	tmp_region_start = region_center - w2/2;
	tmp_region_end = region_center + w2/2;
      }
      if(l==3){ //the region size given by w3;
	tmp_region_start = region_center - w3/2;
	tmp_region_end = region_center + w3/2;
      }
      
      if(tmp_region_start <0) tmp_region_start = 0;
      if(tmp_region_end > contig_sizes[contig_ind]) tmp_region_end = contig_sizes[contig_ind];
      if(tmp_region_end < 0) tmp_region_end = 0;
      if(tmp_region_end > contig_sizes[contig_ind]) tmp_region_end = contig_sizes[contig_ind];
      
      assert(tmp_region_start <= tmp_region_end);
      
      int cur_background_reads = 0;

      for(unsigned int i=0; i<background_coordinates_for[contig_ind].size(); i++){
	int cur_read_coord = background_coordinates_for[contig_ind][i];
	if(cur_read_coord >= tmp_region_start && cur_read_coord <= tmp_region_end){
	  cur_background_reads++;
	  //if(ChIP_calls[k].region.QuEST_id == "R-1725") cout<<"+"<<cur_read_coord<<endl;
	}
      }
      for(unsigned int i=0; i<background_coordinates_rev[contig_ind].size(); i++){
	int cur_read_coord = background_coordinates_rev[contig_ind][i];
	if(cur_read_coord >= tmp_region_start && cur_read_coord <= tmp_region_end){
	  cur_background_reads++;
	  //if(ChIP_calls[k].region.QuEST_id == "R-1725") cout<<"-"<<cur_read_coord<<endl;
	}
      }

      if(l==0){
	if(background_tags > 0){
	  ChIP_calls[k].region.background_tags = cur_background_reads;
	  //if(ChIP_calls[k].region.QuEST_id == "R-1725"){
	  //  cout<<endl<<ChIP_calls[k].region.QuEST_id<<endl;
	  //  cout<<tmp_region_start<<"-"<<tmp_region_end<<" tags: "<<cur_background_reads<<endl;
	  //  exit(0);
	  //}
	}
      }
      double cur_background_lambda = ((double) cur_background_reads) / ((double)(tmp_region_end - tmp_region_start+1));
      ChIP_calls[k].region.region_size_used = (tmp_region_end - tmp_region_start + 1);
      
      if(cur_background_lambda > max_background_lambda) max_background_lambda = cur_background_lambda;
      

      //cout<<"l: "<<l<<" lambda: "<<cur_background_lambda<<" cur_max: "<<max_background_lambda<<endl;
    }
    //cout<<"max_lambda: "<<max_background_lambda<<endl;

    
    if(background_tags == 0){
      ChIP_calls[k].region.background_tags = 0;
    }
    ChIP_calls[k].region.basal_tags = ChIP_lambda_basal;//__int_max(1, (int)ChIP_lambda_basal);
      
      //}
    
    double lambda_renormalized = double_max(max_background_lambda * ChIP_2_background_tag_ratio, ChIP_lambda_basal);

    if(background_tags == 0){
      lambda_renormalized = ChIP_lambda_basal;
    }
    //double cur_p_value = Poisson_pv_approximation(ChIP_calls[k].region.ChIP_tags, lambda_renormalized);
    double log_cur_p_value = log_Poisson_pv_approximation(ChIP_calls[k].region.ChIP_tags, lambda_renormalized);
    //cout<<"lambda_renormalized: "<<lambda_renormalized<<endl;
    ChIP_calls[k].region.log_pv = log_cur_p_value;
    
    //cout<<"w:"<<region_end_coord-region_start_coord<<" expected: "<<lambda_renormalized<<" observed: "<<ChIP_calls[k].region.ChIP_tags<<" basal: "<<ChIP_lambda_basal<<" log_p_v: "<<log_cur_p_value;;
    //    printf("%f",log_cur_p_value);
    //cout<<endl<<endl;
    //exit(1);
  }

  //assigning ranks
  
  vector<double> region_rank_values;
  vector<double> peak_rank_values;

  for(unsigned int i=0; i<ChIP_calls.size(); i++){
    
    double cur_region_rank_value = ChIP_calls[i].region.log_pv;
    region_rank_values.push_back(cur_region_rank_value);
    ChIP_calls[i].region.rank = -1;
    
    for(unsigned int j=0; j<ChIP_calls[i].peaks.size(); j++){
      double cur_peak_rank_value = ChIP_calls[i].peaks[j].log_pv;
      peak_rank_values.push_back(cur_peak_rank_value);
      
      ChIP_calls[i].peaks[j].rank = -1;
    }
  }

  sort(region_rank_values.begin(), region_rank_values.end());
  sort(peak_rank_values.begin(), peak_rank_values.end());
  
  for(unsigned int i=0; i<ChIP_calls.size(); i++){
    double cur_region_rank_value = ChIP_calls[i].region.log_pv;
    bool rank_assigned = false;
    for(unsigned int k=0; k<region_rank_values.size(); k++){
      if(region_rank_values[k] == cur_region_rank_value){
	ChIP_calls[i].region.rank = k+1;
	rank_assigned = true;
      }
    }
    if(!rank_assigned){
      cerr<<"Warning: Failed to assign rank to: "<<endl;
      cerr<<ChIP_calls[i].region.entry_string<<endl;
      cerr<<"rank value: "<<cur_region_rank_value<<endl;
      cerr<<"ChIP tags: "<<ChIP_calls[i].region.ChIP_tags<<endl;
      //exit(1);
      ChIP_calls[i].region.rank = -1;
    }
  }

  for(unsigned int i=0; i<ChIP_calls.size(); i++){
    for(unsigned int j=0; j<ChIP_calls[i].peaks.size(); j++){
      double cur_peak_rank_value = ChIP_calls[i].peaks[j].log_pv;
      bool rank_assigned  = false;
      for(unsigned int k=0; k<peak_rank_values.size(); k++){

	if(peak_rank_values[k] == cur_peak_rank_value){
	    ChIP_calls[i].peaks[j].rank = k+1;
	    rank_assigned = true;
	}
      }
      if(!rank_assigned){
	cerr<<"Failed to assign the rank to: "<<endl;
	cerr<<ChIP_calls[i].peaks[j].entry_string<<endl;
	cerr<<"rank_value: "<<cur_peak_rank_value<<endl;
	exit(1);
      }
    }
  }

  for(unsigned int i=0; i<ChIP_calls.size(); i++){
    assert(ChIP_calls[i].region.rank != -1);
    for(unsigned int j=0; j<ChIP_calls[i].peaks.size(); j++){
      assert(ChIP_calls[i].peaks[j].rank != -1);
    }
  }

  //doing multiple testing correction (Bonferroni)

  for(unsigned int i=0; i<ChIP_calls.size(); i++){
    double cur_region_log_pv = ChIP_calls[i].region.log_pv;
    
    double cur_region_log_qv = 0;
    if( correction_type == "Bonferroni"){      
      cur_region_log_qv = cur_region_log_pv + log((double)genome_size);
    }
    if( correction_type == "Benjamini-Hochberg" ){
      double min_peak_dist = 200;
      cur_region_log_qv = cur_region_log_pv + log((double) genome_size) - log( min_peak_dist ) - log((double)ChIP_calls[i].region.rank);
    }
    ChIP_calls[i].region.log_qv = double_min(cur_region_log_qv,0);

    //cout<<"cur_region_log_pv: "<<cur_region_log_pv<<" log_qv: "<<cur_region_log_qv<<endl;
    
    for(unsigned int j=0; j<ChIP_calls[i].peaks.size(); j++){
      double cur_peak_log_pv = ChIP_calls[i].peaks[j].log_pv;
      double cur_peak_log_qv = 0; // = cur_peak_log_pv + log((double)genome_size);
      
      if( correction_type == "Bonferroni"){      
	cur_peak_log_qv = cur_peak_log_pv + log((double)genome_size);
      }
      if( correction_type == "Benjamini-Hochberg" ){
	double min_peak_dist = 200;
	cur_peak_log_qv = cur_peak_log_pv + log((double) genome_size) - log( min_peak_dist ) - log((double)ChIP_calls[i].region.rank);
      }

      ChIP_calls[i].peaks[j].log_qv = double_min(cur_peak_log_qv,0);
    }
  }

  cout<<endl<<endl<<"calculating region peak shift: "<<endl;

  for(unsigned int k=0; k<ChIP_calls.size(); k++){
    printf("\r%d of %d", k+1, (int)ChIP_calls.size());
    cout.flush();

    int contig_ind = -1;
    for(unsigned int m=0; m<contig_names.size(); m++){
      if(contig_names[m] == ChIP_calls[k].region.chrom){
	contig_ind = m;
      }
    }
    
    assert(contig_ind >= 0);

    ChIP_calls[k].region.ps = -9999;
    ChIP_calls[k].region.ps_cor = -1.0;

    int region_start_coord = ChIP_calls[k].region.start_coord;
    int region_end_coord = ChIP_calls[k].region.end_coord;
    
    int cur_region_width = region_end_coord - region_start_coord + 1;

    if(cur_region_width < region_width){
      int additional_width = (region_width - cur_region_width)/2;
      region_start_coord = region_start_coord - additional_width;
      region_end_coord = region_end_coord + additional_width;
    }

    int start_coord = region_start_coord - bandwidth*3 - dist_threshold;
    int end_coord = region_end_coord + bandwidth*3 + dist_threshold;

    if(start_coord < 0) start_coord = 0;
    if(end_coord> contig_sizes[contig_ind]) end_coord = contig_sizes[contig_ind];

    
    int region_size = end_coord - start_coord + 1;
    
    vector<int> start_counts_for(region_size);
    vector<int> start_counts_rev(region_size);

    vector<double> profile_for(region_size);
    vector<double> profile_rev(region_size);

    for(unsigned int i=0; i<start_counts_for.size(); i++){
      start_counts_for[i] = 0;
      start_counts_rev[i] = 0;
      
      profile_for[i] = 0.0;
      profile_rev[i] = 0.0;
    }
    
    for(unsigned int i=0; i<ChIP_coordinates_for[contig_ind].size(); i++){
      int cur_read_coord = ChIP_coordinates_for[contig_ind][i];
      
      if(cur_read_coord >= start_coord && cur_read_coord <= end_coord){
	start_counts_for[cur_read_coord - start_coord]++;
      }
    }

    for(unsigned int i=0; i<ChIP_coordinates_rev[contig_ind].size(); i++){
      int cur_read_coord = ChIP_coordinates_rev[contig_ind][i];
      
      if(cur_read_coord >= start_coord && cur_read_coord <= end_coord){
	start_counts_rev[cur_read_coord - start_coord]++;
      }
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

      int center_coord = (profile_for.size()-1)/2 + 1;//(region_size-region_size%2)/2+1;

      int offset = (region_end_coord - region_start_coord) / 2;

      int start_for = center_coord - offset - i; //cur_region_width/2 - i;
      int end_for = center_coord + offset - i; //cur_region_width/2 - i + 1;
      
      int start_rev = center_coord - offset + i; //region_width/2 + i;
      int end_rev = center_coord + offset + i; //region_width/2 + i + 1;

      if(start_for < 0){ start_for=0;} // assert(false); }
      if(end_for > (int)profile_for.size()){  end_for=profile_for.size() -1;}// assert(false); }
      if(start_rev < 0){ start_rev =0;} // assert(false); }
      if(end_rev > (int)profile_for.size()){ end_rev=profile_rev.size() -1;} // assert(false); }
 
      double mu_x=0, var_x=0; //forward strand
      double mu_y=0, var_y=0; //reverse strand

      double mu = 0;
      
      for(int j=start_for; j<end_for; j++){
	mu_x += profile_for[j] / ((double) region_width + 1);
      }
      
      for(int j=start_rev; j<end_rev; j++){
	mu_y += profile_rev[j] / ((double) region_width + 1);
      }

      mu = 0.5*(mu_x+mu_y);
      
      /*
      for(int j=start_for; j<end_for; j++){
	var_x += sqr(profile_for[j]-mu_x) / ((double) region_width + 1);
      }
      
      for(int j=start_rev; j<end_rev; j++){
	var_y += sqr(profile_rev[j]-mu_y) / ((double) region_width + 1);
      }
      */ 

      
      for(int j=start_for; j<end_for; j++){
	var_x += sqr(profile_for[j]-mu) / ((double) region_width + 1);
      }
      
      for(int j=start_rev; j<end_rev; j++){
	var_y += sqr(profile_rev[j]-mu) / ((double) region_width + 1);
      }
      

      double s_x = sqrt(var_x);
      double s_y = sqrt(var_y);
      
      double cur_rho = 0;
      
      /*
      for(int j=0; j<=region_width; j++){
	cur_rho += (profile_for[start_for + j] - mu_x) * 
	  (profile_rev[start_rev+j] - mu_y) / (s_x*s_y*((double)region_width+1.0));
	  }	*/
      
      for(int j=0; j<=region_width; j++){
	cur_rho += (profile_for[start_for + j] - mu) * 
	  (profile_rev[start_rev+j] - mu) / (s_x*s_y*((double)region_width+1.0));
      }	
      
      rho[dist_threshold/2 + i] = cur_rho;
      
    }	  
    
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
    //cout<<"ps: "<<peak_shift<<" cor: "<<rho_max<<endl;

    ChIP_calls[k].region.ps = peak_shift;
    ChIP_calls[k].region.ps_cor = rho_max;
  }

  //===================================================================

  cout<<endl<<endl;
  cout<<"calculating peak call peak shift: "<<endl;

  int cur_peak = 0;
  for(unsigned int k=0; k<ChIP_calls.size(); k++){
    for(unsigned int l=0; l<ChIP_calls[k].peaks.size(); l++){
    
      printf("\r%d of %d", cur_peak, peaks_read);
      cout.flush();
      
      cur_peak++;
      
      int contig_ind = -1;
      for(unsigned int m=0; m<contig_names.size(); m++){
	if(contig_names[m] == ChIP_calls[k].peaks[l].chrom){
	  contig_ind = m;
	}
      }
      
      assert(contig_ind >= 0);
      
      ChIP_calls[k].peaks[l].ps = -9999;
      ChIP_calls[k].peaks[l].ps_cor = -1.0;

      int cur_peak_coord = ChIP_calls[k].peaks[l].coord;

      int region_start_coord = cur_peak_coord - region_width / 2;
      int region_end_coord = cur_peak_coord + region_width / 2;
      
      //int cur_region_width = region_end_coord - region_start_coord + 1;
      
      int start_coord = region_start_coord - bandwidth*3 - dist_threshold;
      int end_coord = region_end_coord + bandwidth*3 + dist_threshold;
      

      if(start_coord < 0) start_coord = 0;
      if(end_coord> contig_sizes[contig_ind]) end_coord = contig_sizes[contig_ind];      
      
      int region_size = end_coord - start_coord + 1;
      
      vector<int> start_counts_for(region_size);
      vector<int> start_counts_rev(region_size);
      
      vector<double> profile_for(region_size);
      vector<double> profile_rev(region_size);
      
      for(unsigned int i=0; i<start_counts_for.size(); i++){
	start_counts_for[i] = 0;
	start_counts_rev[i] = 0;
	
	profile_for[i] = 0.0;
	profile_rev[i] = 0.0;
      }
      
      for(unsigned int i=0; i<ChIP_coordinates_for[contig_ind].size(); i++){
	int cur_read_coord = ChIP_coordinates_for[contig_ind][i];
	
	if(cur_read_coord >= start_coord && cur_read_coord <= end_coord){
	  start_counts_for[cur_read_coord - start_coord]++;
	}
      }
      
      for(unsigned int i=0; i<ChIP_coordinates_rev[contig_ind].size(); i++){
	int cur_read_coord = ChIP_coordinates_rev[contig_ind][i];
	
	if(cur_read_coord >= start_coord && cur_read_coord <= end_coord){
	  start_counts_rev[cur_read_coord - start_coord]++;
	}
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
	
	int center_coord = (profile_for.size()-1)/2 + 1;//(region_size-region_size%2)/2+1;
	
	int offset = region_width/2; //(region_end_coord - region_start_coord) / 2;
	
	int start_for = center_coord - offset - i; //cur_region_width/2 - i;
	int end_for = center_coord + offset - i; //cur_region_width/2 - i + 1;
	
	int start_rev = center_coord - offset + i; //region_width/2 + i;
	int end_rev = center_coord + offset + i; //region_width/2 + i + 1;
	
	if(start_for < 0){ start_for=0;} // assert(false); }
	if(end_for > (int)profile_for.size()){  end_for=profile_for.size() -1;}// assert(false); }
	if(start_rev < 0){ start_rev =0;} // assert(false); }
	if(end_rev > (int)profile_for.size()){ end_rev=profile_rev.size() -1;} // assert(false); }
	
	double mu_x=0, var_x=0; //forward strand
	double mu_y=0, var_y=0; //reverse strand
	
	double mu = 0;

	for(int j=start_for; j<end_for; j++){
	  mu_x += profile_for[j] / ((double) region_width + 1);
	}
	
	for(int j=start_rev; j<end_rev; j++){
	  mu_y += profile_rev[j] / ((double) region_width + 1);
	}

	mu = 0.5*(mu_x + mu_y);
	
	/*
	for(int j=start_for; j<end_for; j++){
	  var_x += sqr(profile_for[j]-mu_x) / ((double) region_width + 1);
	}
	
	for(int j=start_rev; j<end_rev; j++){
	  var_y += sqr(profile_rev[j]-mu_y) / ((double) region_width + 1);
	}
	*/
	
	for(int j=start_for; j<end_for; j++){
	  var_x += sqr(profile_for[j]-mu) / ((double) region_width + 1);
	}
	
	for(int j=start_rev; j<end_rev; j++){
	  var_y += sqr(profile_rev[j]-mu) / ((double) region_width + 1);
	}

	double s_x = sqrt(var_x);
	double s_y = sqrt(var_y);
	
	double cur_rho = 0;
	
	/*
	for(int j=0; j<=region_width; j++){
	  cur_rho += (profile_for[start_for + j] - mu_x) * 
	    (profile_rev[start_rev+j] - mu_y) / (s_x*s_y*((double)region_width+1.0));
	}	
	*/

	for(int j=0; j<=region_width; j++){
	  cur_rho += (profile_for[start_for + j] - mu) * 
	    (profile_rev[start_rev+j] - mu) / (s_x*s_y*((double)region_width+1.0));
	}	
	
	rho[dist_threshold/2 + i] = cur_rho;
	
      }	  
      
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
      //cout<<"ps: "<<peak_shift<<" cor: "<<rho_max<<endl;
      
      ChIP_calls[k].peaks[l].ps = peak_shift;
      ChIP_calls[k].peaks[l].ps_cor = double_min(rho_max,1.0);
    }
  }

  cout<<endl;



  for(unsigned int i=0; i<ChIP_calls.size(); i++){
    int cur_ChIP_tags = ChIP_calls[i].region.ChIP_tags;
    ChIP_tags_within_regions += cur_ChIP_tags;
    int cur_background_tags = ChIP_calls[i].region.background_tags;
    background_tags_within_regions += cur_background_tags;
    
    double tag_ef;

    double expected_ChIP_tags_null = ChIP_calls[i].region.basal_tags;// ChIP_tags * ((double)effective_region_size) / genome_size;
    double expected_null_tag_ef = ((double)cur_ChIP_tags) / expected_ChIP_tags_null;//expected_ChIP_tags_null;

    if(cur_background_tags > 0){
      double observed_tag_ef = ( ((double)cur_ChIP_tags) / ((double)cur_background_tags) ) / ChIP_2_background_tag_ratio;
      tag_ef = double_min(observed_tag_ef, expected_null_tag_ef);

      //tag_ef = double_min( ( ((double)cur_ChIP_tags) / ((double)cur_background_tags) ) / ChIP_2_background_tag_ratio, 
      //			   ChIP_tags * ((double)effective_region_size) / genome_size);			    
      //tag_ef = ( ((double)cur_ChIP_tags) / ((double)cur_background_tags) ) / ChIP_2_background_tag_ratio;
    
    }
    else{
      //tag_ef = ChIP_tags * ((double)effective_region_size) / genome_size;			    
      tag_ef = expected_null_tag_ef;
    }
    ChIP_calls[i].region.tag_ef = tag_ef;
  }

  for(unsigned int i=0; i<ChIP_calls.size(); i++){
    //out_str<<ChIP_calls[i].region.entry_string;
    out_str<<ChIP_calls[i].region.QuEST_id<<" "<<ChIP_calls[i].region.chrom<<" ";
    out_str<<ChIP_calls[i].region.start_coord<<"-"<<ChIP_calls[i].region.end_coord<<" ";
    out_str<<"ChIP: "<<ChIP_calls[i].region.ChIP_value / ChIP_basal_level<<" ";
    if(background_used){
      double cur_ef;
      if(ChIP_calls[i].region.background_value <= 0){
	cur_ef = ChIP_calls[i].region.ChIP_value / ChIP_basal_level;
      }
      else{
	cur_ef = double_min( (ChIP_calls[i].region.ChIP_value / ChIP_basal_level) / (ChIP_calls[i].region.background_value / background_basal_level), ChIP_calls[i].region.ChIP_value / ChIP_basal_level);
      }
      out_str<<"control: "<<ChIP_calls[i].region.background_value / background_basal_level<<" ";
      out_str<<"max_pos: "<<ChIP_calls[i].region.max_pos<<" ef: "<<cur_ef; 
    }
    else{
      double cur_ef = ChIP_calls[i].region.ChIP_value / ChIP_basal_level;
      out_str<<"control: 1.0"<<" ";
      out_str<<"max_pos: "<<ChIP_calls[i].region.max_pos<<" ef: "<<cur_ef;
    }
    
    
    out_str<<" ChIP_tags: "<<ChIP_calls[i].region.ChIP_tags;
    out_str<<" background_tags: "<<ChIP_calls[i].region.background_tags<<" tag_ef: "<<ChIP_calls[i].region.tag_ef;
    out_str<<" ps: "<<ChIP_calls[i].region.ps<<" cor: "<<ChIP_calls[i].region.ps_cor;
    out_str<<" -log10_qv: "<<fabs(-ChIP_calls[i].region.log_qv/log((double)10.0));
    out_str<<" -log10_pv: "<<fabs(-ChIP_calls[i].region.log_pv/log((double)10.0));
    out_str<<" qv_rank: "<<ChIP_calls[i].region.rank;
    out_str<<endl;

    for(unsigned int j=0; j<ChIP_calls[i].peaks.size(); j++){
      //out_str<<ChIP_calls[i].peaks[j].entry_string;
      out_str<<ChIP_calls[i].peaks[j].QuEST_id<<" ";
      out_str<<ChIP_calls[i].peaks[j].chrom<<" ";
      out_str<<ChIP_calls[i].peaks[j].coord<<" ";
      out_str<<"ChIP: "<<ChIP_calls[i].peaks[j].ChIP_value/ChIP_basal_level<<" ";
      if(background_used){
	out_str<<"control: "<<ChIP_calls[i].peaks[j].background_value / background_basal_level<<" ";
      }
      else{
	out_str<<"control: 1.0 ";
      }

      out_str<<"region: "<<ChIP_calls[i].peaks[j].region<<" ";
      double cur_ef;
      if(background_used){
	if(ChIP_calls[i].peaks[j].background_value <= 0){
	  cur_ef = ChIP_calls[i].peaks[j].ChIP_value / ChIP_basal_level;
	}
	cur_ef = double_min( (ChIP_calls[i].peaks[j].ChIP_value/ChIP_basal_level) / (ChIP_calls[i].peaks[j].background_value/background_basal_level), ChIP_calls[i].peaks[j].ChIP_value / ChIP_basal_level);
      }
      else{
	cur_ef = ChIP_calls[i].peaks[j].ChIP_value / ChIP_basal_level;
      }
      out_str<<"ef: "<<cur_ef;

      out_str<<" ps: "<<ChIP_calls[i].peaks[j].ps<<" cor: "<<ChIP_calls[i].peaks[j].ps_cor;
      out_str<<" -log10_qv: "<<fabs(-ChIP_calls[i].peaks[j].log_qv/log((double)10.0));
      out_str<<" -log10_pv: "<<fabs(-ChIP_calls[i].peaks[j].log_pv/log((double)10.0));
      out_str<<" qv_rank: "<<ChIP_calls[i].peaks[j].rank<<endl;
    }
    out_str<<endl;
  }
  
  report_str<<"ChIP tags within regions: "<<ChIP_tags_within_regions<<" ( "<<100*((double)ChIP_tags_within_regions) / ((double)ChIP_tags)<<" % )"<<endl;
  report_str<<"background tags within regions: "<<background_tags_within_regions;
  if(background_tags > 0){
    report_str<<" ( "<<100*((double)background_tags_within_regions) / ((double)background_tags)<<" % )"<<endl;
  }


  out_str.close();
  report_str.close();

  cout<<"ChIP tags within regions: "<<ChIP_tags_within_regions<<" ( "<<100*((double)ChIP_tags_within_regions) / ((double)ChIP_tags)<<" % )"<<endl;
  cout<<"background tags within regions: "<<background_tags_within_regions<<" ( "<<100*((double)background_tags_within_regions) / ((double)background_tags)<<" % )"<<endl;

  return 0;

}
