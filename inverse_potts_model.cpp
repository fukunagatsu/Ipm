#include "inverse_potts_model.h"
#include "fmath.hpp"
#include <getopt.h>
#include <stdlib.h>

void InversePottsModel::Run(){
  ReadData();  
  PersistentContrastiveDivergence();
  Output();
}

void InversePottsModel::Output(){  
  ofstream ofs(_output_file_name.c_str());

  ofs << "OG_index_1 OG_index_2 DI" << endl;
  for(int i = 0; i < _item_size; i++){
    for(int j = i+1; j < _item_size; j++){
      
      double temp_sum = 0.0;
      vector<vector<double> > direct_prob(_state_size, vector<double>(_state_size, 0.0));
      
      for(int k = 0; k < _state_size; k++){
	for(int l = 0; l < _state_size; l++){
	  double temp = 0;
	  temp += (k == 0 || l == 0) ? 0 : _j_parameter[GetJID(i,j,k-1,l-1)];
	  temp += (k == 0) ? 0 : _h_parameter[GetHID(i,k-1)];
	  temp += (l == 0) ? 0 : _h_parameter[GetHID(j,l-1)];
	  direct_prob[k][l] = exp(temp);
	  temp_sum += direct_prob[k][l];
	}
      }
      
      vector<double> fi(_state_size, 0.0); vector<double> fj(_state_size, 0.0);
      for(int k = 0; k < _state_size; k++){
	for(int l = 0; l < _state_size; l++){
	  direct_prob[k][l] /= temp_sum;
	  fi[k] += direct_prob[k][l];
	  fj[l] += direct_prob[k][l];
	}
      }
      
      double di = 0.0;
      for(int k = 0; k < _state_size; k++){
	for(int l = 0; l < _state_size; l++){
	  di += direct_prob[k][l] * log(direct_prob[k][l]/(fi[k] * fj[l]));
	}
      }
      ofs << i << " " << j << " " << di << endl;   
    }
  }
  
  ofs.close();
}



void InversePottsModel::DataSampling(vector<int>& data, default_random_engine& engine){
  uniform_real_distribution<> dist(0.0, 1.0);
  
  for(int i = 0; i < _item_size; i++){
    vector<double> hamiltonian(_state_size, 0.0);
   
    for(int j = 1; j < _state_size; j++){      
      hamiltonian[j] = _h_parameter[GetHID(i,j-1)];
	
      for(int k = 0; k < _item_size; k++){
	if(i != k && data[k] != 0){
	  if(i < k){
	    hamiltonian[j] += _j_parameter[GetJID(i,k,j-1,data[k]-1)];
	  }else{
	    hamiltonian[j] += _j_parameter[GetJID(k,i,data[k]-1,j-1)];
	  }
	}
      }
    }
   
    vector<double> prob_vector(_state_size, 1.0);
    double prob = dist(engine);
    bool flag = true;
    for(int j = 0; j < _state_size-1; j++){      
      double temp_sum = 1.0;      
      for(int k = 0; k < _state_size; k++){
	if(k != j){
	  temp_sum += fmath::expd(hamiltonian[k] - hamiltonian[j]);
	}	  
      }
      prob_vector[j] = 1.0 / temp_sum;
      if(j != 0){
	prob_vector[j] += prob_vector[j-1];
      }
      if(prob < prob_vector[j]){
	data[i] = j;
	flag = false;
	break;
      }      
    }
    if(flag){
      data[i] = _state_size-1;
    }
  }
}

void InversePottsModel::DataSamplingQ2(vector<int>& data, default_random_engine& engine){
  for(int i = 0; i < _item_size; i++){
    double hamiltonian = _h_parameter[i];
    for(int j = 0; j <  _item_size; j++){
      if(i != j && data[j] == 1){
	if(i < j){
	  hamiltonian += _j_parameter[GetJID(i,j,0,0)];
	}else{
	  hamiltonian += _j_parameter[GetJID(j,i,0,0)];
	}
      }
    }

    double prob = 1.0 / (1 + fmath::expd(hamiltonian));
    bernoulli_distribution dist(prob);
    data[i] = dist(engine) ? 0 : 1;
  }
}

void InversePottsModel::PersistentContrastiveDivergence(){

  InitiallizeParameters();
  
  random_device seed_gen;
  default_random_engine engine(seed_gen());

  vector<vector<int> > sampled_dataset(_sample_size, vector<int>(_item_size, 0));
  uniform_int_distribution<> ui_dist(0, _data_size-1);
  for(int i = 0; i < _sample_size; i++){    
    int temp_id = ui_dist(engine);
    for(int j = 0; j < _item_size; j++){
      sampled_dataset[i][j] = _dataset[temp_id][j];
    }
  }

  for(int i = 0; i < _loop_size; i++){
    vector<double> sampled_data_single_freq(_item_size*(_state_size-1), 0.0);
    vector<double> sampled_data_pair_freq(_item_size*_item_size*(_state_size-1)*(_state_size-1), 0.0);
    
    for(int j = 0; j < _sample_size; j++){
      if(_state_size == 2){
	DataSamplingQ2(sampled_dataset[j], engine);
      }else{
	DataSampling(sampled_dataset[j], engine);
      }
    }

    for(int j = 0; j < _sample_size; j++){
      for(int k = 0; k < _item_size; k++){
	int k_state = sampled_dataset[j][k];
	if(k_state != 0){
	  sampled_data_single_freq[GetHID(k,k_state-1)] += 1.0/_sample_size;
	  
	  for(int l = k+1;  l < _item_size; l++){
	    int l_state = sampled_dataset[j][l];
	    if(l_state != 0){
	      sampled_data_pair_freq[GetJID(k,l,k_state-1,l_state-1)] += 1.0/_sample_size;
	    }
	  }
	}
      }    
    }
    double h_sum = 0.0; double j_sum = 0.0;
    for(int j = 0; j < _item_size; j++){
      for(int k = 0; k < _state_size-1; k++){
	int h_id = GetHID(j,k);
	double temp = _single_freq[h_id] - sampled_data_single_freq[h_id];	
	
	temp -= _lambda * 2 * _h_parameter[h_id];	
	h_sum += abs(temp);
	_h_parameter[h_id] += _epsilon * temp;
      }
	
      
      for(int k = j+1;  k < _item_size; k++){
	for(int l = 0; l < _state_size-1; l++){
	  for(int m = 0; m < _state_size-1; m++){
	    int j_id = GetJID(j,k,l,m);
	    double temp =  _pair_freq[j_id] - sampled_data_pair_freq[j_id];
	    
	    temp -= _lambda * 2 * _j_parameter[j_id];
	    j_sum += abs(temp);
	    _j_parameter[j_id] += _epsilon * temp;
	  }
	}
      }
    }
    
  }

  _single_freq.shrink_to_fit();
  _pair_freq.shrink_to_fit();
  
}

void InversePottsModel::ReadData(){  
  ifstream ifs(_input_file_name.c_str());
 
  if (!ifs) {
    cout << "Cannot open " + _input_file_name << endl;
    exit(1);
  }
  string buf;
  ifs >> buf; ifs >> buf; ifs >> buf;
  ifs >> _data_size; ifs >> _item_size; ifs >> _state_size;
  _state_size += 1;

  _single_freq.resize(_item_size*(_state_size-1), 0.0);
  _pair_freq.resize(_item_size*_item_size*(_state_size-1)*(_state_size-1),0.0);

  _h_parameter.resize(_item_size*(_state_size-1), 0.0);
  _j_parameter.resize(_item_size*_item_size*(_state_size-1)*(_state_size-1),0.0);

  _dataset.resize(_data_size, vector<int>(_item_size, 0));
  
  for(int i = 0; i < _data_size; i++){
    for(int j = 0; j < _item_size; j++){
      ifs >> _dataset[i][j];
    }

    for(int j = 0; j < _item_size; j++){
      int j_state = _dataset[i][j];
      if(j_state != 0){
	_single_freq[GetHID(j,j_state-1)] += 1.0/_data_size;
	for(int k = j+1;  k < _item_size; k++){
	  int k_state = _dataset[i][k];
	  if(k_state != 0){
	    _pair_freq[GetJID(j,k,j_state-1,k_state-1)] += 1.0/_data_size;
	  }
	}
      }
    }    
  }
  ifs.close();
}

void InversePottsModel::InitiallizeParameters(){
  random_device seed_gen;
  default_random_engine engine(seed_gen());

  normal_distribution<> dist(_init_mu, _init_sigma);
  
  for(int i = 0; i < _item_size; i++){
    for(int j= 0; j < _state_size-1; j++){
      _h_parameter[GetHID(i,j)] = dist(engine);
      
    }
  }
 
  for(int i = 0; i < _item_size; i++){
    for(int j = i+1; j < _item_size; j++){
      for(int k = 0; k < _state_size-1; k++){
	for(int l = 0; l < _state_size-1; l++){
	  _j_parameter[GetJID(i,j,k,l)] = dist(engine);
	}
      }
    }
  }
  
}

double InversePottsModel::GetHID(int i, int j){
  return i*(_state_size - 1)+j;
}
double InversePottsModel::GetJID(int i, int j, int k, int l){
  return i*_item_size*(_state_size-1)*(_state_size-1) + j*(_state_size - 1)*(_state_size - 1) + k*(_state_size - 1) + l;
}

void InversePottsModel::SetParameters(int argc,char* argv[]){
  int c;
  extern char *optarg;
  
  while ((c = getopt(argc, argv, "i:o:e:l:c:s:")) != -1) {
    switch (c) {
    case 'i':
      _input_file_name = optarg;
      break;

    case 'o':
      _output_file_name = optarg;
      break;

    case 'e':
      _epsilon = atof(optarg);  
      break;

    case 'l':
      _lambda = atof(optarg);  
      break;

    case 'c':
      _loop_size = atoi(optarg);  
      break;

    case 's':
      _sample_size = atoi(optarg);  
      break;

    default:
      cerr << "The argument is an invalid command." << endl;
      exit(1); 
    }      
  }
  
  
  if(_input_file_name == "" || _output_file_name == ""){
    cerr << "Please specify your input file name and output file name." << endl;
    exit(1);
  }
}
