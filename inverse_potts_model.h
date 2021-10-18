#ifndef INVERSE_POTTS_MODEL_H
#define INVERSE_POTTS_MODEL_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cfloat>
#include <cmath>

using namespace std;

class InversePottsModel{
 public:
  InversePottsModel(){
    _input_file_name = "";
    _output_file_name = "";
    _data_size = 0;
    _item_size = 0;
    _state_size = 0;
    
    _lambda = 0.0;
    _init_mu = 0.0;
    _init_sigma = 0.0;
    _epsilon = 0.01;
    _loop_size = 3000;
    _sample_size = 200;
  }
  void SetParameters(int argc, char* argv[]);
  void Run();  

 private:
  int _data_size;
  int _item_size;
  int _state_size;
  int _sample_size;
  int _loop_size;
  double _lambda;
  double _init_mu;
  double _init_sigma;
  double _epsilon;
  string _input_file_name;
  string _output_file_name;
  vector<double > _single_freq;
  vector<double > _h_parameter;
  vector<double > _pair_freq;
  vector<double > _j_parameter;
  vector<vector<int> > _dataset;
  vector<int> _pre_calc;

  void ReadData();
  void PersistentContrastiveDivergence();
  void InitiallizeParameters();  
  void DataSampling(vector<int>& data, default_random_engine& engine);
  void DataSamplingQ2(vector<int>& data, default_random_engine& engine);
  void Output();
  
  double GetHID(int i, int j);
  double GetJID(int i, int j, int k, int l);
};

#endif
