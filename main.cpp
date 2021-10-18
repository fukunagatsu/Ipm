#include "inverse_potts_model.h"

int main(int argc, char* argv[]){
  InversePottsModel inverse_potts_model;

  inverse_potts_model.SetParameters(argc,argv);
  inverse_potts_model.Run();


}
