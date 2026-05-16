#include "Garfield/FundamentalConstants.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumGas.hh"

using namespace Garfield;

int main() {

  MediumMagboltz gas;
  gas.LoadGasFile("TIFRH_40000_40500.gas");
  gas.MergeGasFile("TIFRH_40500_41000.gas",false);
  gas.MergeGasFile("TIFRH_41000_41500.gas",false);
  gas.MergeGasFile("TIFRH_41500_42000.gas",false);
  gas.MergeGasFile("TIFRH_42000_42500.gas",false);
  gas.MergeGasFile("TIFRH_42500_43000.gas",false);
  
  gas.WriteGasFile("TIFRH_merged_40000_43000.gas");

  return 0;
}