#include "Garfield/FundamentalConstants.hh"
#include "Garfield/MediumMagboltz.hh"

using namespace Garfield;

int main() {

  // -----------------------------
  // 1. Gas conditions
  // -----------------------------
  const double temperature = 293.15;                // K
  const double pressure = AtmosphericPressure;      // 760 Torr

  MediumMagboltz gas("C2H2F4", 95.5, "iC4H10", 4.2, "SF6", 0.3);
  gas.SetTemperature(temperature);
  gas.SetPressure(pressure);

 
   const int ncoll = 10;  

   gas.SetMaxElectronEnergy(200.);
   gas.SetNumberOfCollisions(10);

  // -----------------------------
  // Save gas table
  // -----------------------------
  gas.Initialise(true);

  gas.WriteGasFile("debug_workhorsemixture.gas");

  return 0;
}
