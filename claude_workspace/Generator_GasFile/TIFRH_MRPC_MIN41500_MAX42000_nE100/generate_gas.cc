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

  // -----------------------------
  // 2. Electric field grid
  // -----------------------------
  const size_t nE = 100;
  const double emin = 41500.;       // V/cm
  const double emax = 42000.;    // V/cm
  constexpr bool useLog = false;

  gas.SetFieldGrid(emin, emax, nE, useLog);

  // -----------------------------
  // 3. Collision statistics
  // -----------------------------
  const int ncoll = 50;   // MUCH better than 10

  gas.GenerateGasTable(ncoll);

  // -----------------------------
  // 4. Save gas table
  // -----------------------------
  gas.WriteGasFile("TIFRH_41500_42000.gas");

  return 0;
}
