// generate_gas.cpp
// Garfield++ gas table generator for RPC-style C2H2F4/iC4H10/SF6 mixture.
//
// E-field grid: 11 linear points, 40000 to 45000 V/cm (step 100 V/cm)
// Composition: 95.5% C2H2F4, 4.2% iC4H10, 0.3% SF6
// Conditions: T = 296.15 K, P = 760 Torr
//
// Usage:
//   ./generate_gas              # default ncoll=10 (production)
//   ./generate_gas 2            # quick sanity test
//   ./generate_gas 5            # rough table
//   ./generate_gas 10           # production
//   ./generate_gas 20           # publication-grade
//
// Output: rpc_95.5_4.2_0.3.gas (in current directory)

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <cstdlib>

#include "Garfield/MediumMagboltz.hh"

int main(int argc, char* argv[]) {
  // ----- Parse ncoll from command line (default = 10) -----
  int ncoll = 20;
  if (argc > 1) {
    ncoll = std::atoi(argv[1]);
    if (ncoll < 1) {
      std::cerr << "ERROR: ncoll must be >= 1 (got " << ncoll << ")\n";
      return 1;
    }
  }

  using clk = std::chrono::steady_clock;
  auto t0 = clk::now();

  std::cout << "==========================================\n";
  std::cout << "  RPC Gas Table Generator (Garfield++)\n";
  std::cout << "==========================================\n";
  std::cout << "  Mixture:   95.5% C2H2F4 / 4.2% iC4H10 / 0.3% SF6\n";
  std::cout << "  T, P:      296.15 K, 760 Torr\n";
  std::cout << "  E-field:   38000 to 39000 V/cm, step 100 V/cm (11 points)\n";
  std::cout << "  ncoll:     " << ncoll << "\n";
  std::cout << "==========================================\n\n";

  // ----- Set up the medium -----
  Garfield::MediumMagboltz gas;
  gas.SetComposition("C2H2F4", 95.5, "iC4H10", 4.2, "SF6", 0.3);
  gas.SetTemperature(296.15);
  gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
  gas.EnableThermalMotion(true);

  // ----- Build the E-field grid: 11 linear points -----
  std::vector<double> eFields;
  eFields.reserve(11);
  for (int i = 0; i <= 50; ++i) {
    eFields.push_back(40000.0 + 100.0 * i);
  }
	
  //for (const auto& item : eFields){
	//  std::cout<<item<<" ";
  //}

  // For E-only tables we use a single B-field point and a single angle.
  std::vector<double> bFields = {0.};
  std::vector<double> angles  = {1.5707963267948966};  // pi/2

  gas.SetFieldGrid(eFields, bFields, angles);

  std::cout << "E-field grid (V/cm): ";
  for (auto e : eFields) std::cout << static_cast<int>(e) << " ";
  std::cout << "\n\n";

  // ----- Run Magboltz -----
  std::cout << "Running Magboltz with ncoll = " << ncoll << " ...\n";
  //std::cout << "Estimated time: " << (ncoll * 11 * 0.5) << " minutes\n";
  //std::cout << "(estimate is rough; SF6-rich mixes can be slower)\n" << std::flush;

  gas.GenerateGasTable(ncoll, /*verbose=*/true);

  // ----- Write the gas file -----
  const std::string outName = "../gas_files/rpc_95.5_4.2_0.3_40-45_ncoll=20.gas";
  gas.WriteGasFile(outName);

  auto t1 = clk::now();
  auto secs = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();

  std::cout << "\n==========================================\n";
  std::cout << "  DONE\n";
  std::cout << "  Output: " << outName << "\n";
  std::cout << "  Time:   " << secs << " s ("
            << std::fixed << std::setprecision(1) << (secs/60.0) << " min)\n";
  std::cout << "==========================================\n";

  return 0;
}
