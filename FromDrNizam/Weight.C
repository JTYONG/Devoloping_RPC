#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/ViewSignal.hh"

#define LOG(x) std::cout << x << std::endl

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  // ──────────────────────────────────────────────────
  // 1. GEOMETRY
  // ──────────────────────────────────────────────────
  const double dMylar = 0.010;   // [cm]
  const double dGlass = 0.300;   // [cm]
  const double dGas   = 0.200;   // [cm]

  const double epsMylar = 3.1;
  const double epsGlass = 7.0;
  const double epsGas   = 1.0;

  std::vector<double> eps       = {epsMylar, epsGlass, epsGas, epsGlass, epsMylar};
  std::vector<double> thickness = {dMylar,   dGlass,   dGas,   dGlass,   dMylar};

  const double dTotal   = std::accumulate(thickness.begin(), thickness.end(), 0.);
  const double gasStart = dMylar + dGlass;        // anode  (lower y) = 0.310
  const double gasEnd   = gasStart + dGas;        // cathode (upper y) = 0.510
  const double margin   = 1.5e-2;                 // 1.5mm injection offset

  const double voltage  = -12000.;                // [V]

  ComponentParallelPlate rpc;
  rpc.Setup(eps.size(), eps, thickness, voltage);

  const std::string electrode = "Readout";
  rpc.AddPlane(electrode, true);

  LOG("Geometry:");
  LOG("  gasStart (anode)  = " << gasStart << " cm");
  LOG("  gasEnd (cathode)  = " << gasEnd   << " cm");
  LOG("  dGas              = " << dGas     << " cm");
  LOG("  dTotal            = " << dTotal   << " cm");
  
  std::cout << "\nWeighting potential scan\n";

for (double y = gasStart; y <= gasEnd; y += 0.02) {

  double wp =
      rpc.WeightingPotential(
          0., y, 0.,
          electrode);

  std::cout
      << y
      << "  "
      << wp
      << std::endl;
}

  // ──────────────────────────────────────────────────
  // 2. GAS
  // ──────────────────────────────────────────────────
  MediumMagboltz gas;
  gas.LoadGasFile("c2h2f4_95-5_ic4h10_4-2_sf6_0-3.gas");

  gas.LoadIonMobility("IonMobility_C8Hn+_iC4H10.txt");

  rpc.SetMedium(&gas);   // called once only

  // ──────────────────────────────────────────────────
  // 3. GEOMETRY BOX — GAS GAP ONLY
  // ──────────────────────────────────────────────────
  // Box covers ONLY the gas gap, not glass or mylar
  // SolidBox(cx, cy, cz, half-lx, half-ly, half-lz)
  const double boxCY  = (gasStart + gasEnd) / 2.;  // centre of gas gap
  const double boxHLY = dGas / 2.;                 // half-length in y

  GeometrySimple geo;
  SolidBox box(0., boxCY, 0.,
               5., boxHLY, 5.);
  geo.AddSolid(&box, &gas);
  rpc.SetGeometry(&geo);

  // Verify box boundaries
  LOG("\nGas box verification:");
  LOG("  Centre y  = " << boxCY           << " cm");
  LOG("  Half-len  = " << boxHLY          << " cm");
  LOG("  y_min     = " << boxCY - boxHLY  << " cm  (should be " << gasStart << ")");
  LOG("  y_max     = " << boxCY + boxHLY  << " cm  (should be " << gasEnd   << ")");
  LOG("  Match     : " << (std::abs(boxCY - boxHLY - gasStart) < 1e-10 &&
                           std::abs(boxCY + boxHLY - gasEnd)   < 1e-10
                           ? " CORRECT" : " MISMATCH"));

  // ──────────────────────────────────────────────────
  // 4. DIAGNOSTIC
  // ──────────────────────────────────────────────────
  double ex, ey, ez, pot;
  Medium* med = nullptr;
  int status  = 0;

  rpc.ElectricField(0., 0.5*(gasStart+gasEnd), 0.,
                    ex, ey, ez, pot, med, status);

  double vx, vy, vz, alpha, eta;
  gas.ElectronVelocity  (ex, ey, ez, 0.,0.,0., vx, vy, vz);
  gas.ElectronTownsend  (ex, ey, ez, 0.,0.,0., alpha);
  gas.ElectronAttachment(ex, ey, ez, 0.,0.,0., eta);

  double ivx, ivy, ivz;
  gas.IonVelocity(ex, ey, ez, 0.,0.,0., ivx, ivy, ivz);

  const double vDrift   = std::abs(vy);
  const double tTransit = dGas / vDrift;
  const double vIon     = std::abs(ivy);
  const double tIon     = (vIon > 0.) ? dGas / vIon : tTransit * 200.;

  LOG("\n=== GAS PROPERTIES ===");
  LOG("Ey              = " << ey          << " V/cm");
  LOG("vy (electron)   = " << vy          << " cm/ns");
  LOG("ivy (ion)       = " << ivy         << " cm/ns");
  LOG("Alpha           = " << alpha       << " /cm");
  LOG("Eta             = " << eta         << " /cm");
  LOG("Electron transit= " << tTransit    << " ns");
  LOG("Ion transit     = " << tIon        << " ns");

  // ──────────────────────────────────────────────────
  // 5. SENSOR
  // Time window covers full ion drift
  // ──────────────────────────────────────────────────
  Sensor sensor(&rpc);
  sensor.SetArea(-1., gasStart, -1., 1., gasEnd, 1.);
  sensor.AddElectrode(&rpc, electrode);

  const double       tWindow = tIon * 1.5;
  const unsigned int nBins   = 10000;
  const double       tMin    = 0.;
  const double       tMax    = tWindow;
  const double       tStep   = (tMax - tMin) / nBins;

  sensor.SetTimeWindow(tMin, tStep, nBins);

  LOG("\nTime window:");
  LOG("  tMax  = " << tMax  << " ns");
  LOG("  tStep = " << tStep << " ns/bin");
  LOG("  nBins = " << nBins);

  // ──────────────────────────────────────────────────
  // 6. AVALANCHE — electrons (microscopic)
  // ──────────────────────────────────────────────────
  AvalancheMicroscopic aval(&sensor);
  aval.EnableSignalCalculation();
  aval.UseWeightingPotential();

  // ──────────────────────────────────────────────────
  // 7. ION DRIFT (MC)
  // ──────────────────────────────────────────────────
  AvalancheMC ionDrift(&sensor);
  ionDrift.EnableSignalCalculation();
  ionDrift.SetDistanceSteps(1.e-3);   // 20 µm steps

  // ──────────────────────────────────────────────────
  // 8. INJECT SINGLE ELECTRON
  // ──────────────────────────────────────────────────
  sensor.ClearSignal();

  const double y0 = gasEnd - margin;   // 1.5mm inside cathode
  LOG("\nInjecting electron at y=" << y0
      << "  (cathode=" << gasEnd
      << "  anode="    << gasStart << ")");

  aval.AvalancheElectron(0., y0, 0., 0., 0.1);

  int ne = 0, ni = 0;
  aval.GetAvalancheSize(ne, ni);
  LOG("\nAvalanche: electrons=" << ne << "  ions=" << ni);

 /* // ──────────────────────────────────────────────────
  // 9. DRIFT IONS — one per avalanche electron, birth point only
  // ──────────────────────────────────────────────────
  int nIonsDrifted = 0;
  int nIonsSkipped = 0;
  const double gapMargin = 1.e-4;   // 1 µm boundary guard

  for (const auto& e : aval.GetElectrons()) {
    if (e.path.empty()) continue;

    // Use birth point (path.front()) — one ion left here per electron
    const auto& birth = e.path.front();

    // Only drift if born inside gas gap
    if (birth.y <= gasStart + gapMargin ||
        birth.y >= gasEnd   - gapMargin) {
      nIonsSkipped++;
      continue;
    }

    ionDrift.DriftIon(birth.x, birth.y, birth.z, birth.t);
    nIonsDrifted++;
  }

  LOG("Ions drifted : " << nIonsDrifted);
  LOG("Ions skipped : " << nIonsSkipped);*/
  
// ─────────────────────────────────────────────────
// ION DRIFT — all ionisation points from avalanche
// ──────────────────────────────────────────────────
int nIonsDrifted = 0;
int nIonsSkipped = 0;
const double gapMargin = 1.e-4;

const auto& ions = aval.GetIons();
LOG("Total ionisation points : " << ions.size());
// Should match ni from GetAvalancheSize()

for (const auto& ion : ions) {
  if (ion.path.empty()) continue;

  const auto& birth = ion.path.front();

  if (birth.y <= gasStart + gapMargin ||
      birth.y >= gasEnd   - gapMargin) {
    nIonsSkipped++;
    continue;
  }

  ionDrift.DriftIon(birth.x, birth.y, birth.z, birth.t);
  nIonsDrifted++;
}

LOG("Ions drifted : " << nIonsDrifted);
LOG("Ions skipped : " << nIonsSkipped);

// Print ion birth position histogram
int posBins[10] = {0};
for (const auto& ion : ions) {
  if (ion.path.empty()) continue;
  const auto& birth = ion.path.front();
  if (birth.y < gasStart || birth.y > gasEnd) continue;
  int bin = int((birth.y - gasStart) / dGas * 10);
  if (bin >= 0 && bin < 10) posBins[bin]++;
}

LOG("Ion birth distribution (0=anode, 9=cathode):");
for (int i = 0; i < 10; i++) {
  LOG("  bin " << i << " [y="
      << gasStart + i*dGas/10. << "-"
      << gasStart + (i+1)*dGas/10. << "]: "
      << posBins[i] << " ions");
}

  // ──────────────────────────────────────────────────
  // 10. SIGNAL SUMMARY
  // ──────────────────────────────────────────────────

  // Read BEFORE IntegrateSignal
  const double totalCharge = sensor.GetTotalInducedCharge(electrode);

  int    nNonZero  = 0;
  double maxSignal = -1e99;
  double minSignal = +1e99;
  double tPeak     = 0.;

  for (unsigned int i = 0; i < nBins; i++) {
    double s = sensor.GetSignal(electrode, i);
    if (std::abs(s) > 1.e-20) nNonZero++;
    if (s > maxSignal) { maxSignal = s; tPeak = tMin + i*tStep; }
    if (s < minSignal)   minSignal = s;
  }

  LOG("\n=== SIGNAL SUMMARY ===");
  LOG("Total charge   : " << totalCharge << " fC");
  LOG("Non-zero bins  : " << nNonZero   << " / " << nBins);
  LOG("Max signal     : " << maxSignal  << " fC/ns  at t=" << tPeak << " ns");
  LOG("Min signal     : " << minSignal  << " fC/ns");
  LOG("Signal duration: " << nNonZero * tStep << " ns");
  
  // Manual charge integration 
double manualCharge = 0.;
double electronCharge = 0.;
double ionCharge = 0.;

for (unsigned int i = 0; i < nBins; i++) {
    double s = sensor.GetSignal(electrode, i);
    manualCharge += s * tStep;
}

// Also get electron and ion contributions separately
for (unsigned int i = 0; i < nBins; i++) {
    double se = sensor.GetElectronSignal(electrode, i);
    double si = sensor.GetIonSignal(electrode, i);
    electronCharge += se * tStep;
    ionCharge      += si * tStep;
}

LOG("=== CHARGE DIAGNOSIS ===");
LOG("Manual total charge   : " << manualCharge    << " fC");
LOG("Manual electron charge: " << electronCharge  << " fC");
LOG("Manual ion charge     : " << ionCharge       << " fC");
LOG("GetTotalInducedCharge : " << sensor.GetTotalInducedCharge(electrode) << " fC");
LOG("Non-zero bins         : " << nNonZero << " / " << nBins);
LOG("Max signal            : " << maxSignal << " fC/ns");
LOG("Min signal            : " << minSignal << " fC/ns");
  

  // ──────────────────────────────────────────────────
  // 11. PLOT
  // ──────────────────────────────────────────────────
  sensor.ExportSignal(electrode, "ElectronIonCurrent");

  TCanvas* c1 = new TCanvas("c1", "Electron + Ion Signal", 900, 700);
  ViewSignal sigView(&sensor);
  sigView.SetCanvas(c1);
  sigView.PlotSignal(electrode);
  c1->SaveAs("rpc_electron_ion_signal.png");
  c1->Update();
  gSystem->ProcessEvents();

  //  Integrate AFTER exporting current
  sensor.IntegrateSignal(electrode);
  sensor.ExportSignal(electrode, "ElectronIonCharge");

  TCanvas* c2 = new TCanvas("c2", "Induced Charge", 900, 700);
  ViewSignal chargeView(&sensor);
  chargeView.SetCanvas(c2);
  chargeView.PlotSignal(electrode);
  c2->SaveAs("rpc_electron_ion_charge.png");
  c2->Update();
  gSystem->ProcessEvents();

  LOG("\nDone.");
  app.Run(true);
  return 0;
}
