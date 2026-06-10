// generate_and_verify_tifr_gas.cc
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ViewMedium.hh"
#include <TApplication.h>
#include <TCanvas.h>
#include <iostream>

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  MediumMagboltz gas;

  // ✅ Set composition explicitly — verify percentages sum to 100
  gas.SetComposition("c2h2f4", 95.5,
                     "isobutane", 4.2,
                     "sf6", 0.3);

  // Verify composition was set correctly
  std::cout << "Number of components: "
            << gas.GetNumberOfComponents() << std::endl;

  for (int i = 0; i < gas.GetNumberOfComponents(); i++) {
    std::string name;
    double frac;
    gas.GetComponent(i, name, frac);
    std::cout << "Component " << i << ": "
              << name << " = " << frac*100 << "%" << std::endl;
  }

  gas.SetTemperature(293.15);
  gas.SetPressure(760.);

  // Field range: 1 kV/cm to 150 kV/cm
  // Your operating point: 60 kV/cm — well within range
  gas.SetFieldGrid(1000., 150000., 20, true);

  std::cout << "Running Magboltz..." << std::endl;
  gas.GenerateGasTable(20, true);

  gas.WriteGasFile("tifr_verified.gas");
  std::cout << "Done." << std::endl;

  // ── Immediately verify at operating point ──────────
  double vx, vy, vz, alpha, eta;
  const double E = 60000.;  // V/cm

  gas.ElectronVelocity  (0., -E, 0., 0., 0., 0., vx, vy, vz);
  gas.ElectronTownsend  (0., -E, 0., 0., 0., 0., alpha);
  gas.ElectronAttachment(0., -E, 0., 0., 0., 0., eta);

  std::cout << "\n=== AT 60 kV/cm ===" << std::endl;
  std::cout << "Drift velocity : " << std::abs(vy)*1e4 << " um/ns" << std::endl;
  std::cout << "Alpha          : " << alpha << " /cm" << std::endl;
  std::cout << "Eta            : " << eta   << " /cm" << std::endl;
  std::cout << "Gain (2mm)     : "
            << std::exp((alpha-eta)*0.2) << std::endl;

  // ── Plot transport coefficients ────────────────────
  TCanvas* cv = new TCanvas("cv", "Drift Velocity", 800, 600);
  ViewMedium view;
  view.SetMedium(&gas);
  view.SetCanvas(cv);
  view.PlotElectronVelocity('e');
  cv->SaveAs("tifr_velocity_verified.png");

  TCanvas* ca = new TCanvas("ca", "Townsend", 800, 600);
  view.SetCanvas(ca);
  view.PlotElectronTownsend('e');
  ca->SaveAs("tifr_townsend.png");

  TCanvas* cat = new TCanvas("cat", "Attachment", 800, 600);
  view.SetCanvas(cat);
  view.PlotElectronAttachment('e');
  cat->SaveAs("tifr_attachment.png");

  app.Run(true);
  return 0;
}
