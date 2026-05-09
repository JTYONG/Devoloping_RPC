#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
#include "TString.h"

// ROOT includes
#include <TFile.h>
#include <TH2D.h>
#include <TGeoManager.h>
#include <TTree.h>

#include <ctime>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
// Standard includes

#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <array>


//#include "Garfield/AvalancheGrid.hh"
//#include "Garfield/AvalancheMicroscopic.hh"
//#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"

// Setup geometry and material
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/GeometrySimple.hh"
#include <Garfield/SolidBox.hh>
#include <Garfield/Solid.hh>

// Medium Classes
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Medium.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/MediumPlastic.hh"
#include "Garfield/FundamentalConstants.hh"

#define LOG(x) std::cout << x << std::endl

using namespace Garfield;

int main(int argc, char *argv[]) {
  // Creates ROOT Application Environtment
  TApplication app("app", &argc, argv);

  const bool debug = true;
  constexpr bool plotSignal = true;
  constexpr bool plotField = true;
  
   static constexpr double fGasGapThickness = 0.025;         // cm - total gas gap thickness            250 micron = 0.25 mm = 0.025 cm
    static constexpr double fGasGapCenterY = 0.0/10.0;           // cm - centered at Y=0                    y = 0 mm =0 cm
    static constexpr double fAnodeCathodeThickness = 0.02/10.0;  // cm - thickness of HV layers             20 micron = 0.02 mm = 0.002 cm
    static constexpr double fStripThickness = 5.00/10.0;         // cm - thickness of copper strips         5 mm = 0.5 cm
    static constexpr double fStripWidthX = 30.0/10;              // cm - readout strip width                28 mm = 2.8 cm
    static constexpr double fDetectorSizeX = 300.0/10.0;         // cm - detector size in X                 300 mm = 30 cm    
    static constexpr double fDetectorSizeZ = 300.0/10.0;         // cm - detector size in Z                 300 mm = 30 cm
    static constexpr double fHoneyCombThickness = 0;     // cm - honey comb layer thickness         0 mm = 0 cm (Remove Honeycomb)
    static constexpr double fMylarThickness = 0.035;         // cm - mylar layer thickness              350 micron = 0.35 mm = 0.035 cm
    static constexpr double fResistiveGlassThickness = 0.12; // cm - resistive glass layer thickness  1.2 mm =0.12 cm
    static constexpr double fAnodeVoltage = 3000.0;               // V - ANODE at +5 kV
    static constexpr double fCathodeVoltage = -3000.0;            // V - CATHODE at -5V
    
    // Calculate the dimension and location of geometry.
    double gasGapTop = fGasGapCenterY + fGasGapThickness/2.0;
    double gasGapBottom = fGasGapCenterY - fGasGapThickness/2.0;
    double fResistiveGlassPosition = fGasGapCenterY + fGasGapThickness/2.0 + fResistiveGlassThickness/2.0;
    double fAnodeCathodePosition = fGasGapCenterY + fGasGapThickness/2.0 + fResistiveGlassThickness + fAnodeCathodeThickness/2.0;
    double fMylarPosition = fGasGapCenterY + fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness/2.0;
    double fStripPosition = fGasGapCenterY + fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness/2.0;

  // Set up the gas (C2H2F4/iC4H10/SF6 90/5/5).
  MediumMagboltz gas;
  gas.LoadGasFile("c2h2f4_ic4h10_sf6.gas");
  gas.Initialise(true);

  // Materials needed.
    MediumPlastic* glass = new Garfield::MediumPlastic();
    MediumConductor* graphite = new Garfield::MediumConductor();
    Medium* mylar = new Garfield::Medium();
    MediumConductor* copper = new Garfield::MediumConductor();

    glass->SetDielectricConstant(10.0);
    mylar->SetDielectricConstant(3.1);

    // Define Geometry
    GeometrySimple* rpc_geometry = new Garfield::GeometrySimple();

    // Define Component
    ComponentNeBem3d *rpc = new ComponentNeBem3d();

    // Readout Strips
    Garfield::SolidBox* Box_TopStrip = new Garfield::SolidBox(0,  fStripPosition ,
                                                               0, fDetectorSizeX/2,
                                                               fStripThickness/2, fDetectorSizeZ/2);
    Box_TopStrip->SetLabel("TopReadout");
    Box_TopStrip->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopStrip,copper);
    const std::string label = "TopReadout";

    Garfield::SolidBox* Box_BottomStrip = new Garfield::SolidBox(0, - fStripPosition,
                                                               0, fDetectorSizeX/2,
                                                               fStripThickness/2, fDetectorSizeZ/2);
    Box_BottomStrip->SetLabel("BottomReadout");
    Box_BottomStrip->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_BottomStrip,copper);
    
    // Mylar Layers
    Garfield::SolidBox* Box_TopMylar = new Garfield::SolidBox(0, fMylarPosition,
                                                               0, fDetectorSizeX/2,
                                                               fMylarThickness/2, fDetectorSizeZ/2);
    Box_TopMylar->SetLabel("TopMylar");
    Box_TopMylar->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopMylar,mylar);

    Garfield::SolidBox* Box_BottomMylar = new Garfield::SolidBox(0, -fMylarPosition,
                                                               0, fDetectorSizeX/2,
                                                               fMylarThickness/2, fDetectorSizeZ/2);
    Box_BottomMylar->SetLabel("BottomMylar");
    Box_BottomMylar->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_BottomMylar,mylar);
    
    // Graphite Layers (Anode and Cathode)
    Garfield::SolidBox* Box_TopGraphite = new Garfield::SolidBox(0, fAnodeCathodePosition,
                                                               0, fDetectorSizeX/2,
                                                               fAnodeCathodeThickness/2, fDetectorSizeZ/2);
    Box_TopGraphite->SetLabel("TopGraphite");
    Box_TopGraphite->SetBoundaryPotential(fAnodeVoltage); // Anode
    rpc_geometry->AddSolid(Box_TopGraphite,graphite);

    Garfield::SolidBox* Box_BottomGraphite = new Garfield::SolidBox(0, -fAnodeCathodePosition,
                                                               0, fDetectorSizeX/2,
                                                               fAnodeCathodeThickness/2, fDetectorSizeZ/2);
    Box_BottomGraphite->SetLabel("BottomGraphite");
    Box_BottomGraphite->SetBoundaryPotential(fCathodeVoltage); // Cathode
    rpc_geometry->AddSolid(Box_BottomGraphite,graphite);

    // Resistive Glass Layers
    Garfield::SolidBox* Box_TopResistiveGlass = new Garfield::SolidBox(0, fResistiveGlassPosition,
                                                               0, fDetectorSizeX/2,
                                                               fResistiveGlassThickness/2, fDetectorSizeZ/2);
    Box_TopResistiveGlass->SetLabel("TopResistiveGlass");
    Box_TopResistiveGlass->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopResistiveGlass,glass);
    Garfield::SolidBox* Box_BottomResistiveGlass = new Garfield::SolidBox(0, -fResistiveGlassPosition,
                                                               0, fDetectorSizeX/2,
                                                               fResistiveGlassThickness/2, fDetectorSizeZ/2);
    Box_BottomResistiveGlass->SetLabel("BottomResistiveGlass");
    Box_BottomResistiveGlass->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_BottomResistiveGlass,glass);

    // Gas Gap
    // Add gas gap volume
    Garfield::SolidBox* Box_GasGap = new Garfield::SolidBox(
        0, fGasGapCenterY, 0,                  // center
        fDetectorSizeX/2,                      // half-width X
        fGasGapThickness/2,                    // half-width Y
        fDetectorSizeZ/2);
    Box_GasGap->SetLabel("GasGap");
    rpc_geometry->AddSolid(Box_GasGap, &gas);           // Use the gas medium

    rpc->SetGeometry(rpc_geometry);
    rpc->SetTargetElementSize(0.01); // 0.01 cm = 0.1mm = 100 microns  
    rpc->SetMinMaxNumberOfElements(1, 10);  //
    rpc->EnableDebugging();
    rpc->Initialise();

  // Create the sensor.
  Sensor sensor(rpc);
  sensor.AddElectrode(rpc, label);

  // Set the time bins.
  const std::size_t nTimeBins = 200;
  const double tmin = 0.;
  const double tmax = 4;
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);

  // Use microscopic tracking for initial stage of the avalanche.
  AvalancheMicroscopic aval(&sensor);
  // Set time until which the calculations will be done microscopically.
  const double tMaxWindow = 4;
  aval.SetTimeWindow(0., tMaxWindow);

  // Use a grid-based method for simulating the avalanche growth 
  // after the initial stage.
  //AvalancheGrid avalgrid(&sensor);
  
  // Start the track in the first gas layer.
  const double dTotal = 2*(fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness);
  const int nX = 5;
  const int nZ = 5;
  const double dY = 1.e-4;
  const int nY = int(dTotal / dY);
  //avalgrid.SetGrid(-0.05, 0.05, nX, -dTotal/2, dTotal/2, nY, -0.05, 0.05, nZ);

  // Prepare the plotting of the induced current and charge.
  ViewSignal *signalView = nullptr;
  TCanvas *cSignal = nullptr;
  if (plotSignal) {
    cSignal = new TCanvas("cSignal", "Current - 1 event", 600, 600);
    signalView = new ViewSignal(&sensor);
    signalView->SetCanvas(cSignal);
  }

  ViewSignal *chargeView = nullptr;
  TCanvas *cCharge = nullptr;
  if (plotSignal) {
    cCharge = new TCanvas("cCharge", "Charge - 1 event", 600, 600);
    chargeView = new ViewSignal(&sensor);
    chargeView->SetCanvas(cCharge);
  }

  // Set up Heed.
  TrackHeed track(&sensor);
  // Set the particle type and momentum [eV/c].
  track.SetParticle("pion");
  track.SetMomentum(7.e9);
  track.CrossInactiveMedia(true);

  std::clock_t start = std::clock();

  auto cD = new TCanvas("cD","Drift",600,600);
  ViewDrift driftView;
  driftView.SetCanvas(cD);
  aval.EnablePlotting(&driftView, 100);
  track.EnablePlotting(&driftView);


  const double y0 = dTotal/2 - fStripThickness - fMylarThickness - fAnodeCathodeThickness - fResistiveGlassThickness;
  // Simulate a charged-particle track.
  track.NewTrack(0, y0-0.00001, 0, 0, 0, -1, 0);
  // Retrieve the clusters along the track.
  for (const auto &cluster : track.GetClusters()) {
    // Loop over the electrons in the cluster.
    for (const auto &electron : cluster.electrons) {
      // Simulate an avalanche (until the set time window).
      aval.AvalancheElectron(electron.x, electron.y, electron.z, electron.t,
                             0., 0., 0., 0.);
      // Transfer the avalanche electrons to the grid.
      //avalgrid.AddElectrons(&aval);
    }
  }
  driftView.Plot();
  // Start grid based avalanche calculation starting from where the microsocpic
  // calculations stopped.
  //LOG("Switching to grid based method.");
  //avalgrid.StartGridAvalanche();
  // Stop timer.
  double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

  LOG("Script: " << "Electrons have drifted. It took " << duration
                 << "s to run.");

  if (plotSignal) {
    // Plot the induced current.
    signalView->PlotSignal(label);
    signalView->EnableLegend();
    cSignal->Update();
    gSystem->ProcessEvents();
    sensor.ExportSignal(label, "Signal");
    // Plot the induced charge.
    sensor.IntegrateSignal(label);
    chargeView->PlotSignal(label);
    cCharge->Update();
    gSystem->ProcessEvents();
    // Export induced current data as an csv file.
    sensor.ExportSignal(label, "Charge");
  }
  /*
  LOG("Script: Total induced charge = " << sensor.GetTotalInducedCharge(label)
                                        << " [fC].");
  LOG("Field Area x = "<< -0.5 <<" to "<<0.5);
  LOG("Field Area y = "<< 0<<" to "<<dTotal);
  */
  // Export and plot Weighting Fields
  ViewField *fieldView = nullptr;
  TCanvas *cField = nullptr;
  
/*
  if (plotField) {
    cField = new TCanvas("cField", "", 600, 600);
    fieldView = new ViewField(&sensor);
    
    fieldView->SetPlane(0,-y0,0,0,-0.001,0);
    fieldView->SetArea(-0.5,0,-0.5,0.5,dTotal,0.5);
    fieldView->SetVoltageRange(-160., 160.);
    
    fieldView->SetCanvas(cField);
    LOG("Calculate Strip Weighting Field - on the GRID");
    rpc->SetWeightingPotentialGrid(-0.5, 0.5, 1, 0, dTotal, 100, -0.5, 0.5, 1, label);
    LOG("Plot Strip Weighting Field - on the GRID");
    cField->SetLeftMargin(0.16);
    fieldView->PlotProfileWeightingField(label,0., 0., 0., 0., dTotal, 0.,"v",true);
    //fieldView->PlotWeightingField(label,"emag");
    //fieldView->PlotContour();  
  }
*/

  LOG("End of Program");
  app.Run(true);
}
