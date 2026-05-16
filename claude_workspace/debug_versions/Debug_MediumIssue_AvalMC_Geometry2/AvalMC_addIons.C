#include <TApplication.h>
#include <TROOT.h>
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
#include <fstream>
#include <string>
// Standard includes

#include <iomanip>
#include <stdexcept>
#include <array>

#include "Garfield/AvalancheGrid.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
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

#include <omp.h>

#define LOG(x) std::cout << x << std::endl

using namespace Garfield;

int main(int argc, char *argv[]) {
  // Creates ROOT Application Environtment
  std::clock_t start = std::clock();
  TApplication app("app", &argc, argv);
  TFile *treefile = new TFile("../result_cache/root_file/aAvalanche.root","recreate");

  TTree *fieldtree = new TTree("ePrimaryTree","Micro: Avalanche Primary Electron Energy Histogram");
  TH1D *ePrimaryEnergy = new TH1D("eAvalancheTree","Micro: Avalanche Primary Electron Energy Distribution",200,0,50);
  ePrimaryEnergy->GetYaxis()->SetTitle("count");
  ePrimaryEnergy->GetXaxis()->SetTitle("Energy (eV)");
  
  TH1D *eSecondaryEnergy = new TH1D("eSecondaryTree","Micro: Avalanche Secondary Electron Energy Distribution",200,0,50);
  eSecondaryEnergy->GetYaxis()->SetTitle("count");
  eSecondaryEnergy->GetXaxis()->SetTitle("Energy (eV)");

  int eventID    = 0;
  int nElectrons = 0;  
  int nSecondaries = 0;
  
  double xp0 = 0., yp0 = 0., zp0 = 0., tp0 = 0., ep0 = 0.;
  double xp1 = 0., yp1 = 0., zp1 = 0., tp1 = 0., ep1 = 0.;
  int    statusPrim = 0;

  std::vector<double> xs0, ys0, zs0, ts0, es0;   // start
  std::vector<double> xs1, ys1, zs1, ts1, es1;   // end
  std::vector<int>    statusSec;

  fieldtree->Branch("eventID",      &eventID,      "eventID/I");
  fieldtree->Branch("nElectrons",   &nElectrons,   "nElectrons/I");
  fieldtree->Branch("nSecondaries", &nSecondaries, "nSecondaries/I");

  // ---- Primary branches (scalar) ----
  fieldtree->Branch("xp0", &xp0, "xp0/D");
  fieldtree->Branch("yp0", &yp0, "yp0/D");
  fieldtree->Branch("zp0", &zp0, "zp0/D");
  fieldtree->Branch("tp0", &tp0, "tp0/D");
  fieldtree->Branch("ep0", &ep0, "ep0/D");
  fieldtree->Branch("xp1", &xp1, "xp1/D");
  fieldtree->Branch("yp1", &yp1, "yp1/D");
  fieldtree->Branch("zp1", &zp1, "zp1/D");
  fieldtree->Branch("tp1", &tp1, "tp1/D");
  fieldtree->Branch("ep1", &ep1, "ep1/D");
  fieldtree->Branch("statusPrim", &statusPrim, "statusPrim/I");
 
  // ---- Secondary branches (vector) ----
  fieldtree->Branch("xs0", &xs0);
  fieldtree->Branch("ys0", &ys0);
  fieldtree->Branch("zs0", &zs0);
  fieldtree->Branch("ts0", &ts0);
  fieldtree->Branch("es0", &es0);
  fieldtree->Branch("xs1", &xs1);
  fieldtree->Branch("ys1", &ys1);
  fieldtree->Branch("zs1", &zs1);
  fieldtree->Branch("ts1", &ts1);
  fieldtree->Branch("es1", &es1);
  fieldtree->Branch("statusSec", &statusSec);

  const bool debug = true;
  constexpr bool plotSignal = true;
  constexpr bool plotField = true;
  
    static constexpr double fGasGapThickness = 0.2;         // cm - total gas gap thickness            2000 micron = 2 mm = 0.2 cm
    static constexpr double fGasGapCenterY = 0.0/10.0;           // cm - centered at Y=0                    y = 0 mm =0 cm
    static constexpr double fAnodeCathodeThickness = 0.02/10.0;  // cm - thickness of HV layers             20 micron = 0.02 mm = 0.002 cm
    static constexpr double fStripThickness = 5.00/10.0;         // cm - thickness of copper strips         5 mm = 0.5 cm
    static constexpr double fStripWidthX = 30.0/10;              // cm - readout strip width                30 mm = 3 cm
    static constexpr double fDetectorSizeX = 300.0/10.0;         // cm - detector size in X                 300 mm = 30 cm    
    static constexpr double fDetectorSizeZ = 300.0/10.0;         // cm - detector size in Z                 300 mm = 30 cm
    static constexpr double fHoneyCombThickness = 0;     // cm - honey comb layer thickness         0 mm = 0 cm (Remove Honeycomb)
    static constexpr double fMylarThickness = 0.035;         // cm - mylar layer thickness              350 micron = 0.35 mm = 0.035 cm
    static constexpr double fResistiveGlassThickness = 0.12; // cm - resistive glass layer thickness  1.2 mm =0.12 cm
    static constexpr double fAnodeVoltage = 5000.0;               // V - ANODE at +5kV
    static constexpr double fCathodeVoltage = -5000.0;            // V - CATHODE at -5V
    
    static constexpr double fReadoutVoltage = 0;		// grounded potential for readout.

    // Calculate the dimension and location of geometry.
    double gasGapTop = fGasGapCenterY + fGasGapThickness/2.0;
    double gasGapBottom = fGasGapCenterY - fGasGapThickness/2.0;
    double fResistiveGlassPosition = fGasGapCenterY + fGasGapThickness/2.0 + fResistiveGlassThickness/2.0;
    double fAnodeCathodePosition = fGasGapCenterY + fGasGapThickness/2.0 + fResistiveGlassThickness + fAnodeCathodeThickness/2.0;
    double fMylarPosition = fGasGapCenterY + fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness/2.0;
    double fStripPosition = fGasGapCenterY + fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness/2.0;
    

  //std::ofstream outFile("configure.txt");

  //if (outFile.is_open()){
	//return 1;
  //}  

  //outFile<< "Configuration"<<std::endl;

  //outFile.close();

  // Set up the gas (C2H2F4/iC4H10/SF6 90/5/5).
  MediumMagboltz gas;
  //gas.LoadGasFile("c2h2f4_ic4h10_sf6.gas");
  gas.SetComposition("C2H2F4", 90, "iC4H10", 5, "SF6", 5);
  gas.SetTemperature(296.15);
  gas.SetPressure(760.0);
  //gas.WriteGasFile("TIFRH_WorkingGas.gas");
  gas.EnableCrossSectionOutput();
  gas.SetMaxElectronEnergy(200.);
  gas.EnablePenningTransfer();
  gas.Initialise(true);
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas.LoadIonMobility(path+"/share/Garfield/Data/IonMobility_SF6-_SF6.txt");


  // Materials needed.
    MediumPlastic* glass = new Garfield::MediumPlastic();
    MediumConductor* graphite = new Garfield::MediumConductor();
    MediumPlastic* mylar = new Garfield::MediumPlastic();
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
    Box_TopStrip->SetBoundaryPotential(fReadoutVoltage);
    rpc_geometry->AddSolid(Box_TopStrip,copper);
    const std::string label1 = "TopReadout";

    Garfield::SolidBox* Box_BottomStrip = new Garfield::SolidBox(0, - fStripPosition,
                                                               0, fDetectorSizeX/2,
                                                               fStripThickness/2, fDetectorSizeZ/2);
    Box_BottomStrip->SetLabel("BottomReadout");
    Box_BottomStrip->SetBoundaryPotential(fReadoutVoltage);
    rpc_geometry->AddSolid(Box_BottomStrip,copper);
    const std::string label2 = "BottomReadout";
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
    Box_TopGraphite->SetBoundaryPotential(fCathodeVoltage); // Cathode
    rpc_geometry->AddSolid(Box_TopGraphite,graphite);
    const std::string label3 = "TopGraphite";

    Garfield::SolidBox* Box_BottomGraphite = new Garfield::SolidBox(0, -fAnodeCathodePosition,
                                                               0, fDetectorSizeX/2,
                                                               fAnodeCathodeThickness/2, fDetectorSizeZ/2);
    Box_BottomGraphite->SetLabel("BottomGraphite");
    Box_BottomGraphite->SetBoundaryPotential(fAnodeVoltage); // Anode
    rpc_geometry->AddSolid(Box_BottomGraphite,graphite);
    const std::string label4 = "BottomGraphite";

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
    rpc->SetNumberOfThreads(16);
    rpc->SetTargetElementSize(0.002); // 0.002 cm = 0.02mm = 20 microns  
    rpc->SetMinMaxNumberOfElements(5, 20);  //
    rpc->EnableDebugging(false);
    //rpc->SetStoreInvMatrix(1);
    rpc->SetReadInvMatrix(1);
    rpc->SetReuseModel();
    //rpc->SetFastVolOptions(1,1,1);
    rpc->Initialise();
    

const std::size_t nx = 50,ny =50,ncont =104;

    ViewField *contourView = nullptr;
  TCanvas *cContour = nullptr;
  cContour = new TCanvas("cContour","Contour",1200,1200);
  contourView = new ViewField();
  contourView->SetCanvas(cContour);
  contourView->SetComponent(rpc);
  contourView->SetPlane(0,-1,0,0,0.,0);
  contourView->SetArea(-14,-0.1,-14,14,0.1,14);
  contourView->SetNumberOfContours(104);
  contourView->PlotContour("emag");
  cContour->SaveAs("../result_cache/png/contouronxz.png");

  // Create the sensor.
  Sensor sensor(rpc);
  sensor.AddElectrode(rpc, label1);
  sensor.AddElectrode(rpc, label2);
  sensor.AddElectrode(rpc, label3);
  sensor.AddElectrode(rpc, label4);


  // Set the time bins.
  const std::size_t nTimeBins = 500;
  const double tmin = 0.;
  const double tmax = 8;
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);
  sensor.ClearSignal();

  /*
  // Set Transfer Function
   auto fT= [](double t){
  	constexpr double tau=25.;
	return (t/tau)*exp(1-t/tau);
  };
  sensor.SetTransferFunction(fT);
  */

   //Shaper shaper(1,25.,1.,"unipolar");
   //sensor.SetTransferFunction(shaper); 

  // Use microscopic tracking for initial stage of the avalanche.
  AvalancheMicroscopic aval(&sensor);
  // Set time until which the calculations will be done microscopically.
  const double tMaxWindow = 8;
  aval.SetTimeWindow(0., tMaxWindow);
  
  aval.EnableSignalCalculation(true);
  aval.UseInducedCharge(true);
  aval.UseWeightingPotential(true);
  aval.EnableElectronEnergyHistogramming(ePrimaryEnergy);
  aval.EnableSecondaryEnergyHistogramming(eSecondaryEnergy);
  aval.EnableDriftLines(true);
  aval.EnableNullCollisionSteps();
  aval.EnableDebugging();
  aval.SetShowProgress(true);
  //aval.SetRunModeOptions(MPRunMode::GPUExclusive,0);
  // Use a grid-based method for simulating the avalanche growth 
  // after the initial stage.
  //AvalancheGrid avalgrid(&sensor);
  
  // For ions
  AvalancheMC avalMCi;
  avalMCi.SetTimeSteps(tstep);
  avalMCi.SetSensor(&sensor);
  //avalMCi.SetAvalancheSteps(0.001);
  avalMCi.SetTimeWindow(0.,tMaxWindow);
  avalMCi.EnableAvalancheSizeLimit(100000);
  avalMCi.EnableSignalCalculation(true);
  avalMCi.EnableDiffusion(true);
  avalMCi.EnableAttachment(true);
  avalMCi.EnableRecombination(true);
  avalMCi.EnableDebugging(true);

  // For negative ions
  AvalancheMC avalMCin;
  avalMCin.SetTimeSteps(tstep);
  avalMCin.SetSensor(&sensor);
  //avalMCin.SetAvalancheSteps(0.001);
  avalMCin.SetTimeWindow(0.,tMaxWindow);
  avalMCin.EnableAvalancheSizeLimit(100000);
  avalMCin.EnableSignalCalculation(true);
  avalMCin.EnableDiffusion(true);
  avalMCin.EnableAttachment(true);
  avalMCin.EnableRecombination(true);
  avalMCin.EnableDebugging(true);

  // Start the track in the first gas layer.
  const double dTotal = 2*(fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness);
  const int nX = 5;
  const int nZ = 5;
  const double dY = 1.e-4;
  const int nY = int(dTotal / dY);
  //avalgrid.SetGrid(-0.05, 0.05, nX, -dTotal/2, dTotal/2, nY, -0.05, 0.05, nZ);
  const double y0 = dTotal/2 - fStripThickness - fMylarThickness - fAnodeCathodeThickness - fResistiveGlassThickness;

  // Prepare the plotting of the induced current and charge.
  /*
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
  */

  // Set up Heed.
  TrackHeed track(&sensor);
  // Set the particle type and momentum [eV/c].
  track.SetParticle("pion");
  track.SetMomentum(7.e9);
  track.CrossInactiveMedia(true);

  auto cD = new TCanvas("cD","Drift",1200,800);
  ViewDrift driftView;
  driftView.SetCanvas(cD);
  aval.EnablePlotting(&driftView, 1000000);
  avalMCi.EnablePlotting(&driftView);
  avalMCin.EnablePlotting(&driftView);
  track.EnablePlotting(&driftView);

  // Simulate a charged-particle track.
  track.NewTrack(0, y0-0.00001, 0, 0, 0, -1, 0);
  // Retrieve the clusters along the track.
  std::size_t nTotal = 0;
  std::size_t nIonMC = 0;
  std::size_t nNegIonMC = 0;
  std::size_t nAvalanche = 0;
  std::size_t nEMicro = 0;
  std::size_t neMicroAval = 0;
  std::size_t nIMicroAval = 0;
  for (const auto &cluster : track.GetClusters()) {
    // Loop over the electrons in the cluster.
    for (const auto &electron : cluster.electrons) {
      // Simulate an avalanche of primary electron(until the set time window).
      ++nEMicro;
      xs0.clear(); ys0.clear(); zs0.clear(); ts0.clear(); es0.clear();
      xs1.clear(); ys1.clear(); zs1.clear(); ts1.clear(); es1.clear();
      statusSec.clear();
      
      xp0 = yp0 = zp0 = tp0 = ep0 = 0.;
      xp1 = yp1 = zp1 = tp1 = ep1 = 0.;
      statusPrim = 0;
      
      // ---- Initial electron ----
      const double xi = 0.05, yi = 0.0, zi = 0.0;
      const double ti = 0.0;
      const double ei = 0.1;
      const double dx = 0.0, dy = 0.0, dz = 0.0;

      aval.AvalancheElectron(electron.x, electron.y, electron.z, electron.t,
                             electron.e, 0., 0., 0.);
      avalMCi.DriftIon(electron.x, electron.y, electron.z, electron.t);
      int ne = 0, ni = 0;
      aval.GetAvalancheSize(ne, ni);
      neMicroAval +=ne;
      nIMicroAval +=ni;
      
      nElectrons = aval.GetNumberOfElectronEndpoints();
      nSecondaries = (nElectrons > 0) ? nElectrons - 1 : 0;
      //eventID      = i;
      
      // ---- Split endpoints into primary (j=0) and secondaries (j>0) ----
    for (int j = 0; j < nElectrons; ++j) {
      double xs, ys, zs, ts, es;
      double xe, ye, ze, te, ee;
      int st;
      aval.GetElectronEndpoint(j,
                               xs, ys, zs, ts, es,
                               xe, ye, ze, te, ee,
                               st);
 
      if (j == 0) {
        // ---- PRIMARY ----
        xp0 = xs; yp0 = ys; zp0 = zs; tp0 = ts; ep0 = es;
        xp1 = xe; yp1 = ye; zp1 = ze; tp1 = te; ep1 = ee;
        statusPrim = st;
      } else {
        // ---- SECONDARY ----
        xs0.push_back(xs); ys0.push_back(ys); zs0.push_back(zs);
        ts0.push_back(ts); es0.push_back(es);
        xs1.push_back(xe); ys1.push_back(ye); zs1.push_back(ze);
        ts1.push_back(te); es1.push_back(ee);
        statusSec.push_back(st);
      }
    }

    fieldtree->Fill();
      // Drift positive ions: AvalancheMicroscopic has no GetIons(); instead,
      // each secondary electron (index > 0) was born from an ionization event,
      // so path[0] of each secondary electron is the ion creation point.
      const auto& allElec = aval.GetElectrons();
      for (std::size_t ie = 1; ie < allElec.size(); ++ie) {
        if (allElec[ie].path.empty()) continue;
        //std::size_t nIe=0,nIi=0;
        const auto& p0 = allElec[ie].path[0];
        //avalMCi.AddIon(p0.x, p0.y, p0.z, p0.t);
        //std::cout<< "Added: Ion MC Avalanche Number: "<<nIonMC<<std::endl;
        //std::cout<< "Added: Number of Electron in AvalMC: " << nIe <<std::endl;
        //std::cout<<  "Added: Number of Ion in AvalMC: " << nIi <<std::endl;
        avalMCi.DriftIon(p0.x, p0.y, p0.z, p0.t);
        ++nIonMC;
        //avalMCi.GetAvalancheSize(nIe,nIi);
        //std::cout<< "Ion MC Avalanche Number: "<<nIonMC<<std::endl;
        //std::cout<< "Number of Electron in AvalMC: " << nIe <<std::endl;
        //std::cout<<  "Number of Ion in AvalMC: " << nIi <<std::endl;
      }
      // Drift negative ions (SF6-): electrons that attached to SF6 molecules.
      // StatusAttached == -7 in Garfield++; verify against your version if needed.
      for (const auto& elec : aval.GetElectrons()) {
        //if (elec.status == -7) {
        //std::size_t nIe=0,nIi=0;
          const auto& pEnd = elec.path.back();
          //avalMCin.AddNegativeIon(pEnd.x, pEnd.y, pEnd.z, pEnd.t);
          //std::cout<< "Added: Ion MC Avalanche Number: "<<nIonMC<<std::endl;
        //std::cout<< "Added: Number of Electron in AvalMC: " << nIe <<std::endl;
        //std::cout<<  "Added: Number of Ion in AvalMC: " << nIi <<std::endl;
          avalMCin.DriftNegativeIon(pEnd.x, pEnd.y, pEnd.z, pEnd.t);
          ++nNegIonMC;
        //avalMCin.GetAvalancheSize(nIe,nIi);
        //std::cout<< "Ion MC Avalanche Number: "<<nNegIonMC<<std::endl;
        //std::cout<< "Number of Electron in AvalMC: " << nIe <<std::endl;
        //std::cout<<  "Number of Ion in AvalMC: " << nIi <<std::endl;
      }
      //Transfer the avalanche electrons to the grid.
      //avalgrid.AddElectrons(&aval);
    }
  }
  driftView.SetColourElectrons(kBlue);
  driftView.SetColourIons(kRed);
  driftView.SetColourNegativeIons(kGreen);
  driftView.Plot();
  cD->Update();
  cD->Write();
  cD->SaveAs("../result_cache/png/DriftPlot.png");
  std::size_t nDriftLine = driftView.GetNumberOfDriftLines();

  // Start grid based avalanche calculation starting from where the microsocpic
  // calculations stopped.
  //LOG("Switching to grid based method.");
  //avalgrid.StartGridAvalanche();
  // Stop timer.
  //double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

  //LOG("Script: " << "Electrons have drifted. It took " << duration<< "s to run.");

  /*               
  if (plotSignal) {
    // Plot the induced current.
    //signalView->PlotSignal(label1);
    //signalView->EnableLegend();
    //cSignal->Update();
    gSystem->ProcessEvents();
    sensor.ExportSignal(label1, "../Signal_"+label1);
    sensor.ExportSignal(label2, "../Signal_"+label2);
    sensor.ExportSignal(label3, "../Signal_"+label3);
    sensor.ExportSignal(label4, "../Signal_"+label4);
  */

  treefile->Write();
  treefile->Close();
  // ---- 2. Save the signals: one tree per electrode, three branches each ----
  const std::vector<std::string> labels = {
    "TopReadout", "BottomReadout", "TopGraphite", "BottomGraphite"
  };
  TFile fout("../result_cache/root_file/signal.root", "RECREATE");

  double TimeStart = 0., TimeStep = 8.;
  std::size_t NumberBins = 0;
  sensor.GetTimeWindow(TimeStart, TimeStep, NumberBins);

  for (const auto& lbl : labels) {
      TTree t(lbl.c_str(),
              ("RPC induced signal at " + lbl).c_str());

      Int_t    timeBin = 0;
      Double_t time    = 0.;   // physical time [ns]
      Double_t sElec   = 0.;
      Double_t sIon    = 0.;
      Double_t sTotal  = 0.;

      t.Branch("TimeBin",        &timeBin, "TimeBin/I");
      t.Branch("Time",           &time,    "Time/D");
      t.Branch("ElectronSignal", &sElec,   "ElectronSignal/D");
      t.Branch("IonSignal",      &sIon,    "IonSignal/D");
      t.Branch("TotalSignal",    &sTotal,  "TotalSignal/D");

      for (unsigned int i = 0; i < NumberBins; ++i) {
          timeBin = static_cast<Int_t>(i);
          time    = TimeStart + (i + 0.5) * TimeStep;   // bin center
          sElec   = sensor.GetElectronSignal(lbl, i);
          sIon    = sensor.GetIonSignal     (lbl, i);
          sTotal  = sensor.GetSignal        (lbl, i);
          t.Fill();
      }

      t.Write();
  }

  fout.Close();

for (const auto& lbl : labels) {
    std::ofstream fout("../result_cache/signal_" + lbl + ".txt");
    fout << "# RPC induced signal at " << lbl << "\n";
    fout << "# TimeBin\tTime[ns]\tElectronSignal\tIonSignal\tTotalSignal\n";
    fout << std::scientific << std::setprecision(6);

    for (unsigned int i = 0; i < NumberBins; ++i) {
        const double time   = TimeStart + (i + 0.5) * TimeStep;
        const double sElec  = sensor.GetElectronSignal(lbl, i);
        const double sIon   = sensor.GetIonSignal     (lbl, i);
        const double sTotal = sensor.GetSignal        (lbl, i);

        fout << i      << '\t'
             << time   << '\t'
             << sElec  << '\t'
             << sIon   << '\t'
             << sTotal << '\n';
    }
}
    /* 
    sensor.ConvoluteSignal(label);
    signalView->PlotSignal(label);
    cSignal->Update();
    gSystem->ProcessEvents();
    sensor.ExportSignal(label,"ConvolutedSignal");
    */
    
    
    // Plot the induced charge.
    //sensor.IntegrateSignal(label);
   // chargeView->PlotSignal(label);
    //cCharge->Update();
    //gSystem->ProcessEvents();
    // Export induced current data as an csv file.
    //sensor.ExportSignal(label, "Charge");
    //gSystem->ProcessEvents();    
  
  /*
  LOG("Script: Total induced charge = " << sensor.GetTotalInducedCharge(label)
                                        << " [fC].");
  LOG("Field Area x = "<< -0.5 <<" to "<<0.5);
  LOG("Field Area y = "<< 0<<" to "<<dTotal);
  */
  // Export and plot Weighting Fields
  //ViewField *fieldView = nullptr;
  //TCanvas *cField = nullptr;
  
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
  std::cout << "Number of Positive Ion Drifts:  " << nIonMC << std::endl;
  std::cout << "Number of Negative Ion Drifts:  " << nNegIonMC << std::endl;
  std::cout << "Number of Primary Electron (Microscopic): " << nEMicro << std::endl;
  std::cout << "Time Bin Size: " << tstep << std::endl;
  std::cout << "Number of Drift Line: " << nDriftLine << std::endl;
  double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

  LOG("Script: " << "Total Run Time: " << duration);
  //std::cout<< path+"share/Garfield/Data/IonMobility_SF6-_SF6.txt" <<std::endl;
  LOG("End of Program");
  app.Run(true);

  
}
