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

struct ClusterPhoton {
  int ClusterId;
  double Photon_x;
  double Photon_y;
  double Photon_z;
  double Photon_t;
  double Photon_energy;
  double Photon_dx;
  double Photon_dy;
  double Photon_dz;
};


struct ClusterElectron {
  int ClusterId;
  double Electron_x;
  double Electron_y;
  double Electron_z;
  double Electron_t;
  double Electron_energy;
  double Electron_dx;
  double Electron_dy;
  double Electron_dz;
};

struct ClusterIon {
  int ClusterId;
  double Ion_x;
  double Ion_y;
  double Ion_z;
  double Ion_t;
  double Ion_energy;
  double Ion_dx;
  double Ion_dy;
  double Ion_dz;
};


int main(int argc, char *argv[]) {
  // Creates ROOT Application Environtment
  std::ofstream OutDebug("../result_cache/output_debug.txt");
  //TApplication app("app", &argc, argv);
  TFile *treefile = new TFile("../result_cache/root_file/aAvalanche.root","recreate");

  TTree *fieldtree = new TTree("ePrimaryTree","Micro: Avalanche Electron energy");
  TH1D *ePrimaryEnergy = new TH1D("eAvalancheTree","Micro: Avalanche Primary Electron Energy Distribution",200,0,50);
  ePrimaryEnergy->GetYaxis()->SetTitle("count");
  ePrimaryEnergy->GetXaxis()->SetTitle("Energy (eV)");
  
  TH1D *eSecondaryEnergy = new TH1D("eSecondaryTree","Micro: Avalanche Secondary Electron Energy Distribution",200,0,50);
  eSecondaryEnergy->GetYaxis()->SetTitle("count");
  eSecondaryEnergy->GetXaxis()->SetTitle("Energy (eV)");

  TFile *clusterfile = new TFile("../result_cache/root_file/cluster_data.root","recreate");
  TTree *treecluster = new TTree("eClusterTree","Heed Ionisation Cluster Data Tree");

  double cluster_x = 0, cluster_y = 0, cluster_z = 0, cluster_t = 0, cluster_energy = 0, cluster_extra= 0;
  struct ClusterElectron cluster_electron;
  struct ClusterIon cluster_ion; 
  struct ClusterPhoton cluster_photon;

  treecluster->Branch("cluster.x",&cluster_x,"cluster.x/D");
  treecluster->Branch("cluster.y",&cluster_y,"cluster.y/D");
  treecluster->Branch("cluster.z",&cluster_z,"cluster.z/D");
  treecluster->Branch("cluster.t",&cluster_t,"cluster.t/D");
  treecluster->Branch("cluster.energy",&cluster_energy,"cluster.energy/D");
  treecluster->Branch("cluster.extra",&cluster_extra,"cluster.extra/D");
  treecluster->Branch("cluster_electron", &cluster_electron,  "ClusterId/I:Electron_x/D:Electron_y/D:Electron_z/D:Electron_t/D:Electron_energy/D:Electron_dx/D:Electron_dy/D:Electron_dz/D");
  treecluster->Branch("cluster_ion", &cluster_ion,"ClusterId/I:Ion_x/D:Ion_y/D:Ion_z/D:Ion_t/D:""Ion_energy/D:Ion_dx/D:Ion_dy/D:Ion_dz/D");
  treecluster->Branch("cluster_photon", &cluster_photon,"ClusterId/I:Photon_x/D:Photon_y/D:Photon_z/D:Photon_t/D:""Photon_energy/D:Photon_dx/D:Photon_dy/D:Photon_dz/D");


  int eventID    = 0;
  int nElectrons = 0;  
  int nSecondaries = 0;
  
  double xp0 = 0., yp0 = 0., zp0 = 0., tp0 = 0., ep0 = 0.;
  double xp1 = 0., yp1 = 0., zp1 = 0., tp1 = 0., ep1 = 0.;
  int    statusPrim = 0;

  std::vector<double> xs0, ys0, zs0, ts0, es0;   // start
  std::vector<double> xs1, ys1, zs1, ts1, es1;   // end
  std::vector<int>    statusSec;

  double bAlpha, bEta, bAlphaEff;

  double bVdElec, bVspeedElec, bMuElec;
  int    bNElecVel;
 
  double Lengthtotal;

  fieldtree->Branch("Ltot",  &Lengthtotal,  "Ltot/D");
  fieldtree->Branch("alpha",     &bAlpha,     "alpha/D");
  fieldtree->Branch("eta",       &bEta,       "eta/D");
  fieldtree->Branch("alphaEff",  &bAlphaEff,  "alphaEff/D");


  //fieldtree->Branch("vdElec",     &bVdElec,     "vdElec/D");        // cm/ns, along drift dir
  //fieldtree->Branch("vspeedElec", &bVspeedElec, "vspeedElec/D");    // cm/ns, |path|/time
  //fieldtree->Branch("muElec",     &bMuElec,     "muElec/D");        // cm^2/(V*ns)
  //fieldtree->Branch("nElecVel",   &bNElecVel,   "nElecVel/I");

  TH1F* hAlpha    = new TH1F("hAlpha",    "Townsend alpha;#alpha [1/cm];events",   200,   0., 200.);
  TH1F* hEta      = new TH1F("hEta",      "Attachment eta;#eta [1/cm];events",     200,   0., 100.);
  TH1F* hAlphaEff = new TH1F("hAlphaEff", "Effective alpha;#alpha_{eff} [1/cm];events", 200, -50., 200.);
  //TH1F* hLogGain  = new TH1F("hLogGain",  "log10(gain);log_{10} n_{e};events",      200,   0.,  10.);

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
    static constexpr double fStripThickness = 0.005;         // cm - thickness of copper strips         50 micron = 0.05 mm = 0.005 cm
    static constexpr double fStripWidthX = 300.0/10;              // cm - readout strip width                30 mm = 3 cm
    static constexpr double fDetectorSizeX = 300.0/10.0;         // cm - detector size in X                 300 mm = 30 cm    
    static constexpr double fDetectorSizeZ = 300.0/10.0;         // cm - detector size in Z                 300 mm = 30 cm
    static constexpr double fHoneyCombThickness = 0;     // cm - honey comb layer thickness         0 mm = 0 cm (Remove Honeycomb)
    static constexpr double fMylarThickness = 0.01;         // cm - mylar layer thickness              100 micron = 0.1 mm = 0.01 cm
    static constexpr double fResistiveGlassThickness = 0.3; // cm - resistive glass layer thickness  3 mm = 0.3 cm
    static constexpr double fAnodeVoltage = 5500.0;               // V - ANODE at +7kV
    static constexpr double fCathodeVoltage = -5500.0;            // V - CATHODE at -7V
    
    static constexpr double fReadoutVoltage = 0;		// grounded potential for readout.

    // Calculate the dimension and location of geometry.
    double gasGapTop = fGasGapCenterY + fGasGapThickness/2.0;
    double gasGapBottom = fGasGapCenterY - fGasGapThickness/2.0;
    double fResistiveGlassPosition = fGasGapCenterY + fGasGapThickness/2.0 + fResistiveGlassThickness/2.0;
    double fAnodeCathodePosition = fGasGapCenterY + fGasGapThickness/2.0 + fResistiveGlassThickness + fAnodeCathodeThickness/2.0;
    double fMylarPosition = fGasGapCenterY + fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness/2.0;
    double fStripPosition = fGasGapCenterY + fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness/2.0;

  // Set up the gas (C2H2F4/iC4H10/SF6 90/5/5).
  MediumMagboltz gas;
  gas.LoadGasFile("rpc_95.5_4.2_0.3_40-45_ncoll=20.gas");
  gas.SetComposition("C2H2F4", 95.5, "iC4H10", 4.2, "SF6", 0.3);
  gas.SetTemperature(296.15);
  gas.SetPressure(760.0);
  gas.EnableCrossSectionOutput();
  gas.SetMaxElectronEnergy(1000.);
  gas.EnablePenningTransfer();
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas.LoadIonMobility(path+"/share/Garfield/Data/IonMobility_SF6-_SF6.txt");
  gas.LoadIonMobility(path+"/share/Garfield/Data/IonMobility_C8Hn+_iC4H10.txt");
  gas.Initialise(true);


  // Materials needed.
    MediumPlastic* glass = new Garfield::MediumPlastic();
   //Medium* graphite = new Garfield::Medium();
   MediumConductor* graphite = new Garfield::MediumConductor();
    MediumPlastic* mylar = new Garfield::MediumPlastic();
    MediumConductor* copper = new Garfield::MediumConductor();
        
   // graphite->SetDielectricConstant(12.0);
    glass->SetDielectricConstant(10.0);
    mylar->SetDielectricConstant(3.1);

    // Define Geometry
    GeometrySimple* rpc_geometry = new Garfield::GeometrySimple();

    // Define Component
    ComponentNeBem3d *rpc = new ComponentNeBem3d();

    // Readout Strips
    Garfield::SolidBox* Box_TopStrip = new Garfield::SolidBox(0, 0,fStripPosition, fDetectorSizeX/2, fDetectorSizeZ/2,fStripThickness/2);
    Box_TopStrip->SetLabel("TopReadout");
    Box_TopStrip->SetBoundaryPotential(fReadoutVoltage);
    rpc_geometry->AddSolid(Box_TopStrip,copper);
    const std::string label1 = "TopReadout";

    Garfield::SolidBox* Box_BottomStrip = new Garfield::SolidBox(0,0,- fStripPosition, fDetectorSizeX/2,fDetectorSizeZ/2,fStripThickness/2);
    Box_BottomStrip->SetLabel("BottomReadout");
    Box_BottomStrip->SetBoundaryPotential(fReadoutVoltage);
    rpc_geometry->AddSolid(Box_BottomStrip,copper);
    const std::string label2 = "BottomReadout";
    // Mylar Layers
    Garfield::SolidBox* Box_TopMylar = new Garfield::SolidBox(0,0, fMylarPosition
                                                               , fDetectorSizeX/2
                                                               , fDetectorSizeZ/2,fMylarThickness/2);
    Box_TopMylar->SetLabel("TopMylar");
    Box_TopMylar->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopMylar,mylar);

    Garfield::SolidBox* Box_BottomMylar = new Garfield::SolidBox(0,0, -fMylarPosition,
                                                                fDetectorSizeX/2, fDetectorSizeZ/2,
                                                               fMylarThickness/2);
    Box_BottomMylar->SetLabel("BottomMylar");
    Box_BottomMylar->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_BottomMylar,mylar);
    
    // Graphite Layers (Anode and Cathode)
    Garfield::SolidBox* Box_TopGraphite = new Garfield::SolidBox(0,0, fAnodeCathodePosition,
                                                              fDetectorSizeX/2, fDetectorSizeZ/2,
                                                               fAnodeCathodeThickness/2);
    Box_TopGraphite->SetLabel("TopGraphite");
    Box_TopGraphite->SetBoundaryPotential(fCathodeVoltage); // Cathode
    //Box_TopGraphite->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopGraphite,graphite);
    const std::string label3 = "TopGraphite";

    Garfield::SolidBox* Box_BottomGraphite = new Garfield::SolidBox(0,0, -fAnodeCathodePosition,
                                                                fDetectorSizeX/2,fDetectorSizeZ/2,
                                                               fAnodeCathodeThickness/2);
    Box_BottomGraphite->SetLabel("BottomGraphite");
    //Box_BottomGraphite->SetBoundaryDielectric();
    Box_BottomGraphite->SetBoundaryPotential(fAnodeVoltage); // Anode
    
    rpc_geometry->AddSolid(Box_BottomGraphite,graphite);
    const std::string label4 = "BottomGraphite";

    // Resistive Glass Layers
    Garfield::SolidBox* Box_TopResistiveGlass = new Garfield::SolidBox(0,0,fResistiveGlassPosition,
                                                              fDetectorSizeX/2, fDetectorSizeZ/2,
                                                               fResistiveGlassThickness/2);
    Box_TopResistiveGlass->SetLabel("TopResistiveGlass");
    Box_TopResistiveGlass->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopResistiveGlass,glass);
    Garfield::SolidBox* Box_BottomResistiveGlass = new Garfield::SolidBox(0,0,  -fResistiveGlassPosition,
                                                               fDetectorSizeX/2,
                                                               fDetectorSizeZ/2,fResistiveGlassThickness/2);
    Box_BottomResistiveGlass->SetLabel("BottomResistiveGlass");
    Box_BottomResistiveGlass->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_BottomResistiveGlass,glass);

    // Gas Gap
    // Add gas gap volume
    Garfield::SolidBox* Box_GasGap = new Garfield::SolidBox(
        0,0, fGasGapCenterY,                   // center
        fDetectorSizeX/2,                      // half-width X
        fDetectorSizeZ/2,fGasGapThickness/2);
    Box_GasGap->SetLabel("GasGap");
    rpc_geometry->AddSolid(Box_GasGap, &gas);           // Use the gas medium

    rpc->SetGeometry(rpc_geometry);
    rpc->SetNumberOfThreads(16);
    rpc->SetTargetElementSize(0.002); // 0.002 cm = 0.02mm = 20 microns  
    rpc->SetMinMaxNumberOfElements(5, 20);  //
    rpc->EnableDebugging();
    rpc->SetStoreInvMatrix(1);
    //rpc->SetReadInvMatrix(1);
    //rpc->SetReuseModel();
    //rpc->SetFastVolOptions(1,1,1);
    rpc->Initialise();


const std::size_t nx = 50,ny =50,ncont =104;

  ViewField *contourView = nullptr;
  TCanvas *cContour = nullptr;
  cContour = new TCanvas("cContour","Contour",1200,1000);
  contourView = new ViewField();
  contourView->SetCanvas(cContour);
  contourView->SetComponent(rpc);
  contourView->SetPlane(0,0,-1,0.,0.,0.);
  contourView->SetArea(-18,-18,-0.1,18,18,0.1);
  contourView->SetNumberOfContours(104);
  contourView->PlotContour("emag");
  cContour->SaveAs("../result_cache/png/contouronxy.png");

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
  
  const double tMaxWindow = 8;
   /*
  // Use microscopic tracking for initial stage of the avalanche.
  AvalancheMicroscopic aval(&sensor);
  //aval.SetRunModeOptions(MPRunMode::GPUExclusive,0);
  const double tMaxWindow = 8;
  aval.SetTimeWindow(0., tMaxWindow);
  
  aval.EnableSignalCalculation(true);
  aval.UseInducedCharge(true);
  aval.UseWeightingPotential(true);
  aval.EnableElectronEnergyHistogramming(ePrimaryEnergy);
  aval.EnableSecondaryEnergyHistogramming(eSecondaryEnergy);
  aval.EnableDriftLines(true);
  aval.EnableNullCollisionSteps();
  //aval.EnableDebugging();
  aval.SetShowProgress(false);
 */
  // Use a grid-based method for simulating the avalanche growth 
  // after the initial stage.
  //AvalancheGrid avalgrid(&sensor);
  
  // For ions
  AvalancheMC avalMCi;
  avalMCi.SetTimeSteps(tstep);
  avalMCi.SetDistanceSteps(0.001);
  avalMCi.SetSensor(&sensor);
  //avalMCi.SetAvalancheSteps(0.001);
  avalMCi.SetTimeWindow(0.,tMaxWindow);
  avalMCi.EnableAvalancheSizeLimit(100000);
  avalMCi.EnableSignalCalculation(true);
  avalMCi.EnableDiffusion(true);
  avalMCi.EnableAttachment(true);
  avalMCi.EnableRecombination(true);
  const unsigned int nThreads = std::thread::hardware_concurrency();
  avalMCi.EnableMultithreading(nThreads);
  avalMCi.EnableDebugging(true);

  // For negative ions
  AvalancheMC avalMCin;
  avalMCin.SetTimeSteps(tstep);
  avalMCin.SetSensor(&sensor);
  avalMCin.SetDistanceSteps(0.001);
  //avalMCin.SetAvalancheSteps(0.001);
  avalMCin.SetTimeWindow(0.,tMaxWindow);
  avalMCin.EnableAvalancheSizeLimit(100000);
  avalMCin.EnableSignalCalculation(true);
  avalMCin.EnableDiffusion(true);
  avalMCin.EnableAttachment(true);
  avalMCin.EnableRecombination(true);
  avalMCin.EnableDebugging(true);
  avalMCin.EnableMultithreading(nThreads);

  // Start the track in the first gas layer.
  const double dTotal = 2*(fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness);
  const int nX = 5;
  const int nZ = 5;
  const double dY = 1.e-4;
  const int nY = int(dTotal / dY);
  //avalgrid.SetGrid(-0.05, 0.05, nX, -dTotal/2, dTotal/2, nY, -0.05, 0.05, nZ);
  const double y0 = dTotal/2 - fStripThickness - fMylarThickness - fAnodeCathodeThickness - fResistiveGlassThickness;

  // Set up Heed.
  TrackHeed track(&sensor);
  // Set the particle type and momentum [eV/c].
  track.SetParticle("Pion");
  track.SetMomentum(7.e9);
  track.CrossInactiveMedia(true);

  std::clock_t start = std::clock();

  TCanvas* cD = new TCanvas("cD","Drift3D",1200,800);
  TCanvas* cD2DXY = new TCanvas("cD2DXY","Drift2DXY",1200,800);
  TCanvas* cD2DXZ = new TCanvas("cD2DXZ","Drift2DXZ",1200,800);

  ViewDrift driftView;
  driftView.SetCanvas(cD);
  driftView.EnableClusterMarkers(true);
  //aval.EnablePlotting(&driftView, 1000000);
  avalMCi.EnablePlotting(&driftView);
  avalMCin.EnablePlotting(&driftView);
  track.EnablePlotting(&driftView);
 
  
  // Simulate a charged-particle track.
track.NewTrack(0, 0, y0-0.00001, 0, 0, 0,-1);
  // Retrieve the clusters along the track.
  std::size_t nIonMC = 0;
  std::size_t nNegIonMC = 0;
  std::size_t nEMicro = 0;
  std::size_t nCluster = 0;
  std::size_t neMicroAval = 0;
  std::size_t nIMicroAval = 0;
  for (const auto &cluster : track.GetClusters()) {
    // Loop over each cluster generated by Track Heed.
    ++nCluster;
    std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    std::cout<<"Cluster: "<<nCluster<<std::endl;
    OutDebug<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
    OutDebug<<"Cluster: "<<nCluster<<"\n";
    cluster_x = cluster.x;
    cluster_y = cluster.y;
    cluster_z = cluster.z;
    cluster_t = cluster.t;
    cluster_energy = cluster.energy;
    cluster_extra = cluster.extra;
    treecluster->Fill();
    int clusterID = 0;
    for (const auto &ion :cluster.ions){
      clusterID = nCluster;
      cluster_ion.Ion_x = ion.x;
      cluster_ion.Ion_y = ion.y;
      cluster_ion.Ion_z = ion.z;
      cluster_ion.Ion_t = ion.t;
      cluster_ion.Ion_energy = ion.e;
      cluster_ion.Ion_dx = ion.dx;
      cluster_ion.Ion_dx = ion.dy;
      cluster_ion.Ion_dy = ion.dz;
      treecluster->Fill();
    }
    for (const auto &photon :cluster.photons){
      clusterID = nCluster;
      cluster_photon.Photon_x = photon.x;
      cluster_photon.Photon_y = photon.y;
      cluster_photon.Photon_z = photon.z;
      cluster_photon.Photon_t = photon.t;
      cluster_photon.Photon_energy = photon.e;
      cluster_photon.Photon_dx = photon.dx;
      cluster_photon.Photon_dx = photon.dy;
      cluster_photon.Photon_dy = photon.dz;
      treecluster->Fill();
    }
    for (const auto &electron : cluster.electrons) {
      clusterID = nCluster;
      cluster_electron.Electron_x = electron.x;
      cluster_electron.Electron_y = electron.y;
      cluster_electron.Electron_z = electron.z;
      cluster_electron.Electron_t = electron.t;
      cluster_electron.Electron_energy = electron.e;
      cluster_electron.Electron_dx = electron.dx;
      cluster_electron.Electron_dx = electron.dy;
      cluster_electron.Electron_dy = electron.dz;
      treecluster->Fill();
    }
    for (const auto &electron : cluster.electrons) {
      //Set Tree Variables, clear value each cluster.
      ++nEMicro; // number of cluster loops
      std::cout<<"Primary Electron: "<<nEMicro<<std::endl;
      OutDebug<<"Primary Electron: "<<nEMicro<<"\n";
      xs0.clear(); ys0.clear(); zs0.clear(); ts0.clear(); es0.clear();
      xs1.clear(); ys1.clear(); zs1.clear(); ts1.clear(); es1.clear();
      statusSec.clear();
      
      xp0 = yp0 = zp0 = tp0 = ep0 = 0.;
      xp1 = yp1 = zp1 = tp1 = ep1 = 0.;
      statusPrim = 0;

      //aval.AvalancheElectron(electron.x, electron.y, electron.z, electron.t,electron.e, 0., 0., 0.);
      avalMCi.AvalancheElectron(electron.x, electron.y, electron.z, electron.t,false,1);
      //avalMCi.DriftIon(electron.x, electron.y, electron.z, electron.t);
      std::size_t ne = 0, ni = 0;
      avalMCi.GetAvalancheSize(ne, ni);
      std::cout<<"ne: "<<ne<<std::endl;
      std::cout<<"ni: "<<ni<<std::endl;
      OutDebug<<"ne: "<<ne<<"\n";
      OutDebug<<"ni: "<<ni<<"\n";
      neMicroAval +=ne; // add number of generated secondaries each cluster loop
      nIMicroAval +=ni; // add numbers of generated ions each cluster loop 
      std::cout<<"Accumulated ne: "<<ne<<std::endl;
      std::cout<<"Accumulated ni: "<<ni<<std::endl;
      OutDebug<<"Accumulated ne: "<<ne<<"\n";
      OutDebug<<"Accumulated ni: "<<ni<<"\n";
      //nElectrons = aval.GetNumberOfElectronEndpoints();
      //nElectrons = aval.GetNumberOfElectronEndpointsGPU();
      //std::cout<<"length of <vector> Electron: "<<nElectrons<<std::endl;
      //OutDebug<<"length of <vector> Electron: "<<nElectrons<<"\n";
      //nSecondaries = (nElectrons > 0) ? nElectrons - 1 : 0;
      // ---- Split endpoints into primary (j=0) and secondaries (j>0) ----
    /*
     for (int j = 0; j < nElectrons; ++j) {
      double xs, ys, zs, ts, es;
      double xe, ye, ze, te, ee;
      int st;
      aval.GetElectronEndpoint(j,xs, ys, zs, ts, es,xe, ye, ze, te, ee,st);
      //aval.GetElectronEndpoint(j,xs, ys, zs, ts, es,xe, ye, ze, te, ee,st);
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
    
      // Drift positive ions: AvalancheMicroscopic has no GetIons(); instead,
      // each secondary electron (index > 0) was born from an ionization event,
      // so path[0] of each secondary electron is the ion creation point.
      const auto& allElec = aval.GetElectrons(); // 
      std::cout<<"length of aval.GetElectron: "<<allElec.size()<<std::endl;
      for (std::size_t ie = 1; ie < allElec.size(); ++ie) {
        if (allElec[ie].path.empty()) continue;
        const auto& p0 = allElec[ie].path[0];
        avalMCi.DriftIon(p0.x, p0.y, p0.z, p0.t);
        ++nIonMC; // Count number of AvalMC loop for secondary ions
        std::cout<<"# of Sec Cation: "<<nIonMC<<std::endl; 
      }
      // Drift negative ions (SF6-): electrons that attached to SF6 molecules.
      // StatusAttached == -7 in Garfield++; verify against your version if needed.
      for (const auto& elec : aval.GetElectrons()) { // loop over all electrons (primary and secondary) by AvalancheElectron
        if (elec.status == -7) {
          const auto& pEnd = elec.path.back();
          avalMCin.DriftNegativeIon(pEnd.x, pEnd.y, pEnd.z, pEnd.t);
          ++nNegIonMC; // Count number of AvalMC loop for secondary negative ions by electron attachment to SF6 
          std::cout<<"# of reattachment (including Primary): "<<nNegIonMC<<std::endl;
	}
      }
    }
    */
     
     /* 
     double Ltot = 0.;
      if (!aval.GetElectrons().empty()) Ltot += aval.GetElectrons()[0].pathLength;
      Lengthtotal = Ltot;
      const double alpha     = (Ltot > 0.) ? double(ni)      / Ltot : 0.;
      const double eta       = (Ltot > 0.) ? double(nNegIonMC) / Ltot : 0.;
      const double alpha_eff = alpha - eta;
      bAlpha    = alpha;
      bEta      = eta;
      bAlphaEff = alpha_eff;
      if (Ltot > 0.) {
          hAlpha   ->Fill(alpha);
          hEta     ->Fill(eta);
          hAlphaEff->Fill(alpha_eff);
        }
      fieldtree->Fill();
     */ 
    }
  }
    fieldtree->Write();
    hAlpha   ->Write();
    hEta     ->Write();
    hAlphaEff->Write();

  //aval.EnableIonisationMarkers(true);
  
  clusterfile->cd();
  clusterfile->Write();
  clusterfile->Close();

  driftView.SetArea(-15, -15, -0.1, 15, 15, 0.1);
  driftView.SetColourElectrons(kBlue);
  driftView.SetColourIons(kRed);
  driftView.SetColourNegativeIons(kGreen);
  driftView.Plot3d(true,false,false);
  cD->Update();
  cD->SaveAs("../result_cache/png/DriftPlot3D.png");
  std::size_t nDriftLine = driftView.GetNumberOfDriftLines();
  std::cout<<"Number of DriftLine: "<<nDriftLine<<std::endl;
  
  driftView.SetCanvas(cD2DXY);
  driftView.SetPlaneXY();
  driftView.SetArea(-15,-15,15,15);
  driftView.SetColourElectrons(kBlue);
  driftView.SetColourIons(kRed);
  driftView.SetColourNegativeIons(kGreen);
  driftView.Plot2d(true,false);
  cD2DXY->Update();
  cD2DXY->SaveAs("../result_cache/png/DriftPlot2DXY.png");
  
  driftView.SetCanvas(cD2DXZ);
  driftView.SetPlaneXY();
  driftView.SetArea(-15,-15,15,15);
  driftView.SetColourElectrons(kBlue);
  driftView.SetColourIons(kRed);
  driftView.SetColourNegativeIons(kGreen);
  driftView.Plot2d(true,false);
  cD2DXZ->Update();
  cD2DXZ->SaveAs("../result_cache/png/DriftPlot2DXZ.png");

  if (plotSignal) {
    gSystem->ProcessEvents();
    sensor.ExportSignal(label1, "../result_cache/sensor_out/Signal_"+label1);
    sensor.ExportSignal(label2, "../result_cache/sensor_out/Signal_"+label2);
    sensor.ExportSignal(label3, "../result_cache/sensor_out/Signal_"+label3);
    sensor.ExportSignal(label4, "../result_cache/sensor_out/Signal_"+label4);
  treefile->cd();  
  cD->Write();         
  cD2DXY->Write();
  cD2DXZ->Write();
  treefile->Write();
  treefile->Close();
  }

    const std::vector<std::string> labels = {
    "TopReadout", "BottomReadout", "TopGraphite", "BottomGraphite"
  };
  TFile fout("../result_cache/root_file/signal.root", "RECREATE");

  for (const auto& lbl : labels) {
      std::vector<double> Timebin, sElec, sIon, sTotal;
      Timebin.reserve(nTimeBins);
      sElec .reserve(nTimeBins);
      sIon  .reserve(nTimeBins);
      sTotal.reserve(nTimeBins);

      for (unsigned int i = 0; i < nTimeBins; ++i) {
          Timebin.push_back(tmin+i*tstep);
          sElec .push_back(sensor.GetElectronSignal(lbl, i));
          sIon  .push_back(sensor.GetIonSignal     (lbl, i));
          sTotal.push_back(sensor.GetSignal        (lbl, i));
      }

      TTree t(lbl.c_str(),
              ("RPC induced signal at " + lbl).c_str());

      // Pointers to the local vectors — they must stay alive until Fill+Write
      std::vector<double>* pTimebin = &Timebin;
      std::vector<double>* pElec  = &sElec;
      std::vector<double>* pIon   = &sIon;
      std::vector<double>* pTotal = &sTotal;
      t.Branch("Time Bin", &Timebin);
      t.Branch("Electron Signal", &pElec);
      t.Branch("Ion Signal",      &pIon);
      t.Branch("Total Signal",    &pTotal);

      t.Fill();    // exactly one entry — the single event
      t.Write();
  }

  fout.Close();

  double TimeStart = 0., TimeStep = 8.;
  std::size_t NumberBins = 0;
  sensor.GetTimeWindow(TimeStart, TimeStep, NumberBins);

  for (const auto& lbl : labels) {
    std::ofstream fout("../result_cache/signal_" + lbl + ".txt");
    fout << "# RPC induced signal at " << lbl << "\n";
    fout << "# TimeBin\tTime[ns]\tElectronSignal\tIonSignal\tTotalSignal\n";
    fout << std::scientific << std::setprecision(6);

    for (unsigned int i = 0; i < NumberBins; ++i) {
        const double time   = TimeStart + (i) * TimeStep;
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
  
  std::cout << "++++ Aval Microscopic ++++"<< std::endl;
  std::cout << "Number of Electron (Avalanche):  " << neMicroAval << std::endl;
  std::cout << "Number of Ion (Avalanche):  " << nIMicroAval << std::endl;
  std::cout << std::endl;
  std::cout << "++++ Aval MC ++++"<< std::endl;
  std::cout << "Number of Positive Ion Drifts:  " << nIonMC << std::endl;
  std::cout << "Number of Negative Ion Drifts:  " << nNegIonMC << std::endl;
  std::cout << "Number of Cluster: " << nEMicro << std::endl;
  std::cout << "Time Bin Size: " << tstep << std::endl;
  std::cout << "Number of Drift Line: " << nDriftLine << std::endl;

  //std::cout<< path+"share/Garfield/Data/IonMobility_SF6-_SF6.txt" <<std::endl;
  LOG("End of Program");
  OutDebug.close();
  //app.Run(true);
}
