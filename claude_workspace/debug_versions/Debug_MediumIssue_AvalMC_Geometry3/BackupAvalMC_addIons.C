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
  TApplication app("app", &argc, argv);

  // ----------------------------------------------------------------------
  // OpenMP run configuration
  // ----------------------------------------------------------------------
  // Set RUN_PARALLEL=false to fall back to a fully serial run (with all
  // debug/plotting features the original code had). Set true to use
  // OpenMP across the per-primary-electron loop.
  const bool RUN_PARALLEL = true;
  const int  nThreads     = RUN_PARALLEL ? omp_get_max_threads() : 1;
  if (RUN_PARALLEL) {
    omp_set_num_threads(nThreads);
    std::cout << "[OpenMP] Running with " << nThreads << " threads.\n";
  } else {
    std::cout << "[OpenMP] Disabled — running serially.\n";
  }

  TFile *treefile = new TFile("../result_cache/root_file/aAvalanche.root","recreate");

  TTree *fieldtree = new TTree("ePrimaryTree","Micro: Avalanche Electron energy");
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

  double bAlpha, bEta, bAlphaEff;

  double bVdElec, bVspeedElec, bMuElec;
  int    bNElecVel;

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

  fieldtree->Branch("alpha",     &bAlpha,     "alpha/D");
  fieldtree->Branch("eta",       &bEta,       "eta/D");
  fieldtree->Branch("alphaEff",  &bAlphaEff,  "alphaEff/D");

  fieldtree->Branch("vdElec",     &bVdElec,     "vdElec/D");        // cm/ns, along drift dir
  fieldtree->Branch("vspeedElec", &bVspeedElec, "vspeedElec/D");    // cm/ns, |path|/time
  fieldtree->Branch("muElec",     &bMuElec,     "muElec/D");        // cm^2/(V*ns)
  fieldtree->Branch("nElecVel",   &bNElecVel,   "nElecVel/I");


  TH1F* hAlpha    = new TH1F("hAlpha",    "Townsend alpha;#alpha [1/cm];events",   200,   0., 200.);
  TH1F* hEta      = new TH1F("hEta",      "Attachment eta;#eta [1/cm];events",     200,   0., 100.);
  TH1F* hAlphaEff = new TH1F("hAlphaEff", "Effective alpha;#alpha_{eff} [1/cm];events", 200, -50., 200.);
  TH1F* hLogGain  = new TH1F("hLogGain",  "log10(gain);log_{10} n_{e};events",      200,   0.,  10.);

  const bool debug = true;
  constexpr bool plotSignal = true;
  constexpr bool plotField  = true;

   static constexpr double fGasGapThickness = 0.2;         // cm - total gas gap thickness            2000 micron = 2 mm = 0.2 cm
    static constexpr double fGasGapCenterY = 0.0/10.0;           // cm - centered at Y=0                    y = 0 mm =0 cm
    static constexpr double fAnodeCathodeThickness = 0.02/10.0;  // cm - thickness of HV layers             20 micron = 0.02 mm = 0.002 cm
    static constexpr double fStripThickness = 5.00/10.0;         // cm - thickness of copper strips         5 mm = 0.5 cm
    static constexpr double fStripWidthX = 300.0/10;              // cm - readout strip width                30 mm = 3 cm
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
    //rpc->SetFastVolOptions(1,1,0);
    rpc->Initialise();


  const std::size_t nx = 50, ny = 50, ncont = 104;

  ViewField *contourView = nullptr;
  TCanvas *cContour = nullptr;
  cContour = new TCanvas("cContour","Contour",1200,1200);
  contourView = new ViewField();
  contourView->SetCanvas(cContour);
  contourView->SetComponent(rpc);
  contourView->SetPlane(0,-1,0,0,0.,0);
  contourView->SetArea(-18,-0.1,-18,18,0.1,18);
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
  const std::size_t nTimeBins = 2000;
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

  // ----------------------------------------------------------------------
  // "Template" AvalancheMicroscopic — kept purely for the SERIAL fallback
  // path and to record the configuration we want.  In the parallel path,
  // we construct one of these per thread (with the same configuration but
  // WITHOUT plotting / debugging / shared histograms) inside the parallel
  // region.  All original setters are preserved here so behaviour matches
  // your current code when RUN_PARALLEL=false.
  // ----------------------------------------------------------------------
  AvalancheMicroscopic aval(&sensor);
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
  //aval.SetRunModeOptions(MPRunMode::GPUExclusive, 0);

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

  // Set up Heed.
  TrackHeed track(&sensor);
  // Set the particle type and momentum [eV/c].
  track.SetParticle("Pion");
  track.SetMomentum(7.e9);
  track.CrossInactiveMedia(true);

  std::clock_t start = std::clock();

  // ----------------------------------------------------------------------
  // Drift-plotting: only enabled in the SERIAL run.  ViewDrift is not
  // thread-safe, so in the parallel path we leave drift plotting off.
  // ----------------------------------------------------------------------
  auto cD = new TCanvas("cD","Drift",1200,800);
  ViewDrift driftView;
  driftView.SetCanvas(cD);

  if (!RUN_PARALLEL) {
    aval.EnablePlotting(&driftView, 1000000);
    avalMCi.EnablePlotting(&driftView);
    avalMCin.EnablePlotting(&driftView);
    track.EnablePlotting(&driftView);
  } else {
    // Track plotting is fine because the track itself runs serially,
    // but the per-thread aval/avalMCi/avalMCin will not be wired to it.
    track.EnablePlotting(&driftView);
  }

  // Simulate a charged-particle track (always serial — TrackHeed is not
  // thread-safe and we only have a single primary track here anyway).
  track.NewTrack(0, y0-0.00001, 0, 0, 0, -1, 0);

  // ----------------------------------------------------------------------
  // Flatten the cluster -> electron iteration into a vector of primaries.
  // This lets us drive the avalanche loop with #pragma omp for.
  // ----------------------------------------------------------------------
  struct Primary {
    double x, y, z, t, e;
  };
  std::vector<Primary> primaries;
  std::size_t nCluster = 0;
  for (const auto& cluster : track.GetClusters()) {
    ++nCluster;
    for (const auto& el : cluster.electrons) {
      primaries.push_back({el.x, el.y, el.z, el.t, el.e});
    }
  }
  const int nPrim = static_cast<int>(primaries.size());

  // Retrieve the clusters along the track.
  std::size_t nTotal = 0;
  std::size_t nIonMC = 0;
  std::size_t nNegIonMC = 0;
  std::size_t nAvalanche = 0;
  std::size_t nEMicro = 0;

  std::size_t neMicroAval = 0;
  std::size_t nIMicroAval = 0;
  std::size_t nDriftLine  = 0;

  // Reduction targets (size_t can't go in OpenMP reduction directly on
  // older compilers, so use long long).
  long long nEMicro_red    = 0;
  long long nIonMC_red     = 0;
  long long nNegIonMC_red  = 0;
  long long neMicroAval_red= 0;
  long long nIMicroAval_red= 0;

  // ======================================================================
  //                           PARALLEL REGION
  // ======================================================================
  #pragma omp parallel if(RUN_PARALLEL) \
      reduction(+:nEMicro_red,nIonMC_red,nNegIonMC_red, \
                  neMicroAval_red,nIMicroAval_red)
  {
    // Decide whether *this thread* uses the shared template object (serial
    // path) or builds its own thread-local transport objects (parallel
    // path).  In the parallel branch we mirror every original setter EXCEPT
    // the ones that touch shared state (energy histograms, drift plotting,
    // debug printing, progress reporting).
    const bool tUseShared = !RUN_PARALLEL;

    // ---------- thread-local transport objects (parallel path) ----------
    AvalancheMicroscopic  avalT_obj(&sensor);
    AvalancheMC           avalMCi_obj;
    AvalancheMC           avalMCin_obj;

    // Per-thread copies of the energy histograms — merged at the end.
    // Names must be unique per thread, otherwise ROOT's global object
    // registry warns "Replacing existing TH1 ... (Potential memory leak)".
    const int  tid       = omp_get_thread_num();
    const std::string hPrimName = "hPrim_t_thr" + std::to_string(tid);
    const std::string hSecName  = "hSec_t_thr"  + std::to_string(tid);
    TH1D hPrim_t(hPrimName.c_str(), "Per-thread primary E",  200, 0, 50);
    TH1D hSec_t (hSecName .c_str(), "Per-thread secondary E",200, 0, 50);
    hPrim_t.SetDirectory(nullptr);
    hSec_t .SetDirectory(nullptr);

    if (!tUseShared) {
      avalT_obj.SetTimeWindow(0., tMaxWindow);
      avalT_obj.EnableSignalCalculation(true);
      avalT_obj.UseInducedCharge(true);
      avalT_obj.UseWeightingPotential(true);
      avalT_obj.EnableElectronEnergyHistogramming(&hPrim_t);
      avalT_obj.EnableSecondaryEnergyHistogramming(&hSec_t);
      avalT_obj.EnableDriftLines(true);
      avalT_obj.EnableNullCollisionSteps();
      // NO EnableDebugging() / SetShowProgress() / EnablePlotting() in parallel.

      avalMCi_obj.SetTimeSteps(tstep);
      avalMCi_obj.SetSensor(&sensor);
      avalMCi_obj.SetTimeWindow(0., tMaxWindow);
      avalMCi_obj.EnableAvalancheSizeLimit(100000);
      avalMCi_obj.EnableSignalCalculation(true);
      avalMCi_obj.EnableDiffusion(true);
      avalMCi_obj.EnableAttachment(true);
      avalMCi_obj.EnableRecombination(true);
      // NO EnableDebugging() in parallel.

      avalMCin_obj.SetTimeSteps(tstep);
      avalMCin_obj.SetSensor(&sensor);
      avalMCin_obj.SetTimeWindow(0., tMaxWindow);
      avalMCin_obj.EnableAvalancheSizeLimit(100000);
      avalMCin_obj.EnableSignalCalculation(true);
      avalMCin_obj.EnableDiffusion(true);
      avalMCin_obj.EnableAttachment(true);
      avalMCin_obj.EnableRecombination(true);
      // NO EnableDebugging() in parallel.
    }

    // Pointers used in the loop body — serial path uses the original
    // shared aval/avalMCi/avalMCin (preserving every feature including
    // drift-line plotting and debug output).
    AvalancheMicroscopic& avalT   = tUseShared ? aval    : avalT_obj;
    AvalancheMC&          avalMCiT  = tUseShared ? avalMCi  : avalMCi_obj;
    AvalancheMC&          avalMCinT = tUseShared ? avalMCin : avalMCin_obj;

    // Thread-local TTree row buffers — avoid clobbering the global ones
    // until we're ready to Fill() under a critical section.
    double xp0_t, yp0_t, zp0_t, tp0_t, ep0_t;
    double xp1_t, yp1_t, zp1_t, tp1_t, ep1_t;
    int    statusPrim_t;
    std::vector<double> xs0_t, ys0_t, zs0_t, ts0_t, es0_t;
    std::vector<double> xs1_t, ys1_t, zs1_t, ts1_t, es1_t;
    std::vector<int>    statusSec_t;
    double bAlpha_t, bEta_t, bAlphaEff_t;
    double bVdElec_t, bVspeedElec_t, bMuElec_t;
    int    bNElecVel_t;
    int    nElectrons_t, nSecondaries_t;
    int    eventID_t;

    // ----------------------------------------------------------------
    //                 PARALLEL FOR over primary electrons
    // ----------------------------------------------------------------
    #pragma omp for schedule(dynamic)
    for (int ip = 0; ip < nPrim; ++ip) {
      const auto& Pe = primaries[ip];

      ++nEMicro_red;

      xs0_t.clear(); ys0_t.clear(); zs0_t.clear(); ts0_t.clear(); es0_t.clear();
      xs1_t.clear(); ys1_t.clear(); zs1_t.clear(); ts1_t.clear(); es1_t.clear();
      statusSec_t.clear();

      xp0_t = yp0_t = zp0_t = tp0_t = ep0_t = 0.;
      xp1_t = yp1_t = zp1_t = tp1_t = ep1_t = 0.;
      statusPrim_t = 0;

      // ---- Initial electron (kept verbatim from your original code) ----
      const double xi = 0.05, yi = 0.0, zi = 0.0;
      const double ti = 0.0;
      const double ei = 0.1;
      const double dx0 = 0.0, dy0 = 0.0, dz0 = 0.0;
      (void)xi; (void)yi; (void)zi; (void)ti; (void)ei;
      (void)dx0; (void)dy0; (void)dz0;

      // Run the microscopic avalanche for this primary.
      avalT.AvalancheElectron(Pe.x, Pe.y, Pe.z, Pe.t, Pe.e, 0., 0., 0.);
      avalMCiT.DriftIon(Pe.x, Pe.y, Pe.z, Pe.t);

      int ne = 0, ni = 0;
      avalT.GetAvalancheSize(ne, ni);
      neMicroAval_red += ne;
      nIMicroAval_red += ni;

      nElectrons_t   = avalT.GetNumberOfElectronEndpoints();
      //nElectrons_t = avalT.GetNumberOfElectronEndpointsGPU();
      nSecondaries_t = (nElectrons_t > 0) ? nElectrons_t - 1 : 0;
      eventID_t      = ip;

      // ---- Split endpoints into primary (j=0) and secondaries (j>0) ----
      for (int j = 0; j < nElectrons_t; ++j) {
        double xs, ys, zs, ts, es;
        double xe, ye, ze, te, ee;
        int st;
        //avalT.GetElectronEndpointGPU(j, xs, ys, zs, ts, es, xe, ye, ze, te, ee, st);
        avalT.GetElectronEndpoint(j,
                                  xs, ys, zs, ts, es,
                                  xe, ye, ze, te, ee,
                                  st);
        if (j == 0) {
          // ---- PRIMARY ----
          xp0_t = xs; yp0_t = ys; zp0_t = zs; tp0_t = ts; ep0_t = es;
          xp1_t = xe; yp1_t = ye; zp1_t = ze; tp1_t = te; ep1_t = ee;
          statusPrim_t = st;
        } else {
          // ---- SECONDARY ----
          xs0_t.push_back(xs); ys0_t.push_back(ys); zs0_t.push_back(zs);
          ts0_t.push_back(ts); es0_t.push_back(es);
          xs1_t.push_back(xe); ys1_t.push_back(ye); zs1_t.push_back(ze);
          ts1_t.push_back(te); es1_t.push_back(ee);
          statusSec_t.push_back(st);
        }
      }

      // ---- Drift positive ions ----
      // AvalancheMicroscopic has no GetIons(); each secondary electron
      // (index > 0) was born from an ionization event, so path[0] of
      // each secondary electron is the ion creation point.
      double Ltot = 0.;
      std::size_t nNegIons = 0;
      const auto& allElec = avalT.GetElectrons();

      int nIo = 0, nEl = 0;
      avalT.GetAvalancheSize(nEl, nIo);
      if (!allElec.empty()) Ltot += allElec[0].pathLength;

      for (std::size_t ie = 1; ie < allElec.size(); ++ie) {
        Ltot += allElec[ie].pathLength;
        if (allElec[ie].path.empty()) continue;
        const auto& p0 = allElec[ie].path[0];
        avalMCiT.DriftIon(p0.x, p0.y, p0.z, p0.t);
        ++nIonMC_red;
      }

      // ---- Drift negative ions (SF6-) ----
      // StatusAttached == -7 in Garfield++; verify against your version if needed.
      for (const auto& elec : allElec) {
        if (elec.status == -7) {
          ++nNegIons;
          if (!elec.path.empty()) {
            const auto& pEnd = elec.path.back();
            avalMCinT.DriftNegativeIon(pEnd.x, pEnd.y, pEnd.z, pEnd.t);
            ++nNegIonMC_red;
          }
        }
      }

      // ---- Field at evaluation point + drift-direction unit vector ----
      double xMid = 0., yMid = 0., zMid = 0.;
      double Ex, Ey, Ez; Medium* m = nullptr; int sst = 0;
      sensor.ElectricField(xMid, yMid, zMid, Ex, Ey, Ez, m, sst);
      const double Emag = std::sqrt(Ex*Ex + Ey*Ey + Ez*Ez);
      const double ux = (Emag > 0.) ? -Ex/Emag : 0.;
      const double uy = (Emag > 0.) ? -Ey/Emag : 0.;
      const double uz = (Emag > 0.) ? -Ez/Emag : 0.;

      // Original code computed v_drift from the primary then overwrote
      // via the loop below.  We keep the (kept-but-unused) variables for
      // parity with the original logic; the real numbers come from the
      // averaging loop.
      if (!allElec.empty() && allElec[0].path.size() >= 2) {
        const auto& a = allElec[0].path.front();
        const auto& b = allElec[0].path.back();
        const double dt0 = b.t - a.t;
        const double dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
        if (dt0 > 0.) {
          const double v_drift_primary = (dx*ux + dy*uy + dz*uz) / dt0;
          (void)v_drift_primary;
        }
      }

      // ---- Average drift / speed over all electrons in this avalanche ----
      double sumVd = 0., sumVspeed = 0.; int nE = 0;
      for (const auto& e : allElec) {
        if (e.path.size() < 2) continue;
        const auto& p0 = e.path.front();
        const auto& p1 = e.path.back();
        const double dt = p1.t - p0.t;
        if (dt <= 0.) continue;

        const double dx = p1.x - p0.x;
        const double dy = p1.y - p0.y;
        const double dz = p1.z - p0.z;

        const double v_drift = (dx*ux + dy*uy + dz*uz) / dt;
        const double v_speed = e.pathLength / dt;

        sumVd     += v_drift;
        sumVspeed += v_speed;
        ++nE;
      }
      const double vdElec_evt     = (nE > 0) ? sumVd     / nE : 0.;
      const double vspeedElec_evt = (nE > 0) ? sumVspeed / nE : 0.;
      const double muElec = (Emag > 0.) ? std::abs(vdElec_evt) / Emag : 0.;

      bVdElec_t     = vdElec_evt;
      bVspeedElec_t = vspeedElec_evt;
      bMuElec_t     = muElec;
      bNElecVel_t   = nE;

      // ---- alpha / eta / alpha_eff ----
      const double alpha     = (Ltot > 0.) ? double(nIo)      / Ltot : 0.;
      const double eta       = (Ltot > 0.) ? double(nNegIons) / Ltot : 0.;
      const double alpha_eff = alpha - eta;
      bAlpha_t    = alpha;
      bEta_t      = eta;
      bAlphaEff_t = alpha_eff;

      // ---- Critical section: write into shared TTree + histograms ----
      #pragma omp critical(treefill)
      {
        // Push thread-local row into the global branch buffers, then Fill().
        eventID      = eventID_t;
        nElectrons   = nElectrons_t;
        nSecondaries = nSecondaries_t;

        xp0 = xp0_t; yp0 = yp0_t; zp0 = zp0_t; tp0 = tp0_t; ep0 = ep0_t;
        xp1 = xp1_t; yp1 = yp1_t; zp1 = zp1_t; tp1 = tp1_t; ep1 = ep1_t;
        statusPrim = statusPrim_t;

        xs0 = xs0_t; ys0 = ys0_t; zs0 = zs0_t; ts0 = ts0_t; es0 = es0_t;
        xs1 = xs1_t; ys1 = ys1_t; zs1 = zs1_t; ts1 = ts1_t; es1 = es1_t;
        statusSec = statusSec_t;

        bAlpha    = bAlpha_t;
        bEta      = bEta_t;
        bAlphaEff = bAlphaEff_t;
        bVdElec     = bVdElec_t;
        bVspeedElec = bVspeedElec_t;
        bMuElec     = bMuElec_t;
        bNElecVel   = bNElecVel_t;

        fieldtree->Fill();

        if (Ltot > 0.) {
          hAlpha   ->Fill(alpha);
          hEta     ->Fill(eta);
          hAlphaEff->Fill(alpha_eff);
        }
        if (ne > 0) hLogGain->Fill(std::log10(double(ne)));
      }
    } // end omp for

    // ---- Merge per-thread energy histograms into the global ones ----
    if (!tUseShared) {
      #pragma omp critical(merge_hist)
      {
        ePrimaryEnergy ->Add(&hPrim_t);
        eSecondaryEnergy->Add(&hSec_t);
      }
    }
  } // end omp parallel

  // Promote reduction totals back into the original size_t variables.
  nEMicro     = static_cast<std::size_t>(nEMicro_red);
  nIonMC      = static_cast<std::size_t>(nIonMC_red);
  nNegIonMC   = static_cast<std::size_t>(nNegIonMC_red);
  neMicroAval = static_cast<std::size_t>(neMicroAval_red);
  nIMicroAval = static_cast<std::size_t>(nIMicroAval_red);

  // ----------------------------------------------------------------------
  // Write tree + histograms into the open TFile.
  // ----------------------------------------------------------------------
  fieldtree->Write();
  hAlpha   ->Write();
  hEta     ->Write();
  hAlphaEff->Write();
  hLogGain ->Write();

  // ----------------------------------------------------------------------
  // Drift-line plot — only meaningful for the serial run, but we always
  // produce the canvas to keep the output set the same.
  // ----------------------------------------------------------------------
  driftView.SetColourElectrons(kBlue);
  driftView.SetColourIons(kRed);
  driftView.SetColourNegativeIons(kGreen);
  driftView.Plot();
  cD->Update();
  cD->Write();
  cD->SaveAs("../result_cache/png/DriftPlot.png");
  nDriftLine = driftView.GetNumberOfDriftLines();

  if (plotSignal) {
    // Plot the induced current.
    //signalView->PlotSignal(label1);
    //signalView->EnableLegend();
    //cSignal->Update();
    gSystem->ProcessEvents();
    sensor.ExportSignal(label1, "../result_cache/sensor_out/Signal_"+label1);
    sensor.ExportSignal(label2, "../result_cache/sensor_out/Signal_"+label2);
    sensor.ExportSignal(label3, "../result_cache/sensor_out/Signal_"+label3);
    sensor.ExportSignal(label4, "../result_cache/sensor_out/Signal_"+label4);

    treefile->cd();
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
          Timebin.push_back(i*nTimeBins);
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
      (void)pTimebin;
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
    std::ofstream fcsv("../result_cache/signal_" + lbl + ".txt");
    fcsv << "# RPC induced signal at " << lbl << "\n";
    fcsv << "# TimeBin\tTime[ns]\tElectronSignal\tIonSignal\tTotalSignal\n";
    fcsv << std::scientific << std::setprecision(6);

    for (unsigned int i = 0; i < NumberBins; ++i) {
        const double time   = TimeStart + (i + 0.5) * TimeStep;
        const double sElec  = sensor.GetElectronSignal(lbl, i);
        const double sIon   = sensor.GetIonSignal     (lbl, i);
        const double sTotal = sensor.GetSignal        (lbl, i);

        fcsv << i      << '\t'
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

  std::cout << "Total Cluster: " << nCluster << std::endl;
  std::cout << "++++ Aval Microscopic ++++"<< std::endl;
  std::cout << "Number of Electron (Avalanche):  " << neMicroAval << std::endl;
  std::cout << "Number of Ion (Avalanche):  " << nIMicroAval << std::endl;
  std::cout << std::endl;
  std::cout << "++++ Aval MC ++++"<< std::endl;
  std::cout << "Number of Positive Ion Drifts:  " << nIonMC << std::endl;
  std::cout << "Number of Negative Ion Drifts:  " << nNegIonMC << std::endl;
  std::cout << "Number of Primary Electron (Microscopic): " << nEMicro << std::endl;
  std::cout << "Time Bin Size: " << tstep << std::endl;
  std::cout << "Number of Drift Line: " << nDriftLine << std::endl;

  // (void) on unused vars from the original code, kept for parity.
  (void)debug; (void)plotField; (void)nx; (void)ny; (void)ncont;
  (void)nX; (void)nZ; (void)dY; (void)nY; (void)y0;
  (void)gasGapTop; (void)gasGapBottom;
  (void)nTotal; (void)nAvalanche; (void)start;

  //std::cout<< path+"share/Garfield/Data/IonMobility_SF6-_SF6.txt" <<std::endl;
  LOG("End of Program");
  app.Run(true);

  return 0;
}