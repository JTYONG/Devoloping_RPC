#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TFolder.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>

using namespace std;

// Forward declaration moved outside the function
void StartPrint();

void Read_NebemPlot(const char* argv) { // FIX 1: char -> const char*

    StartPrint(); // actually call it

    // FIX 4: properly open and check the file
    ifstream file(argv);
    if (!file.is_open()) {
        cerr << "!! File open error: " << argv << endl;
        return;
    }
    file.close(); // TChain will handle reading

    // FIX 2: use std::string for concatenation
    string filename = argv;
    string outfilename = "ROOTPlot_" + filename;

    struct fielddata {
        int pointindex;
        double x, y, z;
        double ex, ey, ez, emag;
        int status;
    };

    fielddata NebemField;

    TChain *NebemFile = new TChain("FieldTree");
    NebemFile->Add(filename.c_str());

    NebemFile->SetBranchAddress("point index", &NebemField.pointindex);
    NebemFile->SetBranchAddress("x",           &NebemField.x);
    NebemFile->SetBranchAddress("y",           &NebemField.y);
    NebemFile->SetBranchAddress("z",           &NebemField.z);
    NebemFile->SetBranchAddress("ex",          &NebemField.ex);
    NebemFile->SetBranchAddress("ey",          &NebemField.ey);
    NebemFile->SetBranchAddress("ez",          &NebemField.ez);
    NebemFile->SetBranchAddress("emag",        &NebemField.emag);
    NebemFile->SetBranchAddress("status",      &NebemField.status);

    int NofPoints = (int)NebemFile->GetEntries();

    // FIX 2: use std::string variable for TFile name
    TFile *treefile = new TFile(outfilename.c_str(), "recreate");
    cout << "Creating ROOTPlot File: " << outfilename << " ....." << endl << endl;

    // FIX 3: add 'new' to create a pointer
    TH2D *EyVSy_TH2D = new TH2D("EyVSy", "E_y (V/cm) vs y (cm)",
                                  100, 0, 300000, 100, -0.0125, 0.0125);
    EyVSy_TH2D->GetYaxis()->SetTitle("y coordinate in gas gap (cm)");
    EyVSy_TH2D->GetXaxis()->SetTitle("Ey (V/cm)");
    EyVSy_TH2D->SetDrawOption("COLZ");

    TH3D *X_Z_Emag_TH3D = new TH3D("X_Z_Emag_TH3D", "E_mag (V/cm) across x,z (all sampled xz)",
                                   100, -1.5, 1.5,100, -1.5, 1.5,100, 0, 300000);
    X_Z_Emag_TH3D->GetYaxis()->SetTitle("z coordinate in gas gap (cm)");
    X_Z_Emag_TH3D->GetXaxis()->SetTitle("x coordinate in gas gap (cm)");
    X_Z_Emag_TH3D->GetZaxis()->SetTitle("Emag (V/cm)");
    X_Z_Emag_TH3D->SetDrawOption("GLBOX");

    TH1D *Ey_profile = new TH1D("Ey_profile", "E_y profile across gas gap",
                              100, -0.0125, 0.0125);
    TProfile *Ey_prof = new TProfile("Ey_prof", "TProfile: Mean E_y vs y",
                                  100, -0.0125, 0.0125, 0, 300000);
    Ey_prof->GetXaxis()->SetTitle("y (cm)");
    Ey_prof->GetYaxis()->SetTitle("Mean E_{y} (V/cm)");
    Ey_prof->SetMarkerStyle(20);     // filled circle (most visible)
    Ey_prof->SetMarkerSize(0.8);     // marker size
    Ey_prof->SetMarkerColor(kBlue);
    Ey_prof->SetLineColor(kBlue);
    Ey_prof->SetLineWidth(2);        // thicker error bar lines
    // Don't hardcode the range, let ROOT find it
Ey_prof->SetMinimum(0);  // remove manual range first
Ey_prof->SetMaximum(0);  // 0 = auto

// Then zoom based on actual mean and RMS
double mean  = Ey_prof->GetMean(2);   // mean of Y values
double sigma = Ey_prof->GetRMS(2);    // RMS of Y values

Ey_prof->SetMinimum(mean - 3*sigma);
Ey_prof->SetMaximum(mean + 3*sigma);

Ey_prof->Draw("E1");             // E1 = error bars with caps
    for (int i = 0; i < NofPoints; i++) { // NofPoints+1 would overflow
        if (i % 1000 == 0) {
            cout << "Loading field point: " << i << endl;
        }
        NebemFile->GetEntry(i);
        // Fill the histogram properly: (x-value, y-value)
        EyVSy_TH2D->Fill(abs(NebemField.ey), NebemField.y);
        X_Z_Emag_TH3D->Fill(NebemField.x, NebemField.z, NebemField.emag);

        // Fill with:
        Ey_profile->Fill(NebemField.y, abs(NebemField.ey));

        Ey_prof->Fill(NebemField.y, abs(NebemField.ey));
        // Increase cap size (default is 0.03)
    gStyle->SetEndErrorSize(5);   // bigger caps — call before Draw()

    Ey_prof->Draw("E1");

    }
cout << "Mean Ey:  " << Ey_prof->GetMean(2)  << endl;
cout << "RMS Ey:   " << Ey_prof->GetRMS(2)   << endl;
cout << "Min bin:  " << Ey_prof->GetMinimum() << endl;
cout << "Max bin:  " << Ey_prof->GetMaximum() << endl;
    treefile->Write();
    treefile->Close();
}

void StartPrint() {
    cout << "*******************************************************************" << endl;
    cout << "*****              Visualization of Nebem Field               *****" << endl;
    cout << "*******************************************************************" << endl;
    cout << endl;
}