#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>

void plot_signals() {
    // Open ROOT file
    TFile* file = TFile::Open("/home/sus/New/Garfieldpp/Examples/RPC/build/rpc_output.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the tree
    TTree* tree = (TTree*)file->Get("SignalTree");
    if (!tree) {
        std::cerr << "Error: SignalTree not found!" << std::endl;
        file->Close();
        return;
    }

    // Set up variables
    double time, signal1, signal2;
    tree->SetBranchAddress("time", &time);
    tree->SetBranchAddress("signal1", &signal1);
    tree->SetBranchAddress("signal2", &signal2);

    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", "RPC Signals with Custom Range", 1000, 600);
    
    // First pass to find data range
    double sig1_min = DBL_MAX, sig1_max = -DBL_MAX;
    double sig2_min = DBL_MAX, sig2_max = -DBL_MAX;
    double time_min = DBL_MAX, time_max = -DBL_MAX;
    
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        sig1_min = TMath::Min(sig1_min, signal1);
        sig1_max = TMath::Max(sig1_max, signal1);
        sig2_min = TMath::Min(sig2_min, signal2);
        sig2_max = TMath::Max(sig2_max, signal2);
        time_min = TMath::Min(time_min, time);
        time_max = TMath::Max(time_max, time);
    }

    // Add small margin (5% of range)
    const double margin_factor = 0.05;
    double sig_range = sig1_max - sig1_min;
    sig1_min -= margin_factor * sig_range;
    sig1_max += margin_factor * sig_range;
    
    sig_range = sig2_max - sig2_min;
    sig2_min -= margin_factor * sig_range;
    sig2_max += margin_factor * sig_range;

    // Create graphs with custom range
    TGraph* graph1 = new TGraph(nEntries);
    TGraph* graph2 = new TGraph(nEntries);
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        graph1->SetPoint(i, time, signal1);
        graph2->SetPoint(i, time, signal2);
    }

    // Configure graphs
    graph1->SetTitle("RPC Signals with Custom Range;Time [ns];Signal [fC/ns]");
    graph1->GetXaxis()->SetLimits(time_min, time_max);
    graph1->GetYaxis()->SetRangeUser(sig1_min-0.4, sig1_max);
    graph1->SetLineColor(kRed);
    graph1->SetMarkerColor(kRed);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(0.5);
    
    graph2->SetLineColor(kBlue);
    graph2->SetMarkerColor(kBlue);
    graph2->SetMarkerStyle(21);
    graph2->SetMarkerSize(0.5);

    // Draw graphs
    graph1->Draw("APL");
    graph2->Draw("PL same");

    // Add legend
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(graph1, "Electrode 1", "lp");
    leg->AddEntry(graph2, "Electrode 2", "lp");
    leg->Draw();

    canvas->Update();
    file->Close();
}

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    plot_signals();
    app.Run();
    return 0;
}
