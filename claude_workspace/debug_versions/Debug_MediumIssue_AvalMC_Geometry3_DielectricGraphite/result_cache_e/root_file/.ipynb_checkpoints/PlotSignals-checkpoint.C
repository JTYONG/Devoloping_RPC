// plotSignals.C
// Read the Electron, Ion, and Total signal vectors from one entry of a tree
// in signal.root and overlay them on a single canvas.
//
// Usage:   root -l 'plotSignals.C("TopReadout")'
//          root -l plotSignals.C            // defaults to "TopReadout"
//
// The file is assumed to have been produced with vector<double> branches
// named "Time Bin", "Electron Signal", "Ion Signal", "Total Signal"
// (the Garfield++ default layout: one waveform per entry).

void PlotSignals(const char* treeName = "BottomReadout",
                 const char* fileName = "signal.root",
                 Long64_t    entry    = 0) {

    // --- Open file and grab the tree ------------------------------------
    TFile *f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
        std::cerr << "Could not open " << fileName << std::endl;
        return;
    }

    TTree *t = (TTree*)f->Get(treeName);
    if (!t) {
        std::cerr << "Tree '" << treeName << "' not found in " << fileName << std::endl;
        return;
    }

    // --- Hook up vector branches ---------------------------------------
    // Pointers MUST be initialised to nullptr; ROOT allocates the vectors.
    std::vector<double> *tbin   = nullptr;
    std::vector<double> *sElec  = nullptr;
    std::vector<double> *sIon   = nullptr;
    std::vector<double> *sTotal = nullptr;

    t->SetBranchAddress("Time Bin",        &tbin);
    t->SetBranchAddress("Electron Signal", &sElec);
    t->SetBranchAddress("Ion Signal",      &sIon);
    t->SetBranchAddress("Total Signal",    &sTotal);

    if (entry >= t->GetEntries()) {
        std::cerr << "Entry " << entry << " out of range (tree has "
                  << t->GetEntries() << " entries)" << std::endl;
        return;
    }
    t->GetEntry(entry);

    const Int_t n = static_cast<Int_t>(tbin->size());
    std::cout << "Loaded " << n << " samples from tree '" << treeName
              << "', entry " << entry << std::endl;

    // --- Build TGraphs --------------------------------------------------
    TGraph *gElec  = new TGraph(n, tbin->data(), sElec->data());
    TGraph *gIon   = new TGraph(n, tbin->data(), sIon->data());
    TGraph *gTotal = new TGraph(n, tbin->data(), sTotal->data());

    gElec ->SetLineColor(kBlue);
    gIon  ->SetLineColor(kRed);
    gTotal->SetLineColor(kBlack);

    gElec ->SetLineWidth(2);
    gIon  ->SetLineWidth(2);
    gTotal->SetLineWidth(2);

    gElec ->SetTitle("Electron signal");
    gIon  ->SetTitle("Ion signal");
    gTotal->SetTitle("Total signal");

    // --- Draw on a single canvas using TMultiGraph ----------------------
    TCanvas *c = new TCanvas("c", treeName, 1680, 1050);
    c->SetGrid();

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(Form("Induced signals at %s (fC/ns);Time bin;Signal", treeName));
    mg->Add(gElec,  "L");
    mg->Add(gIon,   "L");
    mg->Add(gTotal, "L");
    mg->Draw("AL");

    TLegend *leg = new TLegend(0.70, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    //leg->AddEntry(gElec,  "Electron", "l");
    //leg->AddEntry(gIon,   "Ion",      "l");
    leg->AddEntry(gTotal, "Total",    "l");
     leg->AddEntry(gElec,  "Electron", "l");
     leg->AddEntry(gIon,   "Ion",      "l");
    leg->Draw();

    c->Update();

    // Optionally save to file:
    c->SaveAs(Form("../png/Root_signals_%s.png", treeName));
}
