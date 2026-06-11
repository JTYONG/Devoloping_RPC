import ROOT
import csv

# Read CSV: columns are time, i_electron, i_ion_pos, i_total
# (col3 = 0 always, col4 = i_total = col2 in this file)
time       = []
i_electron = []
i_ion      = []
i_total    = []

with open("ElectronIonCurrent.csv") as f:
    for row in csv.reader(f):
        row = [r.strip() for r in row]
        if len(row) < 4:
            continue
        try:
            time.append(float(row[0]))
            i_electron.append(float(row[1]))
            i_ion.append(float(row[2]))
            i_total.append(float(row[3]))
        except ValueError:
            continue  # skip header lines if any

n = len(time)

# Build TGraphs
g_electron = ROOT.TGraph(n)
g_ion      = ROOT.TGraph(n)
g_total    = ROOT.TGraph(n)

for i in range(n):
    g_electron.SetPoint(i, time[i], i_electron[i])
    g_ion.SetPoint(i,      time[i], i_ion[i])
    g_total.SetPoint(i,    time[i], i_total[i])

# Style
g_electron.SetLineColor(ROOT.kBlue)
g_electron.SetLineWidth(2)
g_electron.SetMarkerColor(ROOT.kBlue)
g_electron.SetMarkerStyle(20)
g_electron.SetMarkerSize(0.6)

g_ion.SetLineColor(ROOT.kRed)
g_ion.SetLineWidth(2)
g_ion.SetMarkerColor(ROOT.kRed)
g_ion.SetMarkerStyle(21)
g_ion.SetMarkerSize(0.6)

g_total.SetLineColor(ROOT.kBlack)
g_total.SetLineWidth(2)
g_total.SetLineStyle(2)

# Canvas
c = ROOT.TCanvas("c_signal", "Induced Signal", 900, 600)
c.SetGrid()
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.12)

# Draw electron signal first (sets axis range)
g_electron.SetTitle("Induced Current;Time [ns];Current [fC/ns.]")
g_electron.Draw("APL")
g_electron.GetXaxis().SetRangeUser(time[0], 70.0)
g_ion.Draw("PL SAME")
g_total.Draw("L SAME")

# Legend
leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
leg.SetBorderSize(1)
leg.AddEntry(g_electron, "Electron component", "lp")
leg.AddEntry(g_ion,      "Ion component",      "lp")
leg.AddEntry(g_total,    "Total",              "l")
leg.Draw()

ROOT.gPad.Update()
c.SaveAs("induced_signal.png")
print("Saved induced_signal.png")
