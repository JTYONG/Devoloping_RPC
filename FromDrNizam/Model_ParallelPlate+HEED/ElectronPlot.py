import ROOT
import csv

time       = []
i_electron = []
#i_ion      = []
#i_total    = []

with open("ElectronIonCurrent.csv") as f:
    reader = csv.reader(f)
    next(reader)  # skip header line
    next(reader)  # skip "The induced signal." line
    for row in reader:
        row = [r.strip() for r in row]
        if len(row) < 10:
            continue
        try:
            time.append(float(row[0]))
 #           i_total.append(float(row[3]))
            i_electron.append(float(row[6]))
  #          i_ion.append(float(row[9]))
        except ValueError:
            continue

n = len(time)

g_electron = ROOT.TGraph(n)
#g_ion      = ROOT.TGraph(n)
#g_total    = ROOT.TGraph(n)

for i in range(n):
    g_electron.SetPoint(i, time[i], i_electron[i])
 #   g_ion.SetPoint(i,      time[i], i_ion[i])
  #  g_total.SetPoint(i,    time[i], i_total[i])

# Style

g_electron.SetLineColor(ROOT.kBlue);  g_electron.SetLineWidth(2)
g_electron.SetMarkerColor(ROOT.kBlue); g_electron.SetMarkerStyle(20); g_electron.SetMarkerSize(0.6)
'''
g_ion.SetLineColor(ROOT.kRed);  g_ion.SetLineWidth(2)
g_ion.SetMarkerColor(ROOT.kRed); g_ion.SetMarkerStyle(21); g_ion.SetMarkerSize(0.6)

g_total.SetLineColor(ROOT.kBlack); g_total.SetLineWidth(2); g_total.SetLineStyle(2)
'''
c = ROOT.TCanvas("c_signal", "Induced Signal", 900, 600)
c.SetGrid(); c.SetLeftMargin(0.12); c.SetBottomMargin(0.12)

g_electron.SetTitle("Induced Current;Time [ns];Current [fC/ns]")
g_electron.Draw("APL")
g_electron.GetXaxis().SetRangeUser(time[0], 70.0)
'''
g_ion.Draw("APL")
g_ion.GetXaxis().SetRangeUser(time[0], 70.0)
#g_ion.Draw("PL SAME")
#g_total.Draw("L SAME")
'''
#leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
#leg.SetBorderSize(1)
#leg.AddEntry(g_electron, "Electron signal", "lp")
#leg.AddEntry(g_ion,      "Ion signal",      "lp")
#leg.AddEntry(g_total,    "Total signal",    "l")
#leg.Draw()

ROOT.gPad.Update()
c.SaveAs("induced_Electron_signal.png")
print("Saved induced_Electron_signal.png")
