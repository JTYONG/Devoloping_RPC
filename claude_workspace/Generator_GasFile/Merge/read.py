import ROOT
import Garfield

# Set up the gas.
gas = ROOT.Garfield.MediumMagboltz()
gas.LoadGasFile('build/TIFRH_merged_40000_43000.gas')
#gas.LoadIonMobility('IonMobility_Ar+_Ar.txt')
gas.PrintGas()

view = ROOT.Garfield.ViewMedium(gas)
view.SetMagneticField(2.)
  
cV = ROOT.TCanvas('cV', '', 600, 600)
view.SetCanvas(cV)
view.PlotElectronVelocity()
cV.Update()
cV.SaveAs('plots/PlotElectronVelocity.png')

cD = ROOT.TCanvas('cD', '', 600, 600)
view.SetCanvas(cD)
view.PlotElectronDiffusion()
cV.SaveAs('plots/PlotElectronDiffusion.png')
cV.Update()

cT = ROOT.TCanvas('cT', '', 600, 600)
view.SetCanvas(cT)
view.PlotElectronTownsend()
cV.SaveAs('plots/PlotElectronTownsend.png')
cV.Update()

cA = ROOT.TCanvas('cA', '', 600, 600)
view.SetCanvas(cA)
view.PlotElectronAttachment()
cV.SaveAs('plots/PlotElectronAttachment.png')
cV.Update()
