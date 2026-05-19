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
cD.Update()
cD.SaveAs('plots/PlotElectronDiffusion.png')


cT = ROOT.TCanvas('cT', '', 600, 600)
view.SetCanvas(cT)
view.PlotElectronTownsend()
cT.Update()
cT.SaveAs('plots/PlotElectronTownsend.png')

cA = ROOT.TCanvas('cA', '', 600, 600)
view.SetCanvas(cA)
view.PlotElectronAttachment()
cA.Update()
cA.SaveAs('plots/PlotElectronAttachment.png')

