import ROOT
import Garfield

# Set up the gas.
gas = ROOT.Garfield.MediumMagboltz()
gas.LoadGasFile('rpc_95.5_4.2_0.3.gas')
#gas.LoadIonMobility('IonMobility_Ar+_Ar.txt')
gas.PrintGas()

view = ROOT.Garfield.ViewMedium(gas)
view.SetMagneticField(2.)
  
cV = ROOT.TCanvas('cV', '', 600, 600)
view.SetCanvas(cV)
view.PlotElectronVelocity()
cV.SaveAs('plots/PlotElectronVelocity.png')

cD = ROOT.TCanvas('cD', '', 600, 600)
view.SetCanvas(cD)
view.PlotElectronDiffusion()
cV.SaveAs('plots/PlotElectronDiffusion.png')

cT = ROOT.TCanvas('cT', '', 600, 600)
view.SetCanvas(cT)
view.PlotElectronTownsend()
cV.SaveAs('plots/PlotElectronTownsend.png')


cA = ROOT.TCanvas('cA', '', 600, 600)
view.SetCanvas(cA)
view.PlotElectronAttachment()
cV.SaveAs('plots/PlotElectronAttachment.png')
