import ROOT
import Garfield

# Set up the gas.
gas = ROOT.Garfield.MediumMagboltz()
gas.LoadGasFile('rpc_95.5_4.2_0.3_40-45_ncoll=20.gas')
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

cRN = ROOT.TCanvas('cRN', '', 600, 600)
view.SetCanvas(cRN)
view.PlotElectronReducedTownsendN()
cRN.Update()
cRN.SaveAs('plots/PlotElectronReducedTownsendN.png')

cRP = ROOT.TCanvas('cRP', '', 600, 600)
view.SetCanvas(cRP)
view.PlotElectronReducedTownsendP()
cRP.Update()
cRP.SaveAs('plots/PlotElectronReducedTownsendP.png')

cTOFA = ROOT.TCanvas('cT', '', 600, 600)
view.SetCanvas(cTOFA)
view.PlotElectronTOFAttachment()
cTOFA.Update()
cTOFA.SaveAs('plots/PlotElectronTOFAttachment.png')

cTOFI = ROOT.TCanvas('cT', '', 600, 600)
view.SetCanvas(cTOFI)
view.PlotElectronTOFIonization()
cTOFI.Update()
cTOFI.SaveAs('plots/PlotElectronTOFIonization.png')
