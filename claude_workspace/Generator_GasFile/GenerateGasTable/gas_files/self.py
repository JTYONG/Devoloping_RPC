import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gSystem.Load("libGarfield")
from ROOT import Garfield

gas = Garfield.MediumMagboltz()
gas.LoadGasFile("rpc_95.5_4.2_0.3.gas")

import ctypes
ef, bf, ba = ROOT.std.vector('double')(), ROOT.std.vector('double')(), ROOT.std.vector('double')()
gas.GetFieldGrid(ef, bf, ba)

V = []
alpha = []
eta = []
eff_alpha = []
DL, DT = [], []

ib = ia = 0

for ie in range(ef.size()):
    print(f"Grid Point: {ie}")
    e = gas.UnScaleElectricField(ef[ie])
    v = ctypes.c_double(0.)
    a1 = ctypes.c_double(0.)
    e1 = ctypes.c_double(0.)
    v1 = ctypes.c_double(0.)
    v2 = ctypes.c_double(0.)
    dl1 = ctypes.c_double(0.)
    dt1 = ctypes.c_double(0.)
    if gas.GetElectronVelocityE(ie, ib, ia, v):
        print(f"field: {e}")
        print(f"raw-value: {v.value}")
        print(f"velocity: {gas.ScaleVelocity(v.value)}")
        V.append(gas.ScaleVelocity(v.value))
    ok = gas.GetElectronTownsend(ie, ib, ia, a1)
    print(f"ie={ie} ok={ok} raw={a1.value}")
    if gas.GetElectronTownsend(ie,ib,ia,a1):
        print(f"raw-value: {a1.value-np.log(760)}")
        print(f"townsend: {gas.ScaleTownsend(np.exp(a1.value))}")
        alpha.append(np.exp(a1.value))
    if gas.GetElectronAttachment(ie,ib,ia,e1):
        print(f"raw-value: {e1.value-np.log(760)}")
        print(f"attachment: {gas.ScaleAttachment(np.exp(e1.value))}")
        eta.append(np.exp(e1.value))
    if gas.GetElectronAttachment(ie,ib,ia,v1) and gas.GetElectronTownsend(ie,ib,ia,v2):
        print(f"raw-value: {v2.value-np.log(760),v1.value-np.log(760)}") 
        print(f"effective townsend : {gas.ScaleTownsend(np.exp(v2.value))-gas.ScaleAttachment(np.exp(v1.value))}")
        eff_alpha.append(gas.ScaleTownsend(np.exp(v2.value))-gas.ScaleAttachment(np.exp(v1.value)))    
    if gas.GetElectronLongitudinalDiffusion(ie, ib, ia, dl1) and gas.GetElectronTransverseDiffusion(ie, ib, ia, dt1):
            print(f"raw-value (re-multiplied by sqrt(760)): " f"{dl1.value * np.sqrt(760), dt1.value * np.sqrt(760)}")
            print(f"longitudinal diffusion: {gas.ScaleDiffusion(dl1.value)}")
            print(f"transverse diffusion:   {gas.ScaleDiffusion(dt1.value)}")
            DL.append(gas.ScaleDiffusion(dl1.value))
            DT.append(gas.ScaleDiffusion(dt1.value))
    print("\n")

V1 = np.array(V)
alpha1 = np.array(alpha)
eta1 = np.array(eta)
eff_alpha1 = np.array(eff_alpha)
DL1 = np.array(DL)
DT1 = np.array(DT)

print(f"V: {V1}")
print(f"alpha: {alpha1}")
print(f"eta: {eta1}")
print(f"eff_alpha: {eff_alpha1}")
print(f"length of ef: {len(ef)}")
print(len(V1))
print(len(alpha1))
print(len(eta1))
print(len(eff_alpha1))

#ea = gas.UnScaleElectricField(ef[ie])

fig,ax = plt.subplots(2,2,figsize=(18,18))

ax[0,0].scatter(ef,V1,marker="x",color="blue")
ax[0,0].set_title("Electron Drift Velocity")
ax[0,0].set_xlabel(r"Efield [V/cm]")
ax[0,0].set_ylabel(r"$V_D$ [cm/ns]")

ax[0,1].scatter(ef,alpha1,marker="x",color="blue",label="Townsend")
ax[0,1].scatter(ef,eta1,marker="x",color="red",label="Attachment")
ax[1,1].scatter(ef,eff_alpha1,marker="x",color="black",label="Effective Townsend")
ax[0,1].set_title("Electron Townsend, Attachment and Effective Coefficient")
ax[0,1].set_xlabel(r"Efield [V/cm]")
ax[0,1].set_ylabel(r"$alpha$ $[cm_{-1}]$")

ax[1,0].scatter(ef,DL1,marker="x",color="blue")
ax[1,0].set_title("Electron Longitudinal Difussion")
ax[1,0].set_xlabel(r"Efield [V/cm]")
ax[1,0].set_ylabel(r"$D_L$ $[\sqrt{cm}]$")

ax[1,1].scatter(ef,DT1,marker="x",color="blue")
ax[1,1].set_title("Electron Transverse DIfussion")
ax[1,1].set_xlabel(r"Efield [V/cm]")
ax[1,1].set_ylabel(r"%$D_T$ $[\sqrt{cm}]$")

fig.tight_layout()
fig.savefig("38 -39 kV transport coefficients.png")


