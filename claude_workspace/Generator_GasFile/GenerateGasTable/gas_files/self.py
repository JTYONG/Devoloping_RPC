import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gSystem.Load("libGarfield")
from ROOT import Garfield

gas = Garfield.MediumMagboltz()
gas.LoadGasFile("rpc_95.5_4.2_0.3_40-45_ncoll=20.gas")

import ctypes
ef, bf, ba = ROOT.std.vector('double')(), ROOT.std.vector('double')(), ROOT.std.vector('double')()
gas.GetFieldGrid(ef, bf, ba)

V = []
alpha = []
eta = []
eff_alpha = []

ib = ia = 0
for ie in range(ef.size()):
    print(f"Grid Point: {ie}")
    e = gas.UnScaleElectricField(ef[ie])
    v = ctypes.c_double(0.)
    v1 = ctypes.c_double(0.)
    v2 = ctypes.c_double(0.)

    if gas.GetElectronVelocityE(ie, ib, ia, v):
        print(f"field: {e}")
        print(f"v-value: {v.value}")
        print(f"velocity: {gas.ScaleVelocity(v.value)}")
        V.append(gas.ScaleVelocity(v.value))
    if gas.GetElectronTownsend(ie,ib,ia,v):
        print(f"v-value: {v.value}")
        print(f"townsend: {gas.ScaleTownsend(np.exp(v.value))}")
        alpha.append(np.exp(v.value))
    if gas.GetElectronAttachment(ie,ib,ia,v):
        print(f"v-value: {v.value}")
        print(f"attachment: {gas.ScaleAttachment(np.exp(v.value))}")
        eta.append(np.exp(v.value))
    if gas.GetElectronAttachment(ie,ib,ia,v1) and gas.GetElectronTownsend(ie,ib,ia,v2):
        print(f"v-value: {v2.value,v1.value}") 
        print(f"effective townsend : {gas.ScaleTownsend(np.exp(v2.value))-gas.ScaleAttachment(np.exp(v1.value))}")
        eff_alpha.append(gas.ScaleTownsend(np.exp(v2.value))-gas.ScaleAttachment(np.exp(v1.value)))
    print("\n")

V1 = np.array(V)
alpha1 = np.array(alpha)
eta1 = np.array(eta)
eff_alpha1 = np.array(eff_alpha)

print(f"V: {V1}")
print(f"alpha: {alpha1}")
print(f"eta: {eta1}")
print(f"eff_alpha: {eff_alpha1}")
print(len(ef))
print(len(V1))
print(len(alpha1))
print(len(eta1))
print(len(eff_alpha1))

#ea = gas.UnScaleElectricField(ef[ie])

fig,ax = plt.subplots(2,2,figsize=(18,18))

ax[0,0].plot(ef,V1,marker="x",color="blue")
ax[0,0].set_title("Electron Drift Velocity")
ax[0,0].set_xlabel(r"Efield [V/cm]")
ax[0,0].set_ylabel(r"$V_D$ [cm/ns]")


ax[0,1].plot(ef,alpha1,marker="x",color="blue")
ax[0,1].set_title("Electron Townsend Coefficient")
ax[0,1].set_xlabel(r"Efield [V/cm]")
ax[0,1].set_ylabel(r"$alpha$ $[cm_{-1}]$")


ax[1,0].plot(ef,eta1,marker="x",color="blue")
ax[1,0].set_title("Electron Attachment Coefficient")
ax[1,0].set_xlabel(r"Efield [V/cm]")
ax[1,0].set_ylabel(r"$eta$ $[cm_{-1}]$")

ax[1,1].plot(ef,eff_alpha1,marker="x",color="blue")
ax[1,1].set_title("Electroon Effective Townsend Coefficient")
ax[1,1].set_xlabel(r"Efield [V/cm]")
ax[1,1].set_ylabel(r"$alpha_{eff}$ $[cm^{-1}]$")

fig.tight_layout()
fig.savefig("transport coefficients.png")


