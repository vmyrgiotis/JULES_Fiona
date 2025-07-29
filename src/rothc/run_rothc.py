# run_rothc.py

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == "__main__":
    from simulation import run_soil_model
    from parameters import POOLS, C0_default
else:
    from .simulation import run_soil_model
    from .parameters import POOLS, C0_default

def seasonal_temp(t):   # in °C
    return 10.0 + 10.0*np.sin(2*np.pi*t/365)

def seasonal_moisture(t):
    return np.clip(0.5 + 0.3*np.sin(2*np.pi*t/365), 0, 1)

def seasonal_litter(t):
    return np.clip(0.002 + 0.001*np.sin(2*np.pi*t/365), 0, None)

drivers = {
    "Lambda_c": seasonal_litter,
    "T_soil"  : seasonal_temp,
    "s"       : seasonal_moisture,
    "nu"      : lambda t: 0.7,
}

C0 = [C0_default[p] for p in POOLS]
t_span = (0, 365)

t, C_q10   = run_soil_model(t_span, C0, drivers, temp_fun="Q10")
_, C_rothc = run_soil_model(t_span, C0, drivers, temp_fun="RothC")


plt.figure(figsize=(10,6))
styles = {'Q10':'-', 'RothC':'--'}
for i, pool in enumerate(POOLS):
    plt.plot(t, C_q10[i],   styles['Q10'],   label=f"{pool} (Q10)")
    plt.plot(t, C_rothc[i], styles['RothC'], label=f"{pool} (RothC)")
plt.xlabel("Time (days)")
plt.ylabel("Pool C (kg C m$^{-2}$)")
plt.title("Soil Carbon Pools")
plt.legend(ncol=2, fontsize='small')
plt.show()



pfts = ["tree", "shrub", "grass", "crop"]

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
axes = axes.flatten()

for ax, pft in zip(axes, pfts):
    t, C_all = run_soil_model(t_span, C0, drivers,
                              temp_fun="Q10",
                              pft=pft)
    for i, pool in enumerate(POOLS):
        ax.plot(t, C_all[i], label=pool)
    ax.set_title(f"PFT = {pft}")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("C (kg C m⁻²)")
    ax.legend(fontsize="small")

fig.suptitle("Soil‐C Pool Dynamics for Different PFTs (Q₁₀)", fontsize=16)
fig.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()