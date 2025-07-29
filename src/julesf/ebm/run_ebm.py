# run_ebm.py

import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
    
from julesf.ebm.simulation import run_ebm
from julesf.ebm.parameters import sigma, epsilon_s

def day_fraction(t):
    return (t % 1.0)

# downwelling short and long wave
def Sdn(t, S_max=800.0, phi=0.0):
    θ = 2*np.pi*(day_fraction(t) - 0.25)
    val = S_max * np.sin(θ)
    return np.clip(val, 0, None)

def Ldn(t):
    Ta = Tair(t)
    θ = 2*np.pi*(day_fraction(t) - 0.25)
    ΔT = 20.0 if np.sin(θ) < 0 else 5.0
    T_sky = Ta - ΔT
    return epsilon_s * sigma * T_sky**4


# soil and air temperature

def Tsoil(t, Tsoil_mean=280.0, amp_day=2.0):
    θ = 2*np.pi*(day_fraction(t) - 0.25)
    return Tsoil_mean + amp_day * np.sin(θ)

def Tair(t, T_mean=283.0, amp_day=5.0):
    θ = 2*np.pi*(day_fraction(t) - 0.25)
    return T_mean + amp_day * np.sin(θ)

# specific humidity at fixed relative humidity

def Q1(t, RH=0.6):
    Ta = Tair(t)
    q_star = 0.01 * np.exp(0.07*(Ta - 273.15))
    return RH * q_star

# fractional cover

def nu(t, nu_mean=0.8):
    return nu_mean

drivers = {
    "Sdn":   Sdn,
    "Ldn":   Ldn,
    "Tair":  Tair,
    "Q1":    Q1,
    "nu":    nu,
    "Tsoil": Tsoil,
}

def main():

    t_span = (0, 10)
    T0 = 283.0

    t, Ts = run_ebm(t_span, T0, drivers, dt_out=0.05)

    plt.figure(figsize=(8,4))
    plt.plot(t, Ts, label="T_surf")
    plt.plot(t, [drivers["Tair"](tt) for tt in t], '--', label="T_air")
    plt.plot(t, [drivers["Tsoil"](tt) for tt in t], ':', label="T_soil")
    plt.xlabel("Time (days)")
    plt.ylabel("Temperature (K)")
    plt.legend()
    plt.title("Surface EBM (surface heat capacity = 2.0e3 J m⁻² K⁻¹)")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

