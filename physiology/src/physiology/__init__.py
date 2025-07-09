import numpy as np
import matplotlib.pyplot as plt

# Simulation settings
days = 10
dt = 1  # hourly timestep (hours)
t = np.arange(0, days * 24, dt)  # time in hours

# --- Forcings ---
# 1. Air temperature: sinusoid, mean 20°C, amplitude ±5°C
T_mean = 20.0  # °C
T_amp = 5.0    # °C
T = T_mean + T_amp * np.sin(2 * np.pi * t / 24)

# 2. PAR: zero at night, sinusoidal daytime peak ~1500 µmol m⁻² s⁻¹
I_max = 1500.0  # µmol m⁻² s⁻¹
I_par = np.maximum(0.0, I_max * np.sin(np.pi * t / 24))

# 3. Ambient CO₂ (constant 400 ppm → Pa)
P = 101325.0        # surface pressure, Pa
ca_ppm = 400.0      # ppm
ca = ca_ppm * 1e-6 * P  # convert to Pa

# 4. Leaf internal CO₂ (ci = χ · ca)
chi = 0.7
ci = chi * ca       # Pa

# 5. Ambient O₂
O2 = 0.21 * P       # 21% of total pressure, Pa

# --- PFT parameters for C₃ grass (from Table 2 + eqns) ---
alpha   = 0.12    # quantum efficiency (mol CO₂ mol⁻¹ photon)
omega   = 0.15    # leaf scattering coefficient
f_dr    = 0.015   # dark respiration coefficient
n0      = 0.073   # top-leaf N concentration (kg N per kg C)
n_e     = 0.0008  # conversion factor: mol CO₂ m⁻² s⁻¹ per kg C
Vcmax25 = n_e * n0  # base Vcmax at 25°C (mol CO₂ m⁻² s⁻¹)

# Temperature response Q10 values and RuBisCO specificity
Q10_Kc = 2.1
Q10_Ko = 1.2
Q10_rs = 0.57
T_low  = 0.0     # lower temp parameter (°C)
T_upp  = 36.0    # upper temp parameter (°C)

# --- Temperature-dependent parameter calculations ---
# Arrhenius/Q10 factor
f_T = Q10_Kc ** ((T - 25.0) / 10.0)

# Peaked temperature response for Vcmax
Vcmax = (
    Vcmax25 * f_T
    / ((1.0 + np.exp(0.3 * (T - T_upp)))
       * (1.0 + np.exp(0.3 * (T_low - T))))
)

# Michaelis-Menten constants at temperature T
Kc = 30.0 * Q10_Kc ** ((T - 25.0) / 10.0)
Ko = 3e4  * Q10_Ko ** ((T - 25.0) / 10.0)

# CO₂ compensation point Γ
tau   = 2600.0 * Q10_rs ** ((T - 25.0) / 10.0)  # specificity factor
Gamma = O2 / (2.0 * tau)

# --- Compute rate-limiters ---
Wc = Vcmax * (ci - Gamma) / (ci + Kc * (1.0 + O2 / Ko))        # Rubisco-limited
Wl = alpha * (1.0 - omega) * I_par * (ci - Gamma) / (ci + 2.0 * Gamma)  # Light-limited
We = 0.5 * Vcmax                                            # Export-limited (C₃)

# Simple co-limitation: gross assimilation is the min of the three
W_gross = np.minimum.reduce([Wc, Wl, We])

# --- Respiration and net photosynthesis ---
R_d = f_dr * Vcmax      # dark respiration
A_p = W_gross - R_d     # net potential photosynthesis

# --- Plotting ---
# 1) Rate limiters
plt.figure(figsize=(8, 4))
plt.plot(t, Wc, label='Wc (Rubisco)')
plt.plot(t, Wl, label='Wl (Light)')
plt.plot(t, We, label='We (Export)')
plt.xlabel('Time (h)')
plt.ylabel('Rate (mol CO₂ m⁻² s⁻¹)')
plt.title('Photosynthesis Rate Limiters')
plt.legend()
plt.grid(True)

# 2) Gross, respiration, and net potential
plt.figure(figsize=(8, 4))
plt.plot(t, W_gross, label='W (Gross)')
plt.plot(t, R_d, label='Rd (Dark resp.)')
plt.plot(t, A_p, label='Ap (Net potential)')
plt.xlabel('Time (h)')
plt.ylabel('Rate (mol CO₂ m⁻² s⁻¹)')
plt.title('Gross, Respiration, and Net Photosynthesis')
plt.legend()
plt.grid(True)

# 3) Forcing time series
plt.figure(figsize=(8, 4))
plt.plot(t, T, label='Air Temperature (°C)')
plt.plot(t, I_par, label='PAR (µmol m⁻² s⁻¹)')
plt.axhline(ca, linestyle='--', color='k', label='Ambient CO₂ (Pa)')
plt.xlabel('Time (h)')
plt.title('Forcing Time Series')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
