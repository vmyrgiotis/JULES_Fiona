import numpy as np
import matplotlib.pyplot as plt

# 1.1 Simulation settings
days = 10
dt = 1  # hourly timestep (hours)
t = np.arange(0, days * 24, dt)  # time in hours

# 1.2 Forcings

T_mean = 20.0           # °C - air mean T 20°C
T_amp = 5.0             # °C - T change amplitude ±5°C
T = T_mean + T_amp * np.sin(2 * np.pi * t / 24)
I_max = 1500.0          # PAR peak (during the day) µmol m⁻² s⁻¹
I_par = np.maximum(0.0, I_max * np.sin(np.pi * t / 24))
P = 101325.0            # surface pressure, Pa
ca_ppm = 400.0          # ambient CO2 in ppm
ca = ca_ppm * 1e-6 * P  # ambient CO2 in Pa
chi = 0.7               # ratio of leaf internal CO2 to ambient CO2
ci = chi * ca           # leaf internal CO2 in Pa
O2 = 0.21 * P           # Ambient O2, 21% of total pressure, Pa

# 1.3 Parameters for C₃ grass (from Table 2) 
alpha   = 0.12          # quantum efficiency (mol CO₂ mol⁻¹ photon)
omega   = 0.15          # leaf scattering coefficient
f_dr    = 0.015         # dark respiration coefficient
n0      = 0.073         # top-leaf N concentration (kg N per kg C)
n_e     = 0.0008        # conversion factor: mol CO₂ m⁻² s⁻¹ per kg C
Vcmax25 = n_e * n0      # base Vcmax at 25°C (mol CO₂ m⁻² s⁻¹)

# 1.4 Constants and default parameter values 
Q10_Kc = 2.1
Q10_Ko = 1.2
Q10_rs = 0.57    # RuBisCO specificity
T_low  = 0.0     # lower temp parameter (°C)
T_upp  = 36.0    # upper temp parameter (°C)



# 2.Equations

# 2.1 compute temperature-dependent variables 
# Eq.5, Arrhenius/Q10 factor
f_T = Q10_Kc ** ((T - 25.0) / 10.0)

# Eq.6, Vcmax
Vcmax = (
    Vcmax25 * f_T
    / ((1.0 + np.exp(0.3 * (T - T_upp)))
       * (1.0 + np.exp(0.3 * (T_low - T))))
)

# Eq.7-8 CO₂ compensation point Γ
tau   = 2600.0 * Q10_rs ** ((T - 25.0) / 10.0)  
Gamma = O2 / (2.0 * tau)

# Eq.9-10, Michaelis-Menten constants
Kc = 30.0 * Q10_Kc ** ((T - 25.0) / 10.0)
Ko = 3e4  * Q10_Ko ** ((T - 25.0) / 10.0)

# 2.2 compute limit rates, Eq.1-3 (note that in the future equations for C4 plants are also needed)
Wc = Vcmax * (ci - Gamma) / (ci + Kc * (1.0 + O2 / Ko))                 # Rubisco-limited
Wl = alpha * (1.0 - omega) * I_par * (ci - Gamma) / (ci + 2.0 * Gamma)  # Light-limited
We = 0.5 * Vcmax                                                        # Export-limited 

# 2.3 Gross, respiration and net potential photosynthesis, Eq.11-14
W_gross = np.minimum.reduce([Wc, Wl, We])
R_d = f_dr * Vcmax      # dark respiration
A_p = W_gross - R_d     # net potential photosynthesis





# 3. visualization

# 3.1 Rate limiters
plt.figure(figsize=(8, 4))
plt.plot(t, Wc, label='Wc (Rubisco)')
plt.plot(t, Wl, label='Wl (Light)')
plt.plot(t, We, label='We (Export)')
plt.xlabel('Time (h)')
plt.ylabel('Rate (mol CO₂ m⁻² s⁻¹)')
plt.title('Photosynthesis Rate Limiters')
plt.legend()
plt.grid(True)

# 3.2 Gross, respiration, and net potential
plt.figure(figsize=(8, 4))
plt.plot(t, W_gross, label='W (Gross)')
plt.plot(t, R_d, label='Rd (Dark resp.)')
plt.plot(t, A_p, label='Ap (Net potential)')
plt.xlabel('Time (h)')
plt.ylabel('Rate (mol CO₂ m⁻² s⁻¹)')
plt.title('Gross, Respiration, and Net Photosynthesis')
plt.legend()
plt.grid(True)

# 3.3 Forcing time series
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
