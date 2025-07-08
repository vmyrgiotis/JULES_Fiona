import numpy as np
import matplotlib.pyplot as plt

nx, ny = 5, 5                  # grid dimensions
days = 365                     # simulation period (days)
dt = 24 * 3600                 # daily timestep in seconds
t = np.arange(days)            # time array [0, 1, …, 364]

T0 = 1e-5                      # maximum transmissivity (m²/s)
f = 0.2                        # transmissivity decay coeff. (m⁻¹)
a_s, c_s = 0.3, 1.0            # parameters for Eq 59, f_sat = a_s * exp(-c_s * lambda_c)
theta_eff = 0.2                # effective porosity (m³/m³)

# Generate grid of mean topographic index λ̄ in the range [3, 6]
np.random.seed(0)
lambda_bar = np.random.uniform(3.0, 6.0, size=(nx, ny))

# Forcing W0(t):
# seasonal sinusoid around 1000 mm/yr = 1 m / (365*86400 s)
# a sine wave with peak runoff during mid-year
W0_mean = 1.0 / (365 * 86400)
W0 = W0_mean * (1 + 0.5 * np.sin(2 * np.pi * t / days))

# Initialize
z_bar    = np.full((days, nx, ny), 2.0)   # initial water-table depth = 2 m
Rb       = np.zeros((days, nx, ny))       # subsurface baseflow
lambda_c = np.zeros((days, nx, ny))       # critical topo index
f_sat    = np.zeros((days, nx, ny))       # saturated fraction
Rse      = np.zeros((days, nx, ny))       # excess runoff

# Precompute Rb_max
Rb_max = T0 * np.exp(-lambda_bar)

for i in range(days - 1):
    # Transmissivity and baseflow, Eq. 68
    Tz    = T0 * np.exp(-f * z_bar[i])
    Rb[i] = Tz * np.exp(-lambda_bar)

    # Critical topo-index and saturated fraction, Eq. 69
    lambda_c[i] = np.log(Rb_max / Rb[i])
    f_sat[i]    = a_s * np.exp(-c_s * f * lambda_c[i])

    # Saturation-excess runoff, Eq. 70
    Rse[i] = f_sat[i] * W0[i]

    # Update water-table depth via simple mass balance
    dz = (W0[i] - Rb[i]) * dt / theta_eff
    z_bar[i + 1] = np.maximum(z_bar[i] + dz, 0.0)  # z ≥ 0

# Diagnostics for final day
i = days - 1
Tz    = T0 * np.exp(-f * z_bar[i])
Rb[i] = Tz * np.exp(-lambda_bar)
lambda_c[i] = np.log(Rb_max / Rb[i])
f_sat[i]    = a_s * np.exp(-c_s * lambda_c[i])
Rse[i]      = f_sat[i] * W0[i]

# Visualization

# Rb and Rse
cells = [(0,0), (0,4), (4,0), (4,4)]
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

for (ix, iy) in cells:
    axes[0].plot(t, Rb[:, ix, iy], label=f'Rb cell {ix},{iy}')
    axes[1].plot(t, Rse[:, ix, iy], label=f'Rse cell {ix},{iy}')

axes[0].set_ylabel('Baseflow $R_b$ (kg m$^{-2}$ s$^{-1}$)')
axes[1].set_ylabel('Excess runoff $R_{se}$ (kg m$^{-2}$ s$^{-1}$)')
axes[1].set_xlabel('Time (days)')
axes[0].legend(loc='upper right')
axes[1].legend(loc='upper right')
axes[0].set_title('Subsurface Baseflow and Saturation-Excess Runoff')
plt.tight_layout()
plt.show()

# Catchment-mean intermediate variables
Rb_mean       = Rb.mean(axis=(1,2))
lambda_c_mean = lambda_c.mean(axis=(1,2))
f_sat_mean    = f_sat.mean(axis=(1,2))

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(t, lambda_c_mean, label='Mean $\\lambda_c$')
ax.plot(t, f_sat_mean,    label='Mean $f_{sat}$')
ax.set_ylabel('Value')
ax.set_xlabel('Time (days)')
ax.set_title('Intermediate Variables: Critical Index & Saturated Fraction')
ax.legend()
plt.tight_layout()
plt.show()
