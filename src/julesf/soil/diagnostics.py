import numpy as np
import matplotlib.pyplot as plt

from julesf.soil.simulation import generate_soil_forcings, solve_soil_rhs
from julesf.soil.parameters import SOIL_LAYERS, SOIL_PROPERTIES, TIME_SETTINGS
from julesf.soil.equations_moisture import (
    vertical_water_flux,
    evapotranspiration_extraction,
    lateral_runoff
)
from julesf.soil.equations_thermal import (
    diffusive_heat_flux,
    advective_heat_flux
)

# Thermal properties needed by simple_thermal_conductivity()
THERMAL_PROPERTIES = {
    # soil thermal conductivity parameters (W m⁻¹ K⁻¹)
    'lambda_dry':   0.25,    # dry‐soil conductivity
    'lambda_sat':   2.5,     # saturated‐soil conductivity
    # volumetric heat capacities (J m⁻³ K⁻¹)
    'c_dry':        2.0e6,   # dry‐soil heat capacity
    'c_water':      4.18e6,  # liquid water heat capacity
    # porosity / sat. water content (m³ m⁻³)
    'theta_sat':    0.45,    # your soil’s porosity
}

def diagnose_soil_physics(days=None, dt_hours=None):
    # 1) generate forcings
    days     = days     or TIME_SETTINGS['days']
    dt_hours = dt_hours or TIME_SETTINGS['dt_hours']
    t_forcing, drivers = generate_soil_forcings(days=days, dt_hours=dt_hours)

    # 2) run solver
    t_solve, theta, T_soil = solve_soil_rhs(
        t_span=(0, days*24),
        theta_init=np.full(SOIL_LAYERS['n_layers'], TIME_SETTINGS.get('theta_init', 0.3)),
        T_init=np.full(SOIL_LAYERS['n_layers'], TIME_SETTINGS.get('T_init', 283.15)),
        drivers=drivers,
        method="RK45",
        dt_out=dt_hours
    )

    # 3) align time‐grids: solver often returns one extra point at t=end
    nt_f = len(t_forcing)
    nt_s = len(t_solve)
    nt   = min(nt_f, nt_s)

    t_forcing = t_forcing[:nt]
    t         = t_solve[:nt]
    theta     = theta[:, :nt]
    T_soil    = T_soil[:, :nt]

    # 4) allocate storage with aligned length
    n_layers = SOIL_LAYERS['n_layers']
    depths   = SOIL_LAYERS['depths']
    dz       = SOIL_LAYERS['layer_thickness']

    W = np.zeros((n_layers+1, nt))
    E = np.zeros((n_layers,   nt))
    R = np.zeros((n_layers,   nt))
    G = np.zeros((n_layers+1, nt))
    J = np.zeros((n_layers,   nt))

    # 5) loop over truncated time array
    for i in range(nt):
        th_i = theta[:, i]
        T_i  = T_soil[:, i]

        # moisture
        W[:, i] = vertical_water_flux(th_i, dz, SOIL_PROPERTIES)
        W[0, i] = drivers['precipitation'](t_forcing[i])
        E[:, i] = evapotranspiration_extraction(th_i, SOIL_PROPERTIES, drivers, t_forcing[i])
        R[:, i] = lateral_runoff(th_i, SOIL_PROPERTIES)

        # thermal
        G[:, i] = diffusive_heat_flux(T_i, dz, THERMAL_PROPERTIES, th_i)
        J[:, i] = advective_heat_flux(T_i, W[:, i], dz, THERMAL_PROPERTIES)

    return {
        't':        t,
        'depths':   depths,
        'theta':    theta,
        'W':        W,
        'E':        E,
        'R':        R,
        'T_soil':   T_soil,
        'G':        G,
        'J':        J
    }

def plot_soil_diagnostics(results):
    t      = results['t']
    depths = results['depths']
    theta  = results['theta']
    W      = results['W']
    E      = results['E']
    R      = results['R']
    Tsoil  = results['T_soil']
    G      = results['G']
    J      = results['J']

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle('Soil Physics Diagnostics', fontsize=14)

    # -- Top row: moisture & water
    im0 = axes[0,0].pcolormesh(t, depths, theta, shading='auto')
    fig.colorbar(im0, ax=axes[0,0], label='θ (m³/m³)')
    axes[0,0].set(title='Soil Moisture', xlabel='Time (h)', ylabel='Depth (m)')

    im1 = axes[0,1].pcolormesh(t, np.concatenate(([0], depths)), W, shading='auto')
    fig.colorbar(im1, ax=axes[0,1], label='W (kg/m²/s)')
    axes[0,1].set(title='Vertical Flux W', xlabel='Time (h)', ylabel='Interface depth')

    im2 = axes[0,2].pcolormesh(t, depths, E+R, shading='auto')
    fig.colorbar(im2, ax=axes[0,2], label='E+R (kg/m²/s)')
    axes[0,2].set(title='Evap (E) + Runoff (R)', xlabel='Time (h)')

    # -- Bottom row: temperature & heat flux
    im3 = axes[1,0].pcolormesh(t, depths, Tsoil-273.15, shading='auto')
    fig.colorbar(im3, ax=axes[1,0], label='T (°C)')
    axes[1,0].set(title='Soil Temperature', xlabel='Time (h)', ylabel='Depth (m)')

    im4 = axes[1,1].pcolormesh(t, np.concatenate(([0], depths)), G, shading='auto')
    fig.colorbar(im4, ax=axes[1,1], label='G (W/m²)')
    axes[1,1].set(title='Diffusive Flux G', xlabel='Time (h)')

    im5 = axes[1,2].pcolormesh(t, depths, J, shading='auto')
    fig.colorbar(im5, ax=axes[1,2], label='J (W/m²)')
    axes[1,2].set(title='Advective Flux J', xlabel='Time (h)')

    plt.tight_layout(rect=[0,0,1,0.95])
    plt.show()
    return fig

def plot_layer_fluxes(results):
    """
    For each soil layer, plot T_soil, the downward diffusive flux
    at its top interface G[k], and the advective flux J[k].
    """
    t      = results['t']
    depths = results['depths']
    Tsoil  = results['T_soil']
    G      = results['G']
    J      = results['J']
    nlay   = len(depths)

    fig, axes = plt.subplots(nlay, 1, figsize=(10, 2.5 * nlay), sharex=True)
    for k in range(nlay):
        ax = axes[k]
        ax.plot(t, Tsoil[k] - 273.15, 'k-', label='T (°C)')
        ax2 = ax.twinx()
        ax2.plot(t, G[k],   color='C1', ls='--', label='G (W/m²)')
        ax2.plot(t, J[k],   color='C2', ls=':',  label='J (W/m²)')

        ax.set_ylabel(f'Layer {k+1}\nT (°C)', color='k')
        ax2.set_ylabel('Flux (W/m²)', color='gray')
        ax.set_title(f'Layer {k+1} @ {depths[k]:.2f} m')
        ax.tick_params(axis='y',   labelcolor='k')
        ax2.tick_params(axis='y',  labelcolor='gray')

        # combined legend
        lines, labs = ax.get_legend_handles_labels()
        lines2, labs2 = ax2.get_legend_handles_labels()
        ax.legend(lines+lines2, labs+labs2, loc='upper right', fontsize=8)
        ax.grid(alpha=0.2)

    axes[-1].set_xlabel('Time (h)')
    plt.tight_layout(rect=[0,0,1,0.96])
    plt.show()
    return fig

if __name__ == "__main__":
    res = diagnose_soil_physics()
    plot_soil_diagnostics(res)
    plot_layer_fluxes(res)

    