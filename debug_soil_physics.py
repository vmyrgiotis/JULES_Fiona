"""
Debug script to diagnose soil physics issues using ERA5 forcing
Based on src/julesf/soil/diagnostics.py approach
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.append('src')

from julesf.data.era5_data_reader import load_era5_for_jules
from scripts.era5_simulation_run import create_era5_config
from julesf.soil.simulation import solve_soil_rhs
from julesf.soil.parameters import SOIL_LAYERS, SOIL_PROPERTIES, VAN_GENUCHTEN, THERMAL_PROPERTIES
from julesf.soil.equations_moisture import (
    vertical_water_flux,
    evapotranspiration_extraction,
    lateral_runoff,
    infiltration_rate
)
from julesf.soil.equations_thermal import (
    diffusive_heat_flux,
    advective_heat_flux,
    thermal_rhs
)

def debug_soil_physics_era5(days=2):
    """
    Debug soil physics using actual ERA5 forcing to find what's wrong
    """
    print("=== DEBUGGING SOIL PHYSICS WITH ERA5 ===")
    print("="*60)
    
    # 1) Load ERA5 data and create forcing
    print("Loading ERA5 data...")
    data_file = "data/era5_jules_met_ins_2002-2016.csv"
    metadata_file = "data/era5_jules_met_ins_metadata.csv"
    
    era5_data = load_era5_for_jules(data_file, metadata_file, timestep_hours=0.5)
    config = create_era5_config(era5_data['drivers'], simulation_weeks=1)
    
    # Extract soil drivers
    era5_forcing = config['era5_forcing']
    soil_drivers = {
        'precipitation': era5_forcing['precip'],
        'evapotranspiration': era5_forcing['evapotranspiration'],
        'air_temperature': era5_forcing['Tair'],
        'surface_heat_flux': lambda t: 0.0  # Simplified for now
    }
    
    print(f"✓ ERA5 forcing loaded for {days} days")
    
    # 2) Test forcing functions
    print("\nTesting forcing functions:")
    test_times = [0, 6, 12, 18, 24, 30, 36, 42, 48]
    for t in test_times:
        precip = soil_drivers['precipitation'](t) * 3600 * 1000  # mm/hr
        et = soil_drivers['evapotranspiration'](t) * 3600 * 1000  # mm/hr
        temp = soil_drivers['air_temperature'](t) - 273.15  # °C
        print(f"  t={t:2d}h: P={precip:5.2f} mm/hr, ET={et:5.2f} mm/hr, T={temp:5.1f}°C")
    
    # 3) Set up soil simulation
    dt_hours = 0.5
    t_span = (0, days * 24)
    
    # Initial conditions
    n_layers = SOIL_LAYERS['n_layers']
    theta_init = np.full(n_layers, 0.3)    # Field capacity
    T_init = np.full(n_layers, 283.15)     # 10°C
    
    print(f"\nSoil setup:")
    print(f"  Layers: {n_layers}")
    print(f"  Depths: {SOIL_LAYERS['depths']} m")
    print(f"  Layer thickness: {SOIL_LAYERS['layer_thickness']} m")
    print(f"  Initial θ: {theta_init}")
    print(f"  Initial T: {T_init - 273.15} °C")
    
    # 4) **FORCE BDF METHOD** for solving
    print(f"\nRunning soil solver with BDF method...")
    print(f"  Time span: {t_span} hours")
    print(f"  Timestep: {dt_hours} hours")
    
    t_solve, theta, T_soil = solve_soil_rhs(
        t_span=t_span,
        theta_init=theta_init,
        T_init=T_init,
        drivers=soil_drivers,
        method="BDF",  # ✅ FORCE BDF METHOD
        dt_out=dt_hours
    )
    
    print(f"✓ Solver completed: {len(t_solve)} time points")
    print(f"  θ range: {theta.min():.3f} - {theta.max():.3f} m³/m³")
    print(f"  T range: {T_soil.min():.1f} - {T_soil.max():.1f} K")
    
    # 5) Calculate intermediate fluxes for each time step
    print("\nCalculating intermediate fluxes...")
    
    # Combine parameters
    params = {**SOIL_LAYERS, **SOIL_PROPERTIES, **VAN_GENUCHTEN, **THERMAL_PROPERTIES}
    
    nt = len(t_solve)
    depths = SOIL_LAYERS['depths']
    dz = SOIL_LAYERS['layer_thickness']
    
    # Storage arrays
    W = np.zeros((n_layers+1, nt))  # Vertical water flux
    E = np.zeros((n_layers, nt))    # Evapotranspiration extraction  
    R = np.zeros((n_layers, nt))    # Lateral runoff
    G = np.zeros((n_layers+1, nt))  # Diffusive heat flux
    J = np.zeros((n_layers, nt))    # Advective heat flux
    
    # Calculate fluxes at each time step
    for i in range(nt):
        t_i = t_solve[i]
        theta_i = theta[:, i]
        T_i = T_soil[:, i]
        
        # Water fluxes
        W[:, i] = vertical_water_flux(theta_i, dz, params)
        W[0, i] = infiltration_rate(params, soil_drivers, t_i)  # Surface infiltration
        E[:, i] = evapotranspiration_extraction(theta_i, params, soil_drivers, t_i)
        R[:, i] = lateral_runoff(theta_i, params)
        
        # Heat fluxes  
        G[:, i] = diffusive_heat_flux(T_i, dz, params, theta_i)
        J[:, i] = advective_heat_flux(T_i, W[:, i], dz, params)
    
    # 6) Create comprehensive diagnostic plots
    results = {
        't': t_solve,
        'depths': depths,
        'theta': theta,
        'T_soil': T_soil,
        'W': W,
        'E': E, 
        'R': R,
        'G': G,
        'J': J,
        'soil_drivers': soil_drivers
    }
    
    plot_soil_physics_diagnostics(results)
    plot_layer_by_layer_analysis(results)
    
    return results

def plot_soil_physics_diagnostics(results):
    """
    Comprehensive soil physics diagnostic plots - COMPACT VERSION
    """
    t = results['t']
    depths = results['depths']
    theta = results['theta']
    T_soil = results['T_soil']
    W = results['W']
    E = results['E']
    R = results['R']
    G = results['G']
    J = results['J']
    
    print("\nCreating soil physics diagnostic plots...")
    
    # SMALLER, MORE COMPACT FIGURE
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))  # Reduced from (18, 15)
    fig.suptitle('Soil Physics Diagnostics with ERA5 Forcing', fontsize=14)
    
    # === ROW 1: MOISTURE DYNAMICS ===
    # Soil moisture content
    im0 = axes[0,0].pcolormesh(t/24, depths, theta, shading='auto', cmap='Blues')
    cbar0 = fig.colorbar(im0, ax=axes[0,0], shrink=0.8)
    cbar0.set_label('θ (m³/m³)', fontsize=9)
    cbar0.ax.tick_params(labelsize=8)
    axes[0,0].set_title('Soil Moisture Content', fontsize=10)
    axes[0,0].set_xlabel('Time (days)', fontsize=9)
    axes[0,0].set_ylabel('Depth (m)', fontsize=9)
    axes[0,0].tick_params(labelsize=8)
    axes[0,0].invert_yaxis()
    
    # Vertical water flux
    im1 = axes[0,1].pcolormesh(t/24, np.concatenate(([0], depths)), W, shading='auto', cmap='RdBu_r')
    cbar1 = fig.colorbar(im1, ax=axes[0,1], shrink=0.8)
    cbar1.set_label('W (kg/m²/s)', fontsize=9)
    cbar1.ax.tick_params(labelsize=8)
    axes[0,1].set_title('Vertical Water Flux', fontsize=10)
    axes[0,1].set_xlabel('Time (days)', fontsize=9)
    axes[0,1].set_ylabel('Interface depth (m)', fontsize=9)
    axes[0,1].tick_params(labelsize=8)
    axes[0,1].invert_yaxis()
    
    # Evapotranspiration extraction
    im2 = axes[0,2].pcolormesh(t/24, depths, E, shading='auto', cmap='Reds')
    cbar2 = fig.colorbar(im2, ax=axes[0,2], shrink=0.8)
    cbar2.set_label('E (kg/m²/s)', fontsize=9)
    cbar2.ax.tick_params(labelsize=8)
    axes[0,2].set_title('ET Extraction by Layer', fontsize=10)
    axes[0,2].set_xlabel('Time (days)', fontsize=9)
    axes[0,2].set_ylabel('Depth (m)', fontsize=9)
    axes[0,2].tick_params(labelsize=8)
    axes[0,2].invert_yaxis()
    
    # === ROW 2: THERMAL DYNAMICS ===
    # Soil temperature
    im3 = axes[1,0].pcolormesh(t/24, depths, T_soil-273.15, shading='auto', cmap='coolwarm')
    cbar3 = fig.colorbar(im3, ax=axes[1,0], shrink=0.8)
    cbar3.set_label('T (°C)', fontsize=9)
    cbar3.ax.tick_params(labelsize=8)
    axes[1,0].set_title('Soil Temperature', fontsize=10)
    axes[1,0].set_xlabel('Time (days)', fontsize=9)
    axes[1,0].set_ylabel('Depth (m)', fontsize=9)
    axes[1,0].tick_params(labelsize=8)
    axes[1,0].invert_yaxis()
    
    # Diffusive heat flux
    im4 = axes[1,1].pcolormesh(t/24, np.concatenate(([0], depths)), G, shading='auto', cmap='RdBu_r')
    cbar4 = fig.colorbar(im4, ax=axes[1,1], shrink=0.8)
    cbar4.set_label('G (W/m²)', fontsize=9)
    cbar4.ax.tick_params(labelsize=8)
    axes[1,1].set_title('Diffusive Heat Flux', fontsize=10)
    axes[1,1].set_xlabel('Time (days)', fontsize=9)
    axes[1,1].set_ylabel('Interface depth (m)', fontsize=9)
    axes[1,1].tick_params(labelsize=8)
    axes[1,1].invert_yaxis()
    
    # Advective heat flux
    im5 = axes[1,2].pcolormesh(t/24, depths, J, shading='auto', cmap='RdBu_r')
    cbar5 = fig.colorbar(im5, ax=axes[1,2], shrink=0.8)
    cbar5.set_label('J (W/m²)', fontsize=9)
    cbar5.ax.tick_params(labelsize=8)
    axes[1,2].set_title('Advective Heat Flux', fontsize=10)
    axes[1,2].set_xlabel('Time (days)', fontsize=9)
    axes[1,2].set_ylabel('Depth (m)', fontsize=9)
    axes[1,2].tick_params(labelsize=8)
    axes[1,2].invert_yaxis()
    
    # === ROW 3: FORCING AND BALANCE ===
    # Surface forcing
    drivers = results['soil_drivers']
    t_hours = t
    precip_series = [drivers['precipitation'](tt) * 3600 * 1000 for tt in t_hours]  # mm/hr
    et_total_series = [drivers['evapotranspiration'](tt) * 3600 * 1000 for tt in t_hours]  # mm/hr
    
    axes[2,0].bar(t/24, precip_series, width=0.01, alpha=0.6, color='blue', label='Precipitation')
    ax_twin = axes[2,0].twinx()
    ax_twin.plot(t/24, et_total_series, 'r-', label='Total ET', linewidth=1)
    axes[2,0].set_ylabel('Precipitation (mm/hr)', color='blue', fontsize=9)
    ax_twin.set_ylabel('ET (mm/hr)', color='red', fontsize=9)
    axes[2,0].set_title('Surface Water Forcing', fontsize=10)
    axes[2,0].set_xlabel('Time (days)', fontsize=9)
    axes[2,0].tick_params(labelsize=8)
    ax_twin.tick_params(labelsize=8)
    
    # Water balance per layer
    for layer in range(len(depths)):
        net_water = W[layer, :] - W[layer+1, :] - E[layer, :] - R[layer, :]
        axes[2,1].plot(t/24, net_water * 3600 * 1000, label=f'L{layer+1}', linewidth=1)
    axes[2,1].set_title('Net Water Balance by Layer', fontsize=10)
    axes[2,1].set_ylabel('Net flux (mm/hr)', fontsize=9)
    axes[2,1].set_xlabel('Time (days)', fontsize=9)
    axes[2,1].legend(fontsize=7, ncol=2)
    axes[2,1].grid(True, alpha=0.3)
    axes[2,1].tick_params(labelsize=8)
    
    # Temperature forcing and gradients
    temp_series = [drivers['air_temperature'](tt) - 273.15 for tt in t_hours]
    axes[2,2].plot(t/24, temp_series, 'k-', linewidth=2, label='Air temp')
    for layer in range(len(depths)):
        axes[2,2].plot(t/24, T_soil[layer, :] - 273.15, label=f'Soil L{layer+1}', linewidth=1)
    axes[2,2].set_title('Temperature Profile', fontsize=10)
    axes[2,2].set_ylabel('Temperature (°C)', fontsize=9)
    axes[2,2].set_xlabel('Time (days)', fontsize=9)
    axes[2,2].legend(fontsize=7, ncol=2)
    axes[2,2].grid(True, alpha=0.3)
    axes[2,2].tick_params(labelsize=8)
    
    plt.tight_layout()
    plt.savefig('debug_soil_physics_diagnostics.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✓ Diagnostic plots saved: debug_soil_physics_diagnostics.png")

def plot_layer_by_layer_analysis(results):
    """
    Individual layer analysis plots - COMPACT VERSION
    """
    t = results['t']
    depths = results['depths']
    theta = results['theta']
    T_soil = results['T_soil']
    W = results['W']
    E = results['E']
    G = results['G']
    J = results['J']
    n_layers = len(depths)
    
    print("Creating layer-by-layer analysis plots...")
    
    # SMALLER FIGURE WITH TIGHTER LAYOUT
    fig, axes = plt.subplots(n_layers, 2, figsize=(12, 2.5*n_layers))  # Reduced height per layer
    fig.suptitle('Layer-by-Layer Soil Physics Analysis', fontsize=14)
    
    for k in range(n_layers):
        # Left plot: Moisture and water fluxes
        ax1 = axes[k, 0]
        ax1_twin = ax1.twinx()
        
        # Moisture content (primary axis)
        ax1.plot(t/24, theta[k, :], 'b-', linewidth=1.5, label='θ (m³/m³)')
        ax1.set_ylabel('θ (m³/m³)', color='blue', fontsize=9)
        ax1.tick_params(axis='y', labelcolor='blue', labelsize=8)
        
        # Water fluxes (secondary axis)
        ax1_twin.plot(t/24, W[k, :] * 3600 * 1000, 'g--', label='W_in', linewidth=1)
        ax1_twin.plot(t/24, W[k+1, :] * 3600 * 1000, 'r--', label='W_out', linewidth=1)
        ax1_twin.plot(t/24, E[k, :] * 3600 * 1000, 'm:', label='ET', linewidth=1)
        ax1_twin.set_ylabel('Water flux (mm/hr)', color='gray', fontsize=9)
        ax1_twin.tick_params(axis='y', labelcolor='gray', labelsize=8)
        
        ax1.set_title(f'Layer {k+1} @ {depths[k]:.2f}m - Water Balance', fontsize=10)
        ax1.set_xlabel('Time (days)', fontsize=9)
        ax1.grid(True, alpha=0.3)
        ax1.tick_params(labelsize=8)
        
        # Combined legend - SMALLER
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax1_twin.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=7, ncol=2)
        
        # Right plot: Temperature and heat fluxes
        ax2 = axes[k, 1]
        ax2_twin = ax2.twinx()
        
        # Temperature (primary axis)
        ax2.plot(t/24, T_soil[k, :] - 273.15, 'r-', linewidth=1.5, label='T (°C)')
        ax2.set_ylabel('T (°C)', color='red', fontsize=9)
        ax2.tick_params(axis='y', labelcolor='red', labelsize=8)
        
        # Heat fluxes (secondary axis)
        ax2_twin.plot(t/24, G[k, :], 'orange', linestyle='--', label='G_in', linewidth=1)
        ax2_twin.plot(t/24, G[k+1, :], 'purple', linestyle='--', label='G_out', linewidth=1)
        ax2_twin.plot(t/24, J[k, :], 'brown', linestyle=':', label='J_adv', linewidth=1)
        ax2_twin.set_ylabel('Heat flux (W/m²)', color='gray', fontsize=9)
        ax2_twin.tick_params(axis='y', labelcolor='gray', labelsize=8)
        
        ax2.set_title(f'Layer {k+1} @ {depths[k]:.2f}m - Heat Balance', fontsize=10)
        ax2.set_xlabel('Time (days)', fontsize=9)
        ax2.grid(True, alpha=0.3)
        ax2.tick_params(labelsize=8)
        
        # Combined legend - SMALLER
        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2_twin.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=7, ncol=2)
    
    plt.tight_layout()
    plt.savefig('debug_soil_layer_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✓ Layer analysis plots saved: debug_soil_layer_analysis.png")

def check_solver_method():
    """
    Verify that BDF method is actually being used
    """
    print("\n=== CHECKING SOLVER METHOD ===")
    
    # Look at the solve_soil_rhs function to confirm BDF usage
    import inspect
    from julesf.soil.simulation import solve_soil_rhs
    
    source = inspect.getsource(solve_soil_rhs)
    if 'method=method' in source and 'BDF' in source:
        print("✓ solve_soil_rhs accepts method parameter")
    else:
        print("❌ solve_soil_rhs may not be using the method parameter correctly!")
    
    # Check solve_ivp call
    if 'solve_ivp' in source:
        print("✓ Using scipy.integrate.solve_ivp")
    else:
        print("❌ Not using solve_ivp!")
    
    print(f"Method parameter will be: 'BDF'")

if __name__ == "__main__":
    # Check that BDF method will be used
    check_solver_method()
    
    # Run comprehensive soil physics debugging
    print("\nStarting comprehensive soil physics debug...")
    results = debug_soil_physics_era5(days=2)
    
    print(f"\n=== SUMMARY ===")
    print(f"θ variation: {results['theta'].max() - results['theta'].min():.6f} m³/m³")
    print(f"T variation: {results['T_soil'].max() - results['T_soil'].min():.3f} K")
    print(f"Max W flux: {results['W'].max() * 3600 * 1000:.3f} mm/hr")
    print(f"Max ET flux: {results['E'].max() * 3600 * 1000:.3f} mm/hr")
    print(f"Max heat flux G: {results['G'].max():.1f} W/m²")
    
    if results['theta'].max() - results['theta'].min() < 1e-6:
        print("❌ MOISTURE IS CONSTANT - Problem in moisture equation!")
    
    if results['T_soil'].max() - results['T_soil'].min() < 0.1:
        print("❌ TEMPERATURE IS NEARLY CONSTANT - Problem in thermal equation!")
    
    print("\nCheck the diagnostic plots to identify the specific issues!")