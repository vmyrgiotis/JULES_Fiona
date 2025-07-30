# run_soil.py
import numpy as np
import matplotlib.pyplot as plt

from julesf.soil.simulation import generate_soil_forcings, solve_soil_rhs
from julesf.soil.parameters import SOIL_LAYERS

def run_soil_model():
    """
    Complete soil model workflow:
    1. Set up initial conditions and time settings
    2. Generate forcing data
    3. Run solver
    4. Visualize results
    """
    
    # Time settings
    days = 7
    dt_hours = 0.5
    t_span = (0, days * 24)  # hours
    
    # Initial conditions
    n_layers = SOIL_LAYERS['n_layers']
    theta_init = np.full(n_layers, 0.3)    # m³/m³, field capacity
    T_init = np.full(n_layers, 283.15)     # K, 10°C
    
    print("=== JULES Soil Model Simulation ===")
    print(f"Duration: {days} days ({t_span[1]} hours)")
    print(f"Timestep: {dt_hours} hours")
    print(f"Soil layers: {n_layers}")
    print(f"Layer depths: {SOIL_LAYERS['depths']} m")
    print(f"Initial moisture: {theta_init[0]:.2f} m³/m³")
    print(f"Initial temperature: {T_init[0]-273.15:.1f} °C")
    print()
    
    # Generate forcing data
    print("Generating soil forcings...")
    t_forcing, drivers = generate_soil_forcings(days=days, dt_hours=dt_hours)
    
    # Run solver
    t, theta, T_soil = solve_soil_rhs(
        t_span=t_span,
        theta_init=theta_init,
        T_init=T_init,
        drivers=drivers,
        method="RK45",
        dt_out=dt_hours
    )
    
    # Visualize results
    print("\nGenerating plots...")
    plot_soil_results(t, theta, T_soil, t_forcing, drivers)
    
    return t, theta, T_soil

def plot_soil_results(t, theta, T_soil, t_forcing, drivers):
    """
    Create comprehensive soil model visualization
    
    Parameters:
    - t: time array (hours)
    - theta: moisture content (n_layers × n_times)
    - T_soil: soil temperature (n_layers × n_times)
    - t_forcing: forcing time array
    - drivers: forcing functions
    """
    
    precip_series = [drivers['precipitation'](tt) for tt in t_forcing]
    ET_series = [drivers['evapotranspiration'](tt) for tt in t_forcing]
    T_air_series = [drivers['air_temperature'](tt) for tt in t_forcing]
    t_days = t / 24
    t_forcing_days = t_forcing / 24
    layer_names = [f"Layer {i+1}" for i in range(theta.shape[0])]
    depths = SOIL_LAYERS['depths']
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('JULES Soil Model Results', fontsize=16, y=0.98)  # Move title up slightly
    
    # Plot 1: Soil moisture profiles
    ax1 = axes[0, 0]
    for i in range(theta.shape[0]):
        ax1.plot(t_days, theta[i, :], label=f'{layer_names[i]} ({depths[i]:.1f}m)')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Soil Moisture (m³/m³)')
    ax1.set_title('Soil Moisture Content')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Soil temperature profiles  
    ax2 = axes[0, 1]
    for i in range(T_soil.shape[0]):
        ax2.plot(t_days, T_soil[i, :]-273.15, label=f'{layer_names[i]} ({depths[i]:.1f}m)')
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Soil Temperature (°C)')
    ax2.set_title('Soil Temperature')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Atmospheric forcings
    ax3 = axes[1, 0]
    ax3_twin = ax3.twinx()
    
    # Precipitation (scale up for visibility)
    p1 = ax3.bar(t_forcing_days, np.array(precip_series)*3600*1000, 
                 width=0.02, alpha=0.6, color='blue', label='Precipitation')
    ax3.set_ylabel('Precipitation (mm/hr)', color='blue')
    ax3.tick_params(axis='y', labelcolor='blue')
    
    # Evapotranspiration  
    p2 = ax3_twin.plot(t_forcing_days, np.array(ET_series)*3600*1000, 
                       'r-', label='Evapotranspiration')
    ax3_twin.set_ylabel('ET (mm/hr)', color='red')
    ax3_twin.tick_params(axis='y', labelcolor='red')
    
    ax3.set_xlabel('Time (days)')
    ax3.set_title('Water Fluxes')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Air temperature
    ax4 = axes[1, 1]
    ax4.plot(t_forcing_days, np.array(T_air_series)-273.15, 'k-', linewidth=2)
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Air Temperature (°C)')
    ax4.set_title('Air Temperature')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for suptitle
    plt.subplots_adjust(hspace=0.3, wspace=0.3)  # Add horizontal and vertical spacing
    plt.show()
    
    # Print summary statistics
    print("\n=== Simulation Summary ===")
    print(f"Moisture range: {theta.min():.3f} - {theta.max():.3f} m³/m³")
    print(f"Temperature range: {T_soil.min()-273.15:.1f} - {T_soil.max()-273.15:.1f} °C")
    print(f"Total precipitation: {np.sum(precip_series)*3600*len(t_forcing)/len(t_forcing)*24:.1f} mm")
    print(f"Total ET: {np.sum(ET_series)*3600*len(t_forcing)/len(t_forcing)*24:.1f} mm")

if __name__ == "__main__":
    run_soil_model()
