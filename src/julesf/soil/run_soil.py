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
    Create comprehensive soil model visualization using a professional Morandi/sci-style palette.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from julesf.soil.parameters import SOIL_LAYERS

    # Morandi/sci-style palette for soil layers (deep muted to shallow bright)
    layer_colors = [
        "#D8A0C7", 
        "#5FA434",  
        "#77C2F3",  
        "#2B307A",  
    ]

    # Professional colors for atmospheric forcings
    precip_color = "#164ED1"   # Blue for precipitation
    et_color     = "#CD1F20"   # Red-orange for ET
    air_temp_color = "#000000" # Green for air temperature

    # Evaluate drivers for forcing time series
    # now drivers['precipitation'] and drivers['evapotranspiration']
    # both return kg/m²/s, so ×3600 → kg/m²/hr → mm/hr
    precip_series = [drivers['precipitation'](tt) * 3600.0 for tt in t_forcing]
    ET_series     = [drivers['evapotranspiration'](tt) * 3600.0 for tt in t_forcing]
    T_air_series  = [drivers['air_temperature'](tt) for tt in t_forcing]

    t_days = t / 24
    t_forcing_days = t_forcing / 24
    layer_names = [f"Layer {i+1}" for i in range(theta.shape[0])]
    depths = SOIL_LAYERS['depths']

    # Set up the figure with a professional layout
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    # fig.suptitle('JULES Soil Model Results', fontsize=16, y=0.98)

    # Plot 1: Soil moisture profiles
    ax1 = axes[0, 0]
    for i in range(theta.shape[0]):
        ax1.plot(
            t_days, theta[i, :], color=layer_colors[i],
            label=f'{depths[i]:.1f} m', linewidth=2.2
        )
    ax1.set_xlabel('Time (days)', fontsize=11)
    ax1.set_ylabel('Soil Moisture (m³/m³)', fontsize=11)
    ax1.set_title('Soil Moisture Content', fontsize=13)
    leg1 = ax1.legend(fontsize=10, title='Soil layer')
    leg1.get_title().set_fontsize(11)
    ax1.grid(False)

    # Plot 2: Soil temperature profiles
    ax2 = axes[0, 1]
    for i in range(T_soil.shape[0]):
        ax2.plot(
            t_days, T_soil[i, :] - 273.15, color=layer_colors[i],
            label=f'{depths[i]:.1f} m', linewidth=2
        )
    ax2.set_xlabel('Time (days)', fontsize=11)
    ax2.set_ylabel('Soil Temperature (°C)', fontsize=11)
    ax2.set_title('Soil Temperature', fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(False)

    # Plot 3: Atmospheric forcings
    ax3 = axes[1, 0]
    ax3_twin = ax3.twinx()
    # Precipitation as bar plot
    p1 = ax3.bar(
        t_forcing_days, np.array(precip_series), width=0.02, alpha=0.7,
        color=precip_color, label='Precipitation'
    )
    ax3.set_ylabel('Precipitation (mm/hr)', color=precip_color, fontsize=11)
    ax3.tick_params(axis='y', labelcolor=precip_color)
    # ET as line plot on twin axis
    p2 = ax3_twin.plot(
        t_forcing_days, np.array(ET_series), color=et_color,
        label='Evapotranspiration', linewidth=2
    )
    ax3_twin.set_ylabel('ET (mm/hr)', color=et_color, fontsize=11)
    ax3_twin.tick_params(axis='y', labelcolor=et_color)
    ax3.set_xlabel('Time (days)', fontsize=11)
    ax3.set_title('Water Fluxes', fontsize=13)
    ax3.grid(False)
    # Legends
    ax3.legend([p1], ['Precipitation'], fontsize=10, loc='upper left')
    ax3_twin.legend([p2[0]], ['Evapotranspiration'], fontsize=10, loc='upper right')

    # Plot 4: Air temperature
    ax4 = axes[1, 1]
    ax4.plot(
        t_forcing_days, np.array(T_air_series) - 273.15, color=air_temp_color,
        linewidth=2, label='Air Temperature'
    )
    ax4.set_xlabel('Time (days)', fontsize=11)
    ax4.set_ylabel('Air Temperature (°C)', fontsize=11)
    ax4.set_title('Air Temperature', fontsize=13)
    ax4.grid(False)
    ax4.legend(fontsize=10)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.subplots_adjust(hspace=0.32, wspace=0.28)
    plt.show()

    # Print summary statistics
    print("\n=== Simulation Summary ===")
    print(f"Moisture range: {theta.min():.3f} - {theta.max():.3f} m³/m³")
    print(f"Temperature range: {(T_soil.min()-273.15):.1f} - {(T_soil.max()-273.15):.1f} °C")
    print(f"Total precipitation: {np.sum(precip_series)*3600*24/len(t_forcing):.1f} mm")
    print(f"Total ET: {np.sum(ET_series)*3600*24/len(t_forcing):.1f} mm")

if __name__ == "__main__":
    run_soil_model()
