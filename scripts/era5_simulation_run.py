import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add project root to path
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from src.julesf.data.era5_data_reader import load_era5_for_jules
from src.julesf.coupler.jules_master import jules_master_coupler
from src.julesf.visualization.coupling_plots import visualize_jules_coupling

def create_era5_config(era5_drivers, simulation_weeks=52):
    """
    Create clean ERA5 forcing configuration
    
    Converts ERA5 data drivers into the format expected by fast_models_wrapper:
    external_drivers = {'era5_forcing': {variable: lambda t: value_at_time_t, ...}}
    """
    def estimate_humidity(t):
        """Estimate humidity from air temperature (60% RH assumption)"""
        Ta = era5_drivers['air_temperature'](t)
        # Simple humidity estimate
        q_star = 0.01 * np.exp(0.07 * (Ta - 273.15))
        return 0.6 * q_star
    
    # Map ERA5 data to model variables
    era5_forcing = {
        'Tair': era5_drivers['air_temperature'],        # K
        'Sdn': era5_drivers['shortwave_down'],          # W/m²
        'Ldn': era5_drivers['longwave_down'],           # W/m²
        'Q1': estimate_humidity,                        # kg/kg (estimated)
        'precip': era5_drivers['precipitation'],        # kg/m²/s
        'pressure': era5_drivers['surface_pressure'],   # Pa
        'co2': era5_drivers['co2_concentration'],       # ppm
    }
    
    print(f"✓ ERA5 config created with variables: {list(era5_forcing.keys())}")
    
    return {
        'weeks': simulation_weeks,
        'era5_forcing': era5_forcing,
        'data_source': 'ERA5',
        'time_bounds': era5_drivers.get('time_bounds', (0, simulation_weeks * 7 * 24))
    }

def run_jules_era5_simulation(data_file, metadata_file, weeks=52, output_dir='results', days_to_plot=30):
    """
    Run JULES simulation with ERA5 meteorological forcing
    """
    print("=== JULES Simulation with ERA5 Forcing ===")
    
    # 1. Load ERA5 data
    print("Loading ERA5 meteorological data...")
    era5_data = load_era5_for_jules(data_file, metadata_file, timestep_hours=0.5)
    
    print(f"Available ERA5 variables: {list(era5_data['drivers'].keys())}")
    
    # 2. Create configuration
    print("Creating JULES configuration...")
    config = create_era5_config(era5_data['drivers'], weeks)
    
    # 3. Check data availability
    max_time_hours = era5_data['drivers']['time_bounds'][1]
    max_weeks = int(max_time_hours / (7 * 24))
    
    print(f"ERA5 data covers {max_weeks} weeks ({max_time_hours:.0f} hours)")
    
    if weeks > max_weeks:
        print(f"⚠️  Warning: Requested {weeks} weeks but only {max_weeks} weeks available")
        weeks = max_weeks
        config['weeks'] = weeks
    
    print(f"Running simulation for {weeks} weeks...")
    
    # 4. Test ERA5 drivers before running simulation
    print("Testing ERA5 drivers...")
    era5_forcing = config['era5_forcing']
    test_times = [0, 24, 48]  # Hours
    
    for t in test_times:
        try:
            temp = era5_forcing['Tair'](t)
            sdn = era5_forcing['Sdn'](t)
            precip = era5_forcing['precip'](t)
            print(f"  t={t}h: Tair={temp:.1f}K, Sdn={sdn:.1f}W/m², Precip={precip:.6f}kg/m²/s")
        except Exception as e:
            print(f"  ❌ Driver test failed at t={t}h: {e}")
            raise
    
    # 5. Run JULES simulation
    try:
        print("\n" + "="*50)
        results = jules_master_coupler(
            weeks=weeks,
            triffid_init=None,     # Use defaults
            rothc_init=None,       # Use defaults  
            soil_init=None,        # Use defaults
            external_drivers=config  # Pass ERA5 configuration
        )
        print("="*50)
        
        print(f"✓ Simulation completed successfully!")
        print(f"✓ Forcing used: {results['metadata']['external_forcing']}")
        
        # 6. Generate coupling analysis plots and capture summary
        print("Generating coupling analysis plots…")
        summary = visualize_jules_coupling(
            results,
            save_path=os.path.join(output_dir, 'coupling_plots')
        )

        # 7. Save ERA5 data summary
        save_era5_summary(era5_data['summary'], output_dir)

        # 8. Plot forcing for first N days
        print(f"Plotting first {days_to_plot} days of forcing…")
        plot_era5_forcing(era5_data, output_dir, days_to_plot)

        return {
            'results': results,
            'era5_data': era5_data,
            'config': config,
            'summary': summary   # now defined
        }

    except Exception as e:
        print(f"❌ Simulation failed: {e}")
        raise

def save_era5_summary(era5_summary, output_dir):
    """Save summary of ERA5 meteorological data"""
    summary_file = os.path.join(output_dir, 'era5_summary.txt')
    
    with open(summary_file, 'w') as f:
        f.write("ERA5 Meteorological Data Summary\n")
        f.write("=" * 40 + "\n\n")
        
        for var, stats in era5_summary.items():
            f.write(f"{var}:\n")
            f.write(f"  Mean: {stats['mean']:.3f} {stats['units']}\n")
            f.write(f"  Std:  {stats['std']:.3f} {stats['units']}\n")
            f.write(f"  Min:  {stats['min']:.3f} {stats['units']}\n")
            f.write(f"  Max:  {stats['max']:.3f} {stats['units']}\n\n")
    
    print(f"✓ ERA5 summary saved to: {summary_file}")

def plot_era5_forcing(era5_data, output_dir, days_to_plot=30):
    """Plot ERA5 forcing data for visual inspection"""
    data = era5_data['data']
    time_hours = data['time_hours']
    
    # actual timestep in hours
    dt_h = time_hours[1] - time_hours[0]
    # number of steps for days_to_plot days
    steps = int(days_to_plot * 24 / dt_h)
    end_idx = min(len(time_hours), steps)
    
    time_days = time_hours[:end_idx] / 24
    
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle(f'ERA5 Forcing Data (First {days_to_plot} days)', fontsize=16)
    
    # Temperature
    if 'air_temperature' in data:
        axes[0, 0].plot(time_days, data['air_temperature'][:end_idx] - 273.15)
        axes[0, 0].set_title('Air Temperature')
        axes[0, 0].set_ylabel('Temperature (°C)')
        axes[0, 0].grid(True, alpha=0.3)
    
    # Radiation
    if 'net_shortwave' in data:
        axes[0, 1].plot(time_days, data['net_shortwave'][:end_idx], label='Net SW')
    if 'net_longwave' in data:
        axes[0, 1].plot(time_days, data['net_longwave'][:end_idx], label='Net LW')
    axes[0, 1].set_title('Net Radiation')
    axes[0, 1].set_ylabel('Radiation (W/m²)')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Precipitation
    if 'precipitation' in data:
        axes[1, 0].plot(time_days, data['precipitation'][:end_idx] * 3600)  # Convert to mm/h
        axes[1, 0].set_title('Precipitation')
        axes[1, 0].set_ylabel('Precipitation (mm/h)')
        axes[1, 0].grid(True, alpha=0.3)
    
    # Pressure
    if 'surface_pressure' in data:
        axes[1, 1].plot(time_days, data['surface_pressure'][:end_idx] / 100)  # Convert to hPa
        axes[1, 1].set_title('Surface Pressure')
        axes[1, 1].set_ylabel('Pressure (hPa)')
        axes[1, 1].grid(True, alpha=0.3)
    
    # CO2
    if 'co2_concentration' in data:
        axes[2, 0].plot(time_days, data['co2_concentration'][:end_idx])
        axes[2, 0].set_title('CO₂ Concentration')
        axes[2, 0].set_ylabel('CO₂ (ppm)')
        axes[2, 0].grid(True, alpha=0.3)
    
    # Downward radiation
    if 'shortwave_down' in data:
        axes[2, 1].plot(time_days, data['shortwave_down'][:end_idx], label='SW down')
    if 'longwave_down' in data:
        axes[2, 1].plot(time_days, data['longwave_down'][:end_idx], label='LW down')
    axes[2, 1].set_title('Downward Radiation')
    axes[2, 1].set_ylabel('Radiation (W/m²)')
    axes[2, 1].legend()
    axes[2, 1].grid(True, alpha=0.3)
    
    for ax in axes.flatten():
        ax.set_xlabel('Days')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'era5_forcing_data.png'), dpi=300)
    plt.close()

if __name__ == "__main__":
    # File paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, "data")
    
    data_file = os.path.join(data_dir, "era5_jules_met_ins_2002-2004.csv")
    metadata_file = os.path.join(data_dir, "era5_jules_met_ins_metadata_2002-2004.csv")
    
    # Check if files exist
    if not os.path.exists(data_file):
        print(f"❌ Data file not found: {data_file}")
        exit(1)
    
    if not os.path.exists(metadata_file):
        print(f"❌ Metadata file not found: {metadata_file}")
        exit(1)
    
    print(f"✓ Data file: {data_file}")
    print(f"✓ Metadata file: {metadata_file}")
    
    # Run simulation
    try:
        simulation_results = run_jules_era5_simulation(
            data_file=data_file,
            metadata_file=metadata_file,
            weeks=8,  
            output_dir="results/era5_simulation"
        )
        
        print("\n=== SIMULATION COMPLETE ===")
        print("✓ Results saved to: results/era5_simulation/")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

