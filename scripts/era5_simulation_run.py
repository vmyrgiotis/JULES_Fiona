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
    """
    def estimate_humidity(t):
        Ta = era5_drivers['air_temperature'](t)
        q_star = 0.01 * np.exp(0.07 * (Ta - 273.15))
        return 0.6 * q_star

    era5_forcing = {
        'Tair':            era5_drivers['air_temperature'],        # K
        'Sdn':             era5_drivers['shortwave_down'],         # W/m²
        'Ldn':             era5_drivers['longwave_down'],          # W/m²
        'Q1':              estimate_humidity,                      # 1
        'precip':          era5_drivers['precipitation'],          # kg/m²/s (no change)
        # FIX: convert ET mm/h → m/h → m/s → kg/m²/s
        'evapotranspiration': lambda t: era5_drivers['evapotranspiration'](t) 
                                     / 1000.0  # mm→m
                                     / 3600.0, # h→s
        'PAR':             era5_drivers['PAR'],                    # W/m²
        'pressure':        era5_drivers['surface_pressure'],       # Pa
        'co2':             era5_drivers['co2_concentration'],      # ppm
    }

    print(f"✓ ERA5 config created with variables: {list(era5_forcing.keys())}")

    return {
        'weeks':       simulation_weeks,
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

        print("Creating high-resolution plots...")
        plot_high_resolution_results(results, era5_data, output_dir, weeks_to_plot=2)
    

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

def plot_high_resolution_results(results, era5_data, output_dir, weeks_to_plot=2):
    """Plot high-resolution (0.5-hourly) results from the first few weeks"""
    print(f"Plotting first {weeks_to_plot} weeks of high-resolution results...")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Check what's available in results
    print(f"Available result keys: {list(results.keys())}")
    
    if 'weekly_results' not in results:
        print("❌ No weekly_results found in results!")
        return
    
    weekly_data = results['weekly_results']
    if len(weekly_data) < weeks_to_plot:
        weeks_to_plot = len(weekly_data)
        print(f"⚠️  Only {len(weekly_data)} weeks available, plotting {weeks_to_plot}")
    
    # FIXED: Build time array based on ACTUAL data length, not calculated length
    all_theta = []
    all_T_soil = []
    all_T_surface = []
    
    for week in range(min(weeks_to_plot, len(weekly_data))):
        week_data = weekly_data[week]
        
        # Soil data
        if 'soil_theta' in week_data:
            theta_week = week_data['soil_theta']  # Shape: (n_layers, n_times)
            if len(all_theta) == 0:
                all_theta = theta_week
            else:
                all_theta = np.concatenate([all_theta, theta_week], axis=1)
        
        if 'soil_T' in week_data:
            T_soil_week = week_data['soil_T']
            if len(all_T_soil) == 0:
                all_T_soil = T_soil_week
            else:
                all_T_soil = np.concatenate([all_T_soil, T_soil_week], axis=1)
        
        if 'T_surface' in week_data:
            Ts_week = week_data['T_surface']
            all_T_surface.extend(Ts_week)
    
    # Convert to arrays
    all_theta = np.array(all_theta) if len(all_theta) > 0 else None
    all_T_soil = np.array(all_T_soil) if len(all_T_soil) > 0 else None
    all_T_surface = np.array(all_T_surface) if len(all_T_surface) > 0 else None
    
    # FIXED: Create time array based on ACTUAL data length
    if all_theta is not None:
        n_time_points = all_theta.shape[1]
    elif all_T_surface is not None:
        n_time_points = len(all_T_surface)
    else:
        print("❌ No data available for plotting!")
        return
    
    # Create time array with correct length
    dt_hours = 0.5
    time_hours = np.arange(0, n_time_points * dt_hours, dt_hours)
    time_days = time_hours / 24
    
    # Ensure time and data arrays have same length
    if len(time_days) != n_time_points:
        print(f"⚠️  Adjusting time array: {len(time_days)} -> {n_time_points}")
        time_days = time_days[:n_time_points]
    
    print(f"Plotting data shapes:")
    print(f"  Time points: {len(time_days)}")
    print(f"  Soil moisture: {all_theta.shape if all_theta is not None else 'None'}")
    print(f"  Soil temperature: {all_T_soil.shape if all_T_soil is not None else 'None'}")
    print(f"  Surface temperature: {len(all_T_surface) if all_T_surface is not None else 'None'}")
    
    # CREATE THE PLOTS
    fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    fig.suptitle(f'High-Resolution JULES Results (First {weeks_to_plot} weeks)', fontsize=14)
    
    # ROW 1: SOIL MOISTURE & TEMPERATURE
    if all_theta is not None:
        for layer in range(all_theta.shape[0]):
            axes[0, 0].plot(time_days, all_theta[layer, :], label=f'Layer {layer+1}')
        axes[0, 0].set_title('Soil Moisture (0.5-hourly)')
        axes[0, 0].set_ylabel('θ (m³/m³)')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
    
    if all_T_soil is not None:
        for layer in range(all_T_soil.shape[0]):
            axes[0, 1].plot(time_days, all_T_soil[layer, :] - 273.15, label=f'Layer {layer+1}')
        axes[0, 1].set_title('Soil Temperature (0.5-hourly)')
        axes[0, 1].set_ylabel('T (°C)')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
    
    # ROW 2: SURFACE TEMPERATURE & METEOROLOGY
    if all_T_surface is not None:
        # Make sure lengths match
        surface_time = time_days[:len(all_T_surface)]
        axes[1, 1].plot(surface_time, np.array(all_T_surface) - 273.15, 'r-', linewidth=2)
        axes[1, 1].set_title('Surface Temperature (EBM)')
        axes[1, 1].set_ylabel('T_surface (°C)')
        axes[1, 1].grid(True, alpha=0.3)
    
    # ERA5 forcing data - use same time length as data
    drivers = era5_data['drivers']
    actual_time_hours = time_hours[:len(time_days)]
    
    try:
        # now drivers[...] return kg/m²/s → *3600 → mm/h
        precip_series = [drivers['precipitation'](t) * 3600.0 for t in actual_time_hours]
        et_series     = [drivers['evapotranspiration'](t) * 3600.0 for t in actual_time_hours]
        
        axes[1,0].bar(time_days, precip_series, width=0.01, alpha=0.6, color='blue', label='Precipitation')
        ax_twin = axes[1,0].twinx()
        ax_twin.plot(time_days, et_series, 'r-', label='ET')
        axes[1, 0].set_ylabel('Precipitation (mm/hr)', color='blue')
        ax_twin.set_ylabel('ET (mm/hr)', color='red')
        axes[1, 0].set_title('Water Fluxes')
        
        # ROW 3: NET WATER BALANCE & TEMPERATURE
        net_water = np.array(precip_series) - np.array(et_series)
        axes[2, 0].plot(time_days, net_water, 'purple', linewidth=2)
        axes[2, 0].axhline(y=0, color='black', linestyle='--', alpha=0.5)
        axes[2, 0].set_title('Net Water Input (Precip - ET)')
        axes[2, 0].set_ylabel('Net flux (mm/hr)')
        axes[2, 0].grid(True, alpha=0.3)
        
        temp_series = [drivers['air_temperature'](t) - 273.15 for t in actual_time_hours]
        axes[2, 1].plot(time_days, temp_series, 'k-', linewidth=2)
        axes[2, 1].set_title('Air Temperature')
        axes[2, 1].set_ylabel('T_air (°C)')
        axes[2, 1].grid(True, alpha=0.3)
        
    except Exception as e:
        print(f"⚠️  Warning: Could not plot forcing data: {e}")
        # Fill with placeholder text
        axes[1, 1].text(0.5, 0.5, 'Forcing data\nplotting failed', 
                        ha='center', va='center', transform=axes[1, 1].transAxes)
        axes[2, 0].text(0.5, 0.5, 'Net water balance\nplotting failed', 
                        ha='center', va='center', transform=axes[2, 0].transAxes)
        axes[2, 1].text(0.5, 0.5, 'Temperature\nplotting failed', 
                        ha='center', va='center', transform=axes[2, 1].transAxes)
    
    # Set x-labels
    for i in range(3):
        for j in range(2):
            axes[i, j].set_xlabel('Time (days)')
    
    plt.tight_layout()
    plot_path = os.path.join(output_dir, 'jules_high_resolution_results.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Hourly-resolution soil dynamics plot saved: {plot_path}")
    
    if os.path.exists(plot_path):
        file_size = os.path.getsize(plot_path)
        # print(f"File confirmed: {file_size} bytes")
    else:
        print(f"❌ Plot file NOT created at: {plot_path}")

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, "data")
    
    data_file = os.path.join(data_dir, "era5_jules_met_ins_2002-2016.csv")
    metadata_file = os.path.join(data_dir, "era5_jules_met_ins_metadata.csv")
    
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

