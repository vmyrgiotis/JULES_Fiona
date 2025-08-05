"""
Debug soil water balance - check precipitation vs ET magnitudes
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('src')

from julesf.data.era5_data_reader import load_era5_for_jules
from scripts.era5_simulation_run import create_era5_config
from julesf.soil.equations_moisture import infiltration_rate, evapotranspiration_extraction
from julesf.soil.parameters import SOIL_LAYERS, SOIL_PROPERTIES, VAN_GENUCHTEN, THERMAL_PROPERTIES

def calculate_blaney_criddle_et(temperature_K, latitude_deg=52.0):
    """
    Calculate potential evapotranspiration using Blaney-Criddle method
    PET₀ = p(0.457 T_mean + 8.128) [mm/day]
    """
    # Convert K to °C
    T_celsius = temperature_K - 273.15
    
    # Simplified p factor for mid-latitudes (varies seasonally, but use average)
    # For UK (52°N): p ≈ 0.27 in summer, 0.23 in winter, average ≈ 0.25
    p_factor = 0.25
    
    # Blaney-Criddle formula (mm/day)
    PET_mm_day = p_factor * (0.457 * T_celsius + 8.128)
    
    # Ensure non-negative
    PET_mm_day = max(0, PET_mm_day)
    
    # Convert mm/day to kg/m²/s
    # 1 mm/day = 1 kg/m²/day = 1/(24×3600) kg/m²/s
    PET_kg_m2_s = PET_mm_day / (24 * 3600)
    
    return PET_kg_m2_s

def debug_water_balance():
    """Debug the water balance components"""
    print("="*60)
    print("SOIL WATER BALANCE DIAGNOSTICS")
    print("="*60)
    
    # Load ERA5 data
    data_file = "data/era5_jules_met_ins_2002-2004.csv"
    metadata_file = "data/era5_jules_met_ins_metadata_2002-2004.csv"
    
    era5_data = load_era5_for_jules(data_file, metadata_file, timestep_hours=0.5)
    config = create_era5_config(era5_data['drivers'], simulation_weeks=1)  # FIXED: simulation_weeks
    era5_forcing = config['era5_forcing']
    
    # Setup soil drivers with CORRECT Blaney-Criddle
    soil_drivers = {
        'precipitation': lambda t: max(0, era5_forcing['precip'](t)),
        'evapotranspiration': lambda t: calculate_blaney_criddle_et(era5_forcing['Tair'](t)),  # CORRECT FORMULA
    }
    
    # Get parameters
    params = {**SOIL_LAYERS, **SOIL_PROPERTIES, **VAN_GENUCHTEN, **THERMAL_PROPERTIES}
    
    # Test over 1 week
    times = np.arange(0, 7*24, 1)  # Every hour for 7 days
    precip_values = []
    et_values = []
    infiltration_values = []
    et_extracted_values = []
    
    # Initial soil moisture (moderately dry)
    theta = np.array([0.25, 0.23, 0.21, 0.20])  # 4 layers
    
    print("Testing water balance components...")
    
    for t in times:
        # Test precipitation and infiltration
        precip = soil_drivers['precipitation'](t)
        infiltration = infiltration_rate(params, soil_drivers, t)
        
        # Test ET
        et_potential = soil_drivers['evapotranspiration'](t)
        et_extracted = evapotranspiration_extraction(theta, params, soil_drivers, t)
        et_total = np.sum(et_extracted)
        
        precip_values.append(precip)
        et_values.append(et_potential)
        infiltration_values.append(infiltration)
        et_extracted_values.append(et_total)
    
    # Convert to arrays
    precip_values = np.array(precip_values)
    et_values = np.array(et_values)
    infiltration_values = np.array(infiltration_values)
    et_extracted_values = np.array(et_extracted_values)
    
    # Analysis
    print(f"\nWATER BALANCE ANALYSIS (7 days):")
    print(f"Precipitation:")
    print(f"  Mean: {np.mean(precip_values)*1e6:.2f} μm/s = {np.mean(precip_values)*3600*1000:.3f} mm/hr")
    print(f"  Max:  {np.max(precip_values)*1e6:.2f} μm/s = {np.max(precip_values)*3600*1000:.3f} mm/hr")
    print(f"  Total: {np.sum(precip_values)*3600:.3f} mm/week")
    
    print(f"\nInfiltration (limited by K_sat):")
    print(f"  Mean: {np.mean(infiltration_values)*1e6:.2f} μm/s = {np.mean(infiltration_values)*3600*1000:.3f} mm/hr")
    print(f"  Max:  {np.max(infiltration_values)*1e6:.2f} μm/s = {np.max(infiltration_values)*3600*1000:.3f} mm/hr")
    print(f"  Total: {np.sum(infiltration_values)*3600:.3f} mm/week")
    
    print(f"\nPotential ET:")
    print(f"  Mean: {np.mean(et_values)*1e6:.2f} μm/s = {np.mean(et_values)*3600*1000:.3f} mm/hr")
    print(f"  Max:  {np.max(et_values)*1e6:.2f} μm/s = {np.max(et_values)*3600*1000:.3f} mm/hr")
    print(f"  Total: {np.sum(et_values)*3600:.3f} mm/week")
    
    print(f"\nActual ET extracted:")
    print(f"  Mean: {np.mean(et_extracted_values)*1e6:.2f} μm/s = {np.mean(et_extracted_values)*3600*1000:.3f} mm/hr")
    print(f"  Max:  {np.max(et_extracted_values)*1e6:.2f} μm/s = {np.max(et_extracted_values)*3600*1000:.3f} mm/hr")
    print(f"  Total: {np.sum(et_extracted_values)*3600:.3f} mm/week")
    
    print(f"\nNET WATER BALANCE:")
    net_input = np.sum(infiltration_values) - np.sum(et_extracted_values)
    print(f"  Net input: {net_input*3600:.3f} mm/week")
    print(f"  Net input per day: {net_input*3600/7:.3f} mm/day")
    
    # Check K_sat limit
    K_sat = params['K_sat'] * params.get('rho_water', 1000.0)  # Convert to kg/m²/s
    print(f"\nHYDRAULIC LIMITS:")
    print(f"  K_sat: {K_sat*1e6:.2f} μm/s = {K_sat*3600*1000:.3f} mm/hr")
    
    if np.max(precip_values) > K_sat:
        print(f"  ⚠️  Precipitation sometimes exceeds K_sat!")
    else:
        print(f"  ✓ Precipitation always below K_sat")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    axes[0,0].plot(times/24, precip_values*3600*1000, 'b-', label='Precipitation')
    axes[0,0].plot(times/24, infiltration_values*3600*1000, 'b--', label='Infiltration (actual)')
    axes[0,0].set_ylabel('mm/hr')
    axes[0,0].set_title('Water Input')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    axes[0,1].plot(times/24, et_values*3600*1000, 'r-', label='Potential ET')
    axes[0,1].plot(times/24, et_extracted_values*3600*1000, 'r--', label='Actual ET')
    axes[0,1].set_ylabel('mm/hr')
    axes[0,1].set_title('Evapotranspiration')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # Cumulative
    axes[1,0].plot(times/24, np.cumsum(infiltration_values)*3600, 'b-', label='Cumulative infiltration')
    axes[1,0].plot(times/24, np.cumsum(et_extracted_values)*3600, 'r-', label='Cumulative ET')
    axes[1,0].set_ylabel('mm')
    axes[1,0].set_xlabel('Days')
    axes[1,0].set_title('Cumulative Water Balance')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)
    
    # Net balance
    net_hourly = (infiltration_values - et_extracted_values) * 3600
    axes[1,1].plot(times/24, np.cumsum(net_hourly), 'k-', linewidth=2)
    axes[1,1].set_ylabel('mm')
    axes[1,1].set_xlabel('Days')
    axes[1,1].set_title('Cumulative Net Water Input')
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('debug_water_balance.png', dpi=150, bbox_inches='tight')
    print(f"\n✓ Water balance plot saved: debug_water_balance.png")
    
    # Final diagnosis
    print(f"\n" + "="*60)
    print("DIAGNOSIS:")
    
    if net_input > 0:
        print(f"❌ NET POSITIVE WATER INPUT: {net_input*3600:.3f} mm/week")
        print(f"   This explains the monotonic moisture increase!")
        
        if np.mean(et_extracted_values) < np.mean(infiltration_values) * 0.1:
            print(f"   💡 ET is too low compared to precipitation")
            print(f"   💡 Ratio ET/Precip: {np.mean(et_extracted_values)/np.mean(infiltration_values):.3f}")
        
        if np.max(et_extracted_values) == 0:
            print(f"   💡 ET extraction is zero - check ET function!")
            
    else:
        print(f"✓ Net water balance looks reasonable")
    
    return {
        'precip': precip_values,
        'et': et_extracted_values,
        'infiltration': infiltration_values,
        'net_input': net_input
    }

if __name__ == "__main__":
    results = debug_water_balance()