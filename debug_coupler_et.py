"""
Debug the actual ET being used in the coupler vs our fixed version
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('src')

from julesf.data.era5_data_reader import load_era5_for_jules
from scripts.era5_simulation_run import create_era5_config
from julesf.coupler.fast_models_wrapper import run_fast_models_week, calculate_blaney_criddle_et
from julesf.soil.equations_moisture import infiltration_rate, evapotranspiration_extraction
from julesf.soil.parameters import SOIL_LAYERS, SOIL_PROPERTIES, VAN_GENUCHTEN, THERMAL_PROPERTIES

def debug_coupler_et():
    """Debug what ET is actually being used in the coupler"""
    print("="*60)
    print("COUPLER ET DIAGNOSTIC")
    print("="*60)
    
    # Load ERA5 data
    data_file = "data/era5_jules_met_ins_2002-2004.csv"
    metadata_file = "data/era5_jules_met_ins_metadata_2002-2004.csv"
    
    era5_data = load_era5_for_jules(data_file, metadata_file, timestep_hours=0.5)
    config = create_era5_config(era5_data['drivers'], simulation_weeks=1)
    
    # Test what ET the coupler actually creates
    print("Testing fast_models_wrapper ET setup...")
    
    # Run the fast models wrapper to see what soil drivers it creates
    try:
        fast_results = run_fast_models_week(
            nu_cover=0.8,
            LAI_total=5.0,
            soil_initial=None,
            week_num=0,
            external_drivers=config
        )
        
        # Extract the soil drivers that were actually used
        print("✓ Fast models completed successfully")
        
        # Test the ET function directly from the wrapper
        era5_forcing = config['era5_forcing']
        test_times = [0, 24, 48, 72]  # Hours
        
        print("\nTesting Blaney-Criddle ET function:")
        for t in test_times:
            temp_K = era5_forcing['Tair'](t)
            et_blaney = calculate_blaney_criddle_et(temp_K)
            
            print(f"  t={t}h: T={temp_K:.1f}K → ET={et_blaney*1e6:.2f} μm/s = {et_blaney*3600*1000:.3f} mm/hr")
        
        # Compare with what our debug script calculated
        print("\nComparing with standalone calculation:")
        soil_drivers = {
            'evapotranspiration': lambda t: calculate_blaney_criddle_et(era5_forcing['Tair'](t)),
        }
        
        for t in test_times:
            et_standalone = soil_drivers['evapotranspiration'](t)
            print(f"  t={t}h: Standalone ET={et_standalone*1e6:.2f} μm/s = {et_standalone*3600*1000:.3f} mm/hr")
            
    except Exception as e:
        print(f"❌ Fast models failed: {e}")
        import traceback
        traceback.print_exc()

def check_fast_wrapper_source():
    """Check the actual source code for conflicting ET assignments"""
    print("\n" + "="*60)
    print("CHECKING FAST_MODELS_WRAPPER SOURCE")
    print("="*60)
    
    wrapper_file = "src/julesf/coupler/fast_models_wrapper.py"
    
    try:
        with open(wrapper_file, 'r') as f:
            lines = f.readlines()
        
        et_assignments = []
        for i, line in enumerate(lines):
            if 'evapotranspiration' in line and '=' in line and 'lambda' in line:
                et_assignments.append((i+1, line.strip()))
        
        print(f"Found {len(et_assignments)} ET assignments in fast_models_wrapper.py:")
        
        for line_num, line in et_assignments:
            print(f"  Line {line_num}: {line}")
        
        if len(et_assignments) > 1:
            print("\n❌ MULTIPLE ET ASSIGNMENTS FOUND!")
            print("   The last assignment will override earlier ones.")
            print("   Check which one is actually being used.")
            
    except Exception as e:
        print(f"❌ Could not read wrapper file: {e}")

def test_soil_with_coupler_drivers():
    """Test soil model with the exact drivers from the coupler"""
    print("\n" + "="*60)
    print("TESTING SOIL WITH COUPLER DRIVERS")
    print("="*60)
    
    # Load data
    data_file = "data/era5_jules_met_ins_2002-2004.csv"
    metadata_file = "data/era5_jules_met_ins_metadata_2002-2004.csv"
    
    era5_data = load_era5_for_jules(data_file, metadata_file, timestep_hours=0.5)
    config = create_era5_config(era5_data['drivers'], simulation_weeks=1)
    
    # Create the EXACT same drivers that fast_models_wrapper creates
    era5 = config['era5_forcing']
    
    # This should match what's in fast_models_wrapper.py line ~140
    soil_drivers_coupler = {
        'precipitation': lambda t: max(0, era5['precip'](t)),
        'evapotranspiration': lambda t: calculate_blaney_criddle_et(era5['Tair'](t)),
    }
    
    # Test water balance with coupler drivers
    params = {**SOIL_LAYERS, **SOIL_PROPERTIES, **VAN_GENUCHTEN, **THERMAL_PROPERTIES}
    times = np.arange(0, 24, 1)  # 1 day, hourly
    
    precip_values = []
    et_values = []
    
    for t in times:
        precip = soil_drivers_coupler['precipitation'](t)
        et = soil_drivers_coupler['evapotranspiration'](t)
        
        precip_values.append(precip)
        et_values.append(et)
    
    precip_total = np.sum(precip_values) * 3600  # mm/day
    et_total = np.sum(et_values) * 3600  # mm/day
    
    print(f"Water balance with coupler drivers (1 day):")
    print(f"  Precipitation: {precip_total:.3f} mm/day")
    print(f"  ET:            {et_total:.3f} mm/day")
    print(f"  Net:           {precip_total - et_total:.3f} mm/day")
    
    if precip_total > et_total:
        print("❌ PRECIPITATION > ET → explains monotonic moisture increase!")
    else:
        print("✓ ET > PRECIPITATION → should see moisture decrease")

if __name__ == "__main__":
    debug_coupler_et()
    check_fast_wrapper_source()
    test_soil_with_coupler_drivers()