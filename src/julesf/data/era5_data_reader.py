import pandas as pd
import numpy as np
from datetime import datetime
import os

class ERA5DataReader:
    
    def __init__(self, data_file, metadata_file):
        self.data_file = data_file
        self.metadata_file = metadata_file
        self.data = None
        self.metadata = None
        
        self.var_mapping = {
            # Temperature
            't2m': 'air_temperature',  # K
            
            # Radiation components for EBM
            'ssr': 'net_shortwave',    # W/m²
            'str': 'net_longwave',     # W/m²
            'ssrd': 'shortwave_down',  # W/m²
            'strd': 'longwave_down',   # W/m²
            
            # Hydrology
            'tp': 'precipitation',     # kg/m²/s
            
            # Atmospheric conditions
            'sp': 'surface_pressure',  # Pa
            'co2': 'co2_concentration' # ppm
        }
    
    def load_data(self):
        """Load and process ERA5 data"""
        # Load metadata
        self.metadata = pd.read_csv(self.metadata_file)
        
        # Load main data
        self.data = pd.read_csv(self.data_file)
        
        # Convert time column to datetime
        self.data['datetime'] = pd.to_datetime(self.data['valid_time'], 
                                               format='%d/%m/%Y %H:%M')
        
        # Set datetime as index
        self.data.set_index('datetime', inplace=True)
        
        print(f"Loaded {len(self.data)} records from {self.data.index[0]} to {self.data.index[-1]}")
        
        return self
    
    def extract_jules_variables(self):
        if self.data is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        
        # Extract only the variables we need
        required_vars = list(self.var_mapping.keys())
        
        # Check which variables are available
        available_vars = [var for var in required_vars if var in self.data.columns]
        missing_vars = [var for var in required_vars if var not in self.data.columns]
        
        if missing_vars:
            print(f"Warning: Missing variables: {missing_vars}")
        
        # Extract available variables
        jules_data = {}
        for era5_var, jules_var in self.var_mapping.items():
            if era5_var in available_vars:
                jules_data[jules_var] = self.data[era5_var].values
                print(f"Extracted {era5_var} -> {jules_var}: "
                      f"min={jules_data[jules_var].min():.3f}, "
                      f"max={jules_data[jules_var].max():.3f}")
        
        # Add time information
        jules_data['time'] = self.data.index
        jules_data['time_hours'] = np.arange(len(self.data)) * 1.0  # Hourly data
        
        return jules_data
    
    def create_model_drivers(self, jules_data):
        """Create interpolation functions for model drivers"""
        from scipy.interpolate import interp1d
        
        drivers = {}
        time_hours = jules_data['time_hours']
        
        # Create interpolation functions for each variable
        for var_name, var_data in jules_data.items():
            if var_name not in ['time', 'time_hours']:
                # Create interpolation function
                drivers[var_name] = interp1d(
                    time_hours, var_data, 
                    kind='linear', 
                    bounds_error=False, 
                    fill_value='extrapolate'
                )
        
        # Store time bounds for reference
        drivers['time_bounds'] = (time_hours[0], time_hours[-1])
        drivers['time_step'] = 1.0  # hourly
        
        return drivers
    
    def get_summary_statistics(self):
        if self.data is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        
        summary = {}
        for era5_var, jules_var in self.var_mapping.items():
            if era5_var in self.data.columns:
                data = self.data[era5_var]
                summary[jules_var] = {
                    'min': data.min(),
                    'max': data.max(),
                    'mean': data.mean(),
                    'std': data.std(),
                    'units': self._get_units(era5_var)
                }
        
        return summary
    
    def _get_units(self, variable):
        """Get units for a variable from metadata"""
        if self.metadata is not None:
            unit_row = self.metadata[self.metadata['variable'] == variable]
            if not unit_row.empty:
                return unit_row['units'].iloc[0]
        return 'unknown'
    
    def resample_to_timestep(self, jules_data, target_timestep_hours=0.5):
        """Resample data to target timestep"""
        from scipy.interpolate import interp1d
        
        original_hours = jules_data['time_hours']
        
        # Create new time array
        new_time_hours = np.arange(
            original_hours[0], 
            original_hours[-1] + target_timestep_hours, 
            target_timestep_hours
        )
        
        # Resample each variable
        resampled_data = {
            'time_hours': new_time_hours
        }
        
        for var_name, var_data in jules_data.items():
            if var_name not in ['time', 'time_hours']:
                interp_func = interp1d(
                    original_hours, var_data,
                    kind='linear',
                    bounds_error=False,
                    fill_value='extrapolate'
                )
                resampled_data[var_name] = interp_func(new_time_hours)
        
        print(f"Resampled from {len(original_hours)} to {len(new_time_hours)} points "
              f"(timestep: {target_timestep_hours}h)")
        
        return resampled_data

def load_era5_for_jules(data_file, metadata_file, timestep_hours=0.5):
    """
    Convenience function to load ERA5 data for JULES
    
    Args:
        data_file: Path to ERA5 data CSV
        metadata_file: Path to metadata CSV  
        timestep_hours: Target timestep in hours
        
    Returns:
        Dictionary with model drivers and interpolation functions
    """
    reader = ERA5DataReader(data_file, metadata_file)
    reader.load_data()
    
    # Extract JULES variables
    jules_data = reader.extract_jules_variables()
    
    # Resample to target timestep
    if timestep_hours != 1.0:
        jules_data = reader.resample_to_timestep(jules_data, timestep_hours)
    
    # Create model drivers
    drivers = reader.create_model_drivers(jules_data)
    
    # Add summary statistics
    summary = reader.get_summary_statistics()
    
    return {
        'drivers': drivers,
        'data': jules_data,
        'summary': summary,
        'reader': reader
    }

if __name__ == "__main__":
    # Test the reader - fix file paths
    import os
    
    # Get the correct paths relative to the script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(script_dir, "..", "..", "..")
    data_dir = os.path.join(project_root, "data")
    
    data_file = os.path.join(data_dir, "era5_jules_met_ins_2002-2004.csv")
    metadata_file = os.path.join(data_dir, "era5_jules_met_ins_metadata_2002-2004.csv")
    
    # Check if files exist
    if not os.path.exists(data_file):
        print(f"Error: Data file not found at {data_file}")
        print(f"Looking in directory: {data_dir}")
        if os.path.exists(data_dir):
            print(f"Available files: {os.listdir(data_dir)}")
        exit(1)
    
    if not os.path.exists(metadata_file):
        print(f"Error: Metadata file not found at {metadata_file}")
        exit(1)
    
    print(f"Loading data from: {data_file}")
    print(f"Loading metadata from: {metadata_file}")
    
    result = load_era5_for_jules(data_file, metadata_file)
    
    print("\n=== SUMMARY STATISTICS ===")
    for var, stats in result['summary'].items():
        print(f"{var}: {stats['mean']:.3f} ± {stats['std']:.3f} {stats['units']} "
              f"[{stats['min']:.3f}, {stats['max']:.3f}]")