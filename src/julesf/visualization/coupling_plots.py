import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import os

def create_output_directory(base_path='plots'):
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    return base_path

def plot_vegetation_dynamics(results, save_path=None):
    """Plot vegetation dynamics (TRIFFID model outputs)"""
    weeks = results['triffid'].shape[1] - 1
    time = np.arange(weeks + 1)
    
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('TRIFFID Vegetation Dynamics', fontsize=16)
    
    # Plot leaf biomass
    axs[0, 0].plot(time, results['weekly']['Lb_tree'], 'g-', label='Tree')
    axs[0, 0].plot(time, results['weekly']['Lb_grass'], 'y-', label='Grass')
    axs[0, 0].set_title('Leaf Biomass')
    axs[0, 0].set_ylabel('Leaf Biomass (kg C/m²)')
    axs[0, 0].legend()
    axs[0, 0].grid(True, alpha=0.3)
    
    # Plot vegetation coverage
    axs[0, 1].plot(time, results['weekly']['nu_tree'], 'g-', label='Tree')
    axs[0, 1].plot(time, results['weekly']['nu_grass'], 'y-', label='Grass')
    axs[0, 1].plot(time, results['weekly']['nu_total'], 'k--', label='Total')
    axs[0, 1].set_title('Vegetation Coverage')
    axs[0, 1].set_ylabel('Fraction')
    axs[0, 1].legend()
    axs[0, 1].grid(True, alpha=0.3)
    axs[0, 1].set_ylim(0, 1.05)
    
    # Plot LAI
    axs[1, 0].plot(time, results['weekly']['LAI_total'], 'g-')
    axs[1, 0].set_title('Leaf Area Index')
    axs[1, 0].set_ylabel('LAI (m²/m²)')
    axs[1, 0].grid(True, alpha=0.3)
    
    # Plot NPP
    npp_time = np.arange(weeks)
    axs[1, 1].plot(npp_time, results['weekly']['NPP'], 'r-')
    axs[1, 1].set_title('Net Primary Productivity')
    axs[1, 1].set_ylabel('NPP (kg C/m²/week)')
    axs[1, 1].grid(True, alpha=0.3)
    
    for ax in axs.flatten():
        ax.set_xlabel('Week')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(os.path.join(save_path, 'vegetation_dynamics.png'), dpi=300)
        
    return fig

def plot_soil_carbon_dynamics(results, save_path=None):
    """Plot soil carbon dynamics (RothC model outputs)"""
    weeks = results['rothc'].shape[1] - 1
    time = np.arange(weeks + 1)
    
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('RothC Soil Carbon Dynamics', fontsize=16)
    
    # Plot individual carbon pools
    axs[0, 0].plot(time, results['weekly']['C_dpm'], 'r-', label='DPM')
    axs[0, 0].plot(time, results['weekly']['C_rpm'], 'g-', label='RPM')
    axs[0, 0].plot(time, results['weekly']['C_bio'], 'b-', label='BIO')
    axs[0, 0].plot(time, results['weekly']['C_hum'], 'k-', label='HUM')
    axs[0, 0].set_title('Soil Carbon Pools')
    axs[0, 0].set_ylabel('Carbon (kg C/m²)')
    axs[0, 0].legend()
    axs[0, 0].grid(True, alpha=0.3)
    
    # Plot total soil carbon
    axs[0, 1].plot(time, results['weekly']['soil_C_total'], 'k-')
    axs[0, 1].set_title('Total Soil Carbon')
    axs[0, 1].set_ylabel('Carbon (kg C/m²)')
    axs[0, 1].grid(True, alpha=0.3)
    
    # Plot litterfall
    axs[1, 0].plot(time, results['weekly']['litterfall'], 'm-')
    axs[1, 0].set_title('Litterfall')
    axs[1, 0].set_ylabel('Litterfall (kg C/m²/day)')
    axs[1, 0].grid(True, alpha=0.3)
    
    # Stacked area plot of carbon pools
    carbon_data = np.vstack([
        results['weekly']['C_dpm'],
        results['weekly']['C_rpm'],
        results['weekly']['C_bio'],
        results['weekly']['C_hum']
    ])
    axs[1, 1].stackplot(time, carbon_data, 
                        labels=['DPM', 'RPM', 'BIO', 'HUM'],
                        colors=['r', 'g', 'b', 'k'], alpha=0.6)
    axs[1, 1].set_title('Soil Carbon Pool Composition')
    axs[1, 1].set_ylabel('Carbon (kg C/m²)')
    axs[1, 1].legend(loc='upper left')
    axs[1, 1].grid(True, alpha=0.3)
    
    for ax in axs.flatten():
        ax.set_xlabel('Week')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(os.path.join(save_path, 'soil_carbon_dynamics.png'), dpi=300)
        
    return fig

def plot_soil_hydrothermal(results, save_path=None):
    """Plot soil moisture and temperature dynamics"""
    weeks = len(results['weekly']['soil_moisture_mean'])
    time = np.arange(weeks)
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 8))
    fig.suptitle('Soil Hydrothermal Conditions', fontsize=16)
    
    # Plot soil moisture
    axs[0].plot(time, results['weekly']['soil_moisture_mean'], 'b-')
    axs[0].set_title('Mean Soil Moisture')
    axs[0].set_ylabel('Soil Moisture (m³/m³)')
    axs[0].grid(True, alpha=0.3)
    
    # Plot soil temperature
    axs[1].plot(time, results['weekly']['soil_temp_mean'] - 273.15, 'r-')  # Convert K to °C
    axs[1].set_title('Mean Soil Temperature')
    axs[1].set_ylabel('Temperature (°C)')
    axs[1].grid(True, alpha=0.3)
    
    for ax in axs:
        ax.set_xlabel('Week')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(os.path.join(save_path, 'soil_hydrothermal.png'), dpi=300)
        
    return fig

def plot_photosynthesis_components(results, save_path=None):
    """Plot photosynthesis components (Pi_G and R_p)"""
    weeks = len(results['weekly']['Pi_G'])
    time = np.arange(weeks)
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 8))
    fig.suptitle('Photosynthesis Components', fontsize=16)
    
    # Plot gross photosynthesis
    axs[0].plot(time, results['weekly']['Pi_G'], 'g-')
    axs[0].set_title('Gross Photosynthesis (Pi_G)')
    axs[0].set_ylabel('Pi_G (μmol CO₂/m²/s)')
    axs[0].grid(True, alpha=0.3)
    
    # Plot plant respiration
    axs[1].plot(time, results['weekly']['R_p'], 'r-')
    axs[1].set_title('Plant Respiration (R_p)')
    axs[1].set_ylabel('R_p (μmol CO₂/m²/s)')
    axs[1].grid(True, alpha=0.3)
    
    for ax in axs:
        ax.set_xlabel('Week')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(os.path.join(save_path, 'photosynthesis_components.png'), dpi=300)
        
    return fig

def plot_model_coupling(results, save_path=None):
    """Plot key coupling variables between models"""
    weeks = results['triffid'].shape[1] - 1
    time = np.arange(weeks)
    
    # Create figure with GridSpec for flexible subplot layout
    fig = plt.figure(figsize=(14, 12))
    gs = gridspec.GridSpec(3, 2)
    
    fig.suptitle('Model Coupling Variables', fontsize=16)
    
    # NPP (Physiology → TRIFFID)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(time, results['weekly']['NPP'], 'g-')
    ax1.set_title('NPP (Physiology → TRIFFID)')
    ax1.set_ylabel('NPP (kg C/m²/week)')
    ax1.grid(True, alpha=0.3)
    
    # Litterfall (TRIFFID → RothC)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(time, results['weekly']['litterfall'][1:], 'm-')  # Skip initial value
    ax2.set_title('Litterfall (TRIFFID → RothC)')
    ax2.set_ylabel('Litterfall (kg C/m²/day)')
    ax2.grid(True, alpha=0.3)
    
    # Vegetation Cover (TRIFFID → EBM/Soil)
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(time, results['weekly']['nu_total'][1:], 'b-')  # Skip initial value
    ax3.set_title('Vegetation Cover (TRIFFID → EBM/Soil)')
    ax3.set_ylabel('Fraction')
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 1.05)
    
    # LAI (TRIFFID → Physiology)
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(time, results['weekly']['LAI_total'][1:], 'g-')  # Skip initial value
    ax4.set_title('LAI (TRIFFID → Physiology)')
    ax4.set_ylabel('LAI (m²/m²)')
    ax4.grid(True, alpha=0.3)
    
    # Soil Conditions (Soil → RothC)
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.plot(time, results['weekly']['soil_moisture_mean'], 'b-')
    ax5.set_title('Soil Moisture (Soil → RothC)')
    ax5.set_ylabel('Soil Moisture (m³/m³)')
    ax5.grid(True, alpha=0.3)
    
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.plot(time, results['weekly']['soil_temp_mean'] - 273.15, 'r-')  # Convert K to °C
    ax6.set_title('Soil Temperature (Soil → RothC)')
    ax6.set_ylabel('Temperature (°C)')
    ax6.grid(True, alpha=0.3)
    
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.set_xlabel('Week')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(os.path.join(save_path, 'model_coupling.png'), dpi=300)
        
    return fig

def create_model_summary(results):
    """Create a summary of model state ranges and coupling variables"""
    summary = {
        'Vegetation': {
            'LAI_range': (min(results['weekly']['LAI_total']), max(results['weekly']['LAI_total'])),
            'nu_total_range': (min(results['weekly']['nu_total']), max(results['weekly']['nu_total'])),
            'Leaf_biomass_range': {
                'tree': (min(results['weekly']['Lb_tree']), max(results['weekly']['Lb_tree'])),
                'grass': (min(results['weekly']['Lb_grass']), max(results['weekly']['Lb_grass']))
            }
        },
        'Soil Carbon': {
            'total_range': (min(results['weekly']['soil_C_total']), max(results['weekly']['soil_C_total'])),
            'DPM_range': (min(results['weekly']['C_dpm']), max(results['weekly']['C_dpm'])),
            'RPM_range': (min(results['weekly']['C_rpm']), max(results['weekly']['C_rpm'])),
            'BIO_range': (min(results['weekly']['C_bio']), max(results['weekly']['C_bio'])),
            'HUM_range': (min(results['weekly']['C_hum']), max(results['weekly']['C_hum']))
        },
        'Soil Conditions': {
            'moisture_range': (min(results['weekly']['soil_moisture_mean']), 
                               max(results['weekly']['soil_moisture_mean'])),
            'temp_range_C': (min(results['weekly']['soil_temp_mean']) - 273.15, 
                             max(results['weekly']['soil_temp_mean']) - 273.15)
        },
        'Photosynthesis': {
            'NPP_range': (min(results['weekly']['NPP']), max(results['weekly']['NPP'])),
            'Pi_G_range': (min(results['weekly']['Pi_G']), max(results['weekly']['Pi_G'])),
            'R_p_range': (min(results['weekly']['R_p']), max(results['weekly']['R_p']))
        }
    }
    
    return summary

def print_model_summary(summary):
    """Print summary of model variables to check for reasonable values"""
    print("\n=== MODEL VARIABLE SUMMARY ===")
    
    print("\nVegetation:")
    print(f"  LAI range: {summary['Vegetation']['LAI_range'][0]:.2f} - {summary['Vegetation']['LAI_range'][1]:.2f} m²/m²")
    print(f"  Vegetation cover range: {summary['Vegetation']['nu_total_range'][0]:.2f} - {summary['Vegetation']['nu_total_range'][1]:.2f}")
    print(f"  Leaf biomass range (tree): {summary['Vegetation']['Leaf_biomass_range']['tree'][0]:.2f} - {summary['Vegetation']['Leaf_biomass_range']['tree'][1]:.2f} kg C/m²")
    print(f"  Leaf biomass range (grass): {summary['Vegetation']['Leaf_biomass_range']['grass'][0]:.2f} - {summary['Vegetation']['Leaf_biomass_range']['grass'][1]:.2f} kg C/m²")
    
    print("\nSoil Carbon:")
    print(f"  Total soil carbon range: {summary['Soil Carbon']['total_range'][0]:.2f} - {summary['Soil Carbon']['total_range'][1]:.2f} kg C/m²")
    print(f"  DPM range: {summary['Soil Carbon']['DPM_range'][0]:.2f} - {summary['Soil Carbon']['DPM_range'][1]:.2f} kg C/m²")
    print(f"  RPM range: {summary['Soil Carbon']['RPM_range'][0]:.2f} - {summary['Soil Carbon']['RPM_range'][1]:.2f} kg C/m²")
    print(f"  BIO range: {summary['Soil Carbon']['BIO_range'][0]:.2f} - {summary['Soil Carbon']['BIO_range'][1]:.2f} kg C/m²")
    print(f"  HUM range: {summary['Soil Carbon']['HUM_range'][0]:.2f} - {summary['Soil Carbon']['HUM_range'][1]:.2f} kg C/m²")
    
    print("\nSoil Conditions:")
    print(f"  Soil moisture range: {summary['Soil Conditions']['moisture_range'][0]:.3f} - {summary['Soil Conditions']['moisture_range'][1]:.3f} m³/m³")
    print(f"  Soil temperature range: {summary['Soil Conditions']['temp_range_C'][0]:.2f} - {summary['Soil Conditions']['temp_range_C'][1]:.2f} °C")
    
    print("\nPhotosynthesis:")
    print(f"  NPP range: {summary['Photosynthesis']['NPP_range'][0]:.4f} - {summary['Photosynthesis']['NPP_range'][1]:.4f} kg C/m²/week")
    print(f"  Pi_G range: {summary['Photosynthesis']['Pi_G_range'][0]:.2f} - {summary['Photosynthesis']['Pi_G_range'][1]:.2f} μmol CO₂/m²/s")
    print(f"  R_p range: {summary['Photosynthesis']['R_p_range'][0]:.2f} - {summary['Photosynthesis']['R_p_range'][1]:.2f} μmol CO₂/m²/s")
    
    # Check for any unreasonable values
    warnings = []
    
    if summary['Vegetation']['LAI_range'][1] > 12:
        warnings.append("WARNING: Maximum LAI > 12, which is unusually high")
    if summary['Vegetation']['LAI_range'][0] < 0:
        warnings.append("WARNING: Negative LAI values detected")
    
    if summary['Soil Conditions']['moisture_range'][1] > 0.6:
        warnings.append("WARNING: Soil moisture > 0.6 m³/m³, which is unusually high")
    if summary['Soil Conditions']['moisture_range'][0] < 0.05:
        warnings.append("WARNING: Very low soil moisture detected (< 0.05 m³/m³)")
    
    if summary['Photosynthesis']['NPP_range'][1] <= 0:
        warnings.append("WARNING: No positive NPP detected")
    
    if warnings:
        print("\n=== WARNINGS ===")
        for warning in warnings:
            print(warning)
    else:
        print("\nAll values appear to be within reasonable ranges.")

def visualize_jules_coupling(results, save_path='plots'):
    """Main function to create all plots for JULES coupling visualization"""
    # Create output directory
    output_dir = create_output_directory(save_path)
    
    # Generate plots
    plot_vegetation_dynamics(results, output_dir)
    plot_soil_carbon_dynamics(results, output_dir)
    plot_soil_hydrothermal(results, output_dir)
    plot_photosynthesis_components(results, output_dir)
    plot_model_coupling(results, output_dir)
    
    # Generate and print summary
    summary = create_model_summary(results)
    print_model_summary(summary)
    
    print(f"\nAll plots have been saved to: {output_dir}")
    
    return summary

# Example usage
if __name__ == "__main__":
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..'))
    
    from julesf.coupler.jules_master import run_jules_simulation
    
    # Run simulation
    results = run_jules_simulation()
    
    # Visualize results
    visualize_jules_coupling(results)