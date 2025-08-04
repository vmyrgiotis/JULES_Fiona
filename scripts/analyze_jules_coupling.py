"""
Script to run the JULES model and analyze coupling between components
"""
import os
import sys
import matplotlib.pyplot as plt

# Add the project root to the Python path
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from src.julesf.coupler.jules_master import run_jules_simulation
from src.julesf.visualization.coupling_plots import visualize_jules_coupling

def main():
    print("Running JULES model simulation...")
    
    # Configure simulation
    config = {
        'weeks': 52,  # Run for 52 weeks (1 year)
        # Optional: override default initial conditions
        # 'triffid_init': [6.0, 3.0, 0.8, 0.2],  # [Lb_tree, Lb_grass, nu_tree, nu_grass]
        # 'rothc_init': None,  # Will use defaults
        # 'soil_init': None,  # Will use defaults
    }
    
    # Run simulation
    results = run_jules_simulation(config)
    
    # Set plot style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Visualize and analyze coupling
    print("\nGenerating coupling analysis plots...")
    summary = visualize_jules_coupling(results, save_path='plots/coupling_analysis')
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()