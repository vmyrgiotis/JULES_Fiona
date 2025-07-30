"""
TRIFFID-RothC Coupled Model - Minimal Runner
"""
import numpy as np
import matplotlib.pyplot as plt

from julesf.coupler.triffid_rothc_simulation import solve_triffid_rothc_euler
from julesf.rothc.parameters import POOLS, C0_default

def run_triffid_rothc_model():
    # Settings
    days = 365 * 2
    dt = 7.0
    t_span = (0, days)
    
    # Initial conditions
    triffid_init = [6.0, 3.0, 0.8, 0.2]
    rothc_init = [C0_default[p] for p in POOLS]
    
    print(f"=== TRIFFID-RothC Coupled ({days/365:.1f} years, dt={dt} days) ===")
    
    t, triffid_results, rothc_results, coupling_vars = solve_triffid_rothc_euler(
        t_span=t_span,
        triffid_init=triffid_init,
        rothc_init=rothc_init,
        dt=dt
    )
    
    # Debug: Check array dimensions
    print(f"Debug: t.shape = {t.shape}")
    print(f"Debug: triffid_results.shape = {triffid_results.shape}")
    print(f"Debug: rothc_results.shape = {rothc_results.shape}")
    
    plot_results(t, triffid_results, rothc_results, coupling_vars)
    print_summary(t, triffid_results, rothc_results, coupling_vars)
    
    return t, triffid_results, rothc_results, coupling_vars

def plot_results(t, triffid_results, rothc_results, coupling_vars):
    t_years = t / 365
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('TRIFFID-RothC Coupled Results', fontsize=14, y=0.98)
    
    # Vegetation LAI
    ax1 = axes[0, 0]
    ax1.plot(t_years, triffid_results[0, :], 'g-', label='Tree LAI')
    ax1.plot(t_years, triffid_results[1, :], 'orange', label='Grass LAI')
    ax1.set_ylabel('LAI')
    ax1.set_title('Vegetation LAI')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Soil carbon
    ax2 = axes[0, 1]
    total_soil_c = rothc_results.sum(axis=0)
    ax2.plot(t_years, total_soil_c, 'brown', linewidth=2)
    ax2.set_ylabel('Soil C (kg/m²)')
    ax2.set_title('Total Soil Carbon')
    ax2.grid(True, alpha=0.3)
    
    # Coupling variables
    ax3 = axes[1, 0]
    ax3.plot(t_years, coupling_vars['litterfall'] * 365, 'brown')
    ax3.set_xlabel('Time (years)')
    ax3.set_ylabel('Litterfall (kg C/m²/yr)')
    ax3.set_title('Coupling: Litterfall')
    ax3.grid(True, alpha=0.3)
    
    # Tree vs Grass Competition (CHANGED!)
    ax4 = axes[1, 1]
    ax4.plot(t_years, triffid_results[2, :], 'darkgreen', linewidth=2.5, label='Tree cover')
    ax4.plot(t_years, triffid_results[3, :], 'orange', linewidth=2.5, label='Grass cover')
    ax4.fill_between(t_years, 0, triffid_results[2, :], alpha=0.3, color='darkgreen')
    ax4.fill_between(t_years, triffid_results[2, :], 
                     triffid_results[2, :] + triffid_results[3, :], 
                     alpha=0.3, color='orange')
    ax4.set_xlabel('Time (years)')
    ax4.set_ylabel('Fractional Cover')
    ax4.set_title('Tree vs Grass Competition')
    ax4.set_ylim(0, 1)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

def print_summary(t, triffid_results, rothc_results, coupling_vars):
    print(f"\nFinal LAI: Tree={triffid_results[0,-1]:.2f}, Grass={triffid_results[1,-1]:.2f}")
    print(f"Final soil C: {rothc_results.sum(axis=0)[-1]:.2f} kg C/m²")
    print(f"Average litterfall: {coupling_vars['litterfall'].mean()*365:.3f} kg C/m²/yr")
    
    soil_c_change = rothc_results.sum(axis=0)[-1] - rothc_results.sum(axis=0)[0]
    years = (t[-1] - t[0]) / 365
    print(f"Soil C change: {soil_c_change:+.2f} kg C/m² ({soil_c_change/years:+.3f} kg C/m²/yr)")

def main():
    return run_triffid_rothc_model()

if __name__ == "__main__":
    results = main()