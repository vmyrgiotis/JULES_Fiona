# Create: src/julesf/ebm/diagnostics.py

import numpy as np
import matplotlib.pyplot as plt
from julesf.ebm.simulation import run_ebm
from julesf.ebm.run_ebm import drivers
from julesf.ebm.parameters import sigma, epsilon_s, alpha
from julesf.ebm.equations import H_flux, E_flux, G_flux

def diagnose_ebm_physics(days=2, dt_out=0.5/24):

    t_span = (0, days)
    T0 = 283.0  
    
    # Run EBM
    t, Ts = run_ebm(t_span, T0, drivers, dt_out=dt_out)
    
    # Calculate all energy balance components
    results = {}
    
    for i, time in enumerate(t):
        T_surf = Ts[i]
        
        # Forcing terms
        Sdn_val = drivers["Sdn"](time)
        Ldn_val = drivers["Ldn"](time) 
        Tair_val = drivers["Tair"](time)
        Q1_val = drivers["Q1"](time)
        nu_val = drivers["nu"](time)
        Tsoil_val = drivers["Tsoil"](time)
        
        # Radiation components
        Sw_absorbed = (1 - alpha) * Sdn_val
        Lw_absorbed = epsilon_s * Ldn_val
        Lw_emitted = sigma * epsilon_s * T_surf**4
        
        # Turbulent fluxes
        H_val = H_flux(T_surf, Tair_val)
        E_val = E_flux(T_surf, Q1_val)
        G_val = G_flux(T_surf, Tsoil_val, nu_val)
        
        # Net radiation and energy balance
        Rn = Sw_absorbed + Lw_absorbed - Lw_emitted
        Net_energy = Rn - H_val - E_val - G_val
        
        # Store results
        if i == 0:
            for key in ['time', 'T_surf', 'T_air', 'T_soil', 
                       'Sdn', 'Sw_absorbed', 'Ldn', 'Lw_absorbed', 'Lw_emitted',
                       'H', 'E', 'G', 'Rn', 'Net_energy']:
                results[key] = []
        
        results['time'].append(time)
        results['T_surf'].append(T_surf)
        results['T_air'].append(Tair_val)
        results['T_soil'].append(Tsoil_val)
        results['Sdn'].append(Sdn_val)
        results['Sw_absorbed'].append(Sw_absorbed)
        results['Ldn'].append(Ldn_val)
        results['Lw_absorbed'].append(Lw_absorbed)
        results['Lw_emitted'].append(Lw_emitted)
        results['H'].append(H_val)
        results['E'].append(E_val)
        results['G'].append(G_val)
        results['Rn'].append(Rn)
        results['Net_energy'].append(Net_energy)
    
    # Convert to numpy arrays
    for key in results:
        results[key] = np.array(results[key])
    
    return results

def plot_ebm_diagnostics(results):
    """
    Create comprehensive diagnostic plots
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('EBM Physics Diagnostics', fontsize=16)
    
    t = results['time']
    
    # 1. Temperature evolution
    ax1 = axes[0, 0]
    ax1.plot(t, results['T_surf'], 'r-', label='T_surf', linewidth=2)
    ax1.plot(t, results['T_air'], 'b--', label='T_air')
    ax1.plot(t, results['T_soil'], 'g:', label='T_soil')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('Temperature Components')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Radiation balance
    ax2 = axes[0, 1]
    ax2.plot(t, results['Sdn'], 'orange', label='Sdn (incoming SW)')
    ax2.plot(t, results['Sw_absorbed'], 'red', label='SW absorbed')
    ax2.plot(t, results['Ldn'], 'purple', label='Ldn (incoming LW)')
    ax2.plot(t, results['Lw_absorbed'], 'blue', label='LW absorbed')
    ax2.plot(t, -results['Lw_emitted'], 'black', label='-LW emitted')
    ax2.set_ylabel('Radiation (W/m²)')
    ax2.set_title('Radiation Components')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Turbulent fluxes
    ax3 = axes[0, 2]
    ax3.plot(t, results['H'], 'red', label='H (sensible)')
    ax3.plot(t, results['E'], 'blue', label='E (latent)')
    ax3.plot(t, results['G'], 'green', label='G (ground)')
    ax3.set_ylabel('Heat Flux (W/m²)')
    ax3.set_title('Turbulent Fluxes')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Net energy balance
    ax4 = axes[1, 0]
    ax4.plot(t, results['Rn'], 'black', label='Rn (net radiation)', linewidth=2)
    ax4.plot(t, results['Net_energy'], 'red', label='Net energy balance')
    ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax4.set_ylabel('Energy (W/m²)')
    ax4.set_title('Energy Balance')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. Diurnal patterns
    ax5 = axes[1, 1]
    # Focus on first day for diurnal pattern
    day1_mask = t <= 1.0
    t_day1 = t[day1_mask]
    ax5.plot(t_day1*24, results['Sdn'][day1_mask], 'orange', label='Solar')
    ax5.plot(t_day1*24, results['H'][day1_mask], 'red', label='Sensible')
    ax5.plot(t_day1*24, results['E'][day1_mask], 'blue', label='Latent')
    ax5.set_xlabel('Hour of day')
    ax5.set_ylabel('Flux (W/m²)')
    ax5.set_title('Diurnal Cycle (Day 1)')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. Physics checks
    ax6 = axes[1, 2]
    # Check energy conservation
    energy_imbalance = np.abs(results['Net_energy'])
    ax6.plot(t, energy_imbalance, 'red', label='|Energy imbalance|')
    ax6.set_ylabel('|Energy imbalance| (W/m²)')
    ax6.set_xlabel('Time (days)')
    ax6.set_title('Energy Conservation Check')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    ax6.set_yscale('log')
    
    plt.tight_layout()
    plt.show()
    
    return fig

def check_physics_sanity(results):
    """
    Print physics sanity checks
    """
    print("=== EBM Physics Sanity Checks ===\n")
    
    # Temperature ranges
    T_surf_range = [results['T_surf'].min(), results['T_surf'].max()]
    T_air_range = [results['T_air'].min(), results['T_air'].max()]
    print(f"Surface temp range: {T_surf_range[0]:.1f} - {T_surf_range[1]:.1f} K")
    print(f"Air temp range: {T_air_range[0]:.1f} - {T_air_range[1]:.1f} K")
    
    # Radiation checks
    Sdn_max = results['Sdn'].max()
    Lw_emitted_range = [results['Lw_emitted'].min(), results['Lw_emitted'].max()]
    print(f"\nMax solar radiation: {Sdn_max:.1f} W/m²")
    print(f"LW emission range: {Lw_emitted_range[0]:.1f} - {Lw_emitted_range[1]:.1f} W/m²")
    
    # Flux magnitudes
    H_range = [results['H'].min(), results['H'].max()]
    E_range = [results['E'].min(), results['E'].max()]
    G_range = [results['G'].min(), results['G'].max()]
    print(f"\nSensible heat range: {H_range[0]:.1f} - {H_range[1]:.1f} W/m²")
    print(f"Latent heat range: {E_range[0]:.1f} - {E_range[1]:.1f} W/m²")
    print(f"Ground heat range: {G_range[0]:.1f} - {G_range[1]:.1f} W/m²")
    
    # Energy balance check
    max_imbalance = np.abs(results['Net_energy']).max()
    print(f"\nMax energy imbalance: {max_imbalance:.1e} W/m²")
    
    # Physics warnings
    print("\n=== Potential Issues ===")
    if Sdn_max > 1200:
        print("⚠️  Solar radiation seems high (>1200 W/m²)")
    if T_surf_range[1] - T_surf_range[0] > 50:
        print("⚠️  Large surface temperature swings (>50K)")
    if max_imbalance > 1.0:
        print("⚠️  Poor energy conservation (>1 W/m²)")
    if np.abs(H_range[1]) > 500 or np.abs(E_range[1]) > 500:
        print("⚠️  Very large turbulent fluxes (>500 W/m²)")
    
    print("\n" + "="*50)

def main():
    """Run EBM diagnostics"""
    print("Running EBM physics diagnostics...")
    
    results = diagnose_ebm_physics(days=2, dt_out=0.02)
    check_physics_sanity(results)
    plot_ebm_diagnostics(results)
    
    return results

if __name__ == "__main__":
    results = main()