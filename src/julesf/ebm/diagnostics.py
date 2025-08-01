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
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('EBM Physics Diagnostics', fontsize=14)
    
    t = results['time']
    
    colors = {
        'surf': '#E74C3C',      
        'air': '#3498DB',       
        'soil': '#27AE60',      
        'solar': '#F39C12',     
        'absorbed': '#E67E22',  
        'incoming_lw': '#9B59B6', 
        'absorbed_lw': '#2980B9', 
        'emitted': '#34495E',   
        'sensible': '#C0392B',  
        'latent': '#2E86AB',    
        'ground': '#16A085',    
        'net_rad': '#2C3E50',   
        'balance': '#E74C3C'    
    }
    
    # 1. Temperature evolution
    ax1 = axes[0, 0]
    ax1.plot(t, results['T_surf'], color=colors['surf'], linewidth=2.5, label='T_surf')
    ax1.plot(t, results['T_air'], color=colors['air'], linestyle='--', linewidth=2, label='T_air')
    ax1.plot(t, results['T_soil'], color=colors['soil'], linestyle=':', linewidth=2, label='T_soil')
    ax1.set_ylabel('Temperature (K)', fontsize=10)
    ax1.set_title('Temperature Components', fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.2)
    
    # 2. Radiation balance
    ax2 = axes[0, 1]
    ax2.plot(t, results['Sdn'], color=colors['solar'], linewidth=2, label='Sdn (incoming)')
    ax2.plot(t, results['Sw_absorbed'], color=colors['absorbed'], linewidth=2, label='SW absorbed')
    ax2.plot(t, results['Ldn'], color=colors['incoming_lw'], linewidth=1.5, label='Ldn (incoming LW)')
    ax2.plot(t, results['Lw_absorbed'], color=colors['absorbed_lw'], linewidth=1.5, label='LW absorbed')
    ax2.plot(t, -results['Lw_emitted'], color=colors['emitted'], linewidth=2, label='-LW emitted')
    ax2.set_ylabel('Radiation (W/m²)', fontsize=10)
    ax2.set_title('Radiation Components', fontsize=11)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.2)
    
    # 3. Turbulent fluxes
    ax3 = axes[1, 0]
    ax3.plot(t, results['H'], color=colors['sensible'], linewidth=2.5, label='H (sensible)')
    ax3.plot(t, results['E'], color=colors['latent'], linewidth=2.5, label='E (latent)')
    ax3.plot(t, results['G'], color=colors['ground'], linewidth=2, label='G (ground)')
    ax3.set_ylabel('Heat Flux (W/m²)', fontsize=10)
    ax3.set_title('Turbulent Fluxes', fontsize=11)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.2)
    
    # 4. Net energy balance
    ax4 = axes[1, 1]
    ax4.plot(t, results['Rn'], color=colors['net_rad'], linewidth=3, label='Rn (net radiation)')
    ax4.plot(t, results['Net_energy'], color=colors['balance'], linewidth=2, label='Net energy balance')
    ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.4)
    ax4.set_ylabel('Energy (W/m²)', fontsize=10)
    ax4.set_title('Energy Balance', fontsize=11)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.2)
    
    # # 5. Diurnal patterns
    # ax5 = axes[1, 1]
    # day1_mask = t <= 1.0
    # t_day1 = t[day1_mask]
    # ax5.plot(t_day1*24, results['Sdn'][day1_mask], color=colors['solar'], linewidth=2.5, label='Solar')
    # ax5.plot(t_day1*24, results['H'][day1_mask], color=colors['sensible'], linewidth=2, label='Sensible')
    # ax5.plot(t_day1*24, results['E'][day1_mask], color=colors['latent'], linewidth=2, label='Latent')
    # ax5.set_xlabel('Hour of day', fontsize=10)
    # ax5.set_ylabel('Flux (W/m²)', fontsize=10)
    # ax5.set_title('Diurnal Cycle (Day 1)', fontsize=11)
    # ax5.legend(fontsize=9)
    # ax5.grid(True, alpha=0.2)
    
    # # 6. Physics checks
    # ax6 = axes[1, 2]
    # energy_imbalance = np.abs(results['Net_energy'])
    # ax6.plot(t, energy_imbalance, color='#8E44AD', linewidth=2.5, label='|Energy imbalance|')
    # ax6.set_ylabel('|Energy imbalance| (W/m²)', fontsize=10)
    # ax6.set_xlabel('Time (days)', fontsize=10)
    # ax6.set_title('Energy Conservation Check', fontsize=11)
    # ax6.legend(fontsize=9)
    # ax6.grid(True, alpha=0.2)
    # ax6.set_yscale('log')
    
    # Make plots more compact
    plt.tight_layout(pad=2.0)
    plt.subplots_adjust(top=0.93)
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