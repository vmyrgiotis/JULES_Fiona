import numpy as np
import matplotlib.pyplot as plt
from julesf.physiology.simulation import run_photosynthesis, run_respiration, run_npp, generate_forcings
from julesf.physiology.parameters import SIM_SETTINGS, CANOPY_PARAMS

def diagnose_physiology_physics(pft_key='C3_grass', days=2):
    """
    Comprehensive physiology physics diagnostic
    """
    t_ps, ps_results = run_photosynthesis(pft_key)
    t_resp, resp_results = run_respiration(pft_key)
    t_npp, npp_results = run_npp(pft_key)
    
    min_len = min(len(t_ps), len(t_resp), len(t_npp))
    t = t_ps[:min_len]
    
    settings = SIM_SETTINGS.copy()
    settings['days'] = days
    settings['dt_hours'] = 0.5
    
    t_forc, T_full, I_par_full, ci_full, O2_full = generate_forcings(settings)
    
    if len(T_full) >= min_len:
        T = T_full[:min_len]
        I_par = I_par_full[:min_len]
        ci = ci_full[:min_len]
        O2 = O2_full[:min_len]
    else:
        repeat_factor = int(np.ceil(min_len / len(T_full)))
        T = np.tile(T_full, repeat_factor)[:min_len]
        I_par = np.tile(I_par_full, repeat_factor)[:min_len]
        ci = np.tile(ci_full, repeat_factor)[:min_len]
        O2 = np.tile(O2_full, repeat_factor)[:min_len]
    
    for key in ps_results:
        ps_results[key] = ps_results[key][:min_len]
    for key in resp_results:
        resp_results[key] = resp_results[key][:min_len]
    for key in npp_results:
        npp_results[key] = npp_results[key][:min_len]
    
    print(f"Time grid lengths - PS: {len(t_ps)}, Resp: {len(t_resp)}, NPP: {len(t_npp)}")
    print(f"Using common length: {min_len}")
    
    results = {
        'time': t,
        'T': T,
        'I_par': I_par,
        'ci': ci,
        'O2': O2,
        **ps_results,
        **resp_results,
        **npp_results
    }
    
    return results

def plot_physiology_diagnostics(results):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Physiology Physics Diagnostics', fontsize=14)
    
    t = results['time']
    
    colors = {
        # Rate limiters
        'rubisco': '#E74C3C',      # Red for Rubisco
        'light': '#F39C12',        # Orange for light
        'export': '#27AE60',       # Green for export
        'gross': '#8E44AD',        # Purple for gross
        
        # Photosynthesis components
        'respiration': '#34495E',   # Dark gray for respiration
        'net_potential': '#2980B9', # Blue for net potential
        'canopy': '#16A085',       # Teal for canopy
        
        # Forcings
        'temperature': '#E67E22',   # Orange
        'par': '#F1C40F',          # Yellow
        'co2': '#95A5A6',          # Gray
        
        # NPP components
        'gross_prod': '#9B59B6',   # Light purple
        'maint_resp': '#C0392B',   # Dark red
        'growth_resp': '#D35400',  # Dark orange
        'net_prod': '#2E86AB',     # Ocean blue
    }
    
    # 1. Rate Limiters Comparison (Wc, Wl, We vs W_gross)
    ax1 = axes[0, 0]
    ax1.plot(t, results['Wc'], color=colors['rubisco'], linewidth=2, label='Wc (Rubisco)', alpha=0.8)
    ax1.plot(t, results['Wl'], color=colors['light'], linewidth=2, label='Wl (Light)', alpha=0.8)
    ax1.plot(t, results['We'], color=colors['export'], linewidth=2, label='We (Export)', alpha=0.8)
    ax1.plot(t, results['Wg'], color=colors['gross'], linewidth=3, label='W_gross (min)', linestyle='--')
    ax1.set_ylabel('Rate (mol CO₂ m⁻² s⁻¹)', fontsize=10)
    ax1.set_title('Rate Limiters vs Gross Photosynthesis', fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.2)
    
    # 2. Photosynthesis Balance Check (W_gross = Ap + Rd)
    ax2 = axes[0, 1]
    ax2.plot(t, results['Wg'], color=colors['gross'], linewidth=2.5, label='W_gross')
    ax2.plot(t, results['Rd'], color=colors['respiration'], linewidth=2, label='Rd (dark resp)')
    ax2.plot(t, results['Ap'], color=colors['net_potential'], linewidth=2, label='Ap (net potential)')
    ax2.plot(t, results['Rd'] + results['Ap'], color='black', linewidth=1.5, 
             linestyle=':', label='Rd + Ap (check)', alpha=0.7)
    ax2.set_ylabel('Rate (mol CO₂ m⁻² s⁻¹)', fontsize=10)
    ax2.set_title('Photosynthesis Balance: W = Ap + Rd', fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.2)
    
    # 3. Leaf-to-Canopy Scaling (Ap vs Ac)
    ax3 = axes[0, 2]
    ax3.plot(t, results['Ap'], color=colors['net_potential'], linewidth=2.5, label='Ap (leaf-level)')
    ax3.plot(t, results['Ac'], color=colors['canopy'], linewidth=2.5, label='Ac (canopy)')
    
    # Calculate expected scaling factor
    k_ext = CANOPY_PARAMS['k_ext']
    LAI = CANOPY_PARAMS['LAI']
    scaling_factor = (1.0 - np.exp(-k_ext * LAI)) / k_ext if k_ext > 0 else LAI
    expected_Ac = results['Ap'] * scaling_factor
    ax3.plot(t, expected_Ac, color='black', linewidth=1.5, linestyle=':', 
             label=f'Ap × {scaling_factor:.2f} (expected)', alpha=0.7)
    
    ax3.set_ylabel('Rate (mol CO₂ m⁻² s⁻¹)', fontsize=10)
    ax3.set_title('Leaf-to-Canopy Scaling', fontsize=11)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.2)
    
    # 4. Environmental Forcings
    ax4 = axes[1, 0]
    ax4_twin = ax4.twinx()
    
    # Temperature on left axis
    line1 = ax4.plot(t, results['T'] - 273.15, color=colors['temperature'], linewidth=2, label='Temperature (°C)')
    ax4.set_ylabel('Temperature (°C)', color=colors['temperature'], fontsize=10)
    ax4.tick_params(axis='y', labelcolor=colors['temperature'])
    
    # PAR on right axis
    line2 = ax4_twin.plot(t, results['I_par'], color=colors['par'], linewidth=2, label='PAR (μmol m⁻² s⁻¹)')
    ax4_twin.set_ylabel('PAR (μmol m⁻² s⁻¹)', color=colors['par'], fontsize=10)
    ax4_twin.tick_params(axis='y', labelcolor=colors['par'])
    
    ax4.set_title('Environmental Drivers', fontsize=11)
    ax4.grid(True, alpha=0.2)
    
    # Combined legend
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax4.legend(lines, labels, loc='upper left', fontsize=9)
    
    # 5. NPP Components
    ax5 = axes[1, 1]
    ax5.plot(t, results['Pi_G'], color=colors['gross_prod'], linewidth=2.5, label='Pi_G (gross prod)')
    ax5.plot(t, results['R_pm'], color=colors['maint_resp'], linewidth=2, label='R_pm (maintenance)')
    ax5.plot(t, results['R_pg'], color=colors['growth_resp'], linewidth=2, label='R_pg (growth)')
    ax5.plot(t, results['Pi_net'], color=colors['net_prod'], linewidth=2.5, label='Pi_net (NPP)')
    ax5.plot(t, results['Pi_G'] - results['R_p'], color='black', linewidth=1.5, 
             linestyle=':', label='Pi_G - R_p (check)', alpha=0.7)
    ax5.set_ylabel('Rate (mol CO₂ m⁻² s⁻¹)', fontsize=10)
    ax5.set_xlabel('Time (hours)', fontsize=10)
    ax5.set_title('NPP Balance: Pi_net = Pi_G - R_p', fontsize=11)
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.2)
    
    # 6. Physics Sanity Checks
    ax6 = axes[1, 2]
    
    # Check 1: W_gross should be minimum of limiters
    min_limiters = np.minimum.reduce([results['Wc'], results['Wl'], results['We']])
    gross_error = np.abs(results['Wg'] - min_limiters)
    
    # Check 2: Balance errors
    balance_error = np.abs(results['Wg'] - (results['Ap'] + results['Rd']))
    npp_error = np.abs(results['Pi_net'] - (results['Pi_G'] - results['R_p']))
    
    ax6.semilogy(t, gross_error, color=colors['rubisco'], linewidth=2, label='|W_gross - min(limiters)|')
    ax6.semilogy(t, balance_error, color=colors['net_potential'], linewidth=2, label='|W_gross - (Ap + Rd)|')
    ax6.semilogy(t, npp_error, color=colors['net_prod'], linewidth=2, label='|Pi_net - (Pi_G - R_p)|')
    
    ax6.set_ylabel('Absolute Error', fontsize=10)
    ax6.set_xlabel('Time (hours)', fontsize=10)
    ax6.set_title('Physics Conservation Checks', fontsize=11)
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.2)
    
    plt.tight_layout(pad=2.0)
    plt.subplots_adjust(top=0.93)
    plt.show()
    
    return fig

def check_physiology_sanity(results):
    """
    Print physics sanity checks for physiology
    """
    print("=== Physiology Physics Sanity Checks ===\n")
    
    # Rate limiter checks
    Wc_range = [results['Wc'].min(), results['Wc'].max()]
    Wl_range = [results['Wl'].min(), results['Wl'].max()]
    We_range = [results['We'].min(), results['We'].max()]
    Wg_range = [results['Wg'].min(), results['Wg'].max()]
    
    print("Rate Limiters (mol CO₂ m⁻² s⁻¹):")
    print(f"  Wc (Rubisco): {Wc_range[0]:.2e} - {Wc_range[1]:.2e}")
    print(f"  Wl (Light):   {Wl_range[0]:.2e} - {Wl_range[1]:.2e}")
    print(f"  We (Export):  {We_range[0]:.2e} - {We_range[1]:.2e}")
    print(f"  W_gross:      {Wg_range[0]:.2e} - {Wg_range[1]:.2e}")
    
    # Balance checks
    min_limiters = np.minimum.reduce([results['Wc'], results['Wl'], results['We']])
    gross_error = np.abs(results['Wg'] - min_limiters).max()
    balance_error = np.abs(results['Wg'] - (results['Ap'] + results['Rd'])).max()
    npp_error = np.abs(results['Pi_net'] - (results['Pi_G'] - results['R_p'])).max()
    
    print(f"\nBalance Checks:")
    print(f"  Max |W_gross - min(limiters)|: {gross_error:.2e}")
    print(f"  Max |W_gross - (Ap + Rd)|:    {balance_error:.2e}")
    print(f"  Max |Pi_net - (Pi_G - R_p)|:  {npp_error:.2e}")
    
    # Scaling checks
    k_ext = CANOPY_PARAMS['k_ext']
    LAI = CANOPY_PARAMS['LAI']
    expected_scaling = (1.0 - np.exp(-k_ext * LAI)) / k_ext if k_ext > 0 else LAI
    actual_scaling = np.mean(results['Ac'] / results['Ap']) if np.mean(results['Ap']) > 0 else 0
    
    print(f"\nLeaf-to-Canopy Scaling:")
    print(f"  Expected scaling factor: {expected_scaling:.3f}")
    print(f"  Actual scaling factor:   {actual_scaling:.3f}")
    print(f"  LAI: {LAI}, k_ext: {k_ext}")
    
    # Diurnal patterns
    T_range = [results['T'].min() - 273.15, results['T'].max() - 273.15]
    PAR_range = [results['I_par'].min(), results['I_par'].max()]
    NPP_range = [results['Pi_net'].min(), results['Pi_net'].max()]
    
    print(f"\nDiurnal Ranges:")
    print(f"  Temperature: {T_range[0]:.1f} - {T_range[1]:.1f} °C")
    print(f"  PAR: {PAR_range[0]:.0f} - {PAR_range[1]:.0f} μmol m⁻² s⁻¹")
    print(f"  NPP: {NPP_range[0]:.2e} - {NPP_range[1]:.2e} mol CO₂ m⁻² s⁻¹")
    
    # Physics warnings
    print("\n=== Potential Issues ===")
    if gross_error > 1e-10:
        print("⚠️  W_gross not properly taking minimum of limiters")
    if balance_error > 1e-10:
        print("⚠️  Photosynthesis balance error: W ≠ Ap + Rd")
    if npp_error > 1e-10:
        print("⚠️  NPP balance error: Pi_net ≠ Pi_G - R_p")
    if abs(actual_scaling - expected_scaling) > 0.01:
        print("⚠️  Leaf-to-canopy scaling mismatch")
    if NPP_range[1] < 0:
        print("⚠️  Negative NPP detected - respiration exceeds photosynthesis")
    if any(results['Wg'] < 0):
        print("⚠️  Negative gross photosynthesis detected")
    
    print("\n" + "="*50)

def main():
    """Run physiology diagnostics"""
    print("Running Physiology physics diagnostics...")
    
    # Test different PFTs
    pft = 'C3_grass'  # You can change this to test different plant types
    results = diagnose_physiology_physics(pft_key=pft, days=2)
    
    print(f"Testing PFT: {pft}")
    check_physiology_sanity(results)
    plot_physiology_diagnostics(results)
    
    return results

if __name__ == "__main__":
    results = main()