#visualization.py
import matplotlib.pyplot as plt

def plot_limiters(t, results):
    plt.figure()
    plt.plot(t, results['Wc'], label='Wc')
    plt.plot(t, results['Wl'], label='Wl')
    plt.plot(t, results['We'], label='We')
    plt.title('Rate Limiters')
    plt.xlabel('Time (h)')
    plt.ylabel('Rate')
    plt.legend()
    plt.grid(True)

def plot_forcings(t, T, I_par, ci, ca):
    plt.figure(figsize=(8,4))
    plt.plot(t, T,    label='Air T (°C)')
    plt.plot(t, I_par,label='PAR (µmol m⁻² s⁻¹)')
    plt.plot(t, ci,   label='ci (Pa)')
    #plt.plot(t, O2,   label='O₂ (Pa)')
    plt.axhline(ca, linestyle='--', label='ca (Pa)')
    plt.title('Forcing Time Series')
    plt.xlabel('Time (h)')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True)

def plot_photosynthesis(t, results):
    plt.figure()
    plt.plot(t, results['Wg'], label='Gross W')
    plt.plot(t, results['Rd'], label='Rd')
    plt.plot(t, results['Ap'], label='Ap')
    plt.plot(t, results['Ac'], label='Ac')
    plt.title('Gross, Respiration, Net')
    plt.xlabel('Time (h)')
    plt.ylabel('Rate')
    plt.legend()
    plt.grid(True)

def plot_respiration(t, resp):
    plt.figure()
    plt.plot(t, resp['R_dc'], label='R_dc')
    plt.plot(t, resp['Pi_G'], label='Pi_G')
    plt.plot(t, resp['R_pm'], label='R_pm')
    plt.plot(t, resp['R_pg'], label='R_pg')
    plt.plot(t, resp['R_p'],  label='R_p')
    plt.title('Canopy Respiration & Gross Prod')
    plt.legend(); plt.grid(True)

def plot_npp(t, npp):
    plt.figure()
    plt.plot(t, npp['Pi_net'], label='Pi_net', linewidth=2)
    plt.title('Net Primary Productivity')
    plt.xlabel('Time (h)'); plt.ylabel('Rate')
    plt.legend(); plt.grid(True)
