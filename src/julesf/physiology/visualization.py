#visualization.py
import matplotlib.pyplot as plt

def plot_limiters(t, results, ax=None, show=True):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(t, results['Wc'], label='Wc')
    ax.plot(t, results['Wl'], label='Wl')
    ax.plot(t, results['We'], label='We')
    ax.set_title('Rate Limiters')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Rate')
    ax.legend()
    ax.grid(True)
    if show:
        plt.show()
    return ax

def plot_forcings(t, T, I_par, ci, ca, ax=None, show=True):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(t, T,    label='Air T (°C)')
    ax.plot(t, I_par,label='PAR (µmol m⁻² s⁻¹)')
    ax.plot(t, ci,   label='ci (Pa)')
    ax.axhline(ca, linestyle='--', label='ca (Pa)')
    ax.set_title('Forcing Time Series')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True)
    if show:
        plt.show()
    return ax

def plot_photosynthesis(t, results, ax=None, show=True):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(t, results['Wg'], label='Gross W')
    ax.plot(t, results['Rd'], label='Rd')
    ax.plot(t, results['Ap'], label='Ap')
    ax.plot(t, results['Ac'], label='Ac')
    ax.set_title('Gross, Respiration, Net')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Rate')
    ax.legend()
    ax.grid(True)
    if show:
        plt.show()
    return ax

def plot_respiration(t, resp, ax=None, show=True):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(t, resp['R_dc'], label='R_dc')
    ax.plot(t, resp['Pi_G'], label='Pi_G')
    ax.plot(t, resp['R_pm'], label='R_pm')
    ax.plot(t, resp['R_pg'], label='R_pg')
    ax.plot(t, resp['R_p'],  label='R_p')
    ax.set_title('Canopy Resp & Gross Prod')
    ax.legend()
    ax.grid(True)
    if show:
        plt.show()
    return ax

def plot_npp(t, npp, ax=None, show=True):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(t, npp['Pi_net'], label='Pi_net', linewidth=2)
    ax.set_title('Net Primary Productivity')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Rate')
    ax.legend()
    ax.grid(True)
    if show:
        plt.show()
    return ax
