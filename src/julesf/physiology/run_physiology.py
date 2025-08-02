# run_physiology.py
import matplotlib.pyplot as plt

from julesf.physiology.parameters import SIM_SETTINGS
from julesf.physiology.simulation import run_photosynthesis, run_respiration, run_npp, generate_forcings
from julesf.physiology.visualization import (
    plot_forcings, plot_limiters, plot_photosynthesis, plot_respiration, plot_npp
)

if __name__ == '__main__':
    pft = 'needleleaf_tree'
    t_ps,   ps   = run_photosynthesis(pft)
    t_resp, resp = run_respiration(pft)
    t_npp,  npp  = run_npp(pft)
    t_f,    T, I_par, ci, O2 = generate_forcings(SIM_SETTINGS)
    ca = SIM_SETTINGS['ca_ppm'] * 1e-6 * SIM_SETTINGS['P']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Physiology Diagnostics: {pft}', fontsize=14)

    plot_forcings        (t_f,    T, I_par, ci, ca,    ax=axes[0,0], show=False)
    plot_limiters        (t_ps,   ps,             ax=axes[0,1], show=False)
    plot_photosynthesis  (t_ps,   ps,             ax=axes[0,2], show=False)
    plot_respiration     (t_resp, resp,           ax=axes[1,0], show=False)
    plot_npp             (t_npp,  npp,            ax=axes[1,1], show=False)

    axes[1,2].axis('off')  # blank panel
    plt.tight_layout(pad=2.0, rect=[0,0,1,0.95])
    plt.show()
