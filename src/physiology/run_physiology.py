# run_physiology.py

from simulation import (
    run_photosynthesis,
    run_respiration,
    run_npp,
    generate_forcings
)
from parameters import SIM_SETTINGS
import visualization as viz
import matplotlib.pyplot as plt

if __name__ == '__main__':
    pft = 'needleleaf_tree'

    t, ps = run_photosynthesis(pft)

    _, resp = run_respiration(pft)

    _, npp = run_npp(pft)

    t_forc, T, I_par, ci, O2 = generate_forcings(SIM_SETTINGS)
    ca = SIM_SETTINGS['ca_ppm'] * 1e-6 * SIM_SETTINGS['P']

    viz.plot_forcings(t_forc, T, I_par, ci, ca)
    viz.plot_limiters(t, ps)                             # Wc, Wl, We
    viz.plot_photosynthesis(t, ps)                               # Wg, Rd, Ap, Ac
    viz.plot_respiration(t, resp)                        # R_dc, Pi_G, R_pm, R_pg, R_p
    viz.plot_npp(t, npp)                                 # Pi_net

    plt.show()
