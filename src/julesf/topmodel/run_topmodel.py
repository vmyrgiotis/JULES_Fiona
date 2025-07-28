from julesf.topmodel.simulation import simulate_topmodel
from julesf.topmodel.visualization import plot_runoff, plot_intermediates

if __name__ == "__main__":
    t, lambda_bar, W0, z_bar, Rb, lambda_c, f_sat, Rse = simulate_topmodel()

    # pick some cells to visualize
    corner_cells = [(0,0),(0,4),(4,0),(4,4)]

    plot_runoff(t, Rb, Rse, corner_cells)
    plot_intermediates(t, lambda_c, f_sat)
