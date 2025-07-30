import matplotlib.pyplot as plt

def plot_runoff(t, Rb, Rse, cells):
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True, figsize=(10,6))
    for ix, iy in cells:
        ax1.plot(t, Rb[:,ix,iy], label=f'Rb {ix},{iy}')
        ax2.plot(t, Rse[:,ix,iy], label=f'Rse {ix},{iy}')
    ax1.set_ylabel('Rb (kg m⁻² s⁻¹)')
    ax2.set_ylabel('Rse (kg m⁻² s⁻¹)')
    ax2.set_xlabel('Days')
    ax1.legend(); ax2.legend()
    plt.tight_layout()
    plt.show()

def plot_intermediates(t, lambda_c, f_sat):
    lc = lambda_c.mean(axis=(1,2))
    fs = f_sat.mean(axis=(1,2))
    plt.figure(figsize=(10,4))
    plt.plot(t, lc, label='Mean λc')
    plt.plot(t, fs, label='Mean f_sat')
    plt.xlabel('Days'); plt.ylabel('Value')
    plt.legend(); plt.tight_layout(); plt.show()
