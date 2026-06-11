import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D

mpl.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'lines.linewidth': 2,
    'figure.dpi': 150,
})

gas = ct.Solution('gri30.yaml')
temperatures = np.linspace(300, 2400, 50)
phis = [0.5, 0.8, 1.0, 1.2, 1.4]
fuels = ['CH4', 'H2']
air = 'O2:1.0, N2:3.76'
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

for fuel in fuels:
    fig, ax = plt.subplots(figsize=(8, 5.5))

    for i, phi in enumerate(phis):
        cp_reactants = []
        cp_products = []

        gas.set_equivalence_ratio(phi, fuel, air)
        X_reactants = gas.X

        for T in temperatures:
            gas.TPX = T, ct.one_atm, X_reactants
            cp_reactants.append(gas.cp_mass)

        if fuel == 'CH4':
            if phi <= 1.0:
                X_products = {'CO2': 1.0, 'H2O': 2.0, 'O2': 2.0*(1.0/phi - 1.0), 'N2': 7.52/phi}
            else:
                X_products = {'CO2': 1.0, 'H2O': 2.0, 'N2': 7.52, 'CH4': phi - 1.0}
        elif fuel == 'H2':
            if phi <= 1.0:
                X_products = {'H2O': 1.0, 'O2': 0.5*(1.0/phi - 1.0), 'N2': 1.88/phi}
            else:
                X_products = {'H2O': 1.0, 'N2': 1.88, 'H2': phi - 1.0}

        for T in temperatures:
            gas.TPX = T, ct.one_atm, X_products
            cp_products.append(gas.cp_mass)

        ax.plot(temperatures, cp_reactants, '--', color=colors[i], linewidth=2)
        ax.plot(temperatures, cp_products,  '-',  color=colors[i], linewidth=2)

    fuel_name = r'CH$_4$/air' if fuel == 'CH4' else r'H$_2$/air'
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel(r'$c_{p,\mathrm{mix}}$ [J kg$^{-1}$ K$^{-1}$]')
    ax.set_title(f'Mixture $c_{{p,\\mathrm{{mix}}}}$ — {fuel_name}')
    ax.set_xlim(300, 2400)
    ax.grid(True, alpha=0.3)

    phi_handles   = [Line2D([0],[0], color=c, lw=2, label=rf'$\phi={p}$')
                     for p, c in zip(phis, colors)]
    style_handles = [Line2D([0],[0], color='k', lw=2, ls='--', label='Reactants'),
                     Line2D([0],[0], color='k', lw=2, ls='-',  label='Products')]
    ax.legend(handles=phi_handles + style_handles, ncol=2, framealpha=0.9, fontsize=9)

    fig.tight_layout()
    plt.savefig(f'Part1b_{fuel}_Cp.pdf')
    plt.show()