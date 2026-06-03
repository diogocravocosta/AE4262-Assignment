import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# 1. Setup the gas and standard parameters
gas = ct.Solution('gri30.yaml')
temperatures = np.linspace(300, 2400, 50) # 50 points for smooth curves
phis = [0.5, 0.8, 1.0, 1.2, 1.4]
fuels = ['CH4', 'H2']
air = 'O2:1.0, N2:3.76' 
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

# 2. Loop through each fuel to create SEPARATE plots
for fuel in fuels:
    plt.figure(figsize=(10, 6))
    
    for i, phi in enumerate(phis):
        cp_reactants = []
        cp_products = []
        
        # --- A. REACTANTS ---
        gas.set_equivalence_ratio(phi, fuel, air)
        X_reactants = gas.X 
        
        for T in temperatures:
            gas.TPX = T, ct.one_atm, X_reactants
            cp_reactants.append(gas.cp_mass)
            
        # --- B. PRODUCTS (COMPLETE COMBUSTION) ---
        # Manually define the product moles based on global stoichiometry
        # Cantera will automatically normalize these into mole fractions
        
        if fuel == 'CH4':
            if phi <= 1.0:
                # Eq 10: Lean/Stoich Methane
                X_products = {'CO2': 1.0, 'H2O': 2.0, 'O2': 2.0*(1.0/phi - 1.0), 'N2': 7.52/phi}
            else:
                # Eq 11: Rich Methane
                X_products = {'CO2': 1.0, 'H2O': 2.0, 'N2': 7.52, 'CH4': phi - 1.0}
                
        elif fuel == 'H2':
            if phi <= 1.0:
                # Eq 15: Lean/Stoich Hydrogen
                X_products = {'H2O': 1.0, 'O2': 0.5*(1.0/phi - 1.0), 'N2': 1.88/phi}
            else:
                # Eq 16: Rich Hydrogen
                X_products = {'H2O': 1.0, 'N2': 1.88, 'H2': phi - 1.0}

        # Calculate Cp for the frozen complete combustion products
        for T in temperatures:
            gas.TPX = T, ct.one_atm, X_products
            cp_products.append(gas.cp_mass)
            
        # --- C. PLOTTING ---
        plt.plot(temperatures, cp_reactants, '--', color=colors[i], linewidth=2, 
                 label=f'Reactants ($\phi$={phi})')
        plt.plot(temperatures, cp_products, '-', color=colors[i], linewidth=2, 
                 label=f'Products ($\phi$={phi})')

    if fuel == 'CH4':
        fuel_name = 'Methane ($CH_4$)'
    else:
        fuel_name = 'Hydrogen ($H_2$)'
        
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('Mixture Specific Heat, $C_{p,mix}$ [J / (kg K)]', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc="lower right", title="Mixture State", fontsize=9)
    plt.tight_layout()
    # plt.savefig(f'Part1b_{fuel}_Cp.pdf') # Uncomment to save the new figures
    plt.show()