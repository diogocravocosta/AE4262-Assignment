import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# 1. Setup the gas and standard parameters
gas = ct.Solution('gri30.yaml')
temperatures = np.linspace(300, 2400, 50) # 50 points for smooth curves
phis = [0.5, 0.8, 1.0, 1.2, 1.4]
fuels = ['CH4', 'H2']

# Standard definition of air
air = 'O2:1.0, N2:3.76' 

# Define a list of distinct colors for the 5 equivalence ratios
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

# 2. Loop through each fuel to create SEPARATE plots
for fuel in fuels:
    # Create a new figure for this specific fuel
    plt.figure(figsize=(10, 6))
    
    for i, phi in enumerate(phis):
        cp_reactants = []
        cp_products = []
        
        # --- A. REACTANTS ---
        # Set the unburnt composition
        gas.set_equivalence_ratio(phi, fuel, air)
        X_reactants = gas.X # Save the reactant mole fractions
        
        for T in temperatures:
            # Set state using Temperature, Pressure, and fixed mole fractions (X)
            gas.TPX = T, ct.one_atm, X_reactants
            cp_reactants.append(gas.cp_mass)
            
        # --- B. PRODUCTS (FROZEN COMPOSITION) ---
        # "Ignite" the mixture from 300K to get the product composition
        # HP = Constant Enthalpy & Pressure (Adiabatic Flame)
        gas.TPX = 300, ct.one_atm, X_reactants
        gas.equilibrate('HP') 
        X_products = gas.X # Save the burnt product mole fractions (frozen)
        
        for T in temperatures:
            # Set state using swept Temp, fixed Pressure, and FROZEN product fractions
            gas.TPX = T, ct.one_atm, X_products
            cp_products.append(gas.cp_mass)
            
        # --- C. PLOTTING ---
        # Dashed line for Reactants, Solid line for Products
        plt.plot(temperatures, cp_reactants, '--', color=colors[i], linewidth=2, 
                 label=f'Reactants ($\phi$={phi})')
        plt.plot(temperatures, cp_products, '-', color=colors[i], linewidth=2, 
                 label=f'Products ($\phi$={phi})')

    # Formatting the plot for the current fuel
    if fuel == 'CH4':
        fuel_name = 'Methane ($CH_4$)'
    else:
        fuel_name = 'Hydrogen ($H_2$)'
        
    plt.title(f'Variation of Mixture $C_p$ for {fuel_name}-Air (Frozen Products)')
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('Mixture Specific Heat, $C_{p,mix}$ [J / (kg K)]', fontsize=12)
    
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Place the legend outside the plot so it doesn't cover the lines
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left", title="Mixture State")
    
    # Adjust layout to make room for the external legend and display
    plt.tight_layout()
    plt.show()