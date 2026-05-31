import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the standard GRI-Mech 3.0 mechanism 
gas = ct.Solution('gri30.yaml')

# Define the species list and temperature range
species_list = ['H2', 'O2', 'N2', 'CO2', 'H2O', 'CH4']
temperatures = np.linspace(300, 2400, 20) # 300 K to 2400 K
pressure = ct.one_atm                     # 101325 Pascals

# Retrieve Cp values for each species
cp_data = {'Temperature (K)': temperatures}

for sp in species_list:
    cp_values = []
    for T in temperatures:
        # Set the gas state to pure species (100% mole fraction) at T and P
        gas.TPX = T, pressure, {sp: 1.0}
        
        # Extract the mass-specific heat capacity (J / kg K)
        cp_values.append(gas.cp_mass) 
        
    cp_data[sp] = cp_values

# 4. Generate the Table using Pandas
df_cp = pd.DataFrame(cp_data)
print("Mass Specific Heat Capacity (Cp) in J/(kg K)")
print("-" * 60)
print(df_cp.to_string(index=False, float_format="%.2f"))

# Save the table to a CSV file for Excel
df_cp.to_csv("cantera_cp_table.csv", index=False)

# Plot the data with a logarithmic Y-axis
plt.figure(figsize=(10, 6))
for sp in species_list:
    plt.plot(df_cp['Temperature (K)'], df_cp[sp], marker='o', markersize=4, label=sp)

#plt.title('Variation of Specific Heat ($C_p$) with Temperature Using using Cantera')
plt.xlabel('Temperature (K)')
plt.ylabel('$C_p$ [J / (kg K)] - Log Scale')

# --- THIS IS THE KEY LINE FOR THE LOG SCALE ---
plt.yscale('log') 

# Clean up the grid lines so they look nice on a log scale
plt.grid(True, which="both", linestyle='--', alpha=0.5) 

plt.legend()
plt.tight_layout()
plt.show()