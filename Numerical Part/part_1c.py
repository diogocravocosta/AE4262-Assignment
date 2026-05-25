import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Load all species from GRI 3.0, but filter for ONLY the major ones
# This prevents dissociation (e.g., CO2 cannot break down into CO and O)
all_species = {S.name: S for S in ct.Species.list_from_file('gri30.yaml')}
complete_species_names = ['CH4', 'H2', 'O2', 'N2', 'CO2', 'H2O']
complete_species = [all_species[name] for name in complete_species_names]

# Create a custom non-dissociating gas
gas_complete = ct.Solution(thermo='ideal-gas', species=complete_species)

# Define parameters
# Using the phi range from Part 2a of your assignment
phis = np.linspace(0.5, 1.4, 10) 
T_inlet = 1100 # K
pressure = ct.one_atm
air = 'O2:1.0, N2:3.76'

T_ad_H2 = []
T_ad_CH4 = []

# Compute Adiabatic Flame Temperatures
for phi in phis:
    # --- Hydrogen ---
    gas_complete.TP = T_inlet, pressure
    gas_complete.set_equivalence_ratio(phi, 'H2', air)
    gas_complete.equilibrate('HP') # Constant Enthalpy & Pressure
    T_ad_H2.append(gas_complete.T)
    
    # --- Methane ---
    gas_complete.TP = T_inlet, pressure
    gas_complete.set_equivalence_ratio(phi, 'CH4', air)
    gas_complete.equilibrate('HP')
    T_ad_CH4.append(gas_complete.T)

print(phis)
print(T_ad_H2)
print(T_ad_CH4)
# Plotting the Comparison
plt.figure(figsize=(9, 6))

plt.plot(phis, T_ad_H2, '-o', color='#d62728', linewidth=2, label='Hydrogen ($H_2$)')
plt.plot(phis, T_ad_CH4, '-s', color='#1f77b4', linewidth=2, label='Methane ($CH_4$)')

plt.title('Complete Combustion $T_{ad}$ at $T_{in} = 1100$ K', fontsize=14)
plt.xlabel('Equivalence Ratio, $\phi$', fontsize=12)
plt.ylabel('Adiabatic Flame Temperature (K)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)

# Save and show
plt.savefig('T_ad_Complete_Combustion_1100K.pdf', format='pdf', bbox_inches='tight')
plt.show()