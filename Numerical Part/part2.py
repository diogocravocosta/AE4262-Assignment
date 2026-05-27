import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Setup Parameters & Baselines ---
gas = ct.Solution('gri30.yaml')
air = 'O2:1.0, N2:3.76'

# Baselines
phi_base = 0.8
P_base = 1e5    # 1 bar = 100,000 Pa
T_base = 300.0  # K

fuels = ['CH4', 'H2']
colors = {'CH4': '#1f77b4', 'H2': '#d62728'}
labels = {'CH4': 'Methane ($CH_4$)', 'H2': 'Hydrogen ($H_2$)'}

# =====================================================================
# a) Variation in Tad with Equivalence Ratio (phi)
# =====================================================================
phis = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]

plt.figure(figsize=(8, 5))
for fuel in fuels:
    T_ad_a = []
    for phi in phis:
        gas.TP = T_base, P_base
        gas.set_equivalence_ratio(phi, fuel, air)
        gas.equilibrate('HP') # Constant Enthalpy and Pressure
        T_ad_a.append(gas.T)
    plt.plot(phis, T_ad_a, '-o', color=colors[fuel], label=labels[fuel], linewidth=2)

plt.title('Adiabatic Flame Temp vs Equivalence Ratio\n(Baseline: $T_{in}$=300K, P=1 bar)')
plt.xlabel('Equivalence Ratio, $\phi$')
plt.ylabel('Adiabatic Flame Temperature, $T_{ad}$ (K)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig('Part2a_Tad_vs_Phi.pdf')
plt.show()

# =====================================================================
# b) Variation in Tad with Pressure
# =====================================================================
pressures_bar = [1, 10, 40]
pressures_pa = [p * 1e5 for p in pressures_bar]

plt.figure(figsize=(8, 5))
for fuel in fuels:
    T_ad_b = []
    for P in pressures_pa:
        gas.TP = T_base, P
        gas.set_equivalence_ratio(phi_base, fuel, air)
        gas.equilibrate('HP')
        T_ad_b.append(gas.T)
    plt.plot(pressures_bar, T_ad_b, '-s', color=colors[fuel], label=labels[fuel], linewidth=2, markersize=8)

plt.title('Adiabatic Flame Temp vs Pressure\n(Baseline: $\phi$=0.8, $T_{in}$=300K)')
plt.xlabel('Pressure (bar)')
plt.ylabel('Adiabatic Flame Temperature, $T_{ad}$ (K)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.xticks(pressures_bar) # Force x-axis ticks to show exactly 1, 10, 40
plt.legend()
plt.tight_layout()
plt.savefig('Part2b_Tad_vs_Pressure.pdf')
plt.show()

# =====================================================================
# c) Variation in Tad with Air Inlet Temperature
# =====================================================================
T_ins = [300, 600, 1100]

plt.figure(figsize=(8, 5))
for fuel in fuels:
    T_ad_c = []
    for T_in in T_ins:
        gas.TP = T_in, P_base
        gas.set_equivalence_ratio(phi_base, fuel, air)
        gas.equilibrate('HP')
        T_ad_c.append(gas.T)
    plt.plot(T_ins, T_ad_c, '-^', color=colors[fuel], label=labels[fuel], linewidth=2, markersize=8)

plt.title('Adiabatic Flame Temp vs Inlet Temperature\n(Baseline: $\phi$=0.8, P=1 bar)')
plt.xlabel('Inlet Temperature, $T_{in}$ (K)')
plt.ylabel('Adiabatic Flame Temperature, $T_{ad}$ (K)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.xticks(T_ins)
plt.legend()
plt.tight_layout()
plt.savefig('Part2c_Tad_vs_Tin.pdf')
plt.show()

# =====================================================================
# d) Variation of NO emissions with Equivalence Ratio
# =====================================================================
plt.figure(figsize=(8, 5))

# Find the index of NO in the species array just once
no_index = gas.species_index('NO')

for fuel in fuels:
    NO_ppm_d = []
    for phi in phis:
        gas.TP = T_base, P_base
        gas.set_equivalence_ratio(phi, fuel, air)
        gas.equilibrate('HP')
        
        # Safely extract the NO mole fraction using its index
        no_frac = gas.X[no_index]
        
        # Convert to ppm (parts per million)
        NO_ppm_d.append(no_frac * 1e6) 
        
    plt.plot(phis, NO_ppm_d, '-d', color=colors[fuel], label=labels[fuel], linewidth=2)

plt.title('Equilibrium NO Emissions vs Equivalence Ratio\n(Baseline: $T_{in}$=300K, P=1 bar)')
plt.xlabel('Equivalence Ratio, $\phi$')
plt.ylabel('NO Emissions (ppm)')
plt.yscale('log') # NO varies by orders of magnitude, log scale is best
plt.grid(True, which="both", linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig('Part2d_NO_vs_Phi.pdf')
plt.show()