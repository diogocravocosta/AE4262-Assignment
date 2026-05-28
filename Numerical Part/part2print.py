import cantera as ct

# --- Setup Parameters ---
gas = ct.Solution('gri30.yaml')
air = 'O2:1.0, N2:3.76'

# Baselines
phi_base = 0.8
P_base = 1e5    # 1 bar = 100,000 Pa
T_base = 300.0  # K

phis = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
pressures_bar = [1, 10, 40]
T_ins = [300, 600, 1100]

# --- Initialize Data Storage ---
data = {
    'CH4': {'T_ad_a': [], 'T_ad_b': [], 'T_ad_c': [], 'NO_ppm': []},
    'H2':  {'T_ad_a': [], 'T_ad_b': [], 'T_ad_c': [], 'NO_ppm': []}
}

no_index = gas.species_index('NO')

# --- Computations ---
for fuel in ['CH4', 'H2']:
    
    # Tad vs Equivalence Ratio
    for phi in phis:
        gas.TP = T_base, P_base
        gas.set_equivalence_ratio(phi, fuel, air)
        gas.equilibrate('HP')
        data[fuel]['T_ad_a'].append(gas.T)
        
    # Tad vs Pressure
    for P_bar in pressures_bar:
        gas.TP = T_base, P_bar * 1e5
        gas.set_equivalence_ratio(phi_base, fuel, air)
        gas.equilibrate('HP')
        data[fuel]['T_ad_b'].append(gas.T)
        
    # Tad vs Inlet Temperature
    for T_in in T_ins:
        gas.TP = T_in, P_base
        gas.set_equivalence_ratio(phi_base, fuel, air)
        gas.equilibrate('HP')
        data[fuel]['T_ad_c'].append(gas.T)
        
    # NO Emissions vs Equivalence Ratio
    for phi in phis:
        gas.TP = T_base, P_base
        gas.set_equivalence_ratio(phi, fuel, air)
        gas.equilibrate('HP')
        no_frac = gas.X[no_index]
        data[fuel]['NO_ppm'].append(no_frac * 1e6)

# --- Print Results for Copy/Pasting ---
print("================ DATA DUMP ================\n")

print(f"--- PHIS: {phis} ---\n")

print("1. Part A (T_ad vs Phi)")
print(f"Methane T_ad : {[round(t, 2) for t in data['CH4']['T_ad_a']]}")
print(f"Hydrogen T_ad: {[round(t, 2) for t in data['H2']['T_ad_a']]}\n")

print("2. Part B (T_ad vs Pressure at 1, 10, 40 bar)")
print(f"Methane T_ad : {[round(t, 2) for t in data['CH4']['T_ad_b']]}")
print(f"Hydrogen T_ad: {[round(t, 2) for t in data['H2']['T_ad_b']]}\n")

print("3. Part C (T_ad vs Inlet Temp at 300, 600, 1100 K)")
print(f"Methane T_ad : {[round(t, 2) for t in data['CH4']['T_ad_c']]}")
print(f"Hydrogen T_ad: {[round(t, 2) for t in data['H2']['T_ad_c']]}\n")

print("4. Part D (NO ppm vs Phi)")
print(f"Methane NO ppm : {[round(n, 2) for n in data['CH4']['NO_ppm']]}")
print(f"Hydrogen NO ppm: {[round(n, 2) for n in data['H2']['NO_ppm']]}\n")

print("=====================================================")