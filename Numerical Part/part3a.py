import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# --- Setup ---
gas = ct.Solution('gri30.yaml')
phis = [0.5, 0.8, 1.0]
fuels = ['CH4', 'H2']

# Species indices for standard state lookups
idx_NO = gas.species_index('NO')
idx_N2 = gas.species_index('N2')
idx_O2 = gas.species_index('O2')

# Molecular weights (kg/kmol)
mw_NO = gas.molecular_weights[idx_NO]
mw_N2 = gas.molecular_weights[idx_N2]
mw_O2 = gas.molecular_weights[idx_O2]
mw_CH4 = gas.molecular_weights[gas.species_index('CH4')]
mw_H2 = gas.molecular_weights[gas.species_index('H2')]

data = {'CH4': [], 'H2': []}

for fuel in fuels:
    for phi in phis:
        # ========================================================
        # GET T_EQ FROM PART 2 (Baseline: 300K, 1 bar)
        # ========================================================
        gas.TP = 300, ct.one_atm
        gas.set_equivalence_ratio(phi, fuel, 'O2:1.0, N2:3.76')
        gas.equilibrate('HP')
        T_eq = gas.T 

        # ========================================================
        # GLOBAL REACTION MOLES (per 1 mole of fuel)
        # ========================================================
        if fuel == 'CH4':
            # CH4 + (2/phi)[O2 + 3.76 N2] -> CO2 + 2H2O + 2(1/phi - 1)O2 + (7.52/phi)N2
            n_O2_initial = 2.0 / phi - 2.0
            n_N2_initial = 7.52 / phi
            mass_total = (1.0 * mw_CH4) + (2.0 / phi) * mw_O2 + (7.52 / phi) * mw_N2
        else:
            # H2 + (0.5/phi)[O2 + 3.76 N2] -> H2O + 0.5(1/phi - 1)O2 + (1.88/phi)N2
            n_O2_initial = 0.5 / phi - 0.5
            n_N2_initial = 1.88 / phi
            mass_total = (1.0 * mw_H2) + (0.5 / phi) * mw_O2 + (1.88 / phi) * mw_N2

        # ========================================================
        # EQUILIBRIUM REACTION (N2 + O2 <=> 2NO)
        # ========================================================
        # Set gas to T_eq to extract standard state Gibbs free energies
        gas.TP = T_eq, ct.one_atm
        g_RT = gas.standard_gibbs_RT # G^0 / RT for all species
        
        # Delta G^0 / RT for the reaction: N2 + O2 <=> 2NO
        delta_G_RT = 2 * g_RT[idx_NO] - g_RT[idx_N2] - g_RT[idx_O2]
        
        # Kp = exp(-Delta G^0 / RT)
        Kp = np.exp(-delta_G_RT)
        
        # Solve quadratic: (4 - Kp)x^2 + Kp(n_N2 + n_O2)x - Kp(n_N2)(n_O2) = 0
        a = 4.0 - Kp
        b = Kp * (n_N2_initial + n_O2_initial)
        c = -Kp * n_N2_initial * n_O2_initial
        
        # Extent of reaction (x) - taking the positive root
        x = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
        
        # Moles of NO formed = 2x
        n_NO_formed = 2 * x
        
        # Convert to Mass Fraction (Y_NO)
        Y_NO = (n_NO_formed * mw_NO) / mass_total
        data[fuel].append(Y_NO)

        print(f"[{fuel} | Phi={phi}] T_eq={T_eq:.1f} K | Kp={Kp:.2e} | Y_NO={Y_NO:.2e}")

# ========================================================
# 3. PLOT RESULTS
# ========================================================
plt.figure(figsize=(8, 5))
plt.plot(phis, data['CH4'], '-o', label='Methane ($CH_4$)', linewidth=2, color='#1f77b4', markersize=8)
plt.plot(phis, data['H2'], '-s', label='Hydrogen ($H_2$)', linewidth=2, color='#d62728', markersize=8)

plt.title('NO Mass Fraction as a function of $\phi$', fontsize=14)
plt.xlabel('Equivalence Ratio, $\phi$', fontsize=12)
plt.ylabel('NO Mass Fraction ($Y_{NO}$)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.xticks(phis)
plt.tight_layout()
plt.savefig('Part3a_TwoStep_NO_MassFraction.pdf')
plt.show()