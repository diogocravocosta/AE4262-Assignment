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

data_ppm = {'CH4': [], 'H2': []}

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
        # STEP 1: GLOBAL REACTION MOLES (per 1 mole of fuel)
        # ========================================================
        if fuel == 'CH4':
            # CH4 + (2/phi)[O2 + 3.76 N2] -> 1 CO2 + 2 H2O + 2(1/phi - 1)O2 + (7.52/phi)N2
            n_O2_initial = 2.0 / phi - 2.0
            n_N2_initial = 7.52 / phi
            
            # Total moles = CO2 + H2O + O2 + N2
            n_total = 1.0 + 2.0 + n_O2_initial + n_N2_initial
        else:
            # H2 + (0.5/phi)[O2 + 3.76 N2] -> 1 H2O + 0.5(1/phi - 1)O2 + (1.88/phi)N2
            n_O2_initial = 0.5 / phi - 0.5
            n_N2_initial = 1.88 / phi
            
            # Total moles = H2O + O2 + N2
            n_total = 1.0 + n_O2_initial + n_N2_initial

        # ========================================================
        # STEP 2: EQUILIBRIUM REACTION (N2 + O2 <=> 2NO)
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
        
        # Convert to Mole Fraction (X_NO) and then ppm
        X_NO = n_NO_formed / n_total
        NO_ppm = X_NO * 1e6
        
        # Floor value to avoid log(0) errors at phi=1.0
        plot_ppm = max(NO_ppm, 1e-5)
        data_ppm[fuel].append(plot_ppm)

# ========================================================
# 3. PLOT RESULTS
# ========================================================
plt.figure(figsize=(8, 5))
plt.plot(phis, data_ppm['CH4'], '-o', label='Methane ($CH_4$)', linewidth=2, color='#1f77b4', markersize=8)
plt.plot(phis, data_ppm['H2'], '-s', label='Hydrogen ($H_2$)', linewidth=2, color='#d62728', markersize=8)

#plt.title('NO Emissions vs $\phi$', fontsize=14)
plt.xlabel('Equivalence Ratio, $\phi$', fontsize=12)
plt.ylabel('NO Emissions (ppm)')

# Set to Log Scale
plt.yscale('log')
plt.ylim([1e-2, 1e4])


# --- NEW: Lock the Y-axis bounds to match Part 2d ---
#plt.ylim(1, 10000) 

plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend(fontsize=12)
plt.xticks(phis)
plt.tight_layout()
plt.savefig('Part3a_TwoStep_NO_ppm_log.pdf')
plt.show()