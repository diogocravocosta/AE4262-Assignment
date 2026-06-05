import cantera as ct
import matplotlib.pyplot as plt

phis = [0.6, 0.7, 0.8, 1.0, 1.2, 1.3, 1.4]
speeds = []

for phi in phis:
    gas = ct.Solution('gri30.yaml')
    gas.set_equivalence_ratio(phi, 'CH4', 'O2:1, N2:3.76')
    gas.TP = 300, ct.one_atm

    flame = ct.FreeFlame(gas, width=0.03)
    flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    flame.transport_model = 'mixture-averaged'
    flame.solve(loglevel=0, auto=True)

    sl = flame.velocity[0]
    speeds.append(sl)
    print(f"phi = {phi:.1f}  |  S_L = {sl:.6f} m/s")

# Print copy-paste friendly table
print("\n--- Copy-paste data ---")
print("phi\tS_L [m/s]")
for phi, sl in zip(phis, speeds):
    print(f"{phi:.1f}\t{sl:.6f}")

# Plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(phis, speeds, 'o-', color='#c0392b', markersize=8, linewidth=2)
ax.set_xlabel('Equivalence ratio φ', fontsize=13)
ax.set_ylabel('Laminar flame speed S$_L$ [m/s]', fontsize=13)
ax.set_title('CH$_4$/air – 1 atm, 300 K', fontsize=14)
ax.grid(True, alpha=0.3)