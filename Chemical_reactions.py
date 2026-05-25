import cantera as ct


def chemical_balance(fuel_name: str, oxidizer: str, n2_o2_ratio: float = 3.72):
    """
    Compute the stoichiometric balance for the combustion of a fuel.

    Parameters
    ----------
    fuel_name : str
        Cantera species name of the fuel (e.g. 'CH4', 'H2', 'C2H6').
    oxidizer : str
        'ox'  -> reaction with pure oxygen (O2)
        'air' -> reaction with air (O2 + n2_o2_ratio * N2)
    n2_o2_ratio : float, optional
        Molar N2/O2 ratio for air. Default is 3.72.

    Returns
    -------
    dict with stoichiometric coefficients, enthalpies, LHV, and mass ratios.
    """
    N2_O2_RATIO = n2_o2_ratio

    # --- Fuel composition from Cantera ---
    gas = ct.Solution('gri30.yaml')
    fuel_sp = gas.species(gas.species_index(fuel_name))
    atoms = fuel_sp.composition
    n_C = atoms.get('C', 0)
    n_H = atoms.get('H', 0)

    # --- Stoichiometric balance: CxHy + a O2 -> x CO2 + (y/2) H2O ---
    n_CO2 = n_C
    n_H2O = n_H / 2
    a_O2 = n_CO2 + n_H2O / 2
    

    if oxidizer.lower() == 'air':
        a_N2 = N2_O2_RATIO * a_O2
    elif oxidizer.lower() == 'ox':
        a_N2 = 0.0
    else:
        raise ValueError(f"oxidizer must be 'ox' or 'air', got '{oxidizer}'")

    # --- Enthalpies of formation from Cantera (J/mol, gas phase) ---
    def hf_cantera(species_name):
        sp = gas.species(gas.species_index(species_name))
        return sp.thermo.h(298.15) / 1000  # J/kmol -> J/mol

    Hf_fuel = hf_cantera(fuel_name)
    Hf_CO2 = hf_cantera('CO2')
    Hf_H2O = hf_cantera('H2O')
    Hf_O2 = hf_cantera('O2')
    Hf_N2 = hf_cantera('N2')

    # --- Energy balance (per mole of fuel) ---
    H_reactants = Hf_fuel + a_O2 * Hf_O2 + a_N2 * Hf_N2
    H_products = n_CO2 * Hf_CO2 + n_H2O * Hf_H2O + a_N2 * Hf_N2
    delta_H = H_products - H_reactants

    # --- Molar masses (g/mol) from Cantera ---
    def mw(species_name):
        return gas.molecular_weights[gas.species_index(species_name)]

    MW_fuel = mw(fuel_name)
    MW_O2   = mw('O2')
    MW_N2   = mw('N2')
    MW_CO2  = mw('CO2')
    MW_H2O  = mw('H2O')

    m_fuel = MW_fuel
    m_O2 = MW_O2 * a_O2
    m_N2 = MW_N2 * a_N2
    m_CO2 = MW_CO2 * n_CO2
    m_H2O = MW_H2O * n_H2O

    m_reactants = m_fuel + m_O2 + m_N2
    m_products = m_CO2 + m_H2O + m_N2

    # --- LHV ---
    LHV_kJ_per_g = -delta_H / 1000 / m_fuel
    LHV_J_per_kg = LHV_kJ_per_g * 1e6

    # --- Fuel-oxidizer ratios (by mass, stoichiometric) ---
    F_Ox_stoic = m_fuel / m_O2
    F_Air_stoic = m_fuel / (m_O2 + m_N2)
    Ox_F = 1 / F_Ox_stoic
    Air_F = 1 / F_Air_stoic if F_Air_stoic > 0 else float('inf')

    results = {
        'fuel_formula': fuel_name,
        'n_C': n_C,
        'n_H': n_H,
        'a_O2': a_O2,
        'a_N2': a_N2,
        'n_CO2': n_CO2,
        'n_H2O': n_H2O,
        'Hf_fuel_J_per_mol': Hf_fuel,
        'Hf_CO2_J_per_mol': Hf_CO2,
        'Hf_H2O_J_per_mol': Hf_H2O,
        'delta_H_J_per_mol': delta_H,
        'MW_fuel': MW_fuel,
        'm_fuel_g': m_fuel,
        'm_O2_g': m_O2,
        'm_N2_g': m_N2,
        'm_reactants_g': m_reactants,
        'm_products_g': m_products,
        'LHV_kJ_per_g': LHV_kJ_per_g,
        'LHV_J_per_kg': LHV_J_per_kg,
        'F_Ox_stoic': F_Ox_stoic,
        'F_Air_stoic': F_Air_stoic,
        'Ox_F': Ox_F,
        'Air_F': Air_F,
    }

    return results


if __name__ == '__main__':
    print("=== Methane + Air ===")
    res = chemical_balance('CH4', 'air')
    for k, v in res.items():
        print(f"  {k:25s}: {v}")

    print("\n=== Hydrogen + Air ===")
    res = chemical_balance('H2', 'air')
    for k, v in res.items():
        print(f"  {k:25s}: {v}")