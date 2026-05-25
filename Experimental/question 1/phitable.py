import csv
import cantera as ct
import CoolProp.CoolProp as CP
from Chemical_reactions import chemical_balance

def phi_table(fuel, phi_values, power, D_tube, P_amb=101325, T_amb=298.15, n2_o2_ratio: float = 3.72):
    """
    Power in Watts, P_amb in Pa, T_amb in K, D_tube in mm
    """

    # Hardcoded nominal conditions (from standard reference values)
    T_nom = 273.15  # K
    P_nom = 101325  # Pa
    
    fuel_name = fuel['fuel_formula']
    
    # --- GET MOLECULAR WEIGHT OF FUEL AND AIR ---
    MW_fuel = fuel['MW_fuel']  # Stored from chemical_balance (g/mol)
    
    # Get molecular weight of Air from CoolProp (converted from kg/mol to g/mol)
    MW_air = CP.PropsSI('M', '', 0, '', 0, 'Air') * 1000
    
    # --- CALCULATE FUEL DENSITIES WITH COOLPROP ---
    fuel_density_input = CP.PropsSI('D', 'T', T_amb, 'P', P_amb, fuel_name)
    fuel_density_nominal = CP.PropsSI('D', 'T', T_nom, 'P', P_nom, fuel_name)

    # --- CALCULATE AIR DENSITIES WITH COOLPROP ---
    air_density_input = CP.PropsSI('D', 'T', T_amb, 'P', P_amb, 'Air')
    air_density_nominal = CP.PropsSI('D', 'T', T_nom, 'P', P_nom, 'Air')
    
    # --- Power calculations ---
    react_fuel = power / (fuel['LHV_kJ_per_g'] * 1000)  # Convert LHV to J/kg
    react_air = react_fuel / fuel['F_Air_stoic']
    R = ct.gas_constant / 1000  # J/(mol*K)

    table = [] 

    for phi in phi_values:
        if phi <= 1:
            limit = 'fuel'
            FA = phi * fuel['F_Air_stoic']
            f_f = react_fuel
            air_f = react_fuel / FA
            
        else:
            limit = 'air'
            FA = phi * fuel['F_Air_stoic'] 
            air_f = react_air
            f_f = react_air * FA

        # Calculate mass flow rates.
        f_nlpm = f_f / fuel_density_nominal * 60    # Convert kg/s to NLPM    
        air_nlpm = air_f / air_density_nominal * 60  # Convert kg/s to NLPM
        X_f = 1 / (1 + fuel['a_O2'] * (n2_o2_ratio + 1) / phi)  # Mole fraction of fuel in the mixture
        X_air = 1 - X_f  # Mole fraction of air in the mixture
        MW_mix = X_f * MW_fuel + X_air * MW_air  # g/mol
        rho_mixure = (P_amb * MW_mix) / (R * T_amb) / 1000  # Convert g/mol to kg/mol for density calculation
        v_bulk = (f_f + air_f) / (1000*rho_mixure *3.14159 * (D_tube/1000/2)**2)  # Convert D_tube from mm to m
        reynolds_no = (rho_mixure * v_bulk * (D_tube/1000)) / (CP.PropsSI('V', 'T', T_amb, 'P', P_amb, 'Air'))  # Dynamic viscosity in kg/(m*s)


        # Mapping your loop variables exactly to your layout keys
        row = {
                'Φ': phi,
                'Limit': limit,
                'F/A': FA,
                'Fuel Flow [g/s]': f_f,
                'Air Flow [g/s]': air_f,
                'Fuel flow nlpm': f_nlpm,
                'Air flow nlpm': air_nlpm,
                'Molar fraction Fuel [g/mol]': X_f,   
                'Molar fraction air [g/mol]': X_air,     
                'Molar mass mixture [g/mol]': MW_mix,
                'rho mixture [kg/m3]': rho_mixure,
                'V bulk [m/s]': v_bulk,
                'Re': reynolds_no
            }
        
        table.append(row)

    # --- WRITE CSV FILE (Moved outside the loop) ---
    csv_filename = f"phi_table_{fuel_name}.csv"

    if table:
        headers = list(table[0].keys())
        with open(csv_filename, mode='w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            writer.writerows(table)

    return table


# --- DEFAULT RUN EXECUTION ---
if __name__ == '__main__':
    test_power = 1000.0  # 1 kW target system power

    # --- RUN RUNTIME CASE 1: METHANE (CH4) ---
    methane_phi_values = [0.8, 1.0, 1.2]

    print("Generating Methane data table...")
    methane_profile = chemical_balance('CH4', 'air')
    phi_table(
        fuel=methane_profile, 
        phi_values=methane_phi_values, 
        power=test_power, 
        D_tube=21.0,        # 21 mm tube
        P_amb=101325, 
        T_amb=300.0
    )
    print("-> Successfully saved 'phi_table_CH4.csv'\n")

    # --- RUN RUNTIME CASE 2: HYDROGEN (H2) ---
    import numpy as np
    hydrogen_phi_values = list(np.arange(0.5, 1.31, 0.1))

    print("Generating Hydrogen data table...")
    hydrogen_profile = chemical_balance('H2', 'air')
    phi_table(
        fuel=hydrogen_profile, 
        phi_values=hydrogen_phi_values, 
        power=test_power, 
        D_tube=4.0,         # 4 mm tube
        P_amb=101325, 
        T_amb=300.0
    )
    print("-> Successfully saved 'phi_table_H2.csv'")