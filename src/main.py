import pandas as pd
from constants import *
from thermal_script import *


if __name__ == "__main__":
    # Read the Excel file into a DataFrame
    
    if READ_FROM_EXCEL:
        # ne_melting_data = pd.read_excel(
        #     NE_DATA_FILEPATH, sheet_name=MELTING_SHEETNAME)
        # xe_melting_data = pd.read_excel(
        #     XE_DATA_FILEPATH, sheet_name=MELTING_SHEETNAME)
        # kr_melting_data = pd.read_excel(
        #     KR_DATA_FILEPATH, sheet_name=MELTING_SHEETNAME)
        
        # ne_sublimation_data = read_melting_data(
        #     NE_DATA_FILEPATH, SUBLIMATION_SHEETNAME)
        # xe_sublimation_data = read_melting_data(
        #     XE_DATA_FILEPATH, SUBLIMATION_SHEETNAME)
        # kr_sublimation_data = read_melting_data(
        #     KR_DATA_FILEPATH, SUBLIMATION_SHEETNAME)
        
        # ne_fusion_data = read_fusion_data(
        #     NE_DATA_FILEPATH, FUSION_SHEETNAME)
        # xe_fusion_data = read_fusion_data(
        #     XE_DATA_FILEPATH, FUSION_SHEETNAME)
        # kr_fusion_data = read_fusion_data(
        #     KR_DATA_FILEPATH, FUSION_SHEETNAME)

        # ne_heatsub_data = read_fusion_data(
        #     NE_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
        # xe_heatsub_data = read_fusion_data(
        #     XE_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
        # kr_heatsub_data = read_fusion_data(
        #     KR_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
        
        # ne_thermal_coeff_data = read_thermal_coeff_data(
        #     NE_DATA_FILEPATH, THERMAL_COEFFICIENT_SHEETNAME)
        # xe_thermal_coeff_data = read_thermal_coeff_data(
        #     XE_DATA_FILEPATH, THERMAL_COEFFICIENT_SHEETNAME)
        # kr_thermal_coeff_data = read_thermal_coeff_data(
        #     KR_DATA_FILEPATH, THERMAL_COEFFICIENT_SHEETNAME)
        
        ne_heat_capacity_data = read_data(
            NE_DATA_FILEPATH, HEAT_CAPACITY_SHEETNAME, "Heat Capacity", 2)
        xe_heat_capacity_data = read_data(
            XE_DATA_FILEPATH, HEAT_CAPACITY_SHEETNAME, "Heat Capacity", 2)
        kr_heat_capacity_data = read_data(
            KR_DATA_FILEPATH, HEAT_CAPACITY_SHEETNAME, "Heat Capacity", 2)
    
    else:
        # Read the txt files into DataFrames
        # ne_melting_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\neon_melting_data.txt", sep='\t')
        # xe_melting_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\xenon_melting_data.txt", sep='\t')
        # kr_melting_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\krypton_melting_data.txt", sep='\t')
        
        # ne_sublimation_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\neon_sublimation_data.txt", sep='\t')
        # xe_sublimation_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\xenon_sublimation_data.txt", sep='\t')
        # kr_sublimation_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\krypton_sublimation_data.txt", sep='\t')
        
        # Thermal Coefficients
        # ne_thermal_coeff_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\neon_thermal_coeff_data.txt", sep='\t')
        # xe_thermal_coeff_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\xenon_thermal_coeff_data.txt", sep='\t')
        # kr_thermal_coeff_data = pd.read_csv(
        #     f"{TXT_DATA_FILEPATH}\\krypton_thermal_coeff_data.txt", sep='\t')
        
        # Heat Capacity
        ne_heat_capacity_data = pd.read_csv(
            f"{TXT_DATA_FILEPATH}\\neon_heat_capacity_data.txt", sep='\t')
        xe_heat_capacity_data = pd.read_csv(
            f"{TXT_DATA_FILEPATH}\\xenon_heat_capacity_data.txt", sep='\t')
        kr_heat_capacity_data = pd.read_csv(
            f"{TXT_DATA_FILEPATH}\\krypton_heat_capacity_data.txt", sep='\t')
        
    # Save to txt files 
    
    # Melting Data
    # ne_melting_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\neon_melting_data.txt", sep='\t', index=False)
    # xe_melting_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\xenon_melting_data.txt", sep='\t', index=False)
    # kr_melting_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\krypton_melting_data.txt", sep='\t', index=False)
    
    # Sublimation Data
    # ne_sublimation_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\neon_sublimation_data.txt", sep='\t', index=False)
    # xe_sublimation_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\xenon_sublimation_data.txt", sep='\t', index=False)
    # kr_sublimation_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\krypton_sublimation_data.txt", sep='\t', index=False)
    
    # Fusion Data
    # ne_fusion_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\neon_fusion_data.txt", sep='\t', index=False)
    # xe_fusion_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\xenon_fusion_data.txt", sep='\t', index=False)
    # kr_fusion_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\krypton_fusion_data.txt", sep='\t', index=False)
    
    # Heat of Sublimation Data
    # ne_heatsub_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\neon_heat_of_sublimation_data.txt", sep='\t', index=False)
    # xe_heatsub_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\xenon_heat_of_sublimation_data.txt", sep='\t', index=False)
    # kr_heatsub_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\krypton_heat_of_sublimation_data.txt", sep='\t', index=False)
    
    # Thermal Coefficient Data
    # ne_thermal_coeff_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\neon_thermal_coeff_data.txt", sep='\t', index=False)
    # xe_thermal_coeff_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\xenon_thermal_coeff_data.txt", sep='\t', index=False)
    # kr_thermal_coeff_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\krypton_thermal_coeff_data.txt", sep='\t', index=False)
    
    # Heat Capacity Data
    ne_heat_capacity_data.to_csv(
        f"{TXT_DATA_FILEPATH}\\neon_heat_capacity_data.txt", sep='\t', index=False)
    xe_heat_capacity_data.to_csv(
        f"{TXT_DATA_FILEPATH}\\xenon_heat_capacity_data.txt", sep='\t', index=False)
    kr_heat_capacity_data.to_csv(
        f"{TXT_DATA_FILEPATH}\\krypton_heat_capacity_data.txt", sep='\t', index=False)

    # Fit constants to the data per gas
    # gas_coefficients = fit_melting_pressure_single_gas(
    #     ne_melting_data, 'neon', gas_params)
    # print("Fitted Constants for neon:")
    # print(gas_coefficients)
    
    if DISPLAY_PLOTS:
    #     plot_melting_gas_data(ne_melting_data, 'neon')
    #     plot_melting_pressure_deviation(kr_melting_data, 'krypton')
    
        # Enthalpty of Melting (Fusion)
        # plot_gas_data(xe_fusion_data, 'xenon',
        #               'Change in Enthalpy', 'Heat_of_Fusion', '\Delta H', 'kJ/mol')
        # plot_gas_data(ne_fusion_data, 'neon',
        #               'Change in Enthalpy', 'Heat_of_Fusion', '\Delta H', 'kJ/mol')
        # plot_gas_data(kr_fusion_data, 'krypton',
        #               'Change in Enthalpy', 'Heat_of_Fusion', '\Delta H', 'kJ/mol')
    
        # Enthalpty of Sublimation (Sublimation) 
        # plot_gas_data(xe_heatsub_data, 'xenon',
        #             'Change in Enthalpy', 'Heat_of_Sublimation', '\Delta H', 'kJ/mol')
        # plot_gas_data(ne_heatsub_data, 'neon',
        #             'Change in Enthalpy', 'Heat_of_Sublimation', '\Delta H', 'kJ/mol')
        # plot_gas_data(kr_heatsub_data, 'krypton',
        #             'Change in Enthalpy', 'Heat_of_Sublimation', '\Delta H', 'kJ/mol')
    
        #Enthalpy of Sublimation 
        # plot_gas_data(xe_heatsub_data, 'xenon',
        #             'Change in Enthalpy', 'Heat_of_Sublimation', '\Delta H', 'kJ/mol')
        # plot_gas_data(ne_sublimation_data, 'neon',
        #             'Pressure', 'Sublimation Pressure', 'P', 'MPa')
        # Sublimation Pressure
        # plot_gas_data(xe_sublimation_data, 'xenon',
        #             'Pressure', 'Sublimation Pressure', 'P', 'MPa')
        # plot_gas_data(kr_sublimation_data, 'krypton',
        #             'Pressure', 'Sublimation Pressure', 'P', 'MPa')
        # plot_gas_data(kr_heatsub_data, 'krypton',
        #             'Change in Enthalpy', 'Heat_of_Sublimation', '\Delta H', 'kJ/mol')
        
        # Thermal Expansion Coefficient
        # plot_gas_data(ne_thermal_coeff_data, 'neon','Thermal Expansion Coefficient',
        #               'Thermal Expansion Coefficient', 'alpha', '1/K')
        # plot_gas_data(xe_thermal_coeff_data, 'xenon','Thermal Expansion Coefficient',
        #               'Thermal Expansion Coefficient', 'alpha ', '1/K')
        # plot_gas_data(kr_thermal_coeff_data, 'krypton','Thermal Expansion Coefficient',
        #               'Thermal Expansion Coefficient', 'alpha', '1/K')
        
        # Heat Capacity
        plot_gas_data(ne_heat_capacity_data, 'neon', 'Heat Capacity',
                      'Heat Capacity', 'C_p', 'J/(mol*K)')
        plot_gas_data(xe_heat_capacity_data, 'xenon', 'Heat Capacity',
                      'Heat Capacity', 'C_p', 'J/(mol*K)')
        plot_gas_data(kr_heat_capacity_data, 'krypton', 'Heat Capacity',
                      'Heat Capacity', 'C_p', 'J/(mol*K)')
        