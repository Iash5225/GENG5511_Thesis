import pandas as pd
from constants import *
from thermal_script import *


if __name__ == "__main__":
    # Read the Excel file into a DataFrame
    
    
    # ne_melting_data = read_melting_data(
    #     NE_DATA_FILEPATH, MELTING_SHEETNAME)
    # xe_melting_data = read_melting_data(
    #     XE_DATA_FILEPATH, MELTING_SHEETNAME)
    # kr_melting_data = read_melting_data(
    #     KR_DATA_FILEPATH, MELTING_SHEETNAME)
    
    # ne_fusion_data = read_fusion_data(
    #     NE_DATA_FILEPATH, FUSION_SHEETNAME)
    # xe_fusion_data = read_fusion_data(
    #     XE_DATA_FILEPATH, FUSION_SHEETNAME)
    # kr_fusion_data = read_fusion_data(
    #     KR_DATA_FILEPATH, FUSION_SHEETNAME)
    
    ne_heatsub_data = read_fusion_data(
        NE_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
    xe_heatsub_data = read_fusion_data(
        XE_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
    kr_heatsub_data = read_fusion_data(
        KR_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
    
    
    
    # Fit constants to the data per gas
    # gas_coefficients = fit_melting_pressure_single_gas(
    #     ne_melting_data, 'neon', gas_params)
    # print("Fitted Constants for neon:")
    # print(gas_coefficients)
    # plot_melting_gas_data(ne_melting_data, 'neon')
    # plot_melting_pressure_deviation(kr_melting_data, 'krypton')
    
    # Enthalpty of Melting (Fusion)
    # plot_gas_data(xe_fusion_data, 'xenon',
    #               'Change in Enthalpy', 'Heat_of_Fusion', '\Delta H', 'kJ/mol')
    # plot_gas_data(ne_fusion_data, 'neon',
    #               'Change in Enthalpy', 'Heat_of_Fusion', '\Delta H', 'kJ/mol')
    # plot_gas_data(kr_fusion_data, 'krypton',
    #               'Change in Enthalpy', 'Heat_of_Fusion', '\Delta H', 'kJ/mol')
    
    # Enthalpty of Sublimation (Sublimation) 
    # plot_gas_data(xe_heatsub_data, 'xenon',
    #               'Change in Enthalpy', 'Heat_of_Sublimation', '\Delta H', 'kJ/mol')
    # plot_gas_data(ne_heatsub_data, 'neon',
    #               'Change in Enthalpy', 'Heat_of_Sublimation', '\Delta H', 'kJ/mol')
    # plot_gas_data(kr_heatsub_data, 'krypton',
    #               'Change in Enthalpy', 'Heat_of_Sublimation', '\Delta H', 'kJ/mol')
    
    
    
    