import pandas as pd
from constants import *
from thermal_script import *


if __name__ == "__main__":
    # Read the Excel file into a DataFrame
    ne_melting_data = read_melting_data(
        NE_DATA_FILEPATH, MELTING_SHEETNAME)
    xe_melting_data = read_melting_data(
        XE_DATA_FILEPATH, MELTING_SHEETNAME)
    kr_melting_data = read_melting_data(
        KR_DATA_FILEPATH, MELTING_SHEETNAME)
    
    # Fit constants to the data per gas
    gas_coefficients = fit_melting_pressure_single_gas(
        ne_melting_data, 'neon', gas_params)
    print("Fitted Constants for neon:")
    print(gas_coefficients)
    plot_melting_gas_data(ne_melting_data, 'neon')
    
    
    

