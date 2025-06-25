import pandas as pd
from constants import *
from thermal_script import *


if __name__ == "__main__":
    # Read the Excel file into a DataFrame
    # XENON
    xe_melting_data = pd.read_excel(XE_DATA_FILEPATH, sheet_name=MELTING_SHEETNAME, header=2)
    xe_melting_data = xe_melting_data.filter(
        items=['Year', 'Author', 'Temperature', 'Pressure'])
    xe_melting_data = xe_melting_data.drop([0, 1], axis=0).reset_index(drop=True)
    xe_melting_data['Temperature'] = pd.to_numeric(
        xe_melting_data['Temperature'], errors='coerce')
    xe_melting_data['Pressure'] = pd.to_numeric(
        xe_melting_data['Pressure'], errors='coerce')
    
    # # Fit constants to the data per gas
    # gas_coefficients = fit_melting_pressure_single_gas(
    #     xe_melting_data, 'xenon', gas_params)
    # print("Fitted Constants for xenon:")
    # print(gas_coefficients)

    plot_melting_gas_data(xe_melting_data, 'xenon')
    
    
    

