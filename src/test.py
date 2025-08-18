import pandas as pd
from constants import *
from thermal_script import *
READ_FROM_EXCEL = False

if __name__ == "__main__":
    # Read the Excel file into a DataFrame

    if READ_FROM_EXCEL:
        ne_heatsub_data = read_fusion_data(
            NE_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
        xe_heatsub_data = read_fusion_data(
            XE_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
        kr_heatsub_data = read_fusion_data(
            KR_DATA_FILEPATH, HEAT_OF_SUBLIMATION_SHEETNAME)
    else:
        # Enthalpy of Sublimation
        ne_heatsub_data = pd.read_csv(
            f"{TXT_DATA_FILEPATH}\\neon_heat_of_sublimation_data.txt", sep='\t')
        xe_heatsub_data = pd.read_csv(
            f"{TXT_DATA_FILEPATH}\\xenon_heat_of_sublimation_data.txt", sep='\t')
        kr_heatsub_data = pd.read_csv(
            f"{TXT_DATA_FILEPATH}\\krypton_heat_of_sublimation_data.txt", sep='\t')

    # Save to txt files

    # Heat of Sublimation Data
    # ne_heatsub_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\neon_heat_of_sublimation_data.txt", sep='\t', index=False)
    # xe_heatsub_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\xenon_heat_of_sublimation_data.txt", sep='\t', index=False)
    # kr_heatsub_data.to_csv(
    #     f"{TXT_DATA_FILEPATH}\\krypton_heat_of_sublimation_data.txt", sep='\t', index=False)

