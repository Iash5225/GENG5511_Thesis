import pandas as pd
from constants import *
from thermal_script import *
from pathlib import Path

def load_all_gas_data(gas_name, read_from_excel=True):
    """
    Returns a dictionary of all thermophysical datasets for the specified gas.

    Parameters:
        gas_name (str): 'neon', 'xenon', or 'krypton'
        read_from_excel (bool): If True, reads from Excel; otherwise, from TXT.

    Returns:
        dict[str, pd.DataFrame]: Dictionary of DataFrames for each dataset.
    """
    base_file = {
        'neon': NE_DATA_FILEPATH,
        'xenon': XE_DATA_FILEPATH,
        'krypton': KR_DATA_FILEPATH
    }[gas_name]

    data = {}

    if read_from_excel:
        data['melting'] = pd.read_excel(
            base_file, sheet_name=MELTING_SHEETNAME)
        data['sublimation'] = read_melting_data(
            base_file, SUBLIMATION_SHEETNAME)
        data['fusion'] = read_fusion_data(base_file, FUSION_SHEETNAME)
        data['heatsub'] = read_fusion_data(
            base_file, HEAT_OF_SUBLIMATION_SHEETNAME)
        data['thermal_coeff'] = read_thermal_coeff_data(
            base_file, THERMAL_COEFFICIENT_SHEETNAME)
        data['heat_capacity'] = read_data(
            base_file, HEAT_CAPACITY_SHEETNAME, "Heat Capacity", 2)
        data['bulk_t'], data['bulk_s'] = read_bulk_modulus_data(
            base_file, BULK_MODULUS_SHEETNAME)
        sub_vol, melt_vol = read_cell_volume_data(
            base_file, CELL_VOLUME_SHEETNAME)
        data['cell_volume_sub'] = sub_vol
        data['cell_volume_melt'] = melt_vol
    else:
        path = Path(TXT_DATA_FILEPATH)

        def txt(gas, key):
            return path / f"{gas}_{key}.txt"

        data['melting'] = pd.read_csv(txt(gas_name, 'melting_data_for_fitting'), sep='\t')
        data['sublimation'] = pd.read_csv(
            txt(gas_name, 'sublimation_data_for_fitting'), sep='\t')
        data['fusion'] = pd.read_csv(txt(gas_name, 'fusion_data'), sep='\t')
        data['heatsub'] = pd.read_csv(
            txt(gas_name, 'heat_of_sublimation_data_for_fitting'), sep='\t')
        data['thermal_coeff'] = pd.read_csv(
            txt(gas_name, 'thermal_coeff_data'), sep='\t')
        data['heat_capacity'] = pd.read_csv(
            txt(gas_name, 'heat_capacity_data'), sep='\t')
        data['bulk_s'] = pd.read_csv(txt(gas_name, 'bulk_modulus_s'), sep='\t')
        data['bulk_t'] = pd.read_csv(txt(gas_name, 'bulk_modulus_t'), sep='\t')
        data['cell_volume_sub'] = pd.read_csv(
            txt(gas_name, 'cell_volume_sub'), sep='\t')
        data['cell_volume_melt'] = pd.read_csv(
            txt(gas_name, 'cell_volume_melt'), sep='\t')

    return data
