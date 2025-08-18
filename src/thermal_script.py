import pandas as pd
import matplotlib.pyplot as plt
from pint import UnitRegistry
from constants import *
import numpy as np
from scipy.optimize import curve_fit
ureg = UnitRegistry()


def get_triple_point(gas_name: str, gas_params: dict):
    """
    Extracts the triple point temperature (T_t) and pressure (P_t) for a given gas.

    Parameters:
    - gas_name (str): The name of the gas (e.g., 'xenon', 'krypton', 'neon')
    - gas_params (dict): Dictionary mapping gas names to (T_t, P_t) tuples

    Returns:
    - Tuple (T_t, P_t) if gas is found in gas_params, else raises a KeyError.
    """
    if gas_name not in gas_params:
        raise KeyError(f"Gas '{gas_name}' not found in gas_params dictionary.")

    T_t, P_t = gas_params[gas_name]
    return T_t, P_t


def fit_melting_pressure_single_gas(data, gas_name: str, gas_params: dict):
    """
    Fit the melting pressure equation to experimental data for a single gas.
    Uses known T_t and P_t from gas_params.

    Parameters:
    - data: Pandas DataFrame with 'Temperature' and 'Pressure' columns (Pressure in Pa).
    - gas_name (str): Name of the gas (e.g., 'xenon')
    - gas_params (dict): Dictionary mapping gas names to (T_t, P_t) tuples

    Returns:
    - Dictionary of optimized parameters [e_4, e_5, e_6, e_7]
    """
    # Ensure correct data types
    data['Temperature'] = pd.to_numeric(data['Temperature'], errors='coerce')
    data['Pressure'] = pd.to_numeric(data['Pressure'], errors='coerce')

    # Clean the data
    data = data.dropna(subset=['Temperature', 'Pressure'])
    data = data[
        np.isfinite(data['Temperature']) &
        np.isfinite(data['Pressure'])
    ]

    T_data = data['Temperature'].values
    P_data = data['Pressure'].values / 1e6  # Convert to MPa

    if len(T_data) < 6:
        print("Not enough valid data to fit.")
        return None

    # Get T_t and P_t from the dictionary
    try:
        T_t, P_t = get_triple_point(gas_name, gas_params)
    except KeyError as e:
        print(e)
        return None

    # Initial guess: [e_4, e_5, e_6, e_7]
    initial_guess = [4053.39, 1.291, 2231.27, 0.84562]

    try:
        popt, _ = curve_fit(
            lambda T, e_4, e_5, e_6, e_7: melting_pressure_equation(
                T, e_4, e_5, e_6, e_7, T_t, P_t
            ),
            T_data, P_data,
            p0=initial_guess
        )
        return {
            'e_4': popt[0], 'e_5': popt[1],
            'e_6': popt[2], 'e_7': popt[3],
            'T_t': T_t, 'P_t': P_t
        }
    except RuntimeError:
        print("Curve fitting failed.")
        return None



def melting_pressure_equation(T, e_4, e_5, e_6, e_7, T_t, P_t):
    """
    Equation for melting pressure with safe power handling.
    """
    term1 = np.where(T / T_t - 1 >= 0, (T / T_t - 1) ** e_5, 0)
    term2 = np.where(T / T_t - 1 >= 0, (T / T_t - 1) ** e_7, 0)

    return P_t * (1 + e_4 * term1 + e_6 * term2)


def plot_melting_gas_data(data, gas_name: str):
    """
    plot_melting_gas_data Plotting Melting Pressure Data for a Gas

    Args:
        data (dataframe): DataFrame containing melting pressure data for a specific gas.
        gas_name (str): String representing the name of the gas (e.g., 'xenon', 'krypton', 'neon').
    """
    grouped = data.groupby(['Year', 'Author'])

    # Define markers and colors

    plt.figure(figsize=(10, 6))

    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            # Convert pressure to MPa
            group['Temperature'], group['Pressure'] / 1e6,
            s=50,
            label=f"{year}, {author}",
            edgecolors=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            facecolors='none',
            linewidths=1.5,
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
        )

    if gas_name == 'krypton':
        P_melt = melting_pressure_equation(
            data['Temperature'], KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t)
        X_MIN = MELTING_KRYPTON_X_MIN
        X_MAX = MELTING_KRYPTON_X_MAX
        Y_MIN = MELTING_KRYPTON_Y_MIN
        Y_MAX = MELTING_KRYPTON_Y_MAX
    elif gas_name == "xenon":
        P_melt = melting_pressure_equation(
            data['Temperature'], XENON_E_4, XENON_E_5, XENON_E_6, XENON_E_7, XENON_T_t, XENON_P_t)
        X_MIN = MELTING_XENON_X_MIN
        X_MAX = MELTING_XENON_X_MAX
        Y_MIN = MELTING_XENON_Y_MIN
        Y_MAX = MELTING_XENON_Y_MAX
    else:
        P_melt = melting_pressure_equation(
            data['Temperature'], NEON_E_4, NEON_E_5, NEON_E_6, NEON_E_7, NEON_T_t, NEON_P_t)
        X_MIN = MELTING_NEON_X_MIN
        X_MAX = MELTING_NEON_X_MAX
        Y_MIN = MELTING_NEON_Y_MIN
        Y_MAX = MELTING_NEON_Y_MAX
        
    # plt.scatter(data['Temperature'], P_melt, color='black',
    #          linestyle='-', linewidth=2, label='Calculated Melting Pressure')

    plt.plot(np.sort(data['Temperature']), np.sort(P_melt), color='black')

    # Set labels with italicized text to match LaTeX style in the image
    plt.xlabel(r'$\mathit{T}$ / K', fontsize=14)
    
    plt.ylabel(r'$\mathit{p}$ / MPa', fontsize=14)

    # Set a logarithmic scale for the y-axis
    plt.yscale('log')

    # Set x-axis limits similar to the uploaded image
    # plt.xlim(X_MIN, X_MAX)
    # plt.ylim(Y_MIN, Y_MAX)

    # Place the legend outside the plot
    plt.legend(
        loc='upper left',
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        fontsize=8,
        ncol=1  # Adjust the number of columns in the legend
    )
    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.title(f'Melting Pressure for {gas_name}', fontsize=14)

    output_filepath = f"{OUTPUT_FILEPATH}\{gas_name}_melting_temperatures_plot.png"
   # plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()
    
    
def plot_sublimation_gas_data(data, gas_name: str):
    """
    plot_sublimation_gas_data Plotting Sublimation Pressure Data for a Gas

    Args:
        data (dataframe): DataFrame containing sublimation pressure data for a specific gas.
        gas_name (str): String representing the name of the gas (e.g., 'xenon', 'krypton', 'neon').
    """
    grouped = data.groupby(['Year', 'Author'])

    # Define markers and colors

    plt.figure(figsize=(10, 6))

    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            # Convert pressure to MPa
            group['Temperature'], group['Pressure'] / 1e6,
            s=50,
            label=f"{year}, {author}",
            edgecolors=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            facecolors='none',
            linewidths=1.5,
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
        )

    if gas_name == 'krypton':
        P_sub = sublimation_pressure_equation(
            data['Temperature'], KRYPTON_E_1_SUB, KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t)
    elif gas_name == "xenon":
        P_sub = sublimation_pressure_equation(
            data['Temperature'], XENON_E_1_SUB, XENON_E_2_SUB, XENON_E_3_SUB, XENON_T_t, XENON_P_t)
    else:
        P_sub = sublimation_pressure_equation(
            data['Temperature'], NEON_E_1_SUB, NEON_E_2_SUB, NEON_E_3_SUB,NEON_T_t, NEON_P_t)

    # plt.scatter(data['Temperature'], P_melt, color='black',
    #          linestyle='-', linewidth=2, label='Calculated Melting Pressure')

    plt.plot(np.sort(data['Temperature']), np.sort(P_sub), color='black')

    # Set labels with italicized text to match LaTeX style in the image
    plt.xlabel(r'$\mathit{T}$ / K', fontsize=14)
    plt.ylabel(r'$\mathit{p}$ / MPa', fontsize=14)

    # Set a logarithmic scale for the y-axis
    plt.yscale('log')

    # Set x-axis limits similar to the uploaded image
    # plt.xlim(X_MIN, X_MAX)
    # plt.ylim(Y_MIN, Y_MAX)

    # Place the legend outside the plot
    plt.legend(
        loc='upper left',
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        fontsize=8,
        ncol=1  # Adjust the number of columns in the legend
    )
    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.title(f'Sublimation Pressure for {gas_name}', fontsize=14)

    output_filepath = f"{OUTPUT_FILEPATH}\{gas_name}_sublimation_plot.png"
    #plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()
    
def plot_gas_data(data, gas_name: str,y_axis:str,filename:str = None,y_var: str = 'Pressure',y_units: str = 'MPa'):
    """
    plot_gas_data Plot Gas Data for any variable without fitting a model.

    Args:
        data (_type_):   DataFrame containing gas data.
        gas_name (str):  Gas name to plot (e.g., 'xenon', 'krypton', 'neon').
        y_axis (str):  Variable to plot on the y-axis (e.g., 'Pressure', 'Density', etc.)
    """
    grouped = data.groupby(['Year', 'Author'])

    # Define markers and colors

    plt.figure(figsize=(10, 6))

    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            # Convert pressure to MPa
            group['Temperature'], group[y_axis],
            s=50,
            label=f"{year}, {author}",
            edgecolors=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            facecolors='none',
            linewidths=1.5,
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
        )

    # Set labels with italicized text to match LaTeX style in the image
    plt.xlabel(r'$\mathit{T}$ / K', fontsize=14)
    # plt.ylabel(f'$\\mathit{{{y_var}}}$ / {y_units}', fontsize=14)
    if "$" in y_var:
        label_str = f"{y_var} / {y_units}"
    else:
        ylabel_str = f"$\\mathit{{{y_var}}}$ / {y_units}"

        plt.ylabel(ylabel_str, fontsize=14)



    # Set a logarithmic scale for the y-axis
    plt.yscale('log')

    # Set x-axis limits similar to the uploaded image
    # plt.xlim(X_MIN, X_MAX)
    # plt.ylim(Y_MIN, Y_MAX)

    # Place the legend outside the plot
    plt.legend(
        loc='upper left',
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        fontsize=8,
        ncol=1  # Adjust the number of columns in the legend
    )
    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.title(f'{filename} for {gas_name}', fontsize=14)

    output_filepath = f"{OUTPUT_FILEPATH}\{gas_name}_{filename}_plot.png"
    #plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()


def plot_melting_pressure_deviation(data, gas_name):
    """
    Plot the percentage deviation of calculated melting pressure from experimental data
    for a given gas.
    """
    data['Pressure'] = pd.to_numeric(data['Pressure']/1e6, errors='coerce')

    # Clean and filter numeric data
    # Select constants for the gas
    if gas_name == 'krypton':
        e_4, e_5, e_6, e_7 = KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7
        T_t, P_t = KRYPTON_T_t, KRYPTON_P_t
        X_MIN = MELTING_KRYPTON_X_MIN
        X_MAX = MELTING_KRYPTON_X_MAX
    elif gas_name == 'xenon':
        e_4, e_5, e_6, e_7 = XENON_E_4, XENON_E_5, XENON_E_6, XENON_E_7
        T_t, P_t = XENON_T_t, XENON_P_t
        X_MIN = MELTING_XENON_X_MIN
        X_MAX = MELTING_XENON_X_MAX
    else:
        e_4, e_5, e_6, e_7 = NEON_E_4, NEON_E_5, NEON_E_6, NEON_E_7
        T_t, P_t = NEON_T_t, NEON_P_t
        X_MIN = MELTING_NEON_X_MIN
        X_MAX = MELTING_NEON_X_MAX

    # Calculate theoretical melting pressure
    data['P_calc'] = melting_pressure_equation(
        data['Temperature'], e_4, e_5, e_6, e_7, T_t, P_t)

    # Calculate % deviation
    data['Deviation_percent'] = 100 * \
        (data['Pressure'] - data['P_calc']) / \
        data['Pressure']

    # Plot grouped by year-author
    grouped = data.groupby(['Year', 'Author'])

    plt.figure(figsize=(10, 6))
    plt.tick_params(direction='in', top=True, right=True)
    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            group['Temperature'], group['Deviation_percent'],
            label=f"{year}, {author}",
            s=MARKERSIZE,  # Marker size
            # Assigning custom color
            edgecolor=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            # Assigning custom marker
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
            # marker=markers[i % len(markers)],  # Marker style
            facecolors='none',  # Open symbols
        )
    plt.axhline(0, color='black', linewidth=1)
    plt.xlabel(r'$\mathit{T}$ / K', fontsize=14)
    plt.ylabel(
        r'$100 \cdot (\mathit{p}_{\mathrm{exp}} - \mathit{p}_{\mathrm{calc}})/\mathit{p}_{\mathrm{exp}}$', fontsize=14)
    plt.title(f'Melting Pressure Deviation for {gas_name}', fontsize=14)
    plt.xlim(X_MIN, X_MAX)
    plt.ylim(MELTING_DEVIATION_Y_MIN, MELTING_DEVIATION_Y_MAX)
    # plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    output_filepath = f"{OUTPUT_FILEPATH}\\{gas_name}_melting_pressure_deviation.png"
    #
    # plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()

def read_melting_data(filepath, sheet_name):
    """
    Reads and cleans melting temperature and pressure data from an Excel file.

    Parameters:
    - filepath: Path to the Excel file
    - sheet_name: Sheet name containing the melting data

    Returns:
    - Cleaned Pandas DataFrame with columns: ['Year', 'Author', 'Temperature', 'Pressure']
    """
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=2)
    df = df.filter(items=['Year', 'Author', 'Temperature', 'Pressure'])
    df = df.drop([0, 1], axis=0).reset_index(drop=True)
    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df['Pressure'] = pd.to_numeric(df['Pressure'], errors='coerce')
    return df


def read_fusion_data(filepath, sheet_name):
    """
    Reads and cleans melting temperature and pressure data from an Excel file.

    Parameters:
    - filepath: Path to the Excel file
    - sheet_name: Sheet name containing the melting data

    Returns:
    - Cleaned Pandas DataFrame with columns: ['Year', 'Author', 'Temperature', 'Pressure']
    """
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=2)
    df = df.filter(
        items=['Year', 'Author', 'Temperature', 'Change in Enthalpy'])
    df = df.drop([0, 1], axis=0).reset_index(drop=True)
    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df['Change in Enthalpy'] = pd.to_numeric(df['Change in Enthalpy'], errors='coerce')
    return df

def read_thermal_coeff_data(filepath, sheet_name):
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=3)
    df = df.filter(
        items=['Year', 'Author', 'Temperature', 'Thermal Expansion Coefficient'])
    df = df.drop([0, 1], axis=0).reset_index(drop=True)
    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df['Thermal Expansion Coefficient'] = pd.to_numeric(
        df['Thermal Expansion Coefficient'], errors='coerce')
    return df


def read_data(filepath:str, sheet_name:str, variable:str,header_row=2):
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=header_row)
    df = df.filter(
        items=['Year', 'Author', 'Temperature', variable])
    df = df.drop([0, 1], axis=0).reset_index(drop=True)
    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df[variable] = pd.to_numeric(
        df[variable], errors='coerce')
    return df


def sublimation_pressure_equation(T, e1, e2, e3, T_t, P_t):
    """
    Equation for sublimation pressure.

    Parameters:
    - T: Temperature in Kelvin (array-like)
    - e1, e2, e3: Empirical constants
    - T_t: Triple point temperature (K)
    - P_t: Triple point pressure (MPa)

    Returns:
    - Sublimation pressure (MPa)
    """
    theta = 1 - (T / T_t)
    exponent = (T_t / T) * (e1 * theta + e2 * theta ** (3/2) + e3 * theta ** 5)
    return P_t * np.exp(exponent)


def fit_sublimation_pressure_single_gas(data, gas_name: str, gas_params: dict):
    """
    Fit sublimation pressure data using the sublimation pressure equation,
    emphasizing low-pressure and low-temperature data.
    """
    try:
        T_t, P_t = get_triple_point(gas_name, gas_params)
    except KeyError as e:
        print(e)
        return None

    # Filter data below triple point
    data = data[data['Temperature'] < T_t]
    data = data.dropna(subset=['Temperature', 'Pressure'])

    T_data = data['Temperature'].values
    P_data = data['Pressure'].values / 1e6  # Pa → MPa

    print(
        f"Fitting sublimation pressure for {gas_name} with {len(T_data)} data points.")

    if len(T_data) < 6:
        print("Not enough data points to fit.")
        return None

    # Custom weighting: emphasize low T and P
    weights = 1 / (T_data**2 * P_data**2 + 1e-12)  # Avoid div-by-zero
    
    # weights = 1 / (P_data + 1e-9)  # simpler and stabler
    # weights /= np.max(weights)    # normalize for numerical stability

    initial_guess = [-18.910506919795417, 11.285160377819482, -12.280765945817006]

    try:
        popt, _ = curve_fit(
            lambda T, e1, e2, e3: sublimation_pressure_equation(
                T, e1, e2, e3, T_t, P_t),
            T_data,
            P_data,
            p0=initial_guess,
            sigma=1/weights,
            absolute_sigma=False,
            maxfev=10000
        )
        return {
            'e1': popt[0],
            'e2': popt[1],
            'e3': popt[2],
            'T_t': T_t,
            'P_t': P_t
        }
    except RuntimeError:
        print("Sublimation curve fitting failed.")
        return None


def log_sublimation_pressure_equation(T, e1, e2, e3, T_t, P_t):
    T = np.asarray(T)
    theta = 1 - T / T_t
    safe_theta = np.where(theta >= 0, theta, 0.0)
    exponent = (T_t / T) * (
        e1 * safe_theta +
        e2 * safe_theta ** (3 / 2) +
        e3 * safe_theta ** 5
    )
    return np.log(P_t) + exponent


def fit_log_sublimation_pressure_single_gas(data, gas_name: str, gas_params: dict):
    """
    Fit log of sublimation pressure data using the log-transformed equation.
    """
    try:
        T_t, P_t = get_triple_point(gas_name, gas_params)
    except KeyError as e:
        print(e)
        return None

    # Filter below triple point and clean
    data = data[data['Temperature'] < T_t].dropna(
        subset=['Temperature', 'Pressure'])
    T_data = data['Temperature'].values
    P_data = data['Pressure'].values / 1e6  # Convert to MPa
    log_P_data = np.log(P_data)

    print(
        f"Log-fitting sublimation pressure for {gas_name} with {len(T_data)} points.")

    if len(T_data) < 6:
        print("Not enough data.")
        return None

    # Optional weighting (mild)
    weights = 1 / (T_data + 1e-3)
    weights /= np.max(weights)

    initial_guess = [-10.0, -3.0, 40.0]  # Start mild
    bounds = ([-100, -100, -100], [100, 100, 100])  # Keep reasonable

    try:
        popt, _ = curve_fit(
            lambda T, e1, e2, e3: log_sublimation_pressure_equation(
                T, e1, e2, e3, T_t, P_t),
            T_data,
            log_P_data,
            p0=initial_guess,
            bounds=bounds,
            sigma=1 / weights,
            absolute_sigma=False,
            maxfev=10000
        )

        return {
            'e1': popt[0],
            'e2': popt[1],
            'e3': popt[2],
            'T_t': T_t,
            'P_t': P_t
        }

    except RuntimeError:
        print("Log-sublimation curve fitting failed.")
        return None

def read_bulk_modulus_data(filepath, sheet_name):
    """
    Reads and cleans bulk modulus data from an Excel file.

    Parameters:
    - filepath: Path to the Excel file
    - sheet_name: Sheet name containing the bulk modulus data

    Returns:
    - Cleaned Pandas DataFrame with columns: ['Year', 'Author', 'Temperature', 'Bulk Modulus']
    """
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=2)

    
    # Split the data into two DataFrames
    # DataFrame 1: Left table (up to "Pressure")
    df_left = df.iloc[:, 0:7]  # adjust if more columns needed

    # DataFrame 2: Right table (starts from "Reference" on the right side)
    df_right = df.iloc[:, 7:]  # assumes right table starts from column 8
    
    df_left = df_left.filter(items=['Year', 'Author', 'Temperature', 'Beta T '])
    df_left = df_left.drop([0, 1], axis=0).reset_index(drop=True)
    df_left['Temperature'] = pd.to_numeric(
        df_left['Temperature'], errors='coerce')
    df_left['Beta T '] = pd.to_numeric(
        df_left['Beta T '], errors='coerce')

    df_right = df_right.filter(items=['Year', 'Author', 'Temperature', 'Beta S '])
    df_right = df_right.drop([0, 1], axis=0).reset_index(drop=True)
    df_right['Temperature'] = pd.to_numeric(df_right['Temperature'], errors='coerce')
    df_right['Beta S '] = pd.to_numeric(
        df_right['Beta S '], errors='coerce')

    #rename Beta T and Beta S to remove white space at the end
    df_left.rename(columns={'Beta T ': 'Beta T'}, inplace=True)
    df_right.rename(columns={'Beta S ': 'Beta S'}, inplace=True)
    return df_left, df_right


def read_cell_volume_data(filepath, sheet_name):
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=2)

    # Split into left and right tables
    df_left = df.iloc[:, 0:9]
    df_right = df.iloc[:, 9:]

    # === Clean Left Table ===
    df_left = df_left.filter(
        items=['Year', 'Author', 'Temperature', 'Cell Volume '])
    df_left = df_left.drop([0, 1], axis=0).reset_index(drop=True)
    df_left['Temperature'] = pd.to_numeric(
        df_left['Temperature'], errors='coerce')
    df_left['Cell Volume '] = pd.to_numeric(
        df_left['Cell Volume '], errors='coerce')
    df_left.rename(columns={'Cell Volume ': 'Cell Volume'}, inplace=True)

    # === Clean Right Table ===
    df_right.rename(columns={'Year.1': 'Year'}, inplace=True)
    df_right = df_right.rename(columns={
        'Author.1': 'Author',
        'Temperature.1': 'Temperature',
        'Cell Volume .1': 'Cell Volume '
    })
    df_right = df_right.filter(
        items=['Year', 'Author', 'Temperature', 'Cell Volume '])
    df_right = df_right.drop([0, 1], axis=0).reset_index(drop=True)

    # Drop all rows *at or after* the first one containing "High Pressure"
    high_pressure_idx = df_right[df_right.apply(
        lambda row: row.astype(str).str.contains("high pressure", case=False).any(), axis=1
    )].index

    if not high_pressure_idx.empty:
        df_right = df_right.loc[:high_pressure_idx[0] -
                                1].reset_index(drop=True)

    # Final type conversion and cleanup
    df_right['Temperature'] = pd.to_numeric(
        df_right['Temperature'], errors='coerce')
    df_right['Cell Volume '] = pd.to_numeric(
        df_right['Cell Volume '], errors='coerce')
    df_right.rename(columns={'Cell Volume ': 'Cell Volume'}, inplace=True)

    return df_left, df_right


def plot_all_gas_properties(data, gas_name):
    print(f"Plotting properties for: {gas_name.capitalize()}")

    # Phase-Specific
    if gas_name == 'neon':
        plot_melting_gas_data(data['melting'], gas_name)
    if gas_name == 'krypton':
        plot_melting_pressure_deviation(data['melting'], gas_name)

    # Fusion
    plot_gas_data(data['fusion'], gas_name, 'Change in Enthalpy',
                  'Heat_of_Fusion', r'$\Delta H$', r'kJ/mol')

    # Heat of Sublimation
    plot_gas_data(data['heatsub'], gas_name, 'Change in Enthalpy',
                  'Heat_of_Sublimation', r'$\Delta H$', r'kJ/mol')

    # Sublimation Pressure
    print("Plotting sublimation pressure data...")
    plot_sublimation_gas_data(data['sublimation'], gas_name)

    # Thermal Expansion Coefficient
    plot_gas_data(data['thermal_coeff'], gas_name,
                  'Thermal Expansion Coefficient',
                  'Thermal Expansion Coefficient', r'$\alpha$', r'$1/K$')

    # Heat Capacity
    plot_gas_data(data['heat_capacity'], gas_name,
                  'Heat Capacity', 'Heat Capacity', r'$C_p$', r'J/(mol\cdot K)')

    # Bulk Modulus
    plot_gas_data(data['bulk_s'], gas_name, 'Beta S',
                  'Bulk Modulus S', r'$B_s$', r'$1/\mathrm{MPa}$')
    plot_gas_data(data['bulk_t'], gas_name, 'Beta T',
                  'Bulk Modulus T', r'$B_t$', r'$1/\mathrm{MPa}$')

    # Cell Volume
    plot_gas_data(data['cell_volume_sub'], gas_name, 'Cell Volume',
                  'Sublimation Curve', r'$V_{\mathrm{cell}}$', r'$\mathrm{cm}^3/\mathrm{mol}$')
    plot_gas_data(data['cell_volume_melt'], gas_name, 'Cell Volume',
                  'Melting Curve', r'$V_{\mathrm{cell}}$', r'$\mathrm{cm}^3/\mathrm{mol}$')

    # Combined Plot
    try:
        combined = data['cell_volume_sub'].copy()
        combined['V'] = data['cell_volume_melt']['V'].values
        plot_gas_data(combined, gas_name, 'Cell Volume',
                      'Melting and Sublimation Curve', r'$V$', r'$\mathrm{cm}^3/\mathrm{mol}$')
    except Exception:
        print("Warning: Could not combine melt + sublimation cell volume data.")


def add_psub_column(data: pd.DataFrame, gas_name: str) -> pd.DataFrame:
    """
    Adds sublimation pressure (MPa) column to the given DataFrame
    based on the gas type and its empirical constants.
    
    Args:
        data (pd.DataFrame): must have 'Temperature' column
        gas_name (str): 'krypton', 'xenon', or 'neon'
    
    Returns:
        pd.DataFrame: input data with a new 'P_sub' column (MPa)
    """
    if gas_name == 'krypton':
        P_sub = sublimation_pressure_equation(
            data['Temperature'],
            KRYPTON_E_1_SUB, KRYPTON_E_2_SUB, KRYPTON_E_3_SUB,
            KRYPTON_T_t, KRYPTON_P_t
        )
    elif gas_name == "xenon":
        P_sub = sublimation_pressure_equation(
            data['Temperature'],
            XENON_E_1_SUB, XENON_E_2_SUB, XENON_E_3_SUB,
            XENON_T_t, XENON_P_t
        )
    else:  # default = neon
        P_sub = sublimation_pressure_equation(
            data['Temperature'],
            NEON_E_1_SUB, NEON_E_2_SUB, NEON_E_3_SUB,
            NEON_T_t, NEON_P_t
        )

    # Add new column (MPa)
    data = data.copy()
    data['Pressure'] = P_sub
    return data


def add_pmelt_column(data: pd.DataFrame, gas_name: str) -> pd.DataFrame:
    """
    Adds melting pressure (MPa) column to the given DataFrame
    based on the gas type and its empirical constants.

    Args:
        data (pd.DataFrame): must have 'Temperature' column
        gas_name (str): 'krypton', 'xenon', or 'neon'

    Returns:
        pd.DataFrame: input data with a new 'P_melt' column (MPa)
    """
    if gas_name == 'krypton':
        P_melt = melting_pressure_equation(
            data['Temperature'],
            KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7,
            KRYPTON_T_t, KRYPTON_P_t
        )
    elif gas_name == "xenon":
        P_melt = melting_pressure_equation(
            data['Temperature'],
            XENON_E_4, XENON_E_5, XENON_E_6, XENON_E_7,
            XENON_T_t, XENON_P_t
        )
    else:  # default = neon
        P_melt = melting_pressure_equation(
            data['Temperature'],
            NEON_E_4, NEON_E_5, NEON_E_6, NEON_E_7,
            NEON_T_t, NEON_P_t
        )

    # Add new column (MPa)
    data = data.copy()
    data['Pressure'] = P_melt
    return data
