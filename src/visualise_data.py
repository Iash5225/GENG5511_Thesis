import pandas as pd
import matplotlib.pyplot as plt
from pint import UnitRegistry
from constants import *
import numpy as np
from scipy.optimize import curve_fit
ureg = UnitRegistry()

# Functions
def atm_to_mpa(pressure_atm):
    """
    Convert pressure from atm to MPa using the pint library.

    Parameters:
    pressure_atm (float or pd.Series): Pressure in atm

    Returns:
    float or pd.Series: Pressure in MPa
    """
    # Define the unit
    pressure = pressure_atm * ureg.atm
    # Convert to MPa
    return pressure.to(ureg.megapascal).magnitude

def bar_to_mpa(pressure_bar):
    """
    Convert pressure from bar to MPa using the pint library.
    
    Parameters:
    pressure_bar (float or pd.Series): Pressure in bar
    
    Returns:
    float or pd.Series: Pressure in MPa
    """
    # Define the unit
    pressure = pressure_bar * ureg.bar
    # Convert to MPa
    return pressure.to(ureg.megapascal).magnitude


def cal_to_joules(heat_capacity_cal):
    return heat_capacity_cal*4.184

def convert_pressure(data):
    data['Pressure_MPa_from_atm'] = data['Pressure_atm'].apply(
        lambda x: atm_to_mpa(x) if pd.notnull(x) else None)
    data['Pressure_MPa_from_bar'] = data['Pressure_bar'].apply(
        lambda x: bar_to_mpa(x) if pd.notnull(x) else None)
    data['Pressure_MPa_from_kbar'] = data['Pressure_kbar'].apply(
        lambda x: x * 100 if pd.notnull(x) else None)

    data['Pressure_Mpa'] = data['Pressure_MPa_from_atm'].combine_first(
        data['Pressure_MPa_from_bar']).combine_first(data['Pressure_MPa_from_kbar'])

    return data.drop(columns=['Pressure_MPa_from_atm', 'Pressure_MPa_from_bar', 'Pressure_MPa_from_kbar'])


def convert_heat_capacity(data):
    data['cal_to_joules'] = data['cp_cal/mol/deg'].apply(
        lambda x: cal_to_joules(x) if pd.notnull(x) else None)
    data['millijoules_to_joules'] = data['cp_mJ/mol/K'].apply(
        lambda x: x / 1000 if pd.notnull(x) else None)

    data['cp_J/mol/K'] = data['cp_J/mol/K'].combine_first(
        data['cal_to_joules']).combine_first(data['millijoules_to_joules'].combine_first(data['cp_Ws/mol/deg']))
    return data
    

def plot_melting_gas_data(data, gas_name):
    gas_data = data[data['gas'] == gas_name]
    grouped = gas_data.groupby(['year', 'author'])

    # Define markers and colors

    plt.figure(figsize=(10, 6))

    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            group['Temperature_Kelvin'], group['Pressure_Mpa'],
            s=50,  # Smaller marker size
            label=f"{year}, {author}",
            # Assigning custom color
            edgecolor=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            # Assigning custom marker
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
            # marker=markers[i % len(markers)],  # Unique marker style
            alpha=0.8,  # Slight transparency for overlapping markers
            facecolors='none',  # Open symbols
        )

    if gas_name == 'krypton':
        P_melt = melting_pressure_equation(
            gas_data['Temperature_Kelvin'], KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7,KRYPTON_T_t,KRYPTON_P_t)
        X_MIN = MELTING_KRYPTON_X_MIN
        X_MAX = MELTING_KRYPTON_X_MAX
        Y_MIN = MELTING_KRYPTON_Y_MIN
        Y_MAX= MELTING_KRYPTON_Y_MAX
        # plt.plot(np.sort(gas_data['Temperature_Kelvin']),
        #          np.sort(P_melt), color='black')
    elif gas_name == "xenon":
        P_melt = melting_pressure_equation(
            gas_data['Temperature_Kelvin'], XENON_E_4, XENON_E_5, XENON_E_6, XENON_E_7, XENON_T_t, XENON_P_t)
        X_MIN = MELTING_XENON_X_MIN
        X_MAX = MELTING_XENON_X_MAX
        Y_MIN = MELTING_XENON_Y_MIN
        Y_MAX = MELTING_XENON_Y_MAX
        # plt.plot(np.sort(gas_data['Temperature_Kelvin']),
        #          np.sort(P_melt), color='black')
    else:
        P_melt = melting_pressure_equation(
            gas_data['Temperature_Kelvin'], NEON_E_4, NEON_E_5, NEON_E_6, NEON_E_7, NEON_T_t, NEON_P_t)
        X_MIN = MELTING_NEON_X_MIN
        X_MAX = MELTING_NEON_X_MAX
        Y_MIN = MELTING_NEON_Y_MIN
        Y_MAX = MELTING_NEON_Y_MAX
        # plt.plot(np.sort(gas_data['Temperature_Kelvin']),
        #          np.sort(P_melt), color='black')
    
    # plt.scatter(gas_data['Temperature_Kelvin'], P_melt, color='black',
    #          linestyle='-', linewidth=2, label='Calculated Melting Pressure')
    plt.plot(np.sort(gas_data['Temperature_Kelvin']),
             np.sort(P_melt), color='black')

    # Set labels with italicized text to match LaTeX style in the image
    plt.xlabel(r'$\mathit{T}$ / K', fontsize=14)
    plt.ylabel(r'$\mathit{p}$ / MPa', fontsize=14)

    # Set a logarithmic scale for the y-axis
    plt.yscale('log')

    # Set x-axis limits similar to the uploaded image
    plt.xlim(X_MIN, X_MAX)
    plt.ylim(Y_MIN, Y_MAX)

    # Use minimal grid style
    plt.grid(True, which="both", linestyle='--', alpha=0.3)

    # Place the legend outside the plot
    plt.legend(
        loc='upper left',
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        fontsize=8,
        ncol=1  # Adjust the number of columns in the legend
    )
    plt.grid(True, linestyle='--', alpha=0.6)
    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    plt.title(f'Melting Temperature for {gas_name}', fontsize=14)

    output_filepath = f"{OUTPUT_FILEPATH}\{gas_name}_melting_temperatures_plot.png"
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()


def melting_pressure_equation(T, e_4, e_5, e_6, e_7,T_t,P_t):
    """
    Equation for melting pressure with safe power handling.
    """
    term1 = np.where(T / T_t - 1 >= 0, (T / T_t - 1) ** e_5, 0)
    term2 = np.where(T / T_t - 1 >= 0, (T / T_t - 1) ** e_7, 0)

    return P_t * (1 + e_4 * term1 + e_6 * term2)


def fit_melting_pressure(data, gas_params):
    """
    Fit the melting pressure equation to experimental data to find e_4, e_5, e_6, and e_7.
    Returns a dictionary with gas names as keys and optimized parameters as values.
    
    Parameters:
    - data: Pandas DataFrame containing melting data.
    - gas_params: Dictionary mapping gas names to their (T_t, P_t) values.
    """
    gas_coefficients = {}

    for gas in data['gas'].unique():
        gas_data = data[data['gas'] == gas].copy()

        # Ensure numerical data types
        gas_data['Temperature_Kelvin'] = pd.to_numeric(
            gas_data['Temperature_Kelvin'], errors='coerce')
        gas_data['Pressure_Mpa'] = pd.to_numeric(
            gas_data['Pressure_Mpa'], errors='coerce')

        # Drop NaN and infinite values
        gas_data = gas_data.dropna(
            subset=['Temperature_Kelvin', 'Pressure_Mpa'])
        gas_data = gas_data[(np.isfinite(gas_data['Temperature_Kelvin'].values)) & (
            np.isfinite(gas_data['Pressure_Mpa'].values))]

        # Extract cleaned arrays
        T_data = gas_data['Temperature_Kelvin'].values
        P_data = gas_data['Pressure_Mpa'].values

        # Ensure enough valid data points
        if len(T_data) < 4:
            print(f"Not enough valid data to fit for {gas}. Skipping.")
            gas_coefficients[gas] = None
            continue

        # Get T_t and P_t for this gas
        if gas not in gas_params:
            print(
                f"Skipping {gas} because no T_t and P_t values were provided.")
            continue
        T_t, P_t = gas_params[gas]

        # Initial guesses for e_4, e_5, e_6, e_7
        initial_guess = [1500, 1.7, 4600, 0.98]

        try:
            popt, _ = curve_fit(lambda T, e_4, e_5, e_6, e_7:
                                melting_pressure_equation(
                                    T, e_4, e_5, e_6, e_7, T_t, P_t),
                                T_data, P_data, p0=initial_guess)
            gas_coefficients[gas] = popt
        except RuntimeError:
            print(f"Could not fit curve for {gas}")
            gas_coefficients[gas] = None

    return gas_coefficients


def plot_thermal_coefficient_gas_data(data, gas_name):
    gas_data = data[data['gas'] == gas_name]
    grouped = gas_data.groupby(['year', 'author'])

    # Define markers and colors

    plt.figure(figsize=(10, 6))
    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            group['Temperature_Kelvin'], group['Alpha_Kelvin^-1'],
            label=f"{year}, {author}",
            s=MARKERSIZE,  # Marker size
            # Assigning custom color
            edgecolor=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            # Assigning custom marker
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
            # marker=markers[i % len(markers)],  # Marker style
            facecolors='none',  # Open symbols
        )

    plt.xlabel(r'$\mathit{T}$ / K', fontsize=AXIS_FONT_SIZE)
    plt.ylabel(r'$\mathit{\alpha}$ / $K^{-1}$', fontsize=AXIS_FONT_SIZE)
    plt.title(f'Thermal Coefficient for {gas_name}', fontsize=AXIS_FONT_SIZE)

    # Set a logarithmic scale for the y-axis
    plt.yscale('log')

    # Place the legend outside the plot
    plt.legend(
        loc='upper left',
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        fontsize=8,
        ncol=1  # Adjust the number of columns in the legend
    )
    plt.grid(True, linestyle='--', alpha=0.6)
    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    output_filepath = f"{OUTPUT_FILEPATH}\{gas_name}_thermal_coefficient_plot.png"
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()


def plot_melting_pressure_deviation(data, gas_name):
    """
    Plot the percentage deviation of calculated melting pressure from experimental data
    for a given gas.
    """
    gas_data = data[data['gas'] == gas_name].copy()

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
    gas_data['P_calc'] = melting_pressure_equation(
        gas_data['Temperature_Kelvin'], e_4, e_5, e_6, e_7, T_t, P_t)
    
    print(gas_data['P_calc'])

    # Calculate % deviation
    gas_data['Deviation_percent'] = 100 * \
        (gas_data['Pressure_Mpa'] - gas_data['P_calc']) / \
        gas_data['Pressure_Mpa']

    # Plot grouped by year-author
    grouped = gas_data.groupby(['year', 'author'])

    plt.figure(figsize=(10, 6))
    for i, ((year, author), group) in enumerate(grouped):
        plt.plot(group['Temperature_Kelvin'], group['Deviation_percent'],
                 marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
                 linestyle='None',
                 markersize=10,
                 label=f"{year}, {author}",
                 color=CUSTOMCOLORS[i % len(CUSTOMCOLORS)])

    plt.axhline(0, color='black', linewidth=1)
    plt.xlabel(r'$\mathit{T}$ / K', fontsize=14)
    plt.ylabel(
        r'$100 \cdot (\mathit{p}_{\mathrm{exp}} - \mathit{p}_{\mathrm{calc}})/\mathit{p}_{\mathrm{exp}}$', fontsize=14)
    plt.title(f'Melting Pressure Deviation for {gas_name}', fontsize=14)
    plt.xlim(X_MIN, X_MAX)
    plt.ylim(-10, 20)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    output_filepath = f"{OUTPUT_FILEPATH}\\{gas_name}_melting_pressure_deviation.png"
    plt.savefig(output_filepath, dpi=300)
    plt.show()


def plot_heat_capacity_gas_data(data, gas_name):
    gas_data = data[data['gas'] == gas_name]
    grouped = gas_data.groupby(['year', 'author'])

    # Define markers and colors

    plt.figure(figsize=(10, 6))
    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            group['Temperature_Kelvin'], group['cp_J/mol/K'],
            label=f"{year}, {author}",
            s=MARKERSIZE,  # Marker size
            # Assigning custom color
            edgecolor=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            # Assigning custom marker
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
            # marker=markers[i % len(markers)],  # Marker style
            facecolors='none',  # Open symbols
        )

    plt.xlabel(r'$\mathit{T}$ / K', fontsize=AXIS_FONT_SIZE)
    plt.ylabel(r'$\mathit{c_p}$ / $Jmol^{-1}K^{-1}$', fontsize=AXIS_FONT_SIZE)
    plt.title(f'Heat Capacity for {gas_name}', fontsize=AXIS_FONT_SIZE)

    # Set a logarithmic scale for the y-axis
    plt.yscale('log')

    # Place the legend outside the plot
    plt.legend(
        loc='upper left',
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        fontsize=8,
        ncol=1  # Adjust the number of columns in the legend
    )
    plt.grid(True, linestyle='--', alpha=0.6)
    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    output_filepath = f"{OUTPUT_FILEPATH}\{gas_name}_heat_capacity_plot.png"
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()


# Main Execution
if __name__ == "__main__":
    

    # Read the Excel file into a DataFrame
    melting_data = pd.read_excel(FILEPATH, sheet_name=MELTING_SHEETNAME, header=1)
    thermalcoefficient_data = pd.read_excel(
        FILEPATH, sheet_name=THERMAL_COEFFICIENT_SHEETNAME, header=1)
    
    heat_capacity_data = pd.read_excel(
        FILEPATH, sheet_name=HEAT_CAPACITY_SHEETNAME, header=1)
    
    
    

    melting_data = melting_data.drop(melting_data.columns[0], axis=1)
    # data = data.dropna()  # Drop rows with missing values
    melting_data['Temperature_Kelvin'] = pd.to_numeric(
        melting_data['Temperature_Kelvin'], errors='coerce')
    melting_data['Pressure_atm'] = pd.to_numeric(
        melting_data['Pressure_atm'], errors='coerce')
    thermalcoefficient_data = thermalcoefficient_data.drop(
        thermalcoefficient_data.columns[0], axis=1)
    heat_capacity_data = heat_capacity_data.drop(
        heat_capacity_data.columns[0], axis=1)

    # Convert Pressures
    melting_data = convert_pressure(melting_data)

    heat_capacity_data = convert_heat_capacity(heat_capacity_data)

    # Fit constants to the data per gas
    gas_coefficients = fit_melting_pressure(melting_data, gas_params)
    for gas, params in gas_coefficients.items():
        if params is not None:
            print(f"Optimized Constants for {gas}: e_4={params[0]:.2f}, e_5={params[1]:.3f}, e_6={params[2]:.2f}, e_7={params[3]:.5f}")

    plot_heat_capacity_gas_data(heat_capacity_data, 'krypton')
    plot_heat_capacity_gas_data(heat_capacity_data, 'xenon')
    plot_heat_capacity_gas_data(heat_capacity_data, 'neon')
    # Individual plots for each gas
    plot_melting_gas_data(melting_data, 'krypton')
    plot_melting_gas_data(melting_data, 'xenon')
    plot_melting_gas_data(melting_data, 'neon')

    plot_melting_pressure_deviation(melting_data, 'krypton')
    plot_melting_pressure_deviation(melting_data, 'xenon')
    plot_melting_pressure_deviation(melting_data, 'neon')

    # Individual plots for each gas
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'krypton')
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'xenon')
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'neon')

    
    

