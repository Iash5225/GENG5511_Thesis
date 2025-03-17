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

def melting_auxillary_function(T):
    """
    Compute the melting pressure using the given equation.
    """
    return P_t * (1 + E_4 * (T / T_t - 1) ** E_5 + E_6 * (T / T_t - 1) ** E_7)


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
def plot_melting_gas_data(data, gas_name):
    gas_data = data[data['gas'] == gas_name]
    grouped = gas_data.groupby(['year', 'author'])

    # Define markers and colors

    plt.figure(figsize=(10, 6))

    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            group['Temperature_Kelvin'], group['Pressure_Mpa'],
            s=100,  # Smaller marker size
            label=f"{year}, {author}",
            # Assigning custom color
            edgecolor=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            # Assigning custom marker
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
            # marker=markers[i % len(markers)],  # Unique marker style
            alpha=0.8,  # Slight transparency for overlapping markers
            facecolors='none',  # Open symbols
        )
    P_melt = melting_auxillary_function(gas_data['Temperature_Kelvin'])
    
    # plt.scatter(gas_data['Temperature_Kelvin'], P_melt, color='black',
    #          linestyle='-', linewidth=2, label='Calculated Melting Pressure')
    plt.plot(np.sort(gas_data['Temperature_Kelvin']),
             np.sort(P_melt), color='black')

    # Set labels with italicized text to match LaTeX style in the image
    plt.xlabel(r'$T \,/\, K$', fontsize=14, fontstyle='italic')
    plt.ylabel(r'$p \,/\, MPa$', fontsize=14, fontstyle='italic')

    # Set a logarithmic scale for the y-axis
    plt.yscale('log')

    # Set x-axis limits similar to the uploaded image
    plt.xlim(0, 350)
    plt.ylim(1, 10**4)

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

    output_filepath = f"{OUTPUT_FILEPATH}{gas_name}_melting_temperatures_plot.png"
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()


def melting_pressure_equation(T, e_4, e_5, e_6, e_7):
    """
    Equation for melting pressure with safe power handling.
    """
    term1 = np.where(T / T_t - 1 >= 0, (T / T_t - 1) ** e_5, 0)
    term2 = np.where(T / T_t - 1 >= 0, (T / T_t - 1) ** e_7, 0)

    return P_t * (1 + e_4 * term1 + e_6 * term2)



def fit_melting_pressure(data):
    """
    Fit the melting pressure equation to experimental data to find e_4, e_5, e_6, and e_7.
    Returns a dictionary with gas names as keys and optimized parameters as values.
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

        # print(T_data)
        # print(P_data)

        # Ensure enough valid data points
        if len(T_data) < 4:
            print(f"Not enough valid data to fit for {gas}. Skipping.")
            gas_coefficients[gas] = None
            continue

        # Initial guesses for e_4, e_5, e_6, e_7
        initial_guess = [1500, 1.7, 4600, 0.98]

        try:
            popt, _ = curve_fit(melting_pressure_equation,
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
            s=100,  # Marker size
            # Assigning custom color
            edgecolor=CUSTOMCOLORS[i % len(CUSTOMCOLORS)],
            # Assigning custom marker
            marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
            # marker=markers[i % len(markers)],  # Marker style
            facecolors='none',  # Open symbols
        )

    plt.xlabel('Temperature (Kelvin)', fontsize=12)
    plt.ylabel('Alpha (K^-1)', fontsize=12)
    plt.title(f'Thermal Coefficient for {gas_name}', fontsize=14)

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

    output_filepath = f"{OUTPUT_FILEPATH}{gas_name}_thermal_coefficient_plot.png"
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()


# Main Execution
if __name__ == "__main__":
    

    # Read the Excel file into a DataFrame
    melting_data = pd.read_excel(FILEPATH, sheet_name=MELTING_SHEETNAME, header=1)
    thermalcoefficient_data = pd.read_excel(
        FILEPATH, sheet_name=THERMAL_COEFFICIENT_SHEETNAME, header=1)
    melting_data = melting_data.drop(melting_data.columns[0], axis=1)
    # data = data.dropna()  # Drop rows with missing values
    melting_data['Temperature_Kelvin'] = pd.to_numeric(
        melting_data['Temperature_Kelvin'], errors='coerce')
    melting_data['Pressure_atm'] = pd.to_numeric(
        melting_data['Pressure_atm'], errors='coerce')
    thermalcoefficient_data = thermalcoefficient_data.drop(
        thermalcoefficient_data.columns[0], axis=1)

    # Apply the conversion function
    # Convert Pressure_atm to Pressure_MPa
    # Convert Pressures
    melting_data = convert_pressure(melting_data)

    # Filter data by gas type
    # krypton_data = melting_data[melting_data['gas'] == 'krypton']
    # xenon_data = melting_data[melting_data['gas'] == 'xenon']
    # neon_data = melting_data[melting_data['gas'] == 'neon']

    # krypton_data_thermal_coefficient = thermalcoefficient_data[
    #     thermalcoefficient_data['gas'] == 'krypton']
    # xenon_data_thermal_coefficient = thermalcoefficient_data[thermalcoefficient_data['gas'] == 'xenon']
    # neon_data_thermal_coefficient = thermalcoefficient_data[thermalcoefficient_data['gas'] == 'neon']

    # Fit constants to the data per gas
    gas_coefficients = fit_melting_pressure(melting_data)
    for gas, params in gas_coefficients.items():
        if params is not None:
            print(
                f"Optimized Constants for {gas}: e_4={params[0]:.2f}, e_5={params[1]:.3f}, e_6={params[2]:.2f}, e_7={params[3]:.5f}")

    # Individual plots for each gas
    plot_melting_gas_data(melting_data, 'krypton')
    plot_melting_gas_data(melting_data, 'xenon')
    plot_melting_gas_data(melting_data, 'neon')

    # Fit constants to the data


    # Individual plots for each gas
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'krypton')
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'xenon')
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'neon')

