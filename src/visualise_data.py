import pandas as pd
import matplotlib.pyplot as plt
from pint import UnitRegistry
from constants import *
import numpy as np
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
    krypton_data = melting_data[melting_data['gas'] == 'krypton']
    xenon_data = melting_data[melting_data['gas'] == 'xenon']
    neon_data = melting_data[melting_data['gas'] == 'neon']

    krypton_data_thermal_coefficient = thermalcoefficient_data[
        thermalcoefficient_data['gas'] == 'krypton']
    xenon_data_thermal_coefficient = thermalcoefficient_data[thermalcoefficient_data['gas'] == 'xenon']
    neon_data_thermal_coefficient = thermalcoefficient_data[thermalcoefficient_data['gas'] == 'neon']

    # Individual plots for each gas
    plot_melting_gas_data(melting_data, 'krypton')
    plot_melting_gas_data(melting_data, 'xenon')
    plot_melting_gas_data(melting_data, 'neon')

    # Individual plots for each gas
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'krypton')
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'xenon')
    plot_thermal_coefficient_gas_data(thermalcoefficient_data, 'neon')

