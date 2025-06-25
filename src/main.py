import pandas as pd
import matplotlib.pyplot as plt
from pint import UnitRegistry
from constants import *
import numpy as np
from scipy.optimize import curve_fit
ureg = UnitRegistry()


def melting_pressure_equation(T, e_4, e_5, e_6, e_7, T_t, P_t):
    """
    Equation for melting pressure with safe power handling.
    """
    term1 = np.where(T / T_t - 1 >= 0, (T / T_t - 1) ** e_5, 0)
    term2 = np.where(T / T_t - 1 >= 0, (T / T_t - 1) ** e_7, 0)

    return P_t * (1 + e_4 * term1 + e_6 * term2)

def plot_melting_gas_data(data, gas_name):
    grouped = data.groupby(['Year', 'Author'])

    # Define markers and colors

    plt.figure(figsize=(10, 6))

    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            group['Temperature'], group['Pressure'],
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
            data['Temperature'], KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t)
        X_MIN = MELTING_KRYPTON_X_MIN
        X_MAX = MELTING_KRYPTON_X_MAX
        Y_MIN = MELTING_KRYPTON_Y_MIN
        Y_MAX = MELTING_KRYPTON_Y_MAX
        # plt.plot(np.sort(gas_data['Temperature_Kelvin']),
        #          np.sort(P_melt), color='black')
    elif gas_name == "xenon":
        P_melt = melting_pressure_equation(
            data['Temperature'], XENON_E_4, XENON_E_5, XENON_E_6, XENON_E_7, XENON_T_t, XENON_P_t)
        X_MIN = MELTING_XENON_X_MIN
        X_MAX = MELTING_XENON_X_MAX
        Y_MIN = MELTING_XENON_Y_MIN
        Y_MAX = MELTING_XENON_Y_MAX
        # plt.plot(np.sort(gas_data['Temperature_Kelvin']),
        #          np.sort(P_melt), color='black')
    else:
        P_melt = melting_pressure_equation(
            data['Temperature'], NEON_E_4, NEON_E_5, NEON_E_6, NEON_E_7, NEON_T_t, NEON_P_t)
        X_MIN = MELTING_NEON_X_MIN
        X_MAX = MELTING_NEON_X_MAX
        Y_MIN = MELTING_NEON_Y_MIN
        Y_MAX = MELTING_NEON_Y_MAX
        # plt.plot(np.sort(gas_data['Temperature_Kelvin']),
        #          np.sort(P_melt), color='black')

    # plt.scatter(gas_data['Temperature_Kelvin'], P_melt, color='black',
    #          linestyle='-', linewidth=2, label='Calculated Melting Pressure')
    plt.plot(np.sort(data['Temperature']),
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
    # plt.grid(True, which="both", linestyle='--', alpha=0.3)

    # Place the legend outside the plot
    plt.legend(
        loc='upper left',
        bbox_to_anchor=(1.05, 1),  # Position the legend outside the plot
        fontsize=8,
        ncol=1  # Adjust the number of columns in the legend
    )
    # plt.grid(True, linestyle='--', alpha=0.6)
    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    plt.title(f'Melting Pressure for {gas_name}', fontsize=14)

    # output_filepath = f"{OUTPUT_FILEPATH}\{gas_name}_melting_temperatures_plot.png"
    # plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()
    

# Main Execution
if __name__ == "__main__":

    # Read the Excel file into a DataFrame
    # XENON
    xe_melting_data = pd.read_excel(XE_DATA_FILEPATH, sheet_name=MELTING_SHEETNAME, header=2)
    xe_melting_data = xe_melting_data.filter(
        items=['Year', 'Author', 'Temperature', 'Pressure'])
    xe_melting_data = xe_melting_data.drop([0, 1], axis=0).reset_index(drop=True)
    plot_melting_gas_data(xe_melting_data, 'xenon')
    

