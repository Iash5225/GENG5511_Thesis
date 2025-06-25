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

    # output_filepath = f"{OUTPUT_FILEPATH}\{gas_name}_melting_temperatures_plot.png"
    # plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    plt.show()
