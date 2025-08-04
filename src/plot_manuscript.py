import numpy as np
import matplotlib.pyplot as plt
import os
from GENG5511_Thesis.src.p_functions import psub,pmelt
# from GENG5511_Thesis.src.constants import *
# from constants import *

vle_path = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\VLE.txt"
T_VLE, p_VLE = np.loadtxt(vle_path, skiprows=2, unpack=True)


# Constants
# Plot critical point
T_crit = 150.69
p_crit = 4.863
# Plot triple point
T_triple1 = 83.806
p_triple1 = 0.068891
fontsize = 14
markersize = 8
linewidth = 0.9

# Define RGB color as normalized
mycolor = np.array([0, 0, 0]) / 256

# Placeholder functions for pmelt and psub (they must be defined)

# Generate melting curve
T_melting = np.linspace(83.806, 1000.806, 500)
p_melting = np.array([pmelt(T) for T in T_melting])

# Generate sublimation curve
T_sublimation = np.linspace(0.806, 83.806, 500)
p_sublimation = np.array([psub(T) for T in T_sublimation])

# Begin plotting
plt.figure(figsize=(6, 4))

# Plot VLE curve
plt.plot(p_VLE, T_VLE, color=mycolor, linewidth=linewidth, label='VLE')

# Plot melting curve
plt.plot(p_melting, T_melting, color=mycolor,
         linewidth=linewidth, label='Melting')

# Plot sublimation curve
plt.plot(p_sublimation, T_sublimation, color=mycolor,
         linewidth=linewidth, label='Sublimation')


plt.plot(p_crit, T_crit, 's', color=mycolor,
         markersize=markersize, label='Critical Point')


plt.plot(p_triple1, T_triple1, 'd', color='red',
         markersize=9, label='Triple Point')

# Add text
plt.text(1e-3, 38, 'Solid', fontsize=14)
plt.text(1e-7, 125, 'Vapour', fontsize=14)
plt.text(2, 115, 'Liquid', fontsize=14)

# Labels and formatting
plt.xscale('log')
plt.xlim([1e-10, 1e3])
plt.ylim([0, 250])
plt.xlabel('p / MPa', fontsize=fontsize)
plt.ylabel('T / K', fontsize=fontsize)
plt.tick_params(labelsize=fontsize)
plt.grid(True)
plt.legend(fontsize=10)

plt.show()
