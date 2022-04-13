# TODO Check why integral does not agree with result in NextNano?
# TODO try dx as an array
# for quick view of plot and calculation of electron density integral
# Author: Qi Liu
# Date 26/01/2022


import numpy as np
import nextnanopy as nn
import matplotlib.pyplot as plt
import pandas as pd

# import input file
my_input = nn.InputFile(r'F:\NextNano Data\TLL\TLL_Wafer_Trial.in')

# import output datafile
bandedge_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                            r' - University of Cambridge\Documents\nextnano\Output'
                            r'\TLL_Wafer_Trial\bias_00000\bandedges.dat', product='nextnano++')

edensity_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                            r' - University of Cambridge\Documents\nextnano\Output'
                            r'\TLL_Wafer_Trial\bias_00000\density_electron.dat', product='nextnano++')

wavefunction_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                                r' - University of Cambridge\Documents\nextnano\Output\TLL_Wafer_Trial'
                                r'\bias_00000\Quantum\probabilities_shift_cbr_Gamma_00000.dat', product='nextnano++')

x_coords = bandedge_file.coords['x']
# print(bandedge_file.variables)
Energy_Gamma = bandedge_file.variables['Gamma']
Eenergy_Fermi = bandedge_file.variables['electron_Fermi_level']

x_coords2 = edensity_file.coords['x']
edensity = edensity_file.variables['Electron_density']

# integral interval







# z = function_integral(x, y)
# print(z)
# print(dx)



# print("total charge is " + str(total_charge))

# QW1_Charge = np.sum(edensity.value[indices_QW1] * dx)
# QW2_Charge = np.sum(edensity.value[indices_QW2] * dx)
# print(QW1_Charge)
# print(QW2_Charge)

# test: export edensity to a csv file
# pd.DataFrame(edensity.value).to_csv(r"E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge\Desktop\Research\NextNano\file.csv")


fig, (ax1, ax2) = plt.subplots(1, 2)
axa = ax1.twinx()  # instantiate a second axes that shares the same x-axis
axb = ax2.twinx()  # instantiate a second axes that shares the same x-axis

ax1.plot(x_coords.value, Energy_Gamma.value, label=Energy_Gamma.name, linewidth=0.5)
ax1.plot(x_coords.value, Eenergy_Fermi.value, label=Eenergy_Fermi.name, linewidth=0.8)
ax1.set_xlabel(x_coords.label)
ax1.set_ylabel(Energy_Gamma.label)
# ax1.legend()
# ax1.set_xlim(10,200)
# ax1.set_ylim(-.1,.6)

color = 'tab:green'
axa.set_ylabel(edensity.label, color=color)  # we already handled the x-label with ax1
axa.plot(x_coords2.value, edensity.value, label=edensity.name, color=color, linewidth=0.8)
axa.tick_params(axis='y', labelcolor=color)
# ax2.legend()
# plt.title("Qi's nextnano simulated results")

ax2.plot(x_coords.value, Energy_Gamma.value, label=Energy_Gamma.name, linewidth=0.5)
ax2.plot(x_coords.value, Eenergy_Fermi.value, label=Eenergy_Fermi.name, linewidth=0.8)
# ax1.plot(x_coords2.value,edensity.value, label=edensity.name)
ax2.set_xlabel(x_coords.label)
ax2.set_ylabel(Energy_Gamma.label)
# ax2.legend()
ax2.set_xlim(10, 200)
ax2.set_ylim(-.1, .6)

color = 'tab:cyan'
axb.set_ylabel(edensity.label, color=color)  # we already handled the x-label with ax1
axb.plot(x_coords2.value, edensity.value, label=edensity.name, color=color, linewidth=0.8)
axb.tick_params(axis='y', labelcolor=color)
axb.set_ylim(0, 1)  # unit: 1e18 cm^-3
# ax2.legend()
# plt.title("Qi's nextnano simulated results")
plt.tight_layout()

# fig.savefig('E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge\Desktop\Research\Wafer Growth\Qi_simulation2.png') #save figure to a given file route
plt.show()
