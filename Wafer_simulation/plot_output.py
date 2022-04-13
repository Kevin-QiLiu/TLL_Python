# Plot NextNano output files
# Author: Qi Liu Date: 25/01/2022


import nextnanopy as nn
import matplotlib.pyplot as plt


def get_choice(choices):  # from github
    choice = ""
    while choice not in choices:
        choice = input("Choose one of [%s]:" % ", ".join(choices))
    return choice


def output_file(choice):
    isquantum = 1
    if not isquantum:
        file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                           r' - University of Cambridge\Documents\nextnano\Output'
                           r'\TLL_Wafer_Trial\bias_00000\{choice}.dat'.format(choice=str(choice)),
                           product='nextnano++')  # route of file for non-quantum simulation
    else:
        file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                           r' - University of Cambridge\Documents\nextnano\Output'
                           r'\TLL_Wafer_Trial\bias_00000\Quantum\{choice}.dat'.format(choice=str(choice)),
                           product='nextnano++')  # route of file for quantum simulation, e.g., prob amplitude
    return file


def output_variable(file):
    print(file.variables)
    variable = input("choose one key from above: ")
    return variable


# input parameters
choice1 = get_choice(["bandedges", "density_electron", "amplitudes_shift_cbr_Gamma_00000", "transmission_cbr_Gamma1",
                      "dummy1:variable1,cbr=Q"])
choice1 = str(choice1)
file = output_file(choice1)
variable = output_variable(file)
x_coords = file.coords['x']
variable = file.variables[variable]
xlim_0 = 10
xlim_1 = 200
ylim_0 = -.1
ylim_1 = .6
isVariable2 = 1
isMultVariable = 0

if isVariable2:
    choice2 = get_choice(
        ["bandedges", "density_electron", "amplitudes_shift_cbr_Gamma_00000", "dummy2:variable2,cbr=Q"])
    choice2 = str(choice2)
    file2 = output_file(choice2)
    variable2 = output_variable(file2)
    x_coords2 = file2.coords['x']
    variable2 = file2.variables[variable2]
    isE_fermi = 1
    fig, ax1 = plt.subplots(1)
    ax1.plot(x_coords.value, variable.value, label=variable.name, linewidth=0.5)
    ax1.set_xlabel(x_coords.label)
    ax1.set_ylabel(variable.label)
    ax1.set_xlim(10, 200)
    ax1.set_ylim(-.1, .6)
    ax2 = ax1.twinx()
    color = 'tab:green'
    ax2.set_ylabel(variable2.label, color=color)
    ax2.plot(x_coords2.value, variable2.value, label=variable2.name, color=color, linewidth=0.8)
    ax2.tick_params(axis='y', labelcolor=color)
    if isMultVariable:  # plot multivariables (=6 in this case) often used in plotting amplitudes, wavefunctions, etc.
        var3 = output_variable(file2)
        var4 = output_variable(file2)
        var5 = output_variable(file2)
        var6 = output_variable(file2)
        var3 = file2.variables[var3]
        var4 = file2.variables[var4]
        var5 = file2.variables[var5]
        var6 = file2.variables[var6]
        ax2.plot(x_coords2.value, var3.value, label=var3.name, color='tab:red', linewidth=0.8)
        ax2.plot(x_coords2.value, var4.value, label=var4.name, color='tab:cyan', linewidth=0.8)
        ax2.plot(x_coords2.value, var5.value, label=var5.name, color='tab:blue', linewidth=0.8)
        ax2.plot(x_coords2.value, var6.value, label=var6.name, color=color, linewidth=0.8)
    elif isE_fermi == "True":  # include fermi energy
        bandedge_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                                    r' - University of Cambridge\Documents\nextnano\Output'
                                    r'\TLL_Wafer_Trial\bias_00000\bandedges.dat', product='nextnano++')
        x_coords_fermi = bandedge_file.coords['x']
        Eenergy_Fermi = bandedge_file.variables['electron_Fermi_level']  # define plotting fermi energy
        # define plotting variable
        ax1.plot(x_coords_fermi.value, Eenergy_Fermi.value, label=Eenergy_Fermi.name, linewidth=0.8)
else:  # not include second variable on same graph
    isE_fermi = 1
    if isE_fermi:
        bandedge_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                                    r' - University of Cambridge\Documents\nextnano\Output'
                                    r'\TLL_Wafer_Trial\bias_00000\bandedges.dat', product='nextnano++')
        x_coords_fermi = bandedge_file.coords['x']
        Eenergy_Fermi = bandedge_file.variables['electron_Fermi_level']  # define plotting fermi energy
        # define plotting variable
        x_coords = file.coords['x']
        variable = file.variables[variable]
        fig, (ax1) = plt.subplots(1)
        ax1.plot(x_coords_fermi.value, Eenergy_Fermi.value, label=Eenergy_Fermi.name, linewidth=0.8)
        ax1.plot(x_coords.value, variable.value, label=variable.name, linewidth=0.5)
        ax1.set_xlabel(x_coords.label)
        ax1.set_ylabel(variable.label)
    else:
        x_coords = file.coords['x']
        variable = file.variables[variable]
        fig, (ax1) = plt.subplots(1)
        ax1.plot(x_coords.value, variable.value, label=variable.name, linewidth=0.5)
        ax1.set_xlabel(x_coords.label)
        ax1.set_ylabel(variable.label)

plt.tight_layout()
plt.show()

