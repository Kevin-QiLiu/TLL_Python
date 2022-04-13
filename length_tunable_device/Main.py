# TODO Check why integral does not agree with result in NextNano?
# TODO try dx as an array
# for quick view of plot and calculation of electron density integral
# Author: Qi Liu
# Date 26/01/2022


import numpy as np
import nextnanopy as nn
import matplotlib.pyplot as plt
import pandas as pd


class EDIT_INPUT:

    def __init__(self, input_path, var=None):
        self.path = input_path
        self.file = nn.InputFile(input_path)
        if var is not None:
            self.var = self.file.variables[var].value


class READ_DATA:

    def __init__(self, data_path, var):
        self.path = data_path
        self.file = nn.DataFile(data_path, product='nextnano++')
        self.x = self.file.coords['x'].value
        self.y = self.file.coords['y'].value
        self.z = self.file.coords['z'].value
        self.var = self.file[var].value

    def slice_var(self, vert_z):  # slice 3d array at given z
        idx = (np.abs(self.z - vert_z)).argmin()  # find index in array z
        var_transposed = np.transpose(self.var, (2, 1, 0))  # transpose 3d array so array[z] gives xy plane at given z
        var_sliced = var_transposed[idx]
        return var_transposed, var_sliced, idx


path_input = r'F:\NextNano Data\Length_tunable\trial.in'
input_file = EDIT_INPUT(path_input).file

######################## structual parameters #######################
surface_cap_width = EDIT_INPUT(path_input, 'surface_cap_width').var
upper_AlGaAs_doping_layer_width = EDIT_INPUT(path_input, 'upper_AlGaAs_doping_layer_width').var
upper_AlGaAs_modulation_spacer_width = EDIT_INPUT(path_input, 'upper_AlGaAs_modulation_spacer_width').var
# qw1 parameters
QW1_width = EDIT_INPUT(path_input, 'QW1_width').var
# Starting position of QW1
QW1_start = surface_cap_width + upper_AlGaAs_doping_layer_width + upper_AlGaAs_modulation_spacer_width
QW1_end = QW1_start + QW1_width
Barrier_Width_1 = EDIT_INPUT(path_input, 'Barrier_Width_1').var
Barrier_Width_2 = EDIT_INPUT(path_input, 'Barrier_Width_2').var
Barrier_Period = Barrier_Width_1 + Barrier_Width_2
superlattice_repeat_factor = EDIT_INPUT(path_input, 'superlattice_repeat_factor').var

# Starting position of QW2
QW2_start = QW1_end + Barrier_Period * superlattice_repeat_factor
QW2_width = EDIT_INPUT(path_input, 'QW2_width').var
QW2_end = QW2_start + QW2_width
####################################################################


bandedges_path = r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge\Documents\nextnano\Output\trial(4)\bias_00000\bandedges.fld'
datafile_3d = READ_DATA(bandedges_path, 'Gamma').file
y = READ_DATA(bandedges_path, 'Gamma').y
x = READ_DATA(bandedges_path, 'Gamma').x
z = READ_DATA(bandedges_path, 'Gamma').z

pos_slice = (QW1_start + QW1_end) / 2  # sliced z position
potential = READ_DATA(bandedges_path, 'Gamma').slice_var(pos_slice)[1]

X, Y = np.meshgrid(x, y)

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection='3d')
ax.view_init(30, 150)
map_3D = ax.plot_surface(X, Y, potential, rstride=1, cstride=1, cmap='binary', edgecolor='none')
# pcolor = ax.pcolormesh(x, y, potential, cmap='hot')
cbar = fig.colorbar(map_3D)
cbar.set_label("potential(eV)")

ax.set_xlabel("x(nm)")
ax.set_ylabel("y(nm)")
ax.axis('off')
fig.tight_layout()
plt.show()
