# TO CALCULATE EUCLIDEAN DISTANCE BETWEEN TWO DATA SETS AND MAKE A PLOT
# AUTHOR: QI LIU
# DATE: 07/02/2022

# TODO extract Correlation coefficient between two data sets
# TODO 3D PLOT X ENERGY, Y ARRAY RATIO, Z
import nextnanopy as nn
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from scipy.interpolate import Rbf

# edit input file
my_input = nn.InputFile(r'F:\NextNano Data\TLL\test.in')
Barrier_factor = 20  # import from input file
# output file for comparison
C2617_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                         r'Cambridge\Documents\nextnano\Output\test(8)'
                         r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')
# W0938_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
#                          r'Cambridge\Documents\nextnano\Output\W0938full(1)'
#                          r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')
# W0939_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
#                          r'Cambridge\Documents\nextnano\Output\W0939full(4)'
#                          r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')
# W0940_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
#                          r'Cambridge\Documents\nextnano\Output\W0940full'
#                          r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')
# W0941_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
#                          r'Cambridge\Documents\nextnano\Output\W0941full'
#                          r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')

C2617 = C2617_file.variables['1->2'].value
# W0938 = W0938_file.variables['1->2'].value
# W0939 = W0939_file.variables['1->2'].value
# W0940 = W0940_file.variables['1->2'].value
# W0941 = W0941_file.variables['1->2'].value

var1_name = 'alloy1_x'
# my_input.set_variable(var1_name, value=var)
# Remeber to delete duplicate files after this step
# Sweep variable
# dist_array_C2617 = []
# dist_array_W0938 = []
# dist_array_W0939 = []
# dist_array_W0940 = []
# dist_array_W0941 = []

##countor plpt
# my_input.set_variable(var1_name, value=var)
# Remeber to delete duplicate files after this step
# Sweep variable
dist_array = []

limit = np.linspace(0, 1, 51)
trans_stack = np.zeros((limit.size, C2617.size))  # C2617.SIZE MUST EQUAL TO SIZE OF INPUT FILE!!
count = 0
for var in limit:  # remember the 3rd argument needs to +1
    my_input.set_variable(var1_name, value=var)
    var1 = my_input.variables[var1_name].value
    my_input.save(r'F:\NextNano Data\TLL'
                  r'\C2617_f={factor}_{var1}.in'.format(factor=Barrier_factor, var1=var1_name + '=' + str(var1)),
                  automkdir=False,
                  overwrite=True)
    my_input.execute()
    print(f"New variable: {my_input.get_variable('{inputvar}'.format(inputvar=var1_name)).text}")
    output = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                         r'\Documents\nextnano\Output\C2617_f={factor}_{var1}'
                         r'\bias_00000\transmission_cbr_Gamma1.dat'
                         .format(factor=Barrier_factor, var1=var1_name + '=' + str(var1)), product='nextnano++')
    output_trans = output.variables['1->2'].value
    trans_stack[count, :] = output_trans
    count += 1

Energy = np.array(C2617_file.variables['Energy'].value)
levels = np.linspace(0, 1, 21)  # colour bar levels
X, Y = np.meshgrid(Energy, limit)

fig = plt.figure()
img = plt.contourf(X, Y, trans_stack, 20, levels=levels, cmap=cm.viridis)
plt.xlim(0, 1)
plt.xlabel('Energy (eV)')
plt.ylabel('Alloy Ratio')
plt.title('Barrier factor = {fac}'.format(fac=Barrier_factor))
cbar = plt.colorbar(img)
cbar.set_label('Transmission Probability')
plt.tight_layout()
fig.savefig(
    r'E:\OneDrive--University of Cambridge'
    r'\OneDrive - University of Cambridge\Desktop'
    r'\Research\MBE Wafer Growth\Python figs\TransColourMap_f={factor}.png'.format(factor=Barrier_factor))
plt.show()
