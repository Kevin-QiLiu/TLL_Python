# TO CALCULATE EUCLIDEAN DISTANCE BETWEEN TWO DATA SETS AND MAKE A PLOT
# AUTHOR: QI LIU
# DATE: 07/02/2022

# TODO extract Correlation coefficient between two data sets
# TODO 3D PLOT X ENERGY, Y ARRAY RATIO, Z
import nextnanopy as nn
import numpy as np
import time
import matplotlib.pyplot as plt


def euclidean_dist(array1, array2):
    dummy_array = []
    if array1.size != array2.size:
        return False
    for i in range(len(array1)):
        dist_sqr = (array1[i] - array2[i]) ** 2
        dummy_array = np.append(dummy_array, dist_sqr)
    dist = np.sqrt(np.sum(dummy_array))
    return dist


# edit input file
my_input = nn.InputFile(r'F:\NextNano Data\TLL\C2617_adoptedyiqing.in')
Barrier_factor = 10  # import from input file

C2617_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                         r'Cambridge\Documents\nextnano\Output\C2617_adoptedyiqing(30)'
                         r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')
W0938_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                         r'Cambridge\Documents\nextnano\Output\W0938full(1)'
                         r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')
W0939_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                         r'Cambridge\Documents\nextnano\Output\W0939full(4)'
                         r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')
W0940_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                         r'Cambridge\Documents\nextnano\Output\W0940full'
                         r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')
W0941_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                         r'Cambridge\Documents\nextnano\Output\W0941full'
                         r'\bias_00000\transmission_cbr_Gamma1.dat', product='nextnano++')

C2617 = C2617_file.variables['1->2'].value
W0938 = W0938_file.variables['1->2'].value
W0939 = W0939_file.variables['1->2'].value
W0940 = W0940_file.variables['1->2'].value
W0941 = W0941_file.variables['1->2'].value

var1_name = 'alloy1_x'
# my_input.set_variable(var1_name, value=var)
# Remeber to delete duplicate files after this step
# Sweep variable
dist_array_C2617 = []
dist_array_W0938 = []
dist_array_W0939 = []
dist_array_W0940 = []
dist_array_W0941 = []

limit = np.linspace(0, 1, 51)
for var in limit:  # remember the 3rd argument needs to +1
    my_input.set_variable(var1_name, value=var)
    var1 = my_input.variables[var1_name].value
    my_input.save(r'F:\NextNano Data\TLL'
                  r'\C2617_f={factor}_{var1}.in'.format(factor=Barrier_factor, var1=var1_name + '=' + str(var1)),
                  automkdir=False,
                  overwrite=True)
    my_input.execute()
    time.sleep(0.5)
    print(f"New variable: {my_input.get_variable('{inputvar}'.format(inputvar=var1_name)).text}")
    output = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                         r'\Documents\nextnano\Output\C2617_f={factor}_{var1}'
                         r'\bias_00000\transmission_cbr_Gamma1.dat'
                         .format(factor=Barrier_factor, var1=var1_name + '=' + str(var1)), product='nextnano++')
    output_trans = output.variables['1->2'].value
    dist_C2617 = euclidean_dist(C2617, output_trans)
    dist_W0938 = euclidean_dist(W0938, output_trans)
    dist_W0939 = euclidean_dist(W0939, output_trans)
    dist_W0940 = euclidean_dist(W0940, output_trans)
    dist_W0941 = euclidean_dist(W0941, output_trans)
    dist_array_C2617 = np.append(dist_array_C2617, dist_C2617)
    dist_array_W0938 = np.append(dist_array_W0938, dist_W0938)
    dist_array_W0939 = np.append(dist_array_W0939, dist_W0939)
    dist_array_W0940 = np.append(dist_array_W0940, dist_W0940)
    dist_array_W0941 = np.append(dist_array_W0941, dist_W0941)

fig, axs = plt.subplots(2, 3)
axs[0, 0].plot(limit, dist_array_C2617)
axs[0, 0].set_title(r'R to C2617')
axs[0, 1].plot(limit, dist_array_W0940)
axs[0, 1].set_title(r'R to W0940')
axs[0, 2].plot(limit, dist_array_W0941)
axs[0, 2].set_title(r'R to W0941')
axs[1, 0].plot(limit, dist_array_W0938)
axs[1, 0].set_title(r'R to W0938')
axs[1, 1].plot(limit, dist_array_W0939)
axs[1, 1].set_title(r'R to W0939')

fig.text(0.55, 0.0015, 'Alloy Ratio', ha='center')
fig.text(0.0015, 0.5, 'Euclidean Dist R', va='center', rotation='vertical')
fig.suptitle('Barrier factor = {fac}'.format(fac=Barrier_factor))
# plt.title('Barrier factor = 20')
plt.tight_layout()
fig.savefig(r'E:\OneDrive--University of Cambridge'
            r'\OneDrive - University of Cambridge\Desktop\Research\MBE Wafer Growth'
            r'\Python figs\RVsAlloy_f={factor}2.png'.format(factor=Barrier_factor))
plt.show()

##countor plpt

# print(C2617_file.coords)
# C2617 = C2617_file.variables['1->2'].value
# x = C2617_file.coords['Energy']
# var1_name = 'alloy1_x'
# # my_input.set_variable(var1_name, value=var)
# # Remeber to delete duplicate files after this step
# # Sweep variable
# dist_array = []
#
# limit = np.linspace(0.46, 0.47, 6)
# trans_stack = np.zeros((limit.size, C2617.size))
# count = 0
# for var in limit:  # remember the 3rd argument needs to +1
#     my_input.set_variable(var1_name, value=var)
#     var1 = my_input.variables[var1_name].value
#     my_input.save(r'F:\NextNano Data\TLL'
#                   r'\C2617_adoptedyiqing_{var1}.in'.format(var1=var1_name + '=' + str(var1)), automkdir=False,
#                   overwrite=True)
#     my_input.execute()
#     time.sleep(0.5)
#     print(f"New variable: {my_input.get_variable('{inputvar}'.format(inputvar=var1_name)).text}")
#     output = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
#                          r'\Documents\nextnano\Output\C2617_adoptedyiqing_{var1}'
#                          r'\bias_00000\transmission_cbr_Gamma1.dat'
#                          .format(var1=var1_name + '=' + str(var1)), product='nextnano++')
#     output_trans = output.variables['1->2'].value
#     trans_stack[count, :] = output_trans
#     dist0 = euclidean_dist(C2617, output_trans)
#     dist_array = np.append(dist_array, dist0)
#     count += 1
#
# print(C2617.size, x.size, trans_stack.size)
# plt.plot(limit, dist_array)
# plt.xlabel('Alloy ratio')
# plt.ylabel('R')
# plt.show()
