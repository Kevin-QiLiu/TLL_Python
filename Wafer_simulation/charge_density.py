# TODO Check why integral does not agree with result in NextNano?
# TODO Sweep conc_n and see the integral, electron freezed to avoid parallel conducting??
# TODO Write in OOP Style
# TODO try different ionization energy, what is dx center???? freeze problems in donor region????
# TODO start qw barrier simulation, Try contour plot
# for calculation of charge (electron, ionized donors or acceptors) density integral
# Author: Qi Liu
# Date 26/01/2022


import numpy as np
import nextnanopy as nn
import matplotlib.pyplot as plt
import pandas as pd


def find_nearest(array, value):
    # from stackoverflow
    # return element in array that is nearest to input value and also index of it
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def get_dist_array(array):
    # find distance of element in an array
    # used to calculate differential element dx
    dist_array = []
    for i in range(len(array) - 1):  # -1 for end-point effect, note also len(dist_array) = len(array) - 1
        dist = array[i + 1] - array[i]
        dist_array = np.append(dist_array, dist)
    return dist_array


def dx(coords, factor):
    # to get nonuniform differential element
    # factor convert from nm to whatever
    return (get_dist_array(coords.value)) * factor


def function_integral_global(array1, array2):
    # traverse over the whole array
    # think this as \int f(x) dx, but now dx is not uniform. Integrate over whole system
    # in out case, len(dx)=len(density)-1 so set dx as array2, edensity as array1
    dummy_array = []
    for i in range(len(array2)):
        dummy_array = np.append(dummy_array, (array1[i] * array2[i]))
    integral_global = np.sum(dummy_array)
    return integral_global


def function_integral(array1, array2, index):
    # traverse over specified indices
    # think this as \int f(x) dx, but now dx is not uniform. Integrating limit given in index
    # in out case, len(dx)=len(edensity)-1 so set dx as array2, edensity as array1
    dummy_array = []
    for i in index:
        dummy_array = np.append(dummy_array, (array1[i] * array2[i]))
    integral = np.sum(dummy_array)
    return integral


def total_charge_number(density, diff_x):
    # density and diff_x should be arrays
    return function_integral_global(density, diff_x)


def charge_number(density, diff_x, coords, start, end):
    index_0 = find_nearest(coords.value, start)[1]
    index_1 = find_nearest(coords.value, end)[1]
    indices = [i for i in range(index_0, index_1 + 1)]  # +1 for end-point effect of range()
    return function_integral(density, diff_x, indices)


# import input file
my_input = nn.InputFile(r'F:\NextNano Data\TLL\TLL_Wafer_Trial.in')

# import output electron density datafile
# edensity_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
#                             r' - University of Cambridge\Documents\nextnano\Output'
#                             r'\TLL_Wafer_Trial\bias_00000\density_electron.dat', product='nextnano++')
edensity_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                            r' - University of Cambridge\Documents\nextnano\Output'
                            r'\C2617_adoptedyiqing(9)\bias_00000\density_electron.dat', product='nextnano++')

donor_density_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                                 r' - University of Cambridge\Documents\nextnano\Output'
                                 r'\C2617_adoptedyiqing(9)\bias_00000\density_donor_ionized.dat', product='nextnano++')

acceptor_density_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                                    r' - University of Cambridge\Documents\nextnano\Output'
                                    r'\C2617_adoptedyiqing(9)\bias_00000\density_acceptor_ionized.dat',
                                    product='nextnano++')

# define constants and parameters
convert_factor = 1e-7  # convert nm to cm
convert_factor_nextnano = 1e18  # density has a unit of 1e18 cm^-3 in nextnano as default

# extract output file coordinates
x_coords = edensity_file.coords['x']
edensity = edensity_file.variables['Electron_density'].value * convert_factor_nextnano
x_coords2 = donor_density_file.coords['x']
donor_density = donor_density_file['Ionized_donor_density'].value * convert_factor_nextnano
x_coords3 = acceptor_density_file.coords['x']
acceptor_density = acceptor_density_file['Ionized_acceptor_density'].value * convert_factor_nextnano

######################## structual parameters #######################
surface_cap_width = my_input.variables['surface_cap_width'].value
upper_AlGaAs_doping_layer_width = my_input.variables['upper_AlGaAs_doping_layer_width'].value
upper_AlGaAs_modulation_spacer_width = my_input.variables['upper_AlGaAs_modulation_spacer_width'].value
# qw1 parameters
QW1_width = my_input.variables['QW1_width'].value
# Starting position of QW1
QW1_start = surface_cap_width + upper_AlGaAs_doping_layer_width + upper_AlGaAs_modulation_spacer_width
QW1_end = QW1_start + QW1_width
Barrier_Width_1 = my_input.variables['Barrier_Width_1'].value
Barrier_Width_2 = my_input.variables['Barrier_Width_2'].value
Barrier_Period = Barrier_Width_1 + Barrier_Width_2
superlattice_repeat_factor = my_input.variables['superlattice_repeat_factor'].value

# Starting position of QW2
QW2_start = QW1_end + Barrier_Period * superlattice_repeat_factor
QW2_width = my_input.variables['QW2_width'].value
QW2_end = QW2_start + QW2_width
####################################################################

# set up integrals
QWs = 0  # = 1 if want number in Quantum Well, otherwise can define arbitrary start and end point
diffx = dx(x_coords, convert_factor)
diffx2 = dx(x_coords2, convert_factor)
diffx3 = dx(x_coords3, convert_factor)
region1_start = 10
region1_end = 70
region2_start = 160
region2_end = 220
# for region in QWs:
if QWs:
    region1_start = QW1_start
    region1_end = QW1_end
    region2_start = QW2_start
    region2_end = QW2_end
    print("region1/2 = QW1/2")

# integrating
total_electron = total_charge_number(edensity, diffx)
Region1_electron = charge_number(edensity, diffx, x_coords, region1_start, region1_end)
Region2_electron = charge_number(edensity, diffx, x_coords, region2_start, region2_end)
total_donor = total_charge_number(donor_density, diffx2)
Region1_donor = charge_number(donor_density, diffx2, x_coords2, region1_start, region1_end)
Region2_donor = charge_number(donor_density, diffx2, x_coords2, region2_start, region2_end)
total_acceptor = total_charge_number(acceptor_density, diffx3)
Region1_acceptor = charge_number(acceptor_density, diffx3, x_coords3, region1_start, region1_end)
Region2_acceptor = charge_number(acceptor_density, diffx3, x_coords3, region2_start, region2_end)

print("{Type} = {:.3e}".format(total_electron, Type='total_electron') + " cm^-2")
print("{Type} = {:.3e}".format(total_donor, Type='total_donor') + " cm^-2")
print("{Type} = {:.3e}".format(total_acceptor, Type='total_acceptor') + " cm^-2")
print("{Type} = {:.3e}".format(Region1_electron, Type='Region1_electron') + " cm^-2"
      + ", region1 = [{x0},{x1}]".format(x0=region1_start, x1=region1_end))
print("{Type} = {:.3e}".format(Region2_electron, Type='Region2_electron') + " cm^-2"
      + ", region2 = [{x0},{x1}]".format(x0=region2_start, x1=region2_end))
print("{Type} = {:.3e}".format(Region1_donor, Type='Region1_ionized_donor') + " cm^-2"
      + ", region1 = [{x0},{x1}]".format(x0=region1_start, x1=region1_end))
print("{Type} = {:.3e}".format(Region2_donor, Type='Region2_ionized_donor') + " cm^-2"
      + ", region2 = [{x0},{x1}]".format(x0=region2_start, x1=region2_end))
print("{Type} = {:.3e}".format(Region1_acceptor, Type='Region1_ionized_acceptor') + " cm^-2"
      + ", region1 = [{x0},{x1}]".format(x0=region1_start, x1=region1_end))
print("{Type} = {:.3e}".format(Region2_acceptor, Type='Region2_ionized_acceptor') + " cm^-2"
      + ", region2 = [{x0},{x1}]".format(x0=region2_start, x1=region2_end))
# print(edensity_file.variables['Electron_density'].value)


# test: export edensity to a csv file
# pd.DataFrame(edensity.value).to_csv(r"E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge\Desktop\Research\NextNano\file.csv")
