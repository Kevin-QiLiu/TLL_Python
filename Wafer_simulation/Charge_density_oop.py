# TODO Check why integral does not agree with result in NextNano?
# TODO Sweep conc_n and see the integral, electron freezed to avoid parallel conducting??
# TODO try different ionization energy, what is dx center???? freeze problems in donor region????
# TODO start qw barrier simulation, Try contour plot
# for calculation of charge (electron, ionized donors or acceptors) density integral
# Author: Qi Liu
# Date 26/01/2022
# Updated 01/02/2022


import numpy as np
import nextnanopy as nn
import matplotlib.pyplot as plt
import pandas as pd
from numpy import ndarray


class Distance:
    """Class for one array operation"""

    def __init__(self, array_input):
        self.array = np.asarray(array_input)

    def find_nearest(self, value):
        """ return element in array that is nearest to input value and also index of it,
        adopted from stackoverflow"""
        idx = (np.abs(self.array - value)).argmin()
        return self.array[idx], idx

    def get_dist(self):
        """find distance of element in an array,
        used to calculate differential element dx"""
        dist_array = []
        for i in range(self.array.size - 1):  # -1 for end-point effect, note also len(dist_array) = len(array) - 1
            dist = self.array[i + 1] - self.array[i]
            dist_array = np.append(dist_array, dist)
        return dist_array

    def dx(self):
        """to get nonuniform differential element,
        factor convert from nm to whatever"""
        return self.get_dist()


class Integrate:
    """Class for two arrays"""

    def __init__(self, array1, array2, factor):
        """factor set to be x conversion factor"""
        self.array1 = Distance(array1)  # density obj
        self.array2 = Distance(array2)  # coords obj
        self.arraydx = np.asarray(Distance(array2).dx())  # dx obj, factor is often convert_factor defined latter
        self.factor = factor

    def function_integral(self, index=None):
        """traverse over specified indices,
        think this as \int f(x) dx, but now dx is not uniform. Integrating limit given in index,
        in out case, len(dx)=len(edensity)-1 so set dx as array2, edensity as array1,
        if index = None, traverse over whole array"""
        dummy_array = []
        if index is None:  # The author thanks Dean Xiong for suggesting this way
            index = range(self.arraydx.size)
        for i in index:
            dummy_array = np.append(dummy_array, (self.array1.array[i] * self.arraydx[i]))
        integral = np.sum(dummy_array)
        return integral

    def charge_number(self, start=None, end=None):
        if start and end is not None:
            index_0 = self.array2.find_nearest(start * self.factor)[1]
            index_1 = self.array2.find_nearest(end * self.factor)[1]
            indices = [i for i in range(index_0, index_1 + 1)]  # +1 for end-point effect of range()
            return self.function_integral(indices)
        else:
            return self.function_integral()




if __name__ == '__main__':


    # import input file
    my_input = nn.InputFile(r'F:\NextNano Data\TLL\full_substrate\C2617full.in')

    itr = 2
    # # import output electron density datafile
    # edensity_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
    #                             r' - University of Cambridge\Documents\nextnano\Output'
    #                             r'\C2617full\bias_00000\density_electron.dat'.format(itr=itr),
    #                             product='nextnano++')
    #
    # donor_density_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
    #                                  r' - University of Cambridge\Documents\nextnano\Output'
    #                                  r'\C2617full\bias_00000\density_donor_ionized.dat'.format(itr=itr),
    #                                  product='nextnano++')
    #
    # acceptor_density_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
    #                                     r' - University of Cambridge\Documents\nextnano\Output'
    #                                     r'\C2617full\bias_00000\density_acceptor_ionized.dat'.format(itr=itr),
    #                                     product='nextnano++')
    # import output electron density datafile
    edensity_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                                r' - University of Cambridge\Documents\nextnano\Output'
                                r'\C2617full({itr})\bias_00000\density_electron.dat'.format(itr=itr),
                                product='nextnano++')

    donor_density_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                                     r' - University of Cambridge\Documents\nextnano\Output'
                                     r'\C2617full({itr})\bias_00000\density_donor_ionized.dat'.format(itr=itr),
                                     product='nextnano++')

    acceptor_density_file = nn.DataFile(r'E:\OneDrive--University of Cambridge\OneDrive'
                                        r' - University of Cambridge\Documents\nextnano\Output'
                                        r'\C2617full({itr})\bias_00000\density_acceptor_ionized.dat'.format(itr=itr),
                                        product='nextnano++')
    convert_factor = 1e-7  # x conversion factor
    # convert nm to cm, remeber to add it in x_coords and Integrate()
    convert_factor_nextnano = 1e18  # y conversion factor
    # density has a unit of 1e18 cm^-3 in nextnano as default, remember to add it to density

    # extract output file coordinates
    x_coords = edensity_file.coords['x'].value * convert_factor
    edensity = edensity_file.variables['Electron_density'].value * convert_factor_nextnano
    x_coords2 = donor_density_file.coords['x'].value * convert_factor
    donor_density = donor_density_file['Ionized_donor_density'].value * convert_factor_nextnano
    x_coords3 = acceptor_density_file.coords['x'].value * convert_factor
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
    QWs = 1  # = 1 if we want number in Quantum Well, otherwise can define arbitrary start and end point
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
    total_electron = Integrate(edensity, x_coords, convert_factor).charge_number()
    Region1_electron = Integrate(edensity, x_coords, convert_factor).charge_number(region1_start, region1_end)
    Region2_electron = Integrate(edensity, x_coords, convert_factor).charge_number(region2_start, region2_end)
    total_donor = Integrate(donor_density, x_coords2, convert_factor).charge_number()
    Region1_donor = Integrate(donor_density, x_coords2, convert_factor).charge_number(region1_start, region1_end)
    Region2_donor = Integrate(donor_density, x_coords2, convert_factor).charge_number(region2_start, region2_end)
    total_acceptor = Integrate(acceptor_density, x_coords3, convert_factor).charge_number()
    Region1_acceptor = Integrate(acceptor_density, x_coords3, convert_factor).charge_number(region1_start, region1_end)
    Region2_acceptor = Integrate(acceptor_density, x_coords3, convert_factor).charge_number(region2_start, region2_end)

    print("{Type} = {:.3e}".format(total_electron, Type='total_electron') + " cm^-2")
    print("{Type} = {:.3e}".format(total_donor, Type='total_donor') + " cm^-2")
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

    # test: export edensity to a csv file
    # pd.DataFrame(edensity.value).to_csv(r"E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge\Desktop\Research\NextNano\file.csv")
