# TO CALCULATE EUCLIDEAN DISTANCE BETWEEN TWO DATA SETS AND MAKE A PLOT
# AUTHOR: QI LIU
# DATE: 07/02/2022

# TODO extract Correlation coefficient between two data sets
# TODO 3D PLOT X ENERGY, Y ARRAY RATIO, Z
# TODO DETLTA E TOO SMALL CHANGE IN BOTH OUT AND INPUT FILE!!!

# TODO CHECK COMPATIBILITY
import nextnanopy as nn
import numpy as np
import matplotlib.pyplot as plt
import timeit
from Charge_density_oop import Integrate

d = {}
for i in np.linspace(-0.01, 0.03, 5):
    print(i)

