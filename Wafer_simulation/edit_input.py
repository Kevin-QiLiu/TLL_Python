# edit and sweep NextNano input files
# Author: Qi Liu
# Date: 25/01/2022

# TODO extract Correlation coefficient between two data sets

import nextnanopy as nn
import numpy as np
import time
import matplotlib.pyplot as plt
# edit input file
my_input = nn.InputFile(r'F:\NextNano Data\TLL\C2617_adoptedyiqing.in')
print(f"List of variables: {my_input.variables}")

var1_name = 'alloy1_x'
# my_input.set_variable(var1_name, value=var)
# Remeber to delete duplicate files after this step
# Sweep variable
dist_array = []
for var in np.linspace(0.4, 0.5, 11):  # remember the 3rd argument needs to +1
    my_input.set_variable(var1_name, value=var)
    var1 = my_input.variables[var1_name].value
    my_input.save(r'F:\NextNano Data\TLL'
                  r'\C2617_adoptedyiqing_{var1}.in'.format(var1=var1_name + '=' + str(var1)), automkdir=False,
                  overwrite=True)
    my_input.execute()
    print(f"New variable: {my_input.get_variable('{inputvar}'.format(inputvar=var1_name)).text}")
