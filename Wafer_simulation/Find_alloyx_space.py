# TO CALCULATE EUCLIDEAN DISTANCE BETWEEN TWO DATA SETS AND MAKE A PLOT
# AUTHOR: QI LIU
# DATE: 07/02/2022


# todo check biased cases
# TODO plot delta_E Vs DIFFERENCE
# TODO CHECK MANUALLY IF CALCULATION MAKES SENSE
# TODO PLOT DEVICE TRANSMISSION IN SUBPLOT AS WELL ALTHOUGH WE DO NOT CARE
# TODO TRY OTHER WAFERS
# TODO REPORT TO CHRIS AND CC PEDRO
import nextnanopy as nn
import numpy as np
import matplotlib.pyplot as plt
import timeit
from Charge_density_oop import Integrate

delta_E_limit = np.linspace(-0.01, 0.01, 12)
space = 0.05
delta_E = 0  # offset in ground state energy, could be DC bias or whatever
dict_ratio = {}
dict_ratio_diff = {}
dict_intersect_x = {}
dict_intersect_y = {}
dict_prob = {}
dict_prob_real = {}
dict_diff_min = {}
dict_alloyx_min = {}
start = timeit.default_timer()


for delta_E in delta_E_limit:

    class DEVICE_FILE:

        def __init__(self, Gamma_path, input_path, prob_path=None):
            self.input_file = nn.InputFile(input_path)
            width1 = self.input_file.variables['surface_cap_width'].value
            width2 = self.input_file.variables['upper_AlGaAs_doping_layer_width'].value
            width3 = self.input_file.variables['upper_AlGaAs_modulation_spacer_width'].value
            QW1_start = width1 + width2 + width3
            QW1_width = self.input_file.variables['QW1_width'].value
            self.edge = QW1_start + QW1_width - space  # -0.05 to allow some space
            self.Gamma_file = nn.DataFile(Gamma_path, product='nextnano++')
            self.x = self.Gamma_file.coords['x'].value
            self.Gamma = self.Gamma_file.variables['Band_Edge'].value
            if prob_path is not None:
                self.prob_file = nn.DataFile(prob_path, product='nextnano++')
                self.prob = self.prob_file.variables['1->2'].value
                self.E = self.prob_file.variables['Energy'].value

        def E_1(self):
            idx = (np.abs(self.edge - self.x)).argmin()
            return np.abs(self.Gamma[idx] + delta_E)


    class FILE:

        def __init__(self, prob_path, eigen_path, input_path=None):
            if input_path is not None:
                self.input_file = nn.InputFile(input_path)
            self.prob_file = nn.DataFile(prob_path, product='nextnano++')
            self.eigen_file = nn.DataFile(eigen_path, product='nextnano++')
            self.prob = self.prob_file.variables['1->2'].value
            self.E = self.prob_file.variables['Energy'].value
            # self.E_1 = self.eigen_file.variables['E_1'].value[0] #NOT USING THIS ANY MORE


    def prob(E1, E2, array):
        # E1 IS GROUND STATE ENERGY FROM DEVICE, BY DEFAULT SCALAR
        # E2 IS ENERGY ARRAY IN SIMPLE MODEL TRANS PLOT, VECTOR
        # array IS SIMPLE MODEL TRANS ARRAY
        idx = (np.abs(E1 - E2)).argmin()
        return array[idx]


    def find_nearest(array1, array2):
        array1 = np.asarray(array1)
        array2 = np.asarray(array2)
        idx = (np.abs(array1 - array2)).argmin()
        return array1[idx], idx


    # APPLY BIAS TO QWs?
    is_applied_bias = 0

    # SPECIFY REFERENCE EXISTING WAFER
    if is_applied_bias:
        Model = 'C2617_bias'
        # Model = 'C2618_bias'
        # Model = 'W0938_bias'
        # Model = 'W0939_bias'
        # Model = 'W0940_bias'
        # Model = 'W0941_bias'
        bias_type = '_bias'
    else:
        # Model = 'C2617'
        # Model = 'C2618'
        # Model = 'W0938'
        # Model = 'W0939'
        # Model = 'W0940'
        Model = 'W0941'
        bias_type = ''

    # STANDARD MODEL N=10 FILE, MUST NOT BE ALTERED
    MODEL10_input_path = str(r'F:\NextNano Data\TLL\Existingf=10{bias}.in'.format(bias=bias_type))
    MODEL10_prob_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                            r'Cambridge\Documents\nextnano\Output\Existingf=10{bias}'
                            r'\bias_00000\transmission_cbr_Gamma1.dat'.format(bias=bias_type))
    MODEL10_eigen_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                             r'\Documents\nextnano\Output\Existingf=10{bias}\bias_00000\Quantum'
                             r'\probabilities_shift_cbr_Gamma_00000.dat'.format(bias=bias_type))

    # Model TEST WAFER N=20 STARTING FILE
    ModTEST_input_path = str(r'F:\NextNano Data\TLL\ModTEST{bias}.in'.format(bias=bias_type))
    ModTEST_prob_path0 = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                             r'Cambridge\Documents\nextnano\Output\ModTEST{bias}'
                             r'\bias_00000\transmission_cbr_Gamma1.dat'.format(bias=bias_type))
    ModTEST_eigen_path0 = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                              r'\Documents\nextnano\Output\ModTEST\bias_00000\Quantum\probabilities_shift_cbr_Gamma_00000.dat')

    # DEVICE FILE
    # STANDARD DEVICE N=10 FILE, MUST NOT BE ALTERED

    DEV10_input_path = str(r'F:\NextNano Data\TLL\{type}full.in'.format(type=Model))
    DEV10_Gamma_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                           r'\Documents\nextnano\Output\{type}full\bias_00000\bandedge_Gamma.dat'.format(type=Model))
    DEV10_prob_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                          r'\Documents\nextnano\Output\{type}full\bias_00000\transmission_cbr_Gamma1.dat'.format(
        type=Model))

    # TESTING DEVICE N=20 STRATING FILE

    DevTEST_input_path0 = str(r'F:\NextNano Data\TLL\{type}DEVf=20.in'.format(type=Model))
    DevTEST_Gamma_path0 = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                              r'\Documents\nextnano\Output\{type}DEVf=20\bias_00000\bandedge_Gamma.dat'.format(
        type=Model))
    DevTEST_prob_path0 = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                             r'\Documents\nextnano\Output\{type}DEVf=20\bias_00000\transmission_cbr_Gamma1.dat'.format(
        type=Model))
    # NEW WAFER IN REALITY STARTING FILE (MODEL)
    REAL_input_path = str(r'F:\NextNano Data\TLL\REAL{bias}.in'.format(bias=bias_type))
    REAL_prob_path0 = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of '
                          r'Cambridge\Documents\nextnano\Output\REAL{bias}'
                          r'\bias_00000\transmission_cbr_Gamma1.dat'.format(bias=bias_type))
    REAL_eigen_path0 = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                           r'\Documents\nextnano\Output\REAL{bias}'
                           r'\bias_00000\Quantum\probabilities_shift_cbr_Gamma_00000.dat'.format(bias=bias_type))

    #
    # # edit input file
    MODEL10 = FILE(MODEL10_prob_path, MODEL10_eigen_path, MODEL10_input_path)
    ModTEST = FILE(ModTEST_prob_path0, ModTEST_eigen_path0, ModTEST_input_path)
    REAL = FILE(REAL_prob_path0, REAL_eigen_path0, REAL_input_path)
    DEV10 = DEVICE_FILE(DEV10_Gamma_path, DEV10_input_path, DEV10_prob_path)
    DevTEST = DEVICE_FILE(DevTEST_Gamma_path0, DevTEST_input_path0)

    ## FIXED, DO NOT CHANGE WITH ITERATION, USED LATTER WHEN CHANGING REAL (X=0.33)
    DevTEST0 = DEVICE_FILE(DevTEST_Gamma_path0, DevTEST_input_path0)
    ModTEST0 = FILE(ModTEST_prob_path0, ModTEST_eigen_path0, ModTEST_input_path)
    #

    #
    #
    #
    diff_array = []
    var1_name = 'alloy1_x'
    limit = np.linspace(0.11, 0.21, 21)
    isvarymodeln20 = 1

    for var in limit:
        # SWEEP PARAMETER IN SIMPLE MODEL TEST FILE
        ModTEST.input_file.set_variable(var1_name, value=var)
        # NOTE THIS WILL CHANGE ModTEST input properties although saved in a different path... checked in test!
        var1 = ModTEST.input_file.variables[var1_name].value
        ModTEST.input_file.save(r'F:\NextNano Data\TLL'
                                r'\ModTEST{bias}_{var1}.in'.format(bias=bias_type, var1=var1_name + '=' + str(var1)),
                                automkdir=False,
                                overwrite=True)
        ModTEST.input_file.execute()
        # import data file that corresponds to the input we have just set
        ModTEST_prob_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                                r'\Documents\nextnano\Output\ModTEST{bias}_{var1}'
                                r'\bias_00000\transmission_cbr_Gamma1.dat'
                                .format(bias=bias_type, var1=var1_name + '=' + str(var1)))
        ModTEST_eigen_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                                 r'\Documents\nextnano\Output\ModTEST{bias}_{var1}'
                                 r'\bias_00000\Quantum\probabilities_shift_cbr_Gamma_00000.dat'
                                 .format(bias=bias_type, var1=var1_name + '=' + str(var1)))

        ModTEST_output = FILE(ModTEST_prob_path, ModTEST_eigen_path)  # check input??

        if isvarymodeln20:
            # SWEEP PARAMETER IN DEVICE MODEL TEST FILE
            DevTEST.input_file.set_variable(var1_name, value=var)
            var1 = DevTEST.input_file.variables[var1_name].value

            DevTEST.input_file.save(r'F:\NextNano Data\TLL'
                                    r'\{type}DEVf=20_{var1}.in'.format(type=Model, var1=var1_name + '=' + str(var1)),
                                    automkdir=False,
                                    overwrite=True)
            DevTEST.input_file.execute()
            DevTEST_input_path = str(r'F:\NextNano Data\TLL'
                                     r'\{type}DEVf=20_{var1}.in'.format(type=Model, var1=var1_name + '=' + str(var1))
                                     .format(type=Model, var1=var1_name + '=' + str(var1)))
            DevTEST_Gamma_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                                     r'\Documents\nextnano\Output\{type}DEVf=20_{var1}'
                                     r'\bias_00000\bandedge_Gamma.dat'
                                     .format(type=Model, var1=var1_name + '=' + str(var1)))
            DevTEST_output = DEVICE_FILE(DevTEST_Gamma_path, DevTEST_input_path)

            prob_n10 = prob(DEV10.E_1(), MODEL10.E, MODEL10.prob)  # FIXED
            prob_n20 = prob(DevTEST_output.E_1(), ModTEST_output.E, ModTEST_output.prob)  # DevTEST_ouput.E1()
            diff = np.abs(prob_n10 - prob_n20)
            diff_array = np.append(diff_array, diff)
        else:
            prob_n10 = prob(DEV10.E_1(), MODEL10.E, MODEL10.prob)  # FIXED
            prob_n20 = prob(DevTEST.E_1(), ModTEST_output.E, ModTEST_output.prob)
            diff = np.abs(prob_n10 - prob_n20)
            diff_array = np.append(diff_array, diff)

    diff_min = np.min(diff_array)
    # RESULTANT X VALUE
    alloyx_min = limit[np.where(diff_array == diff_min)][0]

    # RE-TRACK NEW MODEL WAFER WITH RESULTANT X
    MOD_result_prob_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                               r'\Documents\nextnano\Output\ModTEST{bias}_{var1}'
                               r'\bias_00000\transmission_cbr_Gamma1.dat'
                               .format(bias=bias_type, var1=var1_name + '=' + str(alloyx_min)))
    MOD_result_eigen_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                                r'\Documents\nextnano\Output\ModTEST{bias}_{var1}'
                                r'\bias_00000\Quantum\probabilities_shift_cbr_Gamma_00000.dat'
                                .format(bias=bias_type, var1=var1_name + '=' + str(alloyx_min)))

    MOD_result = FILE(MOD_result_prob_path, MOD_result_eigen_path)
    MOD_result_path = MOD_result_prob_path
    # RE-TRACK NEW DEVICE WAFER WITH RESULTANT X

    if isvarymodeln20:
        DEV_result_Gamma_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                                    r'\Documents\nextnano\Output\{type}DEVf=20_{var1}'
                                    r'\bias_00000\bandedge_Gamma.dat'
                                    .format(type=Model, var1=var1_name + '=' + str(alloyx_min)))
        DEV_result_prob_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                                   r'\Documents\nextnano\Output\{type}DEVf=20_{var1}'
                                   r'\bias_00000\transmission_cbr_Gamma1.dat'
                                   .format(type=Model, var1=var1_name + '=' + str(alloyx_min)))
        DEV_result_input_path = str(r'F:\NextNano Data\TLL'
                                    r'\{type}DEVf=20_{var1}.in'.format(type=Model,
                                                                       var1=var1_name + '=' + str(alloyx_min)))

        DEV_result = DEVICE_FILE(DEV_result_Gamma_path, DEV_result_input_path, DEV_result_prob_path)
        Dev_result_path = DEV_result_Gamma_path
    else:  # fixed n = 20 Device E1
        DEV_result = DEVICE_FILE(DevTEST_Gamma_path0, DevTEST_input_path0, DevTEST_prob_path0)
        Dev_result_path = DevTEST_Gamma_path0
    # REAL WAFER SECTION

    # GENERATE NEW WAFER (MODEL) IN REALITY BY CHANGING WIDTH FRACTION + PERIOD LENGTH CONSTRAINT
    x_av = (alloyx_min) * 0.833 / (0.556 + 0.833)
    t_GaAs = lambda t: ((0.33 - x_av) / x_av) * t
    y = lambda x: 1.389 - x  # constraint equation
    t_AlGaAs = np.linspace(0, 1, 100)
    intersect_y = find_nearest(t_GaAs(t_AlGaAs), y(t_AlGaAs))[0]  # t_GaAs
    intersect_x = t_AlGaAs[find_nearest(t_GaAs(t_AlGaAs), y(t_AlGaAs))[1]]  # t_AlGaAs

    var2_name = 'Barrier_Width_2'  # t_AlGaAs
    var3_name = 'Barrier_Width_1'  # t_GaAs

    REAL.input_file.set_variable(var2_name, value=intersect_x)
    REAL.input_file.set_variable(var3_name, value=intersect_y)
    var2 = REAL.input_file.variables[var2_name].value
    var3 = REAL.input_file.variables[var3_name].value
    REAL.input_file.save(r'F:\NextNano Data\TLL\{type}FIN_MODREAL_{var2}_{var3}_dE{delta}.in'
                         .format(type=Model, var2=var2_name + '=' + str(round(var2, 4)),
                                 var3=var3_name + '=' + str(round(var3, 4)), delta=round(delta_E, 4)),
                         automkdir=False, overwrite=True)
    REAL.input_file.execute()

    REAL_prob_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                         r'\Documents\nextnano\Output\{type}FIN_MODREAL_{var2}_{var3}_dE{delta}'
                         r'\bias_00000\transmission_cbr_Gamma1.dat'
                         .format(type=Model, var2=var2_name + '=' + str(round(var2, 4)),
                                 var3=var3_name + '=' + str(round(var3, 4)), delta=round(delta_E, 4)))
    REAL_eigen_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                          r'\Documents\nextnano\Output\{type}FIN_MODREAL_{var2}_{var3}_dE{delta}'
                          r'\bias_00000\Quantum\probabilities_shift_cbr_Gamma_00000.dat'
                          .format(type=Model, var2=var2_name + '=' + str(round(var2, 4)),
                                  var3=var3_name + '=' + str(round(var3, 4)), delta=round(delta_E, 4)))

    MODREAL_result = FILE(REAL_prob_path, REAL_eigen_path)

    # GENERATE NEW WAFER (DEVICE) IN REALITY
    print(DevTEST0.input_file.variables[var1_name].value)
    print('-----------------------------------------------')
    DevTEST0.input_file.set_variable(var2_name, value=intersect_x)
    DevTEST0.input_file.set_variable(var3_name, value=intersect_y)
    var2 = DevTEST0.input_file.variables[var2_name].value
    var3 = DevTEST0.input_file.variables[var3_name].value
    DevTEST0.input_file.save(r'F:\NextNano Data\TLL\{Type}FIN_DEVREAL_{var2}_{var3}_dE{delta}.in'
                             .format(Type=Model, var2=var2_name + '=' + str(round(var2, 4)),
                                     var3=var3_name + '=' + str(round(var3, 4)), delta=round(delta_E, 4))
                             , automkdir=False, overwrite=True)
    DevTEST0.input_file.execute()
    DevTEST0_input_path = str(r'F:\NextNano Data\TLL\{Type}FIN_DEVREAL_{var2}_{var3}_dE{delta}.in'
                              .format(Type=Model, var2=var2_name + '=' + str(round(var2, 4)),
                                      var3=var3_name + '=' + str(round(var3, 4)), delta=round(delta_E, 4)))
    DevTEST0_Gamma_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                              r'\Documents\nextnano\Output\{type}FIN_DEVREAL_{var2}_{var3}_dE{delta}'
                              r'\bias_00000\bandedge_Gamma.dat'
                              .format(type=Model, var2=var2_name + '=' + str(round(var2, 4)),
                                      var3=var3_name + '=' + str(round(var3, 4)), delta=round(delta_E, 4)))
    DevTEST0_prob_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge'
                             r'\Documents\nextnano\Output\{type}FIN_DEVREAL_{var2}_{var3}_dE{delta}'
                             r'\bias_00000\transmission_cbr_Gamma1.dat'
                             .format(type=Model, var2=var2_name + '=' + str(round(var2, 4)),
                                     var3=var3_name + '=' + str(round(var3, 4)), delta=round(delta_E, 4)))
    DEVREAL_result = DEVICE_FILE(DevTEST0_Gamma_path, DevTEST0_input_path, DevTEST0_prob_path)

    prob_wafer = prob(DEV_result.E_1(), MOD_result.E, MOD_result.prob)
    prob_wafer_real = prob(DEVREAL_result.E_1(), MODREAL_result.E, MODREAL_result.prob)

    dict_ratio["delta_E=" + str(delta_E)] = limit
    dict_ratio_diff["delta_E=" + str(delta_E)] = diff_array
    dict_intersect_x["delta_E=" + str(delta_E)] = intersect_x
    dict_intersect_y["delta_E=" + str(delta_E)] = intersect_y
    dict_prob["delta_E=" + str(delta_E)] = prob_wafer
    dict_prob_real["delta_E=" + str(delta_E)] = prob_wafer_real
    dict_diff_min["delta_E=" + str(delta_E)] = diff_min
    dict_alloyx_min["delta_E=" + str(delta_E)] = alloyx_min

    print("---FINISHED CALCULATING---")
    print('x={x}'.format(x=alloyx_min))
    print('x_av={av}'.format(av=x_av))
    print("(Delta_E = {delta})Transmission probability of new wafer is {x:.3e} ".format(x=prob_wafer, delta=delta_E))
    print("(Delta_E = {delta})Transmission probability of new wafer (REAL) is {x:.3e} ".format(x=prob_wafer_real,
                                                                                               delta=delta_E))
    print("PLOTTING GRAPHS...")

    # PLOTTING SECTION
stop = timeit.default_timer()

print('Calculating Time:{time} seconds'.format(time=stop - start))

print(dict_ratio)
print(dict_ratio_diff)
print(dict_intersect_x)
print(dict_intersect_y)
print(dict_prob)
print(dict_prob_real)

fig, ax = plt.subplots(2, 2)

prob = []
prob_real = []
for i in delta_E_limit:
    ax[0, 0].plot(dict_ratio[r'delta_E={x}'.format(x=i)], dict_ratio_diff[r'delta_E={x}'.format(x=i)],
                  label='delta_E={d}'.format(d=round(i, 4)))
    ax[0, 0].scatter(dict_alloyx_min[r'delta_E={x}'.format(x=i)], dict_diff_min[r'delta_E={x}'.format(x=i)],
                     label='Minimum diff ({x},{:.3e})'.format(dict_diff_min[r'delta_E={x}'.format(x=i)],
                                                              x=round(dict_alloyx_min[r'delta_E={x}'.format(x=i)], 4)))
    ax[0, 1].scatter(i, dict_intersect_x[r'delta_E={x}'.format(x=i)])

    ax[1, 0].scatter(i, dict_intersect_y[r'delta_E={x}'.format(x=i)])
    prob = np.append(prob, dict_prob[r'delta_E={x}'.format(x=i)])
    prob_real = np.append(prob_real, dict_prob_real[r'delta_E={x}'.format(x=i)])


ax[0, 0].legend(loc='upper right')
ax[0, 0].set_xlabel("Alloy ratio")
ax[0, 0].set_ylabel("Abs diff in ground state E trans prob")
ax[0, 0].set_title("{type}".format(type=Model))
ax[0, 1].set_xlabel('delta_E')
ax[0, 1].set_ylabel('intersect_x (t_=AlGaAs)')

ax[1, 0].set_xlabel('delta_E')
ax[1, 0].set_ylabel('intersect_y (t_=GaAs)')

ax[1, 1].plot(delta_E_limit, prob, label='probability')
ax[1, 1].plot(delta_E_limit, prob_real, label='probability(REAL)')
ax[1, 1].set_xlabel('delta_E')
ax[1, 1].set_ylabel('Transmission probability')
ax[1, 1].legend(loc='upper right')

fig.set_size_inches(20, 20)
plt.tight_layout()
plt.show()

fig_path = str(r'E:\OneDrive--University of Cambridge'
               r'\OneDrive - University of Cambridge\Desktop\Research\MBE Wafer Growth'
               r'\Python figs\{Type}_model_dE(2).png'.format(Type=Model))

fig.savefig(fig_path, dpi=500)
