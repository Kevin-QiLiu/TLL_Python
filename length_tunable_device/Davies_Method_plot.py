from Davies_SingleGate_Gen import GATE, slice_2dMatrix
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# Define physical constants in SI units
fac_nano = 1e-9
epsilon_0 = 8.854e-12  # unit: Fm^-1
epsilon_r = 12.9  # dimensionless, for GaAs qw
epsilon = epsilon_r * epsilon_0
q = 1.6e-19  # unit: C
depth = 7.66e-8  # unit: m
n_top2DEG = 2.85e15  # unit: m^-2
n_bot2DEG = 1.54e15  # unit: m^-2
V_t = -(q * n_top2DEG * depth) / epsilon
E_F = 0.045  # unit: eV



if __name__ == '__main__':
    # Specify gate geometries
    # Note: analytical solution in Davie's paper is dimensionless, so no need to input SI units
    # Just need to all length units CONSISTENT!
    # We choose nanometers here
    d = depth * (10 ** 9)
    SG_Le = -250
    SG_Ri = 250
    SG_Bo = -200
    SG_To = 200
    num_SG = 6  # number of SG
    MG_Le = -1500
    MG_Ri = SG_Ri
    MG_Bo = -50
    MG_To = 50

    gate_space = -530
    gate_fac = 0.7

    V_SG = V_t * 0.6
    V_gate_1 = V_t*0.6
    V_gate_2 = V_t*0.6
    V_MG = 0

    x_min = -500
    x_max = 3200.0
    y_min = -500
    y_max = 500

    x = np.linspace(x_min, x_max, 700)
    y = np.linspace(y_min, y_max, 700)
    X, Y = np.meshgrid(x, y)

    MidGate = GATE.MG(X, Y, d, MG_Le, MG_Ri, MG_Bo, MG_To, V_MG)
    potential = GATE.SG_array(X, Y, d, SG_Ri, SG_To, num_SG, V_SG, V_gate_1, V_gate_2) + MidGate

    vslice_pos = 1700
    hslice_pos = 0
    vslice = slice_2dMatrix(potential, vslice_pos, x, y)
    hslice = slice_2dMatrix(potential, hslice_pos, y, x)
    vslice_lower_limit = SG_Bo - 50
    vslice_upper_limit = SG_To + 50
    hslice_lower_limit = -1000
    hslice_upper_limit = 4000

    V_fermi = vslice.get_Ef_intersect()[0]
    intersec_V = vslice.get_Ef_intersect()[1]
    y_intersec = vslice.get_Ef_intersect()[2]
    H_fermi = hslice.get_Ef_intersect()[0]
    intersec_H = hslice.get_Ef_intersect()[3]
    x_intersec = hslice.get_Ef_intersect()[4]



    if len(intersec_V) == 0:
        print('Potential vertical slice dose not intersect fermi energy!')
    else:
        print('Potential vertical slice intersects fermi energy at ({x},{y}) '.format(x=y_intersec, y=intersec_V))
    if len(intersec_H) == 0:
        print('Horizontal vertical slice dose not intersect fermi energy!')
    else:
        print('Potential horizontal slice intersects fermi energy at ({x},{y}) '.format(x=x_intersec, y=intersec_H))


    fig = plt.figure(figsize=(25, 25))
    # plt.title("QCP_{num}F device electrostatic potential landscape at z = top_2DEG".format(num=num_SG), fontsize=25)
    ax1 = fig.add_subplot(2, 2, 1, aspect='equal')
    img = ax1.imshow(potential, extent=[x_min, x_max, y_min, y_max], cmap=cm.hot,
                     origin='lower')
    ax1.vlines(vslice_pos, *ax1.get_ylim(), color='b', linestyle='-')
    ax1.hlines(hslice_pos, *ax1.get_xlim(), color='g', linestyle='--')
    ax1.set_xlabel('x (nm)')
    ax1.set_ylabel('y (nm)')

    cax = fig.add_axes([ax1.get_position().x1 + 0.01, ax1.get_position().y0, 0.02, ax1.get_position().height])
    cbar = plt.colorbar(img, cax=cax)  # Similar to fig.colorbar(im, cax = cax)

    cbar.set_label('Potential(eV)')
    img.set_clim(-0.1, 0.32)

    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    map_3D = ax2.plot_surface(X, Y, potential, rstride=1, cstride=1,
                              cmap='jet', edgecolor='none')

    map_3D.set_clim(-0.1, 0.35)
    ax2.set_xlabel("x (nm)")
    ax2.set_ylabel("y (nm)")
    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(y_min, y_max)
    cax2 = fig.add_axes([ax2.get_position().x1 + 0.01, ax2.get_position().y0, 0.02, ax2.get_position().height])
    cbar2 = plt.colorbar(map_3D, cax=cax2)
    cbar2.set_label('Potential (eV)')
    ax2.view_init(50, 185)  # Similar to fig.colorbar(im, cax = cax)
    ax2.axis('off')

    ax3 = fig.add_subplot(2, 2, 3)
    vslice.plot_vslice()
    ax3.set_xlim(vslice_lower_limit, vslice_upper_limit)
    ax3.plot(vslice.slice_axis, V_fermi, 'r--')
    ax3.plot(y_intersec, intersec_V, 'ro')
    ax3.set_xlabel('y (nm)')
    ax3.set_ylabel('Potential(eV)')

    ax4 = fig.add_subplot(2, 2, 4)
    hslice.plot_hslice()
    ax4.set_xlim(hslice_lower_limit, hslice_upper_limit)
    ax4.hlines(E_F, *ax4.get_xlim(), color='r', linestyle='--')
    ax4.plot(hslice.slice_axis, H_fermi, 'r--')
    ax4.plot(x_intersec, intersec_H, 'ro')
    ax4.set_xlabel('x (nm)')
    ax4.set_ylabel('Potential(eV)')

    fig_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge\Desktop\Research\single-wire device\figs\SG_array_test2.png')
    fig.savefig(fig_path, dpi=250, bbox_inches='tight')
    plt.show()
