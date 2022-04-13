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
E_F = 0.045  # unit: eV, import from NextNano++ C2617 wafer simulation
gate_space = -550  # unit: nm
gate_fac = 0.7


class GATE(object):
    """This is a gate class"""

    def __init__(self, x_mesh, y_mesh, z_depth, V_g, dx, dy, width, length):
        self.x_0 = x_mesh
        self.y_0 = y_mesh  # region of gate #1
        self.dx = dx  # unit delta x
        self.dy = dy  # unit delta y
        self.x = self.x_0 + self.dx  # Function translation in x
        self.y = self.y_0 + self.dy  # Function translation in y
        self.d = z_depth  # quantum well depth
        self.V_g = V_g
        self.width = width
        self.length = length

    def g(self, u, v):
        R = np.sqrt(u ** 2 + v ** 2 + self.d ** 2)
        return 1 / (2 * np.pi) * np.arctan2((u * v), (self.d * R))

    def Phi_rectangle(self, L, R, B, T):
        return self.V_g * (self.g(self.x - L, self.y - B) + self.g(self.x - L, T - self.y)
                           + self.g(R - self.x, self.y - B) + self.g(R - self.x, T - self.y))

    def Phi_stripe(self, a):  # a is dimension of stripe
        return self.V_g * 1 / np.pi * (np.arctan2((a + self.x), self.d) + np.arctan2((a - self.x), self.d))

    def Phi_split(self):
        """Davie's method"""
        return self.Phi_stripe(self.width) - self.Phi_rectangle(-self.width, self.width, -self.length, self.length)

    @classmethod
    def Phi_split_channels_single(cls, num, y_unit, X_mesh, Y_mesh, dist, Right, Top, V_SG, gate_space_temp):
        """num is the number of conducting channels, for symmetry num is odd and y_unit is separation of channels, """
        global potential_rec
        num_loop = (num - 1) // 2
        gates_temp = []

        for i in range(-1 * (num_loop), (num_loop) + 1):
            y_space = float(i) * y_unit
            rectangle_temp = cls(X_mesh, Y_mesh, dist, V_SG, gate_space_temp, y_space, Right, Top)
            gates_temp.append(rectangle_temp.Phi_rectangle(-Right, Right, -Top, Top))
            potential_rec = sum(gates_temp)
        split_single = (-1) * (
                cls(X_mesh, Y_mesh, dist, V_SG, gate_space_temp, 0, Right, Top).Phi_stripe(Right) - potential_rec)
        # __init__(self, x_mesh, y_mesh, z_depth, V_g, dx, dy, width, length):
        return split_single

    @classmethod
    def SG_array_channels(cls, X_mesh, Y_mesh, dist, Right, Top, num, V_SG, V_start, V_end,
                          numofchannel, y_unit, p_width, V_PG):
        """def Phi_split_channels_single(cls, num, y_unit, X_mesh, Y_mesh, dist, Right, Top, V_SG):"""
        global potential_gates
        gates = []
        if num == 0:
            print('gates number must be an integer greater than 0')
        for i in range(num):
            gate_space_temp = i * gate_space
            if i == 0:  # specify voltage on first gate
                V_gate = V_start
                GATES_temp = cls(X_mesh, Y_mesh, dist, V_gate, gate_space_temp, 0, Right, Top)
                gates.append(GATES_temp.Phi_split() * (-1))
            elif i == 1:
                V_gate = V_PG  # PG
                GATES_temp = cls(X_mesh, Y_mesh, dist, V_gate, gate_space_temp, 0, Right, Top)
                gates.append(GATES_temp.Phi_stripe(p_width) * (-1))
            elif i == num - 1:
                V_gate = V_end  # on the last gate
                GATES_temp = cls(X_mesh, Y_mesh, dist, V_gate, gate_space_temp, 0, Right, Top)
                gates.append(
                    GATES_temp.Phi_split_channels_single(numofchannel, y_unit, X_mesh, Y_mesh, dist, Right, Top, V_gate,
                                                         gate_space_temp))
            else:
                V_gate = V_SG  # for all others SG
                GATES_temp = cls(X_mesh, Y_mesh, dist, V_gate, gate_space_temp, 0, Right, Top)
                gates.append(
                    GATES_temp.Phi_split_channels_single(numofchannel, y_unit, X_mesh, Y_mesh, dist, Right, Top, V_gate,
                                                         gate_space_temp))
            # (self.Phi_split())
            potential_gates = sum(gates)
        return potential_gates

    @classmethod
    def SG_array(cls, X_mesh, Y_mesh, dist, Right, Top, num, V_SG, V_start, V_end):
        global potential_gates
        gates = []
        if num == 0:
            print('gates number must be an integer greater than 0')
        for i in range(num):
            gate_space_temp = i * gate_space
            if i == 0:  # specify voltage on first gate
                V_gate = V_start
            elif i == num - 1:
                V_gate = V_end  # on the last gate
            else:
                V_gate = V_SG  # for all others SG
            # (self.Phi_split())
            GATES_temp = cls(X_mesh, Y_mesh, dist, V_gate, gate_space_temp, 0, Right, Top)
            gates.append(GATES_temp.Phi_split())
            potential_gates = sum(gates) * (-1)  # in unit eV
        return potential_gates

    @classmethod
    def MG(cls, X_mesh, Y_mesh, dist, Left, Right, Bottom, Top, V_MG):
        MidGate = cls(X_mesh, Y_mesh, dist, V_MG, 0, 0, 0, 0)
        return (MidGate.Phi_rectangle(Left, Right, Bottom, Top)) * (-1)

class slice_2dMatrix(object):
    def __init__(self, matrix, slice_pos, edge_axis, slice_axis):
        self.mat = matrix
        self.slice_pos = slice_pos
        self.edge_axis = edge_axis
        self.slice_axis = slice_axis

    def get_slice(self):
        idx = (np.abs(self.slice_pos - self.edge_axis)).argmin()
        return self.mat[:, idx], self.mat[idx, :]

    def plot_vslice(self):
        plt.plot(self.slice_axis, self.get_slice()[0], 'b-')

    def plot_hslice(self):
        plt.plot(self.slice_axis, self.get_slice()[1], 'g--')

    def get_Ef_intersect(self):
        E_f_line = [E_F] * len(self.slice_axis)
        v_idx = np.argwhere(np.diff(np.sign(self.get_slice()[0] - E_f_line))).flatten()
        h_idx = np.argwhere(np.diff(np.sign(self.get_slice()[1] - E_f_line))).flatten()
        v_intersec = self.get_slice()[0][v_idx]
        v_intersec_y = self.slice_axis[v_idx]
        h_intersec = self.get_slice()[1][h_idx]
        h_intersec_x = self.slice_axis[h_idx]
        return E_f_line, v_intersec, v_intersec_y, h_intersec, h_intersec_x


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

    # gate_space = -530
    # gate_fac = 0.7
    space_y = 550
    V_SG = V_t
    V_gate_1 = V_t
    V_gate_2 = V_t
    V_MG = 0.3
    #
    x_min = -3500
    x_max = 5000.0
    y_min = -7500
    y_max = 7500

    x = np.linspace(x_min, x_max, 700)
    y = np.linspace(y_min, y_max, 700)
    X, Y = np.meshgrid(x, y)
    #
    # MidGate = GATE.MG(X, Y, d, MG_Le, MG_Ri, MG_Bo, MG_To, V_MG)
    # Phi_split_channels_single(cls, num, y_unit, X_mesh, Y_mesh, dist, Right, Top, V_SG):
    # def SG_array_channels(cls, X_mesh, Y_mesh, dist, Right, Top, num, V_SG, V_start, V_end, numofchannel, y_unit):
    # potential = GATE.Phi_split_channels_single(5, space_y, X, Y, d, SG_Ri, SG_To, V_SG, -350)
    MidGate = GATE.MG(X, Y, d, MG_Le, MG_Ri, MG_Bo, MG_To, V_MG)
    potential = GATE.SG_array_channels(X, Y, d, SG_Ri, SG_To, num_SG, V_SG, V_gate_1, V_gate_2, 51, space_y, SG_Ri, 0) + MidGate
    # potential = GATE.Phi_split_channels_single(11, space_y, X, Y, d, SG_Ri, SG_To, V_SG)

    vslice_pos = 1700
    hslice_pos = 0
    vslice = slice_2dMatrix(potential, vslice_pos, x, y)
    hslice = slice_2dMatrix(potential, hslice_pos, y, x)
    vslice_lower_limit = SG_Bo - 3500
    vslice_upper_limit = SG_To + 3500
    hslice_lower_limit = -1000
    hslice_upper_limit = 4000
    #
    V_fermi = vslice.get_Ef_intersect()[0]
    intersec_V = vslice.get_Ef_intersect()[1]
    y_intersec = vslice.get_Ef_intersect()[2]
    H_fermi = hslice.get_Ef_intersect()[0]
    intersec_H = hslice.get_Ef_intersect()[3]
    x_intersec = hslice.get_Ef_intersect()[4]
    #
    #
    #
    if len(intersec_V) == 0:
        print('Potential vertical slice dose not intersect fermi energy!')
    else:
        print('Potential vertical slice intersects fermi energy at ({x},{y}) '.format(x=y_intersec, y=intersec_V))
    if len(intersec_H) == 0:
        print('Horizontal vertical slice dose not intersect fermi energy!')
    else:
        print('Potential horizontal slice intersects fermi energy at ({x},{y}) '.format(x=x_intersec, y=intersec_H))


    fig = plt.figure(figsize=(50, 50))
    # # plt.title("QCP_{num}F device electrostatic potential landscape at z = top_2DEG".format(num=num_SG), fontsize=25)
    ax1 = fig.add_subplot(1, 3, 1, aspect='equal')
    img = ax1.imshow(potential, extent=[x_min, x_max, y_min, y_max], cmap=cm.jet,
                     origin='lower')
    ax1.vlines(vslice_pos, *ax1.get_ylim(), color='b', linestyle='-')
    ax1.hlines(hslice_pos, *ax1.get_xlim(), color='g', linestyle='--')
    ax1.set_xlabel('x (nm)')
    ax1.set_ylabel('y (nm)')

    cax = fig.add_axes([ax1.get_position().x1 + 0.01, ax1.get_position().y0, 0.02, ax1.get_position().height])
    cbar = plt.colorbar(img, cax=cax)  # Similar to fig.colorbar(im, cax = cax)

    cbar.set_label('Potential(eV)')
    img.set_clim(-0.1, 0.32)
    #
    # ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    # map_3D = ax2.plot_surface(X, Y, potential, rstride=1, cstride=1,
    #                           cmap='jet', edgecolor='none')

    # map_3D.set_clim(-0.1, 0.35)
    # ax2.set_xlabel("x (nm)")
    # ax2.set_ylabel("y (nm)")
    # ax2.set_xlim(x_min, x_max)
    # ax2.set_ylim(y_min, y_max)
    # cax2 = fig.add_axes([ax2.get_position().x1 + 0.01, ax2.get_position().y0, 0.02, ax2.get_position().height])
    # cbar2 = plt.colorbar(map_3D, cax=cax2)
    # cbar2.set_label('Potential (eV)')
    # ax2.view_init(50, 185)  # Similar to fig.colorbar(im, cax = cax)
    # ax2.axis('off')

    ax3 = fig.add_subplot(1, 3, 2)
    vslice.plot_vslice()
    ax3.set_xlim(vslice_lower_limit, vslice_upper_limit)
    ax3.plot(vslice.slice_axis, V_fermi, 'r--')
    ax3.plot(y_intersec, intersec_V, 'ro')
    ax3.set_xlabel('y (nm)')
    ax3.set_ylabel('Potential(eV)')

    ax4 = fig.add_subplot(1, 3, 3)
    hslice.plot_hslice()
    ax4.set_xlim(hslice_lower_limit, hslice_upper_limit)
    ax4.hlines(E_F, *ax4.get_xlim(), color='r', linestyle='--')
    ax4.plot(hslice.slice_axis, H_fermi, 'r--')
    ax4.plot(x_intersec, intersec_H, 'ro')
    ax4.set_xlabel('x (nm)')
    ax4.set_ylabel('Potential(eV)')

    fig_path = str(r'E:\OneDrive--University of Cambridge\OneDrive - University of Cambridge\Desktop\Research\single-wire device\figs\SG_array_test3.png')
    fig.savefig(fig_path, dpi=250, bbox_inches='tight')
    plt.tight_layout()
    plt.show()
