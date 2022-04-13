# def SG_array(self, num, V_start, V_end):
#     gates = []
#     for i in range(num):
#         self.dx = i * self.dx
#         if i == 0:
#             self.V_g = V_start
#         elif i == num:
#             self.V_g = V_end
#         else:
#             self.V_g = self.V_fixed
#         gates.append(self.Phi_split())
#     potential_gates = sum(gates) * (-1)  # potential energy in unit eV
#     return potential_gates
num = 3
for i in range(num+1):
    if i == num:
        print('break')
    print(i)
