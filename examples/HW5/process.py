import h5py
import matplotlib.pyplot as plt
import numpy as np
from subprocess import call

mesh_list = [5, 10, 20, 40, 100, 200]
N_list  = [2,4,8,16,32]

pt = []
N_iter = []
rho = []

args = ["./../../SN.exe","."]

for mesh in mesh_list:
    pt_temp = []
    N_iter_temp = []
    rho_temp = []
    for N in N_list:
        with open('input.xml', 'r') as file:
            data = file.readlines()
        data[3] = "<N> %i </N>\n"%N;
        data[18] = "    <mesh> %i </mesh>\n"%mesh;
        with open('input.xml', 'w') as file:
            file.writelines( data )
        call(args)
        f = h5py.File('output.h5', 'r');
        psi = np.array(f['angular_flux'])
        mu = np.array(f['mu_n'])
        w = np.array(f['w_n'])
        r = np.array(f['spectral_radius'])
        N_i = int(np.array(f['N_iter']))
        f.close()
        p = 0
        for i in range(int(N/2),N):
            p = p + psi[-1][i]* mu[i] * w[i]
        pt_temp.append(p)
        N_iter_temp.append(N_i)
        rho_temp.append(r[-1])
    pt.append(pt_temp)
    N_iter.append(N_iter_temp)
    rho.append(rho_temp)

for j in range(len(N_list)):
    out=""
    for i in range(len(mesh_list)):
        out = out + str(pt[i][j]) + " "
    print(out)
print("\n\n")
for j in range(len(N_list)):
    out=""
    for i in range(len(mesh_list)):
        out = out + str(N_iter[i][j]) + " "
    print(out)
print("\n\n")
for j in range(len(N_list)):
    out=""
    for i in range(len(mesh_list)):
        out = out + str(rho[i][j]) + " "
    print(out)
