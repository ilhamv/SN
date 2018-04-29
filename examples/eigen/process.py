from subprocess import call
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Cases run
for zeta in [0.0, 0.9, 0.99]:
    Z = [5., 10., 20., 40., 80.]
    I = len(Z)
    rho = np.zeros(I)
    N = np.zeros(I)
    lamb = np.zeros(I)

    args = ["./../../SN.exe","."]

    for i in range(I):
        with open('input.xml', 'r') as file:
            data = file.readlines()
        data[8] = "    <shift param=\"0.0 0.0 %f\"/>\n"%zeta
        data[21] = "    <space> %f </space>\n"%Z[i]
        data[22] = "    <mesh> %i </mesh>\n"%Z[i]*10
        with open('input.xml', 'w') as file:
            file.writelines( data )
        call(args)
        
        f = h5py.File('output.h5', 'r');
        rho[i] = np.array(f['spectral_radius'])[-1]
        lamb[i] = np.array(f['lambda'])[-1]
        N[i] = np.array(f['N_iter'])
        f.close()

    for i in range(I):
        print(Z[i],N[i],rho[i],lamb[i])

    plt.figure(1)
    plt.plot(Z,N,"o",label="Code")
    plt.xlabel('Z')
    plt.ylabel('# of iterations')
    plt.grid()

    plt.figure(2)
    plt.plot(Z,rho,"o",label="Code")
    plt.xlabel('Z')
    plt.ylabel(r'$\rho$')
    plt.grid()

    plt.figure(3)
    plt.plot(Z,lamb,"o",label="Code")
    plt.xlabel('Z')
    plt.ylabel(r'$\lambda$')
    plt.grid()


plt.figure(1)
plt.show()
plt.figure(2)
plt.show()
plt.figure(3)
plt.show()
