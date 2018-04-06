import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

phi = f['scalar_flux'];
rho = f['spectral_radius'];
z = np.array(f['z']);

Q = 0.1
Z = 10.0
SigmaT = 1.0
D = 1.0 / 3.0 / SigmaT
C1 = (Q/2/D*Z**2+2*Q*Z)/(Z+4*D)
C2 = 2*D*C1
def phi_diffusion(z):
    return -Q/2/D*z**2 + C1*z + C2

plt.plot(z,phi, 'o',label="Code")
plt.plot(z,phi_diffusion(z),label="Analytic Diffusion")
plt.legend()
plt.title("Scalar Flux");
plt.xlabel("z, cm");
plt.ylabel("Scalar flux, /cm^2.s");
plt.show();


