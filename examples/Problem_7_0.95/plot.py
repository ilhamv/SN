import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

phi = f['scalar_flux'];
rho = f['spectral_radius'];
z = f['z'];

plt.plot(z,phi, 'o')
plt.title("Scalar Flux");
plt.xlabel("z, cm");
plt.ylabel("Scalar flux, /cm^2.s");
plt.show();

l = np.arange(1,len(rho)+1);
plt.plot(l,rho, 'o');
plt.title("Spectral Radius");
plt.xlabel("Iteration #");
plt.ylabel("Spectral radius");
plt.show();
