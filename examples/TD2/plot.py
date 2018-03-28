import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

phi = f['scalar_flux_time'];
t = f['time'];
z = f['z'];

f = h5py.File('steady.h5', 'r');
phi0 = f['scalar_flux']

for k in [0,1,2,3,4,5,6,8,10,12,15,20,25]:
    plt.plot(z,phi[k], 'o')
    plt.plot(z,phi0)
    plt.title("t = %s s"%(t[k]));
    plt.xlabel("z, cm");
    plt.ylabel("Scalar flux, /cm^2.s");
    plt.show();
